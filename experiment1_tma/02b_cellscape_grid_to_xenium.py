#!/usr/bin/env python
"""02b_cellscape_grid_to_xenium.py -- put a TMA core grid made on the CellScape scan onto the Xenium cells.

The pathology group dearrayed the CellScape scans in QuPath and exported one table
per slide (<slide>_TMA_cores.csv: core, sector, missing, centroid_x_um / centroid_y_um,
diameter_um). Those coordinates are in the CellScape image (0.18 um pixels, file
orientation); the Xenium section is a different (serial) section on another scanner,
so the grid has to be moved onto it first. This script does that without landmarks:

  1. keep the cores of the sector(s) that are on the Xenium slide (e.g. sector A)
  2. find the orientation: all 8 combinations of 0/90/180/270 degree rotation with and
     without mirroring, plus small tilts (+-4 deg) and scale (+-6 %), each scored by how
     much Xenium tissue (cell density) falls inside the core discs; the best one wins
  3. refine with an affine fit: every core's centre is re-measured on the Xenium cells
     (mean shift), outliers dropped, repeated until stable (-> the transform in the JSON)
  4. follow the section's local warping (serial sections are not identical): every Xenium
     cell goes to its nearest core (so neighbouring cores share no cells), each core moves
     to the centre of its cells, and is kept consistent with a local affine fit to its
     ~12 neighbours (max --max-local-shift x core spacing away from it). Repeated until
     stable. Cores without tissue on Xenium sit where their neighbours put them.
  5. per core: number of Xenium cells -> missing on Xenium, next to the CellScape
     "missing" flag (both / xenium_only / cellscape_only / neither)

Outputs (default folder: <scripts>/config/core_grid/, read by 02_tma_dearray.R)
  <slide_id>.csv               grid for 02 (QuPath mode): core_name, x_um, y_um, diameter_um,
                               missing (= no tissue on Xenium) + cellscape_core, sector,
                               cellscape_missing, n_cells_xenium, shift_um, presence
  <slide_id>_cellscape_to_xenium.json   the transform (CellScape um / px -> Xenium um), also
                               usable to bring CellScape cells into the Xenium frame
  <slide_id>_grid_on_xenium.png         check plot: the grid over the Xenium cells
  <slide_id>_grid_qupath.geojson        the grid as QuPath annotations (one named circle per core,
                               class = presence) on the Xenium morphology image, for manual
                               corrections: open morphology_focus in QuPath, File > Import
                               objects (or drag the file in), move / resize / re-class cores,
                               then run qupath/export_core_edits.groovy -> <slide_id>_edited.csv.
                               02 uses <slide_id>_edited.csv instead of <slide_id>.csv when it
                               exists, so re-running 02b never overwrites manual edits.

core_name drops the sector prefix ("A_1-a" -> "1-a") so 02 matches it to tma_map.csv
by row and column like any QuPath grid; cellscape_core keeps the full name, the
`tma_core` of the CellScape per-cell tables.

Example
  python 02b_cellscape_grid_to_xenium.py --slide-id TMA_ORGAN --sector A \\
      --cores-csv FL0027_TMA_cores.csv --xenium-dir /path/to/output-XETG...
"""

import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import ndimage, signal
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True, help="Xenium slide_id (slides.csv), e.g. TMA_ORGAN")
    p.add_argument("--cores-csv", required=True, help="<slide>_TMA_cores.csv of the CellScape dearray")
    p.add_argument("--xenium-dir", required=True, help="Xenium output folder (output-XETG...) of this slide")
    p.add_argument("--sector", nargs="+", required=True,
                   help="sector letter(s) on the Xenium slide: A for TMA_ORGAN (FL0027), C for TMA_TUMOR (FL0026); "
                        "in the notebook as a list, e.g. ['C']")
    p.add_argument("--out-dir", default=None, help="default: <scripts>/config/core_grid")
    p.add_argument("--bin-um", type=float, default=40.0, help="pixel size of the density image for the search")
    p.add_argument("--max-tilt", type=float, default=4.0, help="largest extra rotation tried, degrees")
    p.add_argument("--max-scale-change", type=float, default=0.06, help="largest scale change tried (0.06 = 6 %%)")
    p.add_argument("--min-cells", type=int, default=30,
                   help="cores with fewer Xenium cells inside are 'missing' on Xenium")
    p.add_argument("--max-local-shift", type=float, default=0.35,
                   help="a core may sit up to this x core spacing away from where its neighbours put it (0 = rigid)")
    p.add_argument("--diameter-um", type=float, default=None,
                   help="core diameter on Xenium; default: median CellScape diameter x scale (one size for all cores)")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------
def read_cellscape_cores(path, sectors):
    g = pd.read_csv(path)
    need = {"core", "sector", "centroid_x_um", "centroid_y_um", "diameter_um"}
    if not need <= set(g.columns):
        raise ValueError(f"{path}: needs columns {sorted(need)}; found {list(g.columns)}")
    g["sector"] = g["sector"].astype(str)
    have = sorted(g["sector"].unique())
    bad = [s for s in sectors if s not in have]
    if bad:
        raise ValueError(f"sector(s) {bad} not in {os.path.basename(path)} (it has {have})")
    g = g[g["sector"].isin(sectors)].reset_index(drop=True)
    g["cellscape_missing"] = (g["missing"].astype(str).str.lower().isin(["1", "true", "t", "yes"])
                              if "missing" in g else False)
    print(f"{len(g)} cores in sector {'+'.join(sectors)} of {os.path.basename(path)} "
          f"({int(g['cellscape_missing'].sum())} flagged missing on CellScape)")
    return g


def read_xenium_cells(xenium_dir):
    for name in ("cells.parquet", "cells.csv.gz"):
        f = os.path.join(xenium_dir, name)
        if os.path.isfile(f):
            cols = ["x_centroid", "y_centroid"]
            c = pd.read_parquet(f, columns=cols) if name.endswith(".parquet") else pd.read_csv(f, usecols=cols)
            xy = c.to_numpy(dtype=np.float64)
            print(f"{len(xy):,} Xenium cells from {f}")
            return xy
    raise FileNotFoundError(f"no cells.parquet / cells.csv.gz in {xenium_dir}")


# ---------------------------------------------------------------------------
# search: orientation x tilt x scale, translation by cross-correlation
# ---------------------------------------------------------------------------
def orientations():
    """The 8 axis-aligned orientations: (name, 2x2 matrix acting on column vectors)."""
    out = []
    for mirror in (False, True):
        for k in range(4):
            a = math.radians(90 * k)
            r = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]]).round(12)
            m = r @ np.diag([-1.0, 1.0]) if mirror else r
            out.append((f"rot{90 * k}{'+mirror' if mirror else ''}", m))
    return out


def rot(deg):
    a = math.radians(deg)
    return np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])


def density_image(xy, b, pad_um):
    """Log cell density, zero mean, padded by pad_um of 'glass' on every side (so a grid placed
    partly off the slide is not favoured over one on the tissue)."""
    x0, y0 = xy.min(0) - 2 * b
    nx = int(math.ceil((xy[:, 0].max() - x0) / b)) + 3
    ny = int(math.ceil((xy[:, 1].max() - y0) / b)) + 3
    h, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[ny, nx],
                             range=[[y0, y0 + ny * b], [x0, x0 + nx * b]])   # rows = y, columns = x
    t = np.log1p(ndimage.gaussian_filter(h, 1.0))
    t /= max(np.percentile(t[t > 0], 99), 1e-9) if (t > 0).any() else 1.0
    t = np.clip(t, 0, 1)
    t -= t.mean()                                  # zero mean: discs on empty glass cost score
    pad = int(math.ceil(pad_um / b))
    t = np.pad(t, pad, constant_values=t.min())
    return t, (x0 - pad * b, y0 - pad * b)


def disc_template(pts, radii, weights, b):
    t0 = pts.min(0) - radii.max() - b
    nx = int(math.ceil((pts[:, 0].max() + radii.max() - t0[0]) / b)) + 2
    ny = int(math.ceil((pts[:, 1].max() + radii.max() - t0[1]) / b)) + 2
    tpl = np.zeros((ny, nx))
    yy, xx = np.mgrid[0:ny, 0:nx]
    for (px, py), r, w in zip(pts, radii, weights):
        cx, cy, rb = (px - t0[0]) / b, (py - t0[1]) / b, r / b
        x_lo, x_hi = int(max(cx - rb - 1, 0)), int(min(cx + rb + 2, nx))
        y_lo, y_hi = int(max(cy - rb - 1, 0)), int(min(cy + rb + 2, ny))
        sub = (xx[y_lo:y_hi, x_lo:x_hi] - cx) ** 2 + (yy[y_lo:y_hi, x_lo:x_hi] - cy) ** 2 <= rb ** 2
        tpl[y_lo:y_hi, x_lo:x_hi][sub] = w
    return tpl, t0


def coarse_search(P, diam, weights, img, img0, b, max_tilt, max_scale):
    c = P.mean(0)
    Pc = P - c
    tilts = np.arange(-max_tilt, max_tilt + 1e-9, 2.0) if max_tilt > 0 else np.array([0.0])
    scales = np.arange(1 - max_scale, 1 + max_scale + 1e-9, 0.02) if max_scale > 0 else np.array([1.0])
    best = None
    results = []
    for name, M in orientations():
        best_o = None
        for tilt in tilts:
            for s in scales:
                A = s * rot(tilt) @ M
                pts = Pc @ A.T
                radii = 0.45 * diam * s
                tpl, t0 = disc_template(pts, radii, weights, b)
                cc = signal.fftconvolve(img, tpl[::-1, ::-1], mode="full")
                iy, ix = np.unravel_index(int(np.argmax(cc)), cc.shape)
                score = cc[iy, ix] / tpl.sum()
                # template pixel (0,0) sits on image pixel (ix - nx + 1, iy - ny + 1)
                off = np.array([ix - tpl.shape[1] + 1, iy - tpl.shape[0] + 1])
                t = np.array(img0) + off * b - t0 - A @ c        # Xenium = A @ P + t
                cand = dict(orientation=name, tilt=float(tilt), scale=float(s), score=float(score), A=A, t=t)
                if best_o is None or score > best_o["score"]:
                    best_o = cand
        results.append(best_o)
        if best is None or best_o["score"] > best["score"]:
            best = best_o
    results.sort(key=lambda r: -r["score"])
    print("orientation search (score = tissue inside the core discs; higher is better):")
    for r in results:
        print(f"  {r['orientation']:<14} tilt {r['tilt']:+.0f} deg  scale {r['scale']:.2f}  score {r['score']:.3f}")
    if len(results) > 1 and results[1]["score"] > 0.9 * results[0]["score"]:
        print("WARNING: the two best orientations score almost the same -- check the plot carefully")
    return best, results


# ---------------------------------------------------------------------------
# refinement on the cells
# ---------------------------------------------------------------------------
def measure_centres(tree, xy, pred, radius, iters=4):
    """Mean shift of every predicted centre onto the Xenium cells around it."""
    c = pred.copy()
    n = np.zeros(len(pred), dtype=int)
    for _ in range(iters):
        for i in range(len(c)):
            idx = tree.query_ball_point(c[i], radius[i])
            n[i] = len(idx)
            if n[i]:
                c[i] = xy[idx].mean(0)
    return c, n


def fit_affine(P, Q, w=None):
    X = np.column_stack([P, np.ones(len(P))])
    W = np.ones(len(P)) if w is None else w
    sol, *_ = np.linalg.lstsq(X * W[:, None], Q * W[:, None], rcond=None)
    return sol[:2].T, sol[2]                                 # A (2x2), t


def refine(P, diam, xy, A, t, min_cells, spacing, rounds=5):
    tree = cKDTree(xy)
    keep = None
    for _ in range(rounds):
        pred = P @ A.T + t
        s = math.sqrt(abs(np.linalg.det(A)))
        obs, n = measure_centres(tree, xy, pred, 0.5 * diam * s)
        ok = n >= min_cells
        res = np.linalg.norm(obs - pred, axis=1)
        if ok.sum() >= 3:
            lim = max(3 * np.median(res[ok]), 0.1 * spacing)
            keep = ok & (res <= lim)
            if keep.sum() >= 3:
                A, t = fit_affine(P[keep], obs[keep])
    pred = P @ A.T + t
    s = math.sqrt(abs(np.linalg.det(A)))
    obs, n = measure_centres(tree, xy, pred, 0.5 * diam * s)
    return A, t, pred, obs, n, keep


def lattice_steps(names, P):
    """The two grid steps (CellScape um) from the core names ('A_3-c' = number 3, letter c), or None."""
    m = pd.Series(names).str.extract(r"(\d+)\W*([A-Za-z]+)$")
    if m.isna().any().any():
        return None
    num = m[0].astype(int).to_numpy()
    let = m[1].str.lower().map(lambda x: sum((ord(ch) - 96) * 26 ** k for k, ch in enumerate(reversed(x)))).to_numpy()
    X = np.column_stack([np.ones(len(P)), num, let])
    if np.linalg.matrix_rank(X) < 3:
        return None
    sol = np.linalg.lstsq(X, P, rcond=None)[0]
    return sol[1], sol[2]                       # um per number step, um per letter step


def grow_grid(P, xy, tree_xy, pred0, spacing_P, sp, min_cells, max_dev):
    """Place the cores outward from a well-fitting seed near the middle: each next core is
    predicted from its already placed neighbours (local affine, so the grid bends with the
    section) and moved onto its own tissue (cells closer to a placed core are not used).
    Returns positions and which cores were found on Xenium tissue."""
    n = len(P)
    near1 = cKDTree(P).query_ball_point(P, 1.5 * spacing_P)
    near3 = cKDTree(P).query_ball_point(P, 3.2 * spacing_P)
    pos = pred0.copy()
    found = np.zeros(n, dtype=bool)
    done = np.zeros(n, dtype=bool)
    r = 0.4 * sp

    def measure(c, mtree):
        k = 0
        for _ in range(6):
            idx = tree_xy.query_ball_point(c, r)
            if not idx:
                return c, 0
            pts = xy[idx]
            if mtree is not None:
                dm, _ = mtree.query(pts)
                pts = pts[dm > np.linalg.norm(pts - c, axis=1)]
            k = len(pts)
            if k == 0:
                return c, 0
            c_new = pts.mean(0)
            if np.linalg.norm(c_new - c) < 1.0:
                return c_new, k
            c = c_new
        return c, k

    def predict(i):
        q = [j for j in near3[i] if found[j] and j != i]
        if len(q) >= 3:
            Pq = P[q] - P[q].mean(0)
            sv = np.linalg.svd(Pq, compute_uv=False)
            if sv[-1] > 0.25 * sv[0]:
                w = 1.0 / (1.0 + np.linalg.norm(P[q] - P[i], axis=1) / spacing_P)
                Al, tl = fit_affine(P[q], pos[q], w)
                return Al @ P[i] + tl
        if q:
            j = min(q, key=lambda j: np.linalg.norm(P[j] - P[i]))
            return pred0[i] + (pos[j] - pred0[j])
        return pred0[i]

    seed = None
    for i in np.argsort(np.linalg.norm(P - P.mean(0), axis=1))[:30]:
        c, k = measure(pred0[i], None)
        if k >= min_cells and np.linalg.norm(c - pred0[i]) <= 0.3 * sp:
            seed = i
            pos[i], found[i], done[i] = c, True, True
            break
    if seed is None:
        return pos, found
    frontier = {j for j in near1[seed] if not done[j]}
    while frontier:
        mtree = cKDTree(pos[found])
        new = sorted(frontier, key=lambda i: -sum(found[j] for j in near1[i]))
        for i in new:
            p = predict(i)
            c, k = measure(p, mtree)
            if k >= min_cells and np.linalg.norm(c - p) <= max_dev * sp:
                pos[i], found[i] = c, True
            else:
                pos[i] = p
            done[i] = True
        frontier = {j for i in new for j in near1[i] if not done[j]}
    for i in np.where(~done)[0]:                # cores not connected to the rest
        pos[i] = predict(i)
    return pos, found


def local_refine(P, xy, start, spacing, min_cells, max_dev, k_nb=12, iters=15):
    """Lloyd-style refinement constrained by the neighbours: cells -> nearest core (within
    0.6 spacing), core -> centre of its cells, unless that is more than max_dev x spacing from
    the position a local affine fit to its k_nb grid neighbours predicts."""
    n_core = len(P)
    nb = cKDTree(P).query(P, k=min(k_nb + 1, n_core))[1][:, 1:]           # neighbours in grid space
    C = start.copy()
    reach = 0.6 * spacing
    for _ in range(iters):
        d, j = cKDTree(C).query(xy, k=1, distance_upper_bound=reach)
        inside = np.isfinite(d)
        lab = j[inside]
        n = np.bincount(lab, minlength=n_core)
        sx = np.bincount(lab, weights=xy[inside, 0], minlength=n_core)
        sy = np.bincount(lab, weights=xy[inside, 1], minlength=n_core)
        good = n >= min_cells
        obs = np.column_stack([sx, sy]) / np.maximum(n, 1)[:, None]
        L = start.copy()
        for i in range(n_core):
            q = nb[i][good[nb[i]]]
            if len(q) >= 4:
                w = 1.0 / (1.0 + np.linalg.norm(P[q] - P[i], axis=1) / spacing)
                Al, tl = fit_affine(P[q], obs[q], w)
                L[i] = Al @ P[i] + tl
        ok = good & (np.linalg.norm(obs - L, axis=1) <= max_dev * spacing)
        newC = np.where(ok[:, None], obs, L)
        moved = float(np.abs(newC - C).max())
        C = newC
        if moved < 1.0:
            break
    return C, ok, L


def xenium_pixel_size(xenium_dir):
    f = os.path.join(xenium_dir, "experiment.xenium")
    try:
        with open(f) as fh:
            return float(json.load(fh)["pixel_size"])
    except Exception:
        return 0.2125


PRESENCE_COLOURS = {"both": [242, 194, 0], "xenium_only": [47, 168, 79], "cellscape_only": [208, 48, 48],
                    "neither": [154, 154, 154]}


def write_qupath_geojson(out, px_um, path, n_vertices=48):
    """Cores as QuPath annotations (circles) in pixels of the Xenium morphology image."""
    th = np.linspace(0, 2 * np.pi, n_vertices, endpoint=False)
    feats = []
    for r in out.itertuples(index=False):
        rad = r.diameter_um / 2 / px_um
        cx, cy = r.x_um / px_um, r.y_um / px_um
        ring = [[float(cx + rad * np.cos(t)), float(cy + rad * np.sin(t))] for t in th]
        ring.append(ring[0])
        feats.append({"type": "Feature",
                      "geometry": {"type": "Polygon", "coordinates": [ring]},
                      "properties": {"objectType": "annotation", "name": str(r.core_name),
                                     "classification": {"name": r.presence,
                                                        "color": PRESENCE_COLOURS.get(r.presence, [200, 200, 200])},
                                     "isLocked": False}})
    with open(path, "w") as fh:
        json.dump({"type": "FeatureCollection", "features": feats}, fh)


def describe(A):
    """Rotation (deg), mirrored, scale x / y of a 2x2 matrix."""
    det = np.linalg.det(A)
    mirrored = det < 0
    B = A @ np.diag([-1.0, 1.0]) if mirrored else A
    angle = math.degrees(math.atan2(B[1, 0], B[0, 0]))
    return dict(rotation_deg=round(angle, 2), mirrored=bool(mirrored),
                scale_x=round(float(np.linalg.norm(A[:, 0])), 4), scale_y=round(float(np.linalg.norm(A[:, 1])), 4))


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
def plot(xy, out, slide_id, sectors, info):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    colours = {"both": "#f2c200", "xenium_only": "#2fa84f", "cellscape_only": "#d03030", "neither": "#9a9a9a"}
    span = xy.max(0) - xy.min(0)
    w = 14
    fig, ax = plt.subplots(figsize=(w, max(4, w * span[1] / max(span[0], 1))))
    b = 20.0
    h, xe, ye = np.histogram2d(xy[:, 0], xy[:, 1], bins=[max(1, int(span[0] / b)), max(1, int(span[1] / b))])
    ax.imshow(np.log1p(h.T), extent=[xe[0], xe[-1], ye[-1], ye[0]], cmap="gray", interpolation="nearest")
    for r in out.itertuples():
        col = colours[r.presence]
        ax.add_patch(Circle((r.x_um, r.y_um), r.diameter_um / 2, fill=False, lw=0.8, ec=col))
        ax.text(r.x_um, r.y_um, r.core_name, color=col, fontsize=5, ha="center", va="center")
    ax.set_xlim(min(xe[0], out.x_um.min() - 600), max(xe[-1], out.x_um.max() + 600))
    ax.set_ylim(max(ye[-1], out.y_um.max() + 600), min(ye[0], out.y_um.min() - 600))     # y down, as Xenium
    ax.set_aspect("equal")
    ax.set_xlabel("x (um, Xenium)")
    ax.set_ylabel("y (um, Xenium)")
    counts = out["presence"].value_counts()
    ax.set_title(f"{slide_id}: CellScape sector {'+'.join(sectors)} on the Xenium cells  "
                 f"[{info['orientation']}, rotation {info['rotation_deg']:.1f} deg, scale {info['scale_x']:.3f}]\n"
                 + "   ".join(f"{k}: {counts.get(k, 0)}" for k in colours)
                 + "   (yellow = tissue on both, green = Xenium only, red = CellScape only, grey = neither)",
                 fontsize=9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def main(argv=None):
    a = parse_args(argv)
    out_dir = a.out_dir or os.path.join(HERE, "config", "core_grid")
    os.makedirs(out_dir, exist_ok=True)

    if isinstance(a.sector, str):          # notebook: sector='C' or 'A, B' instead of a list
        a.sector = [t for t in re.split(r"[\s,;]+", a.sector.strip()) if t]
    g = read_cellscape_cores(a.cores_csv, a.sector)
    xy = read_xenium_cells(a.xenium_dir)
    P = g[["centroid_x_um", "centroid_y_um"]].to_numpy(dtype=np.float64)
    cs_diam = g["diameter_um"].to_numpy(dtype=np.float64)
    # one size for all cores: the CellScape ellipses vary with how each one was drawn
    diam = np.full(len(P), float(np.median(cs_diam)))
    spacing = float(np.median(cKDTree(P).query(P, k=2)[0][:, 1]))
    weights = np.where(g["cellscape_missing"], 0.5, 1.0)

    img, img0 = density_image(xy, a.bin_um, pad_um=float(np.ptp(P, axis=0).max()) * 0.6)
    best, _ = coarse_search(P, diam, weights, img, img0, a.bin_um, a.max_tilt, a.max_scale_change)
    A, t, pred, obs, n, keep = refine(P, diam, xy, best["A"], best["t"], a.min_cells, spacing)

    s = math.sqrt(abs(np.linalg.det(A)))
    d_um = a.diameter_um or float(diam[0] * s)
    tree = cKDTree(xy)
    sp = spacing * s

    def coverage(C):
        return float(np.mean(cKDTree(C).query(xy, k=1)[0] <= 0.5 * d_um))

    if a.max_local_shift > 0:
        # the grid is periodic: also try it one row / column off and keep the placement whose
        # 10 x 20 (or whatever) grid covers most of the Xenium tissue -- the TMA edges decide
        steps = lattice_steps(g["core"].astype(str), P)
        shifts = [(0, 0)] + ([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)] if steps else [])
        cands = []
        for i, j in shifts:
            delta = (i * steps[0] + j * steps[1]) if steps else np.zeros(2)
            pos, found = grow_grid(P, xy, tree, pred + A @ delta, spacing, sp, a.min_cells, a.max_local_shift)
            cands.append(dict(shift=(i, j), pos=pos, found=found, cov=coverage(pos)))
        cands.sort(key=lambda c: -c["cov"])
        print("grid position (number step, letter step) -> Xenium cells covered by the cores:")
        for c in cands:
            print(f"  shift {str(c['shift']):<9} {100 * c['cov']:.1f}%  ({int(c['found'].sum())} cores on tissue)")
        if len(cands) > 1 and cands[1]["cov"] > cands[0]["cov"] - 0.01:
            print("WARNING: two grid positions cover the tissue almost equally -- check the core labels on the plot")
        best_c = cands[0]
        final, local, _ = local_refine(P, xy, best_c["pos"], sp, a.min_cells, 0.25)
        A, t = fit_affine(P[local], final[local]) if local.sum() >= 3 else (A, t)
        pred = P @ A.T + t
    else:
        final, local = pred.copy(), np.zeros(len(P), dtype=bool)
    dist = np.linalg.norm(final - pred, axis=1)
    n_final = np.array([len(tree.query_ball_point(c, 0.5 * d_um)) for c in final])
    xen_missing = n_final < a.min_cells
    cs_missing = g["cellscape_missing"].to_numpy()
    presence = np.select([~xen_missing & ~cs_missing, ~xen_missing & cs_missing, xen_missing & ~cs_missing],
                         ["both", "xenium_only", "cellscape_only"], "neither")

    out = pd.DataFrame({
        "core_name": g["core"].astype(str).str.replace(r"^[A-Za-z]+_", "", regex=True),
        "x_um": final[:, 0].round(2), "y_um": final[:, 1].round(2),
        "diameter_um": round(d_um, 1),
        "missing": xen_missing,
        "cellscape_core": g["core"].astype(str), "sector": g["sector"],
        "cellscape_missing": cs_missing,
        "n_cells_xenium": n_final,
        "shift_um": dist.round(1),
        "cellscape_diameter_um": cs_diam.round(1),
        "presence": presence,
    })

    # how much of the Xenium tissue the grid covers (cells within 0.6 x diameter of a core centre)
    near = cKDTree(final).query(xy, k=1)
    covered = float(np.mean(near[0] <= 0.6 * d_um))
    fit_res = np.linalg.norm(obs[keep] - pred[keep], axis=1) if keep is not None and keep.any() else np.array([np.nan])
    info = dict(orientation=best["orientation"], **describe(A))
    info.update(
        slide_id=a.slide_id, sectors=a.sector, cellscape_cores_csv=os.path.abspath(a.cores_csv),
        xenium_dir=os.path.abspath(a.xenium_dir),
        note="x_xenium_um = A @ [x_cellscape_um, y_cellscape_um] + t ; CellScape um = CellScape file pixels x 0.18 "
             "(file orientation, not the mirrored QuPath view)",
        A=A.round(8).tolist(), t=t.round(4).tolist(),
        cores_used_for_fit=int(keep.sum()) if keep is not None else 0,
        fit_residual_um_median=round(float(np.median(fit_res)), 1),
        fit_residual_um_p90=round(float(np.percentile(fit_res, 90)), 1),
        core_spacing_um=round(spacing * s, 1),
        xenium_cells_inside_grid=round(covered, 3),
        search_score=round(best["score"], 4),
    )

    grid_csv = os.path.join(out_dir, f"{a.slide_id}.csv")
    out.to_csv(grid_csv, index=False)
    with open(os.path.join(out_dir, f"{a.slide_id}_cellscape_to_xenium.json"), "w") as fh:
        json.dump(info, fh, indent=2)
    px_um = xenium_pixel_size(a.xenium_dir)
    gj = os.path.join(out_dir, f"{a.slide_id}_grid_qupath.geojson")
    write_qupath_geojson(out, px_um, gj)
    fig = plot(xy, out, a.slide_id, a.sector, info)
    png = os.path.join(out_dir, f"{a.slide_id}_grid_on_xenium.png")
    fig.savefig(png, dpi=200)

    counts = out["presence"].value_counts()
    print(f"\nCellScape -> Xenium: {info['orientation']} (rotation {info['rotation_deg']} deg, "
          f"mirrored {info['mirrored']}), scale {info['scale_x']:.3f} / {info['scale_y']:.3f}")
    print(f"affine fit on {info['cores_used_for_fit']} cores: residual median {info['fit_residual_um_median']} um, "
          f"90% {info['fit_residual_um_p90']} um (core spacing {info['core_spacing_um']} um)")
    print(f"local refinement: {int(local.sum())} cores placed on their own Xenium tissue; distance from the affine "
          f"grid median {np.median(dist):.0f} um, max {dist.max():.0f} um; core diameter {d_um:.0f} um")
    print("cores: " + ", ".join(f"{k} {counts.get(k, 0)}" for k in ("both", "xenium_only", "cellscape_only", "neither")))
    print(f"{100 * covered:.0f}% of the Xenium cells lie on a grid core")
    if covered < 0.7:
        print("WARNING: much Xenium tissue is outside the grid -- wrong sector, or the fit is off. Check the plot.")
    if info["fit_residual_um_median"] > 0.15 * info["core_spacing_um"]:
        print("WARNING: large fit residuals -- check the plot before using the grid.")
    edited = os.path.join(out_dir, f"{a.slide_id}_edited.csv")
    if os.path.isfile(edited):
        print(f"NOTE: {os.path.basename(edited)} (your manual corrections) exists -- 02 uses it, not this new grid. "
              "Delete it to use the new grid.")
    print(f"\nwritten: {grid_csv}\n         {png}\n         {gj}  (QuPath, for manual corrections)\n"
          f"Next: 02_tma_dearray.R (QuPath grid mode picks up {os.path.basename(grid_csv)}); its core names "
          f"('{out['core_name'].iloc[0]}') are matched to tma_map.csv by row + column.")


if __name__ == "__main__":
    main()
