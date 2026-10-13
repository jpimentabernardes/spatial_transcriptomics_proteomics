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
     (mean shift), outliers dropped, repeated until stable. Each core may then move by
     up to --max-local-shift x core spacing to its own Xenium centre (serial sections
     are not perfectly identical)
  4. per core: number of Xenium cells -> missing on Xenium, next to the CellScape
     "missing" flag (both / xenium_only / cellscape_only / neither)

Outputs (default folder: <scripts>/config/core_grid/, read by 02_tma_dearray.R)
  <slide_id>.csv               grid for 02 (QuPath mode): core_name, x_um, y_um, diameter_um,
                               missing (= no tissue on Xenium) + cellscape_core, sector,
                               cellscape_missing, n_cells_xenium, shift_um, presence
  <slide_id>_cellscape_to_xenium.json   the transform (CellScape um / px -> Xenium um), also
                               usable to bring CellScape cells into the Xenium frame
  <slide_id>_grid_on_xenium.png         check plot: the grid over the Xenium cells

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
    p.add_argument("--max-local-shift", type=float, default=0.2,
                   help="a core may move to its own Xenium centre by up to this x core spacing (0 = rigid grid)")
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
    diam = g["diameter_um"].to_numpy(dtype=np.float64)
    spacing = float(np.median(cKDTree(P).query(P, k=2)[0][:, 1]))
    weights = np.where(g["cellscape_missing"], 0.5, 1.0)

    img, img0 = density_image(xy, a.bin_um, pad_um=float(np.ptp(P, axis=0).max()) * 0.6)
    best, _ = coarse_search(P, diam, weights, img, img0, a.bin_um, a.max_tilt, a.max_scale_change)
    A, t, pred, obs, n, keep = refine(P, diam, xy, best["A"], best["t"], a.min_cells, spacing)

    shift = obs - pred
    dist = np.linalg.norm(shift, axis=1)
    local = (n >= a.min_cells) & (dist <= a.max_local_shift * spacing)
    final = np.where(local[:, None], obs, pred)
    s = math.sqrt(abs(np.linalg.det(A)))
    tree = cKDTree(xy)
    n_final = np.array([len(tree.query_ball_point(c, 0.5 * d * s)) for c, d in zip(final, diam)])
    xen_missing = n_final < a.min_cells
    cs_missing = g["cellscape_missing"].to_numpy()
    presence = np.select([~xen_missing & ~cs_missing, ~xen_missing & cs_missing, xen_missing & ~cs_missing],
                         ["both", "xenium_only", "cellscape_only"], "neither")

    out = pd.DataFrame({
        "core_name": g["core"].astype(str).str.replace(r"^[A-Za-z]+_", "", regex=True),
        "x_um": final[:, 0].round(2), "y_um": final[:, 1].round(2),
        "diameter_um": (diam * s).round(1),
        "missing": xen_missing,
        "cellscape_core": g["core"].astype(str), "sector": g["sector"],
        "cellscape_missing": cs_missing,
        "n_cells_xenium": n_final,
        "shift_um": np.where(local, dist, 0).round(1),
        "presence": presence,
    })

    # how much of the Xenium tissue the grid covers (cells within 0.6 x diameter of a core centre)
    near = cKDTree(final).query(xy, k=1)
    covered = float(np.mean(near[0] <= 0.6 * diam[near[1]] * s))
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
    fig = plot(xy, out, a.slide_id, a.sector, info)
    png = os.path.join(out_dir, f"{a.slide_id}_grid_on_xenium.png")
    fig.savefig(png, dpi=200)

    counts = out["presence"].value_counts()
    print(f"\nCellScape -> Xenium: {info['orientation']} (rotation {info['rotation_deg']} deg, "
          f"mirrored {info['mirrored']}), scale {info['scale_x']:.3f} / {info['scale_y']:.3f}")
    print(f"affine fit on {info['cores_used_for_fit']} cores: residual median {info['fit_residual_um_median']} um, "
          f"90% {info['fit_residual_um_p90']} um (core spacing {info['core_spacing_um']} um)")
    print(f"{int(local.sum())} cores moved to their own Xenium centre (median {np.median(dist[local]) if local.any() else 0:.0f} um)")
    print("cores: " + ", ".join(f"{k} {counts.get(k, 0)}" for k in ("both", "xenium_only", "cellscape_only", "neither")))
    print(f"{100 * covered:.0f}% of the Xenium cells lie on a grid core")
    if covered < 0.7:
        print("WARNING: much Xenium tissue is outside the grid -- wrong sector, or the fit is off. Check the plot.")
    if info["fit_residual_um_median"] > 0.15 * info["core_spacing_um"]:
        print("WARNING: large fit residuals -- check the plot before using the grid.")
    print(f"\nwritten: {grid_csv}\n         {png}\n"
          f"Next: 02_tma_dearray.R (QuPath grid mode picks up {os.path.basename(grid_csv)}); its core names "
          f"('{out['core_name'].iloc[0]}') are matched to tma_map.csv by row + column.")


if __name__ == "__main__":
    main()
