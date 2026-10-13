#!/usr/bin/env python
"""02c_check_cellscape_xenium_match.py -- do the CellScape cores really match the Xenium cores?

Independent check of 02b. The CellScape and Xenium sections are serial sections of the
same block, so the tissue inside every core (its outline, tears, holes, fragments) looks
the same on both -- but only with the right orientation and the right core label. For every
core found on both sections this script compares the CellScape DNA image of the core with
the Xenium cells of the core:

  1. orientation test: the CellScape core in all 8 orientations (0/90/180/270 degrees, with
     and without mirroring; each with a small extra tilt and shift search) against the
     Xenium core -> which orientation fits best, core by core. If 02b is right, nearly all
     cores pick the orientation 02b chose.
  2. label test: the Xenium core against its own CellScape core and against the 8 CellScape
     cores around it (same orientation) -> the own core should fit best. A grid that is one
     row or column off fails this test.

CellScape image: either the <slide>_TMA_preview.png pathology sent (the grid drawn over the
DNA channel; the coloured circles and labels are ignored), or, better, a DNA image of the
CellScape scan as OME-TIFF (only one channel of one low-resolution pyramid level is read; the
DAPI / DNA channel is found by name, or give --cellscape-channel). NOT the *_TMA_labels.ome.tif
(the core label mask has no tissue in it; 02c refuses it). For the PNG the full CellScape image size is needed to convert pixels to
um: known for FL0026 / FL0027 (pathology README), otherwise --cellscape-size-px WxH.

Outputs (next to the 02b grid, default <scripts>/config/core_grid/)
  <slide_id>_match_summary.png   orientation votes, own core vs neighbours, one core in all 8
                                 orientations
  <slide_id>_match_gallery.png   CellScape (magenta) over Xenium (green), white = agreement
  <slide_id>_match_per_core.csv  per core: score per orientation, best orientation, own vs
                                 best neighbour score

Example
  python 02c_check_cellscape_xenium_match.py --slide-id TMA_ORGAN --sector A \\
      --cores-csv FL0027_TMA_cores.csv --cellscape-image FL0027_TMA_preview.png \\
      --xenium-dir /path/to/output-XETG...
"""

import argparse
import json
import math
import os
import re

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
CS_UM_PER_PX = 0.18                                     # CellScape full-resolution pixel size
KNOWN_SIZES = {"FL0026": (174689, 91289), "FL0027": (179811, 91151)}    # pathology README


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True, help="Xenium slide_id, as for 02b")
    p.add_argument("--cores-csv", required=True, help="<slide>_TMA_cores.csv of the CellScape dearray")
    p.add_argument("--sector", nargs="+", required=True, help="sector letter(s), as for 02b, e.g. ['A']")
    p.add_argument("--cellscape-image", required=True,
                   help="<slide>_TMA_preview.png, or a CellScape DNA image (OME-TIFF)")
    p.add_argument("--xenium-dir", required=True, help="Xenium output folder (output-XETG...) of this slide")
    p.add_argument("--grid-dir", default=None,
                   help="folder with the 02b outputs (<slide_id>.csv + _cellscape_to_xenium.json); "
                        "default <scripts>/config/core_grid")
    p.add_argument("--cellscape-size-px", default=None,
                   help="full CellScape image size 'WxH' (only for a PNG preview; known for FL0026/FL0027)")
    p.add_argument("--cellscape-channel", default=None,
                   help="OME-TIFF channel with DNA, as name (e.g. DAPI) or number; default: the first channel "
                        "named DAPI / DNA / Hoechst")
    p.add_argument("--resolution-um", type=float, default=10.0, help="pixel size of the comparison")
    p.add_argument("--max-shift-um", type=float, default=200.0, help="shift search per core")
    p.add_argument("--n-gallery", type=int, default=48, help="cores shown in the gallery")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# images
# ---------------------------------------------------------------------------
def fill_masked(img, valid, sigma):
    """Replace invalid pixels (drawn circles / text) by the blurred valid neighbourhood."""
    num = ndimage.gaussian_filter(np.where(valid, img, 0.0), sigma)
    den = ndimage.gaussian_filter(valid.astype(float), sigma)
    out = img.copy()
    out[~valid] = (num / np.maximum(den, 1e-6))[~valid]
    return out


def ring_score(overlay, cores, um_x, um_y, ox, oy):
    """Fraction of points on the expected core circles that fall on drawn (coloured) pixels."""
    th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    r = cores["diameter_um"].to_numpy()[:, None] / 2
    x = (cores["centroid_x_um"].to_numpy()[:, None] + r * np.cos(th)) / um_x + ox
    y = (cores["centroid_y_um"].to_numpy()[:, None] + r * np.sin(th)) / um_y + oy
    return float(ndimage.map_coordinates(overlay, [y.ravel(), x.ravel()], order=1, mode="constant").mean())


def load_preview(path, size_px, all_cores):
    from PIL import Image
    rgb = np.asarray(Image.open(path).convert("RGB")).astype(np.float64)
    nonwhite = rgb.min(-1) < 235
    rows = np.where(nonwhite.mean(1) > 0.9)[0]
    cols = np.where(nonwhite.mean(0) > 0.9)[0]
    if len(rows) < 10 or len(cols) < 10:
        raise ValueError(f"{path}: could not find the image area inside the figure")
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    rgb = rgb[r0:r1, c0:c1]
    W, H = size_px
    um_x = W * CS_UM_PER_PX / rgb.shape[1]
    um_y = H * CS_UM_PER_PX / rgb.shape[0]
    if abs(um_x / um_y - 1) > 0.03:
        print(f"WARNING: preview aspect does not match the CellScape size {W}x{H} ({um_x:.2f} vs {um_y:.2f} um/px)")
    sat = rgb.max(-1) - rgb.min(-1)
    overlay = sat > 40                                    # yellow / red circles and labels
    gray = rgb.mean(-1)
    # check (and fine-tune) the pixel <-> um calibration on the drawn circles
    ov = ndimage.gaussian_filter(overlay.astype(float), 1.0)
    best = (ring_score(ov, all_cores, um_x, um_y, 0, 0), 1.0, 0.0, 0.0)
    for f in np.linspace(0.985, 1.015, 13):
        for ox in range(-3, 4):
            for oy in range(-3, 4):
                sc = ring_score(ov, all_cores, um_x * f, um_y * f, ox, oy)
                if sc > best[0]:
                    best = (sc, f, ox, oy)
    sc, f, ox, oy = best
    base = float(overlay.mean())
    print(f"CellScape preview: {rgb.shape[1]}x{rgb.shape[0]} px, {um_x * f:.1f} um/px; the drawn core circles "
          f"are found at the expected places ({100 * sc:.0f}% of circle points on drawn pixels, background "
          f"{100 * base:.0f}%)")
    if sc < 3 * base:
        print("WARNING: the drawn circles do not line up with the core table -- is this the preview of the same "
              "slide, and in file orientation?")
    return gray, ~overlay, (um_x * f, um_y * f), (-ox * um_x * f, -oy * um_y * f)


def ome_info(tf):
    """Channel names and pixel size (um) from the OME-XML, if present."""
    xml = tf.ome_metadata or ""
    names = re.findall(r'<Channel[^>]*?Name="([^"]*)"', xml)
    m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', xml)
    unit = re.search(r'PhysicalSizeXUnit="([^"]+)"', xml)
    px = float(m.group(1)) if m else None
    if px is not None and unit and unit.group(1) in ("nm",):
        px /= 1000.0
    return names, px


def load_ome(path, channel, target_um):
    """One channel of a CellScape OME-TIFF, from the pyramid level closest to target_um."""
    import tifffile
    if "label" in os.path.basename(path).lower():
        raise ValueError(f"{os.path.basename(path)} is the core LABEL MASK (0 = background, k = core), not an image "
                         "of the tissue -- it holds only the drawn core ellipses, so it cannot test anything. Use "
                         "the CellScape scan itself (OME-TIFF with a DAPI/DNA channel), or the "
                         "<slide>_TMA_preview.png.")
    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        axes = s.axes
        levels = list(getattr(s, "levels", None) or [s])
        names, px0 = ome_info(tf)
        px0 = px0 or CS_UM_PER_PX
        full_w = s.shape[axes.index("X")]
        k = 0
        for i, lv in enumerate(levels):
            if full_w / lv.shape[axes.index("X")] * px0 <= target_um:
                k = i
        lv = levels[k]
        n_c = lv.shape[axes.index("C")] if "C" in axes else 1
        # which channel: given as index or name, else the first DNA-like name
        if channel in (None, ""):
            hit = [i for i, n in enumerate(names[:n_c]) if re.search(r"dapi|dna|hoechst|nucle", n, re.I)]
            ci = hit[0] if hit else 0
        elif str(channel).isdigit():
            ci = int(channel)
        else:
            hit = [i for i, n in enumerate(names) if n.lower() == str(channel).lower()]
            if not hit:
                raise ValueError(f"channel {channel!r} not in this OME-TIFF; channels: {names}")
            ci = hit[0]
        if ci >= n_c:
            raise ValueError(f"channel {ci} does not exist (the image has {n_c} channels: {names})")
        cname = names[ci] if ci < len(names) else f"#{ci}"
        if "C" in axes and n_c > 1 and len(lv.pages) >= n_c:
            arr = tifffile.imread(path, series=0, level=k, key=ci)       # just this channel's plane
        else:
            arr = lv.asarray()
            if "C" in axes and n_c > 1:
                arr = np.take(arr, ci, axis=axes.index("C"))
        w = lv.shape[axes.index("X")]
    arr = np.squeeze(np.asarray(arr))
    if arr.ndim == 3:                                   # RGB(A) or extra planes: average them
        arr = arr.mean(axis=-1 if arr.shape[-1] <= 4 else 0)
    if np.issubdtype(arr.dtype, np.integer) and n_c == 1:
        u = np.unique(arr)
        if 0 in u and 10 < len(u) < 2000 and u.max() < 5000 and len(u) > 0.5 * (u.max() + 1) \
                and (arr == 0).mean() > 0.3:
            raise ValueError(f"{os.path.basename(path)} looks like a label mask (values 0..{u.max()}), not a DNA "
                             "image of the tissue. Use the CellScape scan itself or the preview PNG.")
    um = full_w / w * px0
    print(f"CellScape image: {os.path.basename(path)}, channel {ci} ({cname}) of {n_c}, pyramid level {k}, "
          f"{arr.shape[1]}x{arr.shape[0]} px, {um:.2f} um/px")
    if n_c > 1 and not re.search(r"dapi|dna|hoechst|nucle", cname, re.I):
        print(f"WARNING: channel {cname!r} does not look like DNA; channels are: {names}. "
              "Set cellscape_channel to the DAPI / DNA channel name.")
    return arr.astype(np.float64), np.ones(arr.shape, bool), (um, um), (0.0, 0.0)


def xenium_density(xy, b):
    x0, y0 = xy.min(0) - 50 * b
    nx = int(math.ceil((xy[:, 0].max() - x0) / b)) + 50
    ny = int(math.ceil((xy[:, 1].max() - y0) / b)) + 50
    h, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[ny, nx], range=[[y0, y0 + ny * b], [x0, x0 + nx * b]])
    return ndimage.gaussian_filter(h, 1.0), (x0, y0)


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def orientations():
    out = []
    for mirror in (False, True):
        for k in range(4):
            a = math.radians(90 * k)
            r = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]]).round(12)
            out.append((f"rot{90 * k}{'+mirror' if mirror else ''}", r @ np.diag([-1.0, 1.0]) if mirror else r))
    return out


def rot(deg):
    a = math.radians(deg)
    return np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])


def sample(img, origin, um, pts):
    """Image value at points given in um (bilinear)."""
    cols = (pts[..., 0] - origin[0]) / um[0] - 0.5
    rows = (pts[..., 1] - origin[1]) / um[1] - 0.5
    return ndimage.map_coordinates(img, [rows, cols], order=1, mode="constant", cval=0.0)


def norm_patch(p):
    p = np.log1p(np.maximum(p, 0))
    sd = p.std()
    return (p - p.mean()) / sd if sd > 0 else p * 0


def best_ncc(big, small):
    """Max normalised cross-correlation of `small` over all positions inside `big`."""
    from skimage.feature import match_template
    if small.std() == 0 or big.std() == 0:
        return -1.0, (0, 0)
    m = match_template(big, small)
    i = np.unravel_index(int(np.argmax(m)), m.shape)
    return float(m[i]), i


def patches(cs, cs_origin, cs_um, xen, xen_origin, b, Pc, Xc, L, half, pad):
    """Xenium patch around Xc, and a function giving the CellScape patch (with search margin)
    around Pc for a given CellScape->Xenium linear map."""
    u = np.arange(-half, half + 1e-9, b)
    gx, gy = np.meshgrid(u, u)
    grid = np.stack([gx, gy], -1)
    xp = norm_patch(sample(xen, xen_origin, (b, b), Xc + grid))
    ub = np.arange(-half - pad, half + pad + 1e-9, b)
    bx, by = np.meshgrid(ub, ub)
    big = np.stack([bx, by], -1)

    def cs_patch(M):
        Minv = np.linalg.inv(M)
        pts = Pc + big @ Minv.T                         # Xenium-frame offsets -> CellScape um
        return norm_patch(sample(cs, cs_origin, cs_um, pts))
    return xp, cs_patch


def overlay_tile(xp, cp):
    def unit(p):
        lo, hi = np.percentile(p, [5, 99.5])
        return np.clip((p - lo) / max(hi - lo, 1e-9), 0, 1)
    c, x = unit(cp), unit(xp)
    return np.dstack([c, x, c])                         # CellScape magenta, Xenium green, both white


# ---------------------------------------------------------------------------
def main(argv=None):
    a = parse_args(argv)
    if isinstance(a.sector, str):
        a.sector = [t for t in re.split(r"[\s,;]+", a.sector.strip()) if t]
    grid_dir = a.grid_dir or os.path.join(HERE, "config", "core_grid")
    grid = pd.read_csv(os.path.join(grid_dir, f"{a.slide_id}.csv"))
    with open(os.path.join(grid_dir, f"{a.slide_id}_cellscape_to_xenium.json")) as fh:
        info = json.load(fh)
    A = np.array(info["A"])
    M_by_name = dict(orientations())
    chosen = info["orientation"]
    L = A @ np.linalg.inv(M_by_name[chosen])            # scale + tilt part of the 02b transform

    all_cores = pd.read_csv(a.cores_csv)
    cs_cores = all_cores[all_cores["sector"].astype(str).isin(a.sector)].set_index("core")

    # images
    if a.cellscape_image.lower().endswith((".png", ".jpg", ".jpeg")):
        size = None
        if a.cellscape_size_px:
            size = tuple(int(v) for v in re.split(r"[x, ]+", a.cellscape_size_px.strip()))
        else:
            for key, v in KNOWN_SIZES.items():
                if key in os.path.basename(a.cellscape_image) or key in os.path.basename(a.cores_csv):
                    size = v
        if size is None:
            raise ValueError("give --cellscape-size-px WxH (full CellScape image size in pixels)")
        cs, valid, cs_um, cs_origin = load_preview(a.cellscape_image, size, all_cores)
    else:
        cs, valid, cs_um, cs_origin = load_ome(a.cellscape_image, a.cellscape_channel, a.resolution_um)
    cs = fill_masked(cs, valid, 2.0)
    cs = ndimage.gaussian_filter(cs, max(0.5, a.resolution_um / max(cs_um) / 2))

    xyf = os.path.join(a.xenium_dir, "cells.parquet")
    xy = pd.read_parquet(xyf, columns=["x_centroid", "y_centroid"]).to_numpy(np.float64) if os.path.isfile(xyf) \
        else pd.read_csv(os.path.join(a.xenium_dir, "cells.csv.gz"), usecols=["x_centroid", "y_centroid"]).to_numpy()
    b = max(a.resolution_um, min(cs_um))
    xen, xen_origin = xenium_density(xy, b)

    both = grid[(grid["presence"] == "both") & grid["cellscape_core"].isin(cs_cores.index)].reset_index(drop=True)
    P_all = cs_cores[["centroid_x_um", "centroid_y_um"]].to_numpy()
    nb_tree = cKDTree(P_all)
    spacing = float(np.median(nb_tree.query(P_all, k=2)[0][:, 1]))
    half = 0.55 * spacing * math.sqrt(abs(np.linalg.det(L)))
    pad = a.max_shift_um
    tilts = (-4.0, 0.0, 4.0)
    names = [n for n, _ in orientations()]
    print(f"{len(both)} cores with tissue on both sections are compared "
          f"(patch {2 * half:.0f} um, {b:.0f} um pixels, shift search +-{pad:.0f} um)")

    rows, tiles = [], []
    for r in both.itertuples(index=False):
        Pc = cs_cores.loc[r.cellscape_core, ["centroid_x_um", "centroid_y_um"]].to_numpy(dtype=float)
        Xc = np.array([r.x_um, r.y_um])
        xp, cs_patch = patches(cs, cs_origin, cs_um, xen, xen_origin, b, Pc, Xc, L, half, pad)
        score = {}
        for name, M in orientations():
            score[name] = max(best_ncc(cs_patch(rot(t) @ L @ M), xp)[0] for t in tilts)
        # label test: own core vs the 8 CellScape cores around it, in the chosen orientation
        M = M_by_name[chosen]
        _, nbr = nb_tree.query(Pc, k=9)
        own = score[chosen]
        others = []
        for j in nbr:
            name_j = cs_cores.index[j]
            if name_j == r.cellscape_core:
                continue
            _, cs_patch_j = patches(cs, cs_origin, cs_um, xen, xen_origin, b, P_all[j], Xc, L, half, pad)
            others.append(max(best_ncc(cs_patch_j(rot(t) @ L @ M), xp)[0] for t in tilts))
        best_name = max(score, key=score.get)
        rows.append(dict(core_name=r.core_name, cellscape_core=r.cellscape_core,
                         best_orientation=best_name, **{f"score_{k}": round(v, 3) for k, v in score.items()},
                         own_score=round(own, 3), best_neighbour_score=round(max(others), 3)))
        # tile for the gallery (chosen orientation, best tilt / shift)
        big = cs_patch(L @ M)
        _, (i, j) = best_ncc(big, xp)
        tiles.append((r.core_name, own, overlay_tile(xp, big[i:i + xp.shape[0], j:j + xp.shape[1]])))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(grid_dir, f"{a.slide_id}_match_per_core.csv"), index=False)

    votes = res["best_orientation"].value_counts().reindex(names, fill_value=0)
    n = len(res)
    own_wins = int((res["own_score"] > res["best_neighbour_score"]).sum())
    print("\norientation test -- cores whose tissue shape fits best in each orientation:")
    for k in names:
        print(f"  {k:<14} {votes[k]:>4}  {'<- 02b chose this' if k == chosen else ''}")
    print(f"label test: the own CellScape core fits better than all 8 neighbours for {own_wins} of {n} cores "
          f"({100 * own_wins / max(n, 1):.0f}%); median score own {res['own_score'].median():.2f}, "
          f"best neighbour {res['best_neighbour_score'].median():.2f}")
    frac = votes[chosen] / max(n, 1)
    if frac >= 0.6 and own_wins / max(n, 1) >= 0.6:
        print(f"VERDICT: orientation and labels confirmed ({100 * frac:.0f}% of cores vote for {chosen}).")
    else:
        print("VERDICT: NOT confirmed -- look at the figures; send them to me before using the grid.")

    plot_summary(a.slide_id, res, votes, chosen, names, tiles, cs, cs_origin, cs_um, xen, xen_origin, b, both,
                 cs_cores, L, half, grid_dir)
    plot_gallery(a.slide_id, tiles, chosen, a.n_gallery, grid_dir)
    print(f"\nwritten to {grid_dir}: {a.slide_id}_match_summary.png, {a.slide_id}_match_gallery.png, "
          f"{a.slide_id}_match_per_core.csv")


def plot_summary(slide_id, res, votes, chosen, names, tiles, cs, cs_origin, cs_um, xen, xen_origin, b, both,
                 cs_cores, L, half, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 8, height_ratios=[1.1, 1])
    ax = fig.add_subplot(gs[0, :4])
    cols = ["#2a78d6" if k == chosen else "#b8b8b8" for k in names]
    ax.bar(range(8), votes.to_numpy(), color=cols)
    ax.set_xticks(range(8))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_ylabel("cores fitting best")
    ax.set_title(f"Orientation test: which orientation makes each core's tissue match? (n = {len(res)})\n"
                 f"blue = the orientation 02b chose ({chosen})", fontsize=10)
    for i, v in enumerate(votes.to_numpy()):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax2 = fig.add_subplot(gs[0, 4:])
    ax2.scatter(res["best_neighbour_score"], res["own_score"], s=12, c="#2a78d6", alpha=0.7)
    lim = [min(res[["own_score", "best_neighbour_score"]].min().min(), 0) - 0.05, 1.0]
    ax2.plot(lim, lim, color="#888", lw=1, ls="--")
    ax2.set_xlim(lim)
    ax2.set_ylim(lim)
    ax2.set_xlabel("best of the 8 neighbouring CellScape cores")
    ax2.set_ylabel("own CellScape core (same name)")
    wins = int((res["own_score"] > res["best_neighbour_score"]).sum())
    ax2.set_title(f"Label test: shape match with own core vs neighbours\nabove the line = own core fits best "
                  f"({wins} of {len(res)})", fontsize=10)
    # one example core in all 8 orientations
    ex = res.loc[(res[[f"score_{k}" for k in names]].max(1) - res[[f"score_{k}" for k in names]].median(1)).idxmax()]
    r = both[both["core_name"] == ex["core_name"]].iloc[0]
    Pc = cs_cores.loc[r["cellscape_core"], ["centroid_x_um", "centroid_y_um"]].to_numpy(dtype=float)
    Xc = np.array([r["x_um"], r["y_um"]])
    xp, cs_patch = patches(cs, cs_origin, cs_um, xen, xen_origin, b, Pc, Xc, L, half, 200.0)
    for i, (name, M) in enumerate(orientations()):
        big = cs_patch(L @ M)
        sc, (p, q) = best_ncc(big, xp)
        axi = fig.add_subplot(gs[1, i])
        axi.imshow(overlay_tile(xp, big[p:p + xp.shape[0], q:q + xp.shape[1]]))
        axi.set_title(f"{name}\nscore {sc:.2f}", fontsize=8, color="#2a78d6" if name == chosen else "black")
        axi.axis("off")
    fig.text(0.5, 0.02, f"core {ex['core_name']} ({ex['cellscape_core']}): CellScape DNA (magenta) over Xenium "
             f"cells (green) in each orientation -- white = tissue on both", ha="center", fontsize=9)
    fig.suptitle(f"{slide_id}: does the CellScape grid match the Xenium cores?", fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(os.path.join(out_dir, f"{slide_id}_match_summary.png"), dpi=150)


def plot_gallery(slide_id, tiles, chosen, n_show, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not tiles:
        return
    idx = np.linspace(0, len(tiles) - 1, min(n_show, len(tiles))).round().astype(int)
    sel = [tiles[i] for i in idx]
    ncol = 8
    nrow = int(math.ceil(len(sel) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2 * ncol, 2.15 * nrow))
    for ax in np.atleast_1d(axs).ravel():
        ax.axis("off")
    for ax, (name, sc, tile) in zip(np.atleast_1d(axs).ravel(), sel):
        ax.imshow(tile)
        ax.set_title(f"{name}  {sc:.2f}", fontsize=8)
    fig.suptitle(f"{slide_id}: CellScape DNA (magenta) over Xenium cells (green), orientation {chosen}; "
                 f"white = tissue on both. Score = shape match (1 = identical)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(out_dir, f"{slide_id}_match_gallery.png"), dpi=150)


if __name__ == "__main__":
    main()
