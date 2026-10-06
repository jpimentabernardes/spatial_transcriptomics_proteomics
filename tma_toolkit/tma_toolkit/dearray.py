"""TMA dearraying: find cores from cell centroids and give them grid positions.

Built for large arrays (e.g. 10 x 20 = 200 cores):
  1. DBSCAN on cell centroids (adaptive radius) -> core seeds; fragments of a torn
     core are merged back.
  2. The array's rotation is estimated from the nearest-neighbour directions
     between core centres and removed, so a TMA mounted a few degrees off still
     gives clean rows/columns.
  3. Rows/columns are found from GAPS in the rotated centres and indexed by the
     grid pitch -- a missing core, or even a whole missing row inside the array,
     does not shift the labels of the others.
  4. Grid positions are matched to core_row/core_col of the metadata sheet.

Always check the labelled plot (`plot_dearray`) against the physical TMA map.
Wrong labels can be fixed with an overrides table (slide_id, detected_core,
core_row, core_col) or by drawing core polygons (`assign_cells_to_polygons`).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import cKDTree


def detect_cores(xy: np.ndarray, eps: float | None = None, min_pts: int = 10, min_cells: int = 200,
                 merge_dist: float = 400.0) -> pd.DataFrame:
    """Cell centroids (N x 2, um) -> one row per core: detected_core, x_center, y_center, radius, n_cells."""
    from sklearn.cluster import DBSCAN
    xy = np.asarray(xy, float)
    if eps is None:
        d, _ = cKDTree(xy).query(xy, k=min_pts + 1)
        eps = 3 * float(np.median(d[:, -1]))
        print(f"DBSCAN eps (adaptive): {eps:.1f} um")
    lab = DBSCAN(eps=eps, min_samples=min_pts).fit_predict(xy)
    ok = lab >= 0
    if not ok.any():
        raise ValueError("DBSCAN found no clusters; set eps explicitly")
    cl = pd.DataFrame({"cl": lab[ok], "x": xy[ok, 0], "y": xy[ok, 1]})
    cent = cl.groupby("cl").agg(x=("x", "median"), y=("y", "median"), n=("x", "size"))

    seeds = cent[cent.n >= min_cells]
    if seeds.empty:
        raise ValueError("no cluster reaches min_cells")
    group = (fcluster(linkage(seeds[["x", "y"]].to_numpy(), "single"), t=merge_dist, criterion="distance")
             if len(seeds) > 1 else np.array([1]))
    lut = dict(zip(seeds.index, group))
    small = cent[cent.n < min_cells]
    if len(small):
        d, i = cKDTree(seeds[["x", "y"]].to_numpy()).query(small[["x", "y"]].to_numpy())
        lut.update({c: group[j] for c, j, dd in zip(small.index, i, d) if dd <= merge_dist})
    core = np.array([lut.get(v, -1) if v >= 0 else -1 for v in lab])

    rows = []
    for g in np.unique(core[core >= 0]):
        p = xy[core == g]
        cx, cy = np.median(p, 0)
        rows.append({"x_center": cx, "y_center": cy,
                     "radius": float(np.quantile(np.hypot(p[:, 0] - cx, p[:, 1] - cy), 0.95)),
                     "n_cells": len(p)})
    cores = pd.DataFrame(rows).sort_values(["y_center", "x_center"]).reset_index(drop=True)
    cores.insert(0, "detected_core", np.arange(1, len(cores) + 1))
    return cores


def estimate_rotation(centers: np.ndarray) -> float:
    """Array rotation in degrees (-45..45) from nearest-neighbour vectors between cores."""
    if len(centers) < 4:
        return 0.0
    _, idx = cKDTree(centers).query(centers, k=2)
    v = centers[idx[:, 1]] - centers
    ang = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
    ang = (ang + 45) % 90 - 45                    # fold onto one grid axis
    return float(np.median(ang))


def _axis_index(v: np.ndarray, tol: float) -> np.ndarray:
    """1D positions -> integer grid index, splitting at gaps > tol and using the pitch
    so an empty row/column in the middle leaves a hole instead of shifting labels."""
    order = np.argsort(v)
    groups = np.zeros(len(v), int)
    g = 0
    for a, b in zip(order[:-1], order[1:]):
        if v[b] - v[a] > tol:
            g += 1
        groups[b] = g
    means = np.array([v[groups == k].mean() for k in range(g + 1)])
    if len(means) == 1:
        return np.zeros(len(v), int)
    pitch = np.median(np.diff(means))
    idx_of_group = np.round((means - means[0]) / pitch).astype(int)
    return idx_of_group[groups]


def assign_grid(cores: pd.DataFrame, row_levels: list, col_levels: list, row_flip: bool = False,
                col_flip: bool = False, rotation_deg: float | None = None) -> pd.DataFrame:
    """Add core_row / core_col to detected cores. Row 1 (e.g. 'A') is at the TOP of the image
    (smallest y) and column 1 at the LEFT, unless flipped."""
    c = cores.copy()
    xy = c[["x_center", "y_center"]].to_numpy()
    rot = estimate_rotation(xy) if rotation_deg is None else rotation_deg
    t = np.radians(-rot)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    uv = (xy - xy.mean(0)) @ R.T
    tol = float(np.median(c["radius"]))
    r = _axis_index(uv[:, 1], tol)
    k = _axis_index(uv[:, 0], tol)
    if row_flip:
        r = r.max() - r
    if col_flip:
        k = k.max() - k
    if r.max() + 1 > len(row_levels) or k.max() + 1 > len(col_levels):
        print(f"WARNING: detected {r.max() + 1} rows x {k.max() + 1} cols but the sheet has "
              f"{len(row_levels)} x {len(col_levels)} -- extra cores are left unlabelled")
    if r.max() + 1 < len(row_levels) or k.max() + 1 < len(col_levels):
        print(f"NOTE: detected {r.max() + 1} rows x {k.max() + 1} cols, sheet has {len(row_levels)} x "
              f"{len(col_levels)}. If a whole EDGE row/column is missing, labels start one position off "
              "-- check the plot and use overrides.")
    c["core_row"] = [row_levels[i] if i < len(row_levels) else None for i in r]
    c["core_col"] = [col_levels[i] if i < len(col_levels) else None for i in k]
    c["grid_rotation_deg"] = rot
    return c


def sort_levels(values) -> list:
    v = pd.unique(pd.Series(values).dropna().astype(str))
    num = pd.to_numeric(pd.Series(v), errors="coerce")
    if num.notna().all():
        return [x for _, x in sorted(zip(num, v))]
    return sorted(v, key=lambda s: (len(s), s))       # A..Z before AA


def apply_overrides(cores: pd.DataFrame, overrides: pd.DataFrame | None, slide_id: str) -> pd.DataFrame:
    if overrides is None or overrides.empty:
        return cores
    c = cores.copy()
    for o in overrides[overrides.slide_id.astype(str) == str(slide_id)].itertuples(index=False):
        m = c.detected_core == int(o.detected_core)
        c.loc[m, "core_row"] = str(o.core_row)
        c.loc[m, "core_col"] = str(o.core_col)
    return c


def assign_cells_to_cores(xy: np.ndarray, cores: pd.DataFrame, factor: float = 1.25) -> np.ndarray:
    """Nearest core centre within factor x radius; returns core_id (object array, None outside)."""
    d, i = cKDTree(cores[["x_center", "y_center"]].to_numpy()).query(np.asarray(xy, float))
    ok = d <= cores["radius"].to_numpy()[i] * factor
    out = np.full(len(xy), None, dtype=object)
    out[ok] = cores["core_id"].to_numpy()[i[ok]]
    return out


def assign_cells_to_polygons(xy: np.ndarray, polygons: dict) -> np.ndarray:
    """polygons: {core_id: (N x 2 array of vertices in um)}; first match wins."""
    from matplotlib.path import Path
    out = np.full(len(xy), None, dtype=object)
    for cid, verts in polygons.items():
        m = Path(np.asarray(verts)).contains_points(xy) & (out == None)  # noqa: E711
        out[m] = cid
    return out


def dearray_slide(adata, core_meta: pd.DataFrame, slide_id: str, row_flip=False, col_flip=False,
                  overrides: pd.DataFrame | None = None, min_counts: int = 10, **detect_kw) -> pd.DataFrame:
    """Detect + label the cores of one slide and write core_id into adata.obs (in place).
    Returns the core table (geometry + core_id)."""
    in_slide = (adata.obs["slide_id"] == slide_id).to_numpy()
    use = in_slide & (adata.obs["n_counts"].to_numpy() >= min_counts)
    cores = detect_cores(adata.obsm["spatial"][use], **detect_kw)
    meta_s = core_meta[core_meta.slide_id == slide_id]
    cores = assign_grid(cores, sort_levels(meta_s.core_row), sort_levels(meta_s.core_col), row_flip, col_flip)
    cores = apply_overrides(cores, overrides, slide_id)
    cores = cores.merge(meta_s[["core_row", "core_col", "core_id"]], on=["core_row", "core_col"], how="left")
    cores["slide_id"] = slide_id
    n_unmatched = cores.core_id.isna().sum()
    if n_unmatched:
        print(f"{slide_id}: {n_unmatched} detected cores have no metadata row (detected_core "
              f"{cores.loc[cores.core_id.isna(), 'detected_core'].tolist()}) -- their cells are dropped")
    missing = sorted(set(meta_s.core_id) - set(cores.core_id.dropna()))
    if missing:
        print(f"{slide_id}: {len(missing)} cores in the sheet were not found on the slide: {missing[:20]}")
    labelled = cores.dropna(subset=["core_id"])
    if "core_id" not in adata.obs:
        adata.obs["core_id"] = None
    col = adata.obs["core_id"].astype(object).to_numpy()
    col[in_slide] = assign_cells_to_cores(adata.obsm["spatial"][in_slide], labelled)
    adata.obs["core_id"] = col
    return cores
