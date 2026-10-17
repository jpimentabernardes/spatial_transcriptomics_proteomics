#!/usr/bin/env python
"""03f_segmentation_disagreement.py -- where do 10x, cellpose and segger disagree?

Scans the cores in small tiles and asks, per tile, whether the three
segmentations group the SAME transcripts into the same cells. The tiles where
they disagree most are then shown close up (default 50 um):

  row 1  DAPI + the three outlines | 10x | cellpose | segger:
         outlines + every transcript, coloured by the cell it was assigned to
         (grey dot = not assigned to any cell)
  row 2  DAPI alone | 10x | cellpose | segger: the same transcripts coloured by
         the lineage of their gene (L1 markers of config/marker_panel_rna.csv),
         so you see e.g. immune transcripts pulled into an epithelial cell

A few tiles where the methods agree best are drawn the same way, as a reference.

Disagreement per tile (transcripts = qv >= 20 gene transcripts in the tile):
  ARI(A, B)     adjusted Rand index between two methods' cell assignments, on the
                transcripts both assigned: 1 = same grouping into cells, ~0 = unrelated
  score         mean over method pairs of (1 - ARI) + 0.5 x (highest - lowest
                fraction of transcripts assigned). High = the methods split/merge
                cells differently AND/OR one leaves many transcripts unassigned.

Inputs: Xenium outs (transcripts.parquet, cell_boundaries.parquet, DAPI),
  02's core table, cellpose output (03a; per-core masks in <cellpose>/masks/ from
  03a --save-masks), segger output (segger_segmentation.parquet, 03b).
  Without saved cellpose masks the scan uses 10x vs segger only and cellpose is
  run on each zoom window with 03a's settings (needs cellpose installed).
  Segger outlines: `segger export boundaries` GeoParquet if --segger-boundaries
  is given, else the convex hull of each segger cell's transcripts (segger
  assigns transcripts, it does not draw cell borders).

Outputs (in --out-dir, default <project>/figures/03f):
  03f_<slide>_disagree_<rank>_<core>.png, 03f_<slide>_agree_<rank>_<core>.png
  03f_<slide>_<...>_transcripts.csv     per zoom: x, y, gene, cell per method
  03f_<slide>_tiles.csv                 every scanned tile with its scores
  03f_<slide>_disagreement_by_organ.png

Example:
  python 03f_segmentation_disagreement.py --slide-id TMA_ORGAN --project-dir /work_ikmb/.../work
  (Xenium outs from config/slides.csv; or pass --xenium-dir /path/output-XETG...)
"""

import argparse
import glob
import importlib
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "py"))
from xenium_io import (PyramidImage, _as_str, feature_type, find_morphology_channels,  # noqa: E402
                       pixel_size_um)

METHODS = ["xenium", "cellpose", "segger"]
METHOD_LABEL = {"xenium": "10x Xenium", "cellpose": "cellpose", "segger": "segger"}
METHOD_COLOURS = {"xenium": "#00a6c7", "cellpose": "#e0a800", "segger": "#b48ee6"}     # as in 03c
LINEAGE_COLOURS = {"Epithelial": "#418ED0", "Immune": "#C17743", "Endothelial": "#D26496",
                   "Mesenchymal": "#00A064", "Neural_NE": "#5D981B"}                     # as in 04
CELL_PALETTE = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3",
                "#808000", "#ffd8b1", "#000075", "#ff69b4"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--xenium-dir", default=None,
                   help="Xenium outs of the slide (default: xenium_dir of the slide in config/slides.csv)")
    p.add_argument("--project-dir", default=os.environ.get("EXP1_PROJECT_DIR"),
                   help="project folder (default: $EXP1_PROJECT_DIR); the paths below default into it")
    p.add_argument("--cores-csv", help="default <project>/tables/02_cores_all.csv")
    p.add_argument("--seg-dir", help="default <project>/data/segmentation")
    p.add_argument("--segger-parquet", help="segger_segmentation.parquet or a folder with it "
                                            "(default <seg-dir>/segger_raw/<slide>)")
    p.add_argument("--segger-boundaries", default=None, help="GeoParquet from `segger export boundaries` (optional)")
    p.add_argument("--tma-map", help="default <code>/config/tma_map.csv (organ per core, for the summary)")
    p.add_argument("--markers-csv", default=os.path.join(HERE, "config", "marker_panel_rna.csv"))
    p.add_argument("--out-dir", help="default <project>/figures/03f")
    p.add_argument("--cores", nargs="*", help="only these cores (default: cores segmented by all methods)")
    p.add_argument("--max-cores", type=int, default=24,
                   help="scan at most this many cores, spread over organs (memory / time)")
    p.add_argument("--tile-um", type=float, default=50.0, help="scan tile size")
    p.add_argument("--window-um", type=float, default=50.0, help="zoom window size")
    p.add_argument("--n-windows", type=int, default=12, help="disagreement windows to draw")
    p.add_argument("--n-agree", type=int, default=2, help="agreement (reference) windows to draw")
    p.add_argument("--max-per-core", type=int, default=2)
    p.add_argument("--min-tx-tile", type=int, default=200, help="tiles with fewer transcripts are ignored")
    p.add_argument("--min-qv", type=float, default=20.0)
    p.add_argument("--cellpose-mode", choices=["nuclei_expand", "cyto"], default="nuclei_expand",
                   help="only for the fallback (no saved 03a masks): same as 03a")
    p.add_argument("--cellpose-expand-um", type=float, default=5.0)
    a = p.parse_args(argv)
    if not a.project_dir and not (a.cores_csv and a.seg_dir and a.out_dir):
        p.error("give --project-dir (or set EXP1_PROJECT_DIR), or all of --cores-csv / --seg-dir / --out-dir")
    a.cores_csv = a.cores_csv or os.path.join(a.project_dir, "tables", "02_cores_all.csv")
    a.seg_dir = a.seg_dir or os.path.join(a.project_dir, "data", "segmentation")
    a.out_dir = a.out_dir or os.path.join(a.project_dir, "figures", "03f")
    a.segger_parquet = a.segger_parquet or os.path.join(a.seg_dir, "segger_raw", a.slide_id)
    a.tma_map = a.tma_map or os.path.join(HERE, "config", "tma_map.csv")
    if not a.xenium_dir:
        sl = pd.read_csv(os.path.join(HERE, "config", "slides.csv"))
        row = sl[sl["slide_id"] == a.slide_id]
        if row.empty:
            p.error(f"{a.slide_id} is not in config/slides.csv -- give --xenium-dir")
        d = str(row["xenium_dir"].iloc[0])
        a.xenium_dir = d if os.path.isabs(d) else os.path.join(a.project_dir or "", d)
        print(f"Xenium outs (config/slides.csv): {a.xenium_dir}")
    if not os.path.isdir(a.xenium_dir):
        p.error(f"Xenium folder not found: {a.xenium_dir}")
    return a


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def read_cells(seg_dir, method, slide):
    f = os.path.join(seg_dir, method, slide, "cells.csv.gz")
    return pd.read_csv(f) if os.path.exists(f) else None


def pick_cores(a, cores, organ_of):
    """Cores in every method's output; at most --max-cores, taken round-robin over organs."""
    ids = list(cores["core_id"])
    for m in ("cellpose", "segger"):
        c = read_cells(a.seg_dir, m, a.slide_id)
        if c is not None and "core_id" in c:
            ids = [i for i in ids if i in set(c["core_id"].dropna())]
        else:
            print(f"no {m} cells.csv.gz for {a.slide_id} -- {m} left out of the scan")
    if a.cores:
        ids = [i for i in ids if i in set(a.cores)]
    if len(ids) > a.max_cores:
        by_org = {}
        for i in ids:
            by_org.setdefault(organ_of.get(i, "?"), []).append(i)
        rr, k = [], 0
        while len(rr) < a.max_cores:
            for o in sorted(by_org):
                if k < len(by_org[o]) and len(rr) < a.max_cores:
                    rr.append(by_org[o][k])
            k += 1
        ids = rr
    return cores[cores["core_id"].isin(ids)].reset_index(drop=True)


def read_xenium_tx(xenium_dir, boxes, min_qv):
    """qv-filtered gene transcripts inside the core boxes, with the 10x cell id."""
    import pyarrow.dataset as ds
    path = os.path.join(xenium_dir, "transcripts.parquet")
    d = ds.dataset(path)
    cols = ["x_location", "y_location", "feature_name", "qv", "cell_id"]
    f = ((ds.field("qv") >= min_qv) &
         (ds.field("x_location") >= float(boxes.xmin.min())) & (ds.field("x_location") <= float(boxes.xmax.max())) &
         (ds.field("y_location") >= float(boxes.ymin.min())) & (ds.field("y_location") <= float(boxes.ymax.max())))
    parts = []
    for batch in d.to_batches(columns=cols, filter=f, batch_size=5_000_000):
        x = batch.column("x_location").to_numpy()
        y = batch.column("y_location").to_numpy()
        hit = np.full(len(x), -1, dtype=np.int64)
        for k, b in enumerate(boxes.itertuples(index=False)):
            m = (hit < 0) & (x >= b.xmin) & (x < b.xmax) & (y >= b.ymin) & (y < b.ymax)
            hit[m] = k
        keep = hit >= 0
        if not keep.any():
            continue
        idx = np.where(keep)[0]
        names = np.asarray(_as_str(batch.column("feature_name")).take(idx).to_pylist(), dtype=object)
        cid = _as_str(batch.column("cell_id")).take(idx).to_pylist()
        parts.append(pd.DataFrame({"core_id": boxes.core_id.to_numpy()[hit[idx]], "x": x[idx], "y": y[idx],
                                   "gene": names, "xenium": [str(v) for v in cid]}))
    tx = pd.concat(parts, ignore_index=True)
    tx = tx[[feature_type(g) == "Gene Expression" for g in tx["gene"]]].reset_index(drop=True)
    tx.loc[tx["xenium"].isin(["UNASSIGNED", "", "-1", "None", "nan"]), "xenium"] = ""
    print(f"{len(tx):,} gene transcripts (qv >= {min_qv}) in {len(boxes)} cores")
    return tx


def tx_key(x, y, gene):
    return pd.Series(np.round(np.asarray(x, float) * 100).astype(np.int64)).astype(str) + "_" + \
           pd.Series(np.round(np.asarray(y, float) * 100).astype(np.int64)).astype(str) + "_" + pd.Series(gene).astype(str)


def add_segger(tx, path, boxes):
    """segger cell per transcript (only `filtered` == True assignments, as in 03b)."""
    import pyarrow.dataset as ds
    s03b = importlib.import_module("03b_segger_to_common")
    try:
        path = s03b.find_segger_parquet(path)
    except FileNotFoundError as e:
        print(f"segger: {e}\n -> segger left out")
        return tx, False
    d = ds.dataset(path)
    xc, yc = ("x", "y") if "x" in d.schema.names else ("x_location", "y_location")
    f = ((ds.field(xc) >= float(boxes.xmin.min())) & (ds.field(xc) <= float(boxes.xmax.max())) &
         (ds.field(yc) >= float(boxes.ymin.min())) & (ds.field(yc) <= float(boxes.ymax.max())))
    s = d.to_table(columns=[xc, yc, "feature_name", "segger_cell_id", "filtered"], filter=f).to_pandas()
    s["feature_name"] = s["feature_name"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
    ok = s["filtered"].fillna(False).astype(bool) & s["segger_cell_id"].notna()
    s["segger"] = np.where(ok, s["segger_cell_id"].astype(str), "")
    s["key"] = tx_key(s[xc], s[yc], s["feature_name"]).to_numpy()
    s = s.drop_duplicates("key")
    tx["key"] = tx_key(tx["x"], tx["y"], tx["gene"]).to_numpy()
    tx = tx.merge(s[["key", "segger"]], on="key", how="left")
    matched = tx["segger"].notna().mean()
    tx["segger"] = tx["segger"].fillna("")
    print(f"segger: {matched:.1%} of the transcripts matched to segger's output")
    if matched < 0.5:
        print("  WARNING: few transcripts matched -- is the segger output from this slide / Xenium run?")
    return tx.drop(columns="key"), True


def add_cellpose(tx, seg_dir, slide):
    """cellpose cell per transcript from 03a's saved masks (masks/<core>.npz); cores without masks -> NaN."""
    cells = read_cells(seg_dir, "cellpose", slide)
    mdir = os.path.join(seg_dir, "cellpose", slide, "masks")
    tx["cellpose"] = pd.Series(np.nan, index=tx.index, dtype=object)
    if cells is None or not os.path.isdir(mdir):
        print(f"cellpose: no saved masks in {mdir} (03a --save-masks) -> scan without cellpose; "
              "cellpose is run on the zoom windows instead")
        return tx, set()
    valid = set(cells["cell_id"].astype(str))
    have = set()
    for core, idx in tx.groupby("core_id").groups.items():
        f = os.path.join(mdir, f"{core}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        masks, (r0, c0), px = z["masks"], z["offset_rc"], float(z["pixel_size"])
        sub = tx.loc[idx]
        cc = np.floor(sub["x"].to_numpy() / px).astype(np.int64) - int(c0)
        rr = np.floor(sub["y"].to_numpy() / px).astype(np.int64) - int(r0)
        inb = (cc >= 0) & (cc < masks.shape[1]) & (rr >= 0) & (rr < masks.shape[0])
        lab = np.zeros(len(sub), dtype=np.int64)
        lab[inb] = masks[rr[inb], cc[inb]]
        ids = np.array([f"{core}-{v}" if v > 0 else "" for v in lab], dtype=object)
        ids[~np.isin(ids, list(valid))] = ""
        tx.loc[idx, "cellpose"] = ids
        have.add(core)
    print(f"cellpose: masks for {len(have)} of {tx.core_id.nunique()} scanned cores")
    return tx, have


# ---------------------------------------------------------------------------
# disagreement
# ---------------------------------------------------------------------------
def ari(a, b):
    """Adjusted Rand index of two labelings (numpy, no sklearn)."""
    if len(a) < 2:
        return np.nan
    ia, ib = pd.factorize(a)[0], pd.factorize(b)[0]
    nb = ib.max() + 1
    nij = np.bincount(ia * nb + ib)
    c2 = lambda n: n * (n - 1) / 2.0  # noqa: E731
    s_ij = c2(nij[nij > 1]).sum()
    s_a = c2(np.bincount(ia)).sum()
    s_b = c2(np.bincount(ib)).sum()
    exp = s_a * s_b / c2(len(a))
    mx = (s_a + s_b) / 2
    return 1.0 if mx == exp else (s_ij - exp) / (mx - exp)


def scan_tiles(tx, boxes, methods, tile, min_tx):
    rows = []
    for b in boxes.itertuples(index=False):
        t = tx[tx.core_id == b.core_id]
        if t.empty:
            continue
        ms = [m for m in methods if t[m].notna().all()]
        if len(ms) < 2:
            continue
        gx = np.floor((t.x.to_numpy() - b.xmin) / tile).astype(int)
        gy = np.floor((t.y.to_numpy() - b.ymin) / tile).astype(int)
        for (i, j), idx in t.groupby([gx, gy]).groups.items():
            s = t.loc[idx]
            if len(s) < min_tx:
                continue
            r = {"core_id": b.core_id, "x0": b.xmin + i * tile, "y0": b.ymin + j * tile, "n_tx": len(s),
                 "methods": "+".join(ms)}
            fr = {}
            for m in ms:
                asg = s[m].to_numpy() != ""
                fr[m] = asg.mean()
                r[f"frac_assigned_{m}"] = round(fr[m], 3)
                r[f"n_cells_{m}"] = int(pd.unique(s[m][asg]).size)
            aris = []
            for k1 in range(len(ms)):
                for k2 in range(k1 + 1, len(ms)):
                    m1, m2 = ms[k1], ms[k2]
                    both = (s[m1].to_numpy() != "") & (s[m2].to_numpy() != "")
                    v = ari(s[m1].to_numpy()[both], s[m2].to_numpy()[both]) if both.sum() >= 20 else np.nan
                    r[f"ari_{m1}_{m2}"] = round(v, 3) if np.isfinite(v) else np.nan
                    aris.append(v)
            aris = [v for v in aris if np.isfinite(v)]
            if not aris:
                continue
            r["mean_ari"] = round(float(np.mean(aris)), 3)
            r["assign_spread"] = round(max(fr.values()) - min(fr.values()), 3)
            r["score"] = round(float(np.mean([1 - v for v in aris])) + 0.5 * r["assign_spread"], 3)
            rows.append(r)
    return pd.DataFrame(rows)


def choose(tiles, n, per_core, ascending, min_dist):
    """Top tiles by score, at most per_core per core and min_dist apart within a core."""
    out, per = [], {}
    if n <= 0:
        return out
    for r in tiles.sort_values("score", ascending=ascending).itertuples(index=False):
        if per.get(r.core_id, 0) >= per_core:
            continue
        if any(o.core_id == r.core_id and np.hypot(o.x0 - r.x0, o.y0 - r.y0) < min_dist for o in out):
            continue
        out.append(r)
        per[r.core_id] = per.get(r.core_id, 0) + 1
        if len(out) >= n:
            break
    return out


# ---------------------------------------------------------------------------
# outlines
# ---------------------------------------------------------------------------
def mask_outlines(crop_region, masks, r0, c0, px, x0, y0, x1, y1, id_fmt, margin=15.0):
    from skimage.measure import regionprops_table
    p = pd.DataFrame(regionprops_table(masks, properties=("label", "centroid")))
    if p.empty:
        return pd.DataFrame(columns=["cell_id", "vertex_x", "vertex_y"])
    p["x"] = (c0 + p["centroid-1"] + 0.5) * px
    p["y"] = (r0 + p["centroid-0"] + 0.5) * px
    p = p[p.x.between(x0 - margin, x1 + margin) & p.y.between(y0 - margin, y1 + margin)]
    ids = {int(v): id_fmt(int(v)) for v in p["label"]}
    return crop_region.outlines(masks, p["label"].to_numpy(), c0, r0, px, ids)


def hull_outlines(t, col):
    from scipy.spatial import ConvexHull
    rows = []
    for cid, g in t[t[col] != ""].groupby(col):
        pts = np.unique(g[["x", "y"]].to_numpy(), axis=0)
        if len(pts) < 3:
            continue
        try:
            h = ConvexHull(pts)
        except Exception:
            continue
        v = pts[h.vertices]
        rows.append(pd.DataFrame({"cell_id": cid, "vertex_x": v[:, 0], "vertex_y": v[:, 1]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["cell_id", "vertex_x", "vertex_y"])


def geoparquet_outlines(path, x0, y0, x1, y1, margin=15.0):
    import geopandas as gpd
    b = gpd.read_parquet(path)
    b = b.cx[x0 - margin:x1 + margin, y0 - margin:y1 + margin]
    id_col = next(c for c in b.columns if c != b.geometry.name)
    rows = []
    for cid, geom in zip(b[id_col].astype(str), b.geometry):
        polys = getattr(geom, "geoms", [geom])
        for poly in polys:
            xy = np.asarray(poly.exterior.coords)
            rows.append(pd.DataFrame({"cell_id": cid, "vertex_x": xy[:, 0], "vertex_y": xy[:, 1]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["cell_id", "vertex_x", "vertex_y"])


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------
def cell_colour(cid):
    """Same colour for the same cell id in every run (python's hash() is salted per process)."""
    import zlib
    return CELL_PALETTE[zlib.crc32(str(cid).encode()) % len(CELL_PALETTE)]


def draw_window(path, title, img, win, tw, outl, avail, gene_lineage, notes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    x0, y0, x1, y1 = win
    lo, hi = np.percentile(img, [1, 99.7]) if img.size else (0, 1)
    g = np.clip((img.astype(float) - lo) / max(hi - lo, 1), 0, 1) ** 0.6
    fig, axes = plt.subplots(2, 4, figsize=(22, 12.8), layout="constrained")
    for ax in axes.ravel():
        ax.imshow(g, cmap="gray", extent=(x0, x1, y1, y0), interpolation="nearest")
        ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
        ax.set_xticks([]); ax.set_yticks([])

    def draw_outl(ax, b, colour, lw=1.1, alpha=1.0):
        if b is None or b.empty:
            return
        for _, v in b.groupby("cell_id", sort=False):
            xy = v[["vertex_x", "vertex_y"]].to_numpy()
            xy = np.vstack([xy, xy[:1]])
            ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=lw, alpha=alpha)

    # row 1, col 1: all outlines
    ax = axes[0, 0]
    for m in METHODS:
        if m in avail:
            draw_outl(ax, outl.get(m), METHOD_COLOURS[m], lw=1.0)
    ax.legend(handles=[Line2D([0], [0], color=METHOD_COLOURS[m], lw=2, label=METHOD_LABEL[m]) for m in METHODS if m in avail],
              loc="lower right", fontsize=9, framealpha=0.8)
    ax.set_title("DAPI + outlines of all methods", fontsize=11)
    axes[1, 0].set_title("DAPI", fontsize=11)
    bar = 10.0 if x1 - x0 >= 30 else 5.0
    axes[1, 0].plot([x0 + 1.5, x0 + 1.5 + bar], [y1 - 2, y1 - 2], color="white", lw=3)
    axes[1, 0].text(x0 + 1.5, y1 - 3, f"{bar:.0f} um", color="white", fontsize=10, va="bottom")

    lin = tw["gene"].map(gene_lineage).fillna("other").to_numpy()
    for k, m in enumerate(METHODS):
        a1, a2 = axes[0, k + 1], axes[1, k + 1]
        if m not in avail:
            for ax in (a1, a2):
                ax.set_title(f"{METHOD_LABEL[m]}: not available", fontsize=11)
            continue
        lab = tw[m].fillna("").to_numpy()
        asg = lab != ""
        for ax in (a1, a2):
            draw_outl(ax, outl.get(m), "white", lw=1.0, alpha=0.9)
            ax.scatter(tw.x[~asg], tw.y[~asg], s=4, c="#9a9a9a", linewidths=0, alpha=0.8)
        a1.scatter(tw.x[asg], tw.y[asg], s=7, c=[cell_colour(c) for c in lab[asg]], linewidths=0)
        cols = np.array([LINEAGE_COLOURS.get(v, "#d9d9d9") for v in lin])
        mk = lin != "other"
        a2.scatter(tw.x[asg & ~mk], tw.y[asg & ~mk], s=3, c="#d9d9d9", linewidths=0, alpha=0.6)
        a2.scatter(tw.x[asg & mk], tw.y[asg & mk], s=12, c=cols[asg & mk], linewidths=0)
        n_cells = len(pd.unique(lab[asg]))
        med = pd.Series(lab[asg]).value_counts().median() if asg.any() else 0
        a1.set_title(f"{METHOD_LABEL[m]}: {n_cells} cells, {100 * asg.mean():.0f}% of transcripts assigned, "
                     f"median {med:.0f}/cell\n{notes.get(m, '')}", fontsize=10)
        a2.set_title(f"{METHOD_LABEL[m]}: transcripts by lineage of their gene", fontsize=10)
    axes[1, 3].legend(handles=[Line2D([0], [0], marker="o", ls="", color=c, label=l) for l, c in LINEAGE_COLOURS.items()] +
                      [Line2D([0], [0], marker="o", ls="", color="#d9d9d9", label="other gene"),
                       Line2D([0], [0], marker="o", ls="", color="#9a9a9a", label="not assigned")],
                      loc="lower right", fontsize=8, framealpha=0.85)
    fig.suptitle(title, fontsize=13)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def summary_plot(path, tiles, organ_of, slide):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = tiles.copy()
    t["organ"] = t["core_id"].map(organ_of).fillna("?")
    order = t.groupby("organ")["score"].median().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(order) + 3), 6))
    ax.boxplot([t.loc[t.organ == o, "score"] for o in order], showfliers=False)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order)
    for k, o in enumerate(order):
        v = t.loc[t.organ == o, "score"].to_numpy()
        ax.scatter(np.random.default_rng(1).normal(k + 1, 0.06, len(v)), v, s=4, alpha=0.3, color="#8e6bbf")
    ax.set_ylabel("disagreement score per tile\n(mean 1 - ARI + 0.5 x assignment spread)")
    ax.set_title(f"{slide}: where do the segmentations disagree? (one dot = one {int(tiles.attrs.get('tile', 50))} um tile)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main(argv=None):
    a = parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    sys.path.insert(0, HERE)
    crop_region = importlib.import_module("crop_region")

    cores = pd.read_csv(a.cores_csv)
    cores = cores[(cores.slide_id == a.slide_id) & cores.core_id.notna()].reset_index(drop=True)
    organ_of = {}
    if os.path.exists(a.tma_map):
        tm = pd.read_csv(a.tma_map, encoding="latin1")
        organ_of = dict(zip(tm.core_id, tm.organ))
    boxes = pick_cores(a, cores, organ_of)
    if boxes.empty:
        sys.exit("No cores to scan (check --cores / the cellpose and segger outputs of this slide).")
    print(f"scanning {len(boxes)} cores: {', '.join(boxes.core_id)}")

    tx = read_xenium_tx(a.xenium_dir, boxes, a.min_qv)
    tx, has_segger = add_segger(tx, a.segger_parquet, boxes)
    if not has_segger:
        tx["segger"] = np.nan
    tx, cp_cores = add_cellpose(tx, a.seg_dir, a.slide_id)

    tiles = scan_tiles(tx, boxes, METHODS, a.tile_um, a.min_tx_tile)
    if tiles.empty:
        sys.exit("No tile could be compared (fewer than two methods, or too few transcripts per tile).")
    tiles["organ"] = tiles.core_id.map(organ_of)
    tiles.attrs["tile"] = a.tile_um
    tiles.to_csv(os.path.join(a.out_dir, f"03f_{a.slide_id}_tiles.csv"), index=False)
    summary_plot(os.path.join(a.out_dir, f"03f_{a.slide_id}_disagreement_by_organ.png"), tiles, organ_of, a.slide_id)
    print(f"{len(tiles)} tiles scored; median score {tiles.score.median():.2f}")

    # gene -> lineage (L1 markers that belong to exactly one lineage)
    gene_lineage = {}
    if os.path.exists(a.markers_csv):
        mk = pd.read_csv(a.markers_csv)
        mk = mk[mk["level"].astype(str) == "1"]
        cnt = {}
        for r in mk.itertuples(index=False):
            for g in str(r.markers).split(";"):
                cnt.setdefault(g.strip(), set()).add(r.cell_type)
        gene_lineage = {g: next(iter(v)) for g, v in cnt.items() if len(v) == 1 and next(iter(v)) in LINEAGE_COLOURS}

    px = pixel_size_um(a.xenium_dir)
    dapi = PyramidImage(find_morphology_channels(a.xenium_dir)["dapi"], level=0)
    picks = [("disagree", r) for r in choose(tiles, a.n_windows, a.max_per_core, False, a.window_um)]
    picks += [("agree", r) for r in choose(tiles[tiles.n_tx >= tiles.n_tx.median()], a.n_agree, 1, True, a.window_um)]

    for rank, (kind, r) in enumerate(picks, 1):
        cx, cy = r.x0 + a.tile_um / 2, r.y0 + a.tile_um / 2
        h = a.window_um / 2
        x0, y0, x1, y1 = cx - h, cy - h, cx + h, cy + h
        c0, r0 = int(np.floor(x0 / px)), int(np.floor(y0 / px))
        c1, r1 = int(np.ceil(x1 / px)), int(np.ceil(y1 / px))
        img = dapi.read(r0, r1, c0, c1)
        if img.ndim == 3:
            img = img[..., 0]
        core_tx = tx[tx.core_id == r.core_id]
        tw = core_tx[core_tx.x.between(x0, x1) & core_tx.y.between(y0, y1)].copy()
        outl, avail, notes = {}, set(), {}

        # 10x
        outl["xenium"] = crop_region.read_boundaries(a.xenium_dir, "cell", x0, y0, x1, y1)
        avail.add("xenium")
        # cellpose: 03a masks, else run on the window
        f = os.path.join(a.seg_dir, "cellpose", a.slide_id, "masks", f"{r.core_id}.npz")
        if r.core_id in cp_cores and os.path.exists(f):
            z = np.load(f)
            m_, (mr0, mc0), mpx = z["masks"], z["offset_rc"], float(z["pixel_size"])
            outl["cellpose"] = mask_outlines(crop_region, m_, int(mr0), int(mc0), mpx, x0, y0, x1, y1,
                                             lambda v, c=r.core_id: f"{c}-{v}")
            avail.add("cellpose"); notes["cellpose"] = "03a result"
        else:
            try:
                ns = SimpleNamespace(x0=x0, y0=y0, x1=x1, y1=y1, cellpose_mode=a.cellpose_mode,
                                     cellpose_expand_um=a.cellpose_expand_um, cellpose_gpu=None,
                                     cellpose_diameter_um=None)
                m_, _, wr0, wc0, _, _, _, _ = crop_region.segment_window(ns, a.xenium_dir, px, 30.0)
                cc = np.floor(tw.x.to_numpy() / px).astype(int) - wc0
                rr = np.floor(tw.y.to_numpy() / px).astype(int) - wr0
                ok = (cc >= 0) & (cc < m_.shape[1]) & (rr >= 0) & (rr < m_.shape[0])
                lab = np.zeros(len(tw), dtype=np.int64)
                lab[ok] = m_[rr[ok], cc[ok]]
                tw["cellpose"] = [f"win-{v}" if v > 0 else "" for v in lab]
                outl["cellpose"] = mask_outlines(crop_region, m_, wr0, wc0, px, x0, y0, x1, y1, lambda v: f"win-{v}")
                avail.add("cellpose"); notes["cellpose"] = "cellpose run on this window (03a settings)"
            except Exception as e:
                print(f"  cellpose on window failed: {type(e).__name__}: {e}")
        # segger
        if has_segger:
            avail.add("segger")
            if a.segger_boundaries and os.path.exists(a.segger_boundaries):
                outl["segger"] = geoparquet_outlines(a.segger_boundaries, x0, y0, x1, y1)
                notes["segger"] = "segger boundaries"
            else:
                near = core_tx[core_tx.segger.isin(set(tw.segger) - {""})]
                outl["segger"] = hull_outlines(near, "segger")
                notes["segger"] = "outline = hull of the cell's transcripts"
        for m in METHODS:
            if m not in tw:
                tw[m] = np.nan

        organ = organ_of.get(r.core_id, "")
        aris = ", ".join(f"{k[4:].replace('_', ' vs ')} {getattr(r, k):.2f}" for k in r._fields
                         if k.startswith("ari_") and pd.notna(getattr(r, k)))
        what = "methods DISAGREE" if kind == "disagree" else "methods AGREE (reference)"
        title = (f"{a.slide_id}  {r.core_id} {organ}  |  x {x0:.0f}-{x1:.0f}, y {y0:.0f}-{y1:.0f} um  |  {what}: "
                 f"score {r.score:.2f}; ARI {aris}")
        stem = f"03f_{a.slide_id}_{kind}_{rank:02d}_{r.core_id}"
        draw_window(os.path.join(a.out_dir, stem + ".png"), title, img, (x0, y0, x1, y1), tw, outl, avail,
                    gene_lineage, notes)
        tw[["x", "y", "gene"] + [m for m in METHODS if m in tw]].to_csv(
            os.path.join(a.out_dir, stem + "_transcripts.csv"), index=False)
        print(f"[{rank}/{len(picks)}] {stem}.png")
    print(f"done -> {a.out_dir}")


if __name__ == "__main__":
    main()
