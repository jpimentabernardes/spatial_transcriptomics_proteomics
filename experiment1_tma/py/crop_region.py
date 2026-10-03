#!/usr/bin/env python
"""crop_region.py -- cut one window out of a Xenium run for close inspection.

Called from R (01_xenium_object_prep.R, section 6) but usable on its own:

    python py/crop_region.py --xenium-dir <outs> --x0 1000 --y0 2000 --x1 1150 --y1 2150 --out-prefix /tmp/win1

Writes
  <prefix>_image.tif           full-resolution crop of a morphology image (uint16, uncompressed,
                               so R's tiff package can read it); DAPI by default
  <prefix>_meta.json           pixel size and the crop's origin/extent in microns
  <prefix>_cell_boundaries.csv.gz     cell outlines (cell_id, vertex_x, vertex_y; microns) in the window
  <prefix>_nucleus_boundaries.csv.gz  nucleus outlines, same format (if present in the outs)

With --cellpose, cellpose is also run on the window (+ a 30 um margin) with the
same model and settings as 03a_segment_cellpose.py, its cells get transcript
counts and the same QC rules as 01 (thresholds passed in), and it writes
  <prefix>_cellpose_cells.csv                  cell_id, x_um, y_um, cell_area, n_counts, n_features,
                                               control_frac, qc_reason
  <prefix>_cellpose_cell_boundaries.csv.gz     outlines, same format as above
  <prefix>_cellpose_nucleus_boundaries.csv.gz  nuclei before expansion (nuclei_expand mode)
If cellpose fails (not installed, no GPU memory, ...) the crop is still written and
the error goes to <prefix>_meta.json ("cellpose_error").

Coordinates are Xenium microns: pixel = micron / pixel_size, origin top-left of the
morphology image -- the same frame as cells.parquet.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xenium_io import PyramidImage, feature_type, find_morphology_channels, pixel_size_um  # noqa: E402

# control types counted in 01's control_frac (nCount_ControlProbe / ControlCodeword / BlankCodeword)
QC_CONTROL_TYPES = {"Negative Control Probe", "Negative Control Codeword", "Blank Codeword"}


def read_boundaries(xenium_dir, kind, x0, y0, x1, y1, margin=15.0):
    """Outline vertices of cells/nuclei that have at least one vertex inside the window."""
    pq = os.path.join(xenium_dir, f"{kind}_boundaries.parquet")
    csv = os.path.join(xenium_dir, f"{kind}_boundaries.csv.gz")
    cols = ["cell_id", "vertex_x", "vertex_y"]
    if os.path.exists(pq):
        import pyarrow.dataset as ds
        f = ((ds.field("vertex_x") >= x0 - margin) & (ds.field("vertex_x") <= x1 + margin) &
             (ds.field("vertex_y") >= y0 - margin) & (ds.field("vertex_y") <= y1 + margin))
        near = ds.dataset(pq).to_table(columns=cols, filter=f).to_pandas()
        if near.empty:
            return near
        ids = near["cell_id"].unique()
        b = ds.dataset(pq).to_table(columns=cols, filter=ds.field("cell_id").isin(ids)).to_pandas()
    elif os.path.exists(csv):
        b = pd.read_csv(csv, usecols=cols)
        inside = b[(b.vertex_x.between(x0 - margin, x1 + margin)) & (b.vertex_y.between(y0 - margin, y1 + margin))]
        b = b[b.cell_id.isin(inside.cell_id.unique())]
    else:
        return None
    b["cell_id"] = b["cell_id"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
    return b


def outlines(masks, labels, c0, r0, px, ids, max_vertices=60):
    """Label image -> outline vertices in microns (one closed polygon per label)."""
    from skimage.measure import find_contours, regionprops
    keep = set(int(v) for v in labels)
    rows = []
    for rp in regionprops(masks):
        if rp.label not in keep:
            continue
        r_min, c_min, _, _ = rp.bbox
        m = np.pad(rp.image, 1).astype(float)
        cs = find_contours(m, 0.5)
        if not cs:
            continue
        c = max(cs, key=len)
        if len(c) > max_vertices:
            c = c[np.linspace(0, len(c) - 1, max_vertices).astype(int)]
        x = (c0 + c_min - 1 + c[:, 1] + 0.5) * px
        y = (r0 + r_min - 1 + c[:, 0] + 0.5) * px
        rows.append(pd.DataFrame({"cell_id": ids[rp.label], "vertex_x": x, "vertex_y": y}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["cell_id", "vertex_x", "vertex_y"])


def cellpose_window(a, xenium_dir, px, prefix, pad=30.0, margin=15.0):
    """Segment the window like 03a does and give each cellpose cell 01's QC flags."""
    import importlib
    import time
    from types import SimpleNamespace

    import pyarrow.compute as pc
    import pyarrow.dataset as ds
    from skimage.measure import regionprops_table
    from skimage.segmentation import expand_labels

    t0 = time.time()
    code_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, code_dir)
    s03a = importlib.import_module("03a_segment_cellpose")      # same model loading / call as 03a

    ch = find_morphology_channels(xenium_dir)
    dapi = PyramidImage(ch["dapi"], level=0)
    c0 = max(0, int(np.floor((a.x0 - pad) / px)))
    r0 = max(0, int(np.floor((a.y0 - pad) / px)))
    c1 = min(dapi.width, int(np.ceil((a.x1 + pad) / px)))
    r1 = min(dapi.height, int(np.ceil((a.y1 + pad) / px)))
    img = dapi.read(r0, r1, c0, c1)

    gpu = a.cellpose_gpu
    if gpu is None:
        try:
            import torch
            gpu = torch.cuda.is_available()
        except ImportError:
            gpu = False
    model, major = s03a.load_cellpose(a.cellpose_mode, gpu)
    opts = SimpleNamespace(flow_threshold=0.4, cellprob_threshold=0.0, force_diameter=False)
    diameter_um = a.cellpose_diameter_um or (8.0 if a.cellpose_mode == "nuclei_expand" else 15.0)
    nuclei = None
    if a.cellpose_mode == "nuclei_expand":
        nuclei = s03a.run_cellpose(model, major, img, diameter_um / px, opts)
        expand_px = int(round(a.cellpose_expand_um / px))
        masks = expand_labels(nuclei, distance=expand_px) if expand_px > 0 else nuclei
    else:
        if ch["boundary"] is None:
            raise RuntimeError("cyto mode needs the multimodal boundary stain image; none found")
        b = PyramidImage(ch["boundary"], level=0).read(r0, r1, c0, c1)
        stack = np.dstack([b, img] + ([np.zeros_like(img)] if major >= 4 else []))
        masks = s03a.run_cellpose(model, major, stack, diameter_um / px, opts)

    props = pd.DataFrame(regionprops_table(masks, properties=("label", "centroid", "area")))
    props["x_um"] = (c0 + props["centroid-1"] + 0.5) * px
    props["y_um"] = (r0 + props["centroid-0"] + 0.5) * px
    props["cell_area"] = props["area"] * px * px
    near = (props.x_um.between(a.x0 - margin, a.x1 + margin) & props.y_um.between(a.y0 - margin, a.y1 + margin))
    props = props[near].reset_index(drop=True)
    props["cell_id"] = ["cp-%d" % v for v in props["label"]]

    # transcripts of the padded crop -> mask label (qv filter as in 01 / 03a)
    f = ((ds.field("x_location") >= c0 * px) & (ds.field("x_location") < c1 * px) &
         (ds.field("y_location") >= r0 * px) & (ds.field("y_location") < r1 * px) &
         (ds.field("qv") >= a.min_qv))
    tx = ds.dataset(os.path.join(xenium_dir, "transcripts.parquet")).to_table(
        columns=["x_location", "y_location", "feature_name"], filter=f)
    names = tx.column("feature_name")
    if names.type == "binary" or str(names.type).startswith("large_binary"):
        names = pc.cast(names, "string")
    names = np.asarray(names.to_pylist(), dtype=object)
    cc = np.floor(tx.column("x_location").to_numpy() / px).astype(np.int64) - c0
    rr = np.floor(tx.column("y_location").to_numpy() / px).astype(np.int64) - r0
    ok = (cc >= 0) & (cc < masks.shape[1]) & (rr >= 0) & (rr < masks.shape[0])
    lab = np.zeros(len(cc), dtype=np.int64)
    lab[ok] = masks[rr[ok], cc[ok]]
    types = np.array([feature_type(n) for n in names], dtype=object)
    is_gene = types == "Gene Expression"
    is_ctrl = np.isin(types, list(QC_CONTROL_TYPES))
    t = pd.DataFrame({"label": lab, "gene": np.where(is_gene, names, None), "is_gene": is_gene, "is_ctrl": is_ctrl})
    t = t[t.label > 0]
    agg = t.groupby("label").agg(n_counts=("is_gene", "sum"), n_ctrl=("is_ctrl", "sum"),
                                 n_features=("gene", "nunique"))
    props = props.merge(agg, left_on="label", right_index=True, how="left").fillna(
        {"n_counts": 0, "n_ctrl": 0, "n_features": 0})
    props["control_frac"] = props.n_ctrl / np.maximum(props.n_counts + props.n_ctrl, 1)

    # 01's QC rules; later assignments win, as in qc_reason() in R/helpers.R
    reason = np.full(len(props), "pass", dtype=object)
    reason[props.control_frac > a.max_control_frac] = "high control"
    if a.max_area is not None:
        reason[props.cell_area > a.max_area] = "too large"
    reason[props.cell_area < a.min_area] = "too small"
    reason[(props.n_counts < a.min_counts) | (props.n_features < a.min_features)] = "low counts"
    props["qc_reason"] = reason
    props[["cell_id", "x_um", "y_um", "cell_area", "n_counts", "n_features", "control_frac", "qc_reason"]].to_csv(
        prefix + "_cellpose_cells.csv", index=False)

    ids = dict(zip(props["label"], props["cell_id"]))
    outlines(masks, props["label"], c0, r0, px, ids).to_csv(prefix + "_cellpose_cell_boundaries.csv.gz", index=False)
    if nuclei is not None:
        outlines(nuclei, props["label"], c0, r0, px, ids).to_csv(
            prefix + "_cellpose_nucleus_boundaries.csv.gz", index=False)
    return {"cellpose_mode": a.cellpose_mode, "cellpose_major": major, "cellpose_gpu": bool(gpu),
            "n_cellpose_cells": int(len(props)), "cellpose_seconds": round(time.time() - t0, 1)}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--xenium-dir", required=True)
    p.add_argument("--x0", type=float, required=True)
    p.add_argument("--y0", type=float, required=True)
    p.add_argument("--x1", type=float, required=True)
    p.add_argument("--y1", type=float, required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--image", default=None,
                   help="morphology image to crop (default: the DAPI file of morphology_focus/)")
    p.add_argument("--channel", type=int, default=0, help="channel index inside a multi-channel file")
    g = p.add_argument_group("cellpose on the window (same settings as 03a)")
    g.add_argument("--cellpose", action="store_true")
    g.add_argument("--cellpose-mode", choices=["nuclei_expand", "cyto"], default="nuclei_expand")
    g.add_argument("--cellpose-expand-um", type=float, default=5.0)
    g.add_argument("--cellpose-diameter-um", type=float, default=None)
    g.add_argument("--cellpose-gpu", type=lambda v: v.lower() in ("1", "true", "yes"), default=None,
                   help="default: use a GPU if PyTorch sees one")
    g.add_argument("--min-qv", type=float, default=20.0)
    g.add_argument("--min-counts", type=float, default=10)
    g.add_argument("--min-features", type=float, default=5)
    g.add_argument("--min-area", type=float, default=6)
    g.add_argument("--max-area", type=float, default=None, help="01's per-slide area ceiling")
    g.add_argument("--max-control-frac", type=float, default=0.05)
    a = p.parse_args()

    px = pixel_size_um(a.xenium_dir)
    path = a.image or find_morphology_channels(a.xenium_dir)["dapi"]
    img = PyramidImage(path, level=0, channel=a.channel)
    c0, r0 = int(np.floor(a.x0 / px)), int(np.floor(a.y0 / px))
    c1, r1 = int(np.ceil(a.x1 / px)), int(np.ceil(a.y1 / px))
    crop = img.read(r0, r1, c0, c1)
    if crop.ndim == 3:                       # several channels -> keep the requested one
        crop = crop[..., 0]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)

    import tifffile
    tifffile.imwrite(a.out_prefix + "_image.tif", crop.astype(np.uint16))
    meta = {"image": path, "pixel_size_um": px,
            "x0_um": c0 * px, "y0_um": r0 * px,
            "x1_um": (c0 + crop.shape[1]) * px, "y1_um": (r0 + crop.shape[0]) * px,
            "width_px": int(crop.shape[1]), "height_px": int(crop.shape[0])}
    for kind in ("cell", "nucleus"):
        b = read_boundaries(a.xenium_dir, kind, a.x0, a.y0, a.x1, a.y1)
        if b is not None:
            b.to_csv(a.out_prefix + f"_{kind}_boundaries.csv.gz", index=False)
            meta[f"n_{kind}_outlines"] = int(b.cell_id.nunique())
    if a.cellpose:
        try:
            meta.update(cellpose_window(a, a.xenium_dir, px, a.out_prefix))
        except Exception as e:                       # keep the crop usable without cellpose
            meta["cellpose_error"] = f"{type(e).__name__}: {e}"
            print("cellpose on the window failed:", meta["cellpose_error"], file=sys.stderr)
    with open(a.out_prefix + "_meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
