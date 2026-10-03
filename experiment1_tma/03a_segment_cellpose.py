#!/usr/bin/env python
"""03a_segment_cellpose.py -- cellpose re-segmentation of Xenium TMA cores.

Runs cellpose on the Xenium morphology image one TMA core at a time (cores
from 02_tma_dearray.R), re-assigns the high-quality transcripts to the new
cells and writes the common segmentation format that
R/helpers.R::load_segmentation() reads (see py/xenium_io.py).

Two modes:
  nuclei_expand  DAPI-only nuclei, then expanded by --expand-um (10x's own
                 fallback is 5 um). Comparable to 10x nucleus-expansion cells.
  cyto           DAPI + the multimodal boundary stain (ATP1A1/CD45/E-cadherin),
                 whole-cell model. Only meaningful if the run used the
                 multimodal cell segmentation kit.

Per core, cells are kept only if their centroid lies inside the core box; the
crop is padded by --pad-um so cells at the edge are segmented whole. The same
logic tiles whole slides if you have no core table (--tile-um).

It also writes assignment_summary.csv with, for the SAME transcripts, how many
the 10x segmentation and cellpose assigned to cells -- one of the metrics in
03c_segmentation_benchmark.R.

ALTERNATIVE: 10x Xenium Ranger can rebuild a complete outs bundle from
cellpose masks (`xeniumranger import-segmentation --nuclei/--cells <mask>`),
which then opens in Xenium Explorer. Use that for visual review; this script
is the lightweight path that does not need Xenium Ranger.

Example (GPU strongly recommended; CPU works but is slow on whole cores):
  python 03a_segment_cellpose.py --slide-id TMA_ORGAN \
      --xenium-dir /path/to/output-XETG... \
      --cores-csv /path/to/tables/02_cores_all.csv \
      --out-dir /path/to/data/segmentation/cellpose/TMA_ORGAN \
      --mode nuclei_expand --gpu
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "py"))
from xenium_io import (PyramidImage, Transcripts, find_morphology_channels, pixel_size_um,  # noqa: E402
                       read_features, regions_from_cores, regions_from_tiles, write_common_format)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--xenium-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cores-csv", help="tables/02_cores_all.csv; omit to tile the whole slide")
    p.add_argument("--cores", nargs="*", help="only these core_ids (quick tests / benchmark subset)")
    p.add_argument("--mode", choices=["nuclei_expand", "cyto"], default="nuclei_expand")
    p.add_argument("--diameter-um", type=float, default=None,
                   help="expected object diameter; default 8 (nuclei) or 15 (cyto). Ignored by cellpose>=4 "
                        "unless --force-diameter")
    p.add_argument("--force-diameter", action="store_true")
    p.add_argument("--expand-um", type=float, default=5.0, help="nuclear expansion (nuclei_expand mode)")
    p.add_argument("--pad-um", type=float, default=30.0)
    p.add_argument("--tile-um", type=float, default=1000.0)
    p.add_argument("--min-qv", type=float, default=20.0)
    p.add_argument("--flow-threshold", type=float, default=0.4)
    p.add_argument("--cellprob-threshold", type=float, default=0.0)
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--save-masks", action="store_true", help="save per-core label masks (.npz)")
    p.add_argument("--qc-png", type=int, default=4, help="write overlay PNGs for the first N regions")
    p.add_argument("--level", type=int, default=0,
                   help="pyramid level of the morphology image: 0 = full resolution (0.2125 um/px), "
                        "1 = 2x downsampled (0.425 um/px, ~4x faster, nuclei still ~19 px). "
                        "Use 1 when running on CPU")
    return p.parse_args()


def load_cellpose(mode, gpu):
    from importlib.metadata import version
    from cellpose import models
    major = int(version("cellpose").split(".")[0])
    if major >= 4:
        # Cellpose-SAM: one generalist model, no model_type / channels
        model = models.CellposeModel(gpu=gpu)
    else:
        model = models.CellposeModel(gpu=gpu, model_type="nuclei" if mode == "nuclei_expand" else "cyto3")
    try:
        import torch
        on_gpu = bool(gpu) and torch.cuda.is_available()
        dev = torch.cuda.get_device_name(0) if on_gpu else f"CPU ({torch.get_num_threads()} threads)"
    except ImportError:
        on_gpu, dev = False, "CPU"
    print(f"cellpose {version('cellpose')} loaded -- running on {dev}")
    if gpu and not on_gpu:
        print("WARNING: gpu=True but PyTorch sees no GPU (job without --gpus-per-node=1, or a CPU-only torch).")
    if not on_gpu and major >= 4:
        print("NOTE: Cellpose-SAM (cellpose 4) on CPU takes many minutes per core. Use a GPU job, "
              "or level=1, or test a few cores first (cores=[...]).")
    return model, major


def run_cellpose(model, major, img, diameter_px, args):
    kw = dict(flow_threshold=args.flow_threshold, cellprob_threshold=args.cellprob_threshold)
    if major < 4:
        kw["diameter"] = diameter_px
        kw["channels"] = [0, 0] if img.ndim == 2 else [1, 2]   # [cytoplasm, nucleus], 1-based
    elif args.force_diameter:
        kw["diameter"] = diameter_px
    out = model.eval(img, **kw)
    return np.asarray(out[0]).astype(np.int32)


def qc_overlay(path, dapi, masks, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries
    step = max(1, int(np.ceil(max(dapi.shape) / 2000)))   # keep PNGs small
    d = dapi[::step, ::step].astype(float)
    lo, hi = np.percentile(d, [1, 99.5])
    d = np.clip((d - lo) / max(hi - lo, 1), 0, 1)
    b = find_boundaries(masks[::step, ::step], mode="inner")
    rgb = np.dstack([d, d, d])
    rgb[b] = [1, 0.2, 0.2]
    plt.figure(figsize=(10, 10))
    plt.imshow(rgb)
    plt.title(title)
    plt.axis("off")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()
    from skimage.measure import regionprops_table
    from skimage.segmentation import expand_labels

    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    px0 = pixel_size_um(args.xenium_dir)
    ch = find_morphology_channels(args.xenium_dir)
    print("Morphology channels:", json.dumps(ch, indent=1))
    dapi = PyramidImage(ch["dapi"], level=args.level)
    px = px0 * dapi.downsample                     # um per pixel at the level that is segmented
    print(f"segmenting pyramid level {args.level}: {px:.4f} um/px ({dapi.width} x {dapi.height} px)")
    if px > 0.5:
        print("WARNING: coarser than ~0.5 um/px, nuclei get too small for reliable segmentation; use level 0 or 1.")
    boundary = None
    if args.mode == "cyto":
        if ch["boundary"] is None:
            sys.exit("cyto mode needs the multimodal boundary stain image; none found.")
        boundary = PyramidImage(ch["boundary"], level=args.level)

    diameter_um = args.diameter_um or (8.0 if args.mode == "nuclei_expand" else 15.0)
    diameter_px = diameter_um / px
    expand_px = int(round(args.expand_um / px))

    features = read_features(args.xenium_dir)
    tx = Transcripts(args.xenium_dir, features, min_qv=args.min_qv)

    if args.cores_csv:
        regions = regions_from_cores(args.cores_csv, args.slide_id, args.cores)
    else:
        regions = regions_from_tiles(dapi.width * px, dapi.height * px, args.tile_um)
    print(f"{len(regions)} regions to segment")

    model, major = load_cellpose(args.mode, args.gpu)
    mask_dir = os.path.join(args.out_dir, "masks")
    qc_dir = os.path.join(args.out_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    if args.save_masks:
        os.makedirs(mask_dir, exist_ok=True)

    coo_rows, coo_cols = [], []
    cell_tables, summaries = [], []
    n_cells_total = 0

    t_loop = time.time()
    for i, reg in enumerate(regions.itertuples(index=False)):
        # padded crop in pixels, clipped to the image
        c0 = max(0, int(np.floor((reg.xmin - args.pad_um) / px)))
        r0 = max(0, int(np.floor((reg.ymin - args.pad_um) / px)))
        c1 = min(dapi.width, int(np.ceil((reg.xmax + args.pad_um) / px)))
        r1 = min(dapi.height, int(np.ceil((reg.ymax + args.pad_um) / px)))
        if c1 <= c0 or r1 <= r0:
            continue
        img_dapi = dapi.read(r0, r1, c0, c1)
        if img_dapi.max() == 0:
            continue

        if args.mode == "nuclei_expand":
            nuclei = run_cellpose(model, major, img_dapi, diameter_px, args)
            masks = expand_labels(nuclei, distance=expand_px) if expand_px > 0 else nuclei
        else:
            img_b = boundary.read(r0, r1, c0, c1)
            stack = np.dstack([img_b, img_dapi])
            if major >= 4:
                stack = np.dstack([img_b, img_dapi, np.zeros_like(img_dapi)])
            masks = run_cellpose(model, major, stack, diameter_px, args)
            nuclei = None

        if masks.max() == 0:
            continue

        props = pd.DataFrame(regionprops_table(masks, properties=("label", "centroid", "area")))
        props["x_centroid"] = (c0 + props["centroid-1"] + 0.5) * px
        props["y_centroid"] = (r0 + props["centroid-0"] + 0.5) * px
        props["cell_area"] = props["area"] * px * px
        if nuclei is not None:
            nuc = pd.DataFrame(regionprops_table(nuclei, properties=("label", "area")))
            props = props.merge(nuc.rename(columns={"area": "nuc_px"}), on="label", how="left")
            props["nucleus_area"] = props["nuc_px"] * px * px
        else:
            props["nucleus_area"] = np.nan
        # keep cells whose centroid is in the (unpadded) region: no duplicates, no cut cells
        inside = ((props.x_centroid >= reg.xmin) & (props.x_centroid < reg.xmax) &
                  (props.y_centroid >= reg.ymin) & (props.y_centroid < reg.ymax))
        props = props[inside].reset_index(drop=True)
        if props.empty:
            continue

        label_to_row = np.full(masks.max() + 1, -1, dtype=np.int64)
        label_to_row[props["label"].to_numpy()] = np.arange(len(props)) + n_cells_total

        # transcripts in the padded crop -> mask label -> kept cell
        idx = tx.window(c0 * px, c1 * px, r0 * px, r1 * px)
        cc = np.floor(tx.x[idx] / px).astype(np.int64) - c0
        rr = np.floor(tx.y[idx] / px).astype(np.int64) - r0
        ok = (cc >= 0) & (cc < masks.shape[1]) & (rr >= 0) & (rr < masks.shape[0])
        idx, cc, rr = idx[ok], cc[ok], rr[ok]
        cell_row = label_to_row[masks[rr, cc]]
        take = (cell_row >= 0) & (tx.feat[idx] >= 0)
        coo_rows.append(cell_row[take])
        coo_cols.append(tx.feat[idx][take])

        # assignment summary on transcripts in the unpadded region, genes only
        interior = ((tx.x[idx] >= reg.xmin) & (tx.x[idx] < reg.xmax) &
                    (tx.y[idx] >= reg.ymin) & (tx.y[idx] < reg.ymax) & tx.is_gene[idx])
        summaries.append({
            "region_id": reg.region_id,
            "n_tx_genes": int(interior.sum()),
            "n_assigned_cellpose": int((interior & (cell_row >= 0)).sum()),
            "n_assigned_xenium": int((interior & tx.xen_assigned[idx]).sum()),
        })

        props["cell_id"] = [f"{reg.region_id}-{lab}" for lab in props["label"]]
        props["core_id"] = reg.region_id
        cell_tables.append(props[["cell_id", "core_id", "x_centroid", "y_centroid", "cell_area", "nucleus_area"]])
        n_cells_total += len(props)

        if args.save_masks:
            os.makedirs(mask_dir, exist_ok=True)       # also when save_masks was switched on after setup
            np.savez_compressed(os.path.join(mask_dir, f"{reg.region_id}.npz"),
                                masks=masks, offset_rc=np.array([r0, c0]), pixel_size=px)
        if i < args.qc_png:
            os.makedirs(qc_dir, exist_ok=True)
            qc_overlay(os.path.join(qc_dir, f"{reg.region_id}_overlay.png"), img_dapi, masks,
                       f"{args.slide_id} {reg.region_id}: {len(props)} cells ({args.mode})")
        print(f"[{i + 1}/{len(regions)}] {reg.region_id}: {len(props):,} cells, "
              f"{summaries[-1]['n_assigned_cellpose'] / max(summaries[-1]['n_tx_genes'], 1):.1%} "
              f"transcripts assigned (10x: "
              f"{summaries[-1]['n_assigned_xenium'] / max(summaries[-1]['n_tx_genes'], 1):.1%}) | "
              f"{(time.time() - t_loop) / (i + 1) / 60:.1f} min/core, "
              f"~{(time.time() - t_loop) / (i + 1) * (len(regions) - i - 1) / 60:.0f} min left", flush=True)

    if n_cells_total == 0:
        sys.exit("No cells segmented -- check the image channel and diameter.")

    rows = np.concatenate(coo_rows)
    cols = np.concatenate(coo_cols)
    counts = sp.coo_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)),
                           shape=(n_cells_total, len(features))).tocsr()
    counts.sum_duplicates()
    cells = pd.concat(cell_tables, ignore_index=True)
    cells["segmentation_method"] = f"cellpose_{args.mode}"
    summary = pd.DataFrame(summaries)
    summary.insert(0, "slide_id", args.slide_id)
    write_common_format(args.out_dir, counts, cells["cell_id"].to_numpy(), features, cells, summary)

    with open(os.path.join(args.out_dir, "run_params.json"), "w") as fh:
        json.dump({**vars(args), "pixel_size_um": px, "pixel_size_level0_um": px0, "cellpose_major": major,
                   "runtime_min": round((time.time() - t0) / 60, 1)}, fh, indent=1)
    print(f"Done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
