#!/usr/bin/env python
"""06a_cellscape_segment.py -- only needed if you have CellScape IMAGES but no
per-cell export (or want the same segmentation logic on both modalities).

If CellScape's own software already exported a per-cell table (mean
intensity per marker + x/y), skip this and point `cellscape_cells` in
config/slides.csv at that export instead.

Segments cells with cellpose on DAPI (+ optionally a summed membrane channel)
tile by tile over the whole slide, then measures the mean intensity of every
channel per cell. Output columns match config.R defaults:
  CellID, X, Y (pixels), Area (pixels), <marker>_Mean ...

Input: either one multi-channel OME-TIFF (channel names read from OME
metadata) or a folder of single-channel TIFFs (file name = marker name).

Example:
  python 06a_cellscape_segment.py --image data/cellscape/CS_TMA_ORGAN/image.ome.tif \
      --dapi DAPI --membrane CD45 PanCK --out data/cellscape/CS_TMA_ORGAN/cells.csv --gpu
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "py"))
from xenium_io import PyramidImage  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--image", required=True, help="multi-channel OME-TIFF or folder of single-channel TIFFs")
    p.add_argument("--out", required=True)
    p.add_argument("--dapi", default="DAPI", help="nuclear channel name (or index)")
    p.add_argument("--membrane", nargs="*", default=[], help="channels summed into a membrane image")
    p.add_argument("--pixel-size-um", type=float, default=None, help="only used to set the default diameter")
    p.add_argument("--diameter-px", type=float, default=None)
    p.add_argument("--expand-px", type=int, default=0, help="expand nuclei by N px when no membrane channel")
    p.add_argument("--tile-px", type=int, default=4096)
    p.add_argument("--overlap-px", type=int, default=128)
    p.add_argument("--min-dapi", type=float, default=None, help="skip tiles whose DAPI max is below this")
    p.add_argument("--gpu", action="store_true")
    return p.parse_args()


def open_channels(path):
    """Returns {marker_name: PyramidImage}."""
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, "*.tif*")))
        return {os.path.basename(f).split(".")[0]: PyramidImage(f) for f in files}
    probe = PyramidImage(path, channel=0)
    names = probe.channel_names()
    n_ch = probe.shape[probe.axes.index("C")] if "C" in probe.axes else 1
    if len(names) != n_ch:
        names = [f"ch{i}" for i in range(n_ch)]
        print("WARNING: no channel names in OME metadata; using ch0, ch1, ... -- pass indices via --dapi")
    return {nm: PyramidImage(path, channel=i) for i, nm in enumerate(names)}


def main():
    a = parse_args()
    from importlib.metadata import version
    from cellpose import models
    from skimage.measure import regionprops_table
    from skimage.segmentation import expand_labels

    ch = open_channels(a.image)
    print("Channels:", list(ch))
    names = list(ch)
    dapi_name = names[int(a.dapi)] if a.dapi.isdigit() else a.dapi
    if dapi_name not in ch:
        sys.exit(f"DAPI channel '{dapi_name}' not found in {names}")
    mem_names = [names[int(m)] if m.isdigit() else m for m in a.membrane]
    H, W = ch[dapi_name].height, ch[dapi_name].width

    major = int(version("cellpose").split(".")[0])
    if major >= 4:
        model = models.CellposeModel(gpu=a.gpu)
    else:
        model = models.CellposeModel(gpu=a.gpu, model_type="cyto3" if mem_names else "nuclei")
    diameter = a.diameter_px or ((10.0 / a.pixel_size_um) if a.pixel_size_um else None)

    step = a.tile_px - a.overlap_px
    rows, next_id = [], 1
    for y0 in range(0, H, step):
        for x0 in range(0, W, step):
            y1, x1 = min(H, y0 + a.tile_px), min(W, x0 + a.tile_px)
            dapi = ch[dapi_name].read(y0, y1, x0, x1).astype(np.float32)
            if a.min_dapi is not None and dapi.max() < a.min_dapi:
                continue
            if mem_names:
                mem = sum(ch[m].read(y0, y1, x0, x1).astype(np.float32) for m in mem_names)
                img = np.dstack([mem, dapi]) if major < 4 else np.dstack([mem, dapi, np.zeros_like(dapi)])
                kw = {"channels": [1, 2]} if major < 4 else {}
            else:
                img = dapi
                kw = {"channels": [0, 0]} if major < 4 else {}
            if major < 4 and diameter:
                kw["diameter"] = diameter
            masks = np.asarray(model.eval(img, **kw)[0]).astype(np.int32)
            if not mem_names and a.expand_px > 0:
                masks = expand_labels(masks, distance=a.expand_px)
            if masks.max() == 0:
                continue
            stack = np.dstack([ch[n].read(y0, y1, x0, x1) for n in names]).astype(np.float32)
            props = pd.DataFrame(regionprops_table(masks, intensity_image=stack,
                                                   properties=("label", "centroid", "area", "mean_intensity")))
            props["X"] = x0 + props["centroid-1"]
            props["Y"] = y0 + props["centroid-0"]
            # keep cells in the tile interior (half the overlap on each shared side)
            half = a.overlap_px / 2
            lo_x = x0 + (half if x0 > 0 else 0); hi_x = x1 - (half if x1 < W else 0)
            lo_y = y0 + (half if y0 > 0 else 0); hi_y = y1 - (half if y1 < H else 0)
            props = props[(props.X >= lo_x) & (props.X < hi_x) & (props.Y >= lo_y) & (props.Y < hi_y)]
            out = pd.DataFrame({"CellID": np.arange(next_id, next_id + len(props)),
                                "X": props["X"].to_numpy(), "Y": props["Y"].to_numpy(),
                                "Area": props["area"].to_numpy()})
            for i, n in enumerate(names):
                out[f"{n}_Mean"] = props[f"mean_intensity-{i}"].to_numpy()
            next_id += len(props)
            rows.append(out)
            print(f"tile y={y0} x={x0}: {len(out)} cells")

    res = pd.concat(rows, ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    res.to_csv(a.out, index=False)
    print(f"Wrote {len(res):,} cells -> {a.out}")


if __name__ == "__main__":
    main()
