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
from xenium_io import PyramidImage, find_morphology_channels, pixel_size_um  # noqa: E402


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
    with open(a.out_prefix + "_meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
