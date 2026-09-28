"""Shared I/O for the experiment1_tma Python scripts.

Everything here converts between Xenium outs, image pyramids and the one
"common segmentation format" the R side reads (R/helpers.R::load_segmentation):

    <out_dir>/cell_feature_matrix/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz}
    <out_dir>/cells.csv.gz          cell_id, core_id, x_centroid, y_centroid, cell_area, nucleus_area
    <out_dir>/assignment_summary.csv  per core: transcripts total / assigned (per method)

Coordinates are always microns in the Xenium frame (x = column, y = row,
y growing downwards), exactly as in cells.parquet / transcripts.parquet.
"""

from __future__ import annotations

import glob
import gzip
import json
import os
import re

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import scipy.io
import scipy.sparse as sp

CONTROL_PREFIXES = {
    "NegControlProbe_": "Negative Control Probe",
    "NegControlCodeword_": "Negative Control Codeword",
    "UnassignedCodeword_": "Unassigned Codeword",
    "BLANK_": "Blank Codeword",
    "DeprecatedCodeword_": "Deprecated Codeword",
    "GenomicControl": "Genomic Control",
    "Intergenic_Region": "Genomic Control",
}


# ---------------------------------------------------------------------------
# Xenium bundle metadata
# ---------------------------------------------------------------------------
def pixel_size_um(xenium_dir: str) -> float:
    """Morphology image pixel size (0.2125 um on current instruments)."""
    f = os.path.join(xenium_dir, "experiment.xenium")
    if os.path.exists(f):
        with open(f) as fh:
            meta = json.load(fh)
        if "pixel_size" in meta:
            return float(meta["pixel_size"])
    print("WARNING: pixel_size not found in experiment.xenium, assuming 0.2125 um")
    return 0.2125


def feature_type(name: str) -> str:
    for prefix, t in CONTROL_PREFIXES.items():
        if name.startswith(prefix):
            return t
    return "Gene Expression"


def read_features(xenium_dir: str) -> pd.DataFrame:
    """Feature list + types from the 10x matrix, so every segmentation uses the
    same feature order and control classes as Xenium's own output."""
    f = os.path.join(xenium_dir, "cell_feature_matrix", "features.tsv.gz")
    if os.path.exists(f):
        feats = pd.read_csv(f, sep="\t", header=None, names=["id", "name", "type"])
    else:
        # fall back to scanning the transcripts table for names
        tx = ds.dataset(os.path.join(xenium_dir, "transcripts.parquet"))
        names = pc.unique(_as_str(tx.to_table(columns=["feature_name"])["feature_name"])).to_pylist()
        feats = pd.DataFrame({"id": names, "name": names})
        feats["type"] = feats["name"].map(feature_type)
    return feats.reset_index(drop=True)


def _as_str(arr):
    """Older XOA versions store strings as binary."""
    t = arr.type
    if pa.types.is_binary(t) or pa.types.is_large_binary(t):
        return pc.cast(arr, pa.string())
    return arr


# ---------------------------------------------------------------------------
# Morphology images
# ---------------------------------------------------------------------------
def find_morphology_channels(xenium_dir: str) -> dict:
    """Returns {'dapi': path, 'boundary': path or None, 'interior_rna': ..., 'interior_protein': ...}.

    XOA >= 2 writes morphology_focus/ with one OME-TIFF per channel
    (names differ between versions: ch0000_dapi.ome.tif or
    morphology_focus_0000.ome.tif). Order is DAPI, boundary
    (ATP1A1/CD45/E-cadherin), interior RNA (18S), interior protein
    (aSMA/vimentin) for the multimodal segmentation kit. XOA 1.x has a single
    DAPI morphology_focus.ome.tif.
    """
    files = sorted(glob.glob(os.path.join(xenium_dir, "morphology_focus", "*.ome.tif")))
    if not files:
        single = os.path.join(xenium_dir, "morphology_focus.ome.tif")
        if os.path.exists(single):
            return {"dapi": single, "boundary": None, "interior_rna": None, "interior_protein": None}
        raise FileNotFoundError(f"No morphology_focus images in {xenium_dir}")
    out = {"dapi": files[0], "boundary": None, "interior_rna": None, "interior_protein": None}
    for f in files:
        b = os.path.basename(f).lower()
        if "dapi" in b:
            out["dapi"] = f
    for key, idx in (("boundary", 1), ("interior_rna", 2), ("interior_protein", 3)):
        if len(files) > idx:
            out[key] = files[idx]
    return out


class PyramidImage:
    """Lazy, region-readable view of one channel of an OME-TIFF pyramid.

    Reads only the requested window through tifffile's zarr store, so a
    whole-slide morphology image never has to fit in memory.
    """

    def __init__(self, path: str, level: int = 0, channel: int = 0):
        import tifffile
        import zarr

        self.path = path
        with tifffile.TiffFile(path) as tif:
            series = tif.series[0]
            self.n_levels = len(series.levels)
            lvl = series.levels[min(level, self.n_levels - 1)]
            self.axes = lvl.axes
            self.shape_full = series.levels[0].shape
            self.shape = lvl.shape
            self.ome_metadata = tif.ome_metadata
        store = tifffile.imread(path, aszarr=True, level=min(level, self.n_levels - 1))
        z = zarr.open(store, mode="r")
        self._z = z["0"] if hasattr(z, "keys") and "0" in list(z.keys()) else z
        self.channel = channel
        yx = [self.axes.index("Y"), self.axes.index("X")]
        full_yx = [self.shape_full[yx[0]], self.shape_full[yx[1]]]
        self.downsample = full_yx[0] / self.shape[yx[0]]
        self.height, self.width = self.shape[yx[0]], self.shape[yx[1]]

    def read(self, y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
        y0, x0 = max(0, int(y0)), max(0, int(x0))
        y1, x1 = min(self.height, int(y1)), min(self.width, int(x1))
        idx = []
        for ax in self.axes:
            if ax == "Y":
                idx.append(slice(y0, y1))
            elif ax == "X":
                idx.append(slice(x0, x1))
            elif ax == "C":
                idx.append(self.channel)
            elif ax in ("S",):          # RGB samples: keep all
                idx.append(slice(None))
            else:                        # Z, T, ... take first plane
                idx.append(0)
        return np.asarray(self._z[tuple(idx)])

    def read_all(self) -> np.ndarray:
        return self.read(0, self.height, 0, self.width)

    def channel_names(self) -> list:
        if not self.ome_metadata:
            return []
        return re.findall(r'<Channel[^>]*Name="([^"]*)"', self.ome_metadata)


# ---------------------------------------------------------------------------
# Transcripts
# ---------------------------------------------------------------------------
class Transcripts:
    """High-quality transcripts of one slide held as compact numpy arrays,
    sorted by x so per-region window queries are a binary search.

    Attributes: x, y (float32 um), feat (int32 index into `features`, -1 =
    not in the feature list), xen_assigned (bool: assigned to a cell by the
    10x segmentation), is_gene (bool).
    """

    def __init__(self, xenium_dir: str, features: pd.DataFrame, min_qv: float = 20.0):
        path = os.path.join(xenium_dir, "transcripts.parquet")
        dataset = ds.dataset(path)
        cols = dataset.schema.names
        want = ["x_location", "y_location", "feature_name", "qv", "cell_id"]
        missing = [c for c in want if c not in cols]
        if missing:
            raise ValueError(f"transcripts.parquet lacks columns {missing}")

        name_to_idx = {n: i for i, n in enumerate(features["name"])}
        gene_mask = (features["type"] == "Gene Expression").to_numpy()
        xs, ys, fs, xa = [], [], [], []
        for batch in dataset.to_batches(columns=want, filter=ds.field("qv") >= min_qv,
                                        batch_size=5_000_000):
            names = _as_str(batch.column("feature_name")).to_numpy(zero_copy_only=False)
            f = pd.Series(names).map(name_to_idx).fillna(-1).to_numpy(np.int32)
            cid = _as_str(batch.column("cell_id"))
            if pa.types.is_integer(cid.type):      # XOA 1.0-1.2: -1 = unassigned
                assigned = pc.greater_equal(cid, 0)
            else:                                   # later: "UNASSIGNED"
                assigned = pc.and_(pc.not_equal(cid, "UNASSIGNED"), pc.not_equal(cid, ""))
            assigned = pc.fill_null(assigned, False).to_numpy(zero_copy_only=False)
            xs.append(batch.column("x_location").to_numpy().astype(np.float32))
            ys.append(batch.column("y_location").to_numpy().astype(np.float32))
            fs.append(f)
            xa.append(np.asarray(assigned, dtype=bool))
        x = np.concatenate(xs); y = np.concatenate(ys)
        feat = np.concatenate(fs); xen = np.concatenate(xa)
        order = np.argsort(x, kind="stable")
        self.x, self.y, self.feat, self.xen_assigned = x[order], y[order], feat[order], xen[order]
        self.is_gene = np.zeros(len(self.feat), dtype=bool)
        valid = self.feat >= 0
        self.is_gene[valid] = gene_mask[self.feat[valid]]
        print(f"Loaded {len(self.x):,} transcripts with qv >= {min_qv}")

    def window(self, xmin, xmax, ymin, ymax) -> np.ndarray:
        """Indices of transcripts inside the window (microns)."""
        lo, hi = np.searchsorted(self.x, [xmin, xmax])
        idx = np.arange(lo, hi)
        yy = self.y[idx]
        return idx[(yy >= ymin) & (yy < ymax)]


# ---------------------------------------------------------------------------
# Regions (TMA cores or tiles)
# ---------------------------------------------------------------------------
def regions_from_cores(cores_csv: str, slide_id: str, only=None) -> pd.DataFrame:
    """Core bounding boxes from 02_tma_dearray.R (tables/02_cores_all.csv)."""
    c = pd.read_csv(cores_csv)
    c = c[(c["slide_id"] == slide_id) & c["core_id"].notna()].copy()
    if only:
        c = c[c["core_id"].isin(only)]
    if c.empty:
        raise ValueError(f"No cores for slide {slide_id} in {cores_csv}")
    return c[["core_id", "xmin", "xmax", "ymin", "ymax"]].rename(columns={"core_id": "region_id"})


def regions_from_tiles(width_um: float, height_um: float, tile_um: float = 1000.0) -> pd.DataFrame:
    """Whole-slide fallback when there is no core table."""
    rows = []
    for yi, y0 in enumerate(np.arange(0, height_um, tile_um)):
        for xi, x0 in enumerate(np.arange(0, width_um, tile_um)):
            rows.append({"region_id": f"tile_{yi}_{xi}", "xmin": x0, "xmax": x0 + tile_um,
                         "ymin": y0, "ymax": y0 + tile_um})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Writer for the common format
# ---------------------------------------------------------------------------
def write_common_format(out_dir: str, counts: sp.csr_matrix, cell_ids, features: pd.DataFrame,
                        cells: pd.DataFrame, summary: pd.DataFrame | None = None):
    """counts: cells x features (csr). Written as features x cells like 10x."""
    mdir = os.path.join(out_dir, "cell_feature_matrix")
    os.makedirs(mdir, exist_ok=True)
    with gzip.open(os.path.join(mdir, "matrix.mtx.gz"), "wb") as fh:
        scipy.io.mmwrite(fh, sp.csc_matrix(counts.T).astype(np.int32), field="integer")
    features[["id", "name", "type"]].to_csv(os.path.join(mdir, "features.tsv.gz"), sep="\t",
                                            header=False, index=False, compression="gzip")
    pd.Series(cell_ids).to_csv(os.path.join(mdir, "barcodes.tsv.gz"), header=False, index=False,
                               compression="gzip")
    cells.to_csv(os.path.join(out_dir, "cells.csv.gz"), index=False, compression="gzip")
    if summary is not None:
        summary.to_csv(os.path.join(out_dir, "assignment_summary.csv"), index=False)
    print(f"Wrote {counts.shape[0]:,} cells x {counts.shape[1]} features -> {out_dir}")
