#!/usr/bin/env python
"""03b_segger_to_common.py -- segger output -> common segmentation format.

Reads segger_segmentation.parquet (one row per transcript), keeps the
recommended assignments (`filtered` == True), and writes the 10x-style matrix
+ cells.csv.gz + assignment_summary.csv that R/helpers.R::load_segmentation()
reads, with core ids from 02_tma_dearray.R.

Cell centroid = mean position of the cell's transcripts. Cell area comes from
`segger export boundaries` (GeoParquet) if available, else it is left NA
(03c then skips the area-based metrics for segger instead of guessing).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "py"))
from xenium_io import Transcripts, read_features, write_common_format  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--xenium-dir", required=True)
    p.add_argument("--segger-parquet", required=True)
    p.add_argument("--boundaries", default=None)
    p.add_argument("--cores-csv", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-transcripts", type=int, default=1)
    p.add_argument("--min-qv", type=float, default=20.0)
    return p.parse_args()


def core_of(x, y, cores):
    """Core id for points by core bounding box (first match), NaN outside."""
    out = np.full(len(x), None, dtype=object)
    for c in cores.itertuples(index=False):
        m = (x >= c.xmin) & (x < c.xmax) & (y >= c.ymin) & (y < c.ymax) & (out == None)  # noqa: E711
        out[m] = c.core_id
    return out


def main():
    a = parse_args()
    features = read_features(a.xenium_dir)
    name_to_idx = {n: i for i, n in enumerate(features["name"])}
    gene_names = set(features.loc[features["type"] == "Gene Expression", "name"])

    seg = pd.read_parquet(a.segger_parquet,
                          columns=["segger_cell_id", "x", "y", "feature_name", "filtered"])
    seg["feature_name"] = seg["feature_name"].astype(str)
    print(f"{len(seg):,} transcripts in segger output")

    cores = pd.read_csv(a.cores_csv)
    cores = cores[(cores["slide_id"] == a.slide_id) & cores["core_id"].notna()]

    keep = seg["filtered"].fillna(False).astype(bool) & seg["segger_cell_id"].notna()
    asg = seg.loc[keep].copy()
    asg["cell"] = asg["segger_cell_id"].astype(str)
    asg["feat"] = asg["feature_name"].map(name_to_idx)
    asg = asg[asg["feat"].notna()]
    asg["feat"] = asg["feat"].astype(np.int64)

    cells = (asg.groupby("cell", sort=False)
                .agg(x_centroid=("x", "mean"), y_centroid=("y", "mean"), n_tx=("x", "size"))
                .reset_index())
    cells = cells[cells["n_tx"] >= a.min_transcripts].reset_index(drop=True)
    cells["core_id"] = core_of(cells["x_centroid"].to_numpy(), cells["y_centroid"].to_numpy(), cores)
    cells = cells[cells["core_id"].notna()].reset_index(drop=True)
    cell_index = pd.Series(np.arange(len(cells)), index=cells["cell"])

    asg = asg[asg["cell"].isin(cell_index.index)]
    rows = cell_index.loc[asg["cell"]].to_numpy()
    counts = sp.coo_matrix((np.ones(len(rows), dtype=np.int32), (rows, asg["feat"].to_numpy())),
                           shape=(len(cells), len(features))).tocsr()
    counts.sum_duplicates()

    cells["cell_area"] = np.nan
    cells["nucleus_area"] = np.nan
    if a.boundaries and os.path.exists(a.boundaries):
        try:
            import geopandas as gpd
            b = gpd.read_parquet(a.boundaries)
            id_col = next(c for c in b.columns if c != b.geometry.name)
            area = pd.Series(b.geometry.area.to_numpy(), index=b[id_col].astype(str))
            cells["cell_area"] = cells["cell"].map(area)
            print(f"cell areas from {a.boundaries} (id column '{id_col}')")
        except Exception as e:  # keep going without areas
            print(f"WARNING: could not read boundaries ({e}); cell_area left NA")

    # Assignment summary on the same transcript definition as 03a: qv-filtered
    # gene transcripts inside each core box. 10x assignment comes from the raw
    # transcripts table; segger's from its own output.
    tx = Transcripts(a.xenium_dir, features, min_qv=a.min_qv)
    seg_gene = seg[seg["feature_name"].isin(gene_names)]
    seg_core = core_of(seg_gene["x"].to_numpy(), seg_gene["y"].to_numpy(), cores)
    seg_ok = keep.loc[seg_gene.index].to_numpy()
    summaries = []
    for c in cores.itertuples(index=False):
        idx = tx.window(c.xmin, c.xmax, c.ymin, c.ymax)
        g = tx.is_gene[idx]
        in_core = seg_core == c.core_id
        summaries.append({
            "slide_id": a.slide_id, "region_id": c.core_id,
            "n_tx_genes": int(g.sum()),
            # segger may keep low-qv transcripts, so its rate uses its own denominator
            "n_tx_genes_segger": int(in_core.sum()),
            "n_assigned_segger": int((in_core & seg_ok).sum()),
            "n_assigned_xenium": int((g & tx.xen_assigned[idx]).sum()),
        })

    cells["cell_id"] = cells["cell"]
    cells["segmentation_method"] = "segger"
    write_common_format(a.out_dir, counts, cells["cell_id"].to_numpy(), features,
                        cells[["cell_id", "core_id", "x_centroid", "y_centroid", "cell_area",
                               "nucleus_area", "segmentation_method"]],
                        pd.DataFrame(summaries))


if __name__ == "__main__":
    main()
