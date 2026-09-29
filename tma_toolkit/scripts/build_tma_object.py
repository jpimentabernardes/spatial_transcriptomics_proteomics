#!/usr/bin/env python
"""Build one AnnData for a TMA project: all slides + core assignment + core metadata + QC.

    python scripts/build_tma_object.py \
        --slides config/slides.csv \
        --core-metadata config/core_metadata.csv \
        --project-dir /work/Spatial_TMA \
        --out /work/Spatial_TMA/data/processed/tma.h5ad

slides.csv          slide_id, xenium_dir [, row_flip, col_flip]  (paths relative to --project-dir)
core_metadata.csv   one row per core (make it with make_core_metadata_template.py)
Optional:
  --overrides        CSV slide_id, detected_core, core_row, core_col   (fix dearray labels)
  --polygons-dir     <dir>/<slide_id>/<core_id>.csv with X,Y columns (um), e.g. Xenium Explorer
                     "Download selection coordinates" -- used instead of automatic dearraying
  --segmentation     xenium | cellpose | segger (common 10x layout under data/segmentation/)

Outputs: the .h5ad, <figures>/dearray_<slide>.png (CHECK THESE against the TMA map),
<tables>/cores.csv, <tables>/core_qc_flags.csv.
Cell QC (flags stored in obs, failing cells removed unless --keep-failing):
  n_counts < --min-counts | n_genes < --min-genes | area outside [--min-area, --max-area]
  | negative-control fraction > --max-ctrl-frac
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tma_toolkit import TMACohort, attach_core_metadata, read_core_metadata  # noqa: E402
from tma_toolkit.dearray import assign_cells_to_polygons, dearray_slide  # noqa: E402
from tma_toolkit.io import load_slides  # noqa: E402
from tma_toolkit.plotting import plot_dearray  # noqa: E402


def truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slides", required=True)
    p.add_argument("--core-metadata", required=True)
    p.add_argument("--project-dir", default=".")
    p.add_argument("--out", required=True)
    p.add_argument("--segmentation", default="xenium")
    p.add_argument("--overrides")
    p.add_argument("--polygons-dir")
    p.add_argument("--figures", default=None, help="default <project-dir>/figures")
    p.add_argument("--tables", default=None, help="default <project-dir>/tables")
    p.add_argument("--eps", type=float, default=None, help="DBSCAN radius (um); default adaptive")
    p.add_argument("--min-cells-per-core", type=int, default=200)
    p.add_argument("--merge-dist", type=float, default=400.0)
    p.add_argument("--min-counts", type=int, default=10)
    p.add_argument("--min-genes", type=int, default=5)
    p.add_argument("--min-area", type=float, default=6.0)
    p.add_argument("--max-area", type=float, default=600.0)
    p.add_argument("--max-ctrl-frac", type=float, default=0.05)
    p.add_argument("--keep-failing", action="store_true", help="keep QC-failing cells (flagged only)")
    return p.parse_args()


def main():
    a = parse_args()
    fig_dir = a.figures or os.path.join(a.project_dir, "figures")
    tab_dir = a.tables or os.path.join(a.project_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True); os.makedirs(tab_dir, exist_ok=True)

    slides = pd.read_csv(a.slides)
    meta = read_core_metadata(a.core_metadata)
    unknown = set(meta.slide_id) - set(slides.slide_id.astype(str))
    if unknown:
        print(f"WARNING: core metadata mentions slides not in slides.csv: {sorted(unknown)}")
    overrides = pd.read_csv(a.overrides) if a.overrides else None

    adata = load_slides(slides, base_dir=a.project_dir, segmentation=a.segmentation)
    adata.obs["core_id"] = None

    geometry = []
    for s in slides.itertuples(index=False):
        sid = str(s.slide_id)
        poly_dir = os.path.join(a.polygons_dir, sid) if a.polygons_dir else None
        if poly_dir and os.path.isdir(poly_dir):
            polys = {}
            for f in glob.glob(os.path.join(poly_dir, "*.csv")):
                d = pd.read_csv(f, comment="#")
                cols = {c.lower(): c for c in d.columns}
                polys[os.path.splitext(os.path.basename(f))[0]] = d[[cols["x"], cols["y"]]].to_numpy(float)
            m = (adata.obs.slide_id == sid).to_numpy()
            col = adata.obs["core_id"].astype(object).to_numpy()
            col[m] = assign_cells_to_polygons(adata.obsm["spatial"][m], polys)
            adata.obs["core_id"] = col
            geometry.append(pd.DataFrame([{"core_id": k, "slide_id": sid, "x_center": v[:, 0].mean(),
                                           "y_center": v[:, 1].mean(),
                                           "radius": np.hypot(*(v - v.mean(0)).T).max()} for k, v in polys.items()]))
            print(f"{sid}: polygon mode, {len(polys)} cores")
        else:
            cores = dearray_slide(adata, meta, sid,
                                  row_flip=truthy(getattr(s, "row_flip", False)),
                                  col_flip=truthy(getattr(s, "col_flip", False)),
                                  overrides=overrides, min_counts=a.min_counts, eps=a.eps,
                                  min_cells=a.min_cells_per_core, merge_dist=a.merge_dist)
            fig = plot_dearray(adata, cores, sid)
            fig.savefig(os.path.join(fig_dir, f"dearray_{sid}.png"), dpi=150, bbox_inches="tight")
            print(f"{sid}: {cores.core_id.notna().sum()} cores labelled -> {fig_dir}/dearray_{sid}.png")
            geometry.append(cores)
    geometry = pd.concat(geometry, ignore_index=True)

    # ---- cell QC flags
    obs = adata.obs
    ctrl = obs[[c for c in obs if c.startswith("n_ctrl_")]].sum(1) if any(c.startswith("n_ctrl_") for c in obs) else 0
    obs["ctrl_frac"] = ctrl / np.maximum(obs["n_counts"] + ctrl, 1)
    obs["qc_low_counts"] = (obs.n_counts < a.min_counts) | (obs.n_genes < a.min_genes)
    obs["qc_area"] = obs.cell_area.notna() & ((obs.cell_area < a.min_area) | (obs.cell_area > a.max_area))
    obs["qc_high_ctrl"] = obs.ctrl_frac > a.max_ctrl_frac
    obs["qc_pass"] = ~(obs.qc_low_counts | obs.qc_area | obs.qc_high_ctrl)
    flags = obs.groupby(obs.core_id.fillna("(outside cores)"), observed=True)[
        ["qc_low_counts", "qc_area", "qc_high_ctrl", "qc_pass"]].mean().mul(100).round(1)
    flags.to_csv(os.path.join(tab_dir, "core_qc_flags.csv"))
    print(f"{obs.core_id.isna().mean():.1%} of cells outside cores; {obs.qc_pass.mean():.1%} pass cell QC")

    keep = obs.core_id.notna().to_numpy() & (a.keep_failing | obs.qc_pass.to_numpy())
    adata = adata[keep].copy()
    adata.obs["core_id"] = adata.obs["core_id"].astype(str)

    # ---- core table = metadata sheet + geometry + QC; metadata copied onto cells
    attach_core_metadata(adata, meta)
    geo = geometry.dropna(subset=["core_id"]).set_index("core_id")
    geo_cols = [c for c in ["detected_core", "x_center", "y_center", "radius", "grid_rotation_deg"] if c in geo]
    cores = meta.join(geo[geo_cols], how="left")
    cores["found_on_slide"] = cores.index.isin(geo.index)
    cohort = TMACohort(adata, cores=cores).refresh_qc()
    cohort.cores.to_csv(os.path.join(tab_dir, "cores.csv"))
    print(cohort)
    print(cohort.summary())
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    cohort.write(a.out)


if __name__ == "__main__":
    main()
