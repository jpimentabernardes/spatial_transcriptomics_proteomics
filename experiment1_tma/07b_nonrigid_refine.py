#!/usr/bin/env python
"""07b_nonrigid_refine.py -- optional non-rigid refinement (GEASO or Spateo),
per TMA core, on top of the rigid/affine result of 07_register_modalities.py.

When to use: only if the 07 QC (registration_qc_cellscape_<slide>.csv and the
overlay figures) shows cores whose SHAPE differs between sections (folds,
stretching, tears) -- a rigid transform cannot fix those. For intact cores the
affine result is usually as good, and non-rigid methods on serial sections
can "over-fit" (warp tissue to make labels match). Every result is therefore
checked and rejected if it moves cells implausibly far or does not improve
cross-modal agreement.

Cross-modality problem: GEASO and Spateo align slices by EXPRESSION
similarity, but RNA and protein measure different features. We give both
methods a SHARED, modality-agnostic representation per cell:
  * lineage composition of the cell's k nearest neighbours (same lineage
    names on both sides -- enforced in 06), and
  * z-scored, neighbourhood-smoothed values of the paired RNA/protein markers
    in config/rna_protein_pairs.csv (e.g. PTPRC <-> CD45).

Methods (APIs checked against the source repositories, 2026):
  geaso   GEASO.alignment.coarse_to_fine_alignment(source, target, rep_layer="latent",
          rep_field="obsm", ...) -> source.obsm["aligned_spatial_nonrigid"]
          Install: git clone https://github.com/xkmaxidian/GEASO (python 3.9 env,
          torch, dgl) and pass --geaso-repo.
  spateo  spateo.align.morpho_align(models=[target, source], rep_layer="latent",
          rep_field="obsm", mode="SN-N") -> source.obsm["align_spatial"]
          Install: pip install spateo-release

Example:
  python 07b_nonrigid_refine.py --project-dir /work/Spatial_TMA --slide-id TMA_ORGAN \
      --method spateo --cores ORG_A1 ORG_B2
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--method", choices=["geaso", "spateo"], default="spateo")
    p.add_argument("--geaso-repo", default=os.path.expanduser("~/GEASO"))
    p.add_argument("--pairs-csv", default=None, help="default: <code>/config/rna_protein_pairs.csv")
    p.add_argument("--cores", nargs="*", help="only these cores (default: all)")
    p.add_argument("--k", type=int, default=15, help="neighbours for the shared representation")
    p.add_argument("--feature-weight", type=float, default=0.5)
    p.add_argument("--max-median-shift-um", type=float, default=20.0)
    p.add_argument("--max-p95-shift-um", type=float, default=60.0)
    p.add_argument("--device", default="cpu", help="cpu or cuda / GPU id")
    return p.parse_args()


def shared_representation(xy, lineage, feats, lineages, k):
    """kNN lineage composition + kNN-smoothed z-scored paired features."""
    k = min(k, len(xy) - 1)
    _, nn = cKDTree(xy).query(xy, k=k + 1)
    lab_idx = np.searchsorted(lineages, lineage)
    onehot = np.zeros((len(xy), len(lineages)), np.float32)
    onehot[np.arange(len(xy)), lab_idx] = 1
    comp = onehot[nn].mean(1)
    if feats.shape[1] == 0:
        return comp
    z = (feats - feats.mean(0)) / (feats.std(0) + 1e-6)
    return comp, z[nn].mean(1).astype(np.float32)


def run_geaso(src_xy, tgt_xy, src_lat, tgt_lat, a):
    sys.path.insert(0, a.geaso_repo)
    import anndata as ad
    import torch
    from GEASO.alignment import coarse_to_fine_alignment
    src = ad.AnnData(X=src_lat.copy()); src.obsm["spatial"] = src_xy.copy(); src.obsm["latent"] = src_lat
    tgt = ad.AnnData(X=tgt_lat.copy()); tgt.obsm["spatial"] = tgt_xy.copy(); tgt.obsm["latent"] = tgt_lat
    dev = torch.device("cuda" if a.device != "cpu" and torch.cuda.is_available() else "cpu")
    n = min(len(src_xy), len(tgt_xy))
    coarse_to_fine_alignment(
        source=src, target=tgt, spatial_key="spatial", rep_layer="latent", rep_field="obsm",
        source_voxel_num=min(3000, n), target_voxel_num=min(3000, n), top_K=10, dis_metric="cos",
        use_latent=True, max_iter=500, non_rigid_start_iter=80, K=100, graph_knn=10, tau=0.5,
        lambda_VF=1e1, device=dev, verbose=False, inplace=True,
    )
    return np.asarray(src.obsm["aligned_spatial_nonrigid"])


def run_spateo(src_xy, tgt_xy, src_lat, tgt_lat, a):
    import anndata as ad
    import spateo as st
    src = ad.AnnData(X=src_lat.copy()); src.obsm["spatial"] = src_xy.copy(); src.obsm["latent"] = src_lat
    tgt = ad.AnnData(X=tgt_lat.copy()); tgt.obsm["spatial"] = tgt_xy.copy(); tgt.obsm["latent"] = tgt_lat
    aligned, _ = st.align.morpho_align(
        models=[tgt, src], rep_layer="latent", rep_field="obsm", spatial_key="spatial",
        key_added="align_spatial", mode="SN-N", dissimilarity="euclidean", max_iter=200,
        device=a.device, verbose=False,
    )
    return np.asarray(aligned[1].obsm["align_spatial"])


def main():
    a = parse_args()
    code_dir = os.path.dirname(os.path.abspath(__file__))
    reg_dir = os.path.join(a.project_dir, "data", "registration")
    from importlib import import_module
    lineage_agreement = import_module("07_register_modalities").lineage_agreement

    pairs = pd.read_csv(a.pairs_csv or os.path.join(code_dir, "config", "rna_protein_pairs.csv"))
    fix = pd.read_csv(os.path.join(reg_dir, f"xenium_cells_{a.slide_id}.csv.gz"))
    mov = pd.read_csv(os.path.join(reg_dir, f"cellscape_aligned_{a.slide_id}.csv.gz"))
    fix["lineage"] = fix["lineage"].astype(str); mov["lineage"] = mov["lineage"].astype(str)
    # R/Seurat turned "_" into "-" in feature names, and the exported columns follow that
    pairs["gene"] = pairs.gene.str.replace("_", "-")
    pairs["protein"] = pairs.protein.str.replace("_", "-")
    pairs = pairs[pairs.gene.isin(fix.columns) & pairs.protein.isin(mov.columns)]
    print(f"{len(pairs)} RNA/protein pairs available for the shared representation")
    lineages = np.array(sorted(set(fix.lineage) | set(mov.lineage)))

    mov["x_final"], mov["y_final"] = mov["x_aligned"], mov["y_aligned"]
    runner = run_geaso if a.method == "geaso" else run_spateo
    cores = a.cores or sorted(set(mov.core_id.dropna()) & set(fix.core_id.dropna()))
    rows = []
    for cid in cores:
        mm, ff = (mov.core_id == cid).to_numpy(), (fix.core_id == cid).to_numpy()
        if mm.sum() < 100 or ff.sum() < 100:
            continue
        mxy = mov.loc[mm, ["x_aligned", "y_aligned"]].to_numpy()
        fxy = fix.loc[ff, ["x_um", "y_um"]].to_numpy()
        mlab, flab = mov.loc[mm, "lineage"].to_numpy(), fix.loc[ff, "lineage"].to_numpy()
        mrep = shared_representation(mxy, mlab, mov.loc[mm, pairs.protein].to_numpy(float), lineages, a.k)
        frep = shared_representation(fxy, flab, fix.loc[ff, pairs.gene].to_numpy(float), lineages, a.k)
        if isinstance(mrep, tuple):
            mlat = np.hstack([mrep[0], a.feature_weight * mrep[1]]).astype(np.float32)
            flat = np.hstack([frep[0], a.feature_weight * frep[1]]).astype(np.float32)
        else:
            mlat, flat = mrep, frep
        try:
            new_xy = runner(mxy, fxy, mlat, flat, a)
        except Exception as e:
            print(f"{cid}: {a.method} failed ({e}) -- keeping affine result")
            rows.append(dict(core_id=cid, accepted=False, reason=f"error: {e}"))
            continue
        shift = np.linalg.norm(new_xy - mxy, axis=1)
        _, lift_before, _ = lineage_agreement(mxy, mlab, fxy, flab)
        _, lift_after, _ = lineage_agreement(new_xy, mlab, fxy, flab)
        ok = (np.median(shift) <= a.max_median_shift_um and np.quantile(shift, 0.95) <= a.max_p95_shift_um
              and not np.isnan(lift_after) and (np.isnan(lift_before) or lift_after >= lift_before))
        reason = "ok" if ok else "rejected: implausible shift or no agreement gain"
        if ok:
            mov.loc[mm, ["x_final", "y_final"]] = new_xy
        rows.append(dict(core_id=cid, accepted=ok, reason=reason, median_shift_um=float(np.median(shift)),
                         p95_shift_um=float(np.quantile(shift, 0.95)), lift_before=lift_before, lift_after=lift_after))
        print(f"{cid}: median shift {np.median(shift):.1f} um, lift {lift_before:.2f} -> {lift_after:.2f} [{reason}]")

    mov["nonrigid_method"] = np.where(mov.core_id.isin([r["core_id"] for r in rows if r.get("accepted")]),
                                      a.method, "none")
    mov.to_csv(os.path.join(reg_dir, f"cellscape_aligned_nonrigid_{a.slide_id}.csv.gz"), index=False)
    pd.DataFrame(rows).to_csv(os.path.join(reg_dir, f"registration_qc_nonrigid_{a.slide_id}.csv"), index=False)
    print("Done. 08_integrate_rna_protein.R uses x_final/y_final from this file when it exists.")


if __name__ == "__main__":
    main()
