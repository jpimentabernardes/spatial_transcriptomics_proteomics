"""Synthetic 200-core Xenium-like TMA (10 x 20, rotated, with missing cores) for tests."""

import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from tma_toolkit.metadata import make_template

ORGANS = ["colon", "lung", "liver", "tonsil", "kidney", "skin", "stomach", "breast"]
TYPES = {"Epithelial": ["EPCAM", "KRT8", "CDH1"], "T_cell": ["PTPRC", "CD3E", "CD8A"],
         "B_cell": ["PTPRC", "MS4A1", "CD79A"], "Macrophage": ["PTPRC", "CD68", "CD163"],
         "Fibroblast": ["COL1A1", "DCN", "LUM"], "Endothelial": ["PECAM1", "VWF", "CDH5"]}
GENES = sorted({g for v in TYPES.values() for g in v} | {"GENE_A", "GENE_B", "MKI67"})
CONTROLS = [f"NegControlProbe_{i}" for i in range(4)]


def make(root: str, slide_id: str = "TMA2", n_rows: int = 10, n_cols: int = 20, cells_per_core: int = 250,
         rotation_deg: float = 2.0, missing_col: int = 12, n_missing_random: int = 5, seed: int = 0):
    rng = np.random.default_rng(seed)
    meta = make_template(slide_id, n_rows, n_cols)
    meta["patient_id"] = [f"P{(i // 2) + 1:03d}" for i in range(len(meta))]      # 2 cores per patient
    meta["organ"] = [ORGANS[(i // 2) % len(ORGANS)] for i in range(len(meta))]
    # blocks of 16 cores (8 organs x 2) alternate normal / tumour -> tissue type balanced across organs
    meta["tissue_type"] = np.where((np.arange(len(meta)) // 16) % 2 == 0, "normal", "tumor")
    meta["site"] = np.where(np.arange(len(meta)) % 2 == 0, "center", "margin")
    meta["age"] = [40 + (i // 2) % 40 for i in range(len(meta))]
    meta["sex"] = ["F" if (i // 2) % 2 else "M" for i in range(len(meta))]
    meta.loc[3, "organ"] = "Colon "                                               # a typo to be caught
    missing = set(meta.index[meta.core_col == missing_col]) | set(rng.choice(len(meta), n_missing_random, replace=False))

    t = np.radians(rotation_deg)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    pitch, radius = 1200.0, 400.0
    rows_xy, truth = [], []
    ri = {r: i for i, r in enumerate(sorted(meta.core_row.unique(), key=lambda s: (len(s), s)))}
    for i, m in meta.iterrows():
        if i in missing:
            continue
        c = np.array([1000 + (m.core_col - 1) * pitch, 1000 + ri[m.core_row] * pitch]) @ R.T
        ang = rng.uniform(0, 2 * np.pi, cells_per_core)
        rad = radius * np.sqrt(rng.uniform(0, 1, cells_per_core))
        xy = c + np.c_[rad * np.cos(ang), rad * np.sin(ang)]
        comp = np.ones(len(TYPES)); comp[ORGANS.index(m.organ.strip().lower()) % len(TYPES)] += 3
        ct = rng.choice(list(TYPES), cells_per_core, p=comp / comp.sum())
        rows_xy.append(xy); truth += [(m.core_id, t_) for t_ in ct]
    xy = np.vstack(rows_xy)
    truth = pd.DataFrame(truth, columns=["core_id", "cell_type"])

    gi = {g: j for j, g in enumerate(GENES)}
    rate = np.full((len(truth), len(GENES)), 0.05)
    for j, t_ in enumerate(truth.cell_type):
        for g in TYPES[t_]:
            rate[j, gi[g]] += 3
    rate[truth.core_id.str.len() > 0, gi["GENE_A"]] += np.where(
        truth.core_id.map(meta.set_index("core_id").tissue_type) == "tumor", 2.0, 0.0)   # tumour-up gene
    X = rng.poisson(rate)
    Cn = rng.poisson(0.02, (len(truth), len(CONTROLS)))
    M = sp.csr_matrix(np.hstack([X, Cn]))

    d = os.path.join(root, "xenium", slide_id)
    os.makedirs(os.path.join(d, "cell_feature_matrix"), exist_ok=True)
    feats = pd.DataFrame({"id": GENES + CONTROLS, "name": GENES + CONTROLS,
                          "type": ["Gene Expression"] * len(GENES) + ["Negative Control Probe"] * len(CONTROLS)})
    ids = [f"c{i}-1" for i in range(len(truth))]
    with gzip.open(os.path.join(d, "cell_feature_matrix", "matrix.mtx.gz"), "wb") as fh:
        scipy.io.mmwrite(fh, sp.csc_matrix(M.T), field="integer")
    feats.to_csv(os.path.join(d, "cell_feature_matrix", "features.tsv.gz"), sep="\t", header=False, index=False)
    pd.Series(ids).to_csv(os.path.join(d, "cell_feature_matrix", "barcodes.tsv.gz"), header=False, index=False)
    pd.DataFrame({"cell_id": ids, "x_centroid": xy[:, 0], "y_centroid": xy[:, 1],
                  "cell_area": rng.uniform(30, 120, len(ids)), "nucleus_area": rng.uniform(10, 40, len(ids))}
                 ).to_csv(os.path.join(d, "cells.csv.gz"), index=False)
    truth.index = [f"{slide_id}_{i}" for i in ids]
    pd.DataFrame({"slide_id": [slide_id], "xenium_dir": [f"xenium/{slide_id}"]}).to_csv(
        os.path.join(root, "slides.csv"), index=False)
    meta.to_csv(os.path.join(root, "core_metadata.csv"), index=False)
    truth.to_csv(os.path.join(root, "truth.csv"))
    return meta, truth, missing
