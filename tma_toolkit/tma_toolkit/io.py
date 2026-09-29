"""Load Xenium slides (or any segmentation in the same 10x layout) into AnnData.

One AnnData for all slides:
  * X / layers["counts"]  genes only (raw counts)
  * obs                   slide_id, cell_id, x_um, y_um, cell_area, nucleus_area, n_counts, n_genes,
                          n_<control class> (negative-control probe / codeword counts per cell)
  * obsm["spatial"]       x_um, y_um (microns, Xenium frame: y grows downwards)
  * uns["tma"]["n_control_features"]  number of features per control class and slide, needed to
                          turn control counts into a per-probe background rate
  * obs_names             "<slide_id>_<cell_id>" (10x cell ids repeat across slides)

Negative controls are kept (as per-cell counts) because they are the background model for
probe-specificity questions; they are never mixed into X.
"""

from __future__ import annotations

import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

CONTROL_TYPES = {
    "Negative Control Probe": "ctrl_probe",
    "Negative Control Codeword": "ctrl_codeword",
    "Unassigned Codeword": "ctrl_unassigned",
    "Blank Codeword": "ctrl_unassigned",
    "Genomic Control": "ctrl_genomic",
    "Deprecated Codeword": "ctrl_deprecated",
}


def _read_10x_mtx(mtx_dir: str):
    m = scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz")).T.tocsr()   # cells x features
    feats = pd.read_csv(os.path.join(mtx_dir, "features.tsv.gz"), sep="\t", header=None)
    feats = feats.iloc[:, :3]
    feats.columns = ["id", "name", "type"][: feats.shape[1]]
    if "type" not in feats:
        feats["type"] = "Gene Expression"
    bc = pd.read_csv(os.path.join(mtx_dir, "barcodes.tsv.gz"), header=None)[0].astype(str).to_numpy()
    return m, feats, bc


def _read_10x_h5(path: str):
    import h5py
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(g["shape"][:])
        m = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).T.tocsr()
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])  # noqa: E731
        feats = pd.DataFrame({"id": dec(g["features/id"][:]), "name": dec(g["features/name"][:]),
                              "type": dec(g["features/feature_type"][:])})
        bc = dec(g["barcodes"][:])
    return m, feats, bc


def read_feature_matrix(seg_dir: str):
    """cells x features CSR, features DataFrame (id, name, type), barcodes."""
    mtx_dir = os.path.join(seg_dir, "cell_feature_matrix")
    h5 = os.path.join(seg_dir, "cell_feature_matrix.h5")
    if os.path.isdir(mtx_dir):
        return _read_10x_mtx(mtx_dir)
    if os.path.exists(h5):
        return _read_10x_h5(h5)
    raise FileNotFoundError(f"No cell_feature_matrix/ or cell_feature_matrix.h5 in {seg_dir}")


def read_cells_table(seg_dir: str) -> pd.DataFrame:
    pq = os.path.join(seg_dir, "cells.parquet")
    csv = os.path.join(seg_dir, "cells.csv.gz")
    if os.path.exists(pq):
        cells = pd.read_parquet(pq)
    elif os.path.exists(csv):
        cells = pd.read_csv(csv)
    else:
        raise FileNotFoundError(f"No cells.parquet or cells.csv.gz in {seg_dir}")
    cells["cell_id"] = cells["cell_id"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
    for col in ("cell_area", "nucleus_area"):
        if col not in cells:
            cells[col] = np.nan
    return cells


def load_xenium(seg_dir: str, slide_id: str, segmentation: str = "xenium") -> ad.AnnData:
    """One slide -> AnnData (see module docstring for the layout)."""
    m, feats, bc = read_feature_matrix(seg_dir)
    cells = read_cells_table(seg_dir).set_index("cell_id")
    keep = np.isin(bc, cells.index)
    if (~keep).any():
        print(f"{slide_id}: {(~keep).sum()} matrix cells without centroid dropped")
    m, bc = m[keep], bc[keep]
    cells = cells.loc[bc]

    is_gene = (feats["type"] == "Gene Expression").to_numpy()
    obs = pd.DataFrame(index=[f"{slide_id}_{b}" for b in bc])
    obs["slide_id"] = slide_id
    obs["cell_id"] = bc
    obs["x_um"] = cells["x_centroid"].to_numpy(float)
    obs["y_um"] = cells["y_centroid"].to_numpy(float)
    obs["cell_area"] = cells["cell_area"].to_numpy(float)
    obs["nucleus_area"] = cells["nucleus_area"].to_numpy(float)
    obs["segmentation"] = segmentation
    n_ctrl = {}
    for ftype, short in CONTROL_TYPES.items():
        sel = (feats["type"] == ftype).to_numpy()
        if sel.any():
            prev = obs[f"n_{short}"].to_numpy() if f"n_{short}" in obs else 0
            obs[f"n_{short}"] = prev + np.asarray(m[:, sel].sum(1)).ravel()
            n_ctrl[short] = n_ctrl.get(short, 0) + int(sel.sum())

    X = m[:, is_gene].astype(np.float32)
    var = pd.DataFrame(index=feats.loc[is_gene, "name"].to_numpy())
    var["feature_id"] = feats.loc[is_gene, "id"].to_numpy()
    obs["n_counts"] = np.asarray(X.sum(1)).ravel()
    obs["n_genes"] = np.asarray((X > 0).sum(1)).ravel()

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = obs[["x_um", "y_um"]].to_numpy()
    a.uns["tma"] = {"n_control_features": {slide_id: n_ctrl}, "n_genes_panel": {slide_id: int(is_gene.sum())}}
    return a


def load_slides(slides: pd.DataFrame, base_dir: str = ".", segmentation: str = "xenium",
                seg_root: str | None = None) -> ad.AnnData:
    """All slides of slides.csv (columns slide_id, xenium_dir) -> one AnnData.

    segmentation != "xenium" reads <seg_root>/<segmentation>/<slide_id>/ (cellpose / segger
    outputs of experiment1_tma 03a/03b, same 10x layout).
    """
    parts, n_ctrl, n_genes = [], {}, {}
    for row in slides.itertuples(index=False):
        if segmentation == "xenium":
            d = row.xenium_dir if os.path.isabs(row.xenium_dir) else os.path.join(base_dir, row.xenium_dir)
        else:
            d = os.path.join(seg_root or os.path.join(base_dir, "data", "segmentation"), segmentation, row.slide_id)
        print(f"Loading {row.slide_id} [{segmentation}] from {d}")
        a = load_xenium(d, row.slide_id, segmentation)
        n_ctrl.update(a.uns["tma"]["n_control_features"])
        n_genes.update(a.uns["tma"]["n_genes_panel"])
        parts.append(a)
    adata = ad.concat(parts, join="outer", merge="first") if len(parts) > 1 else parts[0]
    for col in [c for c in adata.obs if c.startswith("n_ctrl_")]:
        adata.obs[col] = adata.obs[col].fillna(0)
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
    adata.uns["tma"] = {"n_control_features": n_ctrl, "n_genes_panel": n_genes}
    return adata
