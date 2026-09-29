"""Comparisons between cores and groups of cores.

The statistical unit is the CORE (or the patient, via unit="patient_id") --
never the single cell. With millions of cells, cell-level tests make every
difference "significant"; cores and patients are the real replicates.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats


def _counts(adata, layer="counts"):
    X = adata.layers[layer] if layer in adata.layers else adata.X
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    tmp = np.empty(n)
    tmp[order] = ranked
    out[ok] = tmp
    return out


# ---------------------------------------------------------------------------
# per-core tables
# ---------------------------------------------------------------------------
def core_qc(adata, core_key="core_id") -> pd.DataFrame:
    """n_cells, median counts/genes/area and negative-control rate per core."""
    obs = adata.obs
    ctrl = [c for c in obs if c.startswith("n_ctrl_")]
    g = obs.groupby(core_key, observed=True)
    out = pd.DataFrame({
        "n_cells": g.size(),
        "median_counts": g["n_counts"].median(),
        "median_genes": g["n_genes"].median(),
        "median_area": g["cell_area"].median() if "cell_area" in obs else np.nan,
    })
    if ctrl:
        tot_ctrl = obs[ctrl].sum(1).groupby(obs[core_key], observed=True).sum()
        out["ctrl_per_1000_counts"] = 1000 * tot_ctrl / g["n_counts"].sum()
    return out


def composition(adata, label: str, core_key="core_id", normalize=True) -> pd.DataFrame:
    """cores x labels table of fractions (or counts)."""
    t = pd.crosstab(adata.obs[core_key], adata.obs[label])
    return t.div(t.sum(1), axis=0) if normalize else t


def pseudobulk(adata, by: str | list[str] = "core_id", layer="counts", min_cells: int = 20,
               genes: list[str] | None = None):
    """Sum raw counts per group. Returns an AnnData (groups x genes) whose obs holds
    the grouping columns, n_cells, and every core-level column that is constant within
    the group (organ, patient_id, ...)."""
    import anndata as ad
    by = [by] if isinstance(by, str) else list(by)
    obs = adata.obs
    key = obs[by].astype(str).agg("||".join, axis=1)
    n = key.value_counts()
    keep_groups = n[n >= min_cells].index
    keep = key.isin(keep_groups).to_numpy()
    codes, uniq = pd.factorize(key[keep])
    X = _counts(adata)[keep]
    if genes is not None:
        X = X[:, adata.var_names.get_indexer(genes)]
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    pb = G @ X
    meta = obs.loc[keep].assign(_key=key[keep].to_numpy())
    const_cols = [c for c in meta.columns if c != "_key" and meta.groupby("_key", observed=True)[c].nunique(dropna=False).max() == 1]
    pobs = meta.groupby("_key", observed=True)[const_cols].first().reindex(uniq)
    pobs["n_cells"] = n.reindex(uniq).to_numpy()
    pobs.index = uniq.astype(str)
    var = adata.var.iloc[:, :0].copy() if genes is None else pd.DataFrame(index=genes)
    return ad.AnnData(X=sp.csr_matrix(pb), obs=pobs, var=var)


def log_cpm(pb, prior: float = 1.0) -> np.ndarray:
    X = pb.X.toarray() if sp.issparse(pb.X) else np.asarray(pb.X)
    lib = X.sum(1, keepdims=True)
    return np.log2((X + prior) / (lib + 2 * prior) * 1e6)


def gene_by_core(adata, genes: list[str], core_key="core_id", stat="mean", layer=None) -> pd.DataFrame:
    """cores x genes: mean normalised expression ('mean', uses adata.X or `layer`) or
    fraction of cells with >= 1 count ('frac', uses raw counts)."""
    genes = [g for g in genes if g in adata.var_names]
    idx = adata.var_names.get_indexer(genes)
    if stat == "frac":
        X = _counts(adata)[:, idx] > 0
    else:
        X = adata.layers[layer] if layer else adata.X
        X = X[:, idx]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return pd.DataFrame(X.astype(float), columns=genes, index=adata.obs_names).groupby(
        adata.obs[core_key].to_numpy()).mean()


# ---------------------------------------------------------------------------
# group comparisons
# ---------------------------------------------------------------------------
def _unit_table(values: pd.DataFrame, cores: pd.DataFrame, group_col: str, unit: str):
    """Attach the group of each core and optionally average cores -> patients."""
    d = values.join(cores[[group_col] + ([unit] if unit != "core_id" else [])], how="inner")
    if unit != "core_id":
        d = d.groupby([unit, group_col], observed=True).mean(numeric_only=True).reset_index(level=group_col)
    return d


def compare_composition(adata, cores: pd.DataFrame, label: str, group_col: str,
                        unit: str = "core_id", core_key="core_id") -> pd.DataFrame:
    """Per cell type: are its per-core (or per-patient) fractions different between groups?
    Mann-Whitney for 2 groups, Kruskal-Wallis for more; BH across cell types."""
    comp = composition(adata, label, core_key)
    d = _unit_table(comp, cores, group_col, unit)
    groups = d[group_col].dropna().unique()
    rows = []
    for ct in comp.columns:
        samples = [d.loc[d[group_col] == g, ct].to_numpy() for g in groups]
        samples = [s for s in samples if len(s) > 0]
        if len(samples) < 2:
            continue
        try:
            p = (stats.mannwhitneyu(*samples).pvalue if len(samples) == 2 else stats.kruskal(*samples).pvalue)
        except ValueError:
            p = np.nan
        row = {"cell_type": ct, "p": p}
        row.update({f"mean_{g}": d.loc[d[group_col] == g, ct].mean() for g in groups})
        row.update({f"n_{g}": int((d[group_col] == g).sum()) for g in groups})
        rows.append(row)
    res = pd.DataFrame(rows)
    res["padj"] = bh(res["p"].to_numpy()) if len(res) else []
    return res.sort_values("p")


def core_similarity(adata, genes=None, method="spearman", core_key="core_id", min_cells=20) -> pd.DataFrame:
    """Core x core correlation of pseudobulk logCPM profiles."""
    pb = pseudobulk(adata, core_key, min_cells=min_cells, genes=genes)
    lc = pd.DataFrame(log_cpm(pb).T, columns=pb.obs[core_key].to_numpy())
    return lc.corr(method=method)


def similarity_by_relation(corr: pd.DataFrame, cores: pd.DataFrame,
                           same: tuple[str, ...] = ("organ", "patient_id")) -> pd.DataFrame:
    """Long table of core pairs with 'same X / different X' labels, e.g. to ask whether
    two colon cores of different patients are as similar as two cores of one patient."""
    ids = corr.index.to_numpy()
    iu = np.triu_indices(len(ids), 1)
    out = pd.DataFrame({"core_a": ids[iu[0]], "core_b": ids[iu[1]], "r": corr.to_numpy()[iu]})
    for c in same:
        if c in cores:
            a = cores[c].reindex(out.core_a).to_numpy()
            b = cores[c].reindex(out.core_b).to_numpy()
            out[f"same_{c}"] = a == b
    flags = [f"same_{c}" for c in same if f"same_{c}" in out]
    out["relation"] = out[flags].apply(lambda r: ", ".join(("same " if v else "diff ") + k[5:] for k, v in r.items()), axis=1)
    return out


def differential(adata, cores: pd.DataFrame, group_col: str, group_a, group_b, label: str | None = None,
                 label_value=None, unit: str = "core_id", min_cells: int = 20, method: str = "auto") -> pd.DataFrame:
    """Pseudobulk differential expression between two groups of cores (or patients).

    group_a / group_b: value(s) of `group_col` (a value or list). Optionally restrict to
    one cell type (label, label_value), e.g. compare T cells of colon vs lung cores.
    method: 'pydeseq2' (if installed), 'ttest' (Welch on logCPM), or 'auto'.
    """
    sub = adata
    if label is not None:
        sub = adata[(adata.obs[label] == label_value).to_numpy()]
    as_list = lambda v: list(v) if isinstance(v, (list, tuple, set, np.ndarray)) else [v]  # noqa: E731
    A, B = as_list(group_a), as_list(group_b)
    pb = pseudobulk(sub, unit, min_cells=min_cells)
    grp = cores.drop_duplicates(unit).set_index(unit)[group_col] if unit != "core_id" else cores[group_col]
    g = grp.reindex(pb.obs[unit].astype(str).to_numpy()).to_numpy()
    in_a, in_b = np.isin(g, A), np.isin(g, B)
    pb = pb[in_a | in_b].copy()
    cond = np.where(np.isin(g[in_a | in_b], A), "A", "B")
    na, nb = int((cond == "A").sum()), int((cond == "B").sum())
    if min(na, nb) < 2:
        raise ValueError(f"need >= 2 {unit}s per group, have A={na}, B={nb}")
    if min(na, nb) < 3:
        warnings.warn(f"only {na} vs {nb} {unit}s -- results are exploratory", stacklevel=2)

    if method in ("auto", "pydeseq2"):
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            counts = pd.DataFrame(pb.X.toarray().astype(int), index=pb.obs_names, columns=pb.var_names)
            dds = DeseqDataSet(counts=counts, metadata=pd.DataFrame({"cond": cond}, index=pb.obs_names),
                               design="~cond", quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["cond", "A", "B"], quiet=True)
            st.summary()
            r = st.results_df.rename(columns={"log2FoldChange": "log2FC", "pvalue": "p"})
            r["method"] = "pydeseq2"
            return r[["log2FC", "baseMean", "p", "padj", "method"]].assign(n_a=na, n_b=nb).sort_values("p")
        except ImportError:
            if method == "pydeseq2":
                raise
    lc = log_cpm(pb)
    a, b = lc[cond == "A"], lc[cond == "B"]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    r = pd.DataFrame({"log2FC": a.mean(0) - b.mean(0), "mean_logCPM_a": a.mean(0),
                      "mean_logCPM_b": b.mean(0), "t": t, "p": p}, index=pb.var_names)
    r["padj"] = bh(r["p"].to_numpy())
    r["method"] = "welch_t_logCPM"
    return r.assign(n_a=na, n_b=nb).sort_values("p")
