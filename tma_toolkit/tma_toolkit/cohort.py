"""TMACohort: the cells (AnnData) + one row per core (metadata + QC), with core
selection, named core sets and comparisons.

    cohort = TMACohort.read("tma.h5ad")
    cohort.cores                                   # DataFrame, one row per core
    colon = cohort.select(organ="colon", tissue_type=["normal", "tumor"])
    big   = cohort.select(query="n_cells > 1000 and median_counts >= 30")
    both  = colon & big                            # also |, - between selections
    cohort.save_core_set("colon_good", both, "colon cores passing QC")
    sub   = both.adata()                           # AnnData of just those cells
    both.compare_composition("cell_type", group_col="tissue_type")
    both.differential("tissue_type", "tumor", "normal", label="cell_type", label_value="T_cell")
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import compare as C
from . import plotting as P


class TMACohort:
    def __init__(self, adata, cores: pd.DataFrame | None = None, core_key: str = "core_id"):
        self.adata = adata
        self.core_key = core_key
        if cores is None:
            cores = adata.uns.get("tma", {}).get("cores")
            cores = pd.DataFrame(cores) if cores is not None else None
        if cores is None:
            cores = self._cores_from_obs()
        cores = cores.copy()
        cores.index = cores[core_key].astype(str) if core_key in cores else cores.index.astype(str)
        cores.index.name = None
        self.cores = cores
        sets = adata.uns.get("tma", {}).get("core_sets", "{}")
        self.core_sets = json.loads(sets) if isinstance(sets, str) else dict(sets)

    # ---------------------------------------------------------------- build
    def _cores_from_obs(self) -> pd.DataFrame:
        """Core-level columns = obs columns constant within every core."""
        obs = self.adata.obs[self.adata.obs[self.core_key].notna()]
        g = obs.groupby(self.core_key, observed=True)
        const = [c for c in obs.columns if c != self.core_key and g[c].nunique(dropna=False).max() == 1]
        cores = g[const].first()
        cores.insert(0, self.core_key, cores.index.astype(str))
        return cores

    def refresh_qc(self) -> "TMACohort":
        """(Re)compute per-core QC columns from the cells."""
        qc = C.core_qc(self.adata, self.core_key)
        qc.index = qc.index.astype(str)
        self.cores = self.cores.drop(columns=[c for c in qc.columns if c in self.cores]).join(qc)
        self.cores["n_cells"] = self.cores["n_cells"].fillna(0).astype(int)
        return self

    def add_core_columns(self, table: pd.DataFrame, on: str | None = None, to_cells: bool = True) -> "TMACohort":
        """Add/replace core-level columns (e.g. a new spreadsheet column, a score you
        computed per core). `table` indexed by core_id (or give `on`)."""
        t = table.set_index(on) if on else table
        t.index = t.index.astype(str)
        for c in t.columns:
            self.cores[c] = t[c].reindex(self.cores.index)
            if to_cells:
                v = self.cores[c].reindex(self.adata.obs[self.core_key].astype(str)).to_numpy()
                self.adata.obs[c] = pd.Categorical(v) if self.cores[c].dtype == object else v
        return self

    # ---------------------------------------------------------------- io
    def write(self, path: str):
        """h5ad with the core table and core sets stored in uns['tma']."""
        def as_cat(s: pd.Series) -> pd.Categorical:
            # h5ad cannot store mixed-type object columns; keep missing values missing
            out = s.astype(str)
            out[s.isna().to_numpy()] = np.nan
            return pd.Categorical(out)

        tma = dict(self.adata.uns.get("tma", {}))
        cores = self.cores.copy()
        for c in cores.columns:
            if cores[c].dtype == object:
                cores[c] = as_cat(cores[c])
        tma["cores"] = cores
        tma["core_sets"] = json.dumps(self.core_sets)
        self.adata.uns["tma"] = tma
        for c in self.adata.obs.columns:
            if self.adata.obs[c].dtype == object:
                self.adata.obs[c] = as_cat(self.adata.obs[c])
        self.adata.write_h5ad(path, compression="gzip")
        print(f"wrote {path}: {self.adata.n_obs:,} cells, {len(self.cores)} cores, "
              f"{len(self.core_sets)} core sets")

    @classmethod
    def read(cls, path: str, backed: bool = False, core_key: str = "core_id") -> "TMACohort":
        import anndata as ad
        return cls(ad.read_h5ad(path, backed="r" if backed else None), core_key=core_key)

    # ---------------------------------------------------------------- selection
    def all(self) -> "CoreSelection":
        return CoreSelection(self, self.cores.index.tolist(), "all cores")

    def select(self, query: str | None = None, cores=None, exclude=None, include_only: bool = True,
               **filters) -> "CoreSelection":
        """Pick cores by metadata. All conditions are combined with AND.

        query    pandas query string on the core table, e.g. "organ == 'colon' and age > 60"
        cores    explicit list of core_ids
        exclude  core_ids to drop
        include_only  honour the `include` column of the metadata sheet (default True)
        **filters column=value, column=[values] (any of), or column=(low, high) (inclusive range)
        """
        t = self.cores
        m = pd.Series(True, index=t.index)
        desc = []
        if include_only and "include" in t:
            m &= t["include"].astype(str).str.lower().isin(["true", "1"])
        if query:
            m &= t.index.isin(t.query(query).index)
            desc.append(query)
        if cores is not None:
            wanted = [str(c) for c in cores]
            unknown = sorted(set(wanted) - set(t.index))
            if unknown:
                raise KeyError(f"unknown core_ids: {unknown[:10]}")
            m &= t.index.isin(wanted)
            desc.append(f"{len(wanted)} listed cores")
        for col, val in filters.items():
            if col not in t:
                raise KeyError(f"'{col}' is not a core column. Available: {list(t.columns)}")
            if isinstance(val, tuple) and len(val) == 2:
                m &= t[col].between(*val)
            elif isinstance(val, (list, set, np.ndarray, pd.Index)):
                m &= t[col].astype(str).isin([str(v) for v in val])
            else:
                m &= t[col].astype(str) == str(val)
            desc.append(f"{col}={val}")
        if exclude is not None:
            m &= ~t.index.isin([str(c) for c in exclude])
            desc.append(f"minus {len(list(exclude))} excluded")
        sel = CoreSelection(self, t.index[m].tolist(), "; ".join(desc) or "all included cores")
        if sel.n_cores == 0:
            print("WARNING: selection is empty")
        return sel

    # ---------------------------------------------------------------- core sets
    def save_core_set(self, name: str, selection: "CoreSelection", description: str = "") -> None:
        """Name a selection so every analysis on it is reproducible (stored in the h5ad)."""
        self.core_sets[name] = {"cores": list(selection.core_ids), "description": description or selection.description}

    def core_set(self, name: str) -> "CoreSelection":
        s = self.core_sets[name]
        return CoreSelection(self, s["cores"], f"core set '{name}': {s['description']}")

    def export_core_sets(self, path: str) -> None:
        rows = [{"core_set": k, "core_id": c, "description": v["description"]}
                for k, v in self.core_sets.items() for c in v["cores"]]
        pd.DataFrame(rows).to_csv(path, index=False)

    def import_core_sets(self, path: str) -> None:
        """CSV with columns core_set, core_id[, description] -- e.g. curated in Excel."""
        df = pd.read_csv(path)
        for name, g in df.groupby("core_set"):
            desc = g["description"].dropna().iloc[0] if "description" in g and g["description"].notna().any() else ""
            self.core_sets[name] = {"cores": g["core_id"].astype(str).tolist(), "description": desc}

    # ---------------------------------------------------------------- overview
    def summary(self, by=("slide_id", "organ", "tissue_type")) -> pd.DataFrame:
        by = [b for b in by if b in self.cores]
        agg = {"n_cores": (self.core_key, "size")}
        if "n_cells" in self.cores:
            agg["n_cells"] = ("n_cells", "sum")
        if "patient_id" in self.cores:
            agg["n_patients"] = ("patient_id", "nunique")
        return self.cores.groupby(by, observed=True, dropna=False).agg(**agg).reset_index()

    def __repr__(self):
        return (f"TMACohort: {self.adata.n_obs:,} cells, {len(self.cores)} cores, "
                f"{self.cores['slide_id'].nunique() if 'slide_id' in self.cores else '?'} slides, "
                f"core columns: {list(self.cores.columns)}")


class CoreSelection:
    """A set of cores from a cohort. Combine with & (and), | (or), - (minus)."""

    def __init__(self, cohort: TMACohort, core_ids, description: str = ""):
        self.cohort = cohort
        self.core_ids = [c for c in cohort.cores.index if c in set(map(str, core_ids))]
        self.description = description

    # set algebra
    def __and__(self, o): return CoreSelection(self.cohort, set(self.core_ids) & set(o.core_ids), f"({self.description}) AND ({o.description})")
    def __or__(self, o): return CoreSelection(self.cohort, set(self.core_ids) | set(o.core_ids), f"({self.description}) OR ({o.description})")
    def __sub__(self, o): return CoreSelection(self.cohort, set(self.core_ids) - set(o.core_ids), f"({self.description}) MINUS ({o.description})")

    @property
    def n_cores(self): return len(self.core_ids)

    @property
    def cores(self) -> pd.DataFrame:
        return self.cohort.cores.loc[self.core_ids]

    def cell_mask(self) -> np.ndarray:
        return self.cohort.adata.obs[self.cohort.core_key].astype(str).isin(self.core_ids).to_numpy()

    def adata(self, copy: bool = True):
        """AnnData of the cells in these cores (copy by default; views are fragile)."""
        sub = self.cohort.adata[self.cell_mask()]
        return sub.copy() if copy else sub

    def describe(self, by=("organ", "tissue_type")) -> pd.DataFrame:
        by = [b for b in by if b in self.cores]
        print(f"{self.n_cores} cores -- {self.description}")
        return self.cores.groupby(by, observed=True, dropna=False).size().rename("n_cores").reset_index() if by else None

    def sample(self, n_per_group: int, by: str, seed: int = 0) -> "CoreSelection":
        """E.g. 3 random cores per organ for a quick look."""
        picked = []
        for _, d in self.cores.groupby(by, observed=True):
            picked += d.sample(min(n_per_group, len(d)), random_state=seed).index.tolist()
        return CoreSelection(self.cohort, picked, f"{self.description}; {n_per_group} per {by}")

    def to_csv(self, path: str):
        self.cores.to_csv(path)

    def __repr__(self):
        return f"CoreSelection({self.n_cores} cores: {self.description})"

    # comparisons (unit = core unless unit='patient_id')
    def qc(self) -> pd.DataFrame:
        return C.core_qc(self.adata(copy=False), self.cohort.core_key)

    def composition(self, label: str) -> pd.DataFrame:
        return C.composition(self.adata(copy=False), label, self.cohort.core_key)

    def compare_composition(self, label: str, group_col: str, unit: str = "core_id") -> pd.DataFrame:
        return C.compare_composition(self.adata(copy=False), self.cores, label, group_col, unit, self.cohort.core_key)

    def pseudobulk(self, by=None, **kw):
        return C.pseudobulk(self.adata(copy=False), by or self.cohort.core_key, **kw)

    def similarity(self, genes=None, method="spearman", **kw) -> pd.DataFrame:
        return C.core_similarity(self.adata(copy=False), genes, method, self.cohort.core_key, **kw)

    def similarity_by_relation(self, corr=None, same=("organ", "patient_id")) -> pd.DataFrame:
        return C.similarity_by_relation(corr if corr is not None else self.similarity(), self.cores, same)

    def gene_by_core(self, genes, stat="mean", layer=None) -> pd.DataFrame:
        return C.gene_by_core(self.adata(copy=False), genes, self.cohort.core_key, stat, layer)

    def differential(self, group_col, group_a, group_b, **kw) -> pd.DataFrame:
        return C.differential(self.adata(copy=False), self.cohort.cores, group_col, group_a, group_b, **kw)

    # plots
    def plot_cores(self, color: str, max_cores: int = 30, **kw):
        if self.n_cores == 0:
            raise ValueError(f"no cores to plot -- the selection is empty ({self.description})")
        return P.plot_cores(self.cohort.adata, self.core_ids[:max_cores], color, self.cohort.core_key, **kw)

    def plot_composition(self, label: str, group_col: str | None = None, **kw):
        return P.plot_composition(self.composition(label), self.cores, group_col, **kw)

    def plot_similarity(self, corr=None, annotate=("organ", "patient_id"), **kw):
        return P.plot_similarity(corr if corr is not None else self.similarity(), self.cores, annotate, **kw)
