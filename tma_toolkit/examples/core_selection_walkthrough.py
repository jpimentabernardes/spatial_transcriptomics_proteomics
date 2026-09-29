# %% [markdown]
# # Working with 200 TMA cores: metadata, selection, comparison
#
# Run cell by cell in VS Code / Jupyter (the `# %%` markers make cells), after
# building the object once with scripts/build_tma_object.py.
#
# Everything is organised around ONE table with one row per core (`cohort.cores`):
# the metadata sheet you filled in, plus core geometry and QC. Selecting cores is
# a filter on that table; the cells follow.

# %%
import sys
sys.path.insert(0, "..")          # or: pip install -e tma_toolkit
import matplotlib.pyplot as plt
from tma_toolkit import TMACohort

cohort = TMACohort.read("/work/Spatial_TMA/data/processed/tma.h5ad")   # ADAPT
cohort

# %% What is in the TMA?
cohort.summary(by=("organ", "tissue_type"))

# %% One row per core -- sort, filter, export like any DataFrame
cohort.cores.sort_values("median_counts").head(10)       # worst cores first

# %% [markdown]
# ## Selecting cores
# Keyword filters: value, list of values (any of), or (low, high) range.
# `query=` takes any pandas expression. Everything combines with AND.
# The `include` column of the sheet is honoured unless include_only=False.

# %%
colon      = cohort.select(organ="colon")
tumours    = cohort.select(tissue_type="tumor")
gut        = cohort.select(organ=["colon", "small_intestine", "stomach"])
older      = cohort.select(age=(60, 90))
good_qc    = cohort.select(query="n_cells >= 500 and median_counts >= 20")
by_id      = cohort.select(cores=["TMA2_A01", "TMA2_B07"])

colon_tumour_good = colon & tumours & good_qc           # & | - combine selections
colon_tumour_good.describe()

# %% Quick look: 3 random cores per organ
quick = good_qc.sample(3, by="organ")
fig = quick.plot_cores("n_counts")                      # any obs column: cell_type, n_counts, ...

# %% [markdown]
# ## Named core sets: make selections reproducible
# Saved inside the .h5ad (and exportable to CSV), so a figure made in June can be
# regenerated from exactly the same cores in December.

# %%
cohort.save_core_set("colon_tumour_goodQC", colon_tumour_good, "colon tumour cores with >= 500 cells")
cohort.save_core_set("gut_vs_other_immune", good_qc, "all QC-passing cores, for immune comparisons")
cohort.export_core_sets("/work/Spatial_TMA/tables/core_sets.csv")
# cohort.import_core_sets("my_curated_sets.csv")        # core_set, core_id, description
cohort.write("/work/Spatial_TMA/data/processed/tma.h5ad")

sel = cohort.core_set("colon_tumour_goodQC")

# %% [markdown]
# ## Analysing a subset
# `.adata()` returns an AnnData with just those cells -> use scanpy / squidpy as usual.

# %%
sub = sel.adata()
# import scanpy as sc; sc.pp.normalize_total(sub); sc.pp.log1p(sub); sc.tl.pca(sub) ...

# %% [markdown]
# ## Comparing cores and groups of cores
# Unit = core by default; unit="patient_id" averages a patient's cores first
# (use it when patients contribute several cores -- they are not independent).
# Needs a cell-type column in obs (from your annotation) for composition.

# %%
sel = good_qc
comp = sel.compare_composition("cell_type", group_col="organ")             # per cell type test
comp_patient = sel.compare_composition("cell_type", group_col="organ", unit="patient_id")
comp.head()

# %%
sel.plot_composition("cell_type", group_col="organ")

# %% Differential expression between groups of cores (pseudobulk)
de = sel.differential("tissue_type", "tumor", "normal")                      # all cells
de_t = sel.differential("tissue_type", "tumor", "normal",                     # only T cells
                        label="cell_type", label_value="T_cell", unit="patient_id")
de.head(20)

# %% Gene x core tables, e.g. immune targets
imm = sel.gene_by_core(["PTPRC", "CD3E", "CD8A", "CD68"], stat="frac")        # fraction of positive cells
from tma_toolkit.plotting import plot_group_values
plot_group_values(imm, sel.cores, "organ", list(imm.columns), ylabel="fraction of cells > 0")

# %% Core-to-core similarity: do cores group by organ, patient or slide?
corr = sel.similarity()
sel.plot_similarity(corr, annotate=("organ", "tissue_type", "slide_id"))
rel = sel.similarity_by_relation(corr, same=("organ", "patient_id"))
rel.groupby("relation").r.describe()     # same organ+patient vs same organ, different patient, ...

# %% Add a new core-level column later (e.g. pathologist score) -- lands on cells too
# import pandas as pd
# scores = pd.read_csv("pathology_scores.csv", index_col="core_id")      # columns: til_score, ...
# cohort.add_core_columns(scores)
# cohort.select(til_score=(2, 3))
