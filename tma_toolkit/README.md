# tma_toolkit: core metadata, selection and comparison for TMAs (Python)

A Python toolkit for large TMAs, such as the next 200-core dataset. It builds one AnnData
for all slides with one row of metadata per core. You can then select cores by any
metadata, save those selections, and compare cores or groups of cores. It works on Xenium
output, and on cellpose/segger segmentations written in the same 10x layout
(`experiment1_tma` 03a/03b).

```
pip install -e tma_toolkit            # + [de] for pydeseq2, [analysis] for scanpy/squidpy
```

## 1. Metadata sheet (one row per core)

```
python tma_toolkit/scripts/make_core_metadata_template.py --slide-id TMA2 --rows 10 --cols 20 \
       --extra grade stage msi_status --out config/core_metadata.csv
```
This writes 200 rows, `TMA2_A01 … TMA2_J20`, with columns `core_id, slide_id, core_row, core_col,
patient_id, block_id, organ, site, tissue_type, diagnosis, is_control, replicate, include,
exclusion_reason, notes, sex, age, treatment` plus your extra columns. Fill it in Excel and save
as CSV. Any column you add becomes selectable and comparable. Check the sheet with
`--check config/core_metadata.csv`, which also:

- unifies spelling variants (`Colon` / `colon ` / `COLON`) and reports them
- flags duplicate core ids and two cores sharing one grid position
- flags patient-level columns (`sex`, `age`, `treatment`) that disagree between one patient's cores
- converts numbers typed as text into numbers
- checks controlled vocabularies (`validate_core_metadata(df, allowed_values={"organ": [...]})`)

## 2. Build the object

```
python tma_toolkit/scripts/build_tma_object.py --slides config/slides.csv \
       --core-metadata config/core_metadata.csv --project-dir /work/Spatial_TMA \
       --out /work/Spatial_TMA/data/processed/tma.h5ad
```
The build does the following:

1. Loads every slide.
2. Finds cores automatically. It uses DBSCAN, then corrects the array's rotation, then finds
   rows and columns from the gaps between cores. A missing core or a missing interior row or
   column does not shift the labels.
3. Matches each core to the sheet. It reports sheet cores not found on the slide, and slide
   cores missing from the sheet.
4. Flags cell QC and computes per-core QC.
5. Writes `tma.h5ad` plus `tables/cores.csv`.

**Check `figures/dearray_<slide>.png` against the TMA map.** Fix labels with
`--overrides` (`slide_id, detected_core, core_row, core_col`), or draw cores in Xenium Explorer
and pass `--polygons-dir`.

Tested on a synthetic 10 × 20 array rotated by 2°, with one full column and 5 random cores
missing: every present core got its correct label (`pytest tma_toolkit/tests`).

## 3. Select, save, compare

```python
from tma_toolkit import TMACohort
cohort = TMACohort.read("tma.h5ad")
cohort.cores                                   # one row per core: metadata + geometry + QC
sel = cohort.select(organ=["colon", "lung"], tissue_type="tumor", age=(50, 80),
                    query="n_cells >= 500")
sel = sel - cohort.select(cores=["TMA2_C07"])  # &, |, - combine selections
cohort.save_core_set("colon_lung_tumour", sel, "for the immune comparison")   # stored in the h5ad
sub = sel.adata()                              # AnnData of those cells -> scanpy / squidpy

sel.compare_composition("cell_type", group_col="organ")                    # per core
sel.compare_composition("cell_type", group_col="organ", unit="patient_id") # per patient
sel.differential("tissue_type", "tumor", "normal", label="cell_type", label_value="T_cell")
sel.similarity()                               # core x core correlation (pseudobulk)
sel.similarity_by_relation(same=("organ", "patient_id"))
sel.gene_by_core(["PTPRC", "CD3E"], stat="frac")
sel.plot_cores("cell_type"); sel.plot_composition("cell_type", "organ"); sel.plot_similarity()
```

All statistics treat the **core** (or the **patient**) as the replicate, never the cell.
Differential expression uses pydeseq2 if it is installed, otherwise a Welch t-test on logCPM.
`examples/core_selection_walkthrough.py` walks through all of this cell by cell.

## Layout of the AnnData

| where | what |
|---|---|
| `X`, `layers["counts"]` | raw gene counts (controls excluded) |
| `obs` | `slide_id, core_id`, every core-metadata column, `n_counts, n_genes, cell_area`, `n_ctrl_*` (negative-control counts), `qc_*` flags |
| `obsm["spatial"]` | cell centroids in µm |
| `uns["tma"]["cores"]` | the core table (one row per core) |
| `uns["tma"]["core_sets"]` | named core selections |
| `uns["tma"]["n_control_features"]` | control probes per slide (background rate per probe) |
