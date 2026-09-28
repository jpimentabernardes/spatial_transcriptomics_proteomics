# Experiment 1: Multi-organ and multi-tumour TMAs (Xenium + CellScape + H&E)

These scripts cover QC, object prep, segmentation, annotation and first TMA insights on
**Xenium**. They then merge in **CellScape spatial proteomics** from the consecutive section
and **H&E** for histology. The class scripts in the repo root (`00x_*.R`,
MALDI + Xenium) are left unchanged. The patterns from them that proved to work are reused
here: editable parameters at the top, `save_fig()`, explicit sanity checks, and one saved
`.rds` per stage.

## Aims → where they are answered

| Question | Script / output |
|---|---|
| **Aim 1.** Are the colon + custom probes gut-specific, or do immune targets generalise? | `05` (§1.1–1.4) → `tables/05_aim1_*`; checked against protein in `08` → `tables/08_aim1_probe_validation_by_protein.csv` |
| – is a probe above background in organ X? | `05` §1.1: signal/background vs negative-control probes, **per core** |
| – organ-specific, or just more immune cells in that organ? | `05` §1.2: tau on all cells **vs** tau within immune cells |
| – does the immune probe signal land in immune cells in every organ? | `05` §1.3: on-target enrichment + off-target/background ratio (read with MECR from `03c`) |
| – do T cells/macrophages look the same across organs through this panel? | `05` §1.4: between-organ correlation of immune-target profiles per cell type |
| **Aim 2.** How heterogeneous is the data? Can we compare organs, locations, patients? | `05` §2.1–2.5, `09` |
| – organ vs patient vs slide (batch) vs cell type | `05` §2.2: variance partitioning of pseudobulk (core × cell type) |
| – same-organ, different-patient vs same-patient similarity | `05` §2.3 |
| – locations within organ (proximal/distal, tumour centre/margin) | `05` §2.4 (from `location` in the TMA map), `09` (from H&E regions) |

**The statistical unit is the core (or patient), never the cell.** All Aim 1 and Aim 2 statistics use
per-core summaries or pseudobulk.

## Pipeline

```
Xenium outs ─► 01 object + QC ─► 02 TMA dearray + metadata ─┬─► 04 annotation ─► 05 Aim 1 / Aim 2
                                     │                       │                        │
                                     ├─► 03a cellpose ─┐     │                        │
                                     ├─► 03b segger  ──┴─► 03c benchmark              │
                                     │     (choose SEGMENTATION, re-run 01–02)        │
CellScape ─► (06a segment) ─► 06 protein prep + dearray ─► 07 register ─► (07b non-rigid) ─► 08 RNA↔protein
H&E ───────────────────────────────────────────────► 07 --mode he ─► QuPath regions ─► 09 regions → cells
```

| # | Script | Lang | What it does |
|---|---|---|---|
| 01 | `01_xenium_object_prep.R` | R | Loads all slides (any segmentation) into one Seurat v5 object with one FOV per slide. Run metrics, per-cell QC flags, negative-control background. |
| 02 | `02_tma_dearray.R` | R | Finds cores (DBSCAN, or Xenium Explorer / polygon selections) and matches them to `tma_map.csv`. Adds patient/organ/location, core QC, filtering. |
| 03a | `03a_segment_cellpose.py` | Py | Runs cellpose per core (nuclei + expansion, or DAPI + boundary stain) and re-assigns transcripts. |
| 03b | `03b_segment_segger.sh` + `03b_segger_to_common.py` | sh/Py | segger (GNN, transcript-to-cell assignment), converted to the same format. |
| 03c | `03c_segmentation_benchmark.R` | R | 10x vs cellpose vs segger on the same cores: transcripts assigned, counts per cell, **MECR** (contamination), per organ. |
| 04 | `04_annotation_probe_based.R` | R | Hierarchical, probe-based, semi-manual annotation (lineage → immune/stromal/epithelial subsets). Cluster→label decisions go in editable CSVs. |
| 05 | `05_tma_first_insights.R` | R | Aim 1 and Aim 2 analyses (see table above) plus cellular niches. |
| 06a | `06a_cellscape_segment.py` | Py | Only if you have CellScape images without a per-cell export: cellpose + per-channel mean intensities. |
| 06 | `06_cellscape_prep.R` | R | CellScape → Seurat (`PROT` assay, arcsinh). QC, dearray with the same TMA map, annotation with the **same lineage names** as the RNA. |
| 07 | `07_register_modalities.py` | Py | CellScape → Xenium: global affine from **matched core centres**, then per-core refinement (label-aware ICP, or DAPI↔DAPI mutual information). H&E → Xenium from landmarks or a Xenium Explorer alignment matrix. |
| 07b | `07b_nonrigid_refine.py` | Py | Optional non-rigid refinement per core with **GEASO** or **Spateo**, using a shared representation (neighbourhood lineage composition + paired markers). Rejected automatically if it doesn't help. |
| 08 | `08_integrate_rna_protein.R` | R | RNA↔protein at three levels: cell (nearest cell / neighbourhood protein), 50 µm bins, and core × cell type (no registration needed). Also cross-validates Aim 1. |
| 09 | `09_he_regions.R` | R | Pathology regions drawn on H&E (QuPath GeoJSON) → transformed → assigned to cells → composition and immune targets per region. |

## Key design decisions (and why)

- **Object type.** One Seurat v5 object for all slides, one centroid FOV per slide, cells
  named `<slide>_<cell>`. Negative controls are kept as per-cell metadata because they are
  Aim 1's background model. Polygons are not loaded: for multi-slide TMAs they cost many GB,
  no statistic needs them, and they broke `subset()` in the class scripts. Python tools talk
  to R through plain files (10x-style matrix + centroid CSV), with no format converters.
- **The TMA map is the source of truth** for sample identity (`config/tma_map.csv`). Dearraying
  is checked visually (`figures/02_dearray_<slide>.png`) and can be overridden
  (`config/core_overrides.csv`, or per-core polygons in `config/core_selections/<slide>/`).
- **One segmentation for everything.** `03c` compares methods on a subset of cores. You then
  set `SEGMENTATION` in `config/config.R` and re-run 01–02; core definitions are reused.
  Mixing segmentations between organs would make Aim 1/2 uninterpretable.
- **Semi-manual annotation lives in files.** Scores suggest a label per cluster in
  `annotation/*_cluster_map_*.csv`. You edit `label`, set `reviewed=TRUE` and re-run.
- **Serial-section integration is done at the right scale.** A Xenium cell and its nearest
  CellScape cell are not the same cell, so the most important comparison (`08` C) is
  core × cell type, which needs no registration at all. The registration is still useful
  for spatial questions and is quality-controlled per core (lineage-agreement lift).
- **Registration goes coarse → fine, and each step can be rejected.** Global affine from
  matched cores (free landmarks in a TMA), then per-core rigid refinement, then optional
  non-rigid refinement. A refinement is kept only if it improves cross-modal agreement
  within shift/rotation limits.
- **H&E is best done on the Xenium slide after the run.** That makes it the same section, so
  a single affine is essentially exact and H&E regions map onto cells directly (`09`). Both
  omics layers already share DAPI, so DAPI (not H&E) is the anchor for RNA↔protein.

## Getting started

1. Environments: `Rscript env/install_r.R`; `conda env create -f env/py_imaging.yml`
   (cellpose, registration); segger via pixi (see `03b_segment_segger.sh`); Spateo/GEASO via
   `env/py_align.yml`.
2. Fill in `config/slides.csv`, `config/tma_map.csv` (one row per core: patient, organ,
   tissue type, diagnosis, location) and `config/probe_expectations.csv`. That last file lists
   your custom probes and immune targets with the lineage each is expected in, and drives
   Aim 1 §1.3. Check `config/marker_panel_rna.csv` / `marker_panel_protein.csv` against your
   actual panels. Missing markers are reported, not silently used.
3. Set `EXP1_PROJECT_DIR` (data/output root) and run R scripts from inside `experiment1_tma/`:
   ```
   export EXP1_PROJECT_DIR=/work_beegfs/sukmb430/Spatial_TMA
   Rscript 01_xenium_object_prep.R && Rscript 02_tma_dearray.R      # then CHECK figures/02_dearray_*.png
   python 03a_segment_cellpose.py --slide-id ... --cores-csv $EXP1_PROJECT_DIR/tables/02_cores_all.csv ...
   bash 03b_segment_segger.sh <slide> <xenium_outs> $EXP1_PROJECT_DIR
   Rscript 03c_segmentation_benchmark.R                               # choose SEGMENTATION; re-run 01-02
   Rscript 04_annotation_probe_based.R                                # edit annotation/*.csv, re-run
   Rscript 05_tma_first_insights.R
   Rscript 06_cellscape_prep.R                                        # CHECK figures/06_dearray_*.png
   python 07_register_modalities.py --project-dir $EXP1_PROJECT_DIR --slide-id <slide> --mode cellscape
   Rscript 08_integrate_rna_protein.R
   python 07_register_modalities.py ... --mode he ...; (QuPath annotations); Rscript 09_he_regions.R
   ```

## Things to verify with the real data (flagged in the scripts as ADAPT / CHECK)

- **CellScape export columns** (`CS_*` in `config.R`): column names and units differ between
  software versions. `06` stops with the list of columns if they don't match.
- **Species / gene symbols.** The marker files use human symbols.
- **QC thresholds** (`MIN_COUNTS`, areas) are set for a few-hundred-gene panel. Check the
  01/02 figures before accepting them.
- **Xenium Explorer alignment matrix direction** (H&E mode of `07`). The script tests both
  directions and prints which one it used. Confirm on `figures/09_he_regions_check_*.png`.
- **segger / GEASO / Spateo versions.** Commands were checked against the current
  repositories (segger `segger segment`/`segger export`; GEASO `coarse_to_fine_alignment`;
  Spateo `st.align.morpho_align`), but these tools change fast. Run `--help` once on the
  version you installed.
- **Confounding in the TMA design.** If all cores of one organ sit on one slide, organ and
  slide effects cannot be separated (`05` §2.2 says so). Place bridge cores (for example the
  same colon or tonsil block on every TMA) when designing the next TMAs.
