# Experiment 1: Multi-organ and multi-tumour TMAs (Xenium + H&E)

These scripts cover QC, object prep, segmentation, annotation and first TMA insights on
**Xenium**, plus **H&E** regions for histology. Conventions: editable parameters at the top,
`save_fig()`, explicit sanity checks, and one saved `.rds` per stage.

The CellScape (spatial proteomics) steps — protein prep (06), registration onto Xenium (07 via
H&E, 07b) and RNA↔protein integration (08) — are removed for now. They remain in the git
history. `07_register_modalities.py` stays because 09 uses its H&E mode.

## Aims → where they are answered

| Question | Script / output |
|---|---|
| **Aim 1.** Are the colon + custom probes gut-specific, or do immune targets generalise? | `05` (§1.1–1.4) → `tables/05_aim1_*` |
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
                                     ├─► 03d segger from 10x / cellpose ─► 03e         │
                                     │     (choose SEGMENTATION, re-run 01–02)        │
H&E (Xenium slide) ─► 07 register H&E → Xenium (--mode he) ─► QuPath regions ─► 09 regions → cells
```

| # | Script | Lang | What it does |
|---|---|---|---|
| 01 | `01_xenium_object_prep.R` | R | Loads all slides (any segmentation) into one Seurat v5 object with one FOV per slide. Run metrics, per-cell QC flags, negative-control background. |
| 02a | `02a_make_pseudo_tma.R` | R | **Testing only.** Punches virtual TMA cores out of whole sections (e.g. Swiss rolls), with an inner/middle/outer location per core, so the TMA pipeline can be tested before real TMAs exist. Then run 02 as usual. |
| 02 | `02_tma_dearray.R` | R | Finds cores (DBSCAN, or Xenium Explorer / polygon selections) and matches them to `tma_map.csv`. Adds patient/organ/location, core QC, filtering. |
| 02b | `02b_cellscape_grid_to_xenium.py` | Py | TMA grid made on the CellScape scan (pathology's QuPath export `<slide>_TMA_cores.csv`) moved onto the Xenium cells: orientation (rotation / mirroring), tilt, scale and offset found from the tissue itself, refined per core. Writes `config/core_grid/<slide_id>.csv` for 02, the CellScape→Xenium transform (JSON) and a check plot; per core, whether it has tissue on Xenium, on CellScape, both or neither. |
| 02c | `02c_check_cellscape_xenium_match.py` | Py | Independent check of 02b on the tissue itself: each core's CellScape DNA image vs its Xenium cells in all 8 orientations (orientation test) and vs the 8 neighbouring CellScape cores (label test); summary plot, overlay gallery (CellScape magenta, Xenium green) and a verdict. |
| 03a | `03a_segment_cellpose.py` | Py | Runs cellpose per core (nuclei + expansion, or DAPI + boundary stain) and re-assigns transcripts. |
| 03b | `03b_segment_segger.py` (or `.sh`) + `03b_segger_to_common.py` | Py | segger (GNN, transcript-to-cell assignment), run through its pixi environment and converted to the same format. |
| 03c | `03c_segmentation_benchmark.R` | R | 10x vs cellpose vs segger on the same cores: cells called per method and per core, a scorecard (transcripts assigned, counts per cell, **MECR** = contamination, control rate: which method is better in how many cores), and side-by-side DAPI / 10x / cellpose outlines as in 01. |
| 03d | `03d_segger_priors.py` | Py | segger started twice on the same cores: from the 10x cells and from the cellpose cells (03a masks). Builds both inputs in Xenium format, runs segger, converts (`segger_xenium`, `segger_cellpose` for 03c) and writes, per transcript, its source cell in each segmentation and its segger cell (`transcript_attribution.parquet`). |
| 03e | `03e_segger_prior_comparison.R` | R | Does segger re-segment each input and which input is worse: transcript fate (kept / moved to another cell / removed / added), per-cell fate, genes moved most, agreement of the two priors vs the two segger results, MECR and counts before/after, and a printed verdict. |
| 03f | `03f_segmentation_disagreement.py` | Py | Where 10x, cellpose and segger disagree: scans the cores in 50 um tiles (adjusted Rand index of the three cell assignments of the same transcripts + spread in % transcripts assigned), then draws the worst tiles close up: DAPI, each method's outlines, every transcript coloured by its cell and by the lineage of its gene; plus agreement tiles as reference and a per-organ summary. |
| 03g | `03g_cell_size_tumour_vs_organ.R` | R | TMA_TUMOR vs TMA_ORGAN, per segmentation: cell and nucleus area, transcripts and genes per cell, transcripts per um2; per lineage (immune cells = internal control for slide effects), core medians, all cores and organ-matched (`site`). |
| 04 | `04_annotation_probe_based.R` | R | Cell type calling from the transcriptomics: clustering (`RES_CLUSTERS`), FindAllMarkers, dot plot of canonical markers (`CANONICAL_MARKERS`) + top markers per cluster, each cluster in the tissue. You name the clusters in `CLUSTER_TO_CELLTYPE` / `CLUSTER_TO_LINEAGE` (config.R), starting from the template 04 writes to `annotation/`. |
| 04 (TMA) | `04_annotation_hierarchical.R` | R | Multi-organ TMAs: L1 lineage on all cells of both TMAs, then immune (all organs together), stroma (all organs together) and epithelial / neuroendocrine cells per organ family (`config/annotation_families.csv`, normal + tumour together). Per group: clustering, markers, dot plot, UMAP, organ composition and a label template `annotation/04_<seg>_<group>.csv` (pre-filled suggestions; edit `label`, re-run). |
| 04b | `04b_annotation_per_organ.R` | R | After 04: epithelial / neuroendocrine cell types per organ (`site` in `config/annotation_families.csv`: colon, rectum, ileum, stomach, ...), once normal only and once normal + the tumours of that organ (no Harmony, so tumour vs normal is kept). Templates `annotation/04b_<seg>_Organ_<site>.csv` / `OrganTumour_<site>.csv`; columns `cell_type_organ`, `cell_type_organ_tumour`. |
| 05 | `05_tma_first_insights.R` | R | Aim 1 and Aim 2 analyses (see table above) plus cellular niches. |
| 07 | `07_register_modalities.py` | Py | With `--mode he`: transform from the H&E (Xenium slide) to Xenium µm, used by 09. (Its `--mode cellscape` is parked with the CellScape steps.) |
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
- **H&E on the Xenium slide after the run.** Same section as the transcripts, so the H&E aligns to
  the Xenium DAPI almost exactly and pathology regions map straight onto cells (`09`).
- **You name the clusters.** 04 shows markers per cluster and writes a template; the names go in
  `CLUSTER_TO_CELLTYPE` / `CLUSTER_TO_LINEAGE` in `config/config.R`, then 04 is re-run.

## Getting started

1. Environments: `Rscript env/install_r.R`; `conda env create -f env/py_imaging.yml`
   (cellpose, registration); segger via pixi (see `03b_segment_segger.py`).
2. Fill in `config/slides.csv`, `config/tma_map.csv` (one row per core: patient, organ,
   tissue type, diagnosis, location) and `config/probe_expectations.csv`. That last file lists
   your custom probes and immune targets with the lineage each is expected in, and drives
   Aim 1 §1.3. Check `config/marker_panel_rna.csv` against your actual panel. Missing markers are reported, not silently used.
3. Set `EXP1_PROJECT_DIR` (data/output root) and run R scripts from inside `experiment1_tma/`:
   ```
   export EXP1_PROJECT_DIR=/work_beegfs/sukmb430/Spatial_TMA
   Rscript 01_xenium_object_prep.R
   # TMA grid: QuPath TMA dearrayer on the Xenium morphology image + qupath/export_tma_grid.groovy
   #           -> config/core_grid/<slide_id>.csv (otherwise 02 finds the grid itself)
   #   or, from the CellScape dearray (pathology's <slide>_TMA_cores.csv), one call per slide:
   python 02b_cellscape_grid_to_xenium.py --slide-id TMA_ORGAN --sector A --cores-csv FL0027_TMA_cores.csv --xenium-dir <xenium_outs>
   python 02b_cellscape_grid_to_xenium.py --slide-id TMA_TUMOR --sector C --cores-csv FL0026_TMA_cores.csv --xenium-dir <xenium_outs>
   #   manual corrections: import config/core_grid/<slide>_grid_qupath.geojson into QuPath on the Xenium
   #   morphology image, fix cores, run qupath/export_core_edits.groovy -> config/core_grid/<slide>_edited.csv
   #   (02 and 02c use it instead of <slide>.csv)
   python 02c_check_cellscape_xenium_match.py --slide-id TMA_ORGAN --sector A --cores-csv FL0027_TMA_cores.csv --cellscape-image FL0027_TMA_preview.png --xenium-dir <xenium_outs>
   Rscript 02_tma_dearray.R                                           # then CHECK figures/02_dearray_*.png
   python 03a_segment_cellpose.py --slide-id ... --cores-csv $EXP1_PROJECT_DIR/tables/02_cores_all.csv ...
   python 03b_segment_segger.py --slide-id <slide> --xenium-dir <xenium_outs> --project-dir $EXP1_PROJECT_DIR --segger-repo <segger clone>
   python 03d_segger_priors.py --slide-id <slide> --xenium-dir <xenium_outs> --project-dir $EXP1_PROJECT_DIR --segger-repo <segger clone>   # GPU; needs 03a --save-masks
   Rscript 03c_segmentation_benchmark.R                               # choose SEGMENTATION; re-run 01-02
   python 03f_segmentation_disagreement.py --slide-id TMA_ORGAN --xenium-dir <outs> --project-dir <project>
   Rscript 03g_cell_size_tumour_vs_organ.R                            # tumour vs normal cell size / expression
   Rscript 03e_segger_prior_comparison.R                              # segger from 10x vs from cellpose
   Rscript 04_annotation_probe_based.R                                # single-tissue data (e.g. Swiss rolls)
   Rscript 04_annotation_hierarchical.R                               # multi-organ TMAs: edit annotation/04_<seg>_*.csv, re-run
   Rscript 04b_annotation_per_organ.R                                 # per organ, with and without tumours
   Rscript 05_tma_first_insights.R
   python 07_register_modalities.py --mode he --project-dir $EXP1_PROJECT_DIR --slide-id <slide> \
       --xenium-dir <xenium outs> --moving-image <xenium-slide H&E> --refine image
   (QuPath annotations on the Xenium-slide H&E) ; Rscript 09_he_regions.R
   ```

### Python steps as Jupyter notebooks

Every Python step also exists as a notebook in `notebooks/`: 03a, 03b_segment_segger (runs segger
through pixi, then converts), 03b_segger_to_common (conversion only), 03d and 07. Each notebook has the same code as its script, laid out for
interactive use:

- one **Parameters** cell holds the script's command-line options (`--slide-id` → `slide_id`);
- every helper function has its own cell;
- the steps of `main()` are separate cells, so you can look at the images, transforms and
  tables after each step.

Start Jupyter from inside `experiment1_tma/notebooks/`, or anywhere if you open the notebook
from that folder, because the notebooks find the code one folder up. They must run in the
same conda environment as the scripts (`exp1-imaging` from `env/py_imaging.yml`), not Jupyter's
default Python. Register the environment as a
Jupyter kernel once:
```
conda activate exp1-imaging
conda install -c conda-forge ipykernel      # already in the .yml for new environments
python -m ipykernel install --user --name exp1-imaging --display-name exp1-imaging
```
then choose it in Jupyter with **Kernel → Change kernel → exp1-imaging**. The first cell of
each notebook prints which Python it runs on and says this if numpy is missing.

The `.py` scripts remain the source and are what to run as cluster batch jobs. After changing a script, regenerate the
notebooks from the repo root with `python tools/make_notebooks.py`. This **overwrites** the
notebooks, so save your own notebook edits under another name.

## Things to verify with the real data (flagged in the scripts as ADAPT / CHECK)

- **Species / gene symbols.** The marker files use human symbols.
- **QC thresholds** (`MIN_COUNTS`, areas) are set for a few-hundred-gene panel. Check the
  01/02 figures before accepting them.
- **Xenium Explorer alignment matrix direction** (only for `07_register_modalities.py --mode he`).
  The script tests both directions and prints which one it used. Confirm on
  `figures/09_he_regions_check_*.png`.
- **segger version.** Commands were checked against the current repository (`segger segment` /
  `segger export`), but it changes fast. Run `--help` once on the
  version you installed.
- **Confounding in the TMA design.** If all cores of one organ sit on one slide, organ and
  slide effects cannot be separated (`05` §2.2 says so). Place bridge cores (for example the
  same colon or tonsil block on every TMA) when designing the next TMAs.
