## =============================================================================
## config/config.R
##
## Single place for paths and parameters. Every R script in experiment1_tma/
## starts with source("config/config.R"), so run scripts with the working
## directory set to experiment1_tma/ (or set EXP1_CODE_DIR).
##
## Anything marked ADAPT is a judgement call that depends on your data --
## check it the first time real data lands, don't just accept the default.
## =============================================================================

## -----------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------
## ADAPT: where the raw data and outputs live on your cluster. Code stays in
## the git repo; data never goes in the repo.
PROJECT_DIR   <- Sys.getenv("EXP1_PROJECT_DIR", "/work_beegfs/sukmb430/Spatial_TMA")
CODE_DIR      <- Sys.getenv("EXP1_CODE_DIR", getwd())

CONFIG_DIR    <- file.path(CODE_DIR, "config")
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
FIGURES_DIR   <- file.path(PROJECT_DIR, "figures")
TABLES_DIR    <- file.path(PROJECT_DIR, "tables")
ANNOT_DIR     <- file.path(PROJECT_DIR, "annotation")   # editable cluster->label maps
SEG_DIR       <- file.path(PROJECT_DIR, "data", "segmentation")
REG_DIR       <- file.path(PROJECT_DIR, "data", "registration")

for (d in c(PROCESSED_DIR, FIGURES_DIR, TABLES_DIR, ANNOT_DIR, SEG_DIR, REG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

SLIDES_CSV        <- file.path(CONFIG_DIR, "slides.csv")
TMA_MAP_CSV       <- file.path(CONFIG_DIR, "tma_map.csv")
MARKERS_RNA_CSV   <- file.path(CONFIG_DIR, "marker_panel_rna.csv")
MARKERS_PROT_CSV  <- file.path(CONFIG_DIR, "marker_panel_protein.csv")
PROBE_EXPECT_CSV  <- file.path(CONFIG_DIR, "probe_expectations.csv")
RNA_PROT_PAIRS_CSV <- file.path(CONFIG_DIR, "rna_protein_pairs.csv")

## Paths inside slides.csv are relative to PROJECT_DIR unless absolute.
resolve_path <- function(p) {
  if (is.na(p) || !nzchar(p)) return(NA_character_)
  if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(PROJECT_DIR, p)
}

## -----------------------------------------------------------------------
## Segmentation used for everything downstream of 03c.
##   "xenium"   = 10x Xenium Onboard Analysis output (default, always exists)
##   "cellpose" = 03a_segment_cellpose.py output
##   "segger"   = 03b_segment_segger.sh + 03b_segger_to_common.py output
## Start with "xenium", run 03a-03c, then switch if the benchmark says so and
## re-run 01 -> 02.
## -----------------------------------------------------------------------
SEGMENTATION <- Sys.getenv("EXP1_SEGMENTATION", "xenium")

## -----------------------------------------------------------------------
## Cell QC (01). ADAPT to panel size: colon panel + custom add-on is a few
## hundred genes, so per-cell counts are far lower than scRNA-seq.
## -----------------------------------------------------------------------
MIN_QV           <- 20     # 10x high-confidence transcript cutoff
MIN_COUNTS       <- 10
MIN_FEATURES     <- 5
MIN_CELL_AREA    <- 6      # um^2 -- below this is almost always debris
MAX_CELL_AREA    <- 600    # um^2 -- hard ceiling; a per-slide MAD ceiling is also applied
AREA_MAD_UPPER   <- 5      # flag cells > median + 5 MAD of log(area)
MAX_CONTROL_FRAC <- 0.05   # flag cells where >5% of counts are negative controls

## QC maps in 01 (whole-slide plots)
QC_MAP_POINT_SIZE <- 0.3     # bigger = denser-looking tissue
QC_MAP_MAX_CELLS  <- 1e6     # cells drawn per figure (random subsample above this)

## -----------------------------------------------------------------------
## TMA dearraying (02). Cores are usually 0.6-2 mm diameter with >= 200 um
## gaps; if neighbouring cores get merged, set a fixed DBSCAN_EPS_UM below the gap.
## -----------------------------------------------------------------------
DBSCAN_EPS_UM        <- NA    # NA = adaptive (3 x median kNN distance); or a fixed value in um
DBSCAN_MIN_PTS       <- 10
MIN_CELLS_PER_CORE   <- 200   # smaller DBSCAN clusters = debris / tissue fragments
CORE_MERGE_DIST_UM   <- 400   # fragments of one torn core closer than this get merged
CORE_MAX_DIST_FACTOR <- 1.25  # assign cells up to 1.25 x core radius from centre

## Virtual TMA for TESTING on whole sections (02a_make_pseudo_tma.R only).
PSEUDO_TMA_SLIDES       <- NULL       # NULL = all slides; or c("TMA_ORGAN")
PSEUDO_CORE_DIAMETER_UM <- 600        # typical TMA needle: 0.6-1.5 mm
PSEUDO_CORE_PITCH_UM    <- 1000       # centre-to-centre distance of the punch grid
PSEUDO_MIN_CELLS        <- 500        # a punch must contain this many QC-pass cells
PSEUDO_PATIENT_ID       <- NULL       # NULL = slide_id; e.g. "mouse1" if both rolls are one animal
PSEUDO_ORGAN            <- "colon"
PSEUDO_CENTERS          <- NULL       # roll centre per slide, e.g. data.frame(slide_id = "TMA_ORGAN", x = 5000, y = 4800)

## Core-level QC: cores failing these are flagged (not removed) in 02.
CORE_MIN_CELLS         <- 500
CORE_MIN_MEDIAN_COUNTS <- 20

## -----------------------------------------------------------------------
## Clustering / annotation (04)
## -----------------------------------------------------------------------
N_PCS            <- 30
USE_HARMONY      <- TRUE     # correct slide effects, NOT organ/patient effects
RES_CLUSTERS     <- 0.8      # 04: clustering of all cells (ST_clusters); higher = more, finer clusters
RES_LEVEL1       <- 0.3      # 06 (CellScape): lineage clustering
RES_LEVEL2       <- 0.6      # 06 (CellScape): sub-clustering within a lineage
TOP_N_PER_CLUSTER <- 5       # 04: top FindAllMarkers genes per cluster shown in the dot plot
MARKERS_MAX_CELLS <- 2000    # 04: cells per cluster used by FindAllMarkers (speed; Inf = all)
## 04: canonical markers always shown in the dot plot, next to the top markers of
## each cluster. Matched to the panel ignoring case (Epcam = EPCAM); genes not on
## the panel are listed and dropped. ADAPT to tissue / panel.
CANONICAL_MARKERS <- c(
  "Epcam", "Mki67", "Lgr5", "Olfm4", "Stmn1", "Atoh1", "Apoa4", "Fabp1", "Alpi",
  "Slc26a3", "Car4", "Agr2", "Muc2", "Lyz1", "Defa5", "Mmp7", "Dclk1", "Trpm5",
  "Chga", "Cd3d", "Cd4", "Cd8a", "Foxp3", "Il7r", "Ccr5", "Gzmb", "Ncam1",
  "Cd79a", "Ighm", "Xcr1", "Itgax", "Itgam", "Csf1r", "Fcgr1", "Ly6c1", "Ly6g", "Siglecf"
)
## 04: cluster -> cell type, YOUR call, from figures/04_top_markers_dotplot_transcriptomics.png
## and 04_ST_clusters_spatial_<slide>.png. 04 writes a ready-to-fill template with
## each cluster's top markers: annotation/cluster_to_celltype_template_<seg>.R --
## copy it here, fill in the names, re-run 04. NULL = cells are labelled cluster_<n>.
## Cluster numbers change with the data, SEGMENTATION and RES_CLUSTERS: re-check
## the mapping whenever those change (04 stops if a cluster has no name).
CLUSTER_TO_CELLTYPE <- NULL
## Lineage per cluster, one of Epithelial / Immune / Fibroblast / Endothelial /
## Smooth_muscle_pericyte / Neural_glia (level-1 names of marker_panel_rna.csv):
## 05 uses lineage == "Immune" and matches probe_expectations.csv to these names.
## NULL = the marker-signature suggestion written in the template.
CLUSTER_TO_LINEAGE <- NULL
USE_UCELL        <- TRUE     # rank-based scores; robust to small panels. Falls back to AddModuleScore

## -----------------------------------------------------------------------
## Aim 1 / Aim 2 (05)
## -----------------------------------------------------------------------
SBR_LOG2_MIN      <- 1       # gene counted "detected" in a core if >= 2x neg-probe background...
SBR_PADJ_MAX      <- 0.01    # ...and Poisson-test BH-adjusted p below this
MIN_CELLS_PSEUDOBULK <- 20   # min cells per core x cell type for pseudobulk
NICHE_K_NEIGHBOURS <- 20
N_NICHES           <- 8

## -----------------------------------------------------------------------
## CellScape (06). ADAPT: column names differ between CellScape software
## versions/exports -- open your export and set these.
## -----------------------------------------------------------------------
CS_ID_COL        <- "CellID"
CS_X_COL         <- "X"
CS_Y_COL         <- "Y"
CS_AREA_COL      <- "Area"
CS_COORD_UNITS   <- "px"        # "px" or "um"; px are multiplied by cellscape_pixel_size_um
CS_MARKER_REGEX  <- "_Mean$"    # regex picking the per-marker mean-intensity columns
CS_MARKER_STRIP  <- "_Mean$"    # removed from column names to get the marker name
CS_DAPI_MARKER   <- "DAPI"
CS_ARCSINH_COFACTOR <- NA       # NA = per-marker cofactor (5th percentile of positive values)
CS_EXCLUDE_MARKERS <- c("DAPI") # not used for clustering

## -----------------------------------------------------------------------
## RNA <-> protein integration (08)
## -----------------------------------------------------------------------
NN_MAX_DIST_UM   <- 15     # nearest-cell transfer cutoff (serial sections: cells are NOT the same cells)
NBHD_RADIUS_UM   <- 30     # neighbourhood-averaged protein per Xenium cell
BIN_SIZE_UM      <- 50     # spatial bins for RNA-protein correlation

## -----------------------------------------------------------------------
## Zoomed QC inspection (01, section 6)
## -----------------------------------------------------------------------
## Python with tifffile + zarr (env/py_imaging.yml), used to crop the morphology
## image. ADAPT: e.g. "/home/you/miniconda3/envs/exp1-imaging/bin/python"
PYTHON_BIN          <- Sys.getenv("EXP1_PYTHON", "python3")
INSPECT_WINDOW_UM   <- 150     # side of each zoom window
INSPECT_N_PER_SLIDE <- 3       # windows chosen automatically per slide (mix of pass + fail)
## Your own windows instead (centre in um, as read off 01_qc_fail_spatial_by_slide.png):
## INSPECT_CENTERS <- data.frame(slide_id = "TMA_ORGAN", x = 5200, y = 3100)
INSPECT_CENTERS     <- NULL
## Compare with cellpose in each zoom window, same QC rules. Inside a core that 03a
## segmented with save_masks = TRUE, 03a's own result is used (data/segmentation/
## cellpose/<slide>/masks/); elsewhere cellpose runs on the window itself with the
## settings below (keep them equal to 03a's). Running cellpose needs it in
## PYTHON_BIN's environment, ~1-2 min per window on CPU. FALSE = 10x cells only.
INSPECT_CELLPOSE           <- TRUE
INSPECT_CELLPOSE_MODE      <- "nuclei_expand"   # as 03a --mode
INSPECT_CELLPOSE_EXPAND_UM <- 5                 # as 03a --expand-um

set.seed(42)
