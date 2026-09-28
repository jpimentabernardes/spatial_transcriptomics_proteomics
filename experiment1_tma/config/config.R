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

## -----------------------------------------------------------------------
## TMA dearraying (02). Cores are usually 0.6-2 mm diameter with >= 200 um
## gaps; if neighbouring cores get merged, set a fixed DBSCAN_EPS_UM below the gap.
## -----------------------------------------------------------------------
DBSCAN_EPS_UM        <- NA    # NA = adaptive (3 x median kNN distance); or a fixed value in um
DBSCAN_MIN_PTS       <- 10
MIN_CELLS_PER_CORE   <- 200   # smaller DBSCAN clusters = debris / tissue fragments
CORE_MERGE_DIST_UM   <- 400   # fragments of one torn core closer than this get merged
CORE_MAX_DIST_FACTOR <- 1.25  # assign cells up to 1.25 x core radius from centre

## Core-level QC: cores failing these are flagged (not removed) in 02.
CORE_MIN_CELLS         <- 500
CORE_MIN_MEDIAN_COUNTS <- 20

## -----------------------------------------------------------------------
## Clustering / annotation (04)
## -----------------------------------------------------------------------
N_PCS            <- 30
USE_HARMONY      <- TRUE     # correct slide effects, NOT organ/patient effects
RES_LEVEL1       <- 0.3
RES_LEVEL2       <- 0.6
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

set.seed(42)
