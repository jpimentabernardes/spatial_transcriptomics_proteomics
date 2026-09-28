## =============================================================================
## 002_spatial_transcriptomics_prep.R
##
## Spatial transcriptomics (10x Xenium in-situ) analysis using Seurat.
##
## INPUT : Xenium RDS project already created from /outs Xenium output
## outs not included in the class due to size
## OUTPUT: data/processed/seurat_xenium.rds
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(arrow)
})

save_fig <- function(plot_obj, name, width = 9, height = 7, dpi = 150) {
  path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
  message("Saved figure -> ", path)
  invisible(path)
}
FIGURES_DIR <- "figures"


set.seed(42)
setwd('/work_beegfs/sukmb430/Spatial_Meta/data/xenium/')
OUT_DIR <- "/processed"
#ST_PATH <- file.path(OUT_DIR, "seurat_xenium.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


## -----------------------------------------------------------------------
## CONFIG
## -----------------------------------------------------------------------
## Today (public example):
XENIUM_DIR <- "output-XETG00669__0102489__Swiss_Roll_Cut4__20260827__113352/"

## -----------------------------------------------------------------------
## 1. Load the Xenium run
##    segmentations = "cell" uses the multimodal cell-boundary segmentation
##    this particular kit/dataset was generated with (as opposed to nucleus-
##    only segmentation). mols.qv.threshold = 20 is 10x's standard
##    high-confidence transcript cutoff -- don't lower this casually.
## -----------------------------------------------------------------------
xenium_cut4 <- Seurat::LoadXenium(
  XENIUM_DIR,
  fov = "fov",
  segmentations = "cell",
  molecule.coordinates = FALSE)

dim(xenium_cut4)
saveRDS(xenium_cut4, "/work_beegfs/sukmb430/Spatial_Meta/data/processed/xenium_cut4.rds")

## -----------------------------------------------------------------------
## 1. Load the Xenium run
##    segmentations = "cell" uses the multimodal cell-boundary segmentation
##    this particular kit/dataset was generated with (as opposed to nucleus-
##    only segmentation). mols.qv.threshold = 20 is 10x's standard
##    high-confidence transcript cutoff -- don't lower this casually.
## -----------------------------------------------------------------------
xenium.obj <- readRDS(ST_PATH)   
dim(xenium.obj)

## -----------------------------------------------------------------------
## 2. QC
## -----------------------------------------------------------------------
## Drop empty cells (segmented but with zero detected transcripts)
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

p_qc <- VlnPlot(
  xenium.obj,
  features = c("nCount_Xenium", "nFeature_Xenium"),
  pt.size = 0,
  ncol = 2
)
save_fig(p_qc, "Violin_plot_raw")


## ADAPT: pick QC thresholds appropriate to your gene panel size.
## lower per-cell counts than whole-transcriptome scRNA-seq -- do not reuse
## scRNA-seq QC cutoffs unmodified.
p_spatial_qc <- ImageFeaturePlot(xenium.obj, features = "nCount_Xenium", size = 0.75) +
  ggtitle("Total counts per cell, spatially")
print(p_spatial_qc)

save_fig(p_spatial_qc, "Total_cell_raw")

MIN_COUNTS   <- 10
MIN_FEATURES <- 5
xenium.obj <- subset(
  xenium.obj,
  subset = nCount_Xenium >= MIN_COUNTS & nFeature_Xenium >= MIN_FEATURES
)

dim(xenium.obj)
out_path <- file.path(OUT_DIR, "seurat_xenium.rds")

## >>> EXERCISE 1 <<<
## MIN_COUNTS/MIN_FEATURES above were chosen for this specific (small,
## custom) gene panel -- do not reuse whole-transcriptome scRNA-seq
## thresholds unmodified. Try MIN_COUNTS <- 5 and re-run from the subset()
## call above (using the ORIGINAL, unfiltered xenium.obj -- re-load from
## spamtp_transcriptomics_FULL.rds if you already overwrote it). Compare
## dim(xenium.obj) before/after: how many more cells are retained at the
## looser threshold, and does Total_cell_raw.png suggest those extra cells
## are real low-expression cells or likely background/noise?


## NOTE: Xenium is FOV/image-based (like Vizgen/CosMx), not spot-based like
## Visium -- use ImageFeaturePlot()
p_spatial_qc <- ImageFeaturePlot(xenium.obj, features = "nCount_Xenium", size = 0.75) +
  ggtitle("Total counts per cell, spatially")
print(p_spatial_qc)

save_fig(p_spatial_qc, "Total_cell_clean")


## -----------------------------------------------------------------------
## 3. Normalisation, dimensionality reduction, clustering
## -----------------------------------------------------------------------
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", clip.range = c(-10, 10))
xenium.obj <- RunPCA(xenium.obj, npcs = 30, verbose = FALSE)
xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30, verbose = FALSE)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.5, verbose = FALSE, cluster.name = "ST_clusters")
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30, verbose = FALSE)

p_umap <- DimPlot(xenium.obj, group.by = "ST_clusters", label = TRUE) + NoLegend()
p_spatial <- ImageDimPlot(xenium.obj, group.by = "ST_clusters", size = 0.75) + NoLegend()
print(p_umap | p_spatial)

save_fig(p_umap| p_spatial, "Cluster_clean")

## -----------------------------------------------------------------------
## 4. Marker genes per cluster
## -----------------------------------------------------------------------
Idents(xenium.obj) <- "ST_clusters"
markers <- FindAllMarkers(xenium.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top_markers <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 5)
print(top_markers, n = 40)

## Spatial expression of the single most distinctive marker overall, as a
## sanity check that clusters map onto sensible spatial structure (e.g. crypt
## vs. lamina propria vs. muscularis in colon tissue).
top_gene <- top_markers$gene[1]
p_marker <- ImageFeaturePlot(xenium.obj, features = top_gene, size = 0.75) +
  ggtitle(paste0("Top marker: ", top_gene))
print(p_marker)

## >>> EXERCISE 2 <<<
## Only the single most distinctive marker overall was plotted spatially
## above. Pick a DIFFERENT cluster from `top_markers` (any cluster ID other
## than the one top_gene came from) and plot its top gene the same way:
##   my_cluster_gene <- top_markers$gene[top_markers$cluster == "<cluster id>"][1]
##   ImageFeaturePlot(xenium.obj, features = my_cluster_gene, size = 0.75)
## Does the spatial pattern look like a coherent anatomical region (a crypt
## layer, an immune aggregate, an outer tissue edge), or scattered/diffuse?
## A coherent pattern is a good sign that cluster reflects real biology
## rather than a clustering artifact.

## -----------------------------------------------------------------------
## 5. Save for the integration script
## -----------------------------------------------------------------------
out_path <- file.path(OUT_DIR, "seurat_xenium_analyzed.rds")
saveRDS(xenium.obj, out_path)
