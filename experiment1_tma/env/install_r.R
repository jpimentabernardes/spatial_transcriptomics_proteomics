## R packages for experiment1_tma (R >= 4.3). Run once on the cluster.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c(
  "Seurat", "SeuratObject",      # >= 5.0 (v5 assays, FOVs)
  "arrow",                       # cells.parquet / transcripts.parquet
  "hdf5r",                       # only if you read cell_feature_matrix.h5
  "dplyr", "tidyr", "ggplot2", "patchwork", "data.table", "jsonlite",
  "dbscan", "RANN", "sp",        # dearraying, neighbours, point-in-polygon
  "harmony",                     # slide batch correction
  "pheatmap", "Matrix"
)
bioc <- c(
  "UCell",                       # rank-based signature scores (small panels)
  "edgeR", "limma",              # pseudobulk
  "variancePartition",           # Aim 2: variance explained by organ/patient/slide
  "speckle"                      # propeller: differential composition
)
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
