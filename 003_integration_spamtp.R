## =============================================================================
## 003_integration_spamtp.R
##
## Integrates the processed spatial metabolomics object (001) with the
## processed spatial transcriptomics object (002) using SpaMTP.
##
## IMPORTANT CONTEXT: your MALDI slide and the Xenium slide are ADJACENT
## (serial), not identical, tissue sections. There is no shared coordinate
## system to start with -- alignment must be done by matching tissue
## landmarks, and even after alignment you are comparing two different
## physical slices of tissue, not the same cells measured twice. Expect
## noisier correspondence than a same-section case study.
## =============================================================================

suppressPackageStartupMessages({
  library(SpaMTP)
  library(Seurat)
  library(Cardinal)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(mclust)
  library(RANN)
})

set.seed(42)
setwd('/your/own/path/')

FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(plot_obj, name, width = 9, height = 7, dpi = 150) {
  path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
  message("Saved figure -> ", path)
  invisible(path)
}

## Several steps below clustering, subset(),
## and ordinary Seurat processing on the custom SPM/SPT assays 
## drop @images at different points in this pipeline. This
## rebuilds just the centroid coordinates for whichever cells `obj` currently
## has, from a reference object known to still have a valid FOV, and writes
## directly into @images 
reattach_fov <- function(obj, reference_st, fov_name = "fov") {
  if (length(Seurat::Images(obj)) > 0) return(obj)
  obj_cells <- colnames(obj)
  orig_centroids <- reference_st[[fov_name]]@boundaries$centroids
  centroid_idx <- match(obj_cells, orig_centroids@cells)
  stopifnot("Some cells not found in reference_st's centroids" = !anyNA(centroid_idx))
  new_centroids <- orig_centroids
  new_centroids@cells <- orig_centroids@cells[centroid_idx]
  new_centroids@coords <- orig_centroids@coords[centroid_idx, , drop = FALSE]
  new_fov <- reference_st[[fov_name]]
  new_fov@boundaries$centroids <- new_centroids
  obj@images <- list()
  obj@images[[fov_name]] <- new_fov
  message("Reattached FOV '", fov_name, "' for ", length(obj_cells), " cells.")
  obj
}

## ImageDimPlot()'s `cols` argument needs one value per cluster/group level,
## even when split.by == group.by -- confirmed via "Insufficient values in
## manual scale" when passing a single color directly.
plot_clusters_spatial <- function(obj, group_by, reference_st = st, fov = "fov",
                                   size = 0.3, cols = "firebrick", ncol = 6) {
  if (length(Seurat::Images(obj)) == 0) {
    message("No images on this object -- reattaching from reference_st first.")
    obj <- reattach_fov(obj, reference_st = reference_st, fov_name = fov)
  }
  n_levels <- length(unique(obj[[group_by, drop = TRUE]]))
  ImageDimPlot(
    obj,
    fov = fov,
    group.by = group_by,
    split.by = group_by,
    size = size,
    cols = rep(cols, n_levels),
    crop = FALSE,
    dark.background = FALSE
  ) + Seurat::NoLegend() +
    ggtitle(paste0(group_by, " -- spatial location per cluster"))
}

## -----------------------------------------------------------------------
## CONFIG
## -----------------------------------------------------------------------
PROCESSED_DIR <- "data/processed"
SM_PATH <- file.path(PROCESSED_DIR, "spamtp_metabolomics.rds")
ST_PATH <- file.path(PROCESSED_DIR, "seurat_xenium.rds")
OUT_PATH <- file.path(PROCESSED_DIR, "integrated_spamtp.rds")
MAX_DIST <- 100  # microns

sm <- readRDS(SM_PATH)   # SpaMTP metabolomics object (assay "Spatial")
st <- readRDS(ST_PATH)   # Seurat Xenium object (assay "Xenium" / "SCT")

## SpaMTP's mapping functions expect current SeuratObject internals; objects
## saved from an older SeuratObject version need updating first.
if (packageVersion("SeuratObject") >= "5.2.0") {
  sm <- SeuratObject::UpdateSeuratObject(sm)
  st <- SeuratObject::UpdateSeuratObject(st)
}

## -----------------------------------------------------------------------
## 1. Visual sanity check before alignment: are we even looking at
##    comparable tissue outlines? (coordinates will NOT overlap yet)
## -----------------------------------------------------------------------
df_sm <- Seurat::GetTissueCoordinates(sm); df_sm$sample <- "MSI"
df_st <- Seurat::GetTissueCoordinates(st); df_st$sample <- "Xenium"

p_before <- (
  ggplot(df_st, aes(x, y)) + geom_point(size = 0.2, color = "lightblue4") +
    theme_classic() + ggtitle("Xenium (pre-alignment)")
) | (
  ggplot(df_sm, aes(x, y)) + geom_point(size = 0.1, color = "indianred") +
    theme_classic() + ggtitle("MALDI (pre-alignment)")
)
print(p_before)
save_fig(p_before, "03_unimodal_clusters")

## -----------------------------------------------------------------------
## 2. Alignment -- run scripts/003b_manual_alignment.R FIRST if
##    sm_transformed_manual_new.rds doesn't exist yet. AlignSpatialOmics()
##    doesn't work on Xenium/FOV objects (Visium-only @scale.factors slot) --
##    resolved there via manual landmark picking + least-squares affine fit.
## -----------------------------------------------------------------------
sm_transformed_path <- file.path(PROCESSED_DIR, "sm_transformed_manual_new.rds")
stopifnot(
  "Run scripts/003b_manual_alignment.R first" = file.exists(sm_transformed_path)
)
sm_transformed <- readRDS(sm_transformed_path)

## Confirm alignment didn't distort feature patterns (compare a known m/z
## before/after -- ADAPT the mz value to one you annotated in 001 and expect
## to see a clear spatial pattern for).
p_align_check <- SpaMTP::ImageMZPlot(sm, mzs = 861.549, size = 1) |
  SpaMTP::ImageMZPlot(sm_transformed, mzs = 861.549, size = 1)
print(p_align_check)
save_fig(p_align_check, "02_alignment_sanity_check")

## -----------------------------------------------------------------------
## 3. Manual nearest-metabolomics-pixel-per-cell mapping.
##    Assigns each Xenium cell the metabolite intensities of the closest
##    MALDI pixel within MAX_DIST microns; cells with no MALDI pixel that
##    close are dropped.
## -----------------------------------------------------------------------
sm_coords <- Seurat::GetTissueCoordinates(sm_transformed)[, c("x", "y")]
st_coords <- Seurat::GetTissueCoordinates(st)[, c("x", "y")]

nn <- RANN::nn2(data = sm_coords, query = st_coords, k = 1)
keep <- nn$nn.dists[, 1] <= MAX_DIST

sm_counts <- Seurat::GetAssayData(sm_transformed, assay = "Spatial", layer = "counts")
matched_sm_expr <- sm_counts[, nn$nn.idx[keep, 1], drop = FALSE]
kept_cells <- colnames(st)[keep]
colnames(matched_sm_expr) <- kept_cells

## Strip images off BEFORE subsetting -- avoids the heavy FOV/boundary-
## polygon recursion that overflows R's protection stack at this scale.
st_no_img <- st
st_no_img@images <- list()
st_matched <- subset(st_no_img, cells = kept_cells)
st_matched <- reattach_fov(st_matched, reference_st = st)

st_matched[["SPM"]] <- Seurat::CreateAssayObject(counts = matched_sm_expr)
st_matched[["SPT"]] <- st_matched[[DefaultAssay(st_matched)]]
combined <- st_matched
print(combined)

combined <- reattach_fov(combined, reference_st = st)
Seurat::Images(combined)

## -----------------------------------------------------------------------
## 4. Cluster each modality independently
## -----------------------------------------------------------------------
DefaultAssay(combined) <- "SPT"
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 80, reduction.name = "spt.pca", verbose = FALSE)
ElbowPlot(combined, ndims = 80, reduction = 'spt.pca')
combined <- FindNeighbors(combined, dims = 1:40, reduction = "spt.pca", verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:40, reduction = "spt.pca", reduction.name = "spt.umap", verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.3, cluster.name = "ST_clusters", verbose = FALSE)

DefaultAssay(combined) <- "SPM"
combined <- SpaMTP::NormalizeSMData(combined, assay = "SPM")
combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 80, reduction.name = "spm.pca", verbose = FALSE)
ElbowPlot(combined, ndims = 80, reduction = 'spm.pca')
combined <- FindNeighbors(combined, dims = 1:20, reduction = "spm.pca", verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:20, reduction = "spm.pca", reduction.name = "spm.umap", verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.1, cluster.name = "SM_clusters", verbose = FALSE)

## Clustering steps above have been observed to drop @images too, not just
## the mapping step -- reattach defensively before any spatial plot.
combined <- reattach_fov(combined, reference_st = st)

p_unimodal <- (
  DimPlot(combined, group.by = "ST_clusters", reduction = "spt.umap", label = TRUE) + NoLegend() +
    ggtitle("Transcriptomics-only clusters")
) | (
  DimPlot(combined, group.by = "SM_clusters", reduction = "spm.umap") + ggtitle("Metabolomics-only clusters")
)
print(p_unimodal)
save_fig(p_unimodal, "03_clusters_transcriptome_metabolome", width = 18, height = 8)

## >>> EXERCISE 1 <<<
## The number of ST_clusters and SM_clusters found above depends directly on
## the `resolution` argument to FindClusters() (0.3 for ST, 0.1 for SM).
## Try re-running just the SPT block above with resolution = 0.5 instead of
## 0.3. Does the transcriptomics UMAP split into meaningfully more clusters,
## or mostly fragment the same populations into smaller pieces? Compare
## against the spatial facet plot (plot_clusters_spatial(), section 6 below)
## to judge whether the extra clusters correspond to real, spatially
## coherent regions -- or just noise in the clustering.

## -----------------------------------------------------------------------
## 5. Cell type calling, based on transcriptomics.
##    Because Xenium does its own single-cell segmentation, we can use
##    transcriptomics-based clustering for cell typing, then pull the
##    metabolite profile associated with each annotated cell type.
## -----------------------------------------------------------------------
Idents(combined) <- "ST_clusters"
DefaultAssay(combined) <- "SPT"
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, assay = 'SPT')
top_markers <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 5)
print(top_markers, n = 40)

## Visualization: top markers per cluster + canonical cell-type markers
TOP_N_PER_CLUSTER <- 5
top_markers <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = TOP_N_PER_CLUSTER) |>
  ungroup()

canonical_markers <- c(
  "Epcam", "Mki67", "Lgr5", "Olfm4", "Stmn1", "Atoh1", "Apoa4", "Fabp1", "Alpi",
  "Slc26a3", "Car4", "Agr2", "Muc2", "Lyz1", "Defa5", "Mmp7", "Dclk1", "Trpm5",
  "Chga", "Cd3d", "Cd4", "Cd8a", "Foxp3", "Il7r", "Ccr5", "Gzmb", "Ncam1",
  "Cd79a", "Ighm", "Xcr1", "Itgax", "Itgam", "Csf1r", "Fcgr1", "Ly6c1", "Ly6g", "Siglecf"
)
all_features <- unique(c(canonical_markers, top_markers$gene))
present_features <- intersect(all_features, rownames(combined[["SPT"]]))
missing_features <- setdiff(all_features, present_features)
if (length(missing_features) > 0) {
  message("Not present in this Xenium panel, dropped from plot: ",
          paste(missing_features, collapse = ", "))
}

Idents(combined) <- "ST_clusters"
DefaultAssay(combined) <- "SPT"
p_dotplot <- DotPlot(
  combined,
  features = present_features,
  group.by = "ST_clusters"
) +
  scale_color_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))
print(p_dotplot)
save_fig(p_dotplot, "05_top_markers_dotplot_transcriptomics", width = 24, height = 8)


## -----------------------------------------------------------------------
## 6. Spatial location per cluster (see plot_clusters_spatial() above)
## -----------------------------------------------------------------------
p_st <- plot_clusters_spatial(combined, "ST_clusters")
save_fig(p_st, "06_ST_clusters_spatial_facet", width = 20, height = 18)

## -----------------------------------------------------------------------
## 7. Cell type annotation
##    Mapping from ST_clusters -> cell type, built by cross-referencing the
##    dotplot above (top DE genes per cluster) against canonical cell-type
##    markers for mouse colon. See figures/05_top_markers_dotplot_transcriptomics.png.
## -----------------------------------------------------------------------
cluster_to_celltype <- c(
  "0"  = "Colonocytes",
  "1"  = "Fibroblast",
  "2"  = "Goblet cells 1",
  "3"  = "Colonocytes Lipid 1",
  "4"  = "Colonocytes Lipid 2",
  "5"  = "Colonocytes",
  "6"  = "Colonocytes Distal",
  "7"  = "Colonocytes",
  "8"  = "Colonocytes Tryp+",
  "9"  = "Colonocytes Tryp-",
  "10" = "Stromal",
  "11" = "Stem cells",
  "12" = "Colonocytes Lipid 3",
  "13" = "Colonocytes",
  "14" = "Enteric neurons",
  "15" = "Lymphocytes",
  "16" = "Pericytes",
  "17" = "Lymphatic nodes",
  "18" = "Enteroendocrine cells",
  "19" = "Macrophages",
  "20" = "T cells",
  "21" = "Goblet cells 2",
  "22" = "Tuft cells",
  "23" = "DCs",
  "24" = "Macrophages",
  "25" = "Mesenteric lymph nodes",
  "26" = "Paneth cells"
)
combined$cell_type <- factor(unname(cluster_to_celltype[as.character(combined$ST_clusters)]))
stopifnot("Some clusters didn't get a cell_type label" = !anyNA(combined$cell_type))
table(combined$cell_type)

p_celltype <- (DimPlot(combined, group.by = "cell_type", reduction = "spt.umap", label = TRUE) + NoLegend() +
    ggtitle("Transcriptomics-only cell-type"))
print(p_celltype)
save_fig(p_celltype, "06_ST_celltype_spatial_UMAP", width = 12, height = 10)

p_st <- plot_clusters_spatial(combined, "cell_type")
save_fig(p_st, "06_ST_celltype_spatial_facet", width = 20, height = 18)

## >>> EXERCISE 2 <<<
## cluster_to_celltype above is our call, made from the dotplot + canonical
## markers -- not ground truth. Pick 2-3 clusters you find least convincing
## Look up each cluster's full top-5 marker list in
## `top_markers` (or re-inspect 05_top_markers_dotplot_transcriptomics.png)
## and either confirm the existing label or propose a better one. Does your
## revised label change how spatially coherent that cluster looks in
## plot_clusters_spatial(combined, "ST_clusters")?

combined <- reattach_fov(combined, reference_st = st)
p_celltype_spatial <- ImageDimPlot(
  combined,
  fov = "fov",
  group.by = "cell_type",
  size = 0.3,
  dark.background = FALSE
) + ggtitle("Cell-type annotation, spatial")
print(p_celltype_spatial)

save_fig(p_celltype_spatial, "08_celltype_spatial", width = 8, height = 6)

## -----------------------------------------------------------------------
## 8. Differential metabolites per cell type
##    Skip this section first if you're tight on time -- it's the most
##    computationally expensive step and not required to demonstrate the
##    integration itself.
## -----------------------------------------------------------------------
Idents(combined) <- "cell_type"
DefaultAssay(combined) <- "SPM"

sm_meta <- tryCatch(sm_transformed[["Spatial"]][[]], error = function(e) sm_transformed[["Spatial"]]@meta.data)
sm_meta$feature <- rownames(sm_meta)

annotate_and_rank <- function(dem_df) {
  dem_df$feature <- rownames(dem_df)
  dem_df |>
    dplyr::left_join(
      sm_meta[, c("feature", "all_IsomerNames", "Lipid.Maps.Category", "Lipid.Maps.Main.Class")],
      by = "feature"
    ) |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::arrange(dplyr::desc(avg_log2FC))
}


## Mesenteric lymph nodes vs. everything else
dem_tryp_neg <- FindMarkers(
  combined, ident.1 = "Colonocytes Tryp-",
  assay = "SPM", min.pct = 0.1, logfc.threshold = 0.25
) |> annotate_and_rank()

print(head(dem_tryp_neg[, c("feature", "avg_log2FC", "p_val_adj", "all_IsomerNames")], 20))

write.csv(dem_tryp_neg, "figures_dem_tryp_neg_vs_rest.csv", row.names = FALSE)


## -----------------------------------------------------------------------
## plot_metabolite_spatial(): spatial expression of a single metabolite
## (by m/z) in the integrated object. Accepts either a bare number
## (534.037798087828) or the full feature name ("mz-534.037798087828") --
## normalizes internally so you don't have to remember the exact format.
## -----------------------------------------------------------------------
plot_metabolite_spatial <- function(obj, mz, reference_st = st, fov = "fov", size = 0.5) {
  feature_name <- if (grepl("^mz-", as.character(mz))) as.character(mz) else paste0("mz-", mz)
  
  if (length(Seurat::Images(obj)) == 0) {
    message("No images on this object -- reattaching from reference_st first.")
    obj <- reattach_fov(obj, reference_st = reference_st, fov_name = fov)
  }
  
  DefaultAssay(obj) <- "SPM"
  stopifnot(
    "Feature not found in the SPM assay -- check the exact m/z value (see rownames(obj[['SPM']]))" =
      feature_name %in% rownames(obj[["SPM"]])
  )
  
  ImageFeaturePlot(obj, fov = fov, features = feature_name, size = size, dark.background = FALSE) +
    ggtitle(paste0("Spatial expression: ", feature_name))
}

p_mz_478 <- plot_metabolite_spatial(combined, "478.257135603163")
print(p_mz_478)
save_fig(p_mz_478, "09_metabolite_478_spatial")

## >>> EXERCISE 3 <<<
## Repeat the FindMarkers() + annotate_and_rank() pattern above for a
## DIFFERENT cell type of your choice (e.g. "Goblet cells 1", "Macrophages",
## "Stem cells"). What's the top metabolite? Does its LipidMaps category (if
## it has one) make biological sense for that cell type -- e.g. would you
## expect a lipid-storage signature in Colonocytes Lipid clusters, or an
## immunoglobulin/lymphoid-associated signature in Lymphocytes/DCs?

## -----------------------------------------------------------------------
## 9. Gene-metabolite correlation for Mesenteric lymph nodes
## -----------------------------------------------------------------------
Idents(combined) <- "cell_type"

## Top DE metabolite for MLN (from dem_mese, already computed)
top_mln_mz <- as.numeric(sub("^mz-", "", dem_mese$feature[1]))

combined_no_img <- combined
combined_no_img@images <- list()
combined_mln <- subset(combined_no_img, subset = cell_type == "Colonocytes Tryp-")
print(combined_mln)  # sanity check: how many MLN cells?
combined_mln <- reattach_fov(combined_mln, reference_st = st)
Seurat::Images(combined_mln)  # should now show "fov"

genes_correlated_with_top_mz <- SpaMTP::FindCorrelatedFeatures(
  combined_mln,
  mz = top_mln_mz,
  SM.assay = "SPM",
  ST.assay = "SPT",
  nfeatures = 20
)
print(genes_correlated_with_top_mz)

metabolites_correlated_with_Gdf15 <- SpaMTP::FindCorrelatedFeatures(
  combined_mln,
  gene = "Gdf15",
  SM.assay = "SPM",
  ST.assay = "SPT",
  nfeatures = 20
)
print(metabolites_correlated_with_Gdf15)

## -----------------------------------------------------------------------
## 10. Pathway enrichment
##     NOTE: this section requires RaMP-indexed metabolite annotation, which
##     the AnnotateSM() call in 001 (db = HMDB+LipidMaps+ChEBI) does not
##     produce. If this errors with "No metabolite annotation result was
##     found", see scripts/003c_pathway_reannotation.R for the (optional,
##     after-class) fix -- it needs re-annotation with db = chem_props,
##     save.intermediate = TRUE, not a parameter change here.
## -----------------------------------------------------------------------
dem_mese_genes <- FindMarkers(
  combined, ident.1 = "Mesenteric lymph nodes",
  assay = "SPT", min.pct = 0.1, logfc.threshold = 0.25
)
dem_mese_genes$feature <- rownames(dem_mese_genes)

DE_list_mln <- list(genes = dem_mese_genes, metabolites = dem_mese)

pathway_results_mln <- tryCatch({
  SpaMTP::FindRegionalPathways(
    SpaMTP = combined,
    ident = "cell_type",
    DE.list = DE_list_mln,
    analyte_types = c("genes", "metabolites"),
    SM_assay = "SPM",
    ST_assay = "SPT"
  )
}, error = function(e) {
  message("FindRegionalPathways() error: ", conditionMessage(e))
  NULL
})
print(pathway_results_mln)

## -----------------------------------------------------------------------
## 11. Save
## -----------------------------------------------------------------------
saveRDS(combined, OUT_PATH)
