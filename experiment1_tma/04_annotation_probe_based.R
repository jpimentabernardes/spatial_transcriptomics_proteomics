## =============================================================================
## 04_annotation_probe_based.R
##
## Cell type calling from the transcriptomics: cluster all cells, find the
## marker genes of each cluster, look at them next to canonical markers and at
## where each cluster sits in the tissue, then YOU name the clusters.
##
##   1. Panel coverage of the marker signatures (tables/04_panel_coverage.csv)
##   2. Normalise + cluster all cells -> ST_clusters (Harmony on slide only)
##   3. FindAllMarkers per cluster + dot plot of canonical markers (config.R) and
##      the top markers of each cluster   -> figures/04_top_markers_dotplot_transcriptomics.png
##   4. Where is each cluster?             -> figures/04_ST_clusters_spatial_<slide>.png
##   5. Template for your cluster names     -> annotation/cluster_to_celltype_template_<seg>.R
##      (top markers + a lineage suggestion from the marker signatures, per cluster)
##   6. cell_type / lineage from CLUSTER_TO_CELLTYPE / CLUSTER_TO_LINEAGE (config.R);
##      before you have written them, cells are labelled cluster_<n>.
##   7. Cell types on the UMAP, in the tissue and per core.
##
## Workflow: run 04 -> look at the dot plot + spatial plots -> copy the template
## into config.R and fill in the names -> run 04 again.
##
## Gene names in config/ are matched to the panel ignoring case, so human-style
## (EPCAM) and mouse-style (Epcam) names both work.
##
## INPUT : xenium_<seg>_cores.rds (02)
## OUTPUT: xenium_<seg>_annotated.rds, annotation/cluster_to_celltype_template_<seg>.R,
##         tables/04_*.csv, data/registration/xenium_cells_<slide>.csv (for 07/07b)
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen <- readRDS(obj_path("cores"))
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)

## -----------------------------------------------------------------------
## 1. Panel coverage of the marker signatures (used only for the lineage
##    suggestion in the template). A signature with < 2 markers on the
##    panel cannot suggest anything reliably.
## -----------------------------------------------------------------------
l1 <- restrict_to_panel(markers$level1, panel_genes, "Level 1")
save_table(l1$coverage, "04_panel_coverage")
message("Level-1 signature markers on the panel: ", sum(l1$coverage$n_present), " of ", sum(l1$coverage$n_markers))

## -----------------------------------------------------------------------
## 2. Normalise + cluster all cells
##    LogNormalize with scale.factor = median counts (not 1e4: a targeted
##    panel has ~100x fewer counts/cell than scRNA-seq). All panel genes are
##    used -- with a few hundred designed probes, "variable feature"
##    selection mostly throws away informative low-expression immune probes.
## -----------------------------------------------------------------------
DefaultAssay(xen) <- "Xenium"
xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)
VariableFeatures(xen) <- panel_genes
xen <- ScaleData(xen, verbose = FALSE)
xen <- RunPCA(xen, npcs = N_PCS, verbose = FALSE)
red <- "pca"
## Harmony corrects SLIDE (technical batch). Never put organ/patient here:
## those are the biology Aim 1/2 asks about.
if (USE_HARMONY && length(unique(xen$slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
  xen <- harmony::RunHarmony(xen, group.by.vars = "slide_id", reduction.use = "pca",
                             reduction.save = "harmony", verbose = FALSE)
  red <- "harmony"
}
xen <- FindNeighbors(xen, reduction = red, dims = 1:N_PCS, verbose = FALSE)
xen <- FindClusters(xen, resolution = RES_CLUSTERS, cluster.name = "ST_clusters", verbose = FALSE)
xen <- RunUMAP(xen, reduction = red, dims = 1:N_PCS, reduction.name = "umap", verbose = FALSE)
clusters <- levels(xen$ST_clusters)
message(length(clusters), " clusters at RES_CLUSTERS = ", RES_CLUSTERS)

p_umap_cl <- (DimPlot(xen, reduction = "umap", group.by = "ST_clusters", label = TRUE, raster = TRUE) + NoLegend()) |
  DimPlot(xen, reduction = "umap", group.by = "slide_id", raster = TRUE)
save_fig(p_umap_cl, "04_ST_clusters_umap", width = 14, height = 6)

## -----------------------------------------------------------------------
## 3. Cell type calling, based on transcriptomics: marker genes per cluster,
##    shown next to canonical cell-type markers.
##    max.cells.per.ident subsamples big clusters (speed); install presto for
##    a much faster Wilcoxon test.
## -----------------------------------------------------------------------
Idents(xen) <- "ST_clusters"
cl_markers <- FindAllMarkers(xen, assay = "Xenium", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = MARKERS_MAX_CELLS, verbose = FALSE)
save_table(cl_markers, "04_cluster_markers_all")
top_markers <- cl_markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = TOP_N_PER_CLUSTER) |>
  ungroup()
save_table(top_markers, "04_cluster_markers_top")
print(top_markers |> dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2), n = 400)

canonical <- match_case(as_seurat_features(CANONICAL_MARKERS), panel_genes)
all_features <- unique(c(canonical, top_markers$gene))
present_features <- intersect(all_features, panel_genes)
missing_features <- setdiff(all_features, present_features)
if (length(missing_features) > 0) {
  message("Not present in this Xenium panel, dropped from plot: ", paste(missing_features, collapse = ", "))
}
if (length(present_features) == 0) stop("No canonical or marker genes on the panel -- check CANONICAL_MARKERS in config.R")

p_dotplot <- DotPlot(xen, features = present_features, group.by = "ST_clusters", assay = "Xenium") +
  scale_color_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Canonical markers + top markers per cluster (ST_clusters)")
save_fig(p_dotplot, "04_top_markers_dotplot_transcriptomics",
         width = max(16, 0.16 * length(present_features) + 4), height = max(6, 0.3 * length(clusters) + 3))

## -----------------------------------------------------------------------
## 4. Spatial location per cluster, one figure per slide
## -----------------------------------------------------------------------
for (sid in unique(xen$slide_id)) {
  save_fig(plot_clusters_spatial(xen[[]], "ST_clusters", slide = sid),
           paste0("04_ST_clusters_spatial_", sid), width = 20, height = 18)
}

## -----------------------------------------------------------------------
## 5. Template for the cluster names: top markers and a lineage suggestion
##    (level-1 signature of config/marker_panel_rna.csv with the highest
##    mean score; margin = best - second best, small = ambiguous).
## -----------------------------------------------------------------------
top_by_cluster <- top_markers |> group_by(cluster) |> summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
top_by_cluster <- stats::setNames(top_by_cluster$genes, as.character(top_by_cluster$cluster))[clusters]
top_by_cluster[is.na(top_by_cluster)] <- "(no markers)"
suggested <- stats::setNames(rep("Unknown", length(clusters)), clusters)
margin <- stats::setNames(rep(NA_real_, length(clusters)), clusters)
usable <- l1$signatures[lengths(l1$signatures) >= 2]
if (length(usable) >= 2) {
  sc <- score_signatures(xen, usable)[[paste0("score_", names(usable))]]
  agg <- stats::aggregate(sc, by = list(cluster = as.character(xen$ST_clusters)), FUN = mean)
  m <- as.matrix(agg[, -1, drop = FALSE]); dimnames(m) <- list(agg$cluster, names(usable))
  m <- m[clusters, , drop = FALSE]
  suggested[] <- colnames(m)[apply(m, 1, which.max)]
  srt <- t(apply(m, 1, sort, decreasing = TRUE))
  margin[] <- round(srt[, 1] - srt[, 2], 3)
} else {
  message("Fewer than 2 level-1 signatures with >= 2 markers on the panel -- no lineage suggestion.")
}
summary_tab <- data.frame(cluster = clusters, n_cells = as.integer(table(xen$ST_clusters)[clusters]),
                          top_markers = unname(top_by_cluster), suggested_lineage = unname(suggested),
                          lineage_margin = unname(margin))
save_table(summary_tab, "04_cluster_summary")

tmpl_path <- file.path(ANNOT_DIR, sprintf("cluster_to_celltype_template_%s.R", SEGMENTATION))
n <- length(clusters); comma <- ifelse(seq_len(n) < n, ",", " ")
writeLines(c(
  sprintf("## %d clusters (SEGMENTATION = %s, RES_CLUSTERS = %s). Copy into config/config.R,",
          n, SEGMENTATION, RES_CLUSTERS),
  "## fill in the cell type names (see figures/04_top_markers_dotplot_transcriptomics.png),",
  "## check the lineages, and re-run 04. Comments: top markers of the cluster | lineage suggestion margin.",
  "CLUSTER_TO_CELLTYPE <- c(",
  sprintf('  "%s" = ""%s   # %s', clusters, comma, top_by_cluster),
  ")",
  "CLUSTER_TO_LINEAGE <- c(",
  sprintf('  "%s" = "%s"%s   # margin %s', clusters, suggested, comma, margin),
  ")"), tmpl_path)
message("Template for your cluster names -> ", tmpl_path)

## -----------------------------------------------------------------------
## 6. cell_type and lineage from your mapping (config.R)
## -----------------------------------------------------------------------
apply_map <- function(map, what) {
  unlabelled <- setdiff(clusters, names(map)[nzchar(map)])
  if (length(unlabelled) > 0) {
    stop(what, " has no (or an empty) name for cluster(s): ", paste(unlabelled, collapse = ", "),
         "\n  The clusters changed (data / SEGMENTATION / RES_CLUSTERS?) or the template is not filled in yet.",
         "\n  See ", tmpl_path)
  }
  unname(map[as.character(xen$ST_clusters)])
}
if (is.null(CLUSTER_TO_CELLTYPE)) {
  xen$cell_type <- paste0("cluster_", xen$ST_clusters)
  message("CLUSTER_TO_CELLTYPE not set yet -- cells labelled cluster_<n>. Fill in ", tmpl_path)
} else {
  xen$cell_type <- apply_map(CLUSTER_TO_CELLTYPE, "CLUSTER_TO_CELLTYPE")
}
xen$lineage <- if (is.null(CLUSTER_TO_LINEAGE)) unname(suggested[as.character(xen$ST_clusters)]) else
  apply_map(CLUSTER_TO_LINEAGE, "CLUSTER_TO_LINEAGE")
bad_lin <- setdiff(unique(xen$lineage), c(names(markers$level1), "Unknown", "Low_quality"))
if (length(bad_lin) > 0) {
  message("NOTE: lineage values not among the level-1 names of ", basename(MARKERS_RNA_CSV), ": ",
          paste(bad_lin, collapse = ", "), " -- 05 selects immune cells with lineage == \"Immune\".")
}
print(table(cell_type = xen$cell_type, lineage = xen$lineage))

## -----------------------------------------------------------------------
## 7. Cell types on the UMAP, in the tissue and per core
## -----------------------------------------------------------------------
p_celltype <- DimPlot(xen, group.by = "cell_type", reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
  NoLegend() + ggtitle("Transcriptomics-only cell type")
save_fig(p_celltype, "04_celltype_UMAP", width = 12, height = 10)
for (sid in unique(xen$slide_id)) {
  save_fig(plot_clusters_spatial(xen[[]], "cell_type", slide = sid),
           paste0("04_celltype_spatial_facet_", sid), width = 20, height = 18)
}

comp <- xen[[]] |>
  count(organ, core_id, cell_type) |>
  group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
p_comp <- ggplot(comp, aes(core_id, frac, fill = cell_type)) +
  geom_col() + facet_grid(~ organ, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  labs(y = "fraction of cells", title = "Cell-type composition per core")
save_fig(p_comp, "04_composition_per_core", width = 18, height = 7)

## Spatial gallery: a few cores
gallery_cores <- head(unique(xen$core_id[!is.na(xen$core_id)]), 6)
p_gal <- plot_cores_spatial(xen[[]], "cell_type", cores = gallery_cores, size = 0.5) +
  ggtitle("Cell types in a few cores")
save_fig(p_gal, "04_spatial_gallery_celltype", width = 16, height = 12)

## -----------------------------------------------------------------------
## 6. Export for registration / cross-modal alignment (07, 07b):
##    centroids, labels and log-normalised expression of the genes that have
##    a protein partner (config/rna_protein_pairs.csv).
## -----------------------------------------------------------------------
pairs <- read_config_csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pairs$gene <- match_case(as_seurat_features(pairs$gene), panel_genes)
pair_genes <- intersect(pairs$gene, panel_genes)
expr <- GetAssayData(xen, assay = "Xenium", layer = "data")[pair_genes, , drop = FALSE]
for (sid in unique(xen$slide_id)) {
  cells <- colnames(xen)[xen$slide_id == sid]
  out <- data.frame(
    cell = cells,
    slide_id = sid,
    core_id = xen$core_id[cells],
    x_um = xen$x_um[cells],
    y_um = xen$y_um[cells],
    lineage = xen$lineage[cells],
    cell_type = xen$cell_type[cells]
  )
  out <- cbind(out, as.data.frame(as.matrix(Matrix::t(expr[, cells, drop = FALSE]))))
  f <- file.path(REG_DIR, paste0("xenium_cells_", sid, ".csv.gz"))
  data.table_available <- requireNamespace("data.table", quietly = TRUE)
  if (data.table_available) data.table::fwrite(out, f) else utils::write.csv(out, gzfile(f), row.names = FALSE)
  message("Exported ", nrow(out), " cells -> ", f)
}

saveRDS(xen, obj_path("annotated"))
message("Saved -> ", obj_path("annotated"))
