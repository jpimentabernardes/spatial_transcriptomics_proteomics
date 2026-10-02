## =============================================================================
## 04_annotation_probe_based.R
##
## Probe-based, semi-manual, hierarchical cell annotation.
##
##   Level 1 (lineage): Epithelial / Immune / Fibroblast / Endothelial /
##                      Smooth muscle-pericyte / Neural-glia
##   Level 2 (within lineage): e.g. Immune -> T, CD8 T, Treg, B, plasma,
##                      macrophage, DC, mast, neutrophil ...
##
## Why hierarchical: with a few-hundred-gene panel, immune subsets are defined
## by a handful of probes each. Clustering all cells together lets abundant
## epithelium dominate the PCs; re-clustering the immune compartment on its
## own resolves the subsets that Aim 2 compares.
##
## Why "semi-manual": marker scores SUGGEST a label per cluster, written to an
## editable CSV in annotation/. You check the dotplots, edit `label`, set
## reviewed = TRUE and re-run. The decision lives in a versioned file.
##
## Marker definitions: config/marker_panel_rna.csv (HUMAN gene symbols --
## ADAPT if any cores are non-human). Markers missing from the panel are
## reported in tables/04_panel_coverage.csv -- a cell type with <2 markers on
## the panel cannot be annotated reliably, and is worth stating in Aim 1.
##
## INPUT : xenium_<seg>_cores.rds (02)
## OUTPUT: xenium_<seg>_annotated.rds, annotation/*.csv,
##         data/registration/xenium_cells_<slide>.csv (for 07/07b)
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen <- readRDS(obj_path("cores"))
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)

## -----------------------------------------------------------------------
## 1. Panel coverage of the marker definitions
## -----------------------------------------------------------------------
l1 <- restrict_to_panel(markers$level1, panel_genes, "Level 1")
l2 <- lapply(names(markers$level2), function(par) restrict_to_panel(markers$level2[[par]], panel_genes, par))
names(l2) <- names(markers$level2)
coverage <- dplyr::bind_rows(
  cbind(level = 1, parent = "", l1$coverage),
  dplyr::bind_rows(lapply(names(l2), function(par) cbind(level = 2, parent = par, l2[[par]]$coverage)))
)
save_table(coverage, "04_panel_coverage")

## -----------------------------------------------------------------------
## 2. Normalise + embed all cells
##    LogNormalize with scale.factor = median counts (not 1e4: a targeted
##    panel has ~100x fewer counts/cell than scRNA-seq). All panel genes are
##    used -- with a few hundred designed probes, "variable feature"
##    selection mostly throws away informative low-expression immune probes.
## -----------------------------------------------------------------------
DefaultAssay(xen) <- "Xenium"
xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)
VariableFeatures(xen) <- panel_genes

embed_and_cluster <- function(obj, resolution, n_pcs = N_PCS, prefix = "") {
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  red <- "pca"
  ## Harmony corrects SLIDE (technical batch). Never put organ/patient here:
  ## those are the biology Aim 1/2 asks about.
  if (USE_HARMONY && length(unique(obj$slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
    obj <- harmony::RunHarmony(obj, group.by.vars = "slide_id", reduction.use = "pca",
                               reduction.save = "harmony", verbose = FALSE)
    red <- "harmony"
  }
  obj <- FindNeighbors(obj, reduction = red, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, cluster.name = paste0(prefix, "clusters"), verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = 1:n_pcs, reduction.name = paste0(prefix, "umap"), verbose = FALSE)
  obj
}

xen <- embed_and_cluster(xen, RES_LEVEL1, prefix = "l1_")

## -----------------------------------------------------------------------
## 3. Level 1: lineage
## -----------------------------------------------------------------------
xen <- score_signatures(xen, l1$signatures)
score_cols_l1 <- paste0("score_", names(l1$signatures))
top_l1 <- top_markers_per_cluster(xen, "l1_clusters")
map_l1 <- resolve_cluster_map(
  clusters = xen$l1_clusters,
  score_df = xen[[score_cols_l1]],
  map_path = file.path(ANNOT_DIR, sprintf("level1_cluster_map_%s.csv", SEGMENTATION)),
  top_markers = top_l1
)
xen$lineage <- unname(map_l1[as.character(xen$l1_clusters)])

p_dot_l1 <- DotPlot(xen, features = unique(unlist(l1$signatures)), group.by = "l1_clusters") +
  RotatedAxis() + ggtitle("Level-1 markers per cluster (check against annotation/level1_cluster_map)")
save_fig(p_dot_l1, "04_level1_dotplot", width = 16, height = 8)

p_l1 <- (DimPlot(xen, reduction = "l1_umap", group.by = "l1_clusters", label = TRUE, raster = TRUE) + NoLegend()) |
  DimPlot(xen, reduction = "l1_umap", group.by = "lineage", raster = TRUE) |
  DimPlot(xen, reduction = "l1_umap", group.by = "organ", raster = TRUE)
save_fig(p_l1, "04_level1_umap", width = 21, height = 6)

## -----------------------------------------------------------------------
## 4. Level 2: re-cluster each lineage that has sub-types defined.
##    Built as a fresh, image-free object from the counts: faster, and
##    avoids subsetting many FOVs.
## -----------------------------------------------------------------------
xen$cell_type <- xen$lineage
counts_all <- GetAssayData(xen, assay = "Xenium", layer = "counts")

for (par in names(l2)) {
  sigs <- l2[[par]]$signatures
  cells <- colnames(xen)[xen$lineage == par]
  if (length(cells) < 200 || length(sigs) < 2) {
    message("Level 2 skipped for ", par, " (", length(cells), " cells, ", length(sigs), " signatures)")
    next
  }
  message("Level 2: ", par, " (", length(cells), " cells)")
  sub <- CreateSeuratObject(counts = counts_all[, cells], assay = "Xenium",
                            meta.data = xen[[]][cells, c("slide_id", "core_id", "organ")])
  sub <- NormalizeData(sub, scale.factor = median(sub$nCount_Xenium), verbose = FALSE)
  VariableFeatures(sub) <- panel_genes
  sub <- embed_and_cluster(sub, RES_LEVEL2, n_pcs = min(N_PCS, 20), prefix = "l2_")
  sub <- score_signatures(sub, sigs)
  top_l2 <- top_markers_per_cluster(sub, "l2_clusters")
  map_l2 <- resolve_cluster_map(
    clusters = sub$l2_clusters,
    score_df = sub[[paste0("score_", names(sigs))]],
    map_path = file.path(ANNOT_DIR, sprintf("level2_%s_cluster_map_%s.csv", par, SEGMENTATION)),
    top_markers = top_l2
  )
  xen$cell_type[cells] <- unname(map_l2[as.character(sub$l2_clusters)])
  xen[[paste0("l2_cluster_", par)]] <- NA_character_
  xen@meta.data[cells, paste0("l2_cluster_", par)] <- as.character(sub$l2_clusters)

  p_dot <- DotPlot(sub, features = unique(unlist(sigs)), group.by = "l2_clusters") + RotatedAxis() +
    ggtitle(paste0(par, ": level-2 markers per sub-cluster"))
  p_umap <- DimPlot(sub, reduction = "l2_umap", group.by = "l2_clusters", label = TRUE, raster = TRUE) + NoLegend()
  save_fig(p_umap | p_dot, paste0("04_level2_", par), width = 20, height = 7)
  rm(sub); invisible(gc())
}

## >>> CHECK <<<
## Open every annotation/*_cluster_map_*.csv. Clusters with low `margin`,
## clusters whose top_markers contradict the suggested label, and small
## clusters at core edges (cut cells) deserve a label like "Low_quality" or
## "Doublet" rather than forcing a cell type. Then re-run this script.

## -----------------------------------------------------------------------
## 5. Summary figures
## -----------------------------------------------------------------------
comp <- xen[[]] |>
  count(organ, core_id, cell_type) |>
  group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
p_comp <- ggplot(comp, aes(core_id, frac, fill = cell_type)) +
  geom_col() + facet_grid(~ organ, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  labs(y = "fraction of cells", title = "Cell-type composition per core")
save_fig(p_comp, "04_composition_per_core", width = 18, height = 7)

## Spatial gallery: one core per organ
gallery_cores <- xen[[]] |> distinct(organ, core_id) |> group_by(organ) |> slice_head(n = 1) |> pull(core_id)
p_gal <- plot_cores_spatial(xen[[]], "cell_type", cores = gallery_cores, size = 0.3) +
  ggtitle("Cell types in one core per organ")
save_fig(p_gal, "04_spatial_gallery_celltype", width = 16, height = 12)

## -----------------------------------------------------------------------
## 6. Export for registration / cross-modal alignment (07, 07b):
##    centroids, labels and log-normalised expression of the genes that have
##    a protein partner (config/rna_protein_pairs.csv).
## -----------------------------------------------------------------------
pairs <- read_config_csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pairs$gene <- as_seurat_features(pairs$gene)
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
