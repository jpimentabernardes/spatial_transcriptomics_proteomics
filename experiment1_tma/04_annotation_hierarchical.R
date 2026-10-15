## =============================================================================
## 04_annotation_hierarchical.R
##
## Cell type calling for multi-organ TMAs (TMA_ORGAN + TMA_TUMOR together), in
## levels, each clustered on its own and named by YOU:
##
##   L1  lineage, all cells of both TMAs     -> Immune / Epithelial / Mesenchymal /
##                                              Endothelial / Neural_NE (/ Low_quality)
##   L2  Immune, all organs together          -> one consistent set of immune labels
##       Stroma (Mesenchymal + Endothelial), all organs together
##       Epithelial + Neural_NE, per FAMILY (config/annotation_families.csv: organ or
##       diagnosis -> family, normal and tumour together, e.g. GI_glandular, Squamous)
##
## Per group: clustering (Harmony on slide only), FindAllMarkers, a dot plot of the
## group's signature markers (config/marker_panel_rna.csv) + top markers, UMAP,
## which organs each cluster comes from, and a template you fill in:
##   annotation/04_<seg>_<group>.csv   cluster, n_cells, top_markers, top_organs,
##                                     suggestion, margin, label   <- edit `label`
## The label column is pre-filled with the marker-based suggestion; correct it and
## re-run 04. Your labels are kept as long as the clustering is unchanged (same
## clusters with the same cell numbers); if it changed, the old file is kept as
## *.previous.csv and a fresh template is written.
##
## Final columns: lineage (L1), family, annotation_group, cell_type.
## Epithelial cells of cores without a family (lymphoid organs, muscle, fat, ...)
## and families with < MIN_CELLS_GROUP cells get "<lineage>_<family>" labels.
##
## INPUT : xenium_<seg>_cores.rds (02), config/marker_panel_rna.csv,
##         config/annotation_families.csv (or a `family` column in tma_map.csv)
## OUTPUT: xenium_<seg>_annotated.rds, annotation/04_<seg>_*.csv, tables/04_*.csv,
##         figures/04_*.png, data/registration/xenium_cells_<slide>.csv.gz
## =============================================================================

source("config/config.R")
source("R/helpers.R")

ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))
RES_LINEAGE    <- get0("RES_LINEAGE",    ifnotfound = 0.3)
RES_IMMUNE     <- get0("RES_IMMUNE",     ifnotfound = 0.8)
RES_STROMA     <- get0("RES_STROMA",     ifnotfound = 0.5)
RES_EPITHELIAL <- get0("RES_EPITHELIAL", ifnotfound = 0.6)
MIN_CELLS_GROUP <- get0("MIN_CELLS_GROUP", ifnotfound = 500)
LINEAGES <- c("Immune", "Epithelial", "Mesenchymal", "Endothelial", "Neural_NE")

xen <- readRDS(obj_path("cores"))
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)
DefaultAssay(xen) <- "Xenium"

## -----------------------------------------------------------------------
## 0. Family per cell (from the core's organ / diagnosis)
## -----------------------------------------------------------------------
key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
if (!"family" %in% colnames(xen[[]]) || all(is.na(xen$family))) {
  fam <- read_config_csv(ANNOT_FAMILIES_CSV, stringsAsFactors = FALSE, fileEncoding = "latin1")
  xen$family <- fam$family[match(key(xen$organ), key(fam$organ))]
  unmapped <- unique(xen$organ[is.na(xen$family) & !is.na(xen$organ)])
  if (length(unmapped)) message("Organs without a family in ", basename(ANNOT_FAMILIES_CSV), " (-> 'none'): ",
                                paste(unmapped, collapse = "; "))
}
xen$family[is.na(xen$family) | !nzchar(xen$family)] <- "none"
print(table(family = xen$family, slide = xen$slide_id))

## -----------------------------------------------------------------------
## helpers
## -----------------------------------------------------------------------
## cluster a (sub)object: PCA on all panel genes, Harmony on slide, Leiden/Louvain, UMAP
cluster_cells <- function(obj, res, col) {
  set.seed(1)
  VariableFeatures(obj) <- panel_genes
  obj <- ScaleData(obj, verbose = FALSE)
  npc <- min(N_PCS, ncol(obj) - 1, length(panel_genes) - 1)
  obj <- RunPCA(obj, npcs = npc, verbose = FALSE, seed.use = 1)
  red <- "pca"
  if (USE_HARMONY && length(unique(obj$slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
    set.seed(1)
    obj <- harmony::RunHarmony(obj, group.by.vars = "slide_id", reduction.use = "pca",
                               reduction.save = "harmony", verbose = FALSE)
    red <- "harmony"
  }
  obj <- FindNeighbors(obj, reduction = red, dims = 1:npc, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, cluster.name = col, random.seed = 0, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = 1:npc, reduction.name = "umap", seed.use = 42, verbose = FALSE)
  obj
}

## marker-based suggestion per cluster: mean of the scaled average expression of
## each signature's genes (fast on millions of cells); margin = best - second best
suggest_labels <- function(obj, col, sigs) {
  cl <- levels(obj[[col, drop = TRUE]])
  sigs <- lapply(sigs, function(g) intersect(match_case(g, panel_genes), panel_genes))
  sigs <- sigs[lengths(sigs) > 0]
  if (length(sigs) < 1) return(data.frame(cluster = cl, suggestion = "Unknown", margin = NA))
  genes <- unique(unlist(sigs))
  avg <- AverageExpression(obj, assays = "Xenium", features = genes, group.by = col, layer = "data")$Xenium
  avg <- as.matrix(avg)
  colnames(avg) <- sub("^g", "", colnames(avg))            # Seurat prefixes numeric idents with "g"
  z <- t(scale(t(log1p(avg))))
  z[is.na(z)] <- 0
  sc <- sapply(sigs, function(g) colMeans(z[g, , drop = FALSE]))
  sc <- matrix(sc, ncol = length(sigs), dimnames = list(colnames(avg), names(sigs)))[cl, , drop = FALSE]
  best <- colnames(sc)[apply(sc, 1, which.max)]
  srt <- t(apply(sc, 1, sort, decreasing = TRUE))
  data.frame(cluster = cl, suggestion = best,
             margin = if (ncol(sc) > 1) round(srt[, 1] - srt[, 2], 2) else NA)
}

## template: write / read back your labels (kept while the clustering is unchanged)
labels_from_template <- function(summary_tab, group) {
  f <- file.path(ANNOT_DIR, sprintf("04_%s_%s.csv", SEGMENTATION, group))
  summary_tab$label <- summary_tab$suggestion
  if (file.exists(f)) {
    old <- read_config_csv(f, stringsAsFactors = FALSE, colClasses = c(cluster = "character"))
    same <- nrow(old) == nrow(summary_tab) && setequal(old$cluster, summary_tab$cluster) &&
      all(old$n_cells[match(summary_tab$cluster, old$cluster)] == summary_tab$n_cells)
    if (same) {
      lab <- old$label[match(summary_tab$cluster, old$cluster)]
      summary_tab$label <- ifelse(is.na(lab) | !nzchar(trimws(lab)), summary_tab$suggestion, trimws(lab))
    } else {
      file.copy(f, sub("\\.csv$", ".previous.csv", f), overwrite = TRUE)
      message("!! ", group, ": the clustering changed -- your old labels are in ",
              basename(sub("\\.csv$", ".previous.csv", f)), "; a new template was written.")
    }
  }
  utils::write.csv(summary_tab, f, row.names = FALSE)
  message("Labels for ", group, " -> ", f, "  (edit the `label` column, re-run 04)")
  stats::setNames(summary_tab$label, summary_tab$cluster)
}

top_organs <- function(meta, col, by = "organ", n = 3) {
  meta |> filter(!is.na(.data[[by]])) |> count(cl = as.character(.data[[col]]), o = .data[[by]]) |>
    group_by(cl) |> mutate(f = n / sum(n)) |> slice_max(f, n = n, with_ties = FALSE) |>
    summarise(top = paste(sprintf("%s %.0f%%", o, 100 * f), collapse = "; "), .groups = "drop") |>
    (\(d) stats::setNames(d$top, d$cl))()
}

## one annotation group: cluster, markers, figures, template -> named label per cell
annotate_group <- function(obj, group, sigs, res, comp_by) {
  message("\n===== ", group, ": ", format(ncol(obj), big.mark = ","), " cells =====")
  col <- paste0("cl_", group)
  obj <- cluster_cells(obj, res, col)
  cl <- levels(obj[[col, drop = TRUE]])
  Idents(obj) <- col
  mk <- FindAllMarkers(obj, assay = "Xenium", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                       max.cells.per.ident = MARKERS_MAX_CELLS, verbose = FALSE)
  save_table(mk, paste0("04_", group, "_markers"))
  top <- if (nrow(mk)) mk |> group_by(cluster) |> slice_max(avg_log2FC, n = TOP_N_PER_CLUSTER) |> ungroup() else mk
  top_txt <- if (nrow(top)) tapply(top$gene, as.character(top$cluster), paste, collapse = ", ") else character()
  sug <- suggest_labels(obj, col, sigs)
  meta <- obj[[]]
  summ <- data.frame(cluster = cl, n_cells = as.integer(table(obj[[col, drop = TRUE]])[cl]),
                     top_markers = unname(top_txt[cl]), top_organs = unname(top_organs(meta, col)[cl]),
                     suggestion = sug$suggestion, margin = sug$margin)
  summ$top_markers[is.na(summ$top_markers)] <- "(no markers)"
  lab <- labels_from_template(summ, group)

  ## figures: dot plot (signature + top markers), UMAP, organ composition of clusters
  sig_genes <- intersect(unique(unlist(lapply(sigs, match_case, features = panel_genes))), panel_genes)
  feats <- unique(c(sig_genes, as.character(top$gene)))
  if (length(feats)) {
    p_dot <- DotPlot(obj, features = feats, group.by = col, assay = "Xenium") +
      scale_color_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0) +
      RotatedAxis() + theme(axis.text.x = element_text(size = 7)) +
      ggtitle(paste0(group, ": signature markers (config) + top markers per cluster"))
    save_fig(p_dot, paste0("04_", group, "_dotplot"),
             width = max(14, 0.16 * length(feats) + 4), height = max(5, 0.3 * length(cl) + 3))
  }
  obj$label_tmp <- unname(lab[as.character(obj[[col, drop = TRUE]])])
  p_umap <- plot_umap_labelled(obj, col, label_size = 5, legend = FALSE, title = paste(group, "clusters")) |
    plot_umap_labelled(obj, "label_tmp", label_size = 4, title = "current labels") |
    plot_umap_labelled(obj, comp_by, label = FALSE, title = comp_by)
  save_fig(p_umap, paste0("04_", group, "_umap"), width = 30, height = 10, dpi = 160)
  comp <- meta |> filter(!is.na(.data[[comp_by]])) |> count(cl = .data[[col]], o = .data[[comp_by]]) |>
    group_by(o) |> mutate(f = n / sum(n)) |> ungroup()
  p_comp <- ggplot(comp, aes(cl, o, fill = f)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "#2a78d6", name = "fraction of\nrow's cells") +
    theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 7)) +
    labs(x = "cluster", y = comp_by, title = paste0(group, ": where does each cluster come from?"))
  save_fig(p_comp, paste0("04_", group, "_clusters_by_", comp_by),
           width = max(8, 0.35 * length(cl) + 5), height = max(5, 0.16 * length(unique(comp$o)) + 2))
  list(labels = stats::setNames(unname(lab[as.character(obj[[col, drop = TRUE]])]), colnames(obj)),
       umap = Embeddings(obj, "umap"))
}

## -----------------------------------------------------------------------
## 1. Normalise all cells once (log-normalised data reused by every level)
## -----------------------------------------------------------------------
cov <- restrict_to_panel(markers$level1, panel_genes, "Level 1")
save_table(cov$coverage, "04_panel_coverage_level1")
xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)

## -----------------------------------------------------------------------
## 2. L1: lineage on all cells
## -----------------------------------------------------------------------
lin <- annotate_group(xen, "L1_lineage", markers$level1, RES_LINEAGE, comp_by = "family")
xen$lineage <- unname(lin$labels[colnames(xen)])
xen[["umap"]] <- CreateDimReducObject(lin$umap[colnames(xen), ], key = "UMAP_", assay = "Xenium")
rm(lin); gc()
bad <- setdiff(unique(xen$lineage), c(LINEAGES, "Low_quality", "Unknown"))
if (length(bad)) message("NOTE: lineage labels outside ", paste(LINEAGES, collapse = "/"), "/Low_quality: ",
                         paste(bad, collapse = ", "), " -- these cells keep the lineage as cell type.")
print(table(lineage = xen$lineage, slide = xen$slide_id))

## -----------------------------------------------------------------------
## 3. L2 groups
## -----------------------------------------------------------------------
xen$annotation_group <- dplyr::case_when(
  xen$lineage == "Immune" ~ "Immune",
  xen$lineage %in% c("Mesenchymal", "Endothelial") ~ "Stroma",
  xen$lineage %in% c("Epithelial", "Neural_NE") & xen$family != "none" ~ paste0("Epi_", xen$family),
  TRUE ~ "none")
xen$cell_type <- ifelse(xen$lineage %in% c("Epithelial", "Neural_NE"),
                        paste0(xen$lineage, "_", ifelse(xen$family == "none", "other", xen$family)),
                        xen$lineage)
groups <- sort(setdiff(unique(xen$annotation_group), "none"))
group_tab <- as.data.frame(table(annotation_group = xen$annotation_group))
save_table(group_tab, "04_annotation_groups")
for (g in groups) {
  cells <- colnames(xen)[xen$annotation_group == g]
  if (length(cells) < MIN_CELLS_GROUP) {
    message(g, ": only ", length(cells), " cells (< MIN_CELLS_GROUP = ", MIN_CELLS_GROUP, ") -- not sub-clustered")
    next
  }
  if (g == "Immune") {
    sigs <- markers$level2[["Immune"]]; res <- RES_IMMUNE; by <- "family"
  } else if (g == "Stroma") {
    sigs <- markers$level2[["Stroma"]]; res <- RES_STROMA; by <- "family"
  } else {
    f <- sub("^Epi_", "", g)
    sigs <- c(markers$level2[[f]], markers$level2[["Immune"]]["Proliferating"])
    if (is.null(markers$level2[[f]])) message(g, ": no level-2 markers for family '", f, "' in ", basename(MARKERS_RNA_CSV))
    res <- RES_EPITHELIAL; by <- "organ"
  }
  sub <- subset(xen, cells = cells)
  lab <- annotate_group(sub, g, sigs, res, comp_by = by)$labels
  xen$cell_type[match(names(lab), colnames(xen))] <- lab
  rm(sub); gc()
}
xen$celltype <- factor(xen$cell_type)
ct_tab <- xen[[]] |> count(lineage, annotation_group, cell_type, name = "n_cells")
save_table(ct_tab, "04_cell_types")
print(as.data.frame(ct_tab), row.names = FALSE)

## -----------------------------------------------------------------------
## 4. Whole-dataset figures with the final labels
## -----------------------------------------------------------------------
save_fig(plot_umap_labelled(xen, "lineage", label_size = 6, title = "Lineage (L1)"), "04_lineage_umap",
         width = 14, height = 11, dpi = 180)
n_slides <- length(unique(xen$slide_id))
lin_meta <- xen[[]]; lin_meta$cell_type <- lin_meta$lineage
save_fig(plot_celltype_sections(lin_meta, panel_in = 11), "04_lineage_sections",
         width = 11 * n_slides + 4, height = 12, dpi = 220)

## immune cell types per core, cores grouped by family (the Aim 1/2 overview)
imm <- xen[[]][xen$lineage == "Immune" & !is.na(xen$core_id), ]
if (nrow(imm)) {
  pal_imm <- celltype_palette(imm$cell_type, imm$lineage)
  comp <- imm |> count(family, core_id, cell_type) |> group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
  comp$cell_type <- factor(comp$cell_type, levels = names(pal_imm))
  p_imm <- ggplot(comp, aes(core_id, frac, fill = cell_type)) + geom_col(width = 0.9) +
    facet_grid(~ family, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_imm, name = "Immune cell type") + scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 8)) +
    labs(x = "cores", y = "fraction of immune cells", title = "Immune composition per core, by family")
  save_fig(p_imm, "04_immune_composition_per_core", width = 26, height = 8, dpi = 180)
}
lin_comp <- xen[[]] |> filter(!is.na(core_id)) |> count(family, core_id, lineage) |>
  group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
p_lin <- ggplot(lin_comp, aes(core_id, frac, fill = lineage)) + geom_col(width = 0.9) +
  facet_grid(~ family, scales = "free_x", space = "free_x") + scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 8)) +
  labs(x = "cores", y = "fraction of cells", title = "Lineage composition per core, by family")
save_fig(p_lin, "04_lineage_composition_per_core", width = 26, height = 7, dpi = 180)

## -----------------------------------------------------------------------
## 5. Export for registration / cross-modal alignment (07, 07b)
## -----------------------------------------------------------------------
pairs <- read_config_csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pair_genes <- intersect(match_case(as_seurat_features(pairs$gene), panel_genes), panel_genes)
expr <- GetAssayData(xen, assay = "Xenium", layer = "data")[pair_genes, , drop = FALSE]
for (sid in unique(xen$slide_id)) {
  cells <- colnames(xen)[xen$slide_id == sid]
  out <- data.frame(cell = cells, slide_id = sid, core_id = xen$core_id[cells],
                    x_um = xen$x_um[cells], y_um = xen$y_um[cells], family = xen$family[cells],
                    lineage = xen$lineage[cells], cell_type = xen$cell_type[cells])
  if (length(pair_genes)) out <- cbind(out, as.data.frame(as.matrix(Matrix::t(expr[, cells, drop = FALSE]))))
  f <- file.path(REG_DIR, paste0("xenium_cells_", sid, ".csv.gz"))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(out, f) else
    utils::write.csv(out, gzfile(f), row.names = FALSE)
  message("Exported ", nrow(out), " cells -> ", f)
}

saveRDS(xen, obj_path("annotated"))
message("Saved -> ", obj_path("annotated"))
