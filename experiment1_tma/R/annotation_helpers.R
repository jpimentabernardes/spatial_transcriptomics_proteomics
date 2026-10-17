## =============================================================================
## R/annotation_helpers.R -- clustering + marker-suggested labels + editable
## label templates + figures for one annotation group. Used by
## 04_annotation_hierarchical.R and 04b_annotation_per_organ.R.
##
## Expects (set by the calling script): panel_genes, N_PCS, USE_HARMONY,
## MARKERS_MAX_CELLS, TOP_N_PER_CLUSTER, ANNOT_DIR, SEGMENTATION.
## =============================================================================

## -----------------------------------------------------------------------
## cluster a (sub)object: PCA on all panel genes, Harmony on slide, Leiden/Louvain, UMAP
cluster_cells <- function(obj, res, col, harmony = USE_HARMONY) {
  set.seed(1)
  VariableFeatures(obj) <- panel_genes
  obj <- ScaleData(obj, verbose = FALSE)
  npc <- min(N_PCS, ncol(obj) - 1, length(panel_genes) - 1)
  obj <- RunPCA(obj, npcs = npc, verbose = FALSE, seed.use = 1)
  red <- "pca"
  if (harmony && length(unique(obj$slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
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
  if (length(cl) == 1) colnames(avg) <- cl                  # one cluster: Seurat names the column "all"
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
labels_from_template <- function(summary_tab, group, prefix = "04") {
  f <- file.path(ANNOT_DIR, sprintf("%s_%s_%s.csv", prefix, SEGMENTATION, group))
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
  message("Labels for ", group, " -> ", f, "  (edit the `label` column, re-run ", prefix, ")")
  stats::setNames(summary_tab$label, summary_tab$cluster)
}

top_organs <- function(meta, col, by = "organ", n = 3) {
  meta |> filter(!is.na(.data[[by]])) |> count(cl = as.character(.data[[col]]), o = .data[[by]]) |>
    group_by(cl) |> mutate(f = n / sum(n)) |> slice_max(f, n = n, with_ties = FALSE) |>
    summarise(top = paste(sprintf("%s %.0f%%", o, 100 * f), collapse = "; "), .groups = "drop") |>
    (\(d) stats::setNames(d$top, d$cl))()
}

## one annotation group: cluster, markers, figures, template -> named label per cell
annotate_group <- function(obj, group, sigs, res, comp_by, prefix = "04", harmony = USE_HARMONY) {
  message("\n===== ", group, ": ", format(ncol(obj), big.mark = ","), " cells =====")
  col <- paste0("cl_", group)
  obj <- cluster_cells(obj, res, col, harmony = harmony)
  cl <- levels(obj[[col, drop = TRUE]])
  Idents(obj) <- col
  mk <- FindAllMarkers(obj, assay = "Xenium", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                       max.cells.per.ident = MARKERS_MAX_CELLS, verbose = FALSE)
  save_table(mk, paste0(prefix, "_", group, "_markers"))
  top <- if (nrow(mk)) mk |> group_by(cluster) |> slice_max(avg_log2FC, n = TOP_N_PER_CLUSTER) |> ungroup() else mk
  top_txt <- if (nrow(top)) tapply(top$gene, as.character(top$cluster), paste, collapse = ", ") else character()
  sug <- suggest_labels(obj, col, sigs)
  meta <- obj[[]]
  summ <- data.frame(cluster = cl, n_cells = as.integer(table(obj[[col, drop = TRUE]])[cl]),
                     top_markers = unname(top_txt[cl]), top_organs = unname(top_organs(meta, col)[cl]),
                     suggestion = sug$suggestion, margin = sug$margin)
  summ$top_markers[is.na(summ$top_markers)] <- "(no markers)"
  lab <- labels_from_template(summ, group, prefix)

  ## figures: dot plot (signature + top markers), UMAP, organ composition of clusters
  sig_genes <- intersect(unique(unlist(lapply(sigs, match_case, features = panel_genes))), panel_genes)
  feats <- unique(c(sig_genes, as.character(top$gene)))
  if (length(feats)) {
    p_dot <- DotPlot(obj, features = feats, group.by = col, assay = "Xenium") +
      scale_color_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0) +
      RotatedAxis() + theme(axis.text.x = element_text(size = 7)) +
      ggtitle(paste0(group, ": signature markers (config) + top markers per cluster"))
    save_fig(p_dot, paste0(prefix, "_", group, "_dotplot"),
             width = max(14, 0.16 * length(feats) + 4), height = max(5, 0.3 * length(cl) + 3))
  }
  obj$label_tmp <- unname(lab[as.character(obj[[col, drop = TRUE]])])
  p_umap <- plot_umap_labelled(obj, col, label_size = 5, legend = FALSE, title = paste(group, "clusters")) |
    plot_umap_labelled(obj, "label_tmp", label_size = 4, title = "current labels",
                       pal = celltype_palette(obj$label_tmp, if (group == "L1_lineage") obj$label_tmp else obj$lineage)) |
    plot_umap_labelled(obj, comp_by, label = FALSE, title = comp_by)
  save_fig(p_umap, paste0(prefix, "_", group, "_umap"), width = 30, height = 10, dpi = 160)
  comp <- meta |> filter(!is.na(.data[[comp_by]])) |> count(cl = .data[[col]], o = .data[[comp_by]]) |>
    group_by(o) |> mutate(f = n / sum(n)) |> ungroup()
  p_comp <- ggplot(comp, aes(cl, o, fill = f)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "#2a78d6", name = "fraction of\nrow's cells") +
    theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 7)) +
    labs(x = "cluster", y = comp_by, title = paste0(group, ": where does each cluster come from?"))
  save_fig(p_comp, paste0(prefix, "_", group, "_clusters_by_", comp_by),
           width = max(8, 0.35 * length(cl) + 5), height = max(5, 0.16 * length(unique(comp$o)) + 2))
  list(labels = stats::setNames(unname(lab[as.character(obj[[col, drop = TRUE]])]), colnames(obj)),
       umap = Embeddings(obj, "umap"))
}
