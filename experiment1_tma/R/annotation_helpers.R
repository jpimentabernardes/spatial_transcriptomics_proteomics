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
  if ("is_tumour" %in% colnames(meta) && length(unique(meta$is_tumour)) > 1) {
    ## which clusters are shared by normal and tumour cores, which are tumour-only
    pt <- tapply(meta$is_tumour, as.character(meta[[col]]), mean)
    summ$pct_cells_from_tumour_cores <- round(100 * as.numeric(pt[cl]), 1)
  }
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

## =============================================================================
## Per-organ ("site") annotation, shared by 04b (epithelium) and 04c (all lineages)
## =============================================================================

## site (config/annotation_families.csv) and tumour status (slide level) per cell
add_site_and_tumour <- function(xen, families_csv) {
  key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
  fam <- read_config_csv(families_csv, stringsAsFactors = FALSE, fileEncoding = "latin1")
  if (!"site" %in% colnames(fam)) stop(basename(families_csv), " has no `site` column -- update config/.")
  site <- fam$site[match(key(xen$organ), key(fam$organ))]
  site[!is.na(site) & !nzchar(trimws(site))] <- NA
  xen$site <- site
  tum <- tumour_slides()
  xen$is_tumour <- xen$slide_id %in% tum
  message("Tumour slide(s): ", paste(tum, collapse = ", "), "  (from config/slides.csv)")
  print(table(slide = xen$slide_id, tumour = xen$is_tumour))
  xen
}

## sites to run: every site in the csv, or ORGAN_SITES
pick_sites <- function(xen, organ_sites = NULL) {
  sites <- sort(unique(stats::na.omit(xen$site)))
  if (!is.null(organ_sites)) {
    bad <- setdiff(organ_sites, sites)
    if (length(bad)) message("ORGAN_SITES not found in the data: ", paste(bad, collapse = ", "))
    sites <- intersect(sites, organ_sites)
  }
  if (!length(sites)) stop("No sites to run -- fill the `site` column of annotation_families.csv")
  sites
}

## UMAP of one group coloured by core (+ by normal / tumour and by label): do clusters
## follow cell types (mixed cores) or single cores / patients (batch, or a tumour)?
plot_umap_by_core <- function(emb, meta, lab, g, prefix = "04b", pal_lineage = "Epithelial", max_cells = 2e5) {
  d <- data.frame(emb[, 1:2], core = meta[rownames(emb), "core_id"], organ = meta[rownames(emb), "organ"],
                  tissue = ifelse(meta[rownames(emb), "is_tumour"], "tumour", "normal"),
                  label = unname(lab[rownames(emb)]))
  colnames(d)[1:2] <- c("u1", "u2")
  d <- d[!is.na(d$core), ]
  if (!nrow(d)) return(invisible(NULL))
  if (nrow(d) > max_cells) d <- d[sample(nrow(d), max_cells), ]
  d <- d[sample(nrow(d)), ]                                   # no core drawn on top of all others
  d$core_lab <- paste0(d$core, " (", d$organ, ")")
  n_core <- length(unique(d$core_lab))
  pal <- stats::setNames(grDevices::hcl(h = seq(15, 375, length.out = n_core + 1)[seq_len(n_core)],
                                        c = rep(c(85, 55), length.out = n_core),
                                        l = rep(c(55, 72, 42), length.out = n_core)),
                         sort(unique(d$core_lab)))
  pt <- if (nrow(d) > 1e5) 0.25 else if (nrow(d) > 2e4) 0.5 else if (nrow(d) > 2000) 1 else 2
  base <- list(theme_void(base_size = 12), coord_fixed(), theme(plot.title = element_text(size = 13)))
  p_core <- ggplot(d, aes(u1, u2, colour = core_lab)) + geom_point(size = pt, stroke = 0) +
    scale_colour_manual(values = pal, name = "core (organ)") + base +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = if (n_core > 25) 2 else 1)) +
    ggtitle(paste0(g, ": coloured by core (", n_core, " cores)"))
  plots <- list(p_core)
  if (length(unique(d$tissue)) > 1)
    plots <- c(plots, list(ggplot(d, aes(u1, u2, colour = tissue)) + geom_point(size = pt, stroke = 0) +
      scale_colour_manual(values = c(normal = "#3B6FB6", tumour = "#C0392B"), name = NULL) + base +
      guides(colour = guide_legend(override.aes = list(size = 3))) + ggtitle("normal vs tumour cores")))
  pal_l <- celltype_palette(d$label, rep(pal_lineage, nrow(d)))
  plots <- c(plots, list(ggplot(d, aes(u1, u2, colour = label)) + geom_point(size = pt, stroke = 0) +
    scale_colour_manual(values = pal_l, name = NULL) + base +
    guides(colour = guide_legend(override.aes = list(size = 3))) + ggtitle("current labels")))
  save_fig(patchwork::wrap_plots(plots, nrow = 1), paste0(prefix, "_", g, "_umap_by_core"),
           width = 9.5 * length(plots) + 2, height = 9, dpi = 160)
}

## which labels make up each core, normal vs tumour (one joint clustering)
plot_site_composition <- function(meta, lab, site, g, prefix, what, pal_lineage = "Epithelial") {
  m <- meta[names(lab), ]
  m$label <- lab
  comp <- m |> filter(!is.na(core_id)) |> count(is_tumour, organ, core_id, label) |>
    group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
  if (!nrow(comp)) return(invisible(NULL))
  comp$panel <- ifelse(comp$is_tumour, "tumour", "normal")
  pal <- celltype_palette(comp$label, rep(pal_lineage, nrow(comp)))
  p <- ggplot(comp, aes(core_id, frac, fill = factor(label, levels = names(pal)))) +
    geom_col(width = 0.9) +
    facet_grid(~ panel + organ, scales = "free_x", space = "free_x", labeller = label_wrap_gen(18)) +
    scale_fill_manual(values = pal, name = "label") + scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5), panel.grid = element_blank(),
          strip.text = element_text(size = 8)) +
    labs(x = "core", y = paste("fraction of", what, "cells"),
         title = paste0(site, ": ", what, " cell types per core, normal vs tumour (one joint clustering)"))
  save_fig(p, paste0(prefix, "_", g, "_composition_per_core"),
           width = max(10, 0.25 * length(unique(comp$core_id)) + 6), height = 6.5)
}

## Annotate one lineage group per site, in each mode:
##   normal       normal cores of the site      -> group Organ_<site><suffix>
##   with_tumour  normal + tumour cores         -> group OrganTumour_<site><suffix>
## sig_fun(cells) returns the signature list for those cells. Returns a list with
## one named label vector per mode (names = cell names) and the failed groups.
annotate_sites <- function(xen, sites, lineages, sig_fun, res, prefix, suffix = "", modes = c("normal", "with_tumour"),
                           min_cells = 300, what = paste(lineages, collapse = " + "), pal_lineage = "Epithelial") {
  out <- list(normal = character(), with_tumour = character(), failed = list())
  in_lin <- xen$lineage %in% lineages
  for (mode in modes) {
    for (s in sites) {
      sel <- in_lin & xen$site %in% s & (mode != "normal" | !xen$is_tumour)
      cells <- colnames(xen)[sel]
      g <- paste0(if (mode == "normal") "Organ_" else "OrganTumour_", s, suffix)
      if (mode != "normal" && !any(xen$is_tumour[sel])) {
        message(g, ": no tumour cores for this site -- skipped (same as Organ_", s, suffix, ")")
        next
      }
      if (mode == "normal" && !length(cells) && any(xen$is_tumour[in_lin & xen$site %in% s])) {
        message(g, ": tumour cores only (no normal tissue of this site on the TMA) -- see OrganTumour_", s, suffix)
        next
      }
      if (length(cells) < min_cells) {
        message(g, ": only ", length(cells), " ", what, " cells (< ", min_cells, ") -- skipped")
        next
      }
      sub <- subset(xen, cells = cells)
      ann <- tryCatch(annotate_group(sub, g, sig_fun(cells), res, comp_by = "organ", prefix = prefix, harmony = FALSE),
                      error = function(e) {
                        out$failed[[g]] <<- conditionMessage(e)
                        message("!! ", g, " FAILED: ", conditionMessage(e), " -- the other sites continue.")
                        NULL
                      })
      if (!is.null(ann)) {
        plot_umap_by_core(ann$umap, sub[[]], ann$labels, g, prefix, pal_lineage)
        if (mode != "normal") plot_site_composition(sub[[]], ann$labels, s, g, prefix, what, pal_lineage)
        out[[mode]] <- c(out[[mode]], ann$labels)
      }
      rm(sub); invisible(gc())
    }
  }
  out
}

## put named labels into a metadata column (NA elsewhere)
set_labels <- function(xen, col, lab) {
  v <- if (col %in% colnames(xen[[]])) xen[[col, drop = TRUE]] else rep(NA_character_, ncol(xen))
  if (length(lab)) v[match(names(lab), colnames(xen))] <- unname(lab)
  xen[[col]] <- v
  xen
}
