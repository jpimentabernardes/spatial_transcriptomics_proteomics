## =============================================================================
## 03e_segger_prior_comparison.R
##
## segger started from the 10x Xenium cells vs. from the cellpose cells
## (03d_segger_priors.py): does segger re-segment them, which starting
## segmentation is worse, and where did the transcripts go?
##
## segger keeps the cells of the segmentation it starts from and re-assigns the
## transcripts among them. For every transcript 03d knows its source cell in the
## prior and its segger cell, so each prior is described by the FATE of its
## transcripts:
##   same cell              segger agrees with the prior
##   moved to another cell  the prior put it in the wrong cell (spill-over / merged cells)
##   removed from cells     the prior put it in a cell, segger finds it fits no cell
##   added to a cell        the prior missed it (cell too small), segger adds it
##
## How to read it -- which prior is worse?
##   * the prior segger has to CHANGE more agreed less with the transcripts;
##     "moved" is the most telling part (wrong cell, not just a wrong border);
##   * MECR (marker pairs of different lineages in one cell) of the prior vs.
##     its segger result: a bad prior has high MECR that segger lowers;
##   * if both segger runs agree with each other more than the two priors do,
##     segger converges and the choice of prior matters little; if they stay
##     apart, segger inherits the prior's cells (it learns from its nuclei).
##
## INPUT : tables/03d_{transcript_fate,gene_fate,partition_agreement}_<slide>.csv,
##         data/segmentation/segger_priors/<slide>/{cores_used.csv,cell_fate_<prior>.csv.gz},
##         data/segmentation/{cellpose,segger_xenium,segger_cellpose}/<slide>/ + Xenium outs
## OUTPUT: figures/03e_transcript_fate.png, 03e_cell_fate.png, 03e_genes_changed.png,
##         03e_agreement.png, 03e_before_after.png
##         tables/03e_verdict.csv, 03e_metrics_by_core.csv, 03e_gene_changes.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

PRIORS <- c("xenium", "cellpose")
SEGS <- c("xenium", "segger_xenium", "cellpose", "segger_cellpose")       # each prior next to its segger result
SEG_COLOURS <- c(xenium = "#00a6c7", segger_xenium = "#00596b", cellpose = "#e0a800", segger_cellpose = "#8f5f00")
FATES <- c("same cell", "moved to another cell", "removed from cells", "added to a cell", "unassigned in both")
FATE_COLOURS <- c("same cell" = "#c4c4c4", "moved to another cell" = "#d1495b", "removed from cells" = "#edae49",
                  "added to a cell" = "#2a78d6", "unassigned in both" = "#ececec")
MIN_GENE_TX <- 100            # genes with fewer transcripts in cells are not ranked in 03e_genes_changed

slides <- read_slides()
markers <- read_marker_panel(MARKERS_RNA_CSV)

read_03d <- function(stem) {
  fs <- file.path(TABLES_DIR, paste0("03d_", stem, "_", slides$slide_id, ".csv"))
  fs <- fs[file.exists(fs)]
  if (length(fs) == 0) stop("No tables/03d_", stem, "_<slide>.csv -- run 03d_segger_priors.py first.")
  dplyr::bind_rows(lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
}
fate <- read_03d("transcript_fate")
gene <- read_03d("gene_fate")
agr  <- read_03d("partition_agreement")
priors <- intersect(PRIORS, unique(fate$prior))
done_slides <- unique(fate$slide_id)
message("03d results for ", paste(done_slides, collapse = ", "), "; priors: ", paste(priors, collapse = ", "))

## -----------------------------------------------------------------------
## 1. Transcript fate: how much does segger change each prior?
## -----------------------------------------------------------------------
fate$fate <- factor(fate$fate, levels = FATES)
fate$prior <- factor(fate$prior, levels = priors)
tot <- fate |> group_by(slide_id, prior, fate) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(slide_id, prior) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
print(as.data.frame(tidyr::pivot_wider(tot[, c("slide_id", "prior", "fate", "pct")], names_from = fate,
                                       values_from = pct, values_fill = 0)), digits = 3)

p_tot <- ggplot(tot, aes(pct, prior, fill = fate)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(data = tot[tot$pct >= 3, ], aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3.6, colour = "#0b0b0b") +
  facet_wrap(~ slide_id, ncol = 1) +
  scale_fill_manual(values = FATE_COLOURS, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", panel.grid = element_blank()) +
  labs(x = "% of transcripts (core boxes)", y = "segger started from",
       title = "Where segger put each prior's transcripts")

changed <- fate |>
  group_by(slide_id, core_id, prior) |>
  summarise(pct_changed = 100 * sum(n[fate %in% FATES[2:4]]) / sum(n),
            pct_moved = 100 * sum(n[fate == "moved to another cell"]) / sum(n), .groups = "drop")
p_core <- ggplot(changed, aes(prior, pct_changed, group = core_id)) +
  geom_line(colour = "grey75") + geom_point(aes(colour = prior), size = 2.5) +
  scale_colour_manual(values = SEG_COLOURS, guide = "none") +
  facet_wrap(~ slide_id, nrow = 1) + theme_bw(base_size = 13) +
  labs(x = "segger started from", y = "% changed\n(moved + removed + added)",
       title = "Per core (each line = one core)")
save_fig(p_tot / p_core + patchwork::plot_layout(heights = c(length(done_slides), 1.6)),
         "03e_transcript_fate", width = 13, height = 4 + 2.5 * length(done_slides), dpi = 200)

## -----------------------------------------------------------------------
## 2. Per source cell: kept, cut down, or lost?
## -----------------------------------------------------------------------
cf_files <- unlist(lapply(done_slides, function(s)
  file.path(SEG_DIR, "segger_priors", s, paste0("cell_fate_", priors, ".csv.gz"))))
cf <- dplyr::bind_rows(lapply(cf_files[file.exists(cf_files)], utils::read.csv, stringsAsFactors = FALSE))
cf <- cf |> mutate(
  prior = factor(prior, levels = priors),
  pct_kept = n_kept / pmax(n_tx_prior, 1),
  pct_new = (n_gained_from_cells + n_gained_unassigned) / pmax(n_tx_segger, 1),
  cell_fate = factor(dplyr::case_when(
    n_tx_segger == 0 ~ "lost (no transcripts in segger)",
    pct_kept >= 0.8 ~ "kept >= 80% of its transcripts",
    pct_kept >= 0.5 ~ "kept 50-80%",
    TRUE ~ "kept < 50%"),
    levels = c("kept >= 80% of its transcripts", "kept 50-80%", "kept < 50%", "lost (no transcripts in segger)")))
lost <- cf |> group_by(slide_id, prior) |>
  summarise(cells = n(), lost = sum(n_tx_segger == 0),
            lost_without_nuclear_tx = sum(n_tx_segger == 0 & n_nuclear_prior == 0), .groups = "drop")
message("Prior cells segger dropped (segger only keeps cells with transcripts in their nucleus):")
print(as.data.frame(lost))

p_cf <- ggplot(cf, aes(prior, fill = cell_fate)) +
  geom_bar(position = "fill", width = 0.7, colour = "white", linewidth = 0.3) +
  facet_wrap(~ slide_id, nrow = 1) +
  scale_fill_manual(values = c("#7aa974", "#d8c56a", "#e08a4b", "#b23a48"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "segger started from", y = "prior cells", title = "What happened to each prior cell")
dens <- tidyr::pivot_longer(cf[cf$n_tx_segger > 0, c("slide_id", "prior", "pct_kept", "pct_new")],
                            c(pct_kept, pct_new), names_to = "what", values_to = "share")
dens$what <- c(pct_kept = "share of the prior cell's transcripts segger kept",
               pct_new = "share of the segger cell's transcripts that are new")[dens$what]
p_dens <- ggplot(dens, aes(share, colour = prior)) +
  geom_density(linewidth = 1, adjust = 1.2) +
  facet_grid(slide_id ~ what) +
  scale_colour_manual(values = SEG_COLOURS, name = "started from") +
  scale_x_continuous(labels = scales::percent) + theme_bw(base_size = 12) +
  labs(x = NULL, y = "cells (density)")
save_fig(p_cf | p_dens, "03e_cell_fate", width = 17, height = 3 + 3 * length(done_slides), dpi = 200)

## -----------------------------------------------------------------------
## 3. Which genes did segger move? (per prior, genes with >= MIN_GENE_TX in cells)
## -----------------------------------------------------------------------
mecr_genes <- unique(unlist(markers$mecr))
gw <- gene |> group_by(prior, feature_name, fate) |> summarise(n = sum(n), .groups = "drop") |>
  tidyr::pivot_wider(names_from = fate, values_from = n, values_fill = 0)
for (f in FATES) if (!f %in% colnames(gw)) gw[[f]] <- 0
gw <- gw |> mutate(
  n_in_prior_cells = `same cell` + `moved to another cell` + `removed from cells`,
  pct_moved = 100 * `moved to another cell` / pmax(n_in_prior_cells, 1),
  pct_removed = 100 * `removed from cells` / pmax(n_in_prior_cells, 1),
  pct_added = 100 * `added to a cell` / pmax(`same cell` + `moved to another cell` + `added to a cell`, 1),
  mecr_marker = toupper(feature_name) %in% toupper(mecr_genes))
save_table(gw, "03e_gene_changes")
min_tx <- MIN_GENE_TX
if (!any(gw$n_in_prior_cells >= min_tx)) {
  min_tx <- 1
  message("No gene has >= ", MIN_GENE_TX, " transcripts in prior cells (few cores?) -- ranking all genes")
}
top <- gw |> filter(n_in_prior_cells >= min_tx) |> group_by(prior) |>
  slice_max(pct_moved + pct_removed, n = 20, with_ties = FALSE) |> ungroup()
top_long <- tidyr::pivot_longer(top, c(pct_moved, pct_removed), names_to = "fate", values_to = "pct")
top_long$fate <- c(pct_moved = "moved to another cell", pct_removed = "removed from cells")[top_long$fate]
top_long$gene <- paste(top_long$feature_name, top_long$prior, sep = "___")
ord <- top |> arrange(prior, pct_moved + pct_removed) |> mutate(k = paste(feature_name, prior, sep = "___"))
top_long$gene <- factor(top_long$gene, levels = ord$k)
p_genes <- ggplot(top_long, aes(pct, gene, fill = fate)) +
  geom_col(width = 0.75) +
  facet_wrap(~ prior, scales = "free_y", labeller = labeller(prior = function(x) paste("started from", x))) +
  scale_y_discrete(labels = function(x) {
    g <- sub("___.*$", "", x)
    ifelse(toupper(g) %in% toupper(mecr_genes), paste0(g, " *"), g)
  }) +
  scale_fill_manual(values = FATE_COLOURS, name = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  labs(x = "% of the gene's transcripts in prior cells", y = NULL,
       title = "Genes segger re-assigned most (* = MECR lineage marker)")
save_fig(p_genes, "03e_genes_changed", width = 13, height = 8, dpi = 200)

## -----------------------------------------------------------------------
## 4. Do both segger runs converge?
## -----------------------------------------------------------------------
lev <- intersect(SEGS, unique(c(agr$a, agr$b)))
ag_all <- agr[agr$core_id == "all", ]
ag_sym <- dplyr::bind_rows(ag_all, dplyr::rename(ag_all, a = b, b = a),
                           data.frame(slide_id = rep(unique(ag_all$slide_id), each = length(lev)),
                                      a = lev, b = lev, agreement = 1))
ag_sym$a <- factor(ag_sym$a, levels = lev)
ag_sym$b <- factor(ag_sym$b, levels = rev(lev))
p_heat <- ggplot(ag_sym, aes(a, b, fill = agreement)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", agreement)), size = 4.5) +
  facet_wrap(~ slide_id, nrow = 1) +
  scale_fill_gradient(low = "#f3f3f3", high = "#2a78d6", limits = c(min(0.5, ag_sym$agreement, na.rm = TRUE), 1),
                      name = "agreement") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "Same cells? (share of transcripts that stay with the bulk of their cell)")
conv <- agr |> filter(core_id != "all", (a == "xenium" & b == "cellpose") | (a == "segger_xenium" & b == "segger_cellpose")) |>
  mutate(pair = ifelse(a == "xenium", "the two priors\n(10x vs cellpose)", "the two segger results"))
p_conv <- ggplot(conv, aes(pair, agreement, group = core_id)) +
  geom_line(colour = "grey75") + geom_point(size = 2.5, colour = "#2a78d6") +
  facet_wrap(~ slide_id, nrow = 1) + theme_bw(base_size = 12) +
  labs(x = NULL, y = "agreement", title = "Converge? (each line = one core)")
save_fig(p_heat | p_conv, "03e_agreement", width = 8 * length(done_slides) + 6, height = 6, dpi = 200)

## -----------------------------------------------------------------------
## 5. Before / after segger: MECR, counts, cells -- per core, same cores
## -----------------------------------------------------------------------
in_boxes <- function(x, y, cores) {
  out <- rep(NA_character_, length(x))
  for (i in seq_len(nrow(cores))) {
    m <- is.na(out) & x >= cores$xmin[i] & x < cores$xmax[i] & y >= cores$ymin[i] & y < cores$ymax[i]
    out[m] <- cores$core_id[i]
  }
  out
}
metrics <- dplyr::bind_rows(lapply(done_slides, function(sid) {
  cores <- utils::read.csv(file.path(SEG_DIR, "segger_priors", sid, "cores_used.csv"), stringsAsFactors = FALSE)
  srow <- slides[slides$slide_id == sid, ]
  dplyr::bind_rows(lapply(SEGS, function(mth) {
    d <- seg_dir_for(srow, mth)
    if (!dir.exists(d)) { message("  no ", mth, " output for ", sid); return(NULL) }
    counts <- read_feature_matrix(d)[["Gene Expression"]]
    cells <- read_cells_table(d)
    cells$core_id <- in_boxes(cells$x_centroid, cells$y_centroid, cores)
    cells <- cells[!is.na(cells$core_id) & cells$cell_id %in% colnames(counts), ]
    dplyr::bind_rows(lapply(split(cells$cell_id, cells$core_id), function(ids) {
      m <- counts[, ids, drop = FALSE]
      n <- Matrix::colSums(m)
      data.frame(n_cells = length(ids), median_counts = stats::median(n),
                 pct_cells_min_counts = 100 * mean(n >= MIN_COUNTS), mecr = compute_mecr(m, markers$mecr))
    }), .id = "core_id") |> mutate(slide_id = sid, method = mth, .before = 1)
  }))
}))
metrics$method <- factor(metrics$method, levels = SEGS)
metrics$prior <- factor(sub("^segger_", "", metrics$method), levels = priors)
save_table(metrics, "03e_metrics_by_core")

ml <- tidyr::pivot_longer(metrics, c(mecr, median_counts, n_cells, pct_cells_min_counts),
                          names_to = "measure", values_to = "value")
ml$measure <- factor(ml$measure, levels = c("mecr", "median_counts", "n_cells", "pct_cells_min_counts"),
                     labels = c("MECR (lower = cleaner)", "median transcripts per cell",
                                "cells", "% cells >= MIN_COUNTS"))
p_ba <- ggplot(ml, aes(method, value)) +
  geom_line(aes(group = interaction(core_id, prior)), colour = "grey75") +
  geom_point(aes(colour = method), size = 2.3) +
  stat_summary(fun = stats::median, geom = "crossbar", width = 0.5, linewidth = 0.4, colour = "#0b0b0b") +
  facet_wrap(~ measure, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = SEG_COLOURS, guide = "none") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Each prior and its segger result, per core (bar = median)")
save_fig(p_ba, "03e_before_after", width = 18, height = 5.5, dpi = 200)

## -----------------------------------------------------------------------
## 6. Verdict
## -----------------------------------------------------------------------
med <- function(mth, col) stats::median(metrics[[col]][metrics$method == mth], na.rm = TRUE)
verdict <- dplyr::bind_rows(lapply(priors, function(p) {
  f <- tot |> filter(prior == p) |> group_by(fate) |> summarise(n = sum(n), .groups = "drop")
  pct <- function(k) 100 * sum(f$n[f$fate %in% k]) / sum(f$n)
  c_p <- cf[cf$prior == p, ]
  data.frame(
    prior = p,
    pct_tx_changed = pct(FATES[2:4]), pct_tx_moved = pct(FATES[2]), pct_tx_removed = pct(FATES[3]),
    pct_tx_added = pct(FATES[4]),
    pct_cells_kept80 = 100 * mean(c_p$cell_fate == levels(c_p$cell_fate)[1]),
    pct_cells_lost = 100 * mean(c_p$n_tx_segger == 0),
    mecr_prior = med(p, "mecr"), mecr_segger = med(paste0("segger_", p), "mecr"),
    median_counts_prior = med(p, "median_counts"), median_counts_segger = med(paste0("segger_", p), "median_counts"))
}))
verdict$mecr_change_pct <- 100 * (verdict$mecr_segger - verdict$mecr_prior) / verdict$mecr_prior
print(verdict, digits = 3)
save_table(verdict, "03e_verdict")

if (length(priors) == 2) {
  v <- verdict
  more_changed <- v$prior[which.max(v$pct_tx_changed)]
  more_moved <- v$prior[which.max(v$pct_tx_moved)]
  dirtier <- v$prior[which.max(v$mecr_prior)]
  a_pri <- ag_all$agreement[ag_all$a == "xenium" & ag_all$b == "cellpose"]
  a_seg <- ag_all$agreement[ag_all$a == "segger_xenium" & ag_all$b == "segger_cellpose"]
  message("\n==== Verdict ====")
  message(sprintf("segger changed %.1f%% of the transcripts of the 10x cells and %.1f%% of the cellpose cells;",
                  v$pct_tx_changed[v$prior == "xenium"], v$pct_tx_changed[v$prior == "cellpose"]))
  message(sprintf("  moved to ANOTHER cell: %.1f%% (10x) vs %.1f%% (cellpose) -> more wrong-cell transcripts in: %s",
                  v$pct_tx_moved[v$prior == "xenium"], v$pct_tx_moved[v$prior == "cellpose"], more_moved))
  message(sprintf("MECR before segger: %.3f (10x) vs %.3f (cellpose) -> more contaminated: %s",
                  v$mecr_prior[v$prior == "xenium"], v$mecr_prior[v$prior == "cellpose"], dirtier))
  message(sprintf("MECR after segger:  %.3f (from 10x) vs %.3f (from cellpose)",
                  v$mecr_segger[v$prior == "xenium"], v$mecr_segger[v$prior == "cellpose"]))
  if (length(a_pri) && length(a_seg)) {
    message(sprintf("Agreement of the cells: priors %s, segger results %s -> %s",
                    paste(round(a_pri, 3), collapse = "/"), paste(round(a_seg, 3), collapse = "/"),
                    ifelse(mean(a_seg) > mean(a_pri), "segger brings them closer (converges)",
                           "segger keeps the priors' differences (the prior matters)")))
  }
  if (more_changed == dirtier) {
    message("=> The worse starting segmentation is: ", toupper(dirtier),
            " (segger corrects it more AND it is more contaminated before segger).")
  } else {
    message("=> Mixed: segger changes ", more_changed, " more, but ", dirtier, " is more contaminated. ",
            "Weigh 'moved' and MECR over 'removed/added' (those are mostly border / threshold effects).")
  }
}
