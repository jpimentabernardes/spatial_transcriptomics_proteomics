## =============================================================================
## 05_uc_single_gene.R
##
## UC biopsies, after 04: ONE gene at a time -- which cell types express it, and
## how its expression in the cores changes with the Nancy score.
##
## Run:   Rscript 05_uc_single_gene.R IFNG
##        Rscript 05_uc_single_gene.R IFNG IL13 PDCD1      (each gene separately, one after the other)
##        Rscript 05_uc_single_gene.R                      (all genes in UC_SINGLE_GENES)
##
## Per gene <G>, figures/05uc_gene_<G>_*.png and tables/05uc_gene_<G>_*.csv:
##   1_cell_types    which cell types carry the gene's transcripts (share of transcripts,
##                   vs that cell type's share of all transcripts = enrichment), % of each
##                   cell type's cells that are positive, transcripts per positive cell
##   2_nancy         per core, by Nancy grade (one dot = one core): transcripts per 1,000
##                   transcripts of the core, % positive cells, % positive in the main
##                   expressing cell type, transcripts per positive cell
##                   (Spearman over cores; mixed model with patient when patients have
##                   several cores)
##   3_cell_types_by_nancy   per cell type and Nancy grade: % positive cells -- which cell
##                   type's expression goes up or down with Nancy
##   4_example_cores one core per Nancy grade on black (round, with epithelium and immune
##                   cells): cell types + the gene
##
## Low-abundance genes (cytokines such as IL4, IL13, IFNG) may have only a few transcripts
## per core; a cell type then counts as "expressing" only if its share of the gene's
## transcripts is clearly above its share of all transcripts (enrichment > 1.5) and it
## has >= 1% positive cells -- otherwise the signal is close to background.
##
## INPUT : xenium_<seg>_annotated.rds (04 / 04d), tma_map Nancy / Nancy_TMA
## OUTPUT: tables/05uc_gene_<G>_*.csv, figures/05uc_gene_<G>_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_SINGLE_GENES <- get0("UC_SINGLE_GENES", ifnotfound = c("FCRL1", "FCRLA", "CD274", "CX3CR1", "IFNG", "IFNGR1",
                                                         "IFNGR2", "IL13", "IL4", "IL6", "NCAM1", "PDCD1", "SLAMF6"))
args <- commandArgs(trailingOnly = TRUE)
genes_req <- if (length(args)) toupper(trimws(unlist(strsplit(args, ",")))) else UC_SINGLE_GENES

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
rm(xen); invisible(gc())
g_obj <- match_case(genes_req, rownames(counts))
miss <- genes_req[!g_obj %in% rownames(counts)]
if (length(miss)) message("Not on the panel (skipped): ", paste(miss, collapse = ", "))
genes <- stats::setNames(g_obj[g_obj %in% rownames(counts)], genes_req[g_obj %in% rownames(counts)])
if (!length(genes)) stop("None of ", paste(genes_req, collapse = ", "), " is on the panel.")

core_tab <- uc_core_table(m)
core_tab$grade <- factor(ifelse(is.na(core_tab$nancy), "NA", paste("Nancy", core_tab$nancy)),
                         levels = c(paste("Nancy", 0:4), "NA"))
m$grade <- core_tab[m$core_id, "grade"]
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
has_nancy <- any(!is.na(core_tab$nancy))
lab_pal <- celltype_palette(m$label, m$lineage)
lab_pal_dark <- lighten(lab_pal, 0.25)
shape <- core_shape_scores(m)
shape$grade <- core_tab[shape$core_id, "grade"]
tot_core <- tapply(m$nCount_Xenium, m$core_id, sum)
area_core <- tapply(seq_len(nrow(m)), m$core_id, function(i) core_area_mm2(m$x_um[i], m$y_um[i]))
message("Cores: ", nrow(core_tab), ", patients: ", length(unique(core_tab$patient_id)), ", cells: ",
        format(nrow(m), big.mark = ","), " | genes: ", paste(names(genes), collapse = ", "))

stat_txt <- function(st, f) {
  r <- st[st$feature == f, ]
  if (!nrow(r)) return("")
  sprintf("rho %.2f, p %.2g%s", r$spearman_rho, r$p_spearman,
          ifelse(is.na(r$p_mixed_model), "", sprintf(" | with patient p %.2g", r$p_mixed_model)))
}

for (G in names(genes)) {
  message("\n==================== ", G, " ====================")
  x <- counts[genes[[G]], ]
  pos <- x > 0
  pre <- paste0("05uc_gene_", make.names(G))

  ## ---------------------------------------------------------------------
  ## 1. which cell types
  ## ---------------------------------------------------------------------
  ct <- data.frame(label = m$label, lineage = m$lineage, compartment = m$compartment, x = x, tot = m$nCount_Xenium) |>
    group_by(label, lineage, compartment) |>
    summarise(n_cells = n(), n_positive = sum(x > 0), pct_positive = 100 * mean(x > 0),
              transcripts = sum(x), transcripts_per_positive_cell = ifelse(sum(x > 0) > 0, sum(x) / sum(x > 0), NA),
              all_transcripts = sum(tot), .groups = "drop") |>
    mutate(pct_of_gene_transcripts = 100 * transcripts / max(sum(transcripts), 1),
           pct_of_all_transcripts = 100 * all_transcripts / sum(all_transcripts),
           enrichment = pct_of_gene_transcripts / pct_of_all_transcripts,
           expressing = pct_of_gene_transcripts >= 5 & enrichment >= 1.5 & pct_positive >= 1 & n_cells >= 50) |>
    arrange(desc(pct_of_gene_transcripts)) |> as.data.frame()
  save_table(ct, paste0(pre, "_by_cell_type"))
  main <- ct$label[ct$expressing]
  if (!length(main)) main <- ct$label[ct$n_cells >= 50][1]
  n_tx <- sum(x)
  message(sprintf("%s: %s transcripts in %s cells (%.2f%% of cells). Expressing cell types: %s", G,
                  format(n_tx, big.mark = ","), format(sum(pos), big.mark = ","), 100 * mean(pos),
                  if (any(ct$expressing)) paste(sprintf("%s (%.0f%% of transcripts, x%.1f)", ct$label[ct$expressing],
                                                        ct$pct_of_gene_transcripts[ct$expressing], ct$enrichment[ct$expressing]),
                                                collapse = "; ") else "none clearly above the others"))
  if (n_tx < 200) message("!! Only ", n_tx, " transcripts of ", G, " in all cores -- treat the results as indicative.")

  top <- utils::head(ct[ct$n_cells >= 20, ], 15)
  top$label <- factor(top$label, levels = rev(top$label))
  p1 <- ggplot(top, aes(pct_of_gene_transcripts, label, fill = as.character(label))) + geom_col(width = 0.75) +
    geom_point(aes(x = pct_of_all_transcripts), shape = 124, size = 5, colour = "grey15") +
    geom_text(aes(x = pmax(pct_of_gene_transcripts, pct_of_all_transcripts), label = sprintf("x%.1f", enrichment)),
              hjust = -0.3, size = 3.2) +
    scale_fill_manual(values = lab_pal, guide = "none") + scale_x_continuous(expand = expansion(c(0, 0.15))) +
    theme_bw(base_size = 12) +
    labs(x = sprintf("%% of all %s transcripts", G), y = NULL, title = "A. Where are the transcripts?",
         subtitle = "black tick = the cell type's share of ALL transcripts;\nx = enrichment (bar / tick): > 1.5 = preferentially expressed")
  top2 <- top
  p2 <- ggplot(top2, aes(pct_positive, label, fill = as.character(label))) + geom_col(width = 0.75) +
    geom_text(aes(label = ifelse(is.na(transcripts_per_positive_cell), "",
                                 sprintf("%.0f%%  (%.1f per +cell)", pct_positive, transcripts_per_positive_cell))),
              hjust = -0.08, size = 3.2) +
    scale_fill_manual(values = lab_pal, guide = "none") + scale_x_continuous(expand = expansion(c(0, 0.45))) +
    theme_bw(base_size = 12) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "% of the cell type's cells positive", y = NULL, title = "B. How many of its cells express it?",
         subtitle = "(transcripts per positive cell)")
  save_fig((p1 | p2) + plot_annotation(title = sprintf("%s: which cell types express it?", G),
                                       subtitle = sprintf("%s transcripts in %s positive cells (%.2f%% of all cells)",
                                                          format(n_tx, big.mark = ","), format(sum(pos), big.mark = ","),
                                                          100 * mean(pos))),
           paste0(pre, "_1_cell_types"), width = 15, height = max(5.5, 0.38 * nrow(top) + 2.5))

  ## ---------------------------------------------------------------------
  ## 2. per core vs Nancy
  ## ---------------------------------------------------------------------
  main1 <- main[1]
  in_main <- m$label %in% main1
  pc <- data.frame(core_id = core_tab$core_id, row.names = core_tab$core_id)
  pc$transcripts_per_1000 <- 1000 * tapply(x, m$core_id, sum)[pc$core_id] / tot_core[pc$core_id]
  pc$transcripts_per_mm2 <- tapply(x, m$core_id, sum)[pc$core_id] / area_core[pc$core_id]
  pc$pct_positive_all_cells <- 100 * tapply(pos, m$core_id, mean)[pc$core_id]
  n_main <- tapply(in_main, m$core_id, sum)[pc$core_id]
  pc$pct_positive_main <- ifelse(n_main >= UC_MIN_DENOM, 100 * tapply(pos & in_main, m$core_id, sum)[pc$core_id] / n_main, NA)
  npos <- tapply(pos, m$core_id, sum)[pc$core_id]
  pc$transcripts_per_positive_cell <- ifelse(npos >= 3, tapply(x, m$core_id, sum)[pc$core_id] / npos, NA)
  pc <- cbind(core_tab[pc$core_id, c("slide_id", "patient_id", "nancy", "grade", "n_cells")], pc)
  save_table(pc, paste0(pre, "_per_core"))
  metrics <- c(transcripts_per_1000 = sprintf("%s transcripts per 1,000 transcripts", G),
               pct_positive_all_cells = sprintf("%% of all cells positive for %s", G),
               pct_positive_main = sprintf("%% of %s cells positive", main1),
               transcripts_per_positive_cell = "transcripts per positive cell")
  if (has_nancy) {
    st <- assoc_nancy(pc, names(metrics))
    save_table(st, paste0(pre, "_vs_nancy"))
    long <- tidyr::pivot_longer(pc[, c("core_id", "grade", names(metrics))], dplyr::all_of(names(metrics)),
                                names_to = "metric", values_to = "value")
    long <- long[!is.na(long$value), ]
    long$metric <- factor(metrics[long$metric], levels = metrics)
    lab <- data.frame(metric = factor(metrics, levels = metrics), txt = vapply(names(metrics), function(f) stat_txt(st, f), ""))
    p <- ggplot(long, aes(grade, value, colour = grade)) +
      geom_boxplot(outlier.shape = NA, colour = "grey50") +
      geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1.6) +
      geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.03, vjust = 1.3, size = 3.2,
                colour = "grey20", inherit.aes = FALSE) +
      scale_colour_manual(values = grade_cols, guide = "none") + scale_y_continuous(expand = expansion(c(0.05, 0.15))) +
      facet_wrap(~ metric, scales = "free_y", ncol = 2) + theme_bw(base_size = 12) +
      labs(x = NULL, y = NULL, title = sprintf("%s per core, by Nancy grade (one dot = one core)", G),
           subtitle = paste0("main expressing cell type: ", main1,
                             "\n'% positive' up = expressed in more cells; 'per positive cell' up = more per cell"))
    save_fig(p, paste0(pre, "_2_nancy"), width = 12, height = 9)

    ## -------------------------------------------------------------------
    ## 3. per cell type and Nancy grade
    ## -------------------------------------------------------------------
    show_ct <- unique(c(main, utils::head(ct$label[ct$n_cells >= 50], 8)))
    cg <- data.frame(label = m$label, grade = m$grade, core_id = m$core_id, x = x) |> filter(label %in% show_ct) |>
      group_by(label, grade) |>
      summarise(n_cells = n(), pct_positive = 100 * mean(x > 0), mean_counts = mean(x), .groups = "drop") |>
      filter(n_cells >= 20)
    ## per cell type: Spearman over cores of % positive vs Nancy
    rho_ct <- vapply(show_ct, function(l) {
      d <- data.frame(core_id = m$core_id, l = m$label == l, p = pos) |> group_by(core_id) |>
        summarise(n = sum(l), v = ifelse(sum(l) >= UC_MIN_DENOM, 100 * sum(p & l) / sum(l), NA), .groups = "drop")
      d$nancy <- core_tab[d$core_id, "nancy"]
      d <- d[!is.na(d$v) & !is.na(d$nancy), ]
      if (nrow(d) < 6 || stats::sd(d$v) == 0) return(NA_character_)
      t <- suppressWarnings(stats::cor.test(d$v, d$nancy, method = "spearman", exact = FALSE))
      sprintf("rho %+.2f, p %.2g", t$estimate, t$p.value)
    }, character(1))
    save_table(cbind(cg, rho = rho_ct[cg$label]), paste0(pre, "_cell_types_by_nancy"))
    cg$label_txt <- factor(sprintf("%s\n(%s)", cg$label, ifelse(is.na(rho_ct[cg$label]), "too few cores", rho_ct[cg$label])),
                           levels = rev(sprintf("%s\n(%s)", show_ct, ifelse(is.na(rho_ct[show_ct]), "too few cores", rho_ct[show_ct]))))
    p <- ggplot(cg, aes(grade, label_txt, size = pct_positive, colour = mean_counts)) + geom_point() +
      scale_size_area(max_size = 12, name = "% cells\npositive") +
      scale_colour_gradient(low = "#C9C9C9", high = "#B2182B", name = "mean\ncounts / cell") +
      theme_bw(base_size = 11) +
      labs(x = NULL, y = NULL, title = sprintf("%s per cell type and Nancy grade", G),
           subtitle = "rho = Spearman over cores of the cell type's % positive vs Nancy")
    save_fig(p, paste0(pre, "_3_cell_types_by_nancy"), width = 10, height = max(5, 0.62 * length(show_ct) + 2))
  }

  ## ---------------------------------------------------------------------
  ## 4. one example core per Nancy grade, on black
  ## ---------------------------------------------------------------------
  sh <- shape[shape$n_cells >= 300, ]
  sh$pct_pos <- pc[sh$core_id, "pct_positive_all_cells"]
  sh$pick <- sh$score * (0.5 + pmin(sh$pct_pos / max(stats::median(sh$pct_pos, na.rm = TRUE), 0.1), 2))
  grades <- if (has_nancy) intersect(levels(core_tab$grade), as.character(sh$grade)) else "NA"
  grades <- setdiff(grades, "NA")
  pick <- unlist(lapply(grades, function(g) { s <- sh[sh$grade == g, ]; s$core_id[which.max(s$pick)] }))
  if (!length(pick)) pick <- utils::head(sh$core_id[order(-sh$pick)], 4)
  if (length(pick)) {
    cols <- lapply(pick, function(cid) {
      i <- which(m$core_id == cid); d <- m[i, ]
      list(plot_core_black(d, "label", lab_pal_dark, title = sprintf("%s  -  Nancy %s", cid, core_tab[cid, "nancy"]),
                           legend = FALSE, panel_in = 4.2),
           plot_core_gene_black(d, x[i], G, title = G, panel_in = 4.2) +
             (if (cid != utils::tail(pick, 1)) theme(legend.position = "none") else theme()))
    })
    plots <- c(lapply(cols, `[[`, 1), lapply(cols, `[[`, 2))
    fig <- patchwork::wrap_plots(plots, ncol = length(pick), byrow = TRUE)
    leg <- legend_black(plot_core_black(m[m$core_id %in% pick, ], "label", lab_pal_dark))
    if (!is.null(leg)) fig <- patchwork::wrap_plots(fig, leg, widths = c(4.4 * length(pick), 2.6))
    fig <- fig + patchwork::plot_annotation(
      title = sprintf("%s: one core per Nancy grade", G),
      subtitle = "top: cell types; bottom: the gene (dark grey = 0 transcripts); bar = 100 um",
      theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                    plot.title = element_text(colour = "white", face = "bold", size = 17),
                    plot.subtitle = element_text(colour = "grey80", size = 11)))
    save_fig(fig, paste0(pre, "_4_example_cores"), width = 4.4 * length(pick) + 3, height = 10, dpi = 200)
  }
}
message("\nDone.")
