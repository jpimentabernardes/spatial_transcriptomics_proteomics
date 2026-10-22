## =============================================================================
## 05_uc_genes_eosinophils.R
##
## UC biopsy TMAs, after 04 (UC marker panel). Two focused questions:
##
## A. Genes of interest (default KRT7, CXCL1, SLC6A14; EXP1_UC_GENES=KRT7,CXCL1 to change)
##   A1 which cell types express them: % positive cells, mean counts, and where the
##      gene's transcripts go (share per cell type vs that cell type's share of all
##      transcripts) -- overall and in low vs high Nancy cores
##   A2 per core: transcripts per 1,000 transcripts of the core, and % positive cells
##      in the compartment that carries most of the gene -- vs Nancy grade, high vs low
##      (Spearman over cores, Wilcoxon, mixed model with patient)
##   A3 per patient (mean over cores) vs Disease.Control (and the remission outcomes)
##   A4 example cores on a black background: cell types + each gene, low vs high Nancy
##
## B. Eosinophils
##   B1 eosinophils per core (% of cells, % of immune cells, per mm^2), by Nancy and vs
##      the clinical eosinophil columns; cores with a high prevalence (above Q3 + 1.5 IQR
##      of eosinophils per mm^2) are flagged
##   B2 is the eosinophil label right? SIGLEC8 / RNASE3 / ALOX15 (+ other eosinophil genes on
##      the panel) per cell type, and a marker-based eosinophil count per core for comparison
##   B3 close-up of the eosinophil-rich cores (black background): all cell types with the
##      eosinophils ringed, the eosinophils' neighbours, the core's composition vs the cohort,
##      and which cell types sit next to eosinophils (vs the rest of the core)
##
## Eosinophils = cells whose label matches "eosin" (UC_CLASS_PATTERNS in R/uc_helpers.R).
##
## INPUT : xenium_<seg>_annotated.rds (04, UC slides), tma_map Nancy / Nancy_TMA,
##         de-identified clinical table (UC_CLINICAL_CSV, optional)
## OUTPUT: tables/05uc_genes_*.csv, tables/05uc_eos_*.csv,
##         figures/05uc_genes_*.png, figures/05uc_eos_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_GENES <- trimws(strsplit(Sys.getenv("EXP1_UC_GENES", paste(get0("UC_GENES", ifnotfound = c("KRT7", "CXCL1", "SLC6A14")),
                                                             collapse = ",")), ",")[[1]])
UC_EOS_GENES <- get0("UC_EOS_GENES", ifnotfound = c("SIGLEC8", "RNASE3", "ALOX15", "CLC", "CCR3", "IL5RA", "PRG2", "EPX", "RNASE2"))
UC_EOS_MIN_COUNTS <- get0("UC_EOS_MIN_COUNTS", ifnotfound = 2)   # marker-based call: >= 2 eosinophil-gene counts
UC_EOS_TOP_N <- get0("UC_EOS_TOP_N", ifnotfound = 6)             # cores shown close up
UC_EOS_NEIGHBOUR_UM <- get0("UC_EOS_NEIGHBOUR_UM", ifnotfound = 30)
PARTS <- trimws(strsplit(Sys.getenv("EXP1_UC_PARTS", "genes,eos"), ",")[[1]])

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
core_tab <- uc_core_table(m)
clin <- join_uc_clinical(core_tab, read_uc_clinical())
m$nancy <- core_tab[m$core_id, "nancy"]
m$nancy_group <- core_tab[m$core_id, "nancy_group"]
grade_lab <- function(n) factor(ifelse(is.na(n), "NA", paste("Nancy", n)), levels = c(paste("Nancy", 0:4), "NA"))
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
core_tab$grade <- grade_lab(core_tab$nancy)
core_tab$area_mm2 <- vapply(core_tab$core_id, function(cid) { i <- m$core_id == cid; core_area_mm2(m$x_um[i], m$y_um[i]) }, numeric(1))
lab_pal <- celltype_palette(m$label, m$lineage)
lab_pal_dark <- lighten(lab_pal, 0.25)
oc <- if (is.null(clin)) character(0) else intersect(UC_OUTCOME_COLS, colnames(clin))
message("Nancy from ", attr(core_tab, "nancy_col"), "; ", nrow(core_tab), " cores, ",
        length(unique(core_tab$patient_id)), " patients")

## patient table: mean of a per-core value + outcomes + highest Nancy
patient_values <- function(core_vals, cols) {
  d <- cbind(core_tab[core_vals$core_id, c("core_id", "patient_id", "nancy")], core_vals[, cols, drop = FALSE])
  if (length(oc)) d <- cbind(d, clin[d$core_id, oc, drop = FALSE])
  pat <- d |> group_by(patient_id) |>
    summarise(n_cores = n(), nancy_max = suppressWarnings(max(nancy, na.rm = TRUE)),
              dplyr::across(dplyr::all_of(cols), ~ mean(.x, na.rm = TRUE)),
              dplyr::across(dplyr::all_of(oc), ~ dplyr::first(stats::na.omit(.x))), .groups = "drop") |> as.data.frame()
  pat$nancy_max[!is.finite(pat$nancy_max)] <- NA
  for (o in oc) pat[[o]] <- as_outcome(pat[[o]])
  pat
}
## boxplots of per-core values by Nancy grade, with Spearman / high-vs-low p per panel
plot_by_grade <- function(d, cols, title, ncol = 3) {
  d$grade <- core_tab[d$core_id, "grade"]; d$nancy <- core_tab[d$core_id, "nancy"]
  d$patient_id <- core_tab[d$core_id, "patient_id"]
  st <- assoc_nancy(d, cols)
  long <- tidyr::pivot_longer(d[, c("core_id", "grade", cols)], dplyr::all_of(cols), names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric, levels = cols)
  lab <- if (nrow(st)) data.frame(metric = factor(st$feature, levels = cols),
                                  txt = sprintf("rho %.2f (p %.2g)\nhigh vs low p %.2g", st$spearman_rho, st$p_spearman,
                                                st$p_wilcox_high_vs_low)) else NULL
  p <- ggplot(long[!is.na(long$value), ], aes(grade, value, colour = grade)) +
    geom_boxplot(outlier.shape = NA, colour = "grey50") +
    geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1.4) +
    scale_colour_manual(values = grade_cols, guide = "none") +
    facet_wrap(~ metric, scales = "free_y", ncol = ncol) + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = NULL, title = title)
  if (!is.null(lab)) p <- p + geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.2,
                                         size = 2.9, colour = "grey20", inherit.aes = FALSE)
  list(plot = p, stats = st)
}
## patient values vs one outcome, coloured by the patient's highest Nancy
plot_by_outcome <- function(pat, cols, o, title) {
  st <- assoc_outcome(pat, cols, o)
  long <- tidyr::pivot_longer(pat[!is.na(pat[[o]]), c("patient_id", "nancy_max", o, cols)], dplyr::all_of(cols),
                              names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric, levels = cols)
  long$nancy_group <- nancy_group(long$nancy_max)
  p <- ggplot(long[!is.na(long$value), ], aes(.data[[o]], value)) +
    geom_boxplot(outlier.shape = NA, colour = "grey45") +
    geom_point(aes(colour = nancy_group), position = position_jitter(width = 0.15, seed = 1), size = 1.8) +
    scale_colour_manual(values = GROUP_COLOURS(), name = "patient's highest\nNancy", na.value = "grey60") +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) + theme_bw(base_size = 11) +
    labs(x = o, y = "patient mean over cores", title = title)
  if (nrow(st)) p <- p + geom_text(data = data.frame(metric = factor(st$feature, levels = cols),
                                                      txt = sprintf("AUC %.2f, p %.2g", st$auc, st$p)),
                                   aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.3, size = 3, inherit.aes = FALSE)
  list(plot = p, stats = st)
}

## =============================================================================
## A. Genes of interest
## =============================================================================
if ("genes" %in% PARTS) {
g_obj <- match_case(UC_GENES, rownames(counts))
miss <- UC_GENES[!g_obj %in% rownames(counts)]
if (length(miss)) message("Not on the panel (left out): ", paste(miss, collapse = ", "))
G <- stats::setNames(g_obj[g_obj %in% rownames(counts)], UC_GENES[g_obj %in% rownames(counts)])
if (!length(G)) stop("None of ", paste(UC_GENES, collapse = ", "), " is on the panel.")
gc_ <- as.matrix(Matrix::t(counts[G, , drop = FALSE])); colnames(gc_) <- names(G)
gn_ <- t(uc_lognorm(counts[G, , drop = FALSE], m$nCount_Xenium)); colnames(gn_) <- names(G)

## A1 which cell types
by_label <- function(sel, extra = NULL) {
  dplyr::bind_rows(lapply(names(G), function(g) {
    d <- data.frame(label = m$label[sel], lineage = m$lineage[sel], cnt = gc_[sel, g], nrm = gn_[sel, g],
                    tot = m$nCount_Xenium[sel])
    if (!is.null(extra)) d$grp <- extra[sel]
    by <- c("label", "lineage", if (!is.null(extra)) "grp")
    x <- d |> group_by(dplyr::across(dplyr::all_of(by))) |>
      summarise(n_cells = n(), pct_positive = 100 * mean(cnt > 0), mean_counts = mean(cnt), mean_norm = mean(nrm),
                transcripts = sum(cnt), all_transcripts = sum(tot), .groups = "drop")
    if (!is.null(extra)) x <- group_by(x, grp)                # shares within each Nancy group
    x |> mutate(gene = g, pct_of_gene_transcripts = 100 * transcripts / max(sum(transcripts), 1),
                pct_of_all_transcripts = 100 * all_transcripts / sum(all_transcripts),
                enrichment = pct_of_gene_transcripts / pct_of_all_transcripts) |> ungroup()
  }))
}
a1 <- by_label(rep(TRUE, nrow(m)))
save_table(a1 |> arrange(gene, desc(pct_of_gene_transcripts)), "05uc_genes_by_cell_type")
for (g in names(G)) {
  t3 <- a1 |> filter(gene == g, n_cells >= 50) |> arrange(desc(pct_of_gene_transcripts)) |> head(4)
  message(g, ": ", paste(sprintf("%s %.0f%% of transcripts (%.0f%% cells +)", t3$label, t3$pct_of_gene_transcripts,
                                 t3$pct_positive), collapse = "; "))
}
a1$label <- factor(a1$label, levels = rev(names(lab_pal)))
a1$gene <- factor(a1$gene, levels = names(G))
p_dot <- ggplot(a1[a1$n_cells >= 20, ], aes(gene, label, size = pct_positive, colour = mean_counts)) + geom_point() +
  scale_size_area(max_size = 9, name = "% cells\npositive") +
  scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\ncounts / cell") +
  theme_bw(base_size = 11) + labs(x = NULL, y = NULL, title = "Which cell types express them?")
sh <- a1 |> group_by(gene) |> slice_max(pct_of_gene_transcripts, n = 8) |> ungroup()
sh$label_g <- factor(paste(sh$label, sh$gene, sep = "___"),
                     levels = paste(sh$label, sh$gene, sep = "___")[order(sh$gene, sh$pct_of_gene_transcripts)])
p_share <- ggplot(sh, aes(pct_of_gene_transcripts, label_g, fill = as.character(label))) + geom_col(width = 0.75) +
  geom_point(aes(x = pct_of_all_transcripts), shape = 124, size = 4, colour = "grey20") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) + scale_fill_manual(values = lab_pal, guide = "none") +
  facet_wrap(~ gene, scales = "free_y", ncol = 1) + theme_bw(base_size = 11) +
  labs(x = "% of the gene's transcripts", y = NULL, title = "Where do the transcripts go? (top 8)",
       subtitle = "black tick = that cell type's share of ALL transcripts;\nbar beyond the tick = the gene is enriched there")
save_fig(p_dot + p_share + plot_layout(widths = c(1, 1.1)), "05uc_genes_by_cell_type", width = 15,
         height = max(8, 0.22 * length(unique(a1$label)) + 2))

## by Nancy group
if (any(!is.na(m$nancy_group))) {
  a1n <- by_label(!is.na(m$nancy_group), extra = as.character(m$nancy_group))
  save_table(a1n, "05uc_genes_by_cell_type_and_nancy")
  top_lab <- a1 |> filter(n_cells >= 50) |> group_by(gene) |> slice_max(pct_of_gene_transcripts, n = 6) |> pull(label) |>
    as.character() |> unique()
  dd <- a1n[a1n$label %in% top_lab & a1n$n_cells >= 20, ]
  dd$label <- factor(dd$label, levels = rev(intersect(names(lab_pal), top_lab)))
  dd$grp <- factor(dd$grp, levels = levels(nancy_group(0)))
  dd$gene <- factor(dd$gene, levels = names(G))
  p <- ggplot(dd, aes(grp, label, size = pct_positive, colour = mean_counts)) + geom_point() +
    scale_size_area(max_size = 9, name = "% cells\npositive") +
    scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\ncounts / cell") +
    facet_wrap(~ gene, nrow = 1) + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(x = NULL, y = NULL, title = "Same cell types, low vs high Nancy cores",
         subtitle = "more positive cells WITHIN a cell type = higher expression per cell, not just more of that cell type")
  save_fig(p, "05uc_genes_by_cell_type_and_nancy", width = 4 * length(G) + 3, height = max(5, 0.35 * length(top_lab) + 2.5))
}

## A2 per core
main_comp <- vapply(names(G), function(g) {
  s <- tapply(gc_[, g], m$compartment, sum); names(s)[which.max(s)] }, character(1))
message("Compartment with most transcripts: ", paste(names(main_comp), main_comp, sep = " -> ", collapse = ", "))
pc <- data.frame(core_id = core_tab$core_id, row.names = core_tab$core_id)
tot_core <- tapply(m$nCount_Xenium, m$core_id, sum)
for (g in names(G)) {
  pc[[paste0(g, "_per_1000_transcripts")]] <- 1000 * tapply(gc_[, g], m$core_id, sum)[pc$core_id] / tot_core[pc$core_id]
  pc[[paste0(g, "_pct_positive_all_cells")]] <- 100 * tapply(gc_[, g] > 0, m$core_id, mean)[pc$core_id]
  sel <- m$compartment == main_comp[[g]]
  n_c <- tapply(sel, m$core_id, sum)[pc$core_id]
  v <- 100 * tapply(gc_[sel, g] > 0, m$core_id[sel], mean)[pc$core_id]
  pc[[sprintf("%s_pct_positive_%s", g, main_comp[[g]])]] <- ifelse(n_c >= UC_MIN_DENOM, v, NA)
}
cols_core <- setdiff(colnames(pc), "core_id")
out <- cbind(core_tab[, c("slide_id", "core_id", "patient_id", "nancy", "n_cells")], pc[core_tab$core_id, cols_core])
save_table(out, "05uc_genes_per_core")
if (any(!is.na(core_tab$nancy))) {
  show <- c(paste0(names(G), "_per_1000_transcripts"), sprintf("%s_pct_positive_%s", names(G), main_comp))
  r <- plot_by_grade(pc, show, "Genes of interest per core, by Nancy grade (one dot = one core)", ncol = length(G))
  save_table(r$stats, "05uc_genes_vs_nancy")
  save_fig(r$plot + labs(subtitle = paste0("top row: transcripts per 1,000 transcripts of the core; bottom row: % positive ",
                                            "cells in the compartment with most of the gene's transcripts.\nrho = Spearman over cores; ",
                                            "high vs low = Wilcoxon, Nancy >= ", UC_NANCY_HIGH_MIN, " vs below")),
           "05uc_genes_vs_nancy", width = 5 * length(G), height = 9)
}

## A3 per patient vs Disease.Control (and the other outcomes)
if (length(oc)) {
  show <- c(paste0(names(G), "_per_1000_transcripts"), sprintf("%s_pct_positive_%s", names(G), main_comp))
  pat <- patient_values(pc, show)
  save_table(pat, "05uc_genes_per_patient")
  st <- dplyr::bind_rows(lapply(oc, function(o) assoc_outcome(pat, show, o)))
  save_table(st, "05uc_genes_vs_outcomes")
  for (o in intersect(c("Disease.Control", oc[1]), oc)[1]) {
    r <- plot_by_outcome(pat, show, o, sprintf("Genes of interest vs %s (one dot = one patient, mean of its cores)", o))
    save_fig(r$plot, paste0("05uc_genes_vs_", make.names(o)), width = 5 * length(G), height = 8.5)
  }
} else message("A3 skipped: no clinical outcomes (UC_CLINICAL_CSV).")

## A4 example cores, black background: two low-Nancy and two high-Nancy cores
##    (in each group: the core closest to the group's median expression, and the highest)
if (any(!is.na(core_tab$nancy))) {
  z <- scale(as.matrix(pc[, paste0(names(G), "_per_1000_transcripts"), drop = FALSE]))
  z[is.na(z)] <- 0
  sc <- stats::setNames(rowMeans(z), rownames(pc))
  ok <- core_tab$n_cells >= 300
  pick <- unlist(lapply(levels(core_tab$nancy_group), function(gr) {
    cc <- core_tab$core_id[core_tab$nancy_group %in% gr & ok]
    if (!length(cc)) return(NULL)
    s <- sc[cc]
    unique(c(cc[which.min(abs(s - stats::median(s)))], cc[which.max(s)]))
  }))
  if (length(pick)) {
    rows <- list()
    for (cid in pick) {
      d <- m[m$core_id == cid, ]
      ttl <- sprintf("%s  (Nancy %s)", cid, core_tab[cid, "nancy"])
      col_list <- list(plot_core_black(d, "label", lab_pal_dark, title = ttl, legend = FALSE, panel_in = 4))
      for (g in names(G)) {                                   # one colour bar per row: last core only
        col_list[[g]] <- plot_core_gene_black(d, gc_[m$core_id == cid, g], g, title = g, panel_in = 4)
        if (cid != utils::tail(pick, 1)) col_list[[g]] <- col_list[[g]] + theme(legend.position = "none")
      }
      rows[[cid]] <- col_list
    }
    ## columns = cores, rows = cell types + genes
    plots <- unlist(lapply(seq_len(length(G) + 1), function(r) lapply(rows, `[[`, r)), recursive = FALSE)
    leg <- legend_black(plot_core_black(m[m$core_id %in% pick, ], "label", lab_pal_dark))
    fig <- patchwork::wrap_plots(plots, ncol = length(pick), byrow = TRUE)
    if (!is.null(leg)) fig <- patchwork::wrap_plots(fig, leg, widths = c(4 * length(pick), 2.5))
    fig <- fig + patchwork::plot_annotation(
      title = paste("Example cores:", paste(names(G), collapse = ", ")),
      subtitle = "top row: cell types; then one row per gene (dark grey = 0 counts); left: low Nancy, right: high Nancy",
      theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                    plot.title = element_text(colour = "white", face = "bold", size = 16),
                    plot.subtitle = element_text(colour = "grey80")))
    save_fig(fig, "05uc_genes_example_cores", width = 4.3 * length(pick) + 3, height = 4.4 * (length(G) + 1) + 1, dpi = 200)
  }
}
} ## end genes

## =============================================================================
## B. Eosinophils
## =============================================================================
if ("eos" %in% PARTS) {
is_eos <- m$class %in% "Eosinophil"
if (!any(is_eos)) message("!! No cell label matches 'eosin' -- eosinophils are counted with the marker-based call only. ",
                          "Label a cluster 'Eosinophil' in 04 (or adapt UC_CLASS_PATTERNS) for the label-based count.")
eg <- match_case(UC_EOS_GENES, rownames(counts)); eg <- eg[eg %in% rownames(counts)]
message("Eosinophil genes on the panel: ", paste(eg, collapse = ", "))
eos_marker <- if (length(eg)) {
  ec <- Matrix::colSums(counts[eg, , drop = FALSE])
  ec >= UC_EOS_MIN_COUNTS & m$compartment != "Epithelium"
} else rep(FALSE, nrow(m))
if (!any(is_eos)) is_eos <- eos_marker
m$eos <- is_eos
m$eos_marker <- eos_marker

## B1 per core
eos_core <- m |> group_by(core_id) |>
  summarise(n_cells = n(), n_immune = sum(compartment == "Immune"), n_eos = sum(eos), n_eos_marker = sum(eos_marker),
            .groups = "drop") |> as.data.frame()
rownames(eos_core) <- eos_core$core_id
eos_core$area_mm2 <- core_tab[eos_core$core_id, "area_mm2"]
eos_core$eos_pct_cells <- 100 * eos_core$n_eos / eos_core$n_cells
eos_core$eos_pct_immune <- ifelse(eos_core$n_immune >= UC_MIN_DENOM, 100 * eos_core$n_eos / eos_core$n_immune, NA)
eos_core$eos_per_mm2 <- eos_core$n_eos / eos_core$area_mm2
eos_core$eos_marker_per_mm2 <- eos_core$n_eos_marker / eos_core$area_mm2
## intraepithelial: >= half of the 10 nearest neighbours epithelial
ie <- vapply(eos_core$core_id, function(cid) {
  i <- which(m$core_id == cid); e <- which(m$eos[i])
  if (length(e) < 3 || length(i) < 12) return(NA_real_)
  nn <- RANN::nn2(m[i, c("x_um", "y_um")], m[i[e], c("x_um", "y_um")], k = 11)$nn.idx[, -1, drop = FALSE]
  100 * mean(rowMeans(matrix(m$compartment[i][nn] == "Epithelium", nrow = length(e))) >= 0.5)
}, numeric(1))
eos_core$eos_intraepithelial_pct <- ie
q <- stats::quantile(eos_core$eos_per_mm2, c(0.25, 0.75), na.rm = TRUE)
thr <- q[2] + 1.5 * diff(q)
eos_core$high_eosinophils <- eos_core$eos_per_mm2 > thr & eos_core$n_eos >= 10
eos_core <- cbind(core_tab[eos_core$core_id, c("slide_id", "patient_id", "nancy", "grade")], eos_core)
if (!is.null(clin)) for (c in intersect(UC_CLIN_EOS_COLS, colnames(clin))) eos_core[[c]] <- clin[eos_core$core_id, c]
eos_core <- eos_core[order(-eos_core$eos_per_mm2), ]
save_table(eos_core, "05uc_eos_per_core")
message(sprintf("Eosinophils: %s cells (%.2f%% of all); high-eosinophil cores (> %.0f / mm^2, Q3 + 1.5 IQR): %d -- %s",
                format(sum(m$eos), big.mark = ","), 100 * mean(m$eos), thr, sum(eos_core$high_eosinophils),
                paste(eos_core$core_id[eos_core$high_eosinophils], collapse = ", ")))
top_cores <- utils::head(eos_core$core_id[eos_core$n_eos >= 3], UC_EOS_TOP_N)

d <- eos_core; d$core_id <- factor(d$core_id, levels = d$core_id)
p1 <- ggplot(d, aes(core_id, eos_per_mm2, fill = grade)) + geom_col(width = 0.85) +
  geom_hline(yintercept = thr, linetype = 2, colour = "grey30") +
  annotate("text", x = nrow(d), y = thr, label = "Q3 + 1.5 IQR", hjust = 1, vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = grade_cols, name = NULL) + scale_y_continuous(expand = expansion(c(0, 0.05))) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), legend.position = "top") +
  labs(x = "core (ranked)", y = "eosinophils per mm2", title = "Eosinophils per core, ranked (colour = Nancy grade)")
p2 <- ggplot(eos_core, aes(grade, eos_per_mm2, colour = grade)) + geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15, seed = 1)) + scale_colour_manual(values = grade_cols, guide = "none") +
  theme_bw(base_size = 11) + labs(x = NULL, y = "eosinophils per mm2", title = "by Nancy grade")
pl <- p1 / (p2 | patchwork::plot_spacer())
eos_clin <- intersect(UC_CLIN_EOS_COLS, colnames(eos_core))
eos_clin <- eos_clin[vapply(eos_core[eos_clin], function(v) sum(!is.na(v)) >= 3, logical(1))]
if (length(eos_clin)) {
  ec <- eos_clin[1]
  p3 <- ggplot(eos_core[!is.na(eos_core[[ec]]), ], aes(factor(.data[[ec]]), eos_per_mm2)) +
    geom_boxplot(outlier.shape = NA, colour = "grey50") +
    geom_point(aes(colour = grade), position = position_jitter(width = 0.15, seed = 1)) +
    scale_colour_manual(values = grade_cols, name = NULL) + theme_bw(base_size = 11) +
    labs(x = paste("clinical:", ec), y = "eosinophils per mm2 (TMA)", title = "vs the pathology report")
  pl <- p1 / (p2 | p3)
}
save_fig(pl, "05uc_eos_per_core", width = max(14, 0.1 * nrow(eos_core) + 4), height = 10)
if (any(!is.na(core_tab$nancy))) {
  st <- assoc_nancy(eos_core, c("eos_per_mm2", "eos_pct_cells", "eos_pct_immune", "eos_intraepithelial_pct"))
  save_table(st, "05uc_eos_vs_nancy")
}
if (length(oc)) {
  rownames(eos_core) <- eos_core$core_id
  ev <- c("eos_per_mm2", "eos_pct_immune", "eos_intraepithelial_pct")
  pat <- patient_values(eos_core, ev)
  save_table(dplyr::bind_rows(lapply(oc, function(o) assoc_outcome(pat, ev, o))), "05uc_eos_vs_outcomes")
  o <- if ("Disease.Control" %in% oc) "Disease.Control" else oc[1]
  r <- plot_by_outcome(pat, ev, o, sprintf("Eosinophils vs %s (one dot = one patient)", o))
  save_fig(r$plot, paste0("05uc_eos_vs_", make.names(o)), width = 13, height = 5)
}

## B2 eosinophil genes per cell type (is the label right?)
if (length(eg)) {
  dd <- dplyr::bind_rows(lapply(eg, function(g) {
    x <- counts[g, ]
    data.frame(label = m$label, x = x) |> group_by(label) |>
      summarise(n_cells = n(), pct_positive = 100 * mean(x > 0), mean_counts = mean(x), .groups = "drop") |> mutate(gene = g)
  }))
  dd <- dd[dd$n_cells >= 20, ]
  dd$label <- factor(dd$label, levels = rev(names(lab_pal))); dd$gene <- factor(dd$gene, levels = eg)
  p <- ggplot(dd, aes(gene, label, size = pct_positive, colour = mean_counts)) + geom_point() +
    scale_size_area(max_size = 9, name = "% cells\npositive") +
    scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\ncounts / cell") +
    theme_bw(base_size = 11) + labs(x = NULL, y = NULL, title = "Eosinophil genes per cell type",
                                    subtitle = "the eosinophil label should carry SIGLEC8 / RNASE3 (ALOX15 is also epithelial)")
  mk <- data.frame(core_id = eos_core$core_id, label_based = eos_core$eos_per_mm2, marker_based = eos_core$eos_marker_per_mm2)
  r <- suppressWarnings(stats::cor(mk$label_based, mk$marker_based, method = "spearman"))
  p2 <- ggplot(mk, aes(marker_based, label_based)) + geom_abline(linetype = 2, colour = "grey60") +
    geom_point(colour = "#C0392B", alpha = 0.8) + theme_bw(base_size = 11) +
    labs(x = sprintf("marker-based (non-epithelial cells with >= %d counts of %s) per mm2", UC_EOS_MIN_COUNTS,
                     paste(eg, collapse = "/")),
         y = "label-based eosinophils per mm2", title = sprintf("Per core: label vs markers (Spearman %.2f)", r))
  save_fig(p + p2 + plot_layout(widths = c(1, 1)), "05uc_eos_marker_check", width = 15,
           height = max(7, 0.22 * length(unique(dd$label)) + 2))
}

## B3 close-up of the eosinophil-rich cores
if (length(top_cores) && sum(m$eos) > 0) {
  ## cohort reference: median % of each label over cores
  comp <- m |> count(core_id, label) |> group_by(core_id) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
  ref <- comp |> tidyr::complete(core_id, label, fill = list(n = 0, pct = 0)) |> group_by(label) |>
    summarise(cohort_median = stats::median(pct), .groups = "drop")
  over <- list(); nb_tab <- list()
  for (cid in top_cores) {
    i <- which(m$core_id == cid)
    d <- m[i, ]
    e <- which(d$eos)
    ## neighbours of eosinophils within UC_EOS_NEIGHBOUR_UM
    nb <- RANN::nn2(d[, c("x_um", "y_um")], d[e, c("x_um", "y_um")], k = min(30, nrow(d)), searchtype = "radius",
                    radius = UC_EOS_NEIGHBOUR_UM)$nn.idx
    nb_idx <- setdiff(unique(as.vector(nb[nb > 0])), e)
    near <- seq_len(nrow(d)) %in% nb_idx
    nt <- data.frame(label = d$label) |> count(label, name = "n_core") |>
      left_join(data.frame(label = d$label[near]) |> count(label, name = "n_near"), by = "label") |>
      mutate(n_near = dplyr::coalesce(n_near, 0L), pct_core = 100 * n_core / sum(n_core),
             pct_near_eos = 100 * n_near / max(sum(n_near), 1),
             log2_enrichment = log2((pct_near_eos + 0.5) / (pct_core + 0.5)), core_id = cid)
    nb_tab[[cid]] <- nt
    cn <- comp[comp$core_id == cid, ] |> left_join(ref, by = "label")
    ttl <- sprintf("%s  -  Nancy %s  -  %d eosinophils (%.0f / mm2)", cid, core_tab[cid, "nancy"], length(e),
                   eos_core[cid, "eos_per_mm2"])
    pa <- plot_core_black(d, "label", lab_pal_dark, title = "Cell types (eosinophils ringed)", highlight = d$eos,
                          highlight_col = "#FF3EA5", panel_in = 6, legend = FALSE)
    d$show <- ifelse(d$eos, "Eosinophil", ifelse(near, d$label, "other"))
    pal_b <- lab_pal_dark; pal_b["Eosinophil"] <- "#FF3EA5"
    pb <- plot_core_black(d, "show", pal_b, title = sprintf("Eosinophils and their neighbours (< %d um)", UC_EOS_NEIGHBOUR_UM),
                          dim = !(d$eos | near), panel_in = 6) + theme(legend.text = element_text(size = 9))
    cn <- cn[cn$pct >= 1 | cn$cohort_median >= 1, ]
    cn$label <- factor(cn$label, levels = cn$label[order(cn$pct)])
    pc_ <- ggplot(cn, aes(pct, label, fill = as.character(label))) + geom_col(width = 0.75) +
      geom_point(aes(x = cohort_median), shape = 124, size = 4.5, colour = "white") +
      scale_fill_manual(values = lab_pal_dark, guide = "none") + theme_black_chart(11) +
      labs(x = "% of the core's cells", y = NULL, title = "Composition (white tick = median of all cores)")
    nt2 <- nt[nt$n_core >= 5 & nt$label != unique(d$label[e])[1], ]
    nt2$label <- factor(nt2$label, levels = nt2$label[order(nt2$log2_enrichment)])
    pd <- ggplot(nt2, aes(log2_enrichment, label, fill = as.character(label))) + geom_col(width = 0.75) +
      geom_vline(xintercept = 0, colour = "grey70") + scale_fill_manual(values = lab_pal_dark, guide = "none") +
      theme_black_chart(11) +
      labs(x = "log2 (% among eosinophil neighbours / % in the core)", y = NULL,
           title = "Who sits next to the eosinophils?")
    fig <- (pa | pb) / (pc_ | pd) + plot_layout(heights = c(1.6, 1)) +
      patchwork::plot_annotation(title = ttl, theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                                                            plot.title = element_text(colour = "white", face = "bold", size = 16)))
    save_fig(fig, paste0("05uc_eos_core_", make.names(cid)), width = 16, height = 15, dpi = 200)
    over[[cid]] <- plot_core_black(d, "label", lab_pal_dark, title = sprintf("%s (Nancy %s)", cid, core_tab[cid, "nancy"]),
                                   subtitle = sprintf("%d eosinophils, %.0f / mm2", length(e), eos_core[cid, "eos_per_mm2"]),
                                   highlight = d$eos, panel_in = 4.5, legend = FALSE)
  }
  save_table(dplyr::bind_rows(nb_tab), "05uc_eos_neighbours_top_cores")
  leg <- legend_black(plot_core_black(m[m$core_id %in% top_cores, ], "label", lab_pal_dark))
  fig <- patchwork::wrap_plots(over, ncol = 3)
  if (!is.null(leg)) fig <- patchwork::wrap_plots(fig, leg, widths = c(12, 2.5))
  fig <- fig + patchwork::plot_annotation(
    title = sprintf("The %d cores with the most eosinophils per mm2 (eosinophils ringed in pink)", length(over)),
    theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                  plot.title = element_text(colour = "white", face = "bold", size = 16)))
  save_fig(fig, "05uc_eos_top_cores", width = 18, height = 5.2 * ceiling(length(over) / 3) + 1, dpi = 200)
}
} ## end eos
message("Done.")
