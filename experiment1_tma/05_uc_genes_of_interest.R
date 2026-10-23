## =============================================================================
## 05_uc_genes_of_interest.R
##
## UC biopsies (after 04 / 04d): a closer look at a few genes, default CXCL1,
## KRT7 and SLC6A14 (EXP1_UC_GENES=KRT7,CXCL1 or UC_GENES in config.R).
##
## 1. WHICH CELLS express them, and are they the SAME cells?
##    1a  where each gene's transcripts are: share per cell type (one bar per gene)
##        + % positive cells and transcripts per positive cell, per cell type
##    1b  co-expression: how many cells carry one, two or all genes, in which
##        compartment; "of the KRT7+ cells, % also CXCL1+" vs what is expected
##        from the cell types alone (fold > 1 = the genes go together in the same
##        cells); do cores high in one gene also have the others?
##
## 2. HOW DO THEY CHANGE with Nancy and Disease.Control -- more cells, or more per cell?
##    Per core, for each gene and its main cell type(s) (>= 10% of its transcripts):
##      level        transcripts of the gene per 1,000 transcripts of the core
##      abundance    % of the core's cells that are this cell type
##      prevalence   % of this cell type's cells that are positive   -> "more cells"
##      intensity    transcripts per positive cell                    -> "more per cell"
##    vs Nancy (cores: Spearman, mixed model with patient) and vs Disease.Control
##    (patients: Wilcoxon on the patient means; cores: mixed model with patient,
##    with and without Nancy -- Disease.Control and Nancy are linked, so "adjusted
##    for Nancy" says whether the gene tells you something beyond the Nancy score).
##    Summary tile plot: one row per measure, red = higher with Nancy / in
##    Disease.Control yes, blue = lower.
##
## 3. FOUR EXAMPLE CORES on black: one Nancy 0, one Nancy 4, one Disease.Control
##    yes, one no -- picked automatically as round, solid cores with both
##    epithelium and immune cells (tables/05uc_goi_example_core_candidates.csv lists
##    alternatives; UC_GOI_EXAMPLE_CORES in config.R to choose your own, e.g.
##    c(nancy0 = "A_1-a", nancy4 = "C_3-f", dc_yes = "D_2-b", dc_no = "A_5-c")).
##    Rows: cell types, each gene, and which gene(s) each cell carries.
##
## INPUT : xenium_<seg>_annotated.rds (04 / 04d), tma_map Nancy / Nancy_TMA,
##         de-identified clinical table with Disease.Control (UC_CLINICAL_CSV)
## OUTPUT: tables/05uc_goi_*.csv, figures/05uc_goi_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_GOI <- trimws(strsplit(Sys.getenv("EXP1_UC_GENES", paste(get0("UC_GENES", ifnotfound = c("CXCL1", "KRT7", "SLC6A14")),
                                                           collapse = ",")), ",")[[1]])
UC_DC_COL          <- get0("UC_DC_COL", ifnotfound = "Disease.Control")
UC_GOI_MAIN_SHARE  <- get0("UC_GOI_MAIN_SHARE", ifnotfound = 10)   # main cell type: >= 10% of the gene's transcripts
UC_GOI_MAX_TYPES   <- get0("UC_GOI_MAX_TYPES", ifnotfound = 2)     # at most this many main cell types per gene
UC_GOI_EXAMPLE_CORES <- get0("UC_GOI_EXAMPLE_CORES", ifnotfound = NULL)
GENE_COLOURS <- c("#FFD23F", "#3EC1FF", "#FF3EA5", "#7CFC6B", "#FF8C42")   # single genes on black

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
rm(xen); invisible(gc())
core_tab <- uc_core_table(m)
core_tab$grade <- factor(ifelse(is.na(core_tab$nancy), "NA", paste("Nancy", core_tab$nancy)),
                         levels = c(paste("Nancy", 0:4), "NA"))
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
clin <- join_uc_clinical(core_tab, read_uc_clinical())
core_tab$dc <- factor(NA, levels = c("no", "yes"))
if (!is.null(clin) && UC_DC_COL %in% colnames(clin)) {
  v <- as_outcome(clin[core_tab$core_id, UC_DC_COL])
  if (nlevels(v) == 2) core_tab$dc <- factor(as.character(v), levels = levels(v)) else
    message("!! ", UC_DC_COL, " is not a yes/no column (values: ", paste(utils::head(levels(v), 6), collapse = ", "), ")")
} else message("!! No ", UC_DC_COL, " (clinical table missing or column not found) -- Disease.Control parts are skipped.")
has_dc <- any(!is.na(core_tab$dc)) && length(unique(stats::na.omit(core_tab$dc))) == 2
dc_cols <- stats::setNames(c("#8E8E8E", "#2E8B57"), levels(core_tab$dc))
message("Cores: ", nrow(core_tab), ", patients: ", length(unique(core_tab$patient_id)),
        if (has_dc) paste0(" | ", UC_DC_COL, ": ", paste(names(table(core_tab$dc[!duplicated(core_tab$patient_id)])),
                                                         table(core_tab$dc[!duplicated(core_tab$patient_id)]),
                                                         sep = " = ", collapse = ", "), " patients") else "")

g_obj <- match_case(UC_GOI, rownames(counts))
miss <- UC_GOI[!g_obj %in% rownames(counts)]
if (length(miss)) message("Not on the panel (left out): ", paste(miss, collapse = ", "))
G <- stats::setNames(g_obj[g_obj %in% rownames(counts)], UC_GOI[g_obj %in% rownames(counts)])
if (!length(G)) stop("None of ", paste(UC_GOI, collapse = ", "), " is on the panel.")
X <- as.matrix(Matrix::t(counts[G, , drop = FALSE])); colnames(X) <- names(G)
P <- X > 0
lab_pal <- celltype_palette(m$label, m$lineage)
lab_pal_dark <- lighten(lab_pal, 0.25)
gene_pal <- stats::setNames(GENE_COLOURS[seq_along(G)], names(G))

## =============================================================================
## 1. which cells
## =============================================================================
ct_tab <- dplyr::bind_rows(lapply(names(G), function(g) {
  data.frame(label = m$label, lineage = m$lineage, compartment = m$compartment, x = X[, g]) |>
    group_by(label, lineage, compartment) |>
    summarise(n_cells = n(), n_positive = sum(x > 0), pct_positive = 100 * mean(x > 0),
              transcripts_per_positive_cell = ifelse(sum(x > 0) > 0, sum(x) / sum(x > 0), NA),
              transcripts = sum(x), .groups = "drop") |>
    mutate(gene = g, pct_of_gene_transcripts = 100 * transcripts / max(sum(transcripts), 1),
           pct_of_positive_cells = 100 * n_positive / max(sum(n_positive), 1))
}))
save_table(ct_tab |> arrange(gene, desc(pct_of_gene_transcripts)), "05uc_goi_by_cell_type")
main_types <- lapply(stats::setNames(names(G), names(G)), function(g) {
  d <- ct_tab |> filter(gene == g, n_cells >= 50) |> arrange(desc(pct_of_gene_transcripts))
  mt <- d$label[d$pct_of_gene_transcripts >= UC_GOI_MAIN_SHARE]
  utils::head(if (length(mt)) mt else d$label[1], UC_GOI_MAX_TYPES)
})
for (g in names(G)) {
  d <- ct_tab |> filter(gene == g) |> arrange(desc(pct_of_gene_transcripts)) |> head(4)
  message(sprintf("%-8s %s", g, paste(sprintf("%s %.0f%%", d$label, d$pct_of_gene_transcripts), collapse = " | ")))
}

## 1a. one bar per gene: share of its transcripts per cell type
seg <- ct_tab |> group_by(gene) |> arrange(desc(pct_of_gene_transcripts), .by_group = TRUE) |>
  mutate(seg = ifelse(row_number() <= 7 & pct_of_gene_transcripts >= 2, label, "other cell types")) |>
  group_by(gene, seg) |> summarise(pct = sum(pct_of_gene_transcripts), .groups = "drop")
seg_lv <- c(intersect(names(lab_pal), seg$seg), "other cell types")
seg$seg <- factor(seg$seg, levels = rev(seg_lv))
seg$gene <- factor(seg$gene, levels = rev(names(G)))
pal_seg <- c(lab_pal, "other cell types" = "grey85")
p1 <- ggplot(seg, aes(pct, gene, fill = seg)) + geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 12, sprintf("%s\n%.0f%%", vapply(gsub("_", " ", seg), function(z)
              paste(strwrap(z, 16), collapse = "\n"), character(1)), pct), "")), position = position_stack(vjust = 0.5),
            size = 3, lineheight = 0.85) +
  scale_fill_manual(values = pal_seg, name = NULL, breaks = seg_lv) + scale_x_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 13) + theme(panel.grid = element_blank(), legend.position = "right") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "% of the gene's transcripts", y = NULL, title = "A. Which cell types carry each gene's transcripts?")
keep_lab <- ct_tab |> group_by(label) |> filter(max(n_cells) >= 50, max(pct_positive) >= 1) |> pull(label) |> unique()
dd <- ct_tab[ct_tab$label %in% keep_lab, ]
dd$label <- factor(dd$label, levels = rev(intersect(names(lab_pal), keep_lab)))
dd$gene <- factor(dd$gene, levels = names(G))
p2 <- ggplot(dd, aes(gene, label, size = pct_positive, colour = transcripts_per_positive_cell)) + geom_point() +
  scale_size_area(max_size = 11, name = "% cells\npositive") +
  scale_colour_gradient(low = "#C9C9C9", high = "#B2182B", name = "transcripts per\npositive cell") +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = NULL, title = "B. Per cell type",
       subtitle = "dot size = how many of its cells express the gene;\ncolour = how much each expressing cell carries")
save_fig(p1 / p2 + plot_layout(heights = c(1, 2.2)), "05uc_goi_cell_types", width = 13,
         height = max(10, 0.28 * length(keep_lab) + 5))

## 1b. co-expression
anyp <- rowSums(P) > 0
combo <- apply(P[anyp, , drop = FALSE], 1, function(r) paste(colnames(P)[r], collapse = " + "))
cb <- data.frame(combo = combo, compartment = m$compartment[anyp], label = m$label[anyp]) |>
  count(combo, compartment, name = "n_cells")
cb_lab <- data.frame(combo = combo, label = m$label[anyp]) |> count(combo, label) |> group_by(combo) |>
  arrange(desc(n), .by_group = TRUE) |>
  summarise(n_cells = sum(n), top_cell_types = paste(sprintf("%s %.0f%%", head(label, 3), 100 * head(n, 3) / sum(n)), collapse = "; "),
            .groups = "drop") |> arrange(desc(n_cells))
save_table(cb_lab, "05uc_goi_coexpression_combinations")
cb$combo <- factor(cb$combo, levels = cb_lab$combo)
cb$n_genes <- lengths(strsplit(as.character(cb$combo), " + ", fixed = TRUE))
comp_pal <- c(Epithelium = "#418ED0", Immune = "#C0392B", Stroma = "#4DAF4A", Other = "#BDBDBD")
q1 <- ggplot(cb, aes(combo, n_cells, fill = compartment)) + geom_col(width = 0.75) +
  geom_text(data = cb_lab, aes(combo, n_cells, label = format(n_cells, big.mark = ",")), vjust = -0.3, size = 3.2,
            inherit.aes = FALSE) +
  scale_fill_manual(values = comp_pal, name = NULL) + scale_y_continuous(expand = expansion(c(0, 0.08))) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top") +
  labs(x = "genes detected in the cell", y = "cells", title = "C. How many cells carry one, two or all genes?")
## of the A+ cells, % also B+ -- vs expected from their cell types alone
cond <- dplyr::bind_rows(lapply(names(G), function(a) dplyr::bind_rows(lapply(names(G), function(b) {
  if (a == b) return(NULL)
  sa <- P[, a]
  obs <- 100 * mean(P[sa, b])
  pb <- tapply(P[, b], m$label, mean)                           # B+ rate per cell type
  expd <- 100 * mean(pb[m$label[sa]])                           # same cell types, genes independent
  data.frame(gene_a = a, gene_b = b, n_a_positive = sum(sa), pct_a_cells_also_b = obs,
             pct_expected_from_cell_types = expd, fold = obs / max(expd, 1e-6))
}))))
save_table(cond, "05uc_goi_coexpression_conditional")
q2 <- if (nrow(cond)) ggplot(cond, aes(gene_b, gene_a, fill = log2(fold))) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%\n(x%.1f)", pct_a_cells_also_b, fold)), size = 3.6, lineheight = 0.9) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", midpoint = 0, limits = c(-2, 2),
                       oob = scales::squish, name = "log2 fold\nvs expected") +
  theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
  labs(x = "... also positive for", y = "cells positive for ...",
       title = "D. Same cells?",
       subtitle = "% of row-gene+ cells that also carry the column gene\n(x fold = vs expected from their cell types alone)") else NULL
## across cores: are cores high in one gene also high in the others?
lvl <- sapply(names(G), function(g) 1000 * tapply(X[, g], m$core_id, sum) / tapply(m$nCount_Xenium, m$core_id, sum))
lvl <- matrix(lvl, ncol = length(G), dimnames = list(sort(unique(m$core_id)), names(G)))
cc <- if (length(G) > 1) stats::cor(lvl, method = "spearman") else NULL
q3 <- if (!is.null(cc)) {
  cd <- as.data.frame(as.table(cc)); colnames(cd) <- c("a", "b", "rho")
  ggplot(cd, aes(b, a, fill = rho)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 4) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", limits = c(-1, 1), name = "Spearman") +
    theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, title = "E. Across cores", subtitle = "per-core levels (per 1,000 transcripts)")
} else NULL
save_fig(q1 / (patchwork::wrap_plots(Filter(Negate(is.null), list(q2, q3)), nrow = 1)) + plot_layout(heights = c(1, 1)),
         "05uc_goi_coexpression", width = 14, height = 11)

## =============================================================================
## 2. per core: level, abundance, prevalence, intensity -- vs Nancy and Disease.Control
## =============================================================================
mlab <- c(level = "transcripts per 1,000", abundance = "% of core cells", prevalence = "% positive",
          intensity = "transcripts per positive cell")
core_n <- table(m$core_id)
long <- dplyr::bind_rows(lapply(names(G), function(g) {
  d <- data.frame(core_id = m$core_id, label = m$label, x = X[, g], tot = m$nCount_Xenium)
  lv <- d |> group_by(core_id) |> summarise(value = 1000 * sum(x) / sum(tot), .groups = "drop") |>
    mutate(gene = g, cell_type = "all cells", metric = "level")
  per_ct <- dplyr::bind_rows(lapply(main_types[[g]], function(ct) {
    d |> group_by(core_id) |>
      summarise(n_ct = sum(label == ct), n_pos = sum(label == ct & x > 0),
                abundance = 100 * n_ct / n(),
                prevalence = ifelse(n_ct >= UC_MIN_DENOM, 100 * n_pos / n_ct, NA),
                intensity = ifelse(n_pos >= 5, sum(x[label == ct & x > 0]) / n_pos, NA), .groups = "drop") |>
      tidyr::pivot_longer(c(abundance, prevalence, intensity), names_to = "metric", values_to = "value") |>
      transmute(core_id, value, gene = g, cell_type = ct, metric)
  }))
  dplyr::bind_rows(lv, per_ct)
}))
long$measure <- ifelse(long$cell_type == "all cells", sprintf("%s: %s (all cells)", long$gene, mlab[long$metric]),
                       sprintf("%s: %s %s", long$gene, ifelse(long$metric == "abundance", "", paste0(mlab[long$metric], " in")),
                               long$cell_type))
long$measure[long$metric == "abundance"] <- sprintf("%s: %% of core cells that are %s", long$gene[long$metric == "abundance"],
                                                     long$cell_type[long$metric == "abundance"])
long <- cbind(long, core_tab[long$core_id, c("patient_id", "nancy", "grade", "dc")])
measures <- unique(long$measure)
save_table(long, "05uc_goi_per_core_long")
wide <- tidyr::pivot_wider(long[, c("core_id", "measure", "value")], names_from = measure, values_from = value) |> as.data.frame()
wide <- cbind(core_tab[wide$core_id, c("core_id", "patient_id", "nancy", "dc")], wide[, measures, drop = FALSE])
save_table(wide, "05uc_goi_per_core")

## statistics
st_n <- if (any(!is.na(wide$nancy))) assoc_nancy(wide, measures) else data.frame()
lrt <- function(y, full, null, pat) {
  ok <- stats::complete.cases(y, full, pat); y <- y[ok]; full <- full[ok, , drop = FALSE]; pat <- pat[ok]
  null <- if (is.null(null)) NULL else null[ok, , drop = FALSE]
  if (length(y) < 8 || stats::sd(y) == 0) return(c(est = NA, p = NA))
  dfull <- data.frame(y = y, full, pat = pat)
  rhs_f <- paste(colnames(full), collapse = " + "); rhs_n <- if (is.null(null)) "1" else paste(colnames(null), collapse = " + ")
  suppressWarnings(tryCatch({
    if (any(duplicated(pat)) && requireNamespace("lme4", quietly = TRUE)) {
      f1 <- suppressMessages(lme4::lmer(stats::as.formula(paste("y ~", rhs_f, "+ (1 | pat)")), data = dfull, REML = FALSE))
      f0 <- suppressMessages(lme4::lmer(stats::as.formula(paste("y ~", rhs_n, "+ (1 | pat)")), data = dfull, REML = FALSE))
      c(est = unname(lme4::fixef(f1)["dcyes"]), p = stats::anova(f0, f1)[["Pr(>Chisq)"]][2])
    } else {
      f1 <- stats::lm(stats::as.formula(paste("y ~", rhs_f)), data = dfull)
      f0 <- stats::lm(stats::as.formula(paste("y ~", rhs_n)), data = dfull)
      c(est = unname(stats::coef(f1)["dcyes"]), p = stats::anova(f0, f1)[["Pr(>F)"]][2])
    }
  }, error = function(e) c(est = NA, p = NA)))
}
st_dc <- data.frame()
if (has_dc) {
  pat <- wide |> filter(!is.na(dc)) |> group_by(patient_id) |>
    summarise(dc = first(dc), nancy_max = suppressWarnings(max(nancy, na.rm = TRUE)),
              dplyr::across(dplyr::all_of(measures), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |> as.data.frame()
  pat$nancy_max[!is.finite(pat$nancy_max)] <- NA
  for (v in measures) pat[[v]][is.nan(pat[[v]])] <- NA
  save_table(pat, "05uc_goi_per_patient")
  st_dc <- dplyr::bind_rows(lapply(measures, function(v) {
    yes <- pat[[v]][pat$dc == "yes"]; no <- pat[[v]][pat$dc == "no"]
    yes <- yes[!is.na(yes)]; no <- no[!is.na(no)]
    if (length(yes) < 3 || length(no) < 3) return(NULL)
    dcore <- wide[!is.na(wide$dc), ]
    a <- lrt(dcore[[v]], data.frame(dc = dcore$dc), NULL, dcore$patient_id)
    b <- lrt(dcore[[v]], data.frame(dc = dcore$dc, nancy = dcore$nancy), data.frame(nancy = dcore$nancy), dcore$patient_id)
    data.frame(measure = v, n_patients_yes = length(yes), n_patients_no = length(no),
               median_yes = stats::median(yes), median_no = stats::median(no),
               auc_yes_higher = auc_mw(yes, no),
               p_patients_wilcox = suppressWarnings(stats::wilcox.test(yes, no, exact = FALSE))$p.value,
               diff_cores_mixed = a[["est"]], p_cores_mixed = a[["p"]],
               diff_adjusted_for_nancy = b[["est"]], p_adjusted_for_nancy = b[["p"]])
  }))
  if (nrow(st_dc)) for (c in c("p_patients_wilcox", "p_cores_mixed", "p_adjusted_for_nancy"))
    st_dc[[sub("^p_", "padj_", c)]] <- stats::p.adjust(st_dc[[c]], "BH")
  save_table(st_dc, "05uc_goi_vs_disease_control")
}
if (nrow(st_n)) save_table(st_n, "05uc_goi_vs_nancy")

## summary tile plot
tiles <- list()
if (nrow(st_n)) tiles$n <- data.frame(measure = st_n$feature, test = "Nancy\n(cores, Spearman)",
                                      eff = st_n$spearman_rho, p = st_n$p_spearman,
                                      txt = sprintf("rho %+.2f\np %.2g", st_n$spearman_rho, st_n$p_spearman))
if (nrow(st_dc)) {
  tiles$d1 <- data.frame(measure = st_dc$measure, test = "Disease control yes vs no\n(patients, Wilcoxon)",
                         eff = 2 * (st_dc$auc_yes_higher - 0.5), p = st_dc$p_patients_wilcox,
                         txt = sprintf("%.3g vs %.3g\np %.2g", st_dc$median_yes, st_dc$median_no, st_dc$p_patients_wilcox))
  tiles$d2 <- data.frame(measure = st_dc$measure, test = "Disease control yes vs no,\nsame Nancy (cores, mixed model)",
                         eff = sign(st_dc$diff_adjusted_for_nancy), p = st_dc$p_adjusted_for_nancy,
                         txt = sprintf("%+.3g\np %.2g", st_dc$diff_adjusted_for_nancy, st_dc$p_adjusted_for_nancy))
}
tiles <- dplyr::bind_rows(tiles)
if (nrow(tiles)) {
  tiles$score <- sign(tiles$eff) * pmin(-log10(pmax(tiles$p, 1e-10)), 4)
  tiles$measure <- factor(tiles$measure, levels = rev(measures))
  tiles$test <- factor(tiles$test, levels = unique(tiles$test))
  tiles$gene <- sub(":.*$", "", as.character(tiles$measure))
  pt <- ggplot(tiles, aes(test, measure, fill = score)) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = txt), size = 3.1, lineheight = 0.85) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-4, 4),
                         name = "direction x\n-log10 p", breaks = c(-4, -2, 0, 2, 4),
                         labels = c("lower, p<1e-4", "lower, p 0.01", "no change", "higher, p 0.01", "higher, p<1e-4")) +
    scale_x_discrete(position = "top") + facet_grid(gene ~ ., scales = "free_y", space = "free_y") +
    theme_bw(base_size = 12) + theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold")) +
    labs(x = NULL, y = NULL, title = "Genes of interest: what changes with Nancy and with Disease control?",
         subtitle = paste("red = higher with Nancy / in Disease control 'yes', blue = lower.",
                          "\n'% positive' up = expressed in MORE cells; 'transcripts per positive cell' up = MORE per cell;",
                          "\n'% of core cells that are X' up = just more of the expressing cell type."))
  save_fig(pt, "05uc_goi_summary", width = 6 + 3.6 * length(unique(tiles$test)), height = 1.2 + 0.62 * length(measures))
}

## per gene: every measure by Nancy grade (cores) and by Disease control (patients)
for (g in names(G)) {
  ms <- measures[startsWith(measures, paste0(g, ":"))]
  d <- long[long$gene == g & !is.na(long$value), ]
  d$measure <- factor(d$measure, levels = ms)
  sn <- st_n[st_n$feature %in% ms, ]
  pa <- ggplot(d, aes(grade, value, colour = grade)) + geom_boxplot(outlier.shape = NA, colour = "grey50") +
    geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1.3) +
    scale_colour_manual(values = grade_cols, guide = "none") +
    facet_wrap(~ measure, ncol = 1, scales = "free_y", labeller = label_wrap_gen(55)) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = NULL, title = "by Nancy grade (one dot = one core)")
  if (nrow(sn)) pa <- pa + geom_text(data = data.frame(measure = factor(sn$feature, levels = ms),
                                                         txt = sprintf("rho %.2f, p %.2g", sn$spearman_rho, sn$p_spearman)),
                                       aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.3, size = 3,
                                       inherit.aes = FALSE)
  pl <- pa
  if (has_dc) {
    dp <- tidyr::pivot_longer(pat[, c("patient_id", "dc", "nancy_max", ms)], dplyr::all_of(ms), names_to = "measure",
                              values_to = "value")
    dp <- dp[!is.na(dp$value), ]
    dp$measure <- factor(dp$measure, levels = ms)
    dp$nancy_group <- nancy_group(dp$nancy_max)
    sd <- st_dc[st_dc$measure %in% ms, ]
    pb <- ggplot(dp, aes(dc, value)) + geom_boxplot(aes(fill = dc), outlier.shape = NA, alpha = 0.3, colour = "grey45") +
      geom_point(aes(colour = nancy_group), position = position_jitter(width = 0.15, seed = 1), size = 1.8) +
      scale_fill_manual(values = dc_cols, guide = "none") +
      scale_colour_manual(values = GROUP_COLOURS(), name = "patient's\nhighest Nancy", na.value = "grey60") +
      facet_wrap(~ measure, ncol = 1, scales = "free_y", labeller = label_wrap_gen(55)) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(x = UC_DC_COL, y = NULL, title = "by Disease control (one dot = one patient, mean of its cores)")
    if (nrow(sd)) pb <- pb + geom_text(data = data.frame(measure = factor(sd$measure, levels = ms),
                                                           txt = sprintf("patients p %.2g | same Nancy p %.2g",
                                                                         sd$p_patients_wilcox, sd$p_adjusted_for_nancy)),
                                         aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.3, size = 3,
                                         inherit.aes = FALSE)
    pl <- pa | pb
  }
  save_fig(pl + plot_annotation(title = paste(g, "-- more cells, or more per cell?"),
                                subtitle = paste0("main cell type(s): ", paste(main_types[[g]], collapse = ", "),
                                                  "  (>= ", UC_GOI_MAIN_SHARE, "% of the gene's transcripts)")),
           paste0("05uc_goi_", make.names(g), "_nancy_disease_control"), width = if (has_dc) 13 else 7,
           height = 2.6 * length(ms) + 1.5)
}

## =============================================================================
## 3. four example cores on black
## =============================================================================
shape <- dplyr::bind_rows(lapply(split(seq_len(nrow(m)), m$core_id), function(i) {
  xy <- as.matrix(m[i, c("x_um", "y_um")])
  h <- xy[grDevices::chull(xy), , drop = FALSE]
  hx <- c(h[, 1], h[1, 1]); hy <- c(h[, 2], h[1, 2])
  A <- abs(sum(hx[-1] * hy[-length(hy)] - hx[-length(hx)] * hy[-1])) / 2
  Pm <- sum(sqrt(diff(hx)^2 + diff(hy)^2))
  data.frame(core_id = m$core_id[i[1]], n_cells = length(i), roundness = 4 * pi * A / max(Pm^2, 1),
             solidity = min(1, core_area_mm2(xy[, 1], xy[, 2]) * 1e6 / max(A, 1)),
             epithelium_pct = 100 * mean(m$compartment[i] == "Epithelium"),
             immune_pct = 100 * mean(m$compartment[i] == "Immune"),
             goi_positive_pct = 100 * mean(rowSums(P[i, , drop = FALSE]) > 0))
}))
shape$score <- with(shape, roundness * solidity * pmin(epithelium_pct / 35, 1) * pmin(immune_pct / 15, 1) *
                      pmin(n_cells / 1500, 1))
shape <- cbind(shape, core_tab[shape$core_id, c("patient_id", "nancy", "dc")])
shape <- shape[order(-shape$score), ]
lo <- suppressWarnings(min(shape$nancy, na.rm = TRUE)); hi <- suppressWarnings(max(shape$nancy, na.rm = TRUE))
cats <- list(nancy0 = shape$nancy %in% lo, nancy4 = shape$nancy %in% hi,
             dc_yes = shape$dc %in% "yes", dc_no = shape$dc %in% "no")
cand <- dplyr::bind_rows(lapply(names(cats), function(k) if (any(cats[[k]])) cbind(category = k, utils::head(shape[cats[[k]], ], 5))))
save_table(cand, "05uc_goi_example_core_candidates")
pick <- c()
for (k in names(cats)) {
  if (!is.null(UC_GOI_EXAMPLE_CORES) && !is.na(UC_GOI_EXAMPLE_CORES[k])) { pick[k] <- UC_GOI_EXAMPLE_CORES[[k]]; next }
  cc <- shape$core_id[cats[[k]] & !shape$core_id %in% pick]
  if (length(cc)) pick[k] <- cc[1]
}
pick <- pick[!is.na(pick) & pick %in% m$core_id]
if (length(pick)) {
  head_txt <- c(nancy0 = paste("Nancy", lo), nancy4 = paste("Nancy", hi), dc_yes = "Disease control: yes",
                dc_no = "Disease control: no")
  message("Example cores: ", paste(names(pick), pick, sep = " = ", collapse = ", "))
  combo_pal <- c(gene_pal, "2+ genes" = "#FFFFFF")
  cols <- lapply(names(pick), function(k) {
    cid <- pick[[k]]; i <- which(m$core_id == cid); d <- m[i, ]
    sub_t <- sprintf("%s  |  Nancy %s  |  DC %s", cid, core_tab[cid, "nancy"], ifelse(is.na(core_tab[cid, "dc"]), "NA",
                                                                                    as.character(core_tab[cid, "dc"])))
    last <- k == utils::tail(names(pick), 1)
    out <- list(plot_core_black(d, "label", lab_pal_dark, title = head_txt[[k]], subtitle = sub_t, legend = FALSE, panel_in = 4.2))
    for (g in names(G)) {
      p <- plot_core_gene_black(d, X[i, g], g, title = g, panel_in = 4.2)
      out[[g]] <- if (last) p else p + theme(legend.position = "none")
    }
    np <- rowSums(P[i, , drop = FALSE])
    d$combo <- ifelse(np == 0, "none", ifelse(np >= 2, "2+ genes", colnames(P)[max.col(P[i, , drop = FALSE], ties.method = "first")]))
    pc <- plot_core_black(d, "combo", combo_pal, title = "which gene(s)", dim = d$combo == "none", panel_in = 4.2,
                          subtitle = sprintf("%.0f%% of cells carry >= 1 gene", 100 * mean(np > 0)))
    out$combo <- if (last) pc else pc + theme(legend.position = "none")
    out
  })
  nr <- length(cols[[1]])
  plots <- unlist(lapply(seq_len(nr), function(r) lapply(cols, `[[`, r)), recursive = FALSE)
  fig <- patchwork::wrap_plots(plots, ncol = length(pick), byrow = TRUE)
  leg <- legend_black(plot_core_black(m[m$core_id %in% pick, ], "label", lab_pal_dark))
  if (!is.null(leg)) fig <- patchwork::wrap_plots(fig, leg, widths = c(4.4 * length(pick), 2.6))
  fig <- fig + patchwork::plot_annotation(
    title = paste("Example cores:", paste(names(G), collapse = ", ")),
    subtitle = "row 1: cell types; then one row per gene (dark grey = 0 transcripts); last row: which gene(s) each cell carries",
    theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                  plot.title = element_text(colour = "white", face = "bold", size = 18),
                  plot.subtitle = element_text(colour = "grey80", size = 12)))
  save_fig(fig, "05uc_goi_example_cores", width = 4.4 * length(pick) + 3, height = 4.5 * nr + 1.2, dpi = 200)
}
message("Done.")
