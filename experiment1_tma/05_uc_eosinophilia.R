## =============================================================================
## 05_uc_eosinophilia.R
##
## UC biopsies: which cores show eosinophilia, what they look like, and how it
## relates to Nancy -- WITHOUT relying on an eosinophil cluster (04 / 04d).
##
## Eosinophil = a cell that CO-EXPRESSES SIGLEC8 and RNASE3 (>= 1 transcript of each;
## low expression is fine), and
##   - is not epithelial and not a mast cell (mast cells also express SIGLEC8),
##   - does not carry a neutrophil profile: fewer than UC_EOS_MAX_NEUT (default 2)
##     transcripts of CSF3R / S100A9 / S100A12 / TREM1 / CXCL8 / FCGR3B / CXCR2.
## Requiring both genes removes cells with a single stray transcript (background
## decoding, spill-over from a neighbouring eosinophil); the neutrophil filter removes
## co-expressing cells that are really neutrophils. UC_EOS_CALL <- "lenient" goes back to
## "SIGLEC8 or RNASE3" (more sensitive, less specific).
## figures/05uc_eosinophilia_marker_check.png shows what was kept and dropped, and the
## chance rate of the same call in T / B cells (should be ~0).
##
## Per core: eosinophils per mm^2 (and % of immune cells).
## Eosinophilia per core:
##   pathology   the clinical eosinophil column (UC_CLIN_EOS_COLS), if it is a yes / no call
##   TMA         eosinophils per mm^2 above Q3 + 1.5 IQR of all cores (or the top 10%)
## UC_EOSINOPHILIA_BY = "pathology" / "tma" / "auto" (pathology when available) picks the
## call used for the core images.
##
## Figures
##   05uc_eosinophilia_marker_check.png  called eosinophils vs dropped cells vs neutrophils vs
##                                       T / B cells, on eosinophil / neutrophil / mast genes
##   05uc_eosinophilia_vs_nancy.png      eosinophils per mm^2 by Nancy grade (cores; Spearman,
##                                       mixed model with patient), % of cores with eosinophilia per
##                                       Nancy grade, TMA vs the pathology call
##   05uc_eosinophilia_cores.png         the eosinophilia cores on black, ordered by Nancy:
##                                       tissue (epithelium / immune / stroma) + eosinophils in pink
##   05uc_eosinophilia_controls.png      for comparison: low-eosinophil cores of the same Nancy grades
##
## INPUT : xenium_<seg>_annotated.rds (04), tma_map Nancy / Nancy_TMA, clinical table (optional)
## OUTPUT: tables/05uc_eosinophilia_*.csv, figures/05uc_eosinophilia_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_EOS_CALL         <- get0("UC_EOS_CALL", ifnotfound = "strict")      # "strict" (SIGLEC8 AND RNASE3) or "lenient" (OR)
UC_EOS_MAX_NEUT     <- get0("UC_EOS_MAX_NEUT", ifnotfound = 2)         # >= this many neutrophil-gene transcripts = neutrophil
UC_EOS_NEUT_GENES   <- get0("UC_EOS_NEUT_GENES", ifnotfound = c("CSF3R", "S100A9", "S100A12", "TREM1", "CXCL8", "FCGR3B", "CXCR2"))
UC_EOSINOPHILIA_BY  <- get0("UC_EOSINOPHILIA_BY", ifnotfound = "auto")
UC_EOS_MAX_CORES    <- get0("UC_EOS_MAX_CORES", ifnotfound = 12)       # cores shown per figure
COMP_BLACK <- c(Epithelium = "#2F5F8F", Immune = "#8C3B3B", Stroma = "#3C7A4A", Other = "#4A4A4A")
EOS_PINK <- "#FF3EA5"

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
rm(xen); invisible(gc())
on_panel <- function(g) { g <- match_case(g, rownames(counts)); g[g %in% rownames(counts)] }
gs <- on_panel("SIGLEC8"); gr <- on_panel("RNASE3"); ng <- on_panel(UC_EOS_NEUT_GENES)
mg <- on_panel(c("CPA3", "MS4A2", "HDC", "KIT")); ga <- on_panel("ALOX15")
if (!length(gs) && !length(gr)) stop("Neither SIGLEC8 nor RNASE3 is on the panel.")
cnt_of <- function(g) if (length(g)) Matrix::colSums(counts[g, , drop = FALSE]) else rep(0, ncol(counts))
s8 <- cnt_of(gs); r3 <- cnt_of(gr); neut <- cnt_of(ng)
m$eligible <- m$compartment %in% c("Immune", "Stroma", "Other") & !m$class %in% "Mast"
coexp <- if (length(gs) && length(gr)) s8 >= 1 & r3 >= 1 else (s8 + r3) >= 2
if (!(length(gs) && length(gr))) message("Only one of SIGLEC8 / RNASE3 on the panel -> 'co-expression' = >= 2 transcripts of it.")
m$eos_any    <- m$eligible & (s8 + r3) >= 1
m$eos_coexp  <- m$eligible & coexp
m$neut_profile <- neut >= UC_EOS_MAX_NEUT
m$eos_strict <- m$eos_coexp & !m$neut_profile
m$eos <- if (UC_EOS_CALL == "lenient") m$eos_any else m$eos_strict
m$eos_tx <- s8 + r3
message(sprintf("Eosinophil call '%s': %s cells. SIGLEC8 or RNASE3: %s | both: %s | both but neutrophil profile (dropped): %s",
                UC_EOS_CALL, format(sum(m$eos), big.mark = ","), format(sum(m$eos_any), big.mark = ","),
                format(sum(m$eos_coexp), big.mark = ","), format(sum(m$eos_coexp & m$neut_profile), big.mark = ",")))

## chance rate in T / B cells: the same rule where there are no eosinophils
tb <- grepl("t_cell|cd4|cd8|treg|tfh|b_cell|naive|plasma", tolower(m$label)) & m$compartment == "Immune"
bg <- data.frame(group = c("eligible cells", "T / B / plasma cells (background)"),
                 pct_SIGLEC8_or_RNASE3 = c(100 * mean((s8 + r3 >= 1)[m$eligible]), 100 * mean((s8 + r3 >= 1)[tb])),
                 pct_both = c(100 * mean(coexp[m$eligible]), 100 * mean(coexp[tb])),
                 pct_strict_call = c(100 * mean((coexp & !m$neut_profile)[m$eligible]), 100 * mean((coexp & !m$neut_profile)[tb])))
save_table(bg, "05uc_eosinophilia_call_background")
message("Same rules in T / B / plasma cells (chance rate; should be far below the eligible cells):")
print(bg, row.names = FALSE, digits = 3)

core_tab <- uc_core_table(m)
core_tab$grade <- factor(ifelse(is.na(core_tab$nancy), "NA", paste("Nancy", core_tab$nancy)),
                         levels = c(paste("Nancy", 0:4), "NA"))
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
clin <- join_uc_clinical(core_tab, read_uc_clinical())

## -----------------------------------------------------------------------
## 0. marker check: what was kept, what was dropped
## -----------------------------------------------------------------------
grp <- dplyr::case_when(m$eos ~ "eosinophils (called)",
                        m$eos_coexp & m$neut_profile ~ "both genes, but neutrophil profile (dropped)",
                        m$eos_any & !m$eos_coexp ~ "only one of the two genes (dropped)",
                        m$class %in% "Neutrophil" ~ "neutrophil-labelled cells",
                        tb ~ "T / B / plasma cells (background)", TRUE ~ NA_character_)
grp_lv <- c("eosinophils (called)", "both genes, but neutrophil profile (dropped)", "only one of the two genes (dropped)",
            "neutrophil-labelled cells", "T / B / plasma cells (background)")
keep <- which(!is.na(grp))
if (length(keep) > 3e5) keep <- sort(c(which(m$eos), sample(setdiff(keep, which(m$eos)), 3e5)))
genes_show <- c(gs, gr, ga, ng, mg)
dd <- dplyr::bind_rows(lapply(genes_show, function(g) {
  data.frame(group = grp[keep], x = counts[g, keep]) |> group_by(group) |>
    summarise(n = n(), pct = 100 * mean(x > 0), mean_counts = mean(x), .groups = "drop") |> mutate(gene = g)
}))
dd$group <- factor(dd$group, levels = rev(grp_lv))
dd$gene <- factor(dd$gene, levels = genes_show)
dd$set <- factor(ifelse(dd$gene %in% c(gs, gr, ga), "eosinophil", ifelse(dd$gene %in% ng, "neutrophil", "mast")),
                 levels = c("eosinophil", "neutrophil", "mast"))
ng_tab <- data.frame(group = factor(grp[!is.na(grp)], levels = rev(grp_lv))) |> count(group)
p_dot <- ggplot(dd, aes(gene, group, size = pct, colour = mean_counts)) + geom_point() +
  scale_size_area(max_size = 10, name = "% cells\npositive") +
  scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\ncounts") +
  facet_grid(~ set, scales = "free_x", space = "free_x") + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = NULL, title = sprintf("Is the eosinophil call clean? (call: %s)", UC_EOS_CALL),
       subtitle = "called eosinophils should carry SIGLEC8 + RNASE3 and little of the neutrophil / mast genes")
p_n <- ggplot(ng_tab, aes(n, group)) + geom_col(fill = "grey60", width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3.5) +
  scale_x_continuous(expand = expansion(c(0, 0.25))) + theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = "cells", y = NULL, title = "n")
save_fig(p_dot + p_n + plot_layout(widths = c(4, 1)), "05uc_eosinophilia_marker_check", width = 16, height = 5.5)
lab_tab <- data.frame(label_04 = m$label[m$eos]) |> count(label_04, name = "n_called_eosinophils") |>
  mutate(pct = round(100 * n_called_eosinophils / sum(n_called_eosinophils), 1)) |> arrange(desc(n_called_eosinophils))
save_table(lab_tab, "05uc_eosinophilia_called_cells_04_labels")
message("04 labels of the called eosinophils (most sit in the neutrophil / granulocyte cluster if 04 merged them):")
print(utils::head(as.data.frame(lab_tab), 6), row.names = FALSE)

## -----------------------------------------------------------------------
## 1. eosinophils and eosinophilia per core
## -----------------------------------------------------------------------
ec <- m |> group_by(core_id) |>
  summarise(n_cells = n(), area_mm2 = core_area_mm2(x_um, y_um), eosinophils = sum(eos),
            cells_SIGLEC8_or_RNASE3 = sum(eos_any), immune_cells = sum(compartment == "Immune"), .groups = "drop") |>
  mutate(eos_per_mm2 = eosinophils / area_mm2,
         eos_pct_immune = ifelse(immune_cells >= UC_MIN_DENOM, 100 * eosinophils / immune_cells, NA),
         SIGLEC8_or_RNASE3_cells_per_mm2 = cells_SIGLEC8_or_RNASE3 / area_mm2) |>
  as.data.frame()
rownames(ec) <- ec$core_id
ec <- cbind(core_tab[ec$core_id, c("slide_id", "patient_id", "nancy", "grade")], ec)
q <- stats::quantile(ec$eos_per_mm2, c(0.25, 0.75, 0.9), na.rm = TRUE)
thr <- q[[2]] + 1.5 * (q[[2]] - q[[1]])
if (!any(ec$eos_per_mm2 > thr)) thr <- q[[3]]
ec$eosinophilia_tma <- factor(ifelse(ec$eos_per_mm2 > thr & ec$eosinophils >= 3, "yes", "no"), levels = c("no", "yes"))
path_col <- NULL
if (!is.null(clin)) for (c in intersect(UC_CLIN_EOS_COLS, colnames(clin))) {
  v <- clin[ec$core_id, c]
  if (sum(!is.na(v)) < 3) next
  ec[[paste0("pathology_", c)]] <- v
  o <- as_outcome(v)
  if (is.null(path_col) && nlevels(o) == 2) { path_col <- c; ec$eosinophilia_pathology <- factor(as.character(o), levels = c("no", "yes")) }
}
if (!is.null(clin) && is.null(path_col))
  message("No yes / no pathology eosinophil column among ", paste(UC_CLIN_EOS_COLS, collapse = ", "),
          " -- eosinophilia from the TMA only.")
by <- if (UC_EOSINOPHILIA_BY == "auto") (if (is.null(path_col)) "tma" else "pathology") else UC_EOSINOPHILIA_BY
if (by == "pathology" && is.null(path_col)) { message("UC_EOSINOPHILIA_BY = pathology, but no pathology call -> TMA"); by <- "tma" }
ec$eosinophilia <- if (by == "pathology") ec$eosinophilia_pathology else ec$eosinophilia_tma
ec <- ec[order(-ec$eos_per_mm2), ]
save_table(ec, "05uc_eosinophilia_per_core")
message(sprintf("Eosinophilia (%s): %d of %d cores. TMA threshold: %.0f eosinophils / mm2.",
                by, sum(ec$eosinophilia %in% "yes"), nrow(ec), thr))
if (!is.null(path_col)) {
  ct <- table(TMA = ec$eosinophilia_tma, pathology = ec$eosinophilia_pathology)
  message("TMA call vs pathology call (", path_col, "):"); print(ct)
}

## -----------------------------------------------------------------------
## 2. relation to Nancy
## -----------------------------------------------------------------------
st <- if (any(!is.na(ec$nancy))) assoc_nancy(ec, c("eos_per_mm2", "eos_pct_immune")) else data.frame()
save_table(st, "05uc_eosinophilia_vs_nancy")
gt <- ec[!is.na(ec$eosinophilia), ] |> count(grade, eosinophilia, .drop = FALSE) |> group_by(grade) |>
  mutate(n_cores = sum(n), pct = ifelse(n_cores > 0, 100 * n / n_cores, NA)) |> ungroup() |> filter(n_cores > 0)
save_table(gt, "05uc_eosinophilia_by_nancy_grade")
p_trend <- NA
if (any(!is.na(ec$nancy)) && length(unique(stats::na.omit(ec$eosinophilia))) == 2) {
  d <- ec[!is.na(ec$nancy) & !is.na(ec$eosinophilia), ]
  p_trend <- tryCatch(suppressWarnings(stats::cor.test(d$nancy, as.integer(d$eosinophilia == "yes"),
                                                       method = "spearman", exact = FALSE))$p.value, error = function(e) NA)
}
lab <- if (nrow(st)) st[st$feature == "eos_per_mm2", ] else NULL
p1 <- ggplot(ec, aes(grade, eos_per_mm2)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(aes(colour = eosinophilia), position = position_jitter(width = 0.15, seed = 1), size = 2) +
  geom_hline(yintercept = thr, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c(no = "grey55", yes = EOS_PINK), name = paste0("eosinophilia\n(", by, ")"), na.value = "grey80") +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "eosinophils per mm2\n(SIGLEC8 + RNASE3 co-expressing)",
       title = "A. Eosinophils by Nancy grade (one dot = one core)",
       subtitle = if (!is.null(lab) && nrow(lab)) sprintf("Spearman rho %.2f, p %.2g%s\ndashed line = TMA eosinophilia threshold",
                                                        lab$spearman_rho, lab$p_spearman,
                                                        ifelse(is.na(lab$p_mixed_model), "", sprintf("; with patient p %.2g", lab$p_mixed_model))) else NULL)
p2 <- ggplot(gt[gt$eosinophilia == "yes", ], aes(grade, pct, fill = grade)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d / %d", n, n_cores)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = grade_cols, guide = "none") + scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "% of cores with eosinophilia", title = sprintf("B. Cores with eosinophilia (%s call)", by),
       subtitle = if (!is.na(p_trend)) sprintf("trend with Nancy: p %.2g", p_trend) else NULL)
pl <- p1 | p2
if (!is.null(path_col)) {
  p3 <- ggplot(ec[!is.na(ec$eosinophilia_pathology), ], aes(eosinophilia_pathology, eos_per_mm2)) +
    geom_boxplot(outlier.shape = NA, colour = "grey50") +
    geom_point(aes(colour = grade), position = position_jitter(width = 0.15, seed = 1), size = 2) +
    scale_colour_manual(values = grade_cols, name = NULL) + theme_bw(base_size = 12) +
    labs(x = paste("pathology:", path_col), y = "eosinophils per mm2 (TMA)",
         title = "C. Does the TMA count match the pathologist?",
         subtitle = sprintf("AUC %.2f", auc_mw(ec$eos_per_mm2[ec$eosinophilia_pathology %in% "yes"],
                                                ec$eos_per_mm2[ec$eosinophilia_pathology %in% "no"])))
  pl <- p1 | p2 | p3
}
save_fig(pl, "05uc_eosinophilia_vs_nancy", width = if (is.null(path_col)) 14 else 20, height = 6)

## -----------------------------------------------------------------------
## 3. core images on black
## -----------------------------------------------------------------------
plot_eos_core <- function(cid, panel_in = 4.2) {
  d <- m[m$core_id == cid, ]
  ext <- max(diff(range(d$x_um)), diff(range(d$y_um)), 50)
  d$size <- cell_dot_sizes(d$cell_area %||% NA, ext, panel_in)
  e <- d[d$eos, ]
  r <- ec[cid, ]
  ggplot(mapping = aes(x_um, y_um)) +
    geom_point(data = d, aes(size = size, colour = compartment), stroke = 0) +
    geom_point(data = e, aes(size = size * 1.9), colour = EOS_PINK, stroke = 0) +
    scale_colour_manual(values = COMP_BLACK, name = NULL, drop = FALSE) + scale_size_identity() +
    scale_y_reverse() + coord_fixed() + theme_black(10) + theme(legend.position = "none") +
    annotate("segment", x = min(d$x_um), xend = min(d$x_um) + 100, y = max(d$y_um) + 0.03 * ext,
             yend = max(d$y_um) + 0.03 * ext, colour = "white", linewidth = 0.7) +
    labs(title = sprintf("%s  -  Nancy %s", cid, ifelse(is.na(r$nancy), "NA", r$nancy)),
         subtitle = sprintf("%d eosinophils, %.0f / mm2%s", r$eosinophils, r$eos_per_mm2,
                            if (!is.null(path_col)) paste0("  |  pathology: ", ifelse(is.na(r$eosinophilia_pathology), "NA",
                                                                                       as.character(r$eosinophilia_pathology))) else ""))
}
legend_el <- legend_black(
  ggplot(data.frame(k = factor(c("Epithelium", "Immune", "Stroma", "eosinophil"),
                               levels = c("Epithelium", "Immune", "Stroma", "eosinophil")), x = 1:4, y = 1),
         aes(x, y, colour = k)) + geom_point(size = 4) +
    scale_colour_manual(values = c(COMP_BLACK[c("Epithelium", "Immune", "Stroma")], eosinophil = EOS_PINK), name = NULL) +
    theme_black() + guides(colour = guide_legend(ncol = 1)))
eos_txt <- if (UC_EOS_CALL == "lenient") "pink = cells with SIGLEC8 or RNASE3" else
  "pink = eosinophils (SIGLEC8 + RNASE3 co-expressing, no neutrophil profile)"
grid_fig <- function(cores, title, file) {
  if (!length(cores)) return(invisible(NULL))
  ps <- lapply(cores, plot_eos_core)
  nc <- min(4, length(ps))
  fig <- patchwork::wrap_plots(ps, ncol = nc)
  if (!is.null(legend_el)) fig <- patchwork::wrap_plots(fig, legend_el, widths = c(4.3 * nc, 2.2))
  fig <- fig + patchwork::plot_annotation(
    title = title, subtitle = paste0(eos_txt, "; blue = epithelium, red = immune, green = stroma; bar = 100 um"),
    theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                  plot.title = element_text(colour = "white", face = "bold", size = 17),
                  plot.subtitle = element_text(colour = "grey80", size = 11)))
  save_fig(fig, file, width = 4.3 * nc + 2.6, height = 4.6 * ceiling(length(ps) / nc) + 1.2, dpi = 200)
}
ok <- ec$n_cells >= 200
eo <- ec[ok & ec$eosinophilia %in% "yes", ]
eo <- utils::head(eo[order(-eo$eos_per_mm2), ], UC_EOS_MAX_CORES)
eo <- eo[order(eo$nancy, -eo$eos_per_mm2), ]
if (!nrow(eo)) message("No core with eosinophilia (", by, ") -- no core images.")
grid_fig(eo$core_id, sprintf("Cores with eosinophilia (%s call), ordered by Nancy", by), "05uc_eosinophilia_cores")
## controls: per eosinophilia core, the low-eosinophil core of the same grade with the fewest eosinophils
ctrl <- c()
for (i in seq_len(nrow(eo))) {
  cand <- ec[ok & ec$eosinophilia %in% "no" & ec$grade == eo$grade[i] & !ec$core_id %in% ctrl, ]
  if (nrow(cand)) ctrl <- c(ctrl, cand$core_id[which.min(cand$eos_per_mm2)])
}
grid_fig(ctrl, "For comparison: cores WITHOUT eosinophilia, same Nancy grades", "05uc_eosinophilia_controls")
message("Done.")
