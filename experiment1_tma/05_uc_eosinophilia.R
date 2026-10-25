## =============================================================================
## 05_uc_eosinophilia.R
##
## UC biopsies: which cores show eosinophilia, what they look like, and how it
## relates to Nancy -- WITHOUT relying on an eosinophil cluster (04 / 04d).
##
## Eosinophil signal per core, from the transcripts of the eosinophil-specific genes
## (SIGLEC8 + RNASE3 on this panel), counted only in non-epithelial cells that are
## not mast cells (mast cells also express SIGLEC8; ALOX15 is also epithelial):
##   eos_transcripts_per_mm2   SIGLEC8 + RNASE3 transcripts per mm^2 of tissue
##   eos_cells_per_mm2         cells with >= 1 of them, per mm^2
## Eosinophilia per core:
##   pathology   the clinical eosinophil column (UC_CLIN_EOS_COLS), if it is a yes / no call
##   TMA         eos_transcripts_per_mm2 above Q3 + 1.5 IQR of all cores (or the top 10%
##               if no core is that far out)
## UC_EOSINOPHILIA_BY = "pathology" / "tma" / "auto" (pathology when available) picks the
## call used for the images.
##
## Figures
##   05uc_eosinophilia_vs_nancy.png   eosinophil signal by Nancy grade (cores; Spearman,
##                                    mixed model with patient), % of cores with eosinophilia
##                                    per Nancy grade, TMA signal vs the pathology call
##   05uc_eosinophilia_cores.png      the eosinophilia cores on black, ordered by Nancy:
##                                    tissue (epithelium / immune / stroma) and the
##                                    eosinophil-marker cells in pink (bigger = more transcripts)
##   05uc_eosinophilia_controls.png   for comparison: low-eosinophil cores of the same Nancy grades
##
## INPUT : xenium_<seg>_annotated.rds (04), tma_map Nancy / Nancy_TMA, clinical table (optional)
## OUTPUT: tables/05uc_eosinophilia_*.csv, figures/05uc_eosinophilia_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_EOS_SPECIFIC     <- get0("UC_EOS_SPECIFIC", ifnotfound = c("SIGLEC8", "RNASE3", "RNASE2", "CLC", "CCR3", "IL5RA", "PRG2", "EPX"))
UC_EOSINOPHILIA_BY  <- get0("UC_EOSINOPHILIA_BY", ifnotfound = "auto")
UC_EOS_MAX_CORES    <- get0("UC_EOS_MAX_CORES", ifnotfound = 12)      # cores shown per figure
COMP_BLACK <- c(Epithelium = "#2F5F8F", Immune = "#8C3B3B", Stroma = "#3C7A4A", Other = "#4A4A4A")
EOS_PINK <- "#FF3EA5"

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
rm(xen); invisible(gc())
eg <- match_case(UC_EOS_SPECIFIC, rownames(counts)); eg <- eg[eg %in% rownames(counts)]
if (!length(eg)) stop("None of ", paste(UC_EOS_SPECIFIC, collapse = ", "), " is on the panel.")
message("Eosinophil-specific genes used: ", paste(eg, collapse = ", "))
eos_n <- Matrix::colSums(counts[eg, , drop = FALSE])
m$eligible <- m$compartment %in% c("Immune", "Stroma", "Other") & !m$class %in% "Mast"
m$eos_tx <- ifelse(m$eligible, eos_n, 0)
m$eos_pos <- m$eos_tx > 0

core_tab <- uc_core_table(m)
core_tab$grade <- factor(ifelse(is.na(core_tab$nancy), "NA", paste("Nancy", core_tab$nancy)),
                         levels = c(paste("Nancy", 0:4), "NA"))
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
clin <- join_uc_clinical(core_tab, read_uc_clinical())

## -----------------------------------------------------------------------
## 1. eosinophil signal and eosinophilia per core
## -----------------------------------------------------------------------
ec <- m |> group_by(core_id) |>
  summarise(n_cells = n(), area_mm2 = core_area_mm2(x_um, y_um), eos_transcripts = sum(eos_tx),
            eos_cells = sum(eos_pos), immune_cells = sum(compartment == "Immune"), .groups = "drop") |>
  mutate(eos_transcripts_per_mm2 = eos_transcripts / area_mm2, eos_cells_per_mm2 = eos_cells / area_mm2,
         eos_cells_pct_immune = ifelse(immune_cells >= UC_MIN_DENOM, 100 * eos_cells / immune_cells, NA)) |>
  as.data.frame()
rownames(ec) <- ec$core_id
ec <- cbind(core_tab[ec$core_id, c("slide_id", "patient_id", "nancy", "grade")], ec)
q <- stats::quantile(ec$eos_transcripts_per_mm2, c(0.25, 0.75, 0.9), na.rm = TRUE)
thr <- q[[2]] + 1.5 * (q[[2]] - q[[1]])
if (!any(ec$eos_transcripts_per_mm2 > thr)) thr <- q[[3]]
ec$eosinophilia_tma <- factor(ifelse(ec$eos_transcripts_per_mm2 > thr, "yes", "no"), levels = c("no", "yes"))
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
          " -- eosinophilia from the TMA signal only.")
by <- if (UC_EOSINOPHILIA_BY == "auto") (if (is.null(path_col)) "tma" else "pathology") else UC_EOSINOPHILIA_BY
if (by == "pathology" && is.null(path_col)) { message("UC_EOSINOPHILIA_BY = pathology, but no pathology call -> TMA"); by <- "tma" }
ec$eosinophilia <- if (by == "pathology") ec$eosinophilia_pathology else ec$eosinophilia_tma
ec <- ec[order(-ec$eos_transcripts_per_mm2), ]
save_table(ec, "05uc_eosinophilia_per_core")
message(sprintf("Eosinophilia (%s): %d of %d cores. TMA threshold: %.0f eosinophil-gene transcripts / mm2.",
                by, sum(ec$eosinophilia %in% "yes"), nrow(ec), thr))
if (!is.null(path_col)) {
  ct <- table(TMA = ec$eosinophilia_tma, pathology = ec$eosinophilia_pathology)
  message("TMA call vs pathology call (", path_col, "):"); print(ct)
}

## -----------------------------------------------------------------------
## 2. relation to Nancy
## -----------------------------------------------------------------------
st <- if (any(!is.na(ec$nancy))) assoc_nancy(ec, c("eos_transcripts_per_mm2", "eos_cells_per_mm2", "eos_cells_pct_immune")) else data.frame()
save_table(st, "05uc_eosinophilia_vs_nancy")
gt <- ec[!is.na(ec$eosinophilia), ] |> count(grade, eosinophilia) |> group_by(grade) |>
  mutate(n_cores = sum(n), pct = 100 * n / n_cores) |> ungroup()
save_table(gt, "05uc_eosinophilia_by_nancy_grade")
p_trend <- NA
if (any(!is.na(ec$nancy)) && length(unique(stats::na.omit(ec$eosinophilia))) == 2) {
  d <- ec[!is.na(ec$nancy) & !is.na(ec$eosinophilia), ]
  p_trend <- tryCatch(suppressWarnings(stats::cor.test(d$nancy, as.integer(d$eosinophilia == "yes"),
                                                       method = "spearman", exact = FALSE))$p.value, error = function(e) NA)
}
lab <- if (nrow(st)) st[st$feature == "eos_transcripts_per_mm2", ] else NULL
p1 <- ggplot(ec, aes(grade, eos_transcripts_per_mm2)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(aes(colour = eosinophilia), position = position_jitter(width = 0.15, seed = 1), size = 2) +
  geom_hline(yintercept = thr, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c(no = "grey55", yes = EOS_PINK), name = paste0("eosinophilia\n(", by, ")"), na.value = "grey80") +
  theme_bw(base_size = 12) + labs(x = NULL, y = "SIGLEC8 + RNASE3 transcripts per mm2\n(non-epithelial, non-mast cells)",
                                  title = "A. Eosinophil signal by Nancy grade (one dot = one core)",
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
  p3 <- ggplot(ec[!is.na(ec$eosinophilia_pathology), ], aes(eosinophilia_pathology, eos_transcripts_per_mm2)) +
    geom_boxplot(outlier.shape = NA, colour = "grey50") +
    geom_point(aes(colour = grade), position = position_jitter(width = 0.15, seed = 1), size = 2) +
    scale_colour_manual(values = grade_cols, name = NULL) + theme_bw(base_size = 12) +
    labs(x = paste("pathology:", path_col), y = "SIGLEC8 + RNASE3 transcripts per mm2",
         title = "C. Does the TMA signal match the pathologist?",
         subtitle = sprintf("AUC %.2f", auc_mw(ec$eos_transcripts_per_mm2[ec$eosinophilia_pathology %in% "yes"],
                                                ec$eos_transcripts_per_mm2[ec$eosinophilia_pathology %in% "no"])))
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
  e <- d[d$eos_pos, ]
  e <- e[order(e$eos_tx), ]
  r <- ec[cid, ]
  ggplot(mapping = aes(x_um, y_um)) +
    geom_point(data = d, aes(size = size, colour = compartment), stroke = 0) +
    geom_point(data = e, aes(size = size * (1.2 + 0.35 * pmin(eos_tx, 4))), colour = EOS_PINK, stroke = 0) +
    scale_colour_manual(values = COMP_BLACK, name = NULL, drop = FALSE) + scale_size_identity() +
    scale_y_reverse() + coord_fixed() + theme_black(10) + theme(legend.position = "none") +
    annotate("segment", x = min(d$x_um), xend = min(d$x_um) + 100, y = max(d$y_um) + 0.03 * ext,
             yend = max(d$y_um) + 0.03 * ext, colour = "white", linewidth = 0.7) +
    labs(title = sprintf("%s  -  Nancy %s", cid, ifelse(is.na(r$nancy), "NA", r$nancy)),
         subtitle = sprintf("%.0f eos transcripts / mm2%s", r$eos_transcripts_per_mm2,
                            if (!is.null(path_col)) paste0("  |  pathology: ", ifelse(is.na(r$eosinophilia_pathology), "NA",
                                                                                       as.character(r$eosinophilia_pathology))) else ""))
}
legend_el <- legend_black(
  ggplot(data.frame(k = factor(c("Epithelium", "Immune", "Stroma", "eosinophil-marker cell"),
                               levels = c("Epithelium", "Immune", "Stroma", "eosinophil-marker cell")), x = 1:4, y = 1),
         aes(x, y, colour = k)) + geom_point(size = 4) +
    scale_colour_manual(values = c(COMP_BLACK[c("Epithelium", "Immune", "Stroma")], "eosinophil-marker cell" = EOS_PINK), name = NULL) +
    theme_black() + guides(colour = guide_legend(ncol = 1)))
grid_fig <- function(cores, title, subtitle, file) {
  if (!length(cores)) return(invisible(NULL))
  ps <- lapply(cores, plot_eos_core)
  nc <- min(4, length(ps))
  fig <- patchwork::wrap_plots(ps, ncol = nc)
  if (!is.null(legend_el)) fig <- patchwork::wrap_plots(fig, legend_el, widths = c(4.3 * nc, 2.2))
  fig <- fig + patchwork::plot_annotation(title = title, subtitle = subtitle,
    theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                  plot.title = element_text(colour = "white", face = "bold", size = 17),
                  plot.subtitle = element_text(colour = "grey80", size = 11)))
  save_fig(fig, file, width = 4.3 * nc + 2.6, height = 4.6 * ceiling(length(ps) / nc) + 1.2, dpi = 200)
}
ok <- ec$n_cells >= 200
eo <- ec[ok & ec$eosinophilia %in% "yes", ]
eo <- utils::head(eo[order(-eo$eos_transcripts_per_mm2), ], UC_EOS_MAX_CORES)
eo <- eo[order(eo$nancy, -eo$eos_transcripts_per_mm2), ]
if (!nrow(eo)) message("No core with eosinophilia (", by, ") -- no core images.")
grid_fig(eo$core_id, sprintf("Cores with eosinophilia (%s call), ordered by Nancy", by),
         "pink = cells carrying SIGLEC8 / RNASE3 transcripts (bigger = more); blue = epithelium, red = immune, green = stroma; bar = 100 um",
         "05uc_eosinophilia_cores")
## controls: per eosinophilia core, the low-eosinophil core of the same grade with the lowest signal
ctrl <- c()
for (i in seq_len(nrow(eo))) {
  cand <- ec[ok & ec$eosinophilia %in% "no" & ec$grade == eo$grade[i] & !ec$core_id %in% ctrl, ]
  if (nrow(cand)) ctrl <- c(ctrl, cand$core_id[which.min(cand$eos_transcripts_per_mm2)])
}
grid_fig(ctrl, "For comparison: cores WITHOUT eosinophilia, same Nancy grades",
         "pink = cells carrying SIGLEC8 / RNASE3 transcripts; blue = epithelium, red = immune, green = stroma; bar = 100 um",
         "05uc_eosinophilia_controls")
message("Done.")
