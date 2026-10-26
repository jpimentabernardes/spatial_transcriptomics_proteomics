## =============================================================================
## 05_uc_cellscape_concordance.R
##
## UC biopsies: do the Xenium immune cell densities (05_uc_immune_density_table.R)
## show the same patterns as the CellScape (protein) densities of the same cores?
##
## Run:   Rscript 05_uc_cellscape_concordance.R /path/to/cellscape_densities.csv
##        EXP1_SEGMENTATION=xenium Rscript 05_uc_cellscape_concordance.R <csv>   (Xenium-segmentation call)
## (or set UC_CELLSCAPE_DENSITY_CSV in config.R / EXP1_CELLSCAPE_DENSITY)
##
## Cores are matched on TMA_Sector_join + Koordinate_join (A + 1c; or tma_core A_1-c), so
## slide names may differ. Populations = the *_density_per_mm2 columns present in both
## tables, plus sums that do not depend on how the subsets are split:
##   Macrophages_all = M1 + M2 + Special, T_cells_CD4_all = CD4 + CD4FOXP3,
##   T_cells_CD8_all = CD8 + CD8FOXP3, T_cells_all = CD4_all + CD8_all + GammaDelta
## Cores with fewer than UC_CONC_MIN_CELLS cells or UC_CONC_MIN_AREA mm2 in Xenium
## (fragments) are left out.
##
## Per population, over the matched cores:
##   spearman_density     same ranking of cores? (insensitive to a constant offset, e.g. a
##                        different area definition)                -> "similar pattern"
##   spearman_proportion  same ranking of the population's SHARE of the immune cells
##                        (removes area / overall density differences altogether)
##   ccc_log              Lin's concordance on log10(1 + density): same absolute numbers?
##   median_ratio         Xenium / CellScape (cores where both > 0): systematic offset
## Per core: Spearman across populations of the immune composition (does the core's
## immune profile agree?).
##
## OUTPUT: tables/05uc_cellscape_concordance_<seg>_*.csv,
##         figures/05uc_cellscape_concordance_<seg>_scatter.png   one panel per population
##         figures/05uc_cellscape_concordance_<seg>_summary.png   metrics + mean composition
##         figures/05uc_cellscape_concordance_<seg>_per_core.png  agreement of each core's profile
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
CS_CSV <- if (length(args)) args[1] else
  Sys.getenv("EXP1_CELLSCAPE_DENSITY", get0("UC_CELLSCAPE_DENSITY_CSV", ifnotfound = ""))
UC_CONC_MIN_CELLS <- get0("UC_CONC_MIN_CELLS", ifnotfound = 500)
UC_CONC_MIN_AREA  <- get0("UC_CONC_MIN_AREA", ifnotfound = 0.1)
if (!nzchar(CS_CSV) || !file.exists(CS_CSV)) stop("CellScape density table not found: '", CS_CSV,
                                                  "'. Usage: Rscript 05_uc_cellscape_concordance.R /path/to/cellscape.csv")
XEN_CSV <- file.path(PROJECT_DIR, "data", "export", sprintf("05uc_immune_densities_per_core_%s.csv", SEGMENTATION))
if (!file.exists(XEN_CSV)) stop("Xenium density table not found: ", XEN_CSV, " -- run 05_uc_immune_density_table.R with ",
                                "EXP1_SEGMENTATION=", SEGMENTATION, " first.")
pre <- paste0("05uc_cellscape_concordance_", SEGMENTATION)
message("Segmentation: ", SEGMENTATION, "\n  Xenium:    ", XEN_CSV, "\n  CellScape: ", CS_CSV)

## -----------------------------------------------------------------------
## 0. read both tables (CellScape: any separator, decimal comma tolerated)
## -----------------------------------------------------------------------
read_any <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  n_sep <- c("\t" = lengths(regmatches(first, gregexpr("\t", first))), ";" = lengths(regmatches(first, gregexpr(";", first))),
             "," = lengths(regmatches(first, gregexpr(",", first))))
  sep <- names(n_sep)[which.max(n_sep)]
  d <- read_config_csv(path, sep = sep, colClasses = "character", check.names = FALSE, fileEncoding = "latin1")
  d[] <- lapply(d, trimws)
  d
}
to_num <- function(v) { v <- gsub("\\s", "", v); if (!any(grepl("\\.", v)) && any(grepl(",", v))) v <- sub(",", ".", v, fixed = TRUE)
  suppressWarnings(as.numeric(v)) }
add_key <- function(d, what) {
  if (all(c("TMA_Sector_join", "Koordinate_join") %in% colnames(d))) {
    d$key <- core_key(d$TMA_Sector_join, d$Koordinate_join)
  } else if ("tma_core" %in% colnames(d)) {
    d$key <- core_key_from_id(d$tma_core)
  } else stop(what, ": needs TMA_Sector_join + Koordinate_join or tma_core columns; found ", paste(colnames(d), collapse = ", "))
  d
}
xen <- add_key(read_any(XEN_CSV), "Xenium table")
cs  <- add_key(read_any(CS_CSV), "CellScape table")
dcols <- intersect(grep("_density_per_mm2$", colnames(xen), value = TRUE), grep("_density_per_mm2$", colnames(cs), value = TRUE))
if (!length(dcols)) stop("No *_density_per_mm2 columns shared by the two tables.")
only_cs <- setdiff(grep("_density_per_mm2$", colnames(cs), value = TRUE), dcols)
if (length(only_cs)) message("CellScape columns without a Xenium counterpart (left out): ", paste(only_cs, collapse = ", "))
for (c in dcols) { xen[[c]] <- to_num(xen[[c]]); cs[[c]] <- to_num(cs[[c]]) }
if (any(duplicated(cs$key))) {
  message("CellScape: ", sum(duplicated(cs$key)), " duplicated cores -> averaged")
  lab <- cs[!duplicated(cs$key), setdiff(colnames(cs), dcols)]
  cs <- cbind(lab, stats::aggregate(cs[dcols], list(key = cs$key), mean, na.rm = TRUE)[match(lab$key, unique(sort(cs$key))), dcols])
  cs <- cs[, unique(colnames(cs))]
}

## Xenium core size (fragments out): counts table of the same segmentation
cnt_f <- file.path(TABLES_DIR, paste0("05uc_immune_density_counts_per_core_", SEGMENTATION, ".csv"))
if (!file.exists(cnt_f)) {                                                   # older runs: no suffix -- check it matches
  old <- file.path(TABLES_DIR, "05uc_immune_density_counts_per_core.csv")
  if (file.exists(old)) {
    o <- utils::read.csv(old, check.names = FALSE)
    k <- core_key_from_id(o$tma_core); j <- match(k, xen$key)
    same <- !is.na(j) & "B_cells" %in% colnames(o) && isTRUE(all.equal(o$B_cells[!is.na(j)] / o$area_mm2[!is.na(j)],
                                                                       xen$B_cells_density_per_mm2[j[!is.na(j)]], tolerance = 1e-3))
    if (isTRUE(same)) cnt_f <- old else message("tables/05uc_immune_density_counts_per_core.csv is from another run -- ",
                                                 "re-run 05_uc_immune_density_table.R to get the core sizes.")
  }
}
xen$n_cells <- NA_real_; xen$area_mm2 <- NA_real_
if (file.exists(cnt_f)) {
  q <- utils::read.csv(cnt_f, check.names = FALSE)
  j <- match(xen$key, core_key_from_id(q$tma_core))
  xen$n_cells <- q$n_cells[j]; xen$area_mm2 <- q$area_mm2[j]
} else message("No core sizes found -- fragments are NOT filtered.")

## -----------------------------------------------------------------------
## 1. match cores
## -----------------------------------------------------------------------
both <- intersect(xen$key, cs$key)
message(sprintf("Cores: Xenium %d, CellScape %d, matched %d", nrow(xen), nrow(cs), length(both)))
if (length(setdiff(cs$key, xen$key))) message("  CellScape cores not in Xenium (first 10): ",
                                              paste(utils::head(cs$tma_core[!cs$key %in% xen$key] %||% setdiff(cs$key, xen$key), 10), collapse = ", "))
if (length(setdiff(xen$key, cs$key))) message("  Xenium cores not in CellScape (first 10): ",
                                              paste(utils::head(xen$tma_core[!xen$key %in% cs$key], 10), collapse = ", "))
small <- xen$key[!is.na(xen$n_cells) & (xen$n_cells < UC_CONC_MIN_CELLS | xen$area_mm2 < UC_CONC_MIN_AREA)]
if (length(intersect(small, both))) message(sprintf("  left out as fragments (< %d cells or < %.2f mm2): %s", UC_CONC_MIN_CELLS,
                                                    UC_CONC_MIN_AREA, paste(xen$tma_core[xen$key %in% intersect(small, both)], collapse = ", ")))
both <- setdiff(both, small)
if (length(both) < 5) stop("Fewer than 5 matched cores -- check the core names in both tables.")
X <- xen[match(both, xen$key), ]; C <- cs[match(both, cs$key), ]
pop <- sub("_density_per_mm2$", "", dcols)
XM <- as.matrix(X[dcols]); CM <- as.matrix(C[dcols]); colnames(XM) <- colnames(CM) <- pop
rownames(XM) <- rownames(CM) <- X$tma_core
add_sum <- function(M, name, parts) {
  parts <- intersect(parts, colnames(M))
  if (length(parts) < 2) return(M)
  v <- rowSums(M[, parts, drop = FALSE], na.rm = TRUE); v[rowSums(!is.na(M[, parts, drop = FALSE])) == 0] <- NA
  cbind(M, stats::setNames(data.frame(v), name)) |> as.matrix()
}
for (s in list(list("Macrophages_all", c("M1_Macrophages", "M2_Macrophages", "Special_Macrophages")),
               list("T_cells_CD4_all", c("T_cells_CD4", "T_cells_CD4FOXP3")),
               list("T_cells_CD8_all", c("T_cells_CD8", "T_cells_CD8FOXP3")))) {
  XM <- add_sum(XM, s[[1]], s[[2]]); CM <- add_sum(CM, s[[1]], s[[2]])
}
XM <- add_sum(XM, "T_cells_all", c("T_cells_CD4_all", "T_cells_CD8_all", "T_cells_GammaDelta"))
CM <- add_sum(CM, "T_cells_all", c("T_cells_CD4_all", "T_cells_CD8_all", "T_cells_GammaDelta"))
all_pop <- colnames(XM)
## proportions of the immune cells (base populations defined in both platforms)
base <- pop[colSums(!is.na(XM[, pop, drop = FALSE])) > 0 & colSums(!is.na(CM[, pop, drop = FALSE])) > 0]
prop <- function(M) { tot <- rowSums(M[, base, drop = FALSE], na.rm = TRUE); sweep(M, 1, ifelse(tot > 0, tot, NA), "/") }
XP <- prop(XM); CP <- prop(CM)
long <- data.frame(tma_core = rep(rownames(XM), length(all_pop)), population = rep(all_pop, each = nrow(XM)),
                   xenium = as.vector(XM), cellscape = as.vector(CM),
                   xenium_share = as.vector(XP), cellscape_share = as.vector(CP),
                   sector = rep(toupper(sub("_.*$", "", rownames(XM))), length(all_pop)))
save_table(long, paste0(pre, "_matched_cores"))

## -----------------------------------------------------------------------
## 2. concordance per population
## -----------------------------------------------------------------------
lin_ccc <- function(a, b) { if (length(a) < 3 || stats::var(a) + stats::var(b) == 0) return(NA)
  2 * stats::cov(a, b) / (stats::var(a) + stats::var(b) + (mean(a) - mean(b))^2) }
sp <- function(a, b) { ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 5 || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(c(NA, NA))
  t <- suppressWarnings(stats::cor.test(a[ok], b[ok], method = "spearman", exact = FALSE)); c(unname(t$estimate), t$p.value) }
met <- dplyr::bind_rows(lapply(all_pop, function(p) {
  x <- XM[, p]; y <- CM[, p]; ok <- !is.na(x) & !is.na(y)
  if (!any(ok)) return(data.frame(population = p, n_cores = 0, note = "not defined in one of the tables"))
  a <- sp(x, y); b <- sp(XP[, p], CP[, p]); both_pos <- ok & x > 0 & y > 0
  data.frame(population = p, n_cores = sum(ok),
             median_xenium = stats::median(x[ok]), median_cellscape = stats::median(y[ok]),
             pct_cores_zero_xenium = 100 * mean(x[ok] == 0), pct_cores_zero_cellscape = 100 * mean(y[ok] == 0),
             spearman_density = a[1], p_density = a[2], spearman_proportion = b[1], p_proportion = b[2],
             ccc_log = lin_ccc(log10(1 + x[ok]), log10(1 + y[ok])),
             median_ratio_xenium_over_cellscape = if (sum(both_pos) >= 3) stats::median(x[both_pos] / y[both_pos]) else NA,
             note = "")
}))
met$padj_density <- stats::p.adjust(met$p_density, "BH")
met$verdict <- dplyr::case_when(is.na(met$spearman_density) ~ "not comparable (missing / constant)",
                                met$spearman_density >= 0.6 & met$padj_density < 0.05 ~ "similar pattern",
                                met$spearman_density >= 0.3 & met$padj_density < 0.05 ~ "partly similar",
                                TRUE ~ "no clear agreement")
save_table(met, paste0(pre, "_metrics"))
print(met[, c("population", "n_cores", "spearman_density", "spearman_proportion", "ccc_log",
              "median_ratio_xenium_over_cellscape", "verdict")], row.names = FALSE, digits = 2)

## per core: Spearman across base populations of the immune composition
pc <- dplyr::bind_rows(lapply(rownames(XP), function(k) {
  a <- XP[k, base]; b <- CP[k, base]; ok <- !is.na(a) & !is.na(b)
  r <- if (sum(ok) >= 5 && stats::sd(a[ok]) > 0 && stats::sd(b[ok]) > 0) stats::cor(a[ok], b[ok], method = "spearman") else NA
  data.frame(tma_core = k, n_populations = sum(ok), spearman_composition = r,
             immune_xenium_per_mm2 = sum(XM[k, base], na.rm = TRUE), immune_cellscape_per_mm2 = sum(CM[k, base], na.rm = TRUE))
}))
pc <- pc[order(pc$spearman_composition), ]
save_table(pc, paste0(pre, "_per_core"))

## -----------------------------------------------------------------------
## 3. figures
## -----------------------------------------------------------------------
lv <- intersect(c(base, setdiff(all_pop, pop)), met$population[met$n_cores > 0])
lv_sc <- lv[!is.na(met$spearman_density[match(lv, met$population)])]          # panels: comparable populations only
ld <- long[long$population %in% lv_sc & !is.na(long$xenium) & !is.na(long$cellscape), ]
ld$population <- factor(ld$population, levels = lv_sc)
lab <- met[match(lv_sc, met$population), ]
strip <- stats::setNames(sprintf("%s\nrho %.2f | CCC %.2f | Xen/CS x%.2g", lv_sc, lab$spearman_density, lab$ccc_log,
                                 lab$median_ratio_xenium_over_cellscape), lv_sc)
## identical x and y range per panel (so the dashed identity line is the 45-degree diagonal)
lims <- dplyr::bind_rows(lapply(lv_sc, function(p) { v <- c(ld$xenium[ld$population == p], ld$cellscape[ld$population == p])
  data.frame(population = factor(p, levels = lv_sc), xenium = range(v), cellscape = range(v)) }))
sec_cols <- c(A = "#3B6FB6", C = "#C0392B", D = "#2E8B57")
br <- c(0, 1, 10, 100, 1000, 10000)
p_sc <- ggplot(ld, aes(cellscape, xenium)) +
  geom_blank(data = lims) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_point(aes(colour = sector), size = 1.6, alpha = 0.8) +
  scale_x_continuous(trans = "log1p", breaks = br) + scale_y_continuous(trans = "log1p", breaks = br) +
  scale_colour_manual(values = sec_cols, name = "sector", na.value = "grey50") +
  facet_wrap(~ population, scales = "free", ncol = 4, labeller = labeller(population = strip)) +
  theme_bw(base_size = 11) + theme(aspect.ratio = 1, strip.text = element_text(size = 9.5)) +
  labs(x = "CellScape density (cells / mm2)", y = sprintf("Xenium density (%s)", SEGMENTATION),
       title = sprintf("Xenium (%s) vs CellScape, per core (one dot = one core, n = %d)", SEGMENTATION, length(both)),
       subtitle = paste("rho = Spearman (same ranking of cores = similar pattern); CCC = agreement of absolute values (log);",
                        "\nXen/CS = median ratio; dashed = identity (points along a line parallel to it = same pattern, constant offset)"))
save_fig(p_sc, paste0(pre, "_scatter"), width = 16, height = 3.8 * ceiling(length(lv_sc) / 4) + 1.5)

mt <- tidyr::pivot_longer(met[met$population %in% lv, c("population", "spearman_density", "spearman_proportion", "ccc_log")],
                          -population, names_to = "metric", values_to = "value")
mt$metric <- factor(mt$metric, levels = c("spearman_density", "spearman_proportion", "ccc_log"),
                    labels = c("pattern\n(Spearman, density)", "composition\n(Spearman, share)", "absolute\n(CCC, log)"))
mt$population <- factor(mt$population, levels = rev(lv))
p_m <- ggplot(mt, aes(metric, population, fill = value)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(value), "NA", sprintf("%.2f", value))), size = 3.6) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1), na.value = "grey90",
                       name = NULL) +
  scale_x_discrete(position = "top") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "Agreement per population", subtitle = "blue = concordant (1 = perfect)")
mc <- rbind(data.frame(platform = paste0("Xenium (", SEGMENTATION, ")"), population = base, share = colMeans(XP[, base, drop = FALSE], na.rm = TRUE)),
            data.frame(platform = "CellScape", population = base, share = colMeans(CP[, base, drop = FALSE], na.rm = TRUE)))
pal_pop <- stats::setNames(grDevices::hcl.colors(length(base), "Dark 3"), base)
p_c <- ggplot(mc, aes(platform, 100 * share, fill = population)) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(share >= 0.04, sprintf("%.0f%%", 100 * share), "")), position = position_stack(vjust = 0.5), size = 3.2) +
  scale_fill_manual(values = pal_pop, name = NULL) + scale_y_continuous(expand = c(0, 0)) + theme_bw(base_size = 12) +
  labs(x = NULL, y = "% of the immune cells (mean over matched cores)", title = "Mean immune composition")
save_fig(p_m + p_c + plot_layout(widths = c(1.1, 1)), paste0(pre, "_summary"), width = 15, height = max(7, 0.4 * length(lv) + 2.5))

pc$tma_core <- factor(pc$tma_core, levels = pc$tma_core)
p_pc <- ggplot(pc[!is.na(pc$spearman_composition), ], aes(tma_core, spearman_composition)) +
  geom_hline(yintercept = c(0, stats::median(pc$spearman_composition, na.rm = TRUE)), linetype = c(1, 2), colour = c("grey60", "#2166AC")) +
  geom_col(aes(fill = spearman_composition), width = 0.8) +
  scale_fill_gradient2(low = "#B2182B", mid = "grey85", high = "#2166AC", midpoint = 0, limits = c(-1, 1), guide = "none") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
  labs(x = "core (sorted)", y = "Spearman across populations",
       title = "Does each core's immune composition agree between Xenium and CellScape?",
       subtitle = sprintf("dashed = median (%.2f); low bars = cores to look at (tables/%s_per_core.csv)",
                          stats::median(pc$spearman_composition, na.rm = TRUE), pre))
save_fig(p_pc, paste0(pre, "_per_core"), width = max(10, 0.09 * nrow(pc) + 3), height = 5)
message("Done.")
