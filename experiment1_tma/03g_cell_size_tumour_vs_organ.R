## =============================================================================
## 03g_cell_size_tumour_vs_organ.R
##
## Are cells in TMA_TUMOR bigger than in TMA_ORGAN, and do they express more?
## Asked for each segmentation (10x, cellpose, segger), because cell size IS a
## segmentation result: nucleus expansion gives every cell nearly the same
## ring, so real size differences show up mainly in the nucleus area and in
## methods that follow the cell outline.
##
## Per cell: cell area, nucleus area, transcripts (counts) per cell, genes per
## cell, transcripts per um2 of cell (density = expression corrected for size).
## Per cell also a quick lineage from the L1 markers of config/marker_panel_rna.csv
## (the marker set with the most counts per gene, if it clearly wins), the SAME
## rule for every method, so lineages are comparable between methods.
##
## Statistical unit = core (median per core), never the cell.
##   1. all tumour cores vs all normal cores, per method x lineage
##   2. organ-matched: tumour vs normal of the same site (Colon, Stomach, ...;
##      `site` column of config/annotation_families.csv)
## Read immune cells as the internal control: tumour epithelium may be bigger
## for biological reasons, immune cells should not be. If immune cells are also
## bigger / richer on TMA_TUMOR, part of the difference is the slide (run,
## section, staining, background) -- tumour and normal are on different slides,
## so slide and tissue type cannot be fully separated.
##
## INPUT : Xenium outs + data/segmentation/{cellpose,segger}/<slide>/ (whichever exist),
##         tables/02_cores_all.csv, config/tma_map.csv, config/annotation_families.csv
## OUTPUT: tables/03g_cells_per_core.csv        core x method x lineage medians
##         tables/03g_tumour_vs_normal.csv      all tumour vs all normal cores
##         tables/03g_tumour_vs_normal_by_site.csv
##         figures/03g_tumour_vs_normal.png     core medians, per method x lineage
##         figures/03g_tumour_vs_normal_by_site.png  epithelium + immune, matched organs
## =============================================================================

source("config/config.R")
source("R/helpers.R")

METHODS <- c("xenium", "cellpose", "segger")
METHOD_COLOURS <- c(xenium = "#00a6c7", cellpose = "#e0a800", segger = "#8e6bbf")
LINEAGES <- c("Epithelial", "Immune", "Mesenchymal", "Endothelial", "Neural_NE")
MIN_CELLS_CORE <- 20                # core x lineage medians need at least this many cells
ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))

slides <- read_slides()
tma_map <- read_tma_map()
cores_all <- utils::read.csv(file.path(TABLES_DIR, "02_cores_all.csv"), stringsAsFactors = FALSE)
cores_all <- cores_all[!is.na(cores_all$core_id), ]
markers <- read_marker_panel(MARKERS_RNA_CSV)

key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
fam <- read_config_csv(ANNOT_FAMILIES_CSV, stringsAsFactors = FALSE, fileEncoding = "latin1")
core_info <- tma_map[, intersect(c("core_id", "slide_id", "organ", "tissue_type"), colnames(tma_map))]
core_info$site <- if ("site" %in% colnames(fam)) fam$site[match(key(core_info$organ), key(fam$organ))] else NA
core_info$site[!is.na(core_info$site) & !nzchar(core_info$site)] <- NA
tt <- tolower(trimws(as.character(core_info$tissue_type)))
core_info$tumour <- ifelse(is.na(tt) | !nzchar(tt), grepl("tum", tolower(core_info$slide_id)), tt != "normal")

## -----------------------------------------------------------------------
## 1. per-cell table, one method x slide at a time (whole objects are large)
## -----------------------------------------------------------------------
lineage_of <- function(counts) {
  sets <- lapply(markers$level1[intersect(LINEAGES, names(markers$level1))],
                 function(g) intersect(match_case(g, rownames(counts)), rownames(counts)))
  sets <- sets[lengths(sets) > 0]
  raw <- sapply(sets, function(g) Matrix::colSums(counts[g, , drop = FALSE]))
  per_gene <- sweep(raw, 2, lengths(sets), "/")
  o <- t(apply(per_gene, 1, order, decreasing = TRUE))
  best <- per_gene[cbind(seq_len(nrow(per_gene)), o[, 1])]
  second <- per_gene[cbind(seq_len(nrow(per_gene)), o[, 2])]
  top_raw <- raw[cbind(seq_len(nrow(raw)), o[, 1])]
  ifelse(top_raw >= 2 & best >= 2 * second, colnames(per_gene)[o[, 1]], "unclear")
}

cells <- list()
for (mth in METHODS) {
  for (i in seq_len(nrow(slides))) {
    sid <- slides$slide_id[i]
    d <- seg_dir_for(slides[i, ], mth)
    if (!dir.exists(d) || (!dir.exists(file.path(d, "cell_feature_matrix")) &&
                           !file.exists(file.path(d, "cell_feature_matrix.h5")))) {
      message("  no ", mth, " output for ", sid, " -- skipped")
      next
    }
    o <- load_segmentation(slides[i, ], mth)
    m <- o[[]]
    m$core_id <- assign_cells_to_cores(m$x_um, m$y_um, cores_all[cores_all$slide_id == sid, ])
    keep <- !is.na(m$core_id) & m$nCount_Xenium >= MIN_COUNTS
    counts <- GetAssayData(o, assay = "Xenium", layer = "counts")[, keep, drop = FALSE]
    m <- m[keep, ]
    cells[[paste(mth, sid)]] <- data.frame(
      method = mth, slide_id = sid, core_id = m$core_id,
      cell_area = m$cell_area, nucleus_area = m$nucleus_area,
      counts = m$nCount_Xenium, genes = m$nFeature_Xenium,
      lineage = lineage_of(counts))
    message(sprintf("  %s / %s: %s cells in cores (>= %d counts)", mth, sid,
                    format(sum(keep), big.mark = ","), MIN_COUNTS))
    rm(o, counts); invisible(gc())
  }
}
cells <- dplyr::bind_rows(cells)
if (!nrow(cells)) stop("No segmentation output found.")
cells$density <- cells$counts / cells$cell_area
na_area <- cells |> group_by(method) |> summarise(na = mean(is.na(cell_area)), .groups = "drop")
for (k in which(na_area$na > 0.5))
  message("NOTE: ", na_area$method[k], " has no cell areas (segger without `segger export boundaries`) -- ",
          "areas and density are compared for the other methods only.")

## -----------------------------------------------------------------------
## 2. core medians (all cells and per lineage)
## -----------------------------------------------------------------------
METRICS <- c(cell_area = "cell area (um2)", nucleus_area = "nucleus area (um2)", counts = "transcripts per cell",
             genes = "genes per cell", density = "transcripts per um2 of cell")
summarise_cores <- function(d) d |>
  group_by(method, slide_id, core_id, lineage) |>
  summarise(n_cells = n(), across(all_of(names(METRICS)), \(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)),
            .groups = "drop") |>
  filter(n_cells >= MIN_CELLS_CORE)
per_core <- dplyr::bind_rows(summarise_cores(cells),
                             summarise_cores(dplyr::mutate(cells, lineage = "all cells")))
per_core <- per_core |> left_join(core_info[, c("core_id", "organ", "site", "tumour")], by = "core_id")
per_core$group <- ifelse(per_core$tumour, "tumour (TMA_TUMOR)", "normal (TMA_ORGAN)")
save_table(per_core, "03g_cells_per_core")

long <- tidyr::pivot_longer(per_core, all_of(names(METRICS)), names_to = "metric", values_to = "value") |>
  filter(!is.na(value), lineage %in% c("all cells", LINEAGES))

compare <- function(d, by) d |>
  group_by(across(all_of(by))) |>
  summarise(n_normal = sum(!tumour), n_tumour = sum(tumour),
            median_normal = stats::median(value[!tumour]), median_tumour = stats::median(value[tumour]),
            ratio_tumour_vs_normal = median_tumour / median_normal,
            p_wilcoxon = if (n_normal >= 3 && n_tumour >= 3)
              suppressWarnings(stats::wilcox.test(value[tumour], value[!tumour])$p.value) else NA_real_,
            .groups = "drop")
res <- compare(long, c("method", "lineage", "metric"))
res$p_adj <- stats::p.adjust(res$p_wilcoxon, "BH")
res <- res |> mutate(across(where(is.double), \(v) signif(v, 3)))
save_table(res, "03g_tumour_vs_normal")
options(width = 200)
message("\nTumour / normal ratio of core medians (all cells and the two key lineages):")
print(as.data.frame(res |> filter(lineage %in% c("all cells", "Epithelial", "Immune")) |>
                      select(method, lineage, metric, median_normal, median_tumour, ratio_tumour_vs_normal, p_adj)),
      row.names = FALSE)

## -----------------------------------------------------------------------
## 3. figures
## -----------------------------------------------------------------------
lin_show <- c("all cells", "Epithelial", "Immune", "Mesenchymal")
d <- long |> filter(lineage %in% lin_show)
d$lineage <- factor(d$lineage, levels = lin_show)
d$metric <- factor(METRICS[d$metric], levels = METRICS)
d$method <- factor(d$method, levels = METHODS)
lab <- res |> filter(lineage %in% lin_show) |>
  mutate(metric = factor(METRICS[metric], levels = METRICS), lineage = factor(lineage, levels = lin_show),
         method = factor(method, levels = METHODS),
         txt = sprintf("x%.2f%s", ratio_tumour_vs_normal,
                       ifelse(!is.na(p_adj) & p_adj < 0.05, "*", "")))
ymax <- d |> group_by(metric, lineage) |> summarise(y = max(value, na.rm = TRUE) * 1.08, .groups = "drop")
lab <- lab |> left_join(ymax, by = c("metric", "lineage"))
p <- ggplot(d, aes(method, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7) +
  geom_point(aes(colour = group), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.6, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = lab, aes(method, y, label = txt), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c("normal (TMA_ORGAN)" = "#BFD3E6", "tumour (TMA_TUMOR)" = "#F4B6A6"), name = NULL) +
  scale_colour_manual(values = c("normal (TMA_ORGAN)" = "#3B6FB6", "tumour (TMA_TUMOR)" = "#C0392B")) +
  facet_grid(metric ~ lineage, scales = "free_y", switch = "y") + theme_bw(base_size = 11) +
  theme(strip.placement = "outside", legend.position = "top") +
  labs(x = NULL, y = NULL,
       title = "Cell size and expression: TMA_TUMOR vs TMA_ORGAN (one dot = median of one core)",
       subtitle = paste0("label = tumour / normal ratio of the core medians (* BH p < 0.05, Wilcoxon on cores). ",
                         "Immune cells are the internal control: a difference there points to the slide, not the tumour."))
save_fig(p, "03g_tumour_vs_normal", width = 15, height = 16)

## organ-matched: tumour vs normal of the same site, epithelium and immune cells
ms <- long |> filter(!is.na(site), lineage %in% c("Epithelial", "Immune"),
                     metric %in% c("cell_area", "nucleus_area", "counts", "density"))
site_ok <- ms |> group_by(site) |> summarise(ok = any(tumour) && any(!tumour), .groups = "drop") |> filter(ok)
ms <- ms |> filter(site %in% site_ok$site)
if (nrow(ms)) {
  res_site <- compare(ms, c("site", "method", "lineage", "metric")) |>
    mutate(across(where(is.double), \(v) signif(v, 3)))
  save_table(res_site, "03g_tumour_vs_normal_by_site")
  ms$metric <- factor(METRICS[ms$metric], levels = METRICS)
  ms$method <- factor(ms$method, levels = METHODS)
  p2 <- ggplot(ms, aes(site, value, colour = group, shape = method)) +
    geom_point(position = position_dodge(0.7), size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = c("normal (TMA_ORGAN)" = "#3B6FB6", "tumour (TMA_TUMOR)" = "#C0392B"), name = NULL) +
    facet_grid(metric ~ lineage, scales = "free_y", switch = "y") + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.placement = "outside", legend.position = "top") +
    labs(x = NULL, y = NULL, title = "Organ-matched: tumour vs normal cores of the same site (one dot = one core)",
         subtitle = "same site, different slide: epithelium = tumour biology + slide; immune = mostly slide")
  save_fig(p2, "03g_tumour_vs_normal_by_site", width = max(10, 0.9 * nrow(site_ok) + 6), height = 13)
} else message("No site has both tumour and normal cores -- organ-matched comparison skipped ",
               "(check the `site` column of ", basename(ANNOT_FAMILIES_CSV), ").")
