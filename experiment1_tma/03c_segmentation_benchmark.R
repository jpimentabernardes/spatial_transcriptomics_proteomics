## =============================================================================
## 03c_segmentation_benchmark.R
##
## 10x (Xenium Onboard Analysis) vs cellpose vs segger on the SAME cores.
##
## There is no ground truth, so no single number decides. What matters for
## THIS project is that immune probes end up in immune cells and not smeared
## into neighbouring epithelium/stroma -- that is exactly what would fake
## "immune target expressed in organ X" in Aim 1. So the key metric is
## MECR (mutually exclusive co-expression rate: e.g. PTPRC+EPCAM+ cells), read
## together with sensitivity (transcripts assigned, counts/cell):
##   * a method that assigns more transcripts but raises MECR is buying
##     sensitivity with contamination;
##   * lower MECR at similar counts/cell = better.
## Also compare per ORGAN: nuclear-expansion segmentation performs very
## differently in dense lymphoid tissue vs large hepatocytes vs lung.
##
## Run on a subset of cores (BENCH_CORES) -- 2-3 per organ is enough and
## keeps three segmentations of whole TMAs out of memory.
##
## INPUT : Xenium outs + data/segmentation/{cellpose,segger}/<slide>/ (03a/03b)
##         tables/02_cores_all.csv (02)
## OUTPUT: tables/03c_segmentation_metrics_by_core.csv, figures/03c_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

METHODS <- c("xenium", "cellpose", "segger")
## ADAPT: benchmark cores; NULL = all cores present in every method's output.
BENCH_CORES <- NULL
ZOOM_CORE   <- NULL          # core shown side by side; NULL = first benchmark core
ZOOM_BOX_UM <- 250           # side length of the zoom window

slides <- read_slides()
tma_map <- read_tma_map()
cores_all <- utils::read.csv(file.path(TABLES_DIR, "02_cores_all.csv"), stringsAsFactors = FALSE)
cores_all <- cores_all[!is.na(cores_all$core_id), ]
markers <- read_marker_panel(MARKERS_RNA_CSV)

## -----------------------------------------------------------------------
## 1. Load each method, keep benchmark cores only
## -----------------------------------------------------------------------
load_method <- function(method) {
  objs <- list()
  for (i in seq_len(nrow(slides))) {
    d <- seg_dir_for(slides[i, ], method)
    if (!dir.exists(d)) { message("  no ", method, " output for ", slides$slide_id[i]); next }
    o <- load_segmentation(slides[i, ], method)
    cores_s <- cores_all[cores_all$slide_id == slides$slide_id[i], ]
    o$core_id <- assign_cells_to_cores(o$x_um, o$y_um, cores_s)
    keep <- !is.na(o$core_id)
    if (!is.null(BENCH_CORES)) keep <- keep & o$core_id %in% BENCH_CORES
    objs[[slides$slide_id[i]]] <- subset(o, cells = colnames(o)[keep])
  }
  if (length(objs) == 0) return(NULL)
  merge_slides(objs)
}

seg <- lapply(stats::setNames(METHODS, METHODS), load_method)
seg <- seg[!vapply(seg, is.null, logical(1))]
message("Methods with output: ", paste(names(seg), collapse = ", "))
stopifnot("Need at least two segmentations to compare" = length(seg) >= 2)

common_cores <- Reduce(intersect, lapply(seg, function(o) unique(o$core_id)))
message(length(common_cores), " cores present in all methods.")

## -----------------------------------------------------------------------
## 2. Metrics per core x method
## -----------------------------------------------------------------------
core_area_mm2 <- stats::setNames(pi * (cores_all$radius / 1000)^2, cores_all$core_id)

metrics <- dplyr::bind_rows(lapply(names(seg), function(method) {
  o <- seg[[method]]
  counts <- GetAssayData(o, assay = "Xenium", layer = "counts")
  m <- o[[]]
  ctrl_cols <- grep("^nCount_(ControlProbe|ControlCodeword|BlankCodeword)$", colnames(m), value = TRUE)
  m$ctrl <- if (length(ctrl_cols)) rowSums(m[, ctrl_cols, drop = FALSE]) else 0
  dplyr::bind_rows(lapply(common_cores, function(cid) {
    cells <- rownames(m)[m$core_id == cid]
    mm <- m[cells, ]
    data.frame(
      method = method,
      core_id = cid,
      n_cells = length(cells),
      cells_per_mm2 = length(cells) / core_area_mm2[[cid]],
      median_counts = median(mm$nCount_Xenium),
      median_genes = median(mm$nFeature_Xenium),
      median_area = median(mm$cell_area, na.rm = TRUE),
      pct_cells_min_counts = 100 * mean(mm$nCount_Xenium >= MIN_COUNTS),
      ctrl_per_1000 = 1000 * sum(mm$ctrl) / max(sum(mm$nCount_Xenium), 1),
      mecr = compute_mecr(counts[, cells, drop = FALSE], markers$mecr)
    )
  }))
}))

## Transcript assignment rate from the summaries written by 03a/03b. The 10x
## rate is recomputed there on the same transcripts, so all rates share a
## denominator (qv >= 20 gene transcripts in the core box).
read_summary <- function(method) {
  fs <- file.path(SEG_DIR, method, slides$slide_id, "assignment_summary.csv")
  fs <- fs[file.exists(fs)]
  if (length(fs) == 0) return(NULL)
  dplyr::bind_rows(lapply(fs, utils::read.csv))
}
summ <- dplyr::bind_rows(lapply(c("cellpose", "segger"), function(mth) {
  s <- read_summary(mth)
  if (is.null(s)) return(NULL)
  n_col <- paste0("n_assigned_", mth)
  denom <- if (paste0("n_tx_genes_", mth) %in% colnames(s)) s[[paste0("n_tx_genes_", mth)]] else s$n_tx_genes
  dplyr::bind_rows(
    data.frame(method = mth, core_id = s$region_id, pct_tx_assigned = 100 * s[[n_col]] / pmax(denom, 1)),
    data.frame(method = "xenium", core_id = s$region_id, pct_tx_assigned = 100 * s$n_assigned_xenium / pmax(s$n_tx_genes, 1))
  )
}))
if (nrow(summ) > 0) {
  summ <- summ |> distinct(method, core_id, .keep_all = TRUE)
  metrics <- metrics |> left_join(summ, by = c("method", "core_id"))
}

metrics <- metrics |> left_join(tma_map[, c("core_id", "organ", "tissue_type")], by = "core_id")
save_table(metrics, "03c_segmentation_metrics_by_core")

## -----------------------------------------------------------------------
## 3. Figures: paired per core, split by organ
## -----------------------------------------------------------------------
metric_cols <- intersect(c("median_counts", "pct_tx_assigned", "mecr", "median_area",
                           "cells_per_mm2", "ctrl_per_1000"), colnames(metrics))
long <- tidyr::pivot_longer(metrics, all_of(metric_cols), names_to = "metric", values_to = "value")
long$method <- factor(long$method, levels = intersect(METHODS, unique(long$method)))

p_metrics <- ggplot(long, aes(method, value)) +
  geom_line(aes(group = core_id, color = organ), alpha = 0.5) +
  geom_point(aes(color = organ), size = 1.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  ggtitle("Segmentation benchmark -- lines connect the same core across methods",
          subtitle = "Want: higher counts / assigned transcripts WITHOUT higher MECR")
save_fig(p_metrics, "03c_segmentation_metrics", width = 14, height = 8)

p_tradeoff <- ggplot(metrics, aes(median_counts, mecr, color = method)) +
  geom_point(size = 2) +
  geom_path(aes(group = core_id), color = "grey70", linewidth = 0.3) +
  facet_wrap(~ organ) + theme_bw() +
  labs(x = "median transcripts / cell", y = "MECR (lower = less spill-over)",
       title = "Sensitivity vs contamination trade-off, per organ")
save_fig(p_tradeoff, "03c_sensitivity_vs_mecr", width = 12, height = 8)

## -----------------------------------------------------------------------
## 4. Side-by-side zoom: simple lineage call per cell (argmax of mean marker
##    counts) -- enough to SEE mixed cells at epithelium/immune borders.
## -----------------------------------------------------------------------
zoom_core <- ZOOM_CORE %||% common_cores[1]
lin <- lapply(markers$level1, intersect, y = rownames(seg[[1]]))
lin <- lin[lengths(lin) > 0]
zoom <- dplyr::bind_rows(lapply(names(seg), function(method) {
  o <- seg[[method]]
  m <- o[[]]
  cells <- rownames(m)[m$core_id == zoom_core]
  cx <- median(m[cells, "x_um"]); cy <- median(m[cells, "y_um"])
  inside <- cells[abs(m[cells, "x_um"] - cx) < ZOOM_BOX_UM / 2 & abs(m[cells, "y_um"] - cy) < ZOOM_BOX_UM / 2]
  cnt <- GetAssayData(o, assay = "Xenium", layer = "counts")[, inside, drop = FALSE]
  sc <- sapply(lin, function(g) Matrix::colMeans(cnt[g, , drop = FALSE]))
  sc <- matrix(sc, nrow = length(inside), dimnames = list(inside, names(lin)))
  call <- ifelse(rowSums(sc) == 0, "none", colnames(sc)[max.col(sc, ties.method = "first")])
  data.frame(method = method, x_um = m[inside, "x_um"], y_um = m[inside, "y_um"],
             lineage = call, counts = m[inside, "nCount_Xenium"])
}))
p_zoom <- ggplot(zoom, aes(x_um, y_um, color = lineage, size = counts)) +
  geom_point(alpha = 0.8) + scale_size_continuous(range = c(0.3, 2.5)) +
  scale_y_reverse() + coord_fixed() + facet_wrap(~ method) + theme_bw() +
  ggtitle(paste0(zoom_core, ": ", ZOOM_BOX_UM, " um window -- compare cell density and lineage mixing"))
save_fig(p_zoom, "03c_zoom_side_by_side", width = 6 * length(seg), height = 6)

## -----------------------------------------------------------------------
## 5. Summary by method x organ -> decide, then set SEGMENTATION in config.R
## -----------------------------------------------------------------------
decision <- metrics |>
  group_by(method, organ) |>
  summarise(across(all_of(metric_cols), \(v) round(median(v, na.rm = TRUE), 3)), .groups = "drop") |>
  arrange(organ, method)
print(as.data.frame(decision))
save_table(decision, "03c_segmentation_summary")

## >>> DECISION <<<
## Pick ONE segmentation for all slides (mixing methods between organs makes
## Aim 1/2 comparisons uninterpretable). Typical outcome: 10x multimodal is
## fine for epithelium; segger/cellpose help most in dense lymphoid or
## immune-rich cores. Whatever you pick: set SEGMENTATION in config/config.R,
## re-run 01 and 02, and state the choice + these numbers in your methods.
