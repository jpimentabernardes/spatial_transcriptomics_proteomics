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
## INPUT : Xenium outs + data/segmentation/{cellpose,segger,segger_xenium,segger_cellpose}/<slide>/
##         (03a / 03b / 03d; whichever exist)
##         tables/02_cores_all.csv (02)
## OUTPUT: tables/03c_segmentation_metrics_by_core.csv, 03c_scorecard.csv, 03c_segmentation_summary.csv
##         figures/03c_cell_numbers.png      cells called per method (total + per core)
##         figures/03c_scorecard.png         per measure, which method is better in how many cores
##         figures/03c_side_by_side_<core>.png  DAPI | 10x | cellpose | overlay, as in 01
## =============================================================================

source("config/config.R")
source("R/helpers.R")

## methods without output are skipped; segger_xenium / segger_cellpose = segger
## started from the 10x / the cellpose cells (03d_segger_priors.py)
METHODS <- c("xenium", "cellpose", "segger", "segger_xenium", "segger_cellpose")
METHOD_COLOURS <- c(xenium = "#00a6c7", cellpose = "#e0a800", segger = "#8e6bbf",
                    segger_xenium = "#00596b", segger_cellpose = "#8f5f00")
## ADAPT: benchmark cores; NULL = all cores present in every method's output.
BENCH_CORES <- NULL
ZOOM_CORES  <- NULL          # cores shown side by side (window at the core centre); NULL = first 3 common cores
ZOOM_BOX_UM <- 150           # side length of the zoom windows (um)

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
    if (!any(keep)) {
      message("  ", method, " / ", slides$slide_id[i], ": no cells inside the cores of 02_cores_all.csv -- skipped")
      next
    }
    objs[[slides$slide_id[i]]] <- subset(o, cells = colnames(o)[keep])
  }
  if (length(objs) == 0) return(NULL)
  merge_slides(objs)
}

seg <- lapply(stats::setNames(METHODS, METHODS), load_method)
seg <- seg[!vapply(seg, is.null, logical(1))]
message("Methods with output: ", paste(names(seg), collapse = ", "))
stopifnot("Need at least two segmentations to compare" = length(seg) >= 2)

## MECR uses the use_for_mecr genes of marker_panel_rna.csv that are on this panel:
## pairs of genes from DIFFERENT lineages should never be in the same cell
mecr_genes <- lapply(markers$mecr, function(g) intersect(match_case(g, rownames(seg[[1]])), rownames(seg[[1]])))
mecr_genes <- mecr_genes[lengths(mecr_genes) > 0]
message("MECR marker genes on the panel: ",
        paste(sprintf("%s = %s", names(mecr_genes), vapply(mecr_genes, paste, "", collapse = "/")), collapse = "; "))
if (length(mecr_genes) < 2) {
  message("  fewer than 2 lineages with MECR markers on the panel -> MECR cannot be computed; ",
          "add panel genes to use_for_mecr in ", basename(MARKERS_RNA_CSV))
} else if (!"Immune" %in% names(mecr_genes)) {
  message("  no IMMUNE MECR marker on the panel -> MECR does not test immune spill-over; ",
          "add immune panel genes to use_for_mecr in ", basename(MARKERS_RNA_CSV))
}

common_cores <- Reduce(intersect, lapply(seg, function(o) unique(o$core_id)))
message(length(common_cores), " cores present in all methods.")
for (mth in names(seg)) message("  ", mth, ": ", length(unique(seg[[mth]]$core_id)), " cores")
if (length(common_cores) == 0) {
  stop("No core was segmented by all methods -- e.g. 03a was run on other cores (cores = [...]) ",
       "or on another slide. Run the methods on the same cores, or set BENCH_CORES.")
}

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
summ <- dplyr::bind_rows(lapply(setdiff(names(seg), "xenium"), function(mth) {
  s <- read_summary(mth)
  n_col <- paste0("n_assigned_", mth)
  if (is.null(s) || !n_col %in% colnames(s)) return(NULL)
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
## 3. Cell numbers, directly: total per slide, and per core (same cores)
## -----------------------------------------------------------------------
metrics$method <- factor(metrics$method, levels = intersect(METHODS, unique(metrics$method)))
metrics$slide_id <- cores_all$slide_id[match(metrics$core_id, cores_all$core_id)]
tot <- metrics |> group_by(slide_id, method) |> summarise(n_cells = sum(n_cells), .groups = "drop")
print(as.data.frame(tot))
p_tot <- ggplot(tot, aes(slide_id, n_cells, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = format(n_cells, big.mark = ",")), position = position_dodge(0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = METHOD_COLOURS) +
  theme_bw() + labs(x = NULL, y = "cells in the compared cores", title = "Cells called per method")

wide <- tidyr::pivot_wider(metrics[, c("core_id", "slide_id", "method", "n_cells")],
                           names_from = method, values_from = n_cells)
others <- setdiff(levels(metrics$method), "xenium")
p_cores <- lapply(others, function(mth) {
  lims <- range(c(0, wide$xenium, wide[[mth]]), na.rm = TRUE)      # same scale on both axes: diagonal = equal
  ggplot(wide, aes(xenium, .data[[mth]], colour = slide_id)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 2) + coord_equal(xlim = lims, ylim = lims) + theme_bw() +
    labs(x = "cells per core: 10x Xenium", y = paste("cells per core:", mth),
         title = sprintf("%s / Xenium cells: median ratio %.2f", mth, stats::median(wide[[mth]] / wide$xenium)))
})
save_fig(patchwork::wrap_plots(c(list(p_tot), p_cores), nrow = 1), "03c_cell_numbers",
         width = 6 * (1 + length(others)), height = 5.5)

## -----------------------------------------------------------------------
## 4. Which is better? Per measure, paired per core: how many cores each
##    method wins. Higher is better for transcripts assigned to cells,
##    counts per cell and % cells passing the count threshold; lower is
##    better for MECR (marker pairs that should never share a cell = merged
##    or contaminated cells) and negative-control rate. More cells is NOT
##    better by itself: split nuclei and debris also add cells.
## -----------------------------------------------------------------------
better <- c(pct_tx_assigned = "higher", median_counts = "higher", pct_cells_min_counts = "higher",
            mecr = "lower", ctrl_per_1000 = "lower")
better <- better[names(better) %in% colnames(metrics)]
score <- dplyr::bind_rows(lapply(names(better), function(k) {
  w <- tidyr::pivot_wider(metrics[, c("core_id", "method", k)], names_from = method, values_from = all_of(k))
  dplyr::bind_rows(lapply(others, function(mth) {
    d <- w[[mth]] - w$xenium
    ok <- !is.na(d)
    wins <- if (better[[k]] == "higher") sum(d[ok] > 0) else sum(d[ok] < 0)
    data.frame(measure = k, better_if = better[[k]], method = mth, n_cores = sum(ok),
               xenium_median = round(stats::median(w$xenium, na.rm = TRUE), 3),
               method_median = round(stats::median(w[[mth]], na.rm = TRUE), 3),
               method_wins = wins, xenium_wins = sum(ok) - wins - sum(d[ok] == 0),
               winner = ifelse(wins > sum(ok) / 2, mth, ifelse(sum(ok) - wins - sum(d[ok] == 0) > sum(ok) / 2, "xenium", "tie")))
  }))
}))
print(score)
save_table(score, "03c_scorecard")

long <- tidyr::pivot_longer(metrics, all_of(c("n_cells", names(better))), names_to = "measure", values_to = "value")
lab <- score |> group_by(measure) |>
  summarise(lab = paste0(measure[1], " (", better_if[1], " = better): ",
                         paste(sprintf("xenium better in %d, %s in %d of %d cores", xenium_wins, method,
                                       method_wins, n_cores), collapse = "; ")), .groups = "drop")
lab <- rbind(lab, data.frame(measure = "n_cells", lab = "n_cells (descriptive: more is not better by itself)"))
long$panel <- lab$lab[match(long$measure, lab$measure)]
p_score <- ggplot(long, aes(method, value, group = core_id)) +
  geom_line(colour = "grey70") + geom_point(aes(colour = method), size = 2) +
  scale_colour_manual(values = METHOD_COLOURS, guide = "none") +
  facet_wrap(~ panel, scales = "free_y", ncol = 2, labeller = label_wrap_gen(55)) +
  theme_bw() + labs(x = NULL, y = NULL, title = "Segmentation scorecard -- each line is one core")
save_fig(p_score, "03c_scorecard", width = 12, height = 4 * ceiling((length(better) + 1) / 2))

## -----------------------------------------------------------------------
## 5. Side by side on the image, as in 01: DAPI | 10x | cellpose | overlay,
##    for a window at the centre of a few compared cores. cellpose = 03a's
##    own result where its masks were saved (03a save_masks = TRUE), else
##    cellpose run on the window with the 01 inspection settings.
## -----------------------------------------------------------------------
if (!requireNamespace("tiff", quietly = TRUE)) {
  message("Side-by-side images skipped: install.packages('tiff').")
} else {
  zoom_cores <- ZOOM_CORES %||% head(common_cores, 3)
  for (cid in zoom_cores) {
    cr <- cores_all[cores_all$core_id == cid, ][1, ]
    h <- ZOOM_BOX_UM / 2
    win <- data.frame(x0 = cr$x_center - h, y0 = cr$y_center - h, x1 = cr$x_center + h, y1 = cr$y_center + h)
    sid <- cr$slide_id
    cp_opts <- if ("cellpose" %in% names(seg)) {
      list(mode = INSPECT_CELLPOSE_MODE, expand_um = INSPECT_CELLPOSE_EXPAND_UM, max_area = MAX_CELL_AREA,
           masks_dir = file.path(SEG_DIR, "cellpose", sid, "masks"))
    } else NULL
    crop <- tryCatch(
      crop_morphology(resolve_path(slides$xenium_dir[slides$slide_id == sid]), win$x0, win$y0, win$x1, win$y1,
                      paste0("03c_", cid), cellpose = cp_opts),
      error = function(e) { message("Window ", cid, ": ", conditionMessage(e)); NULL })
    if (is.null(crop)) next
    xm <- seg$xenium[[]]
    p_side <- plot_segmentation_compare(crop, win, xm[xm$slide_id == sid, ],
                                        title = sprintf("%s: %d um at the core centre", cid, ZOOM_BOX_UM))
    save_fig(p_side, paste0("03c_side_by_side_", cid), width = 20, height = 5.8, dpi = 200)
  }
}

## -----------------------------------------------------------------------
## 5. Summary by method x organ -> decide, then set SEGMENTATION in config.R
## -----------------------------------------------------------------------
metric_cols <- intersect(c("n_cells", "median_counts", "pct_tx_assigned", "pct_cells_min_counts", "mecr",
                           "median_area", "ctrl_per_1000"), colnames(metrics))
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
