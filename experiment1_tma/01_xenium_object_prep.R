## =============================================================================
## 01_xenium_object_prep.R
##
## QC and object preparation for all Xenium TMA slides.
##
## OBJECT TYPE -- why this structure (decided once here, used everywhere):
##   * One Seurat v5 object for ALL slides, one centroid FOV per slide
##     (FOV name = slide_id). Continuity with the class scripts, and Seurat is
##     where the annotation/statistics tooling you already know lives.
##   * Cells renamed <slide_id>_<cell_id>: 10x cell ids repeat across slides.
##   * Genes in assay "Xenium"; negative-control probes/codewords kept as
##     per-cell metadata (nCount_ControlProbe, ...). They are the background
##     model for Aim 1 (is a probe's signal above noise in organ X?), so they
##     are never thrown away.
##   * Centroids only, no cell polygons. Polygons for multi-slide TMAs cost
##     many GB and are not needed for any statistic here; plot them in Xenium
##     Explorer when you need to look at actual cell shapes.
##   * TMA metadata (core, patient, organ, tissue type, location) is added per
##     cell in 02 from config/tma_map.csv -- the TMA map, not the file name, is
##     the source of truth for sample identity.
##   * Python tools (cellpose, segger, Spateo, GEASO) exchange data with this
##     object through plain files (10x-style matrices, CSV of centroids), not
##     through format converters -- fewer moving parts to break.
##
## The same script loads any segmentation (config: SEGMENTATION), so after the
## benchmark in 03c you switch segmentation and re-run 01 -> 02 unchanged.
##
## INPUT : config/slides.csv, Xenium outs (or 03a/03b segmentation outputs)
## OUTPUT: data/processed/xenium_<seg>_raw.rds  (all cells, QC flags, not filtered)
##         tables/01_run_metrics.csv, tables/01_slide_background.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

slides <- read_slides()
print(slides[, c("slide_id", "tma_type", "xenium_dir")])

## -----------------------------------------------------------------------
## 1. Run-level metrics from Xenium Onboard Analysis (metrics_summary.csv)
##    Look at these BEFORE any cell-level QC: a slide with low decoded
##    transcripts / high negative-control rate is a run problem, not
##    something cell filtering can fix.
## -----------------------------------------------------------------------
run_metrics <- lapply(seq_len(nrow(slides)), function(i) {
  f <- file.path(resolve_path(slides$xenium_dir[i]), "metrics_summary.csv")
  if (!file.exists(f)) return(NULL)
  m <- utils::read.csv(f, check.names = FALSE, colClasses = "character")
  m$slide_id <- slides$slide_id[i]
  m
})
run_metrics <- dplyr::bind_rows(run_metrics)
if (nrow(run_metrics) > 0) {
  run_metrics <- utils::type.convert(run_metrics, as.is = TRUE)
  save_table(run_metrics, "01_run_metrics")
  ## Columns differ between XOA versions -- print what is there, pick by name.
  key_cols <- intersect(
    c("slide_id", "num_cells_detected", "fraction_transcripts_assigned",
      "median_transcripts_per_cell", "decoded_transcripts_per_100um2",
      "estimated_number_of_false_positive_transcripts_per_cell",
      "fraction_transcripts_decoded_q20", "negative_control_probe_counts_per_control_per_cell",
      "negative_decoder_codeword_counts_per_control_per_cell"),
    colnames(run_metrics)
  )
  print(run_metrics[, key_cols, drop = FALSE])
}

## -----------------------------------------------------------------------
## 2. Load every slide with an identical structure and merge
## -----------------------------------------------------------------------
objs <- lapply(seq_len(nrow(slides)), function(i) load_segmentation(slides[i, ], SEGMENTATION))
names(objs) <- slides$slide_id
xen <- merge_slides(objs)
rm(objs); invisible(gc())

xen$tma_type <- slides$tma_type[match(xen$slide_id, slides$slide_id)]
print(xen)
print(Images(xen))   # one FOV per slide
print(table(xen$slide_id))

## -----------------------------------------------------------------------
## 3. Per-cell QC metrics
## -----------------------------------------------------------------------
meta <- xen[[]]
ctrl_cols <- grep("^nCount_(ControlProbe|ControlCodeword|BlankCodeword)$", colnames(meta), value = TRUE)
meta$nCount_controls <- if (length(ctrl_cols) > 0) rowSums(meta[, ctrl_cols, drop = FALSE]) else 0
meta$control_frac <- meta$nCount_controls / pmax(meta$nCount_Xenium + meta$nCount_controls, 1)
meta$counts_per_um2 <- meta$nCount_Xenium / meta$cell_area

## Per-slide upper area cut: merged cells/doublets from segmentation. MAD on
## log(area) because area is right-skewed.
meta <- meta |>
  group_by(slide_id) |>
  mutate(
    area_upper = exp(median(log(cell_area), na.rm = TRUE) +
                       AREA_MAD_UPPER * mad(log(cell_area), na.rm = TRUE)),
    area_upper = pmin(area_upper, MAX_CELL_AREA)
  ) |>
  ungroup() |>
  as.data.frame()
rownames(meta) <- colnames(xen)

## Flags, not deletions: each reason kept separately so you can see WHY cells
## go, per core, in 02. Cells with missing area (some segmentations) are not
## penalised on area.
meta$qc_low_counts   <- meta$nCount_Xenium < MIN_COUNTS | meta$nFeature_Xenium < MIN_FEATURES
meta$qc_small        <- !is.na(meta$cell_area) & meta$cell_area < MIN_CELL_AREA
meta$qc_large        <- !is.na(meta$cell_area) & meta$cell_area > meta$area_upper
meta$qc_high_control <- meta$control_frac > MAX_CONTROL_FRAC
meta$qc_pass <- !(meta$qc_low_counts | meta$qc_small | meta$qc_large | meta$qc_high_control)
meta$qc_reason <- qc_reason(meta)          # first failing criterion, or "pass"

for (col in c("nCount_controls", "control_frac", "counts_per_um2", "area_upper",
              "qc_low_counts", "qc_small", "qc_large", "qc_high_control", "qc_pass", "qc_reason")) {
  xen[[col]] <- meta[[col]]
}

qc_summary <- meta |>
  group_by(slide_id) |>
  summarise(
    n_cells = n(),
    median_counts = median(nCount_Xenium),
    median_genes = median(nFeature_Xenium),
    median_area = median(cell_area, na.rm = TRUE),
    pct_low_counts = 100 * mean(qc_low_counts),
    pct_small = 100 * mean(qc_small),
    pct_large = 100 * mean(qc_large),
    pct_high_control = 100 * mean(qc_high_control),
    pct_pass = 100 * mean(qc_pass),
    .groups = "drop"
  )
print(qc_summary)
save_table(qc_summary, "01_cell_qc_by_slide")

## -----------------------------------------------------------------------
## 4. Slide-level background: counts per negative-control probe vs counts per
##    gene. The ratio is a rough false-discovery estimate for the whole slide;
##    05 repeats this per core and per gene, which is what Aim 1 needs.
## -----------------------------------------------------------------------
bg <- meta |>
  group_by(slide_id) |>
  summarise(
    gene_counts = sum(nCount_Xenium),
    n_genes = first(nFeatTotal_Gene),
    ctrl_probe_counts = if ("nCount_ControlProbe" %in% colnames(meta)) sum(nCount_ControlProbe) else NA,
    n_ctrl_probes = if ("nFeatTotal_ControlProbe" %in% colnames(meta)) first(nFeatTotal_ControlProbe) else NA,
    ctrl_codeword_counts = if ("nCount_ControlCodeword" %in% colnames(meta)) sum(nCount_ControlCodeword) else NA,
    n_ctrl_codewords = if ("nFeatTotal_ControlCodeword" %in% colnames(meta)) first(nFeatTotal_ControlCodeword) else NA,
    .groups = "drop"
  ) |>
  mutate(
    counts_per_gene = gene_counts / n_genes,
    counts_per_ctrl_probe = ctrl_probe_counts / n_ctrl_probes,
    counts_per_ctrl_codeword = ctrl_codeword_counts / n_ctrl_codewords,
    background_ratio_probe = counts_per_ctrl_probe / counts_per_gene,
    background_ratio_codeword = counts_per_ctrl_codeword / counts_per_gene
  )
print(bg)
save_table(bg, "01_slide_background")

## -----------------------------------------------------------------------
## 5. Figures
## -----------------------------------------------------------------------
## Violins on a log scale (a handful of huge or very bright cells otherwise
## squash everything), with the QC thresholds as dashed lines and the % of
## cells each threshold removes.
qc_violin <- function(col, label, trans, breaks, lo = NULL, hi = NULL) {
  pct_lo <- if (!is.null(lo)) 100 * mean(meta[[col]] < lo, na.rm = TRUE) else NA
  pct_hi <- if (!is.null(hi)) 100 * mean(meta[[col]] > hi, na.rm = TRUE) else NA
  sub <- paste(c(if (!is.na(pct_lo)) sprintf("%.1f%% below %g", pct_lo, lo),
                 if (!is.na(pct_hi)) sprintf("%.1f%% above %g", pct_hi, hi)), collapse = ", ")
  ggplot(meta, aes(slide_id, .data[[col]], fill = slide_id)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    { if (!is.null(lo)) geom_hline(yintercept = lo, linetype = 2) } +
    { if (!is.null(hi)) geom_hline(yintercept = hi, linetype = 2) } +
    scale_y_continuous(trans = trans, breaks = breaks) +
    scale_fill_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")[seq_len(nrow(slides))],
                      guide = "none") +
    theme_bw(base_size = 10) + labs(x = NULL, y = NULL, title = label, subtitle = sub)
}
p_vln <- qc_violin("nCount_Xenium", "transcripts per cell", "log1p", c(0, 3, 10, 30, 100, 300, 1000, 3000), lo = MIN_COUNTS) |
  qc_violin("nFeature_Xenium", "genes per cell", "log1p", c(0, 2, 5, 10, 20, 50, 100), lo = MIN_FEATURES) |
  qc_violin("cell_area", "cell area (um2)", "log1p", c(0, 5, 10, 30, 100, 300, 1000, 3000),
            lo = MIN_CELL_AREA, hi = MAX_CELL_AREA) |
  qc_violin("control_frac", "negative-control fraction", "sqrt", c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1),
            hi = MAX_CONTROL_FRAC)
save_fig(p_vln, "01_qc_violin_by_slide", width = 16, height = 5)

## Spatial maps with micron axes: read coordinates here for INSPECT_CENTERS
## (section 6). QC_MAP_POINT_SIZE / QC_MAP_MAX_CELLS are set in config.R.
p_spatial <- plot_cores_spatial(
  meta |> mutate(log10_counts = log10(nCount_Xenium + 1)),
  color_by = "log10_counts", facet = FALSE, axes = TRUE,
  size = QC_MAP_POINT_SIZE, max_cells = QC_MAP_MAX_CELLS
) + facet_wrap(~ slide_id) +
  ggtitle("log10 transcripts per cell (look for whole cores that are dim = tissue/run problem)")
save_fig(p_spatial, "01_counts_spatial_by_slide", width = 8 * nrow(slides), height = 8)

## Coloured by WHY a cell fails -- tells you which threshold to look at.
p_pass <- plot_cores_spatial(meta[order(meta$qc_reason != "pass"), ],      # failing cells drawn on top
                             color_by = "qc_reason", facet = FALSE, axes = TRUE,
                             size = QC_MAP_POINT_SIZE, max_cells = QC_MAP_MAX_CELLS) +
  facet_wrap(~ slide_id) +
  scale_color_manual(values = c(QC_REASON_COLORS[-1], pass = "grey80"), drop = FALSE, name = "cell QC") +
  ggtitle("Cell QC: pass (grey) and the first criterion each failing cell fails")
save_fig(p_pass, "01_qc_fail_spatial_by_slide", width = 8 * nrow(slides), height = 8)

## >>> CHECK <<<
## If failing cells cluster at core EDGES that is expected (cut cells). If a
## WHOLE core fails, that core is a tissue-quality issue (necrosis, fixation,
## detachment) -- keep it flagged in 02 rather than silently letting the
## cell filter delete it; for Aim 2 a missing core is a missing sample.

## -----------------------------------------------------------------------
## 6. Zoomed inspection: QC pass / fail cells on the morphology image
##    For each window: where it is on the slide | the DAPI image alone | DAPI
##    with segmentation outlines coloured by QC result. Look for: failing
##    cells in tissue gaps / folds (debris = correctly removed), large cells
##    spanning several nuclei (merged segmentation), small cells without a
##    clear nucleus, or good-looking nuclei that fail only on counts (a
##    threshold that is too strict for this tissue).
##    Windows are picked automatically where pass and fail cells mix; set
##    INSPECT_CENTERS in config.R to look at a specific place instead.
##    With INSPECT_CELLPOSE = TRUE a second row shows the same window
##    segmented by cellpose (03a settings, same QC rules): QC counts of both,
##    10x vs cellpose outlines on the image, and cellpose cells by QC result.
## -----------------------------------------------------------------------
if (!requireNamespace("tiff", quietly = TRUE)) {
  message("Section 6 skipped: install.packages('tiff') to read the image crops.")
} else {
  xm <- xen[[]]
  xm$cell_id <- as.character(xm$cell_id)
  windows <- if (!is.null(INSPECT_CENTERS)) {
    h <- INSPECT_WINDOW_UM / 2
    with(INSPECT_CENTERS, data.frame(slide_id = slide_id, x0 = x - h, y0 = y - h, x1 = x + h, y1 = y + h))
  } else pick_inspection_windows(xm, INSPECT_WINDOW_UM, INSPECT_N_PER_SLIDE)
  print(windows)
  save_table(windows, "01_inspection_windows")

  for (k in seq_len(nrow(windows))) {
    w <- windows[k, ]
    sid <- w$slide_id
    tag <- sprintf("%s_x%d_y%d", sid, round(w$x0), round(w$y0))
    cp_opts <- if (isTRUE(INSPECT_CELLPOSE)) {
      list(mode = INSPECT_CELLPOSE_MODE, expand_um = INSPECT_CELLPOSE_EXPAND_UM,
           max_area = xm$area_upper[xm$slide_id == sid][1])
    } else NULL
    crop <- tryCatch(
      crop_morphology(resolve_path(slides$xenium_dir[slides$slide_id == sid]),
                      w$x0, w$y0, w$x1, w$y1, tag, cellpose = cp_opts),
      error = function(e) { message("Window ", tag, ": ", conditionMessage(e)); NULL })
    if (is.null(crop)) next
    p_zoom <- plot_inspection(xm[xm$slide_id == sid, ], crop, w,
                              title = sprintf("%s  x %.0f-%.0f, y %.0f-%.0f um", sid, w$x0, w$x1, w$y0, w$y1))
    save_fig(p_zoom, paste0("01_qc_zoom_", tag), width = 15, height = if (is.null(crop$cp_cells)) 5.5 else 11,
             dpi = 200)
  }
}

saveRDS(xen, obj_path("raw"))
message("Saved -> ", obj_path("raw"))
