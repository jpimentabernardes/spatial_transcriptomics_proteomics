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

for (col in c("nCount_controls", "control_frac", "counts_per_um2", "area_upper",
              "qc_low_counts", "qc_small", "qc_large", "qc_high_control", "qc_pass")) {
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
p_vln <- VlnPlot(xen, features = c("nCount_Xenium", "nFeature_Xenium", "cell_area", "control_frac"),
                 group.by = "slide_id", pt.size = 0, ncol = 4, log = FALSE)
save_fig(p_vln, "01_qc_violin_by_slide", width = 16, height = 5)

p_spatial <- plot_cores_spatial(
  meta |> mutate(log10_counts = log10(nCount_Xenium + 1)),
  color_by = "log10_counts", facet = FALSE
) + facet_wrap(~ slide_id) + theme(aspect.ratio = NULL) +
  ggtitle("log10 transcripts per cell (look for whole cores that are dim = tissue/run problem)")
save_fig(p_spatial, "01_counts_spatial_by_slide", width = 8 * nrow(slides), height = 8)

p_pass <- plot_cores_spatial(meta, color_by = "qc_pass", facet = FALSE) +
  facet_wrap(~ slide_id) + scale_color_manual(values = c(`TRUE` = "grey60", `FALSE` = "red"))
save_fig(p_pass, "01_qc_fail_spatial_by_slide", width = 8 * nrow(slides), height = 8)

## >>> CHECK <<<
## If failing cells cluster at core EDGES that is expected (cut cells). If a
## WHOLE core fails, that core is a tissue-quality issue (necrosis, fixation,
## detachment) -- keep it flagged in 02 rather than silently letting the
## cell filter delete it; for Aim 2 a missing core is a missing sample.

saveRDS(xen, obj_path("raw"))
message("Saved -> ", obj_path("raw"))
