## =============================================================================
## 02_tma_dearray.R
##
## Assigns every cell to a TMA core and attaches the sample metadata
## (patient, organ, tissue type, diagnosis, location) from config/tma_map.csv.
## Then does core-level QC and applies the cell filter from 01.
##
## Two ways to define cores (per slide, picked automatically):
##   A. Polygon mode -- if config/core_selections/<slide_id>/ exists: one CSV
##      per core (file name = core_id) with X,Y vertices in microns, e.g.
##      drawn with Xenium Explorer's lasso and "Download selection
##      coordinates". Most accurate; use it for torn/irregular tumour cores.
##   B. Automatic mode -- DBSCAN on cell centroids finds cores, rows/columns
##      are inferred from the core positions and matched to core_row/core_col
##      in tma_map.csv. Fast, but the grid inference MUST be checked on the
##      labelled plot (missing cores or a rotated TMA shift labels). Fix any
##      mistakes with config/core_overrides.csv:
##         slide_id,detected_core,core_row,core_col
##
## Core definitions are saved once (tables/02_cores_<slide>.csv) and reused
## when you re-run with another segmentation, so core identity never depends
## on the segmentation method. DELETE that file to re-dearray (e.g. after
## editing core_overrides.csv or the DBSCAN parameters). The same file drives per-core cropping in 03a
## and per-core registration in 07.
##
## INPUT : xenium_<seg>_raw.rds (01), config/tma_map.csv
## OUTPUT: xenium_<seg>_cores.rds (cells in included cores, QC-passing)
##         tables/02_cores_<slide>.csv, tables/02_core_qc.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen <- readRDS(obj_path("raw"))
slides <- read_slides()
tma_map <- read_tma_map()
meta <- xen[[]]

overrides_path <- file.path(CONFIG_DIR, "core_overrides.csv")
overrides <- if (file.exists(overrides_path)) {
  read_config_csv(overrides_path, stringsAsFactors = FALSE, colClasses = "character")
} else NULL

meta$core_id <- NA_character_
all_cores <- list()

for (i in seq_len(nrow(slides))) {
  sid <- slides$slide_id[i]
  in_slide <- meta$slide_id == sid
  cores_file <- file.path(TABLES_DIR, paste0("02_cores_", sid, ".csv"))
  sel_dir <- file.path(CONFIG_DIR, "core_selections", sid)
  map_s <- tma_map[tma_map$slide_id == sid, ]
  if (nrow(map_s) == 0) stop("No rows for slide ", sid, " in tma_map.csv")

  if (dir.exists(sel_dir)) {
    ## ---- A. polygon mode --------------------------------------------------
    polys <- read_core_selections(sel_dir)
    unknown <- setdiff(names(polys), map_s$core_id)
    if (length(unknown) > 0) warning(sid, ": selection files not in tma_map: ", paste(unknown, collapse = ", "))
    meta$core_id[in_slide] <- assign_cells_to_polygons(meta$x_um[in_slide], meta$y_um[in_slide], polys)
    cores <- do.call(rbind, lapply(names(polys), function(nm) {
      p <- polys[[nm]]
      cx <- mean(p$x); cy <- mean(p$y)
      data.frame(core_id = nm, x_center = cx, y_center = cy,
                 radius = max(sqrt((p$x - cx)^2 + (p$y - cy)^2)),
                 xmin = min(p$x), xmax = max(p$x), ymin = min(p$y), ymax = max(p$y))
    }))
    message(sid, ": polygon mode, ", nrow(cores), " cores.")

  } else if (file.exists(cores_file)) {
    ## ---- previously defined cores (e.g. re-run with another segmentation) --
    cores <- utils::read.csv(cores_file, stringsAsFactors = FALSE)
    cores <- cores[!is.na(cores$core_id), ]
    meta$core_id[in_slide] <- assign_cells_to_cores(meta$x_um[in_slide], meta$y_um[in_slide], cores)
    message(sid, ": reusing core definitions from ", cores_file)

  } else {
    ## ---- B. automatic dearraying ------------------------------------------
    ## Use cells with some signal so empty debris segments don't bridge cores.
    use <- in_slide & meta$nCount_Xenium >= MIN_COUNTS
    cores <- detect_cores(meta$x_um[use], meta$y_um[use])
    cores <- assign_grid(
      cores,
      row_levels = sort_levels(map_s$core_row),
      col_levels = sort_levels(map_s$core_col),
      row_flip = isTRUE(slides$row_flip[i]),
      col_flip = isTRUE(slides$col_flip[i])
    )
    if (!is.null(overrides)) {
      ov <- overrides[overrides$slide_id == sid, ]
      for (j in seq_len(nrow(ov))) {
        k <- cores$detected_core == as.integer(ov$detected_core[j])
        cores$core_row[k] <- ov$core_row[j]
        cores$core_col[k] <- ov$core_col[j]
      }
    }
    cores <- cores |>
      left_join(map_s[, c("core_row", "core_col", "core_id")], by = c("core_row", "core_col"))
    dups <- cores$core_id[!is.na(cores$core_id) & duplicated(cores$core_id)]
    if (length(dups) > 0) {
      warning(sid, ": several detected cores map to the same grid position: ",
              paste(unique(dups), collapse = ", "), " -- check the plot and add overrides.")
    }
    ## Bounding boxes (used for per-core cropping in 03a/03b/07) cover the same
    ## area as the cell assignment: radius x CORE_MAX_DIST_FACTOR.
    half <- cores$radius * CORE_MAX_DIST_FACTOR
    cores$xmin <- cores$x_center - half
    cores$xmax <- cores$x_center + half
    cores$ymin <- cores$y_center - half
    cores$ymax <- cores$y_center + half
    unmatched <- cores[is.na(cores$core_id), ]
    if (nrow(unmatched) > 0) {
      message(sid, ": ", nrow(unmatched), " detected cores have no tma_map entry (detected_core ",
              paste(unmatched$detected_core, collapse = ", "), ") -- cells there are dropped.")
    }
    utils::write.csv(cores, cores_file, row.names = FALSE)
    cores <- cores[!is.na(cores$core_id), ]
    meta$core_id[in_slide] <- assign_cells_to_cores(meta$x_um[in_slide], meta$y_um[in_slide], cores)
  }

  cores$slide_id <- sid
  all_cores[[sid]] <- cores

  ## ---- the plot you MUST look at: detected cores + labels vs the TMA map --
  d <- meta[in_slide, ]
  d$core_label <- ifelse(is.na(d$core_id), "unassigned", "core")
  if (nrow(d) > 3e5) d <- d[sample(nrow(d), 3e5), ]
  lab <- cores
  lab$text <- if ("detected_core" %in% colnames(lab)) {
    paste0(lab$core_id, "\n(", lab$core_row, lab$core_col, ", #", lab$detected_core, ")")
  } else lab$core_id
  p_dearray <- ggplot(d, aes(x_um, y_um)) +
    geom_point(aes(color = core_label), size = 0.05, stroke = 0) +
    geom_text(data = lab, aes(x_center, y_center, label = text), size = 2.5, fontface = "bold") +
    scale_color_manual(values = c(core = "grey60", unassigned = "red")) +
    scale_y_reverse() + coord_fixed() + theme_bw() +
    ggtitle(paste0(sid, ": core assignment -- compare every label with your TMA map"))
  save_fig(p_dearray, paste0("02_dearray_", sid), width = 12, height = 12)
}

cores_all <- dplyr::bind_rows(all_cores)
save_table(cores_all, "02_cores_all")

## -----------------------------------------------------------------------
## Attach TMA metadata to cells
## -----------------------------------------------------------------------
map_cols <- setdiff(colnames(tma_map), c("slide_id", "core_row", "core_col", "notes"))
meta_join <- meta[, c("slide_id", "core_id")] |>
  left_join(tma_map[, map_cols], by = "core_id")
for (col in setdiff(map_cols, "core_id")) xen[[col]] <- meta_join[[col]]
xen$core_id <- meta$core_id

message("Cells outside any core: ", sum(is.na(xen$core_id)),
        " (", round(100 * mean(is.na(xen$core_id)), 1), "%)")

## -----------------------------------------------------------------------
## Core-level QC. Flag, don't delete: whether a poor core is excluded is an
## explicit decision (set include = FALSE in tma_map.csv), not a side effect.
## -----------------------------------------------------------------------
m <- xen[[]]
core_qc <- m |>
  filter(!is.na(core_id)) |>
  group_by(slide_id, core_id, patient_id, organ, tissue_type) |>
  summarise(
    n_cells = n(),
    n_pass = sum(qc_pass),
    pct_pass = 100 * mean(qc_pass),
    median_counts = median(nCount_Xenium[qc_pass]),
    median_genes = median(nFeature_Xenium[qc_pass]),
    median_area = median(cell_area[qc_pass], na.rm = TRUE),
    control_frac = sum(nCount_controls) / sum(nCount_Xenium + nCount_controls),
    .groups = "drop"
  ) |>
  mutate(core_flag = case_when(
    n_pass < CORE_MIN_CELLS ~ "few_cells",
    median_counts < CORE_MIN_MEDIAN_COUNTS ~ "low_counts",
    TRUE ~ "ok"
  ))
print(as.data.frame(core_qc))
save_table(core_qc, "02_core_qc")

missing_cores <- setdiff(tma_map$core_id[tma_map$slide_id %in% slides$slide_id], core_qc$core_id)
if (length(missing_cores) > 0) {
  message("Cores in tma_map.csv with NO cells (lost during sectioning?): ",
          paste(missing_cores, collapse = ", "))
}

p_core_qc <- ggplot(core_qc, aes(reorder(core_id, median_counts), median_counts, fill = organ)) +
  geom_col() +
  geom_hline(yintercept = CORE_MIN_MEDIAN_COUNTS, linetype = 2) +
  facet_grid(~ slide_id, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  labs(x = "core", y = "median transcripts / cell (QC-pass cells)",
       title = "Core quality -- differences here confound every organ comparison in Aim 1/2")
save_fig(p_core_qc, "02_core_qc_median_counts", width = 16, height = 6)

## -----------------------------------------------------------------------
## Filter: QC-pass cells in cores that are in the map and included
## -----------------------------------------------------------------------
xen$core_flag <- core_qc$core_flag[match(xen$core_id, core_qc$core_id)]
keep <- colnames(xen)[!is.na(xen$core_id) & xen$qc_pass & xen$include %in% TRUE]
message("Keeping ", length(keep), " of ", ncol(xen), " cells.")
xen <- subset(xen, cells = keep)

saveRDS(xen, obj_path("cores"))
message("Saved -> ", obj_path("cores"))
