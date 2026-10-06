## =============================================================================
## 02a_make_pseudo_tma.R   (TESTING ONLY -- skip for real TMAs)
##
## Turns whole tissue sections (e.g. Swiss rolls) into a "virtual TMA" so the
## TMA pipeline can be tested before real TMAs exist: circular cores of
## PSEUDO_CORE_DIAMETER_UM are punched on a regular grid (pitch
## PSEUDO_CORE_PITCH_UM) wherever there is enough tissue -- exactly what a TMA
## needle does. Cells outside the punches are dropped in 02, as on a real TMA.
##
## Why: with one sample per slide, sample and slide are the same thing, so
## nothing can be compared within a slide and slide effects cannot be
## separated from biology. With virtual cores every slide holds many cores,
## and each core gets a LOCATION (inner / middle / outer ring of the roll,
## i.e. along the proximal-distal axis of a Swiss roll) that varies WITHIN a
## slide -- so location effects can be tested separately from slide effects,
## like organ/location on the real TMAs.
##
## Writes, for every slide in PSEUDO_TMA_SLIDES:
##   config/core_selections/<slide_id>/<core_id>.csv   circle outline (X, Y in um)
##   config/tma_map.csv                                  one row per virtual core
##                                                       (old file backed up first)
## then run 02_tma_dearray.R as usual: it finds core_selections/ and uses
## polygon mode. Edit the generated tma_map.csv freely (patient, organ, ...).
##
## INPUT : xenium_<seg>_raw.rds (01)
## OUTPUT: the files above + figures/02a_pseudo_tma_<slide>.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen <- readRDS(obj_path("raw"))
slides <- read_slides()
meta <- xen[[]]
pseudo_slides <- if (is.null(PSEUDO_TMA_SLIDES)) slides$slide_id else PSEUDO_TMA_SLIDES

r <- PSEUDO_CORE_DIAMETER_UM / 2
circle <- function(cx, cy, n = 64) {
  a <- seq(0, 2 * pi, length.out = n + 1)[-1]
  data.frame(X = cx + r * cos(a), Y = cy + r * sin(a))
}

new_rows <- list()
for (sid in pseudo_slides) {
  m <- meta[meta$slide_id == sid & meta$qc_pass, ]
  stopifnot("slide not found in the object" = nrow(m) > 0)

  ## candidate punch centres on a regular grid over the tissue
  gx <- seq(min(m$x_um) + r, max(m$x_um) - r, by = PSEUDO_CORE_PITCH_UM)
  gy <- seq(min(m$y_um) + r, max(m$y_um) - r, by = PSEUDO_CORE_PITCH_UM)
  cand <- expand.grid(col = seq_along(gx), row = seq_along(gy))
  cand$x <- gx[cand$col]; cand$y <- gy[cand$row]

  ## keep punches that land on tissue: enough QC-pass cells inside the circle
  nn <- RANN::nn2(cbind(m$x_um, m$y_um), cbind(cand$x, cand$y),
                  k = min(nrow(m), 5000), searchtype = "radius", radius = r)
  cand$n_cells <- rowSums(nn$nn.idx > 0)
  cores <- cand[cand$n_cells >= PSEUDO_MIN_CELLS, ]
  if (nrow(cores) == 0) stop(sid, ": no punch reaches PSEUDO_MIN_CELLS -- lower it or the diameter")

  ## grid labels like a real TMA: rows A, B, ... from the top, columns 1, 2, ...
  rows_used <- sort(unique(cores$row)); cols_used <- sort(unique(cores$col))
  cores$core_row <- LETTERS_EXT(match(cores$row, rows_used))
  cores$core_col <- match(cores$col, cols_used)
  cores$core_id <- sprintf("%s_%s%02d", sid, cores$core_row, cores$core_col)

  ## location within the section: radial distance from the roll's centre, split
  ## into thirds. Centre = middle of the tissue's bounding box (the median cell
  ## would be pulled outwards, as outer turns of a roll hold more cells), or
  ## your own centre from PSEUDO_CENTERS (read it off the 01 maps' um axes).
  own <- if (!is.null(PSEUDO_CENTERS)) PSEUDO_CENTERS[PSEUDO_CENTERS$slide_id == sid, ] else NULL
  if (!is.null(own) && nrow(own) == 1) {
    cx <- own$x; cy <- own$y
  } else {
    cx <- (min(m$x_um) + max(m$x_um)) / 2; cy <- (min(m$y_um) + max(m$y_um)) / 2
  }
  d <- sqrt((cores$x - cx)^2 + (cores$y - cy)^2)
  cores$location <- as.character(cut(d, stats::quantile(d, c(0, 1/3, 2/3, 1)),
                                     labels = c("inner", "middle", "outer"), include.lowest = TRUE))

  sel_dir <- file.path(CONFIG_DIR, "core_selections", sid)
  unlink(sel_dir, recursive = TRUE)
  dir.create(sel_dir, recursive = TRUE)
  for (k in seq_len(nrow(cores))) {
    utils::write.csv(circle(cores$x[k], cores$y[k]), file.path(sel_dir, paste0(cores$core_id[k], ".csv")),
                     row.names = FALSE)
  }
  ## stale automatic core definitions would otherwise shadow nothing, but keep things clean
  unlink(file.path(TABLES_DIR, paste0("02_cores_", sid, ".csv")))

  new_rows[[sid]] <- data.frame(
    slide_id = sid, core_row = cores$core_row, core_col = cores$core_col, core_id = cores$core_id,
    patient_id = PSEUDO_PATIENT_ID %||% sid, organ = PSEUDO_ORGAN, tissue_type = "normal",
    diagnosis = "", location = cores$location, include = TRUE,
    notes = sprintf("virtual core, %d QC-pass cells", cores$n_cells)
  )
  message(sid, ": ", nrow(cores), " virtual cores (", paste(names(table(cores$location)),
          table(cores$location), sep = "=", collapse = ", "), ")")

  ## figure: tissue + punches, labelled like the real dearray plot
  circ <- do.call(rbind, lapply(seq_len(nrow(cores)), function(k)
    cbind(circle(cores$x[k], cores$y[k]), core_id = cores$core_id[k], location = cores$location[k])))
  bg <- m[sample(nrow(m), min(nrow(m), 2e5)), ]
  p <- ggplot() +
    geom_point(data = bg, aes(x_um, y_um), size = 0.05, colour = "grey75", stroke = 0) +
    geom_polygon(data = circ, aes(X, Y, group = core_id, colour = location), fill = NA, linewidth = 0.4) +
    geom_text(data = cores, aes(x, y, label = paste0(core_row, core_col)), size = 2) +
    scale_colour_manual(values = c(inner = "#2a78d6", middle = "#1baf7a", outer = "#eb6834")) +
    scale_y_reverse() + coord_fixed() + theme_bw() +
    annotate("point", x = cx, y = cy, shape = 3, size = 4) +
    labs(x = "x (um)", y = "y (um)", caption = "+ = centre used for inner / middle / outer",
         title = paste0(sid, ": virtual TMA, ", nrow(cores), " cores of ", PSEUDO_CORE_DIAMETER_UM, " um"))
  save_fig(p, paste0("02a_pseudo_tma_", sid), width = 11, height = 11)
}

## ---- tma_map.csv: replace the rows of the virtual-TMA slides, keep the rest --
old <- if (file.exists(TMA_MAP_CSV)) read_config_csv(TMA_MAP_CSV, stringsAsFactors = FALSE) else NULL
if (!is.null(old)) {
  bak <- paste0(TMA_MAP_CSV, ".", format(Sys.time(), "%Y%m%d%H%M%S"), ".bak")
  file.copy(TMA_MAP_CSV, bak)
  message("Backed up old tma_map.csv -> ", bak)
  old <- old[!old$slide_id %in% pseudo_slides, , drop = FALSE]
}
## all columns as text so old (hand-edited) and new rows always combine
as_chr <- function(d) { if (!is.null(d)) d[] <- lapply(d, as.character); d }
tma_map_new <- dplyr::bind_rows(as_chr(old), as_chr(dplyr::bind_rows(new_rows)))
utils::write.csv(tma_map_new, TMA_MAP_CSV, row.names = FALSE)
message("Wrote ", nrow(tma_map_new), " rows -> ", TMA_MAP_CSV, "\nNow run 02_tma_dearray.R (polygon mode).")
