## =============================================================================
## 09_he_regions.R
##
## H&E as the histology anchor: pathology regions drawn on the H&E are
## transferred onto Xenium (and CellScape) cells, so Aim 2 can compare immune
## targets between LOCATIONS WITHIN a core (tumour vs stroma vs invasive
## margin vs lymphoid aggregate vs normal mucosa), not only between cores.
##
## Recommended set-up: H&E on the Xenium slide after the run (10x post-Xenium
## H&E protocol) -> same section as the transcripts, so a global affine (07,
## --mode he) is essentially exact. If the H&E is on another section, treat
## it like CellScape (per-core refinement in 07).
##
## Workflow:
##   1. Register: python 07_register_modalities.py --mode he ...
##        -> data/registration/transforms_he_<slide>.json (H&E px -> Xenium um)
##   2. Annotate regions on the H&E in QuPath (classes e.g. Tumor, Stroma,
##      Immune_aggregate, Normal_epithelium, Necrosis, Muscle). Export:
##      File > Export objects as GeoJSON (pixel coordinates, full resolution)
##        -> data/he/annotations/<slide_id>.geojson
##   3. Run this script.
##
## Overlapping annotations: the SMALLEST polygon wins (e.g. an immune
## aggregate drawn inside a stroma region). Holes in polygons are ignored.
##
## OUTPUT: xenium_<seg>_he.rds (adds `he_region`), tables/09_*.csv, figures/09_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

ANNOT_GEOJSON_DIR <- file.path(PROJECT_DIR, "data", "he", "annotations")

in_path <- c(obj_path("integrated"), obj_path("insights"), obj_path("annotated"))
xen <- readRDS(in_path[file.exists(in_path)][1])
cores_all <- utils::read.csv(file.path(TABLES_DIR, "02_cores_all.csv"), stringsAsFactors = FALSE)

read_qupath_geojson <- function(path) {
  gj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  feats <- if (!is.null(gj$features)) gj$features else gj
  out <- list()
  for (f in feats) {
    cls <- f$properties$classification$name %||% f$properties$name %||% "unclassified"
    g <- f$geometry
    rings <- switch(g$type,
                    Polygon = list(g$coordinates[[1]]),
                    MultiPolygon = lapply(g$coordinates, function(p) p[[1]]),
                    NULL)
    for (r in rings) {
      xy <- do.call(rbind, lapply(r, function(p) as.numeric(unlist(p)[1:2])))
      out[[length(out) + 1]] <- list(class = cls, xy = xy)
    }
  }
  out
}

apply_tf <- function(M, xy) {
  h <- cbind(xy, 1) %*% t(M)
  h[, 1:2, drop = FALSE]
}

polygon_area <- function(xy) abs(sum(xy[, 1] * c(xy[-1, 2], xy[1, 2]) - c(xy[-1, 1], xy[1, 1]) * xy[, 2])) / 2

xen$he_region <- NA_character_
region <- stats::setNames(rep(NA_character_, ncol(xen)), colnames(xen))

for (sid in unique(xen$slide_id)) {
  tf_path <- file.path(REG_DIR, paste0("transforms_he_", sid, ".json"))
  gj_path <- file.path(ANNOT_GEOJSON_DIR, paste0(sid, ".geojson"))
  if (!file.exists(tf_path) || !file.exists(gj_path)) {
    message(sid, ": missing ", if (!file.exists(tf_path)) tf_path else gj_path, " -- skipped")
    next
  }
  tf <- jsonlite::fromJSON(tf_path)
  G <- as.matrix(tf$global)
  polys <- read_qupath_geojson(gj_path)
  cores_s <- cores_all[cores_all$slide_id == sid & !is.na(cores_all$core_id), ]

  ## H&E px -> Xenium um (global), then the per-core correction of the core
  ## the polygon falls in, if 07 refined that core.
  for (k in seq_along(polys)) {
    xy <- apply_tf(G, polys[[k]]$xy)
    cid <- assign_cells_to_cores(mean(xy[, 1]), mean(xy[, 2]), cores_s, factor = 2)
    if (!is.na(cid) && !is.null(tf$cores[[cid]])) xy <- apply_tf(as.matrix(tf$cores[[cid]]), xy)
    polys[[k]]$xy_um <- xy
    polys[[k]]$area <- polygon_area(xy)
  }
  polys <- polys[order(vapply(polys, `[[`, numeric(1), "area"))]   # smallest first

  cells <- colnames(xen)[xen$slide_id == sid]
  x <- xen$x_um[match(cells, colnames(xen))]
  y <- xen$y_um[match(cells, colnames(xen))]
  reg_s <- rep(NA_character_, length(cells))
  for (p in polys) {
    inside <- is.na(reg_s) & sp::point.in.polygon(x, y, p$xy_um[, 1], p$xy_um[, 2]) > 0
    reg_s[inside] <- p$class
  }
  region[cells] <- reg_s
  message(sid, ": ", length(polys), " polygons; ", round(100 * mean(!is.na(reg_s)), 1),
          "% of cells inside an annotation")

  ## Sanity plot: polygons over cells for a few cores. Misplaced polygons =
  ## wrong transform direction/landmarks, not a biology result.
  show <- head(unique(xen$core_id[xen$slide_id == sid]), 6)
  ## each polygon is drawn whole, in the facet of the core its centroid falls in
  poly_df <- dplyr::bind_rows(lapply(seq_along(polys), function(k) {
    xy <- polys[[k]]$xy_um
    data.frame(id = k, class = polys[[k]]$class, x_um = xy[, 1], y_um = xy[, 2],
               core_id = assign_cells_to_cores(mean(xy[, 1]), mean(xy[, 2]), cores_s, factor = 2))
  }))
  d <- data.frame(core_id = xen$core_id[match(cells, colnames(xen))], x_um = x, y_um = y,
                  he_region = ifelse(is.na(reg_s), "none", reg_s))
  d <- d[d$core_id %in% show, ]
  poly_df <- poly_df[poly_df$core_id %in% show, ]
  p <- ggplot(d, aes(x_um, y_um)) +
    geom_point(aes(color = he_region), size = 0.2) +
    geom_polygon(data = poly_df, aes(group = id), fill = NA, color = "black", linewidth = 0.3) +
    facet_wrap(~ core_id, scales = "free") + scale_y_reverse() + theme_bw() + theme(aspect.ratio = 1) +
    ggtitle(paste0(sid, ": H&E annotations on Xenium cells"))
  save_fig(p, paste0("09_he_regions_check_", sid), width = 14, height = 10)
}
xen$he_region <- unname(region)

## -----------------------------------------------------------------------
## Region-level summaries for Aim 2: composition and immune targets by
## region within organ (unit = core x region).
## -----------------------------------------------------------------------
m <- xen[[]] |> filter(!is.na(he_region))
if (nrow(m) > 0) {
  comp <- m |> count(organ, core_id, he_region, cell_type) |>
    group_by(core_id, he_region) |> mutate(frac = n / sum(n)) |> ungroup()
  save_table(comp, "09_composition_by_he_region")
  p_comp <- comp |> group_by(organ, he_region, cell_type) |> summarise(frac = mean(frac), .groups = "drop") |>
    ggplot(aes(he_region, frac, fill = cell_type)) + geom_col() + facet_wrap(~ organ) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "mean fraction (over cores)", title = "Cell-type composition per H&E region, by organ")
  save_fig(p_comp, "09_composition_by_he_region", width = 16, height = 10)

  probes <- read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE)
  probes$gene <- match_case(as_seurat_features(probes$gene), rownames(xen[["Xenium"]]))
  targets <- intersect(probes$gene[probes$category == "immune_target"], rownames(xen[["Xenium"]]))
  if (length(targets) > 0) {
    cnt <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m), drop = FALSE]
    grp <- paste(m$core_id, m$he_region, sep = "||")
    ok <- names(which(table(grp) >= MIN_CELLS_PSEUDOBULK))
    pb <- pseudobulk(cnt[, grp %in% ok, drop = FALSE], grp[grp %in% ok])
    lcpm <- edgeR::cpm(as.matrix(pb), log = TRUE, prior.count = 1)[targets, , drop = FALSE]
    long <- as.data.frame(as.table(lcpm), stringsAsFactors = FALSE)
    colnames(long) <- c("gene", "group", "logCPM")
    long$core_id <- sub("\\|\\|.*$", "", long$group)
    long$he_region <- sub("^.*\\|\\|", "", long$group)
    long$organ <- m$organ[match(long$core_id, m$core_id)]
    save_table(long, "09_immune_targets_by_he_region")
    p_t <- ggplot(long, aes(he_region, logCPM)) + geom_boxplot(outlier.size = 0.4) +
      facet_grid(gene ~ organ, scales = "free_y") + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(size = 6)) +
      ggtitle("Immune targets by H&E region (pseudobulk per core x region)")
    save_fig(p_t, "09_immune_targets_by_he_region", width = 16, height = max(6, 1.2 * length(targets)))
  }
}

saveRDS(xen, obj_path("he"))
message("Saved -> ", obj_path("he"))
