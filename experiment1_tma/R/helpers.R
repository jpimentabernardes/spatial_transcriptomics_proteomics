## =============================================================================
## R/helpers.R
##
## Shared functions for all experiment1_tma R scripts. Sourced right after
## config/config.R. Nothing in here reads data on its own -- only functions.
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## -----------------------------------------------------------------------
## Output helpers (same pattern as the class scripts)
## -----------------------------------------------------------------------
save_fig <- function(plot_obj, name, width = 9, height = 7, dpi = 150) {
  path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
  message("Saved figure -> ", path)
  invisible(path)
}

save_table <- function(df, name) {
  path <- file.path(TABLES_DIR, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  message("Saved table  -> ", path)
  invisible(path)
}

obj_path <- function(stage, segmentation = SEGMENTATION) {
  file.path(PROCESSED_DIR, sprintf("xenium_%s_%s.rds", segmentation, stage))
}

## -----------------------------------------------------------------------
## Config readers
## -----------------------------------------------------------------------
## read.csv() for sheets you edit by hand. Excel and most editors save the last
## row without a trailing newline; R then warns "incomplete final line", which is
## harmless (the row is read) -- so that one warning is muffled here.
read_config_csv <- function(path, ...) {
  withCallingHandlers(
    utils::read.csv(path, ...),
    warning = function(w) {
      if (grepl("incomplete final line", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

read_slides <- function(path = SLIDES_CSV) {
  s <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("slides.csv needs slide_id and xenium_dir" =
              all(c("slide_id", "xenium_dir") %in% colnames(s)))
  bad <- s$slide_id[s$slide_id != make.names(s$slide_id)]
  if (length(bad) > 0) {
    stop("slide_id must be a syntactically valid R name (used as FOV name): ",
         paste(bad, collapse = ", "), ". Use letters, digits, _ and . only.")
  }
  for (col in c("row_flip", "col_flip", "cs_row_flip", "cs_col_flip")) {
    if (col %in% colnames(s)) s[[col]] <- as.logical(s[[col]]) %in% TRUE
  }
  s
}

read_tma_map <- function(path = TMA_MAP_CSV) {
  m <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("slide_id", "core_row", "core_col", "core_id", "patient_id", "organ")
  stopifnot("tma_map.csv is missing required columns" = all(need %in% colnames(m)))
  for (opt in c("tissue_type", "diagnosis", "location")) if (!opt %in% colnames(m)) m[[opt]] <- NA_character_
  if (!"include" %in% colnames(m)) m$include <- TRUE
  m$include <- as.logical(m$include) %in% TRUE
  m$core_row <- as.character(m$core_row)
  m$core_col <- as.character(m$core_col)
  dup <- m$core_id[duplicated(m$core_id)]
  if (length(dup) > 0) stop("core_id must be unique across slides: ", paste(unique(dup), collapse = ", "))
  m
}

## Seurat turns "_" into "-" in feature names. Apply the same to every gene /
## protein name read from config, or custom probes like "IL23_A" silently
## never match.
as_seurat_features <- function(x) gsub("_", "-", x, fixed = TRUE)

## Config gene names -> the spelling used in the object, ignoring upper / lower case,
## so the human-style names in config/ (EPCAM, PTPRC, CD3E) also find mouse genes
## (Epcam, Ptprc, Cd3e). Exact matches win; names not in `features` come back
## unchanged (and are then reported as missing by the caller).
match_case <- function(x, features) {
  hit <- match(toupper(x), toupper(features))
  ifelse(x %in% features | is.na(hit), x, features[hit])
}

split_markers <- function(x) {
  x <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ";")))
  as_seurat_features(x[nzchar(x)])
}

fov_key <- function(id) paste0("fov", tolower(gsub("[^[:alnum:]]", "", id)), "_")

## Returns list(level1 = named list, level2 = list(parent -> named list), mecr = named list)
read_marker_panel <- function(path) {
  p <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  to_list <- function(df) stats::setNames(lapply(df$markers, split_markers), df$cell_type)
  out <- list(level1 = to_list(p[p$level == 1, ]), level2 = list(), mecr = list())
  for (par in unique(p$parent[p$level == 2])) {
    out$level2[[par]] <- to_list(p[p$level == 2 & p$parent == par, ])
  }
  if ("use_for_mecr" %in% colnames(p)) {
    l1 <- p[p$level == 1, ]
    out$mecr <- stats::setNames(lapply(l1$use_for_mecr, split_markers), l1$cell_type)
    out$mecr <- out$mecr[lengths(out$mecr) > 0]
  }
  out
}

## Keep only markers present in the panel; report what is missing. For a
## targeted panel this is not cosmetic -- a cell type with 0-1 present
## markers cannot be annotated reliably, and that is itself an Aim 1 finding.
restrict_to_panel <- function(sig_list, panel_genes, label = "", min_reliable = 2) {
  sig_list <- lapply(sig_list, match_case, features = panel_genes)
  cov <- data.frame(
    cell_type = names(sig_list),
    n_markers = lengths(sig_list),
    n_present = vapply(sig_list, function(g) sum(g %in% panel_genes), integer(1)),
    missing = vapply(sig_list, function(g) paste(setdiff(g, panel_genes), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  weak <- cov$cell_type[cov$n_present < min_reliable]
  if (length(weak) > 0) {
    message(label, " signatures with <", min_reliable, " markers in panel (unreliable): ", paste(weak, collapse = ", "))
  }
  sig <- lapply(sig_list, intersect, panel_genes)
  list(signatures = sig[lengths(sig) > 0], coverage = cov)
}

## -----------------------------------------------------------------------
## Xenium feature types. 10x ships several control-feature classes next to
## the genes; they are the basis of every background / specificity estimate
## in this project, so keep them (as per-cell metadata), don't discard them.
## -----------------------------------------------------------------------
CONTROL_SHORT <- c(
  "Negative Control Probe"    = "ControlProbe",
  "Negative Control Codeword" = "ControlCodeword",
  "Unassigned Codeword"       = "BlankCodeword",
  "Blank Codeword"            = "BlankCodeword",
  "Genomic Control"           = "GenomicControl",
  "Deprecated Codeword"       = "DeprecatedCodeword"
)

control_short_name <- function(type) {
  ifelse(type %in% names(CONTROL_SHORT), CONTROL_SHORT[type], make.names(type))
}

## Reads a 10x-style cell_feature_matrix (directory or .h5). Always returns a
## named list by feature type, even if there is only gene expression.
read_feature_matrix <- function(seg_dir) {
  mtx_dir <- file.path(seg_dir, "cell_feature_matrix")
  h5 <- file.path(seg_dir, "cell_feature_matrix.h5")
  counts <- if (dir.exists(mtx_dir)) {
    Seurat::Read10X(mtx_dir)
  } else if (file.exists(h5)) {
    Seurat::Read10X_h5(h5)
  } else {
    stop("No cell_feature_matrix/ or cell_feature_matrix.h5 in ", seg_dir)
  }
  if (!is.list(counts)) counts <- list("Gene Expression" = counts)
  counts
}

read_cells_table <- function(seg_dir) {
  pq <- file.path(seg_dir, "cells.parquet")
  cells <- if (file.exists(pq)) {
    as.data.frame(arrow::read_parquet(pq))
  } else if (file.exists(file.path(seg_dir, "cells.csv.gz"))) {
    utils::read.csv(file.path(seg_dir, "cells.csv.gz"), stringsAsFactors = FALSE)
  } else {
    stop("No cells.parquet or cells.csv.gz in ", seg_dir)
  }
  cells$cell_id <- as.character(cells$cell_id)
  for (col in c("cell_area", "nucleus_area")) if (!col %in% colnames(cells)) cells[[col]] <- NA_real_
  if (!"segmentation_method" %in% colnames(cells)) cells$segmentation_method <- NA_character_
  cells
}

## Builds one Seurat object per slide with an identical structure for every
## segmentation method, so downstream code never branches on method:
##   - assay "Xenium" = genes only
##   - controls summed into metadata: nCount_ControlProbe, nCount_BlankCodeword, ...
##     plus nFeatTotal_<type> (number of control features on the panel), needed
##     to turn control counts into a per-probe background rate
##   - x_um / y_um centroids in metadata AND as a centroid FOV named slide_id
##   - cells named <slide_id>_<cell_id> (10x cell ids are NOT unique across slides)
## Only centroids are loaded (no polygons): polygons for multi-slide TMAs cost
## many GB and are what broke subset() in the class integration script.
build_xenium_object <- function(slide_id, counts_list, cells, segmentation) {
  gex <- counts_list[["Gene Expression"]]
  stopifnot("No 'Gene Expression' features found" = !is.null(gex))

  common <- intersect(colnames(gex), cells$cell_id)
  if (length(common) < ncol(gex)) {
    message(slide_id, ": ", ncol(gex) - length(common),
            " cells in the matrix have no centroid and are dropped.")
  }
  gex <- gex[, common, drop = FALSE]
  cells <- cells[match(common, cells$cell_id), , drop = FALSE]
  new_names <- paste0(slide_id, "_", common)
  colnames(gex) <- new_names

  meta <- data.frame(
    slide_id = slide_id,
    cell_id = common,
    x_um = cells$x_centroid,
    y_um = cells$y_centroid,
    cell_area = cells$cell_area,
    nucleus_area = cells$nucleus_area,
    seg_boundary_source = cells$segmentation_method,
    segmentation = segmentation,
    row.names = new_names,
    stringsAsFactors = FALSE
  )
  for (type in setdiff(names(counts_list), "Gene Expression")) {
    short <- control_short_name(type)
    m <- counts_list[[type]][, common, drop = FALSE]
    prev <- if (paste0("nCount_", short) %in% colnames(meta)) meta[[paste0("nCount_", short)]] else 0
    meta[[paste0("nCount_", short)]] <- prev + Matrix::colSums(m)
    prev_n <- if (paste0("nFeatTotal_", short) %in% colnames(meta)) meta[[paste0("nFeatTotal_", short)]] else 0
    meta[[paste0("nFeatTotal_", short)]] <- prev_n + nrow(m)
  }
  meta$nFeatTotal_Gene <- nrow(gex)

  obj <- CreateSeuratObject(counts = gex, assay = "Xenium", meta.data = meta)
  fov <- CreateFOV(
    coords = data.frame(x = meta$x_um, y = meta$y_um, cell = rownames(meta)),
    type = "centroids",
    assay = "Xenium",
    key = fov_key(slide_id)
  )
  obj[[slide_id]] <- fov
  obj
}

seg_dir_for <- function(slide_row, segmentation) {
  if (segmentation == "xenium") return(resolve_path(slide_row$xenium_dir))
  file.path(SEG_DIR, segmentation, slide_row$slide_id)
}

load_segmentation <- function(slide_row, segmentation = SEGMENTATION) {
  d <- seg_dir_for(slide_row, segmentation)
  stopifnot("Segmentation output folder not found" = dir.exists(d))
  message("Loading ", slide_row$slide_id, " [", segmentation, "] from ", d)
  build_xenium_object(
    slide_id = slide_row$slide_id,
    counts_list = read_feature_matrix(d),
    cells = read_cells_table(d),
    segmentation = segmentation
  )
}

## merge() keeps one FOV per slide (named by slide_id). JoinLayers() because
## v5 merge splits counts into counts.1, counts.2, ... per input object.
merge_slides <- function(objs) {
  if (length(objs) == 1) return(objs[[1]])
  obj <- merge(objs[[1]], y = objs[-1])
  obj <- JoinLayers(obj)
  ## control feature totals are per-slide constants; merge may leave NA for
  ## slides whose panel lacked a control class -- set those to 0
  for (col in grep("^(nCount|nFeatTotal)_", colnames(obj[[]]), value = TRUE)) {
    v <- obj[[col, drop = TRUE]]
    v[is.na(v)] <- 0
    obj[[col]] <- v
  }
  obj
}

## -----------------------------------------------------------------------
## TMA dearraying
## -----------------------------------------------------------------------

## Grid-based clustering: bins of `cell` um, bins with >= 2 cells are tissue,
## 8-connected tissue bins form one cluster. Same result as DBSCAN for
## well-separated TMA cores, with no compiled dependency (fallback when the
## dbscan package cannot be loaded). Returns a cluster id per point, 0 = noise.
grid_components <- function(xy, cell) {
  gx <- floor((xy[, 1] - min(xy[, 1])) / cell) + 2      # +2: one empty border bin
  gy <- floor((xy[, 2] - min(xy[, 2])) / cell) + 2
  nx <- max(gx) + 1; ny <- max(gy) + 1
  lin <- gx + (gy - 1) * nx
  occ <- matrix(tabulate(lin, nbins = nx * ny) >= 2, nx, ny)
  lab <- matrix(NA_real_, nx, ny)
  lab[occ] <- seq_len(sum(occ))
  shift <- function(m, dx, dy) {                        # m[i + dx, j + dy], NA outside
    out <- matrix(NA_real_, nx, ny)
    xs <- max(1, 1 - dx):min(nx, nx - dx); ys <- max(1, 1 - dy):min(ny, ny - dy)
    out[xs, ys] <- m[xs + dx, ys + dy]
    out
  }
  repeat {                                              # propagate the smallest label
    old <- lab
    for (d in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))) {
      nb <- shift(lab, d[1], d[2])
      lab[occ] <- pmin(lab[occ], nb[occ], na.rm = TRUE)
    }
    if (identical(old, lab)) break
  }
  cl <- lab[lin]
  out <- integer(length(cl))                           # 0 = noise
  ok <- !is.na(cl)
  out[ok] <- as.integer(factor(cl[ok]))
  out
}

## DBSCAN on cell centroids -> core seeds -> merge fragments -> core table.
## Returns data.frame(core_idx, x_center, y_center, radius, n_cells_dbscan, detected_core)
detect_cores <- function(x, y,
                         eps = DBSCAN_EPS_UM, min_pts = DBSCAN_MIN_PTS,
                         min_cells = MIN_CELLS_PER_CORE, merge_dist = CORE_MERGE_DIST_UM) {
  xy <- cbind(x, y)
  if (is.na(eps)) {
    ## Adaptive: 3 x the median distance to the min_pts-th neighbour. Works for
    ## dense lymphoid cores and sparse lung/adipose alike, while staying far
    ## below the >= 200 um gaps between cores.
    eps <- 3 * stats::median(RANN::nn2(xy, k = min_pts + 1)$nn.dists[, min_pts + 1])
    message("Clustering radius eps (adaptive): ", round(eps, 1), " um")
  }
  cl <- if (requireNamespace("dbscan", quietly = TRUE)) {
    tryCatch(dbscan::dbscan(xy, eps = eps, minPts = min_pts)$cluster, error = function(e) NULL)
  }
  if (is.null(cl)) {
    message("dbscan package not usable -- using grid-based clustering instead")
    cl <- grid_components(xy, eps)
  }
  cl[cl == 0] <- NA
  if (all(is.na(cl))) stop("DBSCAN found no clusters -- increase DBSCAN_EPS_UM.")

  cent <- data.frame(cl = cl, x = x, y = y) |>
    filter(!is.na(cl)) |>
    group_by(cl) |>
    summarise(x = median(x), y = median(y), n = n(), .groups = "drop")

  ## Seeds = big clusters. Seeds closer than merge_dist are one torn core.
  seeds <- cent[cent$n >= min_cells, ]
  if (nrow(seeds) == 0) stop("No cluster reaches MIN_CELLS_PER_CORE cells.")
  seed_group <- if (nrow(seeds) > 1) {
    stats::cutree(stats::hclust(stats::dist(seeds[, c("x", "y")]), method = "single"), h = merge_dist)
  } else 1L
  lut <- stats::setNames(seed_group, seeds$cl)

  ## Small fragments attach to the nearest seed if within merge_dist.
  small <- cent[cent$n < min_cells, ]
  if (nrow(small) > 0) {
    nn <- RANN::nn2(seeds[, c("x", "y")], small[, c("x", "y")], k = 1)
    ok <- nn$nn.dists[, 1] <= merge_dist
    lut <- c(lut, stats::setNames(seed_group[nn$nn.idx[ok, 1]], small$cl[ok]))
  }
  core_idx <- unname(lut[as.character(cl)])

  cores <- data.frame(core_idx = core_idx, x = x, y = y) |>
    filter(!is.na(core_idx)) |>
    group_by(core_idx) |>
    summarise(
      x_center = median(x), y_center = median(y),
      radius = stats::quantile(sqrt((x - median(x))^2 + (y - median(y))^2), 0.95),
      n_cells_dbscan = n(),
      .groups = "drop"
    ) |>
    arrange(y_center, x_center) |>
    mutate(detected_core = seq_len(n())) |>
    as.data.frame()
  cores
}

## Nearest core centre, but only within CORE_MAX_DIST_FACTOR x radius. Works
## for any segmentation's centroids once cores are defined.
assign_cells_to_cores <- function(x, y, cores, factor = CORE_MAX_DIST_FACTOR) {
  nn <- RANN::nn2(cores[, c("x_center", "y_center")], cbind(x, y), k = 1)
  idx <- nn$nn.idx[, 1]
  ok <- nn$nn.dists[, 1] <= cores$radius[idx] * factor
  out <- cores$core_id[idx]
  out[!ok] <- NA
  out
}

## Rows = clusters of y centres, columns = clusters of x centres. Row A is at
## the TOP of the image (smallest y) unless row_flip -- set *_flip in
## slides.csv if the TMA was mounted rotated. ALWAYS check the labelled plot.
assign_grid <- function(cores, row_levels, col_levels, row_flip = FALSE, col_flip = FALSE) {
  cluster_axis <- function(v, k) {
    k <- min(k, length(unique(round(v))))
    if (k <= 1) return(rep(1L, length(v)))
    g <- stats::cutree(stats::hclust(stats::dist(v), method = "average"), k = k)
    rank_of_group <- rank(tapply(v, g, mean))
    as.integer(rank_of_group[as.character(g)])
  }
  r <- cluster_axis(cores$y_center, length(row_levels))
  c <- cluster_axis(cores$x_center, length(col_levels))
  if (max(r) < length(row_levels) || max(c) < length(col_levels)) {
    warning("Fewer detected rows/columns than in tma_map.csv (missing cores?). ",
            "Row/column labels may be shifted -- verify the dearray plot and use ",
            "config/core_overrides.csv if needed.")
  }
  if (row_flip) r <- max(r) + 1L - r
  if (col_flip) c <- max(c) + 1L - c
  cores$core_row <- row_levels[r]
  cores$core_col <- col_levels[c]
  cores
}

## TMA row names: 1 -> A, 26 -> Z, 27 -> AA, ...
LETTERS_EXT <- function(i) {
  vapply(i, function(k) {
    s <- ""
    while (k > 0) { s <- paste0(LETTERS[(k - 1) %% 26 + 1], s); k <- (k - 1) %/% 26 }
    s
  }, character(1))
}

sort_levels <- function(v) {
  v <- unique(as.character(v))
  num <- suppressWarnings(as.numeric(v))
  if (!anyNA(num)) v[order(num)] else sort(v)
}

## Optional polygon mode: one CSV per core in config/core_selections/<slide_id>/,
## file name = core_id, with X/Y columns in microns (e.g. Xenium Explorer
## "Download selection coordinates"; lines starting with # are skipped).
read_core_selections <- function(sel_dir) {
  files <- list.files(sel_dir, pattern = "\\.csv$", full.names = TRUE)
  lapply(stats::setNames(files, tools::file_path_sans_ext(basename(files))), function(f) {
    d <- read_config_csv(f, comment.char = "#", check.names = FALSE)
    xcol <- grep("^x$", colnames(d), ignore.case = TRUE, value = TRUE)[1]
    ycol <- grep("^y$", colnames(d), ignore.case = TRUE, value = TRUE)[1]
    stopifnot("Selection CSV needs X and Y columns" = !is.na(xcol) && !is.na(ycol))
    data.frame(x = d[[xcol]], y = d[[ycol]])
  })
}

assign_cells_to_polygons <- function(x, y, polygons) {
  out <- rep(NA_character_, length(x))
  for (nm in names(polygons)) {
    p <- polygons[[nm]]
    inside <- sp::point.in.polygon(x, y, p$x, p$y) > 0
    out[inside & is.na(out)] <- nm
  }
  out
}

## -----------------------------------------------------------------------
## Scoring and annotation
## -----------------------------------------------------------------------

## Adds one score column per signature named score_<signature>.
score_signatures <- function(obj, sigs, use_ucell = USE_UCELL, assay = "Xenium") {
  sigs <- sigs[lengths(sigs) > 0]
  DefaultAssay(obj) <- assay
  if (use_ucell && requireNamespace("UCell", quietly = TRUE)) {
    ## serial (ncores = 1): no BiocParallel workers; if UCell still fails (version /
    ## object-specific problems), say why and fall back to AddModuleScore below
    res <- tryCatch(
      UCell::AddModuleScore_UCell(obj, features = sigs, name = "_UCell", assay = assay, ncores = 1),
      error = function(e) {
        message("UCell failed (", conditionMessage(e), ") -- falling back to Seurat::AddModuleScore.")
        NULL
      })
    if (!is.null(res)) {
      for (s in names(sigs)) {
        res[[paste0("score_", s)]] <- res[[paste0(s, "_UCell"), drop = TRUE]]
        res[[paste0(s, "_UCell")]] <- NULL
      }
      return(res)
    }
  } else if (use_ucell) {
    message("UCell not installed -- falling back to AddModuleScore.")
  }
  ## control genes are sampled per expression bin, so ctrl must stay well
  ## below the bin size (a few hundred genes / nbin) on a targeted panel
  n_genes <- nrow(obj[[assay]])
  nbin <- 12
  obj <- AddModuleScore(obj, features = sigs, name = "tmpscore", assay = assay,
                        ctrl = max(2, min(50, floor(n_genes / (2 * nbin)))), nbin = nbin)
  for (i in seq_along(sigs)) {
    obj[[paste0("score_", names(sigs)[i])]] <- obj[[paste0("tmpscore", i), drop = TRUE]]
    obj[[paste0("tmpscore", i)]] <- NULL
  }
  obj
}

## Dot plot of marker genes per group (colour = scaled mean expression, size = % of
## cells with counts), computed straight from the assay's data layer. Replaces
## Seurat::DotPlot, which goes through FetchData() and fails ("non-numeric argument
## to mathematical function") when a metadata column is returned instead of a gene.
marker_dotplot <- function(obj, features, group.by, assay = "Xenium", layer = "data", log_data = TRUE) {
  features <- intersect(unique(features), rownames(obj[[assay]]))
  if (length(features) == 0) stop("marker_dotplot: none of the features are in assay ", assay)
  m <- GetAssayData(obj, assay = assay, layer = layer)[features, , drop = FALSE]
  if (!methods::is(m, "CsparseMatrix")) m <- methods::as(as.matrix(m), "CsparseMatrix")
  grp <- droplevels(as.factor(obj[[group.by, drop = TRUE]]))
  keep <- !is.na(grp)
  m <- m[, keep, drop = FALSE]; grp <- grp[keep]
  G <- Matrix::sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                            dims = c(length(grp), nlevels(grp)), dimnames = list(NULL, levels(grp)))
  n <- Matrix::colSums(G)
  lin <- m
  if (log_data) lin@x <- expm1(lin@x)                       # log-normalised -> linear, like DotPlot
  avg <- as.matrix(lin %*% G) / rep(n, each = nrow(m))
  pct <- as.matrix((m > 0) %*% G) / rep(n, each = nrow(m)) * 100
  sc <- t(scale(t(avg)))
  sc[is.na(sc)] <- 0
  sc <- pmin(pmax(sc, -2.5), 2.5)
  d <- data.frame(feature = factor(rep(features, times = ncol(avg)), levels = features),
                  group = factor(rep(colnames(avg), each = nrow(avg)), levels = levels(grp)),
                  avg_scaled = as.vector(sc), pct = as.vector(pct))
  ggplot(d, aes(feature, group)) +
    geom_point(aes(size = pct, colour = avg_scaled)) +
    scale_colour_gradient(low = "lightgrey", high = "blue", name = "Average\nexpression\n(scaled)") +
    scale_size(range = c(0, 6), limits = c(0, 100), name = "% expressed") +
    labs(x = "Features", y = group.by) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## Semi-manual annotation. Writes an editable CSV (cluster, suggested label,
## top markers, label, reviewed). You edit `label`, set reviewed = TRUE, and
## re-run: the file is then used as-is. This keeps the manual decision in a
## versionable file instead of a hard-coded vector in the script.
resolve_cluster_map <- function(clusters, score_df, map_path, top_markers = NULL) {
  clusters <- as.character(clusters)
  mean_scores <- stats::aggregate(score_df, by = list(cluster = clusters), FUN = mean)
  sig_names <- sub("^score_", "", colnames(score_df))
  m <- as.matrix(mean_scores[, -1, drop = FALSE])
  best <- apply(m, 1, which.max)
  srt <- t(apply(m, 1, sort, decreasing = TRUE))
  margin <- if (ncol(m) > 1) srt[, 1] - srt[, 2] else srt[, 1]
  template <- data.frame(
    cluster = mean_scores$cluster,
    n_cells = as.integer(table(clusters)[mean_scores$cluster]),
    suggested = sig_names[best],
    margin = round(margin, 3),
    top_markers = if (is.null(top_markers)) "" else unname(top_markers[mean_scores$cluster]),
    label = sig_names[best],
    reviewed = FALSE,
    stringsAsFactors = FALSE
  )
  template$top_markers[is.na(template$top_markers)] <- ""

  if (file.exists(map_path)) {
    existing <- read_config_csv(map_path, stringsAsFactors = FALSE, colClasses = c(cluster = "character"))
    if (setequal(existing$cluster, template$cluster)) {
      if (!any(as.logical(existing$reviewed) %in% TRUE)) {
        message("Using UNREVIEWED labels from ", map_path, " -- edit `label`, set reviewed=TRUE.")
      }
      return(stats::setNames(existing$label, existing$cluster))
    }
    bak <- paste0(map_path, ".", format(Sys.time(), "%Y%m%d%H%M%S"), ".bak")
    file.copy(map_path, bak)
    warning("Clusters changed since ", map_path, " was written (backup: ", bak, "). Rewriting template.")
  }
  utils::write.csv(template, map_path, row.names = FALSE)
  message("Wrote annotation template -> ", map_path,
          "\n  Low `margin` = ambiguous cluster; check it on the dotplot before accepting.")
  stats::setNames(template$label, template$cluster)
}

## Top-N positive markers per cluster as "G1;G2;G3" strings, on downsampled
## cells (FindAllMarkers on 1M cells is not a good use of a day).
top_markers_per_cluster <- function(obj, ident, n = 5, max_cells = 500) {
  Idents(obj) <- ident
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                       max.cells.per.ident = max_cells, verbose = FALSE)
  if (nrow(mk) == 0) return(NULL)
  mk |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = n) |>
    summarise(genes = paste(gene, collapse = ";"), .groups = "drop") |>
    (\(d) stats::setNames(d$genes, as.character(d$cluster)))()
}

## -----------------------------------------------------------------------
## Quantitative helpers
## -----------------------------------------------------------------------

## Sum counts by group: returns features x groups (sparse).
pseudobulk <- function(counts, group) {
  g <- factor(group)
  counts %*% Matrix::t(Matrix::fac2sparse(g))
}

## Mutually exclusive co-expression rate (Hartman & Satija 2024 idea): for
## marker pairs from DIFFERENT lineages, fraction of cells expressing both
## among cells expressing either. Lower = cleaner segmentation (less
## transcript spill-over between neighbouring cells).
compute_mecr <- function(counts, lineage_markers) {
  lineage_markers <- lapply(lineage_markers, function(g) intersect(match_case(g, rownames(counts)), rownames(counts)))
  lineage_markers <- lineage_markers[lengths(lineage_markers) > 0]
  genes <- unique(unlist(lineage_markers))
  if (length(genes) < 2) return(NA_real_)
  B <- counts[genes, , drop = FALSE] > 0
  lin_of <- stats::setNames(rep(names(lineage_markers), lengths(lineage_markers)), unlist(lineage_markers))
  pairs <- utils::combn(genes, 2)
  pairs <- pairs[, lin_of[pairs[1, ]] != lin_of[pairs[2, ]], drop = FALSE]
  if (ncol(pairs) == 0) return(NA_real_)
  rates <- apply(pairs, 2, function(p) {
    a <- B[p[1], ]; b <- B[p[2], ]
    either <- sum(a | b)
    if (either == 0) NA_real_ else sum(a & b) / either
  })
  mean(rates, na.rm = TRUE)
}

## Tissue-specificity index tau (Yanai 2005): 0 = uniform, 1 = one organ only.
tau_index <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || max(x) <= 0) return(NA_real_)
  sum(1 - x / max(x)) / (length(x) - 1)
}

## -----------------------------------------------------------------------
## Plotting
## -----------------------------------------------------------------------

## ggplot on metadata centroids instead of ImageDimPlot(): fast on millions of
## cells, never depends on @images surviving a subset, and can facet by core.
## axes = TRUE draws micron axes + grid, so you can read off coordinates (e.g.
## for INSPECT_CENTERS in config.R); theme_void otherwise.
plot_cores_spatial <- function(meta, color_by, cores = NULL, facet = TRUE,
                               max_cells = 3e5, size = 0.15, title = NULL, axes = FALSE) {
  d <- meta
  if (!is.null(cores)) d <- d[d$core_id %in% cores, , drop = FALSE]
  ## sort() keeps the caller's row order, so rows placed last are drawn on top
  if (nrow(d) > max_cells) d <- d[sort(sample(nrow(d), max_cells)), , drop = FALSE]
  p <- ggplot(d, aes(x_um, y_um, color = .data[[color_by]])) +
    geom_point(size = size, stroke = 0) +
    ggtitle(title %||% color_by)
  p <- if (axes) {
    p + scale_x_continuous(breaks = scales::breaks_pretty(10)) +
      scale_y_reverse(breaks = scales::breaks_pretty(10)) +   # image convention: y grows downwards
      theme_minimal(base_size = 9) +
      theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank()) +
      labs(x = "x (um)", y = "y (um)")
  } else {
    p + scale_y_reverse() + theme_void()
  }
  ## coord_fixed() cannot be combined with free facet scales; cores are ~round,
  ## so a square panel keeps them undistorted enough.
  do_facet <- facet && "core_id" %in% colnames(d)
  p <- p + if (do_facet) theme(aspect.ratio = 1) else coord_fixed()
  if (is.numeric(d[[color_by]])) {
    p <- p + scale_color_viridis_c()
  } else {
    p <- p + guides(color = guide_legend(override.aes = list(size = 3)))
  }
  if (do_facet) {
    p <- p + facet_wrap(~ core_id, scales = "free") + theme(strip.text = element_text(size = 7))
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## -----------------------------------------------------------------------
## Zoomed QC inspection: cells (pass / fail) on top of the morphology image
## -----------------------------------------------------------------------

## First failing QC reason per cell ("pass" if none) -- for colouring.
qc_reason <- function(meta) {
  r <- rep("pass", nrow(meta))
  r[meta$qc_high_control] <- "high control"
  r[meta$qc_large] <- "too large"
  r[meta$qc_small] <- "too small"
  r[meta$qc_low_counts] <- "low counts"
  factor(r, levels = c("pass", "low counts", "too small", "too large", "high control"))
}
QC_REASON_COLORS <- c(`pass` = "#9ec5f4", `low counts` = "#eda100", `too small` = "#e34948",
                      `too large` = "#eb6834", `high control` = "#e87ba4")

## Windows (square, window_um wide) that contain BOTH passing and failing cells:
## the slide is cut into a grid and the windows with the most balanced mix
## (largest min(n_pass, n_fail)) are returned, n_per_slide per slide.
pick_inspection_windows <- function(meta, window_um = 150, n_per_slide = 3, min_cells = 30) {
  d <- meta[!is.na(meta$x_um), c("slide_id", "x_um", "y_um", "qc_pass")]
  d$gx <- floor(d$x_um / window_um)
  d$gy <- floor(d$y_um / window_um)
  w <- d |>
    group_by(slide_id, gx, gy) |>
    summarise(n = n(), n_pass = sum(qc_pass), n_fail = sum(!qc_pass), .groups = "drop") |>
    filter(n >= min_cells) |>
    mutate(score = pmin(n_pass, n_fail)) |>
    group_by(slide_id) |>
    slice_max(score, n = n_per_slide, with_ties = FALSE) |>
    ungroup()
  data.frame(slide_id = w$slide_id,
             x0 = w$gx * window_um, y0 = w$gy * window_um,
             x1 = (w$gx + 1) * window_um, y1 = (w$gy + 1) * window_um,
             n_cells = w$n, n_pass = w$n_pass, n_fail = w$n_fail)
}

## Crop the morphology image (and cell / nucleus outlines) for one window via
## py/crop_region.py. Python does the reading because Xenium OME-TIFFs are tiled,
## pyramidal and compressed in ways R's TIFF readers do not handle reliably.
crop_morphology <- function(xenium_dir, x0, y0, x1, y1, tag, image = NULL,
                            python = PYTHON_BIN, out_dir = file.path(PROCESSED_DIR, "inspection")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(out_dir, tag)
  args <- c(file.path(CODE_DIR, "py", "crop_region.py"), "--xenium-dir", shQuote(xenium_dir),
            "--x0", x0, "--y0", y0, "--x1", x1, "--y1", y1, "--out-prefix", shQuote(prefix))
  if (!is.null(image)) args <- c(args, "--image", shQuote(image))
  out <- system2(python, args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("crop_region.py failed (set PYTHON_BIN in config.R to a Python with tifffile + zarr):\n",
         paste(utils::tail(out, 10), collapse = "\n"))
  }
  read_b <- function(kind) {
    f <- paste0(prefix, "_", kind, "_boundaries.csv.gz")
    if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  }
  list(image = tiff::readTIFF(paste0(prefix, "_image.tif"), as.is = TRUE),
       meta = jsonlite::fromJSON(paste0(prefix, "_meta.json")),
       cells = read_b("cell"), nuclei = read_b("nucleus"))
}

## Three panels: where the window is on the slide | morphology alone |
## morphology + cell outlines (or centroids) coloured by QC result.
## y is plotted as -y so the image keeps its orientation (y grows downwards).
plot_inspection <- function(meta_slide, crop, win, title = "", point_size = 1.2) {
  m <- crop$meta
  img <- crop$image
  lim <- stats::quantile(img, c(0.01, 0.997), na.rm = TRUE)
  g <- pmin(pmax((img - lim[1]) / max(lim[2] - lim[1], 1), 0), 1)^0.6      # contrast stretch
  ras <- grDevices::as.raster(g)
  bg <- annotation_raster(ras, xmin = m$x0_um, xmax = m$x1_um, ymin = -m$y1_um, ymax = -m$y0_um)
  frame <- list(coord_fixed(xlim = c(win$x0, win$x1), ylim = c(-win$y1, -win$y0), expand = FALSE),
                theme_void(), theme(plot.title = element_text(size = 9)))

  cells <- meta_slide[meta_slide$x_um >= win$x0 - 20 & meta_slide$x_um <= win$x1 + 20 &
                        meta_slide$y_um >= win$y0 - 20 & meta_slide$y_um <= win$y1 + 20, ]
  cells$qc <- qc_reason(cells)

  p_img <- ggplot() + bg + frame + ggtitle("morphology (DAPI)")

  p_cells <- ggplot() + bg
  if (!is.null(crop$nuclei) && nrow(crop$nuclei) > 0) {
    p_cells <- p_cells + geom_polygon(data = crop$nuclei, aes(vertex_x, -vertex_y, group = cell_id),
                                      fill = NA, colour = "white", linewidth = 0.15, alpha = 0.6)
  }
  if (!is.null(crop$cells) && nrow(crop$cells) > 0) {
    b <- crop$cells
    b$qc <- cells$qc[match(b$cell_id, cells$cell_id)]
    b <- b[!is.na(b$qc), ]
    p_cells <- p_cells + geom_polygon(data = b, aes(vertex_x, -vertex_y, group = cell_id, colour = qc),
                                      fill = NA, linewidth = 0.35)
  }
  p_cells <- p_cells +
    geom_point(data = cells, aes(x_um, -y_um, colour = qc), size = point_size) +
    scale_colour_manual(values = QC_REASON_COLORS, drop = FALSE, name = "cell QC") +
    frame + ggtitle({
      inw <- cells$x_um >= win$x0 & cells$x_um <= win$x1 & cells$y_um >= win$y0 & cells$y_um <= win$y1
      sprintf("cells in window: %d pass, %d fail (outlines = segmentation)",
              sum(cells$qc[inw] == "pass"), sum(cells$qc[inw] != "pass"))
    })

  ov <- meta_slide[sample(nrow(meta_slide), min(nrow(meta_slide), 1e5)), ]
  p_where <- ggplot(ov, aes(x_um, -y_um)) +
    geom_point(size = 0.05, colour = "grey60", stroke = 0) +
    annotate("rect", xmin = win$x0, xmax = win$x1, ymin = -win$y1, ymax = -win$y0,
             fill = NA, colour = "#e34948", linewidth = 0.8) +
    coord_fixed() + theme_void() + ggtitle("window position")

  (p_where | p_img | p_cells) + plot_layout(widths = c(1, 1.4, 1.4)) +
    plot_annotation(title = title)
}
