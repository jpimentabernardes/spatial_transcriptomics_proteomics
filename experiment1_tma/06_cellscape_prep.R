## =============================================================================
## 06_cellscape_prep.R
##
## Spatial proteomics (CellScape, consecutive section) -> QC'd, dearrayed,
## annotated Seurat object with the SAME core ids and the SAME lineage names
## as the Xenium object. Those two shared vocabularies are what make the
## cross-modal steps (07, 08) possible:
##   * core ids   -> matched cores give an automatic initial alignment (07)
##   * lineages   -> label-aware registration refinement (07/07b) and
##                   cross-modal agreement checks (08)
##
## INPUT : per-cell CellScape table (CellScape export, or 06a output), path in
##         slides.csv `cellscape_cells`. Column names -> config.R CS_* settings.
## OUTPUT: data/processed/cellscape_annotated.rds
##         data/registration/cellscape_cells_<xenium slide_id>.csv.gz (for 07)
##         tables/06_cs_cores_all.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

slides <- read_slides()
slides <- slides[!is.na(slides$cellscape_cells) & nzchar(slides$cellscape_cells), ]
stopifnot("No slide has `cellscape_cells` set in slides.csv" = nrow(slides) > 0)
tma_map <- read_tma_map()
markers_p <- read_marker_panel(MARKERS_PROT_CSV)
CS_DAPI_MARKER <- as_seurat_features(CS_DAPI_MARKER)
CS_EXCLUDE_MARKERS <- as_seurat_features(CS_EXCLUDE_MARKERS)

read_cs_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path))
  } else utils::read.csv(path, check.names = FALSE)
}

## -----------------------------------------------------------------------
## 1. Load each CellScape slide -> Seurat (assay PROT)
##    counts layer = raw mean intensities; data layer = arcsinh transform.
## -----------------------------------------------------------------------
per_marker_cofactor <- function(v) {
  pos <- v[v > 0]
  if (length(pos) == 0) return(1)
  max(stats::quantile(pos, 0.05), 1e-3)
}

objs <- list()
for (i in seq_len(nrow(slides))) {
  cs_id <- slides$cellscape_slide_id[i]
  tab <- read_cs_table(resolve_path(slides$cellscape_cells[i]))
  need <- c(CS_ID_COL, CS_X_COL, CS_Y_COL)
  if (!all(need %in% colnames(tab))) {
    stop("CellScape table for ", cs_id, " lacks ", paste(setdiff(need, colnames(tab)), collapse = ", "),
         ". Columns are: ", paste(head(colnames(tab), 30), collapse = ", "), " ... -> set CS_* in config.R")
  }
  mcols <- grep(CS_MARKER_REGEX, colnames(tab), value = TRUE)
  stopifnot("No marker columns matched CS_MARKER_REGEX" = length(mcols) > 0)
  mnames <- as_seurat_features(sub(CS_MARKER_STRIP, "", mcols))
  scale_xy <- if (CS_COORD_UNITS == "px") as.numeric(slides$cellscape_pixel_size_um[i]) else 1
  stopifnot("cellscape_pixel_size_um missing in slides.csv" = !is.na(scale_xy))

  cells <- paste0(cs_id, "_", tab[[CS_ID_COL]])
  raw <- t(as.matrix(tab[, mcols, drop = FALSE]))
  raw[is.na(raw)] <- 0
  dimnames(raw) <- list(mnames, cells)

  cof <- if (is.na(CS_ARCSINH_COFACTOR)) apply(raw, 1, per_marker_cofactor) else rep(CS_ARCSINH_COFACTOR, nrow(raw))
  trans <- asinh(sweep(pmax(raw, 0), 1, cof, "/"))

  meta <- data.frame(
    cs_slide_id = cs_id,
    slide_id = slides$slide_id[i],          # the Xenium slide this section pairs with
    x_um = tab[[CS_X_COL]] * scale_xy,
    y_um = tab[[CS_Y_COL]] * scale_xy,
    cell_area = if (CS_AREA_COL %in% colnames(tab)) tab[[CS_AREA_COL]] * scale_xy^2 else NA_real_,
    row.names = cells
  )
  o <- CreateSeuratObject(counts = Matrix::Matrix(raw, sparse = TRUE), assay = "PROT", meta.data = meta)
  o <- SetAssayData(o, assay = "PROT", layer = "data", new.data = Matrix::Matrix(trans, sparse = TRUE))
  o[[cs_id]] <- CreateFOV(data.frame(x = meta$x_um, y = meta$y_um, cell = cells),
                          type = "centroids", assay = "PROT", key = fov_key(cs_id))
  objs[[cs_id]] <- o
  message(cs_id, ": ", ncol(o), " cells, ", nrow(o), " markers (", paste(mnames, collapse = ", "), ")")
}
cs <- merge_slides(objs)
rm(objs); invisible(gc())

## -----------------------------------------------------------------------
## 2. Cell QC
##    - no/low nucleus signal (segmentation artefact)
##    - area outliers (merged cells)
##    - "everything bright" cells: autofluorescent debris, RBCs, edge
##      artefacts -- they are positive for every marker and would otherwise
##      form a fake "multi-positive" cluster
## -----------------------------------------------------------------------
dat <- GetAssayData(cs, assay = "PROT", layer = "data")
m <- cs[[]]
if (CS_DAPI_MARKER %in% rownames(dat)) {
  ld <- log1p(as.numeric(dat[CS_DAPI_MARKER, ]))
  m$qc_low_dapi <- as.logical(ave(ld, m$cs_slide_id, FUN = function(v) v < stats::median(v) - 3 * stats::mad(v)))
} else {
  message("No '", CS_DAPI_MARKER, "' marker -- skipping DAPI QC.")
  m$qc_low_dapi <- FALSE
}
la <- log(m$cell_area)
m$qc_area <- if (all(is.na(la))) FALSE else
  !is.na(la) & abs(la - stats::median(la, na.rm = TRUE)) > 4 * stats::mad(la, na.rm = TRUE)
use_m <- setdiff(rownames(dat), CS_EXCLUDE_MARKERS)
p99 <- apply(dat[use_m, , drop = FALSE], 1, stats::quantile, 0.99)
m$frac_markers_high <- Matrix::colMeans(dat[use_m, , drop = FALSE] > p99)
m$qc_allbright <- m$frac_markers_high >= 0.5
m$qc_pass <- !(m$qc_low_dapi | m$qc_area | m$qc_allbright)
for (col in c("qc_low_dapi", "qc_area", "frac_markers_high", "qc_allbright", "qc_pass")) cs[[col]] <- m[[col]]
print(m |> group_by(cs_slide_id) |>
        summarise(n = n(), pct_low_dapi = 100 * mean(qc_low_dapi), pct_area = 100 * mean(qc_area),
                  pct_allbright = 100 * mean(qc_allbright), pct_pass = 100 * mean(qc_pass)))

## -----------------------------------------------------------------------
## 3. Dearray with the same TMA map (the consecutive section carries the
##    same cores). A section can be mounted flipped relative to the Xenium
##    one -> cs_row_flip / cs_col_flip in slides.csv. CHECK THE PLOT.
## -----------------------------------------------------------------------
m <- cs[[]]
m$core_id <- NA_character_
cs_cores <- list()
for (i in seq_len(nrow(slides))) {
  cs_id <- slides$cellscape_slide_id[i]
  sel <- m$cs_slide_id == cs_id
  map_s <- tma_map[tma_map$slide_id == slides$slide_id[i], ]
  sel_dir <- file.path(CONFIG_DIR, "core_selections", cs_id)
  if (dir.exists(sel_dir)) {
    polys <- read_core_selections(sel_dir)
    m$core_id[sel] <- assign_cells_to_polygons(m$x_um[sel], m$y_um[sel], polys)
    cores <- do.call(rbind, lapply(names(polys), function(nm) {
      p <- polys[[nm]]
      data.frame(core_id = nm, x_center = mean(p$x), y_center = mean(p$y),
                 radius = max(sqrt((p$x - mean(p$x))^2 + (p$y - mean(p$y))^2)))
    }))
  } else {
    use <- sel & m$qc_pass
    cores <- detect_cores(m$x_um[use], m$y_um[use])
    cores <- assign_grid(cores, sort_levels(map_s$core_row), sort_levels(map_s$core_col),
                         row_flip = isTRUE(slides$cs_row_flip[i]), col_flip = isTRUE(slides$cs_col_flip[i]))
    cores <- cores |> left_join(map_s[, c("core_row", "core_col", "core_id")], by = c("core_row", "core_col"))
    cores <- cores[!is.na(cores$core_id), ]
    m$core_id[sel] <- assign_cells_to_cores(m$x_um[sel], m$y_um[sel], cores)
  }
  cores$cs_slide_id <- cs_id
  cores$slide_id <- slides$slide_id[i]
  cs_cores[[cs_id]] <- cores

  d <- m[sel, ]
  if (nrow(d) > 3e5) d <- d[sample(nrow(d), 3e5), ]
  p <- ggplot(d, aes(x_um, y_um)) +
    geom_point(aes(color = is.na(core_id)), size = 0.05, stroke = 0) +
    geom_text(data = cores, aes(x_center, y_center, label = core_id), size = 2.5, fontface = "bold") +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "red"), name = "unassigned") +
    scale_y_reverse() + coord_fixed() + theme_bw() +
    ggtitle(paste0(cs_id, " (CellScape): core assignment -- must match the Xenium labels in 02"))
  save_fig(p, paste0("06_dearray_", cs_id), width = 12, height = 12)
}
cs_cores_all <- dplyr::bind_rows(cs_cores)
save_table(cs_cores_all, "06_cs_cores_all")

cs$core_id <- m$core_id
jm <- tma_map[match(cs$core_id, tma_map$core_id), ]
for (col in c("patient_id", "organ", "tissue_type", "diagnosis", "location", "include")) cs[[col]] <- jm[[col]]
keep <- colnames(cs)[!is.na(cs$core_id) & cs$qc_pass & cs$include %in% TRUE]
message("Keeping ", length(keep), " of ", ncol(cs), " CellScape cells.")
cs <- subset(cs, cells = keep)

## -----------------------------------------------------------------------
## 4. Embed + cluster on scaled arcsinh intensities (all markers except
##    DAPI etc.). Harmony on slide if >1 slide.
## -----------------------------------------------------------------------
use_m <- setdiff(rownames(cs), CS_EXCLUDE_MARKERS)
VariableFeatures(cs) <- use_m
cs <- ScaleData(cs, features = use_m, verbose = FALSE)
n_pc <- min(20, length(use_m) - 1)
cs <- RunPCA(cs, features = use_m, npcs = n_pc, verbose = FALSE, approx = FALSE)
red <- "pca"
if (USE_HARMONY && length(unique(cs$cs_slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
  cs <- harmony::RunHarmony(cs, group.by.vars = "cs_slide_id", reduction.save = "harmony", verbose = FALSE)
  red <- "harmony"
}
cs <- FindNeighbors(cs, reduction = red, dims = 1:n_pc, verbose = FALSE)
cs <- FindClusters(cs, resolution = RES_LEVEL1, cluster.name = "p1_clusters", verbose = FALSE)
cs <- RunUMAP(cs, reduction = red, dims = 1:n_pc, reduction.name = "p_umap", verbose = FALSE)

## Protein signature score = mean of scaled marker values (UCell ranks make
## little sense on a 20-40 marker panel).
score_protein <- function(obj, sigs) {
  z <- GetAssayData(obj, assay = "PROT", layer = "scale.data")
  as.data.frame(sapply(sigs, function(g) {
    g <- intersect(g, rownames(z))
    if (length(g) == 0) rep(NA_real_, ncol(z)) else Matrix::colMeans(z[g, , drop = FALSE])
  }, simplify = TRUE), row.names = colnames(z))
}

## -----------------------------------------------------------------------
## 5. Level 1 + level 2 annotation (same editable-map workflow as 04)
## -----------------------------------------------------------------------
## one good antibody is enough to define a protein lineage (unlike RNA probes)
p1 <- restrict_to_panel(markers_p$level1, rownames(cs), "Protein level 1", min_reliable = 1)
save_table(p1$coverage, "06_protein_panel_coverage_level1")
sc1 <- score_protein(cs, p1$signatures)
colnames(sc1) <- paste0("score_", colnames(sc1))
map1 <- resolve_cluster_map(cs$p1_clusters, sc1,
                            file.path(ANNOT_DIR, "cellscape_level1_cluster_map.csv"),
                            top_markers = top_markers_per_cluster(cs, "p1_clusters", n = 4))
cs$lineage <- unname(map1[as.character(cs$p1_clusters)])
cs$cell_type <- cs$lineage
save_fig(DotPlot(cs, features = use_m, group.by = "p1_clusters", assay = "PROT") + RotatedAxis() +
           ggtitle("CellScape level-1 clusters"), "06_level1_dotplot", width = 14, height = 7)

for (par in names(markers_p$level2)) {
  sigs <- restrict_to_panel(markers_p$level2[[par]], rownames(cs), par, min_reliable = 1)$signatures
  cells <- colnames(cs)[cs$lineage == par]
  if (length(cells) < 200 || length(sigs) < 2) next
  sub <- CreateSeuratObject(counts = GetAssayData(cs, assay = "PROT", layer = "counts")[, cells],
                            assay = "PROT", meta.data = cs[[]][cells, c("cs_slide_id", "core_id")])
  sub <- SetAssayData(sub, assay = "PROT", layer = "data",
                      new.data = GetAssayData(cs, assay = "PROT", layer = "data")[, cells])
  sub <- ScaleData(sub, features = use_m, verbose = FALSE)
  sub <- RunPCA(sub, features = use_m, npcs = n_pc, verbose = FALSE, approx = FALSE)
  sub <- FindNeighbors(sub, dims = 1:n_pc, verbose = FALSE)
  sub <- FindClusters(sub, resolution = RES_LEVEL2, cluster.name = "p2_clusters", verbose = FALSE)
  sc2 <- score_protein(sub, sigs)
  colnames(sc2) <- paste0("score_", colnames(sc2))
  map2 <- resolve_cluster_map(sub$p2_clusters, sc2,
                              file.path(ANNOT_DIR, sprintf("cellscape_level2_%s_cluster_map.csv", par)),
                              top_markers = top_markers_per_cluster(sub, "p2_clusters", n = 4))
  cs$cell_type[cells] <- unname(map2[as.character(sub$p2_clusters)])
  save_fig(DotPlot(sub, features = use_m, group.by = "p2_clusters", assay = "PROT") + RotatedAxis() +
             ggtitle(paste0("CellScape level-2: ", par)), paste0("06_level2_", par), width = 14, height = 6)
}

p_umap <- DimPlot(cs, reduction = "p_umap", group.by = "lineage", raster = TRUE) |
  DimPlot(cs, reduction = "p_umap", group.by = "cell_type", raster = TRUE)
save_fig(p_umap, "06_cellscape_umap", width = 16, height = 7)

## Lineage names must match the RNA side, otherwise 07/08 silently compare nothing.
rna_lineages <- names(read_marker_panel(MARKERS_RNA_CSV)$level1)
odd <- setdiff(unique(cs$lineage), rna_lineages)
if (length(odd) > 0) warning("CellScape lineages not used on the RNA side: ", paste(odd, collapse = ", "),
                             " -- rename in the cluster map so 07/08 can match them.")

## -----------------------------------------------------------------------
## 6. Export for registration (07/07b)
## -----------------------------------------------------------------------
pairs <- read_config_csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pairs$protein <- as_seurat_features(pairs$protein)
pair_prot <- intersect(pairs$protein, rownames(cs))
dat <- GetAssayData(cs, assay = "PROT", layer = "data")
for (sid in unique(cs$slide_id)) {
  cells <- colnames(cs)[cs$slide_id == sid]
  out <- data.frame(cell = cells, cs_slide_id = cs$cs_slide_id[cells], slide_id = sid,
                    core_id = cs$core_id[cells], x_um = cs$x_um[cells], y_um = cs$y_um[cells],
                    lineage = cs$lineage[cells], cell_type = cs$cell_type[cells])
  out <- cbind(out, as.data.frame(as.matrix(Matrix::t(dat[pair_prot, cells, drop = FALSE])), check.names = FALSE))
  f <- file.path(REG_DIR, paste0("cellscape_cells_", sid, ".csv.gz"))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(out, f) else utils::write.csv(out, gzfile(f), row.names = FALSE)
  message("Exported ", nrow(out), " CellScape cells -> ", f)
}

saveRDS(cs, file.path(PROCESSED_DIR, "cellscape_annotated.rds"))
message("Saved -> ", file.path(PROCESSED_DIR, "cellscape_annotated.rds"))
