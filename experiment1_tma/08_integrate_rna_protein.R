## =============================================================================
## 08_integrate_rna_protein.R
##
## Spatial transcriptomics (Xenium) <-> spatial proteomics (CellScape) on the
## registered coordinates from 07/07b.
##
## The sections are CONSECUTIVE, so a Xenium cell and its nearest CellScape
## cell are generally NOT the same cell (~5 um apart in z, different nuclei
## cut). Integration is therefore done at three levels, from most to least
## dependent on registration accuracy:
##
##   A. cell level   -- nearest CellScape cell within NN_MAX_DIST_UM (label
##                      transfer, confusion matrix) and neighbourhood-averaged
##                      protein within NBHD_RADIUS_UM added to each Xenium cell
##                      as assay "PROTnbhd". Use for spatial visualisation and
##                      niche-level questions, not single-cell correlations.
##   B. spatial bins -- BIN_SIZE_UM grid inside each core: RNA vs protein of
##                      each marker pair, Spearman per core.
##   C. core x cell type -- RNA mean in Xenium T cells of core X vs protein
##                      mean in CellScape T cells of core X. Needs NO
##                      registration, only matched cores + shared labels. This
##                      is the most robust comparison, and it is the one that
##                      answers the Aim 1 question "is the RNA probe signal
##                      real in organ X?" with an orthogonal readout.
##
## INPUT : xenium_<seg>_insights.rds (05) or _annotated.rds (04)
##         cellscape_annotated.rds (06)
##         data/registration/cellscape_aligned[_nonrigid]_<slide>.csv.gz (07/07b)
##         tables/05_aim1_sbr_gene_by_organ.csv (05, optional)
## OUTPUT: xenium_<seg>_integrated.rds, tables/08_*.csv, figures/08_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen_path <- if (file.exists(obj_path("insights"))) obj_path("insights") else obj_path("annotated")
xen <- readRDS(xen_path)
cs <- readRDS(file.path(PROCESSED_DIR, "cellscape_annotated.rds"))
pairs <- utils::read.csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pairs$gene <- as_seurat_features(pairs$gene)
pairs$protein <- as_seurat_features(pairs$protein)
pairs <- pairs[pairs$gene %in% rownames(xen[["Xenium"]]) & pairs$protein %in% rownames(cs[["PROT"]]), ]
message(nrow(pairs), " RNA/protein pairs on both panels: ", paste(pairs$pair_id, collapse = ", "))

## -----------------------------------------------------------------------
## 0. Registered CellScape coordinates (non-rigid if 07b was run)
## -----------------------------------------------------------------------
aligned <- dplyr::bind_rows(lapply(unique(xen$slide_id), function(sid) {
  f_nr <- file.path(REG_DIR, paste0("cellscape_aligned_nonrigid_", sid, ".csv.gz"))
  f_af <- file.path(REG_DIR, paste0("cellscape_aligned_", sid, ".csv.gz"))
  if (file.exists(f_nr)) {
    d <- utils::read.csv(f_nr); d$x_reg <- d$x_final; d$y_reg <- d$y_final; d$reg <- "nonrigid"
  } else if (file.exists(f_af)) {
    d <- utils::read.csv(f_af); d$x_reg <- d$x_aligned; d$y_reg <- d$y_aligned; d$reg <- "affine"
  } else { message("No registration for ", sid, " -- run 07 first"); return(NULL) }
  d[, c("cell", "slide_id", "core_id", "x_reg", "y_reg", "reg")]
}))
aligned <- aligned[aligned$cell %in% colnames(cs), ]
cs_meta <- cs[[]][aligned$cell, ]
cs_meta$x_reg <- aligned$x_reg
cs_meta$y_reg <- aligned$y_reg
prot <- GetAssayData(cs, assay = "PROT", layer = "data")

## =============================================================================
## A. Cell level: nearest-cell transfer + neighbourhood protein
## =============================================================================
xm <- xen[[]]
xen_cells <- colnames(xen)
nn_dist <- rep(NA_real_, length(xen_cells)); names(nn_dist) <- xen_cells
nn_cell <- rep(NA_character_, length(xen_cells)); names(nn_cell) <- xen_cells
n_nbhd <- rep(0L, length(xen_cells)); names(n_nbhd) <- xen_cells
nbhd <- matrix(0, nrow = nrow(prot), ncol = length(xen_cells), dimnames = list(rownames(prot), xen_cells))

for (cid in intersect(unique(xm$core_id), unique(cs_meta$core_id))) {
  xc <- xen_cells[xm$core_id == cid]
  cc <- rownames(cs_meta)[cs_meta$core_id == cid]
  if (length(xc) < 10 || length(cc) < 10) next
  cxy <- as.matrix(cs_meta[cc, c("x_reg", "y_reg")])
  qxy <- as.matrix(xm[xc, c("x_um", "y_um")])
  nn1 <- RANN::nn2(cxy, qxy, k = 1)
  nn_dist[xc] <- nn1$nn.dists[, 1]
  nn_cell[xc] <- cc[nn1$nn.idx[, 1]]
  k <- min(50, length(cc))
  nr <- RANN::nn2(cxy, qxy, k = k, searchtype = "radius", radius = NBHD_RADIUS_UM)
  idx <- nr$nn.idx
  i <- rep(seq_along(xc), k)[as.vector(idx) > 0]
  j <- as.vector(idx)[as.vector(idx) > 0]
  if (length(i) == 0) next
  W <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(length(xc), length(cc)))
  cnt <- Matrix::rowSums(W)
  W <- Matrix::Diagonal(x = 1 / pmax(cnt, 1)) %*% W
  nbhd[, xc] <- as.matrix(prot[, cc, drop = FALSE] %*% Matrix::t(W))
  n_nbhd[xc] <- as.integer(cnt)
}

matched <- !is.na(nn_dist) & nn_dist <= NN_MAX_DIST_UM
xen$cs_nn_dist <- nn_dist
cs_lineage <- stats::setNames(cs[["lineage", drop = TRUE]], colnames(cs))
cs_type <- stats::setNames(cs[["cell_type", drop = TRUE]], colnames(cs))
xen$cs_nn_lineage <- ifelse(matched, unname(cs_lineage[nn_cell]), NA)
xen$cs_nn_cell_type <- ifelse(matched, unname(cs_type[nn_cell]), NA)
xen$n_cs_nbhd <- n_nbhd
xen[["PROTnbhd"]] <- CreateAssayObject(data = Matrix::Matrix(nbhd, sparse = TRUE))
message(round(100 * mean(matched), 1), "% of Xenium cells have a CellScape cell within ", NN_MAX_DIST_UM, " um")
rm(nbhd); invisible(gc())

## Cross-modal lineage agreement vs chance (random pairing within the core).
chance_agreement <- function(a, b) {
  pa <- prop.table(table(a)); pb <- prop.table(table(b))
  l <- intersect(names(pa), names(pb))
  sum(pa[l] * pb[l])
}
agree <- xen[[]] |>
  filter(!is.na(cs_nn_lineage)) |>
  group_by(core_id, organ) |>
  summarise(n = n(), agreement = mean(lineage == cs_nn_lineage),
            expected = chance_agreement(lineage, cs_nn_lineage), .groups = "drop") |>
  mutate(lift = agreement / expected)
save_table(agree, "08_crossmodal_lineage_agreement_by_core")

conf <- xen[[]] |> filter(!is.na(cs_nn_lineage)) |> count(lineage, cs_nn_lineage) |>
  group_by(lineage) |> mutate(frac = n / sum(n)) |> ungroup()
p_conf <- ggplot(conf, aes(cs_nn_lineage, lineage, fill = frac)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
  scale_fill_viridis_c() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "nearest CellScape cell (protein)", y = "Xenium cell (RNA)",
       title = "Cross-modal lineage agreement (row-normalised). Diagonal = consistent annotation + registration")
save_fig(p_conf, "08_crossmodal_lineage_confusion", width = 9, height = 7)

## =============================================================================
## B. Spatial bins: RNA vs protein per marker pair, per core
## =============================================================================
rna <- GetAssayData(xen, assay = "Xenium", layer = "data")
bin_id <- function(core, x, y) paste(core, floor(x / BIN_SIZE_UM), floor(y / BIN_SIZE_UM), sep = "|")
xm <- xen[[]]
xb <- bin_id(xm$core_id, xm$x_um, xm$y_um)
cb <- bin_id(cs_meta$core_id, cs_meta$x_reg, cs_meta$y_reg)
bins_ok <- intersect(names(which(table(xb) >= 5)), names(which(table(cb) >= 5)))

bin_mean <- function(mat, bins, keep) {
  sel <- bins %in% keep
  f <- factor(bins[sel], levels = keep)
  s <- mat[, sel, drop = FALSE] %*% Matrix::t(Matrix::fac2sparse(f, drop.unused.levels = FALSE))
  sweep(as.matrix(s), 2, as.numeric(table(f)), "/")
}
rna_bins <- bin_mean(rna[pairs$gene, , drop = FALSE], xb, bins_ok)
prot_bins <- bin_mean(prot[pairs$protein, rownames(cs_meta), drop = FALSE], cb, bins_ok)
bin_core <- sub("\\|.*$", "", bins_ok)

bin_cor <- dplyr::bind_rows(lapply(unique(bin_core), function(cid) {
  sel <- bin_core == cid
  if (sum(sel) < 8) return(NULL)
  data.frame(core_id = cid, pair_id = pairs$pair_id, n_bins = sum(sel),
             spearman = vapply(seq_len(nrow(pairs)), function(k) {
               suppressWarnings(stats::cor(rna_bins[pairs$gene[k], sel], prot_bins[pairs$protein[k], sel],
                                           method = "spearman"))
             }, numeric(1)))
}))
bin_cor$organ <- xm$organ[match(bin_cor$core_id, xm$core_id)]
save_table(bin_cor, "08_rna_protein_bin_correlation_by_core")

p_bin <- bin_cor |> group_by(pair_id, organ) |> summarise(rho = median(spearman, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(organ, pair_id, fill = rho)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("RNA vs protein, ", BIN_SIZE_UM, " um bins (median Spearman across cores)"))
save_fig(p_bin, "08_rna_protein_bin_correlation", width = 10, height = 7)

## =============================================================================
## C. Core x cell type: registration-free RNA vs protein comparison
## =============================================================================
## Harmonised label: the cell type where BOTH modalities use that name (e.g.
## CD8_T), otherwise the lineage. Without this, a compartment sub-typed
## differently on the two sides (RNA: goblet/absorptive; protein: just
## "Epithelial") silently drops out and the correlation loses its range.
shared_types <- intersect(unique(xen$cell_type), unique(cs$cell_type))
harmonised_label <- function(meta) ifelse(meta$cell_type %in% shared_types, meta$cell_type, meta$lineage)
message("Core x cell type comparison; shared cell types: ", paste(shared_types, collapse = ", "),
        " (all other cells compared at lineage level)")

group_means <- function(mat, meta, feats) {
  g <- paste(meta$core_id, harmonised_label(meta), sep = "||")
  keep <- names(which(table(g) >= MIN_CELLS_PSEUDOBULK))
  bin_mean(mat[feats, rownames(meta), drop = FALSE], g, keep)
}
rna_ct <- group_means(rna, xen[[]], pairs$gene)
prot_ct <- group_means(prot, cs[[]], pairs$protein)
common_g <- intersect(colnames(rna_ct), colnames(prot_ct))
ct_long <- dplyr::bind_rows(lapply(seq_len(nrow(pairs)), function(k) {
  data.frame(pair_id = pairs$pair_id[k], group = common_g,
             rna = rna_ct[pairs$gene[k], common_g], protein = prot_ct[pairs$protein[k], common_g])
}))
ct_long$core_id <- sub("\\|\\|.*$", "", ct_long$group)
ct_long$label <- sub("^.*\\|\\|", "", ct_long$group)
ct_long$organ <- xm$organ[match(ct_long$core_id, xm$core_id)]
save_table(ct_long, "08_rna_protein_core_celltype_means")

## Pearson on the (already log / arcsinh) values is the headline metric: a
## marker specific to one compartment has one high group and many ~0 groups,
## and Spearman would mostly rank noise among the ~0 groups.
ct_cor <- ct_long |>
  group_by(pair_id, organ) |>
  summarise(n_groups = n(),
            pearson = if (n() >= 5) suppressWarnings(stats::cor(rna, protein)) else NA_real_,
            spearman = if (n() >= 5) suppressWarnings(stats::cor(rna, protein, method = "spearman")) else NA_real_,
            .groups = "drop")
save_table(ct_cor, "08_rna_protein_core_celltype_correlation")

p_ct <- ggplot(ct_long, aes(rna, protein, color = label)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~ pair_id, scales = "free") + theme_bw() +
  labs(x = "RNA (mean log-norm expression)", y = "protein (mean arcsinh intensity)",
       title = "RNA vs protein per core x cell type (no registration needed)")
save_fig(p_ct, "08_rna_protein_core_celltype_scatter", width = 16, height = 12)

## Composition agreement per core (lineage level)
comp_r <- xen[[]] |> count(core_id, lineage) |> group_by(core_id) |> mutate(rna = n / sum(n)) |> select(-n)
comp_p <- cs[[]] |> count(core_id, lineage) |> group_by(core_id) |> mutate(protein = n / sum(n)) |> select(-n)
comp <- inner_join(comp_r, comp_p, by = c("core_id", "lineage"))
p_comp <- ggplot(comp, aes(rna, protein, color = lineage)) +
  geom_abline(linetype = 2) + geom_point() + theme_bw() + coord_equal() +
  labs(x = "fraction of cells (Xenium)", y = "fraction of cells (CellScape)",
       title = "Per-core lineage composition: RNA vs protein section")
save_fig(p_comp, "08_composition_rna_vs_protein", width = 8, height = 7)

## =============================================================================
## Aim 1 tie-in: RNA detection (05) x protein concordance (C) per organ
## =============================================================================
sbr_f <- file.path(TABLES_DIR, "05_aim1_sbr_gene_by_organ.csv")
if (file.exists(sbr_f)) {
  sbr_org <- utils::read.csv(sbr_f, stringsAsFactors = FALSE)
  val <- ct_cor |>
    left_join(pairs[, c("pair_id", "gene", "protein")], by = "pair_id") |>
    left_join(sbr_org[, c("gene", "organ", "frac_cores_detected", "median_log2_sbr")], by = c("gene", "organ")) |>
    mutate(call = case_when(
      is.na(pearson) ~ "too few groups",
      frac_cores_detected >= 0.5 & pearson >= 0.5 ~ "RNA detected, protein-concordant (validated)",
      frac_cores_detected >= 0.5 & pearson < 0.2 ~ "RNA detected, protein-discordant (check specificity)",
      frac_cores_detected < 0.5 & pearson >= 0.5 ~ "RNA weak but tracks protein (sensitivity-limited)",
      TRUE ~ "inconclusive"
    ))
  save_table(val, "08_aim1_probe_validation_by_protein")
  print(as.data.frame(val[, c("pair_id", "organ", "frac_cores_detected", "pearson", "call")]))
}

saveRDS(xen, obj_path("integrated"))
message("Saved -> ", obj_path("integrated"))
