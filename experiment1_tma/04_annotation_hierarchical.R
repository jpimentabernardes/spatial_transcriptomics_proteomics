## =============================================================================
## 04_annotation_hierarchical.R
##
## Cell type calling for multi-organ TMAs (TMA_ORGAN + TMA_TUMOR together), in
## levels, each clustered on its own and named by YOU:
##
##   L1  lineage, all cells of both TMAs     -> Immune / Epithelial / Mesenchymal /
##                                              Endothelial / Neural_NE (/ Low_quality)
##   L2  Immune, all organs together          -> one consistent set of immune labels
##       Stroma (Mesenchymal + Endothelial), all organs together
##       Epithelial + Neural_NE, per FAMILY (config/annotation_families.csv: organ or
##       diagnosis -> family, normal and tumour together, e.g. GI_glandular, Squamous)
##
## Per group: clustering (Harmony on slide only), FindAllMarkers, a dot plot of the
## group's signature markers (config/marker_panel_rna.csv) + top markers, UMAP,
## which organs each cluster comes from, and a template you fill in:
##   annotation/04_<seg>_<group>.csv   cluster, n_cells, top_markers, top_organs,
##                                     suggestion, margin, label   <- edit `label`
## The label column is pre-filled with the marker-based suggestion; correct it and
## re-run 04. Your labels are kept as long as the clustering is unchanged (same
## clusters with the same cell numbers); if it changed, the old file is kept as
## *.previous.csv and a fresh template is written.
##
## Final columns: lineage (L1), family, annotation_group, cell_type.
## Epithelial cells of cores without a family (lymphoid organs, muscle, fat, ...)
## and families with < MIN_CELLS_GROUP cells get "<lineage>_<family>" labels.
##
## INPUT : xenium_<seg>_cores.rds (02), config/marker_panel_rna.csv,
##         config/annotation_families.csv (or a `family` column in tma_map.csv)
## OUTPUT: xenium_<seg>_annotated.rds, annotation/04_<seg>_*.csv, tables/04_*.csv,
##         figures/04_*.png, data/registration/xenium_cells_<slide>.csv.gz
## =============================================================================

source("config/config.R")
source("R/helpers.R")

ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))
RES_LINEAGE    <- get0("RES_LINEAGE",    ifnotfound = 0.3)
RES_IMMUNE     <- get0("RES_IMMUNE",     ifnotfound = 0.8)
RES_STROMA     <- get0("RES_STROMA",     ifnotfound = 0.5)
RES_EPITHELIAL <- get0("RES_EPITHELIAL", ifnotfound = 0.6)
MIN_CELLS_GROUP <- get0("MIN_CELLS_GROUP", ifnotfound = 500)
## families whose epithelium is sub-clustered (NULL = all); others keep "<lineage>_<family>"
ANNOT_FAMILIES_DETAIL <- get0("ANNOT_FAMILIES_DETAIL", ifnotfound = NULL)
LINEAGES <- c("Immune", "Epithelial", "Mesenchymal", "Endothelial", "Neural_NE")

xen <- readRDS(obj_path("cores"))
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)
DefaultAssay(xen) <- "Xenium"

## -----------------------------------------------------------------------
## 0. Family per cell (from the core's organ / diagnosis)
## -----------------------------------------------------------------------
key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
if (!"family" %in% colnames(xen[[]]) || all(is.na(xen$family))) {
  fam <- read_config_csv(ANNOT_FAMILIES_CSV, stringsAsFactors = FALSE, fileEncoding = "latin1")
  xen$family <- fam$family[match(key(xen$organ), key(fam$organ))]
  unmapped <- unique(xen$organ[is.na(xen$family) & !is.na(xen$organ)])
  if (length(unmapped)) message("Organs without a family in ", basename(ANNOT_FAMILIES_CSV), " (-> 'none'): ",
                                paste(unmapped, collapse = "; "))
}
xen$family[is.na(xen$family) | !nzchar(xen$family)] <- "none"
print(table(family = xen$family, slide = xen$slide_id))

## clustering / suggestion / template / figure helpers: R/annotation_helpers.R (shared with 04b)
source("R/annotation_helpers.R")

## -----------------------------------------------------------------------
## 1. Normalise all cells once (log-normalised data reused by every level)
## -----------------------------------------------------------------------
cov <- restrict_to_panel(markers$level1, panel_genes, "Level 1")
save_table(cov$coverage, "04_panel_coverage_level1")
xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)
ANNOT_GLOBAL_MEDIAN <- stats::median(xen$nCount_Xenium)   # low-information flag: clusters vs all cells

## -----------------------------------------------------------------------
## 2. L1: lineage on all cells
## -----------------------------------------------------------------------
lin <- annotate_group(xen, "L1_lineage", markers$level1, RES_LINEAGE, comp_by = "family")
xen$lineage <- unname(lin$labels[colnames(xen)])
xen[["umap"]] <- CreateDimReducObject(lin$umap[colnames(xen), ], key = "UMAP_", assay = "Xenium")
rm(lin); gc()
bad <- setdiff(unique(xen$lineage), c(LINEAGES, "Low_quality", "Unknown"))
if (length(bad)) message("NOTE: lineage labels outside ", paste(LINEAGES, collapse = "/"), "/Low_quality: ",
                         paste(bad, collapse = ", "), " -- these cells keep the lineage as cell type.")
print(table(lineage = xen$lineage, slide = xen$slide_id))

## -----------------------------------------------------------------------
## 3. L2 groups
## -----------------------------------------------------------------------
xen$annotation_group <- dplyr::case_when(
  xen$lineage == "Immune" ~ "Immune",
  xen$lineage %in% c("Mesenchymal", "Endothelial") ~ "Stroma",
  xen$lineage %in% c("Epithelial", "Neural_NE") & xen$family != "none" &
    (is.null(ANNOT_FAMILIES_DETAIL) | xen$family %in% ANNOT_FAMILIES_DETAIL) ~ paste0("Epi_", xen$family),
  TRUE ~ "none")
xen$cell_type <- ifelse(xen$lineage %in% c("Epithelial", "Neural_NE"),
                        paste0(xen$lineage, "_", ifelse(xen$family == "none", "other", xen$family)),
                        xen$lineage)
groups <- sort(setdiff(unique(xen$annotation_group), "none"))
failed <- list()
group_tab <- as.data.frame(table(annotation_group = xen$annotation_group))
save_table(group_tab, "04_annotation_groups")
for (g in groups) {
  cells <- colnames(xen)[xen$annotation_group == g]
  if (length(cells) < MIN_CELLS_GROUP) {
    message(g, ": only ", length(cells), " cells (< MIN_CELLS_GROUP = ", MIN_CELLS_GROUP, ") -- not sub-clustered")
    next
  }
  if (g == "Immune") {
    sigs <- markers$level2[["Immune"]]; res <- RES_IMMUNE; by <- "family"
  } else if (g == "Stroma") {
    sigs <- markers$level2[["Stroma"]]; res <- RES_STROMA; by <- "family"
  } else {
    f <- sub("^Epi_", "", g)
    sigs <- c(markers$level2[[f]], markers$level2[["Immune"]]["Proliferating"])
    if (is.null(markers$level2[[f]])) message(g, ": no level-2 markers for family '", f, "' in ", basename(MARKERS_RNA_CSV))
    res <- RES_EPITHELIAL; by <- "organ"
  }
  sub <- subset(xen, cells = cells)
  lab <- tryCatch(annotate_group(sub, g, sigs, res, comp_by = by)$labels, error = function(e) {
    failed[[g]] <<- conditionMessage(e)
    message("!! ", g, " FAILED: ", conditionMessage(e), "\n   -> its cells keep '", xen$cell_type[cells[1]],
            "'-type labels; the other groups continue.")
    NULL
  })
  if (!is.null(lab)) xen$cell_type[match(names(lab), colnames(xen))] <- lab
  rm(sub); gc()
}
if (length(failed)) {
  message("\n!! Groups that failed (send this to the analyst):")
  for (g in names(failed)) message("   ", g, ": ", failed[[g]])
}
xen$celltype <- factor(xen$cell_type)
ct_tab <- xen[[]] |> count(lineage, annotation_group, cell_type, name = "n_cells")
save_table(ct_tab, "04_cell_types")
print(as.data.frame(ct_tab), row.names = FALSE)

## -----------------------------------------------------------------------
## 4. Whole-dataset figures with the final labels
## -----------------------------------------------------------------------
pal_lin <- celltype_palette(xen$lineage, xen$lineage)       # same colours as the sections / bar plots
save_fig(plot_umap_labelled(xen, "lineage", label_size = 6, title = "Lineage (L1)", pal = pal_lin), "04_lineage_umap",
         width = 14, height = 11, dpi = 180)
n_slides <- length(unique(xen$slide_id))
lin_meta <- xen[[]]; lin_meta$cell_type <- lin_meta$lineage
save_fig(plot_celltype_sections(lin_meta, panel_in = 11), "04_lineage_sections",
         width = 11 * n_slides + 4, height = 12, dpi = 220)

## immune cell types per core, cores grouped by family (the Aim 1/2 overview)
imm <- xen[[]][xen$lineage == "Immune" & !is.na(xen$core_id), ]
if (nrow(imm)) {
  pal_imm <- celltype_palette(imm$cell_type, imm$lineage)
  comp <- imm |> count(family, core_id, cell_type) |> group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
  comp$cell_type <- factor(comp$cell_type, levels = names(pal_imm))
  p_imm <- ggplot(comp, aes(core_id, frac, fill = cell_type)) + geom_col(width = 0.9) +
    facet_grid(~ family, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_imm, name = "Immune cell type") + scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11, angle = 90, hjust = 0)) +
    labs(x = "cores", y = "fraction of immune cells", title = "Immune composition per core, by family")
  save_fig(p_imm, "04_immune_composition_per_core", width = 26, height = 10, dpi = 180)
}
lin_comp <- xen[[]] |> filter(!is.na(core_id)) |> count(family, core_id, lineage) |>
  group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
lin_comp$lineage <- factor(lin_comp$lineage, levels = names(pal_lin))
p_lin <- ggplot(lin_comp, aes(core_id, frac, fill = lineage)) + geom_col(width = 0.9) +
  scale_fill_manual(values = pal_lin, name = "Lineage") +
  facet_grid(~ family, scales = "free_x", space = "free_x") + scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11, angle = 90, hjust = 0)) +
  labs(x = "cores", y = "fraction of cells", title = "Lineage composition per core, by family")
save_fig(p_lin, "04_lineage_composition_per_core", width = 26, height = 9, dpi = 180)

## -----------------------------------------------------------------------
## 5. Export for registration / cross-modal alignment (07, 07b)
## -----------------------------------------------------------------------
pairs <- read_config_csv(RNA_PROT_PAIRS_CSV, stringsAsFactors = FALSE)
pair_genes <- intersect(match_case(as_seurat_features(pairs$gene), panel_genes), panel_genes)
expr <- GetAssayData(xen, assay = "Xenium", layer = "data")[pair_genes, , drop = FALSE]
for (sid in unique(xen$slide_id)) {
  cells <- colnames(xen)[xen$slide_id == sid]
  out <- data.frame(cell = cells, slide_id = sid, core_id = xen$core_id[cells],
                    x_um = xen$x_um[cells], y_um = xen$y_um[cells], family = xen$family[cells],
                    lineage = xen$lineage[cells], cell_type = xen$cell_type[cells])
  if (length(pair_genes)) out <- cbind(out, as.data.frame(as.matrix(Matrix::t(expr[, cells, drop = FALSE]))))
  f <- file.path(REG_DIR, paste0("xenium_cells_", sid, ".csv.gz"))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(out, f) else
    utils::write.csv(out, gzfile(f), row.names = FALSE)
  message("Exported ", nrow(out), " cells -> ", f)
}

saveRDS(xen, obj_path("annotated"))
message("Saved -> ", obj_path("annotated"))
