## =============================================================================
## 04c_annotation_per_organ_all_lineages.R
##
## Per-organ annotation of ALL lineages, after 04 (and optionally 04b): for every
## site (config/annotation_families.csv, `site`) three separate clusterings, each
## with its own marker set, normal only and normal + matched tumours:
##
##   Epithelium  Epithelial + Neural_NE    organ-family markers   -> SAME groups and label
##                                                                   files as 04b (04b_<seg>_Organ_<site>.csv),
##                                                                   so your 04b labels are reused
##   Immune      Immune                    immune markers (L2)    -> 04c_<seg>_Organ_<site>_Immune.csv
##   Stroma      Mesenchymal + Endothelial stroma markers (L2)    -> 04c_<seg>_Organ_<site>_Stroma.csv
##
## (and OrganTumour_<site>[_Immune|_Stroma] for normal + tumour). Like 04: edit the
## `label` column of a template and re-run; labels are kept while the clustering is
## unchanged. Groups below MIN_CELLS_ORGAN cells are skipped (their cells keep the 04
## label). No Harmony (tumour and normal are on different slides).
##
## Final column cell_type_organ_all: ONE organ-level label per cell =
##   normal cores  -> the normal-only clustering of its site and lineage group
##   tumour cores  -> the normal + tumour clustering
##   otherwise     -> the 04 label (cell_type); cell_type_organ_source says which.
## The 04 labels stay in `cell_type` for comparisons ACROSS organs (one shared
## immune / stroma clustering); the organ-level labels are for within-organ detail.
##
## Run:  Rscript 04c_annotation_per_organ_all_lineages.R
##       EXP1_04C_GROUPS=Immune,Stroma Rscript ...   (skip the epithelium, e.g. after 04b)
##       EXP1_04B_MODES=normal Rscript ...           (normal only)
##
## INPUT : xenium_<seg>_annotated.rds (04 / 04b), config/annotation_families.csv,
##         config/marker_panel_rna.csv
## OUTPUT: columns cell_type_organ_<Epithelium|Immune|Stroma>, cell_type_organ_tumour_<...>,
##         cell_type_organ_all, cell_type_organ_source in xenium_<seg>_annotated.rds;
##         annotation/04c_<seg>_cell_types_per_organ_all.csv.gz, tables/04c_*.csv,
##         figures/04b_* (epithelium) and 04c_* (immune, stroma): dot plot, UMAP,
##         UMAP by core, clusters by organ, composition per core (with tumour)
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/annotation_helpers.R")

ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))
ORGAN_SITES     <- get0("ORGAN_SITES", ifnotfound = NULL)
RES_ORGAN       <- get0("RES_ORGAN", ifnotfound = 0.6)
RES_IMMUNE      <- get0("RES_IMMUNE", ifnotfound = 0.8)
RES_STROMA      <- get0("RES_STROMA", ifnotfound = 0.5)
MIN_CELLS_ORGAN <- get0("MIN_CELLS_ORGAN", ifnotfound = 300)
MODES  <- trimws(strsplit(Sys.getenv("EXP1_04B_MODES", "normal,with_tumour"), ",")[[1]])
GROUPS <- trimws(strsplit(Sys.getenv("EXP1_04C_GROUPS", "Epithelium,Immune,Stroma"), ",")[[1]])

xen <- readRDS(obj_path("annotated"))
if (!all(c("lineage", "cell_type") %in% colnames(xen[[]])))
  stop("Run 04_annotation_hierarchical.R to the end first (lineage / cell_type missing).")
DefaultAssay(xen) <- "Xenium"
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)
if (!"data" %in% Layers(xen[["Xenium"]]))
  xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)

xen <- add_site_and_tumour(xen, ANNOT_FAMILIES_CSV)
sites <- pick_sites(xen, ORGAN_SITES)

## lineage groups: cells, markers, resolution, file prefix / suffix
epi_sigs <- function(cells) {
  fams <- unique(xen$family[match(cells, colnames(xen))])
  s <- do.call(c, unname(markers$level2[intersect(fams, names(markers$level2))]))
  s <- c(s, markers$level2[["Immune"]]["Proliferating"])
  s[!duplicated(names(s))]
}
GROUP_DEF <- list(
  Epithelium = list(lineages = c("Epithelial", "Neural_NE"), sigs = epi_sigs, res = RES_ORGAN,
                    prefix = "04b", suffix = "", pal = "Epithelial"),
  Immune     = list(lineages = "Immune", sigs = function(cells) markers$level2[["Immune"]], res = RES_IMMUNE,
                    prefix = "04c", suffix = "_Immune", pal = "Immune"),
  Stroma     = list(lineages = c("Mesenchymal", "Endothelial"), sigs = function(cells) markers$level2[["Stroma"]],
                    res = RES_STROMA, prefix = "04c", suffix = "_Stroma", pal = "Mesenchymal"))
bad <- setdiff(GROUPS, names(GROUP_DEF))
if (length(bad)) stop("EXP1_04C_GROUPS: unknown group(s) ", paste(bad, collapse = ", "),
                      " -- use ", paste(names(GROUP_DEF), collapse = ", "))
if (is.null(markers$level2[["Stroma"]])) message("NOTE: no level-2 'Stroma' markers in ", basename(MARKERS_RNA_CSV))

cells_tab <- xen[[]][xen$site %in% sites, ] |>
  mutate(group = dplyr::case_when(lineage %in% GROUP_DEF$Epithelium$lineages ~ "Epithelium",
                                  lineage %in% GROUP_DEF$Immune$lineages ~ "Immune",
                                  lineage %in% GROUP_DEF$Stroma$lineages ~ "Stroma", TRUE ~ "other")) |>
  count(site, tumour = ifelse(is_tumour, "tumour", "normal"), group, name = "n_cells") |>
  tidyr::pivot_wider(names_from = group, values_from = n_cells, values_fill = 0)
save_table(cells_tab, "04c_cells_per_site_and_group")
message("Cells per site and lineage group (groups < ", MIN_CELLS_ORGAN, " cells are skipped):")
print(as.data.frame(cells_tab), row.names = FALSE)

## -----------------------------------------------------------------------
## 1. per group: annotate every site (normal / normal + tumour)
## -----------------------------------------------------------------------
failed <- list()
for (gname in GROUPS) {
  gd <- GROUP_DEF[[gname]]
  message("\n################ ", gname, " (", paste(gd$lineages, collapse = " + "), ") ################")
  res <- annotate_sites(xen, sites, gd$lineages, gd$sigs, gd$res, prefix = gd$prefix, suffix = gd$suffix,
                        modes = MODES, min_cells = MIN_CELLS_ORGAN, what = gname, pal_lineage = gd$pal)
  xen <- set_labels(xen, paste0("cell_type_organ_", gname), res$normal)
  xen <- set_labels(xen, paste0("cell_type_organ_tumour_", gname), res$with_tumour)
  failed <- c(failed, res$failed)
}
if (length(failed)) {
  message("\n!! Groups that failed (send this to the analyst):")
  for (g in names(failed)) message("   ", g, ": ", failed[[g]])
}

## -----------------------------------------------------------------------
## 2. one organ-level label per cell
## -----------------------------------------------------------------------
m <- xen[[]]
grp_of <- dplyr::case_when(m$lineage %in% GROUP_DEF$Epithelium$lineages ~ "Epithelium",
                           m$lineage %in% GROUP_DEF$Immune$lineages ~ "Immune",
                           m$lineage %in% GROUP_DEF$Stroma$lineages ~ "Stroma", TRUE ~ NA_character_)
pick <- function(prefix) {
  v <- rep(NA_character_, nrow(m))
  for (gname in GROUPS) {
    col <- paste0(prefix, gname)
    if (col %in% colnames(m)) v[grp_of %in% gname] <- m[[col]][grp_of %in% gname]
  }
  v
}
lab_normal <- pick("cell_type_organ_")
lab_tumour <- pick("cell_type_organ_tumour_")
organ_lab <- ifelse(m$is_tumour, lab_tumour, lab_normal)
xen$cell_type_organ_source <- ifelse(is.na(organ_lab), "04 (shared)", "organ")
xen$cell_type_organ_all <- ifelse(is.na(organ_lab), as.character(m$cell_type), organ_lab)

m <- xen[[]]
tab <- m |> filter(site %in% sites) |>
  count(site, tumour = ifelse(is_tumour, "tumour", "normal"), lineage, cell_type_organ_all, cell_type_organ_source,
        name = "n_cells")
save_table(tab, "04c_cell_types_per_organ_all")
keep <- c("slide_id", "core_id", "organ", "site", "is_tumour", "x_um", "y_um", "lineage", "cell_type",
          paste0("cell_type_organ_", GROUPS), paste0("cell_type_organ_tumour_", GROUPS),
          "cell_type_organ_all", "cell_type_organ_source")
out <- data.frame(cell = rownames(m), m[, intersect(keep, colnames(m))])
f <- file.path(ANNOT_DIR, sprintf("04c_%s_cell_types_per_organ_all.csv.gz", SEGMENTATION))
utils::write.csv(out, gzfile(f), row.names = FALSE)
message("Per-cell table -> ", f)
message("Organ-level label for ", round(100 * mean(xen$cell_type_organ_source == "organ"), 1),
        "% of cells; the rest keep the 04 label.")

saveRDS(xen, obj_path("annotated"))
message("Saved (cell_type_organ_* columns added) -> ", obj_path("annotated"))
