## =============================================================================
## 04b_annotation_per_organ.R
##
## Epithelial (+ neuroendocrine) cell types annotated PER ORGAN, after 04:
## colon separate from rectum, ileum, stomach, ... Each organ ("site") is run
## twice:
##   normal       only the normal cores of that organ (TMA_ORGAN)
##   with_tumour  the normal cores + the tumours of that organ (TMA_TUMOR),
##                clustered together -> which tumour clusters resemble which
##                normal cell types, and which are tumour-only
## Sites come from the `site` column of config/annotation_families.csv, which
## pairs each organ with its tumours, for all organ systems, e.g.
##   Colon  = colon transversum + sigmoid mucosa + colon / sigmoid adenocarcinoma + colorectal NEC
##   Lung   = bronchus + lung + NSCLC adeno, lung SCC, SCLC, lung NET
##   Skin   = skin + extramammary Paget, nodular melanoma, Merkel cell carcinoma
##   Brain, Oral_cavity, Salivary_gland, Tonsil, Thymus, Breast, Stomach, Liver, Pancreas, ...
## Normal only (no matched tumour on the TMA): cervix, endometrium, ovary, placenta,
## kidney, prostate, testis, thyroid, ...; tumour only: Larynx, Nasal_cavity, Pleura.
## Edit the column to split or merge sites (e.g. Colon_sigmoid vs Colon_transversum);
## organs left empty (muscle, fat, vessels, spleen, lymph node: no epithelium) are not run.
## ORGAN_SITES in config.R picks a subset, e.g. ORGAN_SITES <- c("Colon", "Lung").
##
## Immune and stromal cells keep their 04 labels (one consistent set across all
## organs, so organs stay comparable); set ORGAN_LINEAGES to sub-cluster them
## per organ as well.
##
## No Harmony here: normal and tumour cores are on different slides, so
## integrating on slide would remove exactly the tumour-vs-normal difference.
##
## Per site and mode, like 04: clustering, FindAllMarkers, dot plot, UMAP,
## clusters by organ / diagnosis, and a label template you edit and re-run:
##   annotation/04b_<seg>_Organ_<site>.csv         (normal)
##   annotation/04b_<seg>_OrganTumour_<site>.csv   (normal + tumour)
## plus figures/04b_<group>_umap_by_core.png (UMAP coloured by core, normal/tumour and label)
## and figures/04b_OrganTumour_<site>_composition_per_core.png.
##
## Run:  Rscript 04b_annotation_per_organ.R           (both modes)
##       EXP1_04B_MODES=normal Rscript 04b_annotation_per_organ.R
##
## INPUT : xenium_<seg>_annotated.rds (04), config/annotation_families.csv (site),
##         config/marker_panel_rna.csv
## OUTPUT: new columns site, is_tumour, cell_type_organ, cell_type_organ_tumour in
##         xenium_<seg>_annotated.rds (re-run 04b after re-running 04),
##         annotation/04b_<seg>_cell_types_per_organ.csv.gz, tables/04b_*.csv,
##         figures/04b_*.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/annotation_helpers.R")

ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))
ORGAN_SITES     <- get0("ORGAN_SITES", ifnotfound = NULL)                        # NULL = every site in the csv
ORGAN_LINEAGES  <- get0("ORGAN_LINEAGES", ifnotfound = c("Epithelial", "Neural_NE"))
RES_ORGAN       <- get0("RES_ORGAN", ifnotfound = 0.6)
MIN_CELLS_ORGAN <- get0("MIN_CELLS_ORGAN", ifnotfound = 300)
MODES <- trimws(strsplit(Sys.getenv("EXP1_04B_MODES", "normal,with_tumour"), ",")[[1]])

xen <- readRDS(obj_path("annotated"))
if (!all(c("lineage", "cell_type") %in% colnames(xen[[]])))
  stop("Run 04_annotation_hierarchical.R to the end first (lineage / cell_type missing).")
DefaultAssay(xen) <- "Xenium"
panel_genes <- rownames(xen[["Xenium"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)
if (!"data" %in% Layers(xen[["Xenium"]]))                 # 04 saves the object normalised
  xen <- NormalizeData(xen, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)

## -----------------------------------------------------------------------
## 0. site and tumour status per cell
## -----------------------------------------------------------------------
xen <- add_site_and_tumour(xen, ANNOT_FAMILIES_CSV)
sites <- pick_sites(xen, ORGAN_SITES)
in_lin <- xen$lineage %in% ORGAN_LINEAGES
site_tab <- xen[[]][in_lin & xen$site %in% sites, ] |>
  count(site, tumour = ifelse(is_tumour, "tumour", "normal"), organ, name = "n_cells") |>
  arrange(site, tumour, organ)
save_table(site_tab, "04b_cells_per_site")
message(paste(ORGAN_LINEAGES, collapse = " + "), " cells per site:")
print(as.data.frame(site_tab |> group_by(site, tumour) |> summarise(n_cells = sum(n_cells), .groups = "drop")),
      row.names = FALSE)

## signatures: level-2 markers of the families present (normal + tumour), + Proliferating
site_sigs <- function(cells) {
  fams <- unique(xen$family[match(cells, colnames(xen))])
  s <- do.call(c, unname(markers$level2[intersect(fams, names(markers$level2))]))
  s <- c(s, markers$level2[["Immune"]]["Proliferating"])
  s[!duplicated(names(s))]
}

## -----------------------------------------------------------------------
## 1. normal / normal + tumour, per site (R/annotation_helpers.R: annotate_sites)
## -----------------------------------------------------------------------
res <- annotate_sites(xen, sites, ORGAN_LINEAGES, site_sigs, RES_ORGAN, prefix = "04b", modes = MODES,
                      min_cells = MIN_CELLS_ORGAN)
xen <- set_labels(xen, "cell_type_organ", res$normal)
xen <- set_labels(xen, "cell_type_organ_tumour", res$with_tumour)
if (length(res$failed)) {
  message("\n!! Groups that failed (send this to the analyst):")
  for (g in names(res$failed)) message("   ", g, ": ", res$failed[[g]])
}

## -----------------------------------------------------------------------
## 2. tables + save
## -----------------------------------------------------------------------
m <- xen[[]]
done <- !is.na(m$cell_type_organ) | !is.na(m$cell_type_organ_tumour)
tab <- m[done, ] |> count(site, tumour = ifelse(is_tumour, "tumour", "normal"),
                          cell_type_organ, cell_type_organ_tumour, name = "n_cells")
save_table(tab, "04b_cell_types_per_organ")
out <- data.frame(cell = rownames(m)[done], m[done, intersect(c("slide_id", "core_id", "organ", "site", "is_tumour",
                                                                  "lineage", "cell_type", "cell_type_organ",
                                                                  "cell_type_organ_tumour"), colnames(m))])
f <- file.path(ANNOT_DIR, sprintf("04b_%s_cell_types_per_organ.csv.gz", SEGMENTATION))
utils::write.csv(out, gzfile(f), row.names = FALSE)
message("Per-cell table -> ", f)

saveRDS(xen, obj_path("annotated"))
message("Saved (columns site, is_tumour, cell_type_organ, cell_type_organ_tumour added) -> ", obj_path("annotated"))
