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
key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
fam <- read_config_csv(ANNOT_FAMILIES_CSV, stringsAsFactors = FALSE, fileEncoding = "latin1")
if (!"site" %in% colnames(fam)) stop(basename(ANNOT_FAMILIES_CSV), " has no `site` column -- update config/.")
xen$site <- fam$site[match(key(xen$organ), key(fam$organ))]
xen$site[!is.na(xen$site) & !nzchar(trimws(xen$site))] <- NA
## tumour = the cells of the tumour TMA (slide level, see tumour_slides() in R/helpers.R)
tum_slides <- tumour_slides()
xen$is_tumour <- xen$slide_id %in% tum_slides
message("Tumour slide(s): ", paste(tum_slides, collapse = ", "), "  (from config/slides.csv)")
print(table(slide = xen$slide_id, tumour = xen$is_tumour))

sites <- sort(unique(stats::na.omit(xen$site)))
if (!is.null(ORGAN_SITES)) {
  bad <- setdiff(ORGAN_SITES, sites)
  if (length(bad)) message("ORGAN_SITES not found in the data: ", paste(bad, collapse = ", "))
  sites <- intersect(sites, ORGAN_SITES)
}
if (!length(sites)) stop("No sites to run -- fill the `site` column of ", basename(ANNOT_FAMILIES_CSV))
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

## UMAP of one site coloured by core (+ by normal / tumour and by label): do clusters
## follow cell types (mixed cores) or single cores / patients (batch, or a tumour)?
plot_umap_by_core <- function(emb, meta, lab, g, max_cells = 2e5) {
  d <- data.frame(emb[, 1:2], core = meta[rownames(emb), "core_id"], organ = meta[rownames(emb), "organ"],
                  tissue = ifelse(meta[rownames(emb), "is_tumour"], "tumour", "normal"),
                  label = unname(lab[rownames(emb)]))
  colnames(d)[1:2] <- c("u1", "u2")
  d <- d[!is.na(d$core), ]
  if (nrow(d) > max_cells) d <- d[sample(nrow(d), max_cells), ]
  d <- d[sample(nrow(d)), ]                                   # no core drawn on top of all others
  d$core_lab <- paste0(d$core, " (", d$organ, ")")
  n_core <- length(unique(d$core_lab))
  pal <- stats::setNames(grDevices::hcl(h = seq(15, 375, length.out = n_core + 1)[seq_len(n_core)],
                                        c = rep(c(85, 55), length.out = n_core),
                                        l = rep(c(55, 72, 42), length.out = n_core)),
                         sort(unique(d$core_lab)))
  pt <- if (nrow(d) > 1e5) 0.25 else if (nrow(d) > 2e4) 0.5 else if (nrow(d) > 2000) 1 else 2
  base <- list(theme_void(base_size = 12), coord_fixed(), theme(plot.title = element_text(size = 13)))
  p_core <- ggplot(d, aes(u1, u2, colour = core_lab)) + geom_point(size = pt, stroke = 0) +
    scale_colour_manual(values = pal, name = "core (organ)") + base +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = if (n_core > 25) 2 else 1)) +
    ggtitle(paste0(g, ": coloured by core (", n_core, " cores)"))
  plots <- list(p_core)
  if (length(unique(d$tissue)) > 1)
    plots <- c(plots, list(ggplot(d, aes(u1, u2, colour = tissue)) + geom_point(size = pt, stroke = 0) +
      scale_colour_manual(values = c(normal = "#3B6FB6", tumour = "#C0392B"), name = NULL) + base +
      guides(colour = guide_legend(override.aes = list(size = 3))) + ggtitle("normal vs tumour cores")))
  pal_l <- celltype_palette(d$label, rep("Epithelial", nrow(d)))
  plots <- c(plots, list(ggplot(d, aes(u1, u2, colour = label)) + geom_point(size = pt, stroke = 0) +
    scale_colour_manual(values = pal_l, name = NULL) + base +
    guides(colour = guide_legend(override.aes = list(size = 3))) + ggtitle("current labels")))
  save_fig(patchwork::wrap_plots(plots, nrow = 1), paste0("04b_", g, "_umap_by_core"),
           width = 9.5 * length(plots) + 2, height = 9, dpi = 160)
}

## -----------------------------------------------------------------------
## 1. normal / normal + tumour, per site
## -----------------------------------------------------------------------
xen$cell_type_organ <- NA_character_
xen$cell_type_organ_tumour <- NA_character_
failed <- list()
for (mode in MODES) {
  out_col <- if (mode == "normal") "cell_type_organ" else "cell_type_organ_tumour"
  for (s in sites) {
    sel <- in_lin & xen$site %in% s & (mode != "normal" | !xen$is_tumour)
    cells <- colnames(xen)[sel]
    g <- paste0(if (mode == "normal") "Organ_" else "OrganTumour_", s)
    if (mode != "normal" && !any(xen$is_tumour[sel])) {
      message(g, ": no tumour cores for this site -- skipped (same as Organ_", s, ")")
      next
    }
    if (mode == "normal" && !length(cells) && any(xen$is_tumour[in_lin & xen$site %in% s])) {
      message(g, ": tumour cores only (no normal tissue of this site on the TMA) -- see OrganTumour_", s)
      next
    }
    if (length(cells) < MIN_CELLS_ORGAN) {
      message(g, ": only ", length(cells), " cells (< MIN_CELLS_ORGAN = ", MIN_CELLS_ORGAN, ") -- skipped")
      next
    }
    sub <- subset(xen, cells = cells)
    ann <- tryCatch(annotate_group(sub, g, site_sigs(cells), RES_ORGAN, comp_by = "organ",
                                   prefix = "04b", harmony = FALSE),
                    error = function(e) {
                      failed[[g]] <<- conditionMessage(e)
                      message("!! ", g, " FAILED: ", conditionMessage(e), " -- the other sites continue.")
                      NULL
                    })
    lab <- ann$labels
    if (!is.null(ann)) plot_umap_by_core(ann$umap, sub[[]], lab, g)
    if (!is.null(lab)) {
      v <- xen[[out_col, drop = TRUE]]
      v[match(names(lab), colnames(xen))] <- lab
      xen[[out_col]] <- v
      if (mode != "normal") {
        ## which labels make up each core, normal vs tumour (same cell types, side by side)
        m <- xen[[]][names(lab), ]
        m$label <- lab
        comp <- m |> filter(!is.na(core_id)) |> count(is_tumour, organ, core_id, label) |>
          group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
        comp$panel <- ifelse(comp$is_tumour, "tumour", "normal")
        pal <- celltype_palette(comp$label, rep("Epithelial", nrow(comp)))
        p <- ggplot(comp, aes(core_id, frac, fill = factor(label, levels = names(pal)))) +
          geom_col(width = 0.9) +
          facet_grid(~ panel + organ, scales = "free_x", space = "free_x", labeller = label_wrap_gen(18)) +
          scale_fill_manual(values = pal, name = "label") + scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 11) +
          theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5), panel.grid = element_blank(),
                strip.text = element_text(size = 8)) +
          labs(x = "core", y = paste("fraction of", paste(ORGAN_LINEAGES, collapse = " + "), "cells"),
               title = paste0(s, ": cell types per core, normal vs tumour (one joint clustering)"))
        save_fig(p, paste0("04b_", g, "_composition_per_core"),
                 width = max(10, 0.25 * length(unique(comp$core_id)) + 6), height = 6.5)
      }
    }
    rm(sub); invisible(gc())
  }
}
if (length(failed)) {
  message("\n!! Groups that failed (send this to the analyst):")
  for (g in names(failed)) message("   ", g, ": ", failed[[g]])
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
