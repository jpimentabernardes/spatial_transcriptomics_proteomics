## =============================================================================
## 05c_checkpoint_genes_by_compartment.R
##
## PDCD1 (PD-1), CD274 (PD-L1), NCAM1 (CD56), CD8A and CD8B per cell, with
## coordinates, core, organ and tissue compartment, for TMA_ORGAN and TMA_TUMOR.
##
## Compartment per cell from its lineage + cell type label, rules in
## config/compartments.csv (first matching row wins; edit freely):
##   Epithelial/Tumour   epithelial + neuroendocrine cells (tumour cells on TMA_TUMOR)
##   T cells             T-cell labels (T_cell, CD4_T, CD8_T, Treg, Tfh, gd_T, ...)
##   Macrophage/DC       macrophages, monocytes, dendritic cells
##   NK/ILC, B/Plasma, Other immune, Stroma
## Label used: cell_type_organ_all (04c, organ-level) when present, else cell_type (04);
## set CHECKPOINT_LABEL_COL in config.R to force one.
##
## A cell is "positive" for a gene at >= CHECKPOINT_MIN_COUNTS transcripts (default 1).
## One transcript can be background or spill-over from a neighbour (see 03c MECR / 03f);
## set 2 for a stricter call -- the counts are in the table, so you can re-threshold.
##
## INPUT : xenium_<seg>_annotated.rds (04, ideally after 04c), config/compartments.csv
## OUTPUT: data/export/05c_punch_file_<NTA|MTA>_<seg>.csv     ONE ROW PER CORE (patient): per marker
##           (PD1, PDL1, CD56, CD8 = CD8A + CD8B, CD8A, CD8B) transcripts and transcripts per cell,
##           for all cells and per compartment (EpiTumour, MacDC, T) -- section 3b
##         data/export/05c_checkpoint_cells_<slide>.csv.gz    one row per cell (in a core):
##           cell, slide_id, is_tumour, core_id, organ, site, x_um, y_um, lineage, label,
##           compartment, nCount_Xenium, <gene> (counts), <gene>_norm (log-normalised)
##         tables/05c_checkpoint_by_core_compartment.csv   per core x compartment:
##           n_cells, pct positive, mean counts, mean normalised expression per gene
##         tables/05c_checkpoint_by_site_compartment.csv   per site x normal/tumour x compartment
##         figures/05c_dotplot_compartments.png        genes x compartment, TMA_ORGAN vs TMA_TUMOR
##         figures/05c_dotplot_by_site.png             genes x organ site, per compartment
##         figures/05c_per_core_fractions.png          % PD-1+ T, % CD8+ T, % PD-L1+ macrophage/DC,
##                                                     % PD-L1+ epithelium/tumour, % CD56+ ... per core
##         figures/05c_spatial_<slide>.png             where the positive cells are, per gene
##         figures/05c_example_cores_<slide>.png       6 cores close up: compartments + PD-1 / PD-L1 cells
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/annotation_helpers.R")

CHECKPOINT_GENES <- get0("CHECKPOINT_GENES", ifnotfound = c(PDCD1 = "PDCD1", CD274 = "CD274", CD56 = "NCAM1",
                                                            CD8A = "CD8A", CD8B = "CD8B"))
## extra genes: EXP1_05C_GENES=CX3CR1,IFNG,IL6 (or CHECKPOINT_EXTRA_GENES in config.R), added to the
## list above (genes already in it are skipped); each gets its own columns in the punch file
extra <- trimws(unlist(strsplit(Sys.getenv("EXP1_05C_GENES", paste(get0("CHECKPOINT_EXTRA_GENES", ifnotfound = ""),
                                                                   collapse = ",")), ",")))
extra <- unique(toupper(extra[nzchar(extra)]))
extra <- extra[!extra %in% toupper(CHECKPOINT_GENES)]
if (length(extra)) CHECKPOINT_GENES <- c(CHECKPOINT_GENES, stats::setNames(extra, extra))
CHECKPOINT_MIN_COUNTS <- get0("CHECKPOINT_MIN_COUNTS", ifnotfound = 1)
CHECKPOINT_LABEL_COL  <- get0("CHECKPOINT_LABEL_COL", ifnotfound = NULL)
COMPARTMENTS_CSV <- get0("COMPARTMENTS_CSV", ifnotfound = file.path(CONFIG_DIR, "compartments.csv"))
ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))
MIN_CELLS_CORE <- 20                       # per-core fractions need this many cells of the compartment
EXPORT_DIR <- file.path(PROJECT_DIR, "data", "export")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

COMP_LEVELS <- c("Epithelial/Tumour", "T cells", "Macrophage/DC", "NK/ILC", "B/Plasma", "Other immune", "Stroma", "Other")
COMP_COLOURS <- c("Epithelial/Tumour" = "#418ED0", "T cells" = "#C0392B", "Macrophage/DC" = "#E69F00",
                  "NK/ILC" = "#8E44AD", "B/Plasma" = "#16A085", "Other immune" = "#C17743",
                  "Stroma" = "#7F8C8D", "Other" = "#BDBDBD")

xen <- readRDS(obj_path("annotated"))
if (!all(c("site", "is_tumour") %in% colnames(xen[[]]))) xen <- add_site_and_tumour(xen, ANNOT_FAMILIES_CSV)
genes <- match_case(CHECKPOINT_GENES, rownames(xen[["Xenium"]]))
names(genes) <- names(CHECKPOINT_GENES)
miss <- names(genes)[!genes %in% rownames(xen[["Xenium"]])]
if (length(miss)) message("Not on the panel (left out): ", paste(sprintf("%s (%s)", miss, genes[miss]), collapse = ", "))
genes <- genes[!names(genes) %in% miss]
if (!length(genes)) stop("None of the checkpoint genes is on the panel.")
message("Genes in this run (", length(genes), "): ", paste(names(genes), collapse = ", "),
        if (!length(extra)) "  -- no extra genes (set EXP1_05C_GENES=IFNG,IL6,...)" else "")

label_col <- CHECKPOINT_LABEL_COL %||%
  (if ("cell_type_organ_all" %in% colnames(xen[[]])) "cell_type_organ_all" else "cell_type")
message("Cell type labels from: ", label_col,
        if (label_col == "cell_type") "  (run 04c for organ-level immune / stroma labels)" else "")

## -----------------------------------------------------------------------
## 1. compartment per cell
## -----------------------------------------------------------------------
rules <- read_config_csv(COMPARTMENTS_CSV, stringsAsFactors = FALSE)
m <- xen[[]]
m$label <- as.character(m[[label_col]])
m$label[is.na(m$label)] <- as.character(m$cell_type[is.na(m$label)])
m$compartment <- NA_character_
for (k in seq_len(nrow(rules))) {
  hit <- is.na(m$compartment) & m$lineage %in% rules$lineage[k] &
    grepl(rules$pattern[k], ifelse(is.na(m$label), "", m$label), ignore.case = TRUE)
  m$compartment[hit] <- rules$compartment[k]
}
m$compartment[is.na(m$compartment)] <- "Other"
m$compartment <- factor(m$compartment, levels = unique(c(intersect(COMP_LEVELS, m$compartment),
                                                         setdiff(unique(m$compartment), COMP_LEVELS))))
for (cp in setdiff(levels(m$compartment), names(COMP_COLOURS)))
  COMP_COLOURS[cp] <- grDevices::hcl(h = 360 * stats::runif(1), c = 60, l = 55)
comp_tab <- m |> count(slide_id, compartment, label, name = "n_cells") |> arrange(slide_id, compartment, desc(n_cells))
save_table(comp_tab, "05c_labels_per_compartment")
message("Compartments (check that every label landed where you expect: tables/05c_labels_per_compartment.csv):")
print(table(compartment = m$compartment, slide = m$slide_id))

## -----------------------------------------------------------------------
## 2. per-cell export (cells in cores)
## -----------------------------------------------------------------------
in_core <- !is.na(m$core_id)
cnt <- GetAssayData(xen, assay = "Xenium", layer = "counts")[genes, , drop = FALSE]
nrm <- if ("data" %in% Layers(xen[["Xenium"]])) GetAssayData(xen, assay = "Xenium", layer = "data")[genes, , drop = FALSE] else
  log1p(sweep(cnt, 2, pmax(m$nCount_Xenium, 1), "/") * stats::median(m$nCount_Xenium))
cnt <- as.matrix(Matrix::t(cnt)); colnames(cnt) <- names(genes)
nrm <- round(as.matrix(Matrix::t(nrm)), 4); colnames(nrm) <- paste0(names(genes), "_norm")
if (!"patient_id" %in% colnames(m)) m$patient_id <- NA
cells <- data.frame(cell = rownames(m), slide_id = m$slide_id, is_tumour = m$is_tumour, core_id = m$core_id,
                    patient_id = m$patient_id, organ = m$organ, site = m$site, x_um = round(m$x_um, 2), y_um = round(m$y_um, 2),
                    lineage = m$lineage, label = m$label, compartment = m$compartment,
                    nCount_Xenium = m$nCount_Xenium, cnt, nrm, check.names = FALSE)[in_core, ]
for (sid in unique(cells$slide_id)) {
  f <- file.path(EXPORT_DIR, sprintf("05c_checkpoint_cells_%s_%s.csv.gz", SEGMENTATION, sid))
  d <- cells[cells$slide_id == sid, ]
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(d, f) else
    utils::write.csv(d, gzfile(f), row.names = FALSE)
  message("Exported ", format(nrow(d), big.mark = ","), " cells -> ", f)
}
gene_names <- names(genes)
pos <- cells[, gene_names, drop = FALSE] >= CHECKPOINT_MIN_COUNTS

## -----------------------------------------------------------------------
## 3. summaries per core x compartment and per site x compartment
## -----------------------------------------------------------------------
summ_by <- function(by, keep = rep(TRUE, nrow(cells))) {
  d <- cells[keep, ]; pz <- pos[keep, , drop = FALSE]
  grp <- interaction(lapply(d[, by, drop = FALSE], function(v) ifelse(is.na(v), "NA", as.character(v))),
                     drop = TRUE, sep = "\r")
  idx <- split(seq_len(nrow(d)), grp)
  dplyr::bind_rows(lapply(idx, function(i) {
    row <- d[i[1], by, drop = FALSE]
    row$n_cells <- length(i)
    for (g in gene_names) {
      row[[paste0("pct_pos_", g)]] <- round(100 * mean(pz[i, g]), 2)
      row[[paste0("mean_counts_", g)]] <- round(mean(d[[g]][i]), 4)
      row[[paste0("mean_norm_", g)]] <- round(mean(d[[paste0(g, "_norm")]][i]), 4)
    }
    row
  }))
}
by_core <- summ_by(c("slide_id", "is_tumour", "site", "organ", "core_id", "compartment"))
save_table(by_core, "05c_checkpoint_by_core_compartment")
by_site <- summ_by(c("slide_id", "is_tumour", "site", "compartment"))
save_table(by_site, "05c_checkpoint_by_site_compartment")

## -----------------------------------------------------------------------
## 3b. "Punch file" for the comparison with IHC / RNA-FISH: ONE ROW PER CORE.
##     For each marker (PD1 = PDCD1, PDL1 = CD274, CD56 = NCAM1, CD8 = CD8A + CD8B,
##     and CD8A / CD8B separately):
##       <M>_transcripts            transcripts of the marker in all cells of the core
##       <M>_per_cell               ... divided by the number of cells in the core
##       <M>_transcripts_<comp>     transcripts in the cells of one compartment
##       <M>_per_cell_<comp>        ... divided by the number of cells of that compartment
##     compartments: EpiTumour = Epithelial/Tumour, MacDC = Macrophage/DC, T = T cells
##     (from the cell type annotation, see config/compartments.csv), plus n_cells_<comp>.
##     Cells = QC-passed cells whose centre lies in the core; transcripts that no cell
##     was assigned to are not counted.
## -----------------------------------------------------------------------
PUNCH_COMPARTMENTS <- get0("PUNCH_COMPARTMENTS", ifnotfound = c(EpiTumour = "Epithelial/Tumour",
                                                              MacDC = "Macrophage/DC", T = "T cells"))
pm <- cells[, gene_names, drop = FALSE]
short <- c(PDCD1 = "PD1", CD274 = "PDL1", CD56 = "CD56", CD8A = "CD8A", CD8B = "CD8B")[gene_names]
names(pm) <- ifelse(is.na(short), gene_names, short)
if (all(c("CD8A", "CD8B") %in% names(pm))) pm <- cbind(pm[, setdiff(names(pm), c("CD8A", "CD8B")), drop = FALSE],
                                                    CD8 = pm$CD8A + pm$CD8B, pm[, c("CD8A", "CD8B")])
pm <- as.matrix(pm)
key <- paste(cells$slide_id, cells$core_id, sep = "\r")
core_tab <- cells[!duplicated(key), c("slide_id", "is_tumour", "core_id", "patient_id", "organ", "site")]
core_tab$tma <- ifelse(core_tab$is_tumour, "MTA", "NTA")
rownames(core_tab) <- key[!duplicated(key)]
add_block <- function(tab, sel, suffix) {
  n <- as.numeric(table(factor(key[sel], levels = rownames(tab)))[rownames(tab)])
  s <- matrix(0, nrow(tab), ncol(pm), dimnames = list(rownames(tab), colnames(pm)))
  if (any(sel)) {                                  # rowsum() drops cores without cells of the compartment
    r <- rowsum(pm[sel, , drop = FALSE], key[sel])
    s[rownames(r), ] <- r
  }
  tab[[paste0("n_cells", suffix)]] <- as.integer(n)
  if (suffix == "") tab$total_transcripts <- as.numeric(rowsum(cells$nCount_Xenium, factor(key, levels = rownames(tab)))[rownames(tab), 1])
  for (mk in colnames(pm)) {
    tab[[paste0(mk, "_transcripts", suffix)]] <- s[, mk]
    tab[[paste0(mk, "_per_cell", suffix)]] <- round(ifelse(n > 0, s[, mk] / n, NA), 5)
  }
  tab
}
punch <- add_block(core_tab, rep(TRUE, nrow(cells)), "")
for (cn in names(PUNCH_COMPARTMENTS)) punch <- add_block(punch, cells$compartment %in% PUNCH_COMPARTMENTS[[cn]], paste0("_", cn))
## matching keys to the pathology grids, where 02 kept them
cores_file <- file.path(TABLES_DIR, "02_cores_all.csv")
if (file.exists(cores_file)) {
  ca <- utils::read.csv(cores_file, stringsAsFactors = FALSE)
  for (col in intersect(c("core_name", "cellscape_core"), colnames(ca)))
    punch[[col]] <- ca[[col]][match(paste(punch$slide_id, punch$core_id), paste(ca$slide_id, ca$core_id))]
}
front <- intersect(c("tma", "slide_id", "core_id", "core_name", "cellscape_core", "patient_id", "organ", "site"),
                   colnames(punch))
punch <- punch[order(punch$slide_id, punch$core_id), c(front, setdiff(colnames(punch), c(front, "is_tumour")))]
rownames(punch) <- NULL
save_table(punch, "05c_punch_file_per_core")
for (t in unique(punch$tma)) {
  f <- file.path(EXPORT_DIR, sprintf("05c_punch_file_%s_%s.csv", t, SEGMENTATION))
  utils::write.csv(punch[punch$tma == t, ], f, row.names = FALSE)
  message("Punch file (", t, ", ", sum(punch$tma == t), " cores) -> ", f)
}

## -----------------------------------------------------------------------
## 4. figures
## -----------------------------------------------------------------------
main_comp <- intersect(c("Epithelial/Tumour", "T cells", "Macrophage/DC", "NK/ILC"), levels(cells$compartment))
long_dot <- function(keep, by) {
  dd <- summ_by(by, keep)
  tidyr::pivot_longer(dd, dplyr::starts_with(c("pct_pos_", "mean_norm_")),
                      names_to = c(".value", "gene"), names_pattern = "(pct_pos|mean_norm)_(.*)")
}
## 4a. dot plot genes x compartment, per slide
dd <- long_dot(rep(TRUE, nrow(cells)), c("slide_id", "compartment"))
dd$gene <- factor(dd$gene, levels = gene_names)
p1 <- ggplot(dd, aes(gene, compartment, size = pct_pos, colour = mean_norm)) + geom_point() +
  scale_size_area(max_size = 11, name = sprintf("%% cells\n>= %d count", CHECKPOINT_MIN_COUNTS)) +
  scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\nnormalised") +
  scale_y_discrete(limits = rev) + facet_wrap(~ slide_id) + theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, title = "Checkpoint / cytotoxicity genes per compartment",
       subtitle = paste("label column:", label_col))
save_fig(p1, "05c_dotplot_compartments", width = 12, height = 6)

## 4b. dot plot per site (main compartments), normal vs tumour
ds <- long_dot(cells$compartment %in% main_comp & !is.na(cells$site),
               c("site", "is_tumour", "compartment")) |> filter(n_cells >= 50)
if (nrow(ds)) {
  ds$gene <- factor(ds$gene, levels = gene_names)
  ds$tissue <- ifelse(ds$is_tumour, "tumour", "normal")
  ds$compartment <- factor(ds$compartment, levels = main_comp)
  p2 <- ggplot(ds, aes(gene, site, size = pct_pos, colour = mean_norm)) + geom_point() +
    scale_size_area(max_size = 7, name = sprintf("%% cells\n>= %d count", CHECKPOINT_MIN_COUNTS)) +
    scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\nnormalised") +
    scale_y_discrete(limits = rev) + facet_grid(tissue ~ compartment, scales = "free_y", space = "free_y") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, title = "Per organ site: normal (TMA_ORGAN) vs tumour (TMA_TUMOR)",
         subtitle = "site x compartment with >= 50 cells")
  save_fig(p2, "05c_dotplot_by_site", width = 15, height = max(7, 0.22 * length(unique(ds$site)) * 2 + 3))
}

## 4c. per-core fractions of the key populations
frac_def <- list(
  "% PD-1+ of T cells"                 = list(comp = "T cells", genes = "PDCD1"),
  "% CD8A/B+ of T cells"               = list(comp = "T cells", genes = c("CD8A", "CD8B")),
  "% PD-L1+ of macrophage/DC"          = list(comp = "Macrophage/DC", genes = "CD274"),
  "% PD-L1+ of epithelium/tumour"      = list(comp = "Epithelial/Tumour", genes = "CD274"),
  "% CD56+ of NK/ILC"                  = list(comp = "NK/ILC", genes = "CD56"),
  "% PD-1+ of CD8A/B+ T cells"         = list(comp = "T cells", genes = "PDCD1", within = c("CD8A", "CD8B")))
fr <- dplyr::bind_rows(lapply(names(frac_def), function(nm) {
  f <- frac_def[[nm]]
  if (!all(c(f$genes, f$within) %in% gene_names)) return(NULL)
  sel <- cells$compartment == f$comp & !is.na(cells$site)
  if (!is.null(f$within)) sel <- sel & rowSums(pos[, f$within, drop = FALSE]) > 0
  if (!any(sel)) return(NULL)
  d <- cells[sel, c("slide_id", "is_tumour", "site", "core_id")]
  d$pos <- rowSums(pos[sel, f$genes, drop = FALSE]) > 0
  d |> group_by(slide_id, is_tumour, site, core_id) |>
    summarise(n = n(), pct = 100 * mean(pos), .groups = "drop") |>
    filter(n >= MIN_CELLS_CORE) |> mutate(metric = nm)
}))
if (nrow(fr)) {
  save_table(fr, "05c_per_core_fractions")
  fr$metric <- factor(fr$metric, levels = names(frac_def))
  fr$tissue <- ifelse(fr$is_tumour, "tumour (TMA_TUMOR)", "normal (TMA_ORGAN)")
  site_order <- fr |> group_by(site) |> summarise(n = n_distinct(core_id)) |> arrange(site) |> pull(site)
  p3 <- ggplot(fr, aes(factor(site, levels = site_order), pct, colour = tissue)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.9, alpha = 0.7) +
    scale_colour_manual(values = c("normal (TMA_ORGAN)" = "#3B6FB6", "tumour (TMA_TUMOR)" = "#C0392B"), name = NULL) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    labs(x = NULL, y = "% of the compartment's cells in the core",
         title = "Per core (one dot = one core with >= 20 cells of that compartment)",
         subtitle = sprintf("positive = >= %d transcript(s)", CHECKPOINT_MIN_COUNTS))
  save_fig(p3, "05c_per_core_fractions", width = max(12, 0.45 * length(site_order) + 4),
           height = 3.2 * length(unique(fr$metric)) + 1.5)
}

## 4d. whole slide: where are the positive cells? one panel per gene
for (sid in unique(cells$slide_id)) {
  d <- cells[cells$slide_id == sid, ]
  bg <- d[sample(nrow(d), min(nrow(d), 3e5)), c("x_um", "y_um")]
  pp <- lapply(gene_names, function(g) {
    p <- d[pos[cells$slide_id == sid, g], ]
    ggplot() + geom_point(data = bg, aes(x_um, -y_um), size = 0.05, colour = "grey85", stroke = 0) +
      geom_point(data = p, aes(x_um, -y_um, colour = compartment), size = 0.25, stroke = 0) +
      scale_colour_manual(values = COMP_COLOURS, drop = FALSE, name = NULL) + coord_fixed() + theme_void(base_size = 12) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      ggtitle(sprintf("%s (%s): %s positive cells", g, genes[[g]], format(nrow(p), big.mark = ",")))
  })
  save_fig(patchwork::wrap_plots(pp, ncol = 3, guides = "collect") +
             patchwork::plot_annotation(title = paste0(sid, ": cells with >= ", CHECKPOINT_MIN_COUNTS,
                                                       " transcript(s), coloured by compartment")),
           paste0("05c_spatial_", sid), width = 24, height = 8 * ceiling(length(gene_names) / 3) + 1, dpi = 200)
}

## 4e. six cores close up: all cells by compartment, PD-1+ (red ring) and PD-L1+ (black ring) cells
for (sid in unique(cells$slide_id)) {
  d <- cells[cells$slide_id == sid, ]
  pd <- pos[cells$slide_id == sid, , drop = FALSE]
  sc <- d |> mutate(pdl1 = if ("CD274" %in% gene_names) pd[, "CD274"] else FALSE,
                    pd1 = if ("PDCD1" %in% gene_names) pd[, "PDCD1"] else FALSE) |>
    group_by(core_id, organ) |> summarise(n = n(), score = mean(pdl1) + mean(pd1), .groups = "drop") |>
    filter(n >= min(500, max(n) / 2)) |> arrange(desc(score))
  top <- head(sc$core_id, 6)
  if (!length(top)) next
  pp <- lapply(top, function(cid) {
    dc <- d[d$core_id == cid, ]
    pc <- pd[d$core_id == cid, , drop = FALSE]
    sz <- if (nrow(dc) < 2000) 1.6 else if (nrow(dc) < 10000) 0.9 else 0.5
    p <- ggplot(dc, aes(x_um, -y_um)) + geom_point(aes(colour = compartment), size = sz, stroke = 0) +
      scale_colour_manual(values = COMP_COLOURS, drop = FALSE, name = NULL)
    if ("PDCD1" %in% gene_names)
      p <- p + geom_point(data = dc[pc[, "PDCD1"], ], shape = 21, size = 2.5 * sz + 0.6, colour = "red", stroke = 0.4)
    if ("CD274" %in% gene_names)
      p <- p + geom_point(data = dc[pc[, "CD274"], ], shape = 21, size = 3.2 * sz + 0.8, colour = "black", stroke = 0.3)
    p + coord_fixed() + theme_void(base_size = 11) + guides(colour = guide_legend(override.aes = list(size = 3))) +
      ggtitle(sprintf("%s  %s", cid, sc$organ[match(cid, sc$core_id)]))
  })
  save_fig(patchwork::wrap_plots(pp, ncol = 3, guides = "collect") +
             patchwork::plot_annotation(title = paste0(sid, ": cores with the most PD-1+ / PD-L1+ cells"),
                                        subtitle = "colour = compartment; red ring = PD-1 (PDCD1)+, black ring = PD-L1 (CD274)+"),
           paste0("05c_example_cores_", sid), width = 21, height = 14, dpi = 200)
}
message("Done.")
