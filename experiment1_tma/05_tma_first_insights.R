## =============================================================================
## 05_tma_first_insights.R
##
## First insights from the multi-organ / multi-tumour TMAs, organised by aim.
##
## AIM 1 -- probe specificity: are the colon + custom probes gut-specific, or
##          do the immune targets generalise to other organs?
##   1.1 Signal over background per gene x core (vs negative-control probes).
##       "Detected" = clearly above the non-binding-probe rate in that core.
##       (segger has no control probes -> the Xenium-segmentation object of
##       the same cores is used for this step.)
##   1.1b GI specificity per probe: detected in GI glandular cores (colon,
##       rectum, appendix, small intestine, stomach) vs non-GI cores, and
##       GI/non-GI expression on all cells vs within epithelium / immune cells.
##   1.2 Family specificity (tau) -- computed twice: on all cells AND within
##       immune cells only. An immune probe with high tau on all cells but low
##       tau within immune cells is NOT organ-specific; organs just differ in
##       how many immune cells they contain. Separating these two is the key
##       step for answering Aim 1 honestly.
##   1.3 Lineage specificity: does an immune probe's signal land in immune
##       cells in every organ, or does it appear in epithelium/stroma in some
##       organs (non-specific binding, or segmentation spill-over -> compare
##       with MECR in 03c)?
##   1.4 Generalisability: is the immune-target profile of e.g. T cells in
##       lung the same as in GI cores (correlation per cell type and family)?
##
## AIM 2 -- heterogeneity and comparability between organs, locations, patients
##   2.1 Composition per core + entropy + differential composition (propeller)
##   2.2 Variance partitioning of pseudobulk expression: how much variance is
##       organ vs patient vs slide (batch) vs cell type vs residual?
##   2.3 Core-to-core similarity: do cores cluster by organ, patient or slide?
##       Same-patient vs different-patient correlation within an organ.
##   2.4 Location within organ (e.g. proximal/distal, tumour centre/margin):
##       pseudobulk limma-voom, patient as blocking factor.
##   2.5 Cellular niches (neighbourhood composition) per organ.
##
## Statistical unit = CORE (or patient), never cell. Millions of cells give
## tiny p-values for trivial differences; cores/patients are the replicates.
##
## Organs are grouped by the families of 04 (config/annotation_families.csv).
## Quick run, probe specificity only:  EXP1_05_AIMS=1 Rscript 05_tma_first_insights.R
##
## INPUT : xenium_<seg>_annotated.rds (04), config/probe_expectations.csv
##         (panel = colon / custom, category, expected_lineage per gene)
## OUTPUT: tables/05_*.csv, figures/05_*.png, xenium_<seg>_insights.rds
## =============================================================================

source("config/config.R")
source("R/helpers.R")

## which aims to run: EXP1_05_AIMS=1 (probe specificity only, quick), 2, or 1,2 (default)
AIMS <- trimws(strsplit(Sys.getenv("EXP1_05_AIMS", "1,2"), ",")[[1]])
## families counted as "GI" (colon, rectum, appendix, small intestine, stomach)
GI_FAMILIES <- get0("GI_FAMILIES", ifnotfound = "GI_glandular")
ANNOT_FAMILIES_CSV <- get0("ANNOT_FAMILIES_CSV", ifnotfound = file.path(CONFIG_DIR, "annotation_families.csv"))

xen <- readRDS(obj_path("annotated"))
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")
meta <- xen[[]]
genes <- rownames(counts)

probes <- read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE)
probes$gene <- match_case(as_seurat_features(probes$gene), genes)
missing_probes <- setdiff(probes$gene, genes)
if (length(missing_probes) > 0) message("In probe_expectations.csv but not on panel: ", paste(missing_probes, collapse = ", "))
probes <- probes[probes$gene %in% genes, ]
for (v in c("expected_lineage", "category")) probes[[v]][is.na(probes[[v]])] <- ""
rownames(probes) <- probes$gene
custom_genes <- probes$gene[probes$panel == "custom"]
## immune targets of the custom panel: immune-cell genes + cytokine / receptor / signalling genes
immune_targets <- probes$gene[probes$panel == "custom" & probes$category %in% c("immune_cell", "immune_signalling")]
gene_panel <- function(g) ifelse(g %in% probes$gene, probes[g, "panel"], "other")
gene_cat <- function(g) ifelse(g %in% probes$gene, probes[g, "category"], "")

## family (from 04) and region: GI_glandular / Gastro_other (oesophagus, anus, liver, pancreas, ...) / non_GI
key <- function(x) tolower(trimws(gsub("[^[:print:]]", "", as.character(x))))
fam <- read_config_csv(ANNOT_FAMILIES_CSV, stringsAsFactors = FALSE, fileEncoding = "latin1")
if (!"family" %in% colnames(meta)) meta$family <- fam$family[match(key(meta$organ), key(fam$organ))]
meta$family[is.na(meta$family)] <- "none"
organ_group <- fam$group[match(key(meta$organ), key(fam$organ))]
meta$region <- dplyr::case_when(meta$family %in% GI_FAMILIES ~ "GI_glandular",
                                organ_group %in% "Gastro" ~ "Gastro_other",
                                TRUE ~ "non_GI")
REGIONS <- c("GI_glandular", "Gastro_other", "non_GI")

meta <- meta[!is.na(meta$core_id), ]
counts <- counts[, rownames(meta)]
core_info <- meta |>
  distinct(core_id, slide_id, patient_id, organ, tissue_type, diagnosis, location, family, region) |>
  distinct(core_id, .keep_all = TRUE) |>
  as.data.frame()
rownames(core_info) <- core_info$core_id
imm_types <- setdiff(unique(meta$cell_type[meta$lineage == "Immune"]), NA)
message("Cores per region: ", paste(names(table(core_info$region)), table(core_info$region), sep = "=", collapse = ", "))

## =============================================================================
## AIM 1
## =============================================================================
if ("1" %in% AIMS) {

## -----------------------------------------------------------------------
## 1.1 Signal-to-background per gene x core
##     lambda = counts per negative-control probe in this core = what a
##     non-binding probe collects here (autofluorescence, off-target
##     hybridisation, mis-decoding). Background differs a lot between organs
##     (e.g. autofluorescent lung/liver), so it must be per core, not global.
##     Whether a probe is above background is a property of the probe, not of
##     the segmentation: when the segmentation output has no control probes
##     (segger), the Xenium-segmentation object of the same cores is used.
##     CAVEAT: gene probe sets and control probes are not identical in probe
##     number/design, so treat SBR as a ranking, not an absolute FDR.
## -----------------------------------------------------------------------
ctrl_name <- function(m) {
  for (s in c("ControlProbe", "ControlCodeword")) {
    cc <- paste0("nCount_", s)
    if (cc %in% colnames(m) && sum(m[[cc]], na.rm = TRUE) > 0) return(s)
  }
  NULL
}
ctrl_col <- ctrl_name(meta)
bg_counts <- counts; bg_meta <- meta
if (is.null(ctrl_col) && SEGMENTATION != "xenium") {
  f <- obj_path("cores", "xenium")
  if (file.exists(f)) {
    message(SEGMENTATION, " output has no negative-control probe counts -> probe detection (1.1) and background ",
            "use the Xenium-segmentation object of the same cores: ", f)
    bg <- readRDS(f)
    bg_meta <- bg[[]]
    bg_meta <- bg_meta[!is.na(bg_meta$core_id) & bg_meta$core_id %in% core_info$core_id, ]
    bg_counts <- GetAssayData(bg, assay = "Xenium", layer = "counts")[, rownames(bg_meta)]
    rm(bg); invisible(gc())
    ctrl_col <- ctrl_name(bg_meta)
  }
}
if (is.null(ctrl_col)) stop("No negative-control probe counts found. For segger, first make the Xenium-segmentation ",
                            "object of the same cores: EXP1_SEGMENTATION=xenium Rscript 01_xenium_object_prep.R && ",
                            "EXP1_SEGMENTATION=xenium Rscript 02_tma_dearray.R")
message("Background model uses: ", ctrl_col)
bg_core <- bg_meta |>
  group_by(core_id) |>
  summarise(ctrl = sum(.data[[paste0("nCount_", ctrl_col)]]),
            n_ctrl = first(.data[[paste0("nFeatTotal_", ctrl_col)]]),
            n_cells = n(), .groups = "drop") |>
  mutate(lambda = ctrl / n_ctrl, bg_per_cell = lambda / n_cells)
save_table(bg_core, "05_aim1_background_per_core")
sbr_genes <- intersect(genes, rownames(bg_counts))
pb_bg <- pseudobulk(bg_counts[sbr_genes, , drop = FALSE], bg_meta$core_id)          # genes x cores
lambda <- stats::setNames(bg_core$lambda, bg_core$core_id)[colnames(pb_bg)]
rm(bg_counts); invisible(gc())

sbr <- as.data.frame(as.table(as.matrix(pb_bg)), stringsAsFactors = FALSE)
colnames(sbr) <- c("gene", "core_id", "counts")
sbr$lambda <- lambda[sbr$core_id]
sbr$log2_sbr <- log2((sbr$counts + 0.5) / (sbr$lambda + 0.5))
sbr$p <- stats::ppois(sbr$counts - 1, sbr$lambda, lower.tail = FALSE)
sbr$padj <- stats::p.adjust(sbr$p, method = "BH")
sbr$detected <- sbr$log2_sbr >= SBR_LOG2_MIN & sbr$padj < SBR_PADJ_MAX
sbr <- cbind(sbr, core_info[sbr$core_id, c("slide_id", "patient_id", "organ", "tissue_type", "family", "region")])
sbr$panel <- gene_panel(sbr$gene)
sbr$category <- gene_cat(sbr$gene)
save_table(sbr, "05_aim1_sbr_gene_by_core")

sbr_family <- sbr |>
  group_by(gene, panel, category, family) |>
  summarise(frac_cores_detected = mean(detected), median_log2_sbr = median(log2_sbr),
            n_cores = n(), .groups = "drop")
save_table(sbr_family, "05_aim1_sbr_gene_by_family")
save_table(sbr |> group_by(gene, panel, category, region, tissue_type, organ) |>
             summarise(frac_cores_detected = mean(detected), median_log2_sbr = median(log2_sbr),
                       n_cores = n(), .groups = "drop"), "05_aim1_sbr_gene_by_organ")

## Panel-wide view: how much of each panel is above background per family?
fam_order <- sbr_family |> group_by(family) |> summarise(m = mean(frac_cores_detected >= 0.5)) |>
  arrange(desc(m)) |> pull(family)
p_det <- sbr_family |>
  group_by(family, panel) |>
  summarise(pct_genes_detected = 100 * mean(frac_cores_detected >= 0.5), n_cores = first(n_cores), .groups = "drop") |>
  mutate(family = factor(family, levels = fam_order)) |>
  ggplot(aes(family, pct_genes_detected, fill = panel)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c(colon = "#3B6FB6", custom = "#D1495B")) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = "% of panel genes above background\nin >= half the cores",
       title = "Aim 1: how much of each panel 'works' per organ family")
save_fig(p_det, "05_aim1_pct_genes_detected_by_family", width = 11, height = 5.5)

## -----------------------------------------------------------------------
## 1.1b GI specificity of every probe.
##   Detection: fraction of GI_glandular vs non_GI cores in which the probe is
##   above background. Expression: CPM per core (each core weighted equally),
##   GI vs non-GI, on all cells AND within epithelium / within immune cells --
##   a gene that is "GI-specific" on all cells but not within immune cells is
##   driven by how many immune cells the GI cores contain, not by the probe.
## -----------------------------------------------------------------------
region_cpm <- function(cells, min_cells = MIN_CELLS_PSEUDOBULK) {
  n <- table(meta[cells, "core_id"])
  cells <- cells[meta[cells, "core_id"] %in% names(n)[n >= min_cells]]
  if (length(cells) == 0) return(NULL)
  pb <- pseudobulk(counts[, cells, drop = FALSE], meta[cells, "core_id"])
  cpm <- sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e6
  reg <- core_info[colnames(cpm), "region"]
  out <- sapply(REGIONS, function(r) if (any(reg == r)) rowMeans(cpm[, reg == r, drop = FALSE]) else rep(NA_real_, nrow(cpm)))
  rownames(out) <- rownames(cpm)
  out
}
lfc <- function(m) if (is.null(m)) NA_real_ else log2((m[, "GI_glandular"] + 1) / (m[, "non_GI"] + 1))
cpm_reg_all <- region_cpm(colnames(counts))
cpm_reg_epi <- region_cpm(colnames(counts)[meta$lineage %in% "Epithelial"])
cpm_reg_imm <- region_cpm(colnames(counts)[meta$lineage %in% "Immune"])

gi_spec <- sbr |>
  group_by(gene, panel, category) |>
  summarise(det_GI = mean(detected[region == "GI_glandular"]),
            det_gastro_other = mean(detected[region == "Gastro_other"]),
            det_nonGI = mean(detected[region == "non_GI"]),
            sbr_GI = median(log2_sbr[region == "GI_glandular"]),
            sbr_nonGI = median(log2_sbr[region == "non_GI"]), .groups = "drop") |>
  as.data.frame()
gi_spec$cpm_GI_all <- cpm_reg_all[gi_spec$gene, "GI_glandular"]
gi_spec$cpm_nonGI_all <- cpm_reg_all[gi_spec$gene, "non_GI"]
gi_spec$log2fc_GI_all <- lfc(cpm_reg_all)[gi_spec$gene]
gi_spec$log2fc_GI_epithelium <- if (is.null(cpm_reg_epi)) NA else lfc(cpm_reg_epi)[gi_spec$gene]
gi_spec$log2fc_GI_immune <- if (is.null(cpm_reg_imm)) NA else lfc(cpm_reg_imm)[gi_spec$gene]
gi_spec$class <- dplyr::case_when(
  gi_spec$det_GI < 0.5 & gi_spec$det_nonGI < 0.5 ~ "not above background",
  gi_spec$det_GI >= 0.5 & gi_spec$det_nonGI <= 0.2 ~ "GI-restricted",
  gi_spec$det_GI >= 0.5 & gi_spec$det_nonGI >= 0.5 ~ "broad (GI and non-GI)",
  gi_spec$det_GI >= 0.5 ~ "GI-enriched",
  TRUE ~ "mainly non-GI")
gi_spec <- gi_spec[order(gi_spec$panel, gi_spec$class, -gi_spec$det_GI), ]
save_table(gi_spec, "05_aim1_gi_specificity_per_gene")
class_tab <- gi_spec |> count(panel, category, class, name = "n_genes")
save_table(class_tab, "05_aim1_gi_specificity_summary")
print(as.data.frame(class_tab), row.names = FALSE)

class_cols <- c("GI-restricted" = "#1B7837", "GI-enriched" = "#7FBF7B", "broad (GI and non-GI)" = "#3B6FB6",
                "mainly non-GI" = "#E08214", "not above background" = "#BDBDBD")
p_cls <- ggplot(class_tab, aes(category, n_genes, fill = factor(class, levels = names(class_cols)))) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = class_cols, name = NULL) + scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "genes", title = "Aim 1: GI specificity of the probes",
       subtitle = "GI-restricted = above background in >= 50% of GI glandular cores and <= 20% of non-GI cores")
save_fig(p_cls, "05_aim1_gi_specificity_summary", width = 10, height = 6)

p_sc <- ggplot(gi_spec, aes(det_nonGI, det_GI)) +
  annotate("rect", xmin = -Inf, xmax = 0.2, ymin = 0.5, ymax = Inf, fill = "#1B7837", alpha = 0.08) +
  geom_abline(linetype = 2, color = "grey60") +
  geom_point(aes(color = class), size = 1.8, alpha = 0.8, position = position_jitter(0.008, 0.008, seed = 1)) +
  scale_color_manual(values = class_cols, name = NULL) +
  facet_wrap(~ panel) + coord_equal() + theme_bw(base_size = 12) +
  labs(x = "fraction of non-GI cores with the probe above background",
       y = "fraction of GI glandular cores above background",
       title = "Aim 1: probe detection, GI vs non-GI cores (shaded = GI-restricted)")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_sc <- p_sc + ggrepel::geom_text_repel(data = gi_spec[gi_spec$panel == "custom" & gi_spec$class != "not above background", ], aes(label = gene),
                                          size = 2.3, max.overlaps = 30, segment.size = 0.2)
}
save_fig(p_sc, "05_aim1_gi_detection_scatter", width = 14, height = 7.5)

## composition check: GI enrichment on all cells vs within the expected lineage
p_comp <- ggplot(gi_spec[gi_spec$category %in% c("immune_cell", "immune_signalling"), ],
                 aes(log2fc_GI_all, log2fc_GI_immune, color = panel)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(linetype = 2, color = "grey60") + geom_point(alpha = 0.7) +
  scale_color_manual(values = c(colon = "#3B6FB6", custom = "#D1495B")) + theme_bw(base_size = 12) +
  labs(x = "log2 GI / non-GI, all cells", y = "log2 GI / non-GI, within immune cells",
       title = "Aim 1: are immune probes GI-specific, or are there just more immune cells in GI cores?",
       subtitle = "near y = 0: same expression in immune cells of GI and non-GI organs (= the probe generalises)")
save_fig(p_comp, "05_aim1_gi_enrichment_all_vs_immune", width = 9, height = 7.5)

## custom panel: detection per family
if (length(custom_genes) > 0) {
  hm <- sbr_family |> filter(gene %in% custom_genes) |>
    mutate(family = factor(family, levels = fam_order),
           category = factor(category, levels = c("immune_cell", "immune_signalling", "non_immune")))
  p_hm <- ggplot(hm, aes(family, gene, fill = 100 * frac_cores_detected)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "#F2F2F2", high = "#B2182B", limits = c(0, 100), name = "% cores\nabove bg") +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.text.y = element_text(angle = 0),
          panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, title = "Aim 1: custom panel -- % of cores per family with the probe above background")
  save_fig(p_hm, "05_aim1_custom_detection_by_family", width = 10, height = max(6, 0.16 * length(custom_genes) + 2))
}

## -----------------------------------------------------------------------
## 1.2 Family specificity (tau): all cells vs within immune cells.
##     Each core weighted equally: CPM per core, then mean over the family's
##     cores (otherwise one huge core dominates its family).
## -----------------------------------------------------------------------
family_mean_cpm <- function(cells) {
  pb <- pseudobulk(counts[, cells, drop = FALSE], meta[cells, "core_id"])
  cpm <- sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e6
  f <- core_info[colnames(cpm), "family"]
  sapply(split(colnames(cpm), f), function(cc) rowMeans(cpm[, cc, drop = FALSE]))
}
cpm_all <- family_mean_cpm(colnames(counts))
imm_cells <- colnames(counts)[meta$lineage %in% "Immune"]
cpm_imm <- if (length(imm_cells) > 100) family_mean_cpm(imm_cells) else NULL

tau_tab <- data.frame(gene = genes, tau_all = apply(cpm_all, 1, tau_index))
tau_tab$top_family_all <- colnames(cpm_all)[apply(cpm_all, 1, which.max)]
if (!is.null(cpm_imm)) {
  tau_tab$tau_immune <- apply(cpm_imm, 1, tau_index)
  tau_tab$top_family_immune <- colnames(cpm_imm)[apply(cpm_imm, 1, which.max)]
}
tau_tab$panel <- gene_panel(tau_tab$gene)
tau_tab$category <- gene_cat(tau_tab$gene)
save_table(tau_tab, "05_aim1_family_specificity_tau")

if (!is.null(cpm_imm)) {
  lab <- tau_tab[tau_tab$gene %in% immune_targets, ]
  p_tau <- ggplot(tau_tab, aes(tau_all, tau_immune)) +
    geom_abline(linetype = 2, color = "grey50") +
    geom_point(aes(color = panel), alpha = 0.6) +
    scale_color_manual(values = c(colon = "#3B6FB6", custom = "#D1495B", other = "grey60")) +
    theme_bw(base_size = 12) +
    labs(x = "tau across families, all cells", y = "tau across families, immune cells only",
         title = "Aim 1: family specificity -- below the diagonal = driven by immune cell abundance, not the probe")
  if (requireNamespace("ggrepel", quietly = TRUE))
    p_tau <- p_tau + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.3, max.overlaps = 30)
  save_fig(p_tau, "05_aim1_tau_all_vs_immune", width = 10, height = 8)
}

## -----------------------------------------------------------------------
## 1.3 Lineage specificity per probe x family.
##     enrichment = (share of the gene's counts in its expected lineage) /
##                  (share of ALL transcripts in that lineage). 1 = no
##                  specificity; >>1 = on target.
##     offtarget_vs_bg = mean counts per NON-expected-lineage cell / what a
##                  negative-control probe collects per cell. ~1 = the
##                  off-target signal is just background; >>1 = real signal
##                  (or spill-over) outside the expected cells.
## -----------------------------------------------------------------------
lin_probes <- probes[nzchar(probes$expected_lineage), ]
if (nrow(lin_probes) > 0) {
  bg_cell <- stats::setNames(bg_core$bg_per_cell, bg_core$core_id)[meta$core_id]
  bg_cell[is.na(bg_cell)] <- 0
  grp <- paste(meta$family, meta$lineage, sep = "||")
  pb_fl <- as.matrix(pseudobulk(counts[lin_probes$gene, , drop = FALSE], grp))      # genes x family||lineage
  fl <- data.frame(group = colnames(pb_fl), family = sub("\\|\\|.*$", "", colnames(pb_fl)),
                   lineage = sub("^.*\\|\\|", "", colnames(pb_fl)))
  fl$n_cells <- as.numeric(table(grp)[fl$group])
  fl$tot <- as.numeric(tapply(meta$nCount_Xenium, grp, sum)[fl$group])
  fl$bg <- as.numeric(tapply(bg_cell, grp, sum)[fl$group])
  lin_spec <- dplyr::bind_rows(lapply(unique(fl$family), function(f) {
    cols <- which(fl$family == f)
    dplyr::bind_rows(lapply(seq_len(nrow(lin_probes)), function(j) {
      g <- lin_probes$gene[j]
      on <- cols[fl$lineage[cols] == lin_probes$expected_lineage[j]]
      off <- setdiff(cols, on)
      g_on <- sum(pb_fl[g, on]); g_all <- sum(pb_fl[g, cols])
      n_on <- sum(fl$n_cells[on]); n_off <- sum(fl$n_cells[off])
      data.frame(family = f, gene = g, panel = lin_probes$panel[j], category = lin_probes$category[j],
                 expected_lineage = lin_probes$expected_lineage[j], n_expected_cells = n_on,
                 enrichment = if (n_on > 0 && g_all > 0) (g_on / g_all) / (sum(fl$tot[on]) / sum(fl$tot[cols])) else NA,
                 mean_expected = if (n_on > 0) g_on / n_on else NA,
                 mean_offtarget = if (n_off > 0) (g_all - g_on) / n_off else NA,
                 offtarget_vs_bg = if (n_off > 0) ((g_all - g_on) / n_off) / max(sum(fl$bg[off]) / n_off, 1e-6) else NA)
    }))
  }))
  lin_spec$region <- ifelse(lin_spec$family %in% GI_FAMILIES, "GI_glandular", "other")
  save_table(lin_spec, "05_aim1_lineage_specificity")
  hm <- lin_spec[lin_spec$panel == "custom" & lin_spec$n_expected_cells >= 50, ]
  if (nrow(hm) > 0) {
    hm$family <- factor(hm$family, levels = fam_order)
    p_lin <- ggplot(hm, aes(family, gene, fill = pmax(pmin(log2(enrichment), 4), -4))) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_text(aes(label = sprintf("%.0f", offtarget_vs_bg)), size = 1.9) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2\nenrichment") +
      facet_grid(expected_lineage ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.text.y = element_text(angle = 0),
            panel.grid = element_blank()) +
      labs(x = NULL, y = NULL,
           title = "Aim 1: custom probes -- enrichment in their expected lineage, per family",
           subtitle = "red = on target in that family; label = off-target signal / negative-control background (per cell)")
    save_fig(p_lin, "05_aim1_lineage_specificity", width = 11, height = max(6, 0.16 * length(unique(hm$gene)) + 2))
  }
}

## -----------------------------------------------------------------------
## 1.4 Generalisability of the immune targets to immune cells of other organs.
##     For each immune cell type: pseudobulk per family over the custom
##     immune targets, correlated with the same cell type in GI glandular
##     cores. High correlation = the probes read out that cell type the same
##     way in GI and non-GI organs. Plus per gene: % of the cell type's cells
##     with >= 1 count, GI vs non-GI.
## -----------------------------------------------------------------------
gen_genes <- if (length(immune_targets) >= 3) immune_targets else genes
gen_tab <- list(); det_tab <- list()
for (ct in imm_types) {
  cells <- colnames(counts)[meta$cell_type %in% ct & meta$lineage %in% "Immune"]
  f <- meta[cells, "family"]
  ok <- names(which(table(f) >= MIN_CELLS_PSEUDOBULK))
  if (!any(ok %in% GI_FAMILIES) || length(ok) < 2) next
  cells <- cells[f %in% ok]
  reg <- ifelse(meta[cells, "family"] %in% GI_FAMILIES, "GI_glandular", "non_GI")
  pos <- counts[gen_genes, cells, drop = FALSE] > 0
  pct_gi <- 100 * Matrix::rowMeans(pos[, reg == "GI_glandular", drop = FALSE])
  det_tab[[ct]] <- data.frame(cell_type = ct, gene = gen_genes, pct_cells_GI = pct_gi,
                              pct_cells_nonGI = if (any(reg == "non_GI")) 100 * Matrix::rowMeans(pos[, reg == "non_GI", drop = FALSE]) else NA,
                              n_GI = sum(reg == "GI_glandular"), n_nonGI = sum(reg == "non_GI"))
  ## correlate only over targets this cell type actually expresses (>= 2% of its GI cells);
  ## genes at zero everywhere would only add ties and noise
  use <- gen_genes[pct_gi >= 2]
  if (length(use) < 5) next
  pb <- pseudobulk(counts[use, cells, drop = FALSE], meta[cells, "family"])
  lcpm <- log1p(sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e4)
  ref <- rowMeans(lcpm[, colnames(lcpm) %in% GI_FAMILIES, drop = FALSE])
  others <- setdiff(colnames(lcpm), GI_FAMILIES)
  if (length(others))
    gen_tab[[ct]] <- data.frame(cell_type = ct, family = others, n_genes = length(use),
                                n_cells = as.numeric(table(meta[cells, "family"])[others]),
                                spearman_vs_GI = sapply(others, function(o) stats::cor(lcpm[, o], ref, method = "spearman")))
}
gen_tab <- dplyr::bind_rows(gen_tab); det_tab <- dplyr::bind_rows(det_tab)
if (nrow(gen_tab) > 0) {
  save_table(gen_tab, "05_aim1_generalisability_immune_profiles")
  ct_order <- gen_tab |> group_by(cell_type) |> summarise(m = median(spearman_vs_GI, na.rm = TRUE)) |>
    arrange(desc(m)) |> pull(cell_type)
  p_gen <- ggplot(gen_tab, aes(factor(cell_type, levels = ct_order), spearman_vs_GI)) +
    geom_boxplot(outlier.shape = NA, color = "grey40") +
    geom_point(aes(size = n_cells), position = position_jitter(width = 0.15, seed = 1), alpha = 0.6, color = "#3B6FB6") +
    scale_size_area(max_size = 3.5, name = "cells") + coord_cartesian(ylim = c(min(0, gen_tab$spearman_vs_GI, na.rm = TRUE), 1)) +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = "Spearman r vs the same cell type in GI cores",
         title = "Aim 1: do the custom immune targets read out immune cells the same way outside the GI tract?",
         subtitle = "one dot = one non-GI organ family; genes = custom immune targets expressed by >= 2% of that cell type in GI cores")
  save_fig(p_gen, "05_aim1_generalisability", width = 12, height = 6.5)
}
if (nrow(det_tab) > 0) {
  save_table(det_tab, "05_aim1_immune_target_detection_GI_vs_nonGI")
  p_det2 <- ggplot(det_tab, aes(pct_cells_GI, pct_cells_nonGI)) +
    geom_abline(linetype = 2, color = "grey60") + geom_point(alpha = 0.5, size = 1, color = "#D1495B") +
    facet_wrap(~ cell_type) + coord_equal() + theme_bw(base_size = 10) +
    labs(x = "% cells with >= 1 count, GI glandular cores", y = "% cells, non-GI cores",
         title = "Aim 1: custom immune targets per immune cell type, GI vs non-GI (one dot = one gene)")
  save_fig(p_det2, "05_aim1_immune_target_detection_GI_vs_nonGI", width = 14, height = 12)
}

if (length(immune_targets) > 0 && length(imm_cells) > 0) {
  sub_imm <- CreateSeuratObject(counts = counts[, imm_cells], assay = "Xenium", meta.data = meta[imm_cells, ])
  sub_imm <- NormalizeData(sub_imm, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)
  sub_imm$type_region <- paste(sub_imm$cell_type, ifelse(sub_imm$region == "GI_glandular", "GI", "non-GI"), sep = " | ")
  p_dot_imm <- marker_dotplot(sub_imm, features = immune_targets, group.by = "type_region") +
    RotatedAxis() + ggtitle("Aim 1: custom immune targets within immune cell types, GI vs non-GI cores")
  save_fig(p_dot_imm, "05_aim1_immune_targets_dotplot_GI_vs_nonGI", width = 22,
           height = max(6, 0.22 * length(unique(sub_imm$type_region)) + 3))
  rm(sub_imm); invisible(gc())
}
} ## end AIM 1

if ("2" %in% AIMS) {
  ## =============================================================================
  ## AIM 2
  ## =============================================================================

  ## -----------------------------------------------------------------------
  ## 2.1 Composition per core, heterogeneity, differential composition
  ## -----------------------------------------------------------------------
  comp <- meta |>
    count(core_id, cell_type) |>
    group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup() |>
    left_join(core_info, by = "core_id")
  save_table(comp, "05_aim2_composition_per_core")

  entropy <- comp |>
    group_by(core_id, organ, patient_id, tissue_type) |>
    summarise(shannon = -sum(frac * log(frac)), immune_frac = sum(frac[cell_type %in% imm_types]),
              .groups = "drop")
  save_table(entropy, "05_aim2_core_entropy")

  p_imm <- comp |> filter(cell_type %in% imm_types) |>
    ggplot(aes(organ, frac)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(color = patient_id), position = position_jitter(width = 0.2), size = 1.2) +
    facet_wrap(~ cell_type, scales = "free_y") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(y = "fraction of all cells in core", title = "Aim 2: immune subset abundance per core, by organ (dot colour = patient)")
  save_fig(p_imm, "05_aim2_immune_fractions_by_organ", width = 16, height = 10)

  ## Heterogeneity within organ: CV of each cell type's fraction across cores
  cv_tab <- comp |>
    group_by(organ, cell_type) |>
    summarise(n_cores = n_distinct(core_id), mean_frac = mean(frac),
              cv = stats::sd(frac) / mean(frac), .groups = "drop")
  save_table(cv_tab, "05_aim2_composition_cv_within_organ")

  if (requireNamespace("speckle", quietly = TRUE)) {
    cores_per_organ <- table(core_info$organ)
    ok <- names(cores_per_organ)[cores_per_organ >= 2]
    if (length(ok) >= 2) {
      sel <- meta$organ %in% ok
      prop <- speckle::propeller(clusters = meta$cell_type[sel], sample = meta$core_id[sel],
                                 group = meta$organ[sel], transform = "logit")
      save_table(cbind(cell_type = rownames(prop), prop), "05_aim2_propeller_organ")
      message("propeller treats cores as independent -- cores from one patient are not; ",
              "read p-values as screening, confirm at patient level.")
    }
  } else message("speckle not installed -- skipping propeller.")

  ## -----------------------------------------------------------------------
  ## 2.2 Pseudobulk per core x cell type, variance partitioning
  ## -----------------------------------------------------------------------
  grp <- paste(meta$core_id, meta$cell_type, sep = "||")
  n_grp <- table(grp)
  keep_grp <- names(n_grp)[n_grp >= MIN_CELLS_PSEUDOBULK]
  pb <- pseudobulk(counts[, grp %in% keep_grp, drop = FALSE], grp[grp %in% keep_grp])
  pb_meta <- data.frame(group = colnames(pb),
                        core_id = sub("\\|\\|.*$", "", colnames(pb)),
                        cell_type = sub("^.*\\|\\|", "", colnames(pb)),
                        n_cells = as.integer(n_grp[colnames(pb)]))
  pb_meta <- cbind(pb_meta, core_info[pb_meta$core_id, c("slide_id", "patient_id", "organ", "tissue_type", "location")])
  rownames(pb_meta) <- pb_meta$group
  lcpm <- edgeR::cpm(as.matrix(pb), log = TRUE, prior.count = 1)
  saveRDS(list(counts = pb, meta = pb_meta), file.path(PROCESSED_DIR, "05_pseudobulk_core_celltype.rds"))

  vp_terms <- c("cell_type", "organ", "patient_id", "slide_id", "location", "tissue_type")
  vp_terms <- vp_terms[vapply(vp_terms, function(v) length(unique(stats::na.omit(pb_meta[[v]]))) > 1, logical(1))]
  vp_genes <- if (length(immune_targets) >= 3) unique(c(immune_targets, custom_genes)) else rownames(lcpm)
  if (requireNamespace("variancePartition", quietly = TRUE) && length(vp_terms) >= 2) {
    pbm <- pb_meta
    for (v in vp_terms) pbm[[v]] <- ifelse(is.na(pbm[[v]]) | pbm[[v]] == "", "NA", as.character(pbm[[v]]))
    form <- stats::as.formula(paste("~", paste0("(1|", vp_terms, ")", collapse = " + ")))
    message("variancePartition formula: ", deparse(form))
    vp <- tryCatch(
      variancePartition::fitExtractVarPartModel(lcpm[vp_genes, , drop = FALSE], form, pbm),
      error = function(e) { message("variancePartition failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(vp)) {
      vp_df <- cbind(gene = rownames(vp), as.data.frame(vp))
      save_table(vp_df, "05_aim2_variance_partition")
      save_fig(variancePartition::plotVarPart(vp) +
                 ggtitle("Aim 2: variance explained per factor (pseudobulk core x cell type)"),
               "05_aim2_variance_partition", width = 10, height = 6)
      message("If slide_id explains more than patient_id, correct batch before cross-patient comparisons; ",
              "if organ and slide are confounded (all cores of an organ on one slide) they cannot be separated.")
    }
  } else message("variancePartition not installed or too few factors -- skipping 2.2.")

  ## -----------------------------------------------------------------------
  ## 2.3 Core-to-core similarity (whole-core pseudobulk, and immune-only)
  ## -----------------------------------------------------------------------
  pb_core <- pseudobulk(counts, meta$core_id)
  core_lcpm <- edgeR::cpm(as.matrix(pb_core), log = TRUE, prior.count = 1)
  cc <- stats::cor(core_lcpm, method = "spearman")
  ann <- core_info[colnames(cc), c("organ", "patient_id", "slide_id")]
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    grDevices::png(file.path(FIGURES_DIR, "05_aim2_core_similarity_heatmap.png"), width = 12, height = 11,
                   units = "in", res = 150)
    pheatmap::pheatmap(cc, annotation_col = ann, annotation_row = ann, show_rownames = FALSE,
                       show_colnames = FALSE, main = "Aim 2: core similarity (Spearman, pseudobulk) -- clusters by organ, patient or slide?")
    grDevices::dev.off()
  }
  pairs_cc <- which(upper.tri(cc), arr.ind = TRUE)
  pair_df <- data.frame(a = rownames(cc)[pairs_cc[, 1]], b = colnames(cc)[pairs_cc[, 2]], r = cc[pairs_cc])
  pair_df <- pair_df |>
    mutate(same_organ = core_info[a, "organ"] == core_info[b, "organ"],
           same_patient = core_info[a, "patient_id"] == core_info[b, "patient_id"],
           same_slide = core_info[a, "slide_id"] == core_info[b, "slide_id"],
           relation = case_when(
             same_organ & same_patient ~ "same organ, same patient",
             same_organ & !same_patient ~ "same organ, different patient",
             !same_organ & same_patient ~ "different organ, same patient",
             TRUE ~ "different organ, different patient"))
  save_table(pair_df, "05_aim2_core_pair_similarity")
  p_pairs <- ggplot(pair_df, aes(relation, r, fill = same_slide)) +
    geom_boxplot(outlier.size = 0.5) + theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(y = "Spearman r between cores", fill = "same slide",
         title = "Aim 2: can we compare between patients? (same-organ, different-patient vs same-patient)")
  save_fig(p_pairs, "05_aim2_core_pair_similarity", width = 11, height = 6)

  ## -----------------------------------------------------------------------
  ## 2.4 Location within organ: pseudobulk limma-voom per organ x cell type,
  ##     patient as blocking factor (duplicateCorrelation). Runs only where
  ##     the design has >= 2 locations with >= 2 cores each.
  ## -----------------------------------------------------------------------
  loc_results <- list()
  for (org in unique(pb_meta$organ)) {
    for (ct in unique(pb_meta$cell_type)) {
      sel <- pb_meta$organ == org & pb_meta$cell_type == ct & !is.na(pb_meta$location) & pb_meta$location != ""
      d <- pb_meta[sel, ]
      loc_n <- table(d$location)
      if (sum(loc_n >= 2) < 2) next
      d <- d[d$location %in% names(loc_n)[loc_n >= 2], ]
      y <- edgeR::DGEList(as.matrix(pb[, d$group, drop = FALSE]))
      y <- y[edgeR::filterByExpr(y, group = d$location), , keep.lib.sizes = FALSE]
      if (nrow(y) < 5) next
      y <- edgeR::calcNormFactors(y)
      design <- stats::model.matrix(~ 0 + location, data = d)
      colnames(design) <- make.names(sub("^location", "", colnames(design)))
      v <- limma::voom(y, design)
      ## Blocking on patient only helps if a patient contributes cores to more
      ## than one location (paired design); otherwise patient is nested in location.
      block <- d$patient_id
      use_block <- any(tapply(d$location, block, function(v) length(unique(v))) > 1)
      fit <- if (use_block) {
        corfit <- limma::duplicateCorrelation(v, design, block = block)
        limma::lmFit(v, design, block = block, correlation = corfit$consensus)
      } else limma::lmFit(v, design)
      locs <- colnames(design)
      for (k in 2:length(locs)) {
        con <- limma::makeContrasts(contrasts = paste0(locs[k], "-", locs[1]), levels = design)
        tt <- limma::topTable(limma::eBayes(limma::contrasts.fit(fit, con)), number = Inf)
        tt$gene <- rownames(tt); tt$organ <- org; tt$cell_type <- ct
        tt$contrast <- paste0(locs[k], " vs ", locs[1]); tt$patient_blocked <- use_block
        loc_results[[paste(org, ct, k)]] <- tt
      }
    }
  }
  if (length(loc_results) > 0) {
    save_table(dplyr::bind_rows(loc_results), "05_aim2_location_de")
  } else message("No organ had >= 2 locations with >= 2 cores each -- 2.4 skipped (check `location` in tma_map.csv).")

  ## -----------------------------------------------------------------------
  ## 2.5 Cellular niches: composition of each cell's K nearest neighbours
  ##     (within its core), k-means into N_NICHES niches.
  ## -----------------------------------------------------------------------
  types <- sort(unique(meta$cell_type))
  nbhd <- matrix(0, nrow = nrow(meta), ncol = length(types), dimnames = list(rownames(meta), types))
  for (cid in unique(meta$core_id)) {
    idx <- which(meta$core_id == cid)
    if (length(idx) <= NICHE_K_NEIGHBOURS) next
    nn <- RANN::nn2(meta[idx, c("x_um", "y_um")], k = NICHE_K_NEIGHBOURS + 1)$nn.idx[, -1]
    lab <- matrix(match(meta$cell_type[idx][nn], types), nrow = length(idx))
    for (t in seq_along(types)) nbhd[idx, t] <- rowMeans(lab == t)
  }
  km <- stats::kmeans(nbhd, centers = N_NICHES, nstart = 5, iter.max = 50)
  meta$niche <- paste0("N", km$cluster)
  xen <- AddMetaData(xen, stats::setNames(meta$niche, rownames(meta)), col.name = "niche")
  niche_comp <- km$centers
  rownames(niche_comp) <- paste0("N", seq_len(N_NICHES))
  save_table(cbind(niche = rownames(niche_comp), as.data.frame(niche_comp)), "05_aim2_niche_composition")
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    grDevices::png(file.path(FIGURES_DIR, "05_aim2_niche_composition.png"), width = 10, height = 6, units = "in", res = 150)
    pheatmap::pheatmap(niche_comp, main = "Niche = mean neighbourhood composition")
    grDevices::dev.off()
  }
  niche_by_organ <- meta |> count(organ, core_id, niche) |> group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
  p_niche <- ggplot(niche_by_organ, aes(core_id, frac, fill = niche)) +
    geom_col() + facet_grid(~ organ, scales = "free_x", space = "free_x") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
    labs(title = "Aim 2: niche composition per core")
  save_fig(p_niche, "05_aim2_niches_per_core", width = 18, height = 6)
  saveRDS(xen, obj_path("insights"))
  message("Saved -> ", obj_path("insights"))
} ## end AIM 2
