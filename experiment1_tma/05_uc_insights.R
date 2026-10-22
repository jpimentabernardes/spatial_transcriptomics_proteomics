## =============================================================================
## 05_uc_insights.R
##
## Ulcerative colitis biopsy TMAs (TMA_UC_A, TMA_UC_CD), after 04 with the UC
## marker panel (config/marker_panel_rna_uc.csv). One script per question block:
##
## AIM 1 -- does the panel work in UC biopsies, and which cell types are missing?
##   1.1 Probe signal over background per gene x core (negative-control probes;
##       segger has none -> the Xenium-segmentation object of the same cores),
##       colon vs custom panel, low vs high Nancy cores.
##   1.2 Custom probes: share of each probe's transcripts in its expected lineage
##       (epithelium / immune / stroma), low vs high Nancy -- on target in
##       inflamed tissue too?
##   1.3 Missing cell types: for every level-2 signature of the marker panel, cells
##       co-expressing >= 2 of its genes, and the labels they received. Many
##       co-expressing cells but no matching label = a cell type the annotation
##       missed (merged into another cluster); very few = not detectable with this
##       panel / segmentation. Plus cells left unlabelled / unclear per core.
##
## AIM 2 -- heterogeneity within the same Nancy score
##   Per-core features (R/uc_helpers.R): compartment and cell-class fractions,
##   densities per mm^2, intraepithelial neutrophils / eosinophils / T cells,
##   neutrophil aggregates (crypt-abscess proxy), plasma cells near the crypt
##   base (basal plasmacytosis proxy), goblet ratio, niche fractions.
##   2.1 composition per core, by Nancy grade
##   2.2 cellular niches (neighbourhood composition) per core
##   2.3 feature heatmap, cores ordered by Nancy grade
##   2.4 spread within each grade: PCA, distance to the grade centre, which features
##       make a core deviate; same patient vs different patients
##
## AIM 3 -- can spatial features sharpen the Nancy call, and do they relate to
##          Comprehensive Disease Control (Disease.Control)?
##   3.1 agreement Nancy_TMA vs Nancy (tma_map) vs the clinical Nancy columns
##   3.2 each feature vs Nancy (core level; mixed model with patient)
##   3.3 ordinal model Nancy ~ a-priori spatial features, leave-one-patient-out:
##       a continuous "spatial Nancy" per core; cores far from their grade
##   3.4 patient level: features vs Disease.Control and the remission outcomes,
##       and whether a feature adds to Nancy (cross-validated AUC)
##
## Statistical unit = core (Nancy_TMA is scored per core) or patient (outcomes);
## never cell. With few patients per group read p-values as screening.
##
## Run:   Rscript 05_uc_insights.R                 (all aims)
##        EXP1_UC_AIMS=1 Rscript 05_uc_insights.R  (Aim 1 only; 2,3 = Aims 2 and 3)
##
## INPUT : xenium_<seg>_annotated.rds (04, UC slides), tma_map columns Nancy / Nancy_TMA
##         (copied into the object by 02), config/marker_panel_rna_uc.csv,
##         config/probe_expectations.csv, de-identified clinical table
##         (UC_CLINICAL_CSV, default config/clinical_uc.csv; optional)
## OUTPUT: tables/05uc_*.csv, figures/05uc_*.png,
##         data/processed/05uc_core_features.rds (input of 05_uc_genes_eosinophils.R)
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

AIMS <- trimws(strsplit(Sys.getenv("EXP1_UC_AIMS", "1,2,3"), ",")[[1]])
## Aim 3.3 model features (chosen before looking at the data: the Nancy items --
## acute infiltrate, chronic infiltrate, epithelial damage); missing ones are skipped
UC_MODEL_FEATURES <- get0("UC_MODEL_FEATURES", ifnotfound = c(
  "Neutrophil_per_mm2", "Neutrophil_intraepithelial_pct", "Neutrophil_aggregated_pct", "Plasma_per_mm2",
  "Immune_pct", "Eosinophil_per_mm2", "Inflamed_epi_pct_epithelium", "Goblet_ratio", "Epithelium_pct"))

xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
genes <- rownames(counts)
core_tab <- uc_core_table(m)
clin <- join_uc_clinical(core_tab, read_uc_clinical())
m$nancy <- core_tab[m$core_id, "nancy"]
m$nancy_group <- core_tab[m$core_id, "nancy_group"]
nancy_col <- attr(core_tab, "nancy_col")
if (all(is.na(core_tab$nancy))) message("!! No Nancy score in the object (columns Nancy_TMA / Nancy of tma_map.csv) -- ",
                                        "re-run 02 after adding them; Nancy comparisons are skipped.")
message("Nancy per core from: ", nancy_col, "  |  cores: ", nrow(core_tab), ", patients: ",
        length(unique(core_tab$patient_id)), ", cells: ", format(nrow(m), big.mark = ","))
cohort <- core_tab |> count(nancy = ifelse(is.na(nancy), "NA", nancy), name = "n_cores") |>
  left_join(core_tab |> group_by(nancy = ifelse(is.na(nancy), "NA", nancy)) |>
              summarise(n_patients = n_distinct(patient_id), median_cells = stats::median(n_cells), .groups = "drop"),
            by = "nancy")
save_table(cohort, "05uc_cores_per_nancy")
print(as.data.frame(cohort), row.names = FALSE)
lab_pal <- celltype_palette(m$label, m$lineage)

probes <- if (file.exists(PROBE_EXPECT_CSV)) read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE) else
  data.frame(gene = character(0), panel = character(0), category = character(0), expected_lineage = character(0))
probes$gene <- match_case(as_seurat_features(probes$gene), genes)
probes <- probes[probes$gene %in% genes & !duplicated(probes$gene), ]
for (v in c("panel", "category", "expected_lineage")) probes[[v]][is.na(probes[[v]])] <- ""
rownames(probes) <- probes$gene
gene_panel <- function(g) ifelse(g %in% probes$gene, probes[g, "panel"], "other")

## =============================================================================
## AIM 1
## =============================================================================
if ("1" %in% AIMS) {
## -----------------------------------------------------------------------
## 1.1 signal over background per gene x core
## -----------------------------------------------------------------------
ctrl_name <- function(mm) {
  for (s in c("ControlProbe", "ControlCodeword")) {
    cc <- paste0("nCount_", s)
    if (cc %in% colnames(mm) && sum(mm[[cc]], na.rm = TRUE) > 0) return(s)
  }
  NULL
}
bg_meta <- m; bg_counts <- counts
ctrl_col <- ctrl_name(m)
if (is.null(ctrl_col) && SEGMENTATION != "xenium" && file.exists(obj_path("cores", "xenium"))) {
  message(SEGMENTATION, " has no negative-control probes -> background from the Xenium-segmentation object ",
          "of the same cores: ", obj_path("cores", "xenium"))
  bg <- readRDS(obj_path("cores", "xenium"))
  bg_meta <- bg[[]]; bg_meta <- bg_meta[!is.na(bg_meta$core_id) & bg_meta$core_id %in% core_tab$core_id, ]
  bg_counts <- GetAssayData(bg, assay = "Xenium", layer = "counts")[, rownames(bg_meta)]
  rm(bg); invisible(gc())
  ctrl_col <- ctrl_name(bg_meta)
}
bg_core <- NULL
if (is.null(ctrl_col)) {
  message("!! 1.1 skipped: no negative-control probe counts. For segger first run 01 + 02 with ",
          "EXP1_SEGMENTATION=xenium (makes xenium_xenium_cores.rds).")
} else {
  bg_core <- bg_meta |> group_by(core_id) |>
    summarise(ctrl = sum(.data[[paste0("nCount_", ctrl_col)]]), n_ctrl = first(.data[[paste0("nFeatTotal_", ctrl_col)]]),
              n_cells_bg = n(), .groups = "drop") |>
    mutate(lambda = ctrl / n_ctrl, bg_per_cell = lambda / n_cells_bg)
  pb_bg <- as.matrix(pseudobulk(bg_counts[intersect(genes, rownames(bg_counts)), , drop = FALSE], bg_meta$core_id))
  lam <- stats::setNames(bg_core$lambda, bg_core$core_id)[colnames(pb_bg)]
  sbr <- data.frame(gene = rep(rownames(pb_bg), ncol(pb_bg)), core_id = rep(colnames(pb_bg), each = nrow(pb_bg)),
                    counts = as.vector(pb_bg))
  sbr$lambda <- lam[sbr$core_id]
  sbr$log2_sbr <- log2((sbr$counts + 0.5) / (sbr$lambda + 0.5))
  sbr$padj <- stats::p.adjust(stats::ppois(sbr$counts - 1, sbr$lambda, lower.tail = FALSE), "BH")
  sbr$detected <- sbr$log2_sbr >= SBR_LOG2_MIN & sbr$padj < SBR_PADJ_MAX
  sbr$nancy_group <- core_tab[sbr$core_id, "nancy_group"]
  sbr$panel <- gene_panel(sbr$gene)
  per_gene <- sbr |> group_by(gene, panel) |>
    summarise(pct_cores_detected = 100 * mean(detected),
              pct_cores_detected_low = 100 * mean(detected[nancy_group %in% levels(nancy_group)[1]]),
              pct_cores_detected_high = 100 * mean(detected[nancy_group %in% levels(nancy_group)[2]]),
              median_log2_sbr = stats::median(log2_sbr), median_counts_per_core = stats::median(counts), .groups = "drop") |>
    mutate(category = ifelse(gene %in% probes$gene, probes[gene, "category"], ""),
           usable = ifelse(pct_cores_detected >= 50, "above background in >= 50% of cores", "mostly at background")) |>
    arrange(panel, desc(pct_cores_detected))
  save_table(per_gene, "05uc_aim1_probe_detection_per_gene")
  message("Genes above background in >= 50% of UC cores: ",
          paste(sprintf("%s %d/%d", names(table(per_gene$panel)), tapply(per_gene$pct_cores_detected >= 50, per_gene$panel, sum),
                        as.integer(table(per_gene$panel))), collapse = ", "))
  summ <- sbr |> filter(!is.na(nancy_group)) |> group_by(gene, panel, nancy_group) |>
    summarise(det = mean(detected) >= 0.5, .groups = "drop") |>
    group_by(panel, nancy_group) |> summarise(pct_genes = 100 * mean(det), n_genes = n(), .groups = "drop")
  if (nrow(summ)) {
    p <- ggplot(summ, aes(panel, pct_genes, fill = nancy_group)) +
      geom_col(position = position_dodge(0.8), width = 0.75) +
      geom_text(aes(label = sprintf("%.0f%%", pct_genes)), position = position_dodge(0.8), vjust = -0.4, size = 3.5) +
      scale_fill_manual(values = GROUP_COLOURS(), name = NULL) + scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
      theme_bw(base_size = 13) + theme(legend.position = "top") +
      labs(x = NULL, y = "% of the panel's genes above background\nin >= half of the cores",
           title = "Aim 1.1: how much of each panel works in UC biopsies")
    save_fig(p, "05uc_aim1_panel_detection", width = 7, height = 5.5)
  }
  cg <- per_gene[per_gene$panel == "custom", ]
  if (nrow(cg)) {
    cg$gene <- factor(cg$gene, levels = cg$gene[order(cg$pct_cores_detected)])
    d <- tidyr::pivot_longer(cg, c(pct_cores_detected_low, pct_cores_detected_high), names_to = "grp", values_to = "pct")
    d$grp <- factor(ifelse(d$grp == "pct_cores_detected_low", levels(nancy_group(0))[1], levels(nancy_group(0))[2]),
                    levels = levels(nancy_group(0)))
    p <- ggplot(d, aes(pct, gene)) +
      geom_line(aes(group = gene), colour = "grey70") + geom_point(aes(colour = grp), size = 2.2) +
      geom_vline(xintercept = 50, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = GROUP_COLOURS(), name = NULL) +
      facet_grid(category ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 10) +
      theme(legend.position = "top", strip.text.y = element_text(angle = 0)) +
      labs(x = "% of cores with the probe above background", y = NULL,
           title = "Aim 1.1: custom panel in UC biopsies, low vs high Nancy cores")
    save_fig(p, "05uc_aim1_custom_detection", width = 9, height = max(6, 0.14 * nrow(cg) + 2))
  }
}

## -----------------------------------------------------------------------
## 1.2 lineage specificity of the custom probes, low vs high Nancy
##     on_target_pct  = % of the probe's transcripts in cells of its expected lineage
##     lineage_share  = % of ALL transcripts in that lineage (what a non-specific
##                      probe would give); enrichment = on_target_pct / lineage_share
## -----------------------------------------------------------------------
lin_comp <- function(l) dplyr::case_when(l %in% EPI_LINEAGES ~ "Epithelium", l %in% IMMUNE_LINEAGES ~ "Immune",
                                        l %in% STROMA_LINEAGES ~ "Stroma", TRUE ~ NA_character_)
lp <- probes[nzchar(probes$expected_lineage) & probes$panel == "custom", ]
lp$expected <- lin_comp(lp$expected_lineage)
lp <- lp[!is.na(lp$expected), ]
if (nrow(lp) && any(!is.na(m$nancy_group))) {
  grp <- paste(m$nancy_group, m$compartment, sep = "||")
  pb <- as.matrix(pseudobulk(counts[lp$gene, , drop = FALSE], grp))
  tot <- tapply(m$nCount_Xenium, grp, sum)
  lin_spec <- dplyr::bind_rows(lapply(levels(m$nancy_group), function(g) {
    cols <- colnames(pb)[startsWith(colnames(pb), paste0(g, "||"))]
    comp <- sub("^.*\\|\\|", "", cols)
    dplyr::bind_rows(lapply(seq_len(nrow(lp)), function(j) {
      on <- cols[comp == lp$expected[j]]
      g_all <- sum(pb[lp$gene[j], cols])
      data.frame(nancy_group = g, gene = lp$gene[j], category = lp$category[j], expected = lp$expected[j],
                 transcripts = g_all,
                 on_target_pct = if (g_all > 0) 100 * sum(pb[lp$gene[j], on]) / g_all else NA,
                 lineage_share = 100 * sum(tot[on]) / sum(tot[cols]))
    }))
  }))
  lin_spec$enrichment <- lin_spec$on_target_pct / lin_spec$lineage_share
  save_table(lin_spec, "05uc_aim1_custom_lineage_specificity")
  d <- lin_spec[lin_spec$transcripts >= 20, ]
  if (nrow(d)) {
    ord <- d |> group_by(gene) |> summarise(v = mean(on_target_pct, na.rm = TRUE)) |> arrange(v) |> pull(gene)
    d$gene <- factor(d$gene, levels = ord)
    d$nancy_group <- factor(d$nancy_group, levels = levels(nancy_group(0)))
    p <- ggplot(d, aes(on_target_pct, gene)) +
      geom_point(aes(x = lineage_share), shape = 124, size = 4, colour = "grey45") +
      geom_line(aes(group = gene), colour = "grey75") +
      geom_point(aes(colour = nancy_group), size = 2.3) +
      scale_colour_manual(values = GROUP_COLOURS(), name = NULL) + scale_x_continuous(limits = c(0, 100)) +
      facet_grid(expected ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 10) +
      theme(legend.position = "top", strip.text.y = element_text(angle = 0)) +
      labs(x = "% of the probe's transcripts in its expected compartment", y = NULL,
           title = "Aim 1.2: are the custom probes on target in UC biopsies?",
           subtitle = "dot = probe; grey tick = that compartment's share of all transcripts (what a non-specific probe would give).\nDot far right of the tick = on target; dot at the tick = no specificity")
    save_fig(p, "05uc_aim1_custom_lineage_specificity", width = 9, height = max(6, 0.15 * length(ord) + 2.5))
  }
}

## -----------------------------------------------------------------------
## 1.3 missing cell types: signature co-expression vs assigned labels
## -----------------------------------------------------------------------
mk <- read_marker_panel(MARKERS_RNA_CSV)
sigs <- unlist(unname(lapply(names(mk$level2), function(par)
  stats::setNames(mk$level2[[par]], paste0(par, ": ", names(mk$level2[[par]]))))), recursive = FALSE)
sigs <- lapply(sigs, function(g) intersect(match_case(g, genes), genes))
sigs <- sigs[lengths(sigs) >= 2 & !grepl(": Proliferating$", names(sigs))]
## names as "_word_word_" (plural s dropped), so whole words are compared: T_cell does not match Mast_cell
norm_name <- function(x) vapply(strsplit(tolower(x), "[^[:alnum:]]+"), function(w)
  paste0("_", paste(sub("(.{3,})s$", "\\1", w[nzchar(w)]), collapse = "_"), "_"), character(1))
labels_n <- table(m$label)
miss_tab <- list(); miss_hm <- list()
for (s in names(sigs)) {
  ct <- sub("^.*: ", "", s)
  B <- counts[sigs[[s]], , drop = FALSE] > 0
  pos <- Matrix::colSums(B) >= 2
  n_pos <- sum(pos)
  ## a label "matches" when one name contains the other (Plasma_cell ~ Plasma cells), or when both fall in
  ## the same narrow class (Neutrophil, Goblet, ...); broad classes (T cells, macrophages, ...) don't count,
  ## so Treg is not "found" just because CD4_T exists
  sig_class <- uc_classes(ct)
  ln <- norm_name(names(labels_n)); sn <- norm_name(ct)
  by_name <- vapply(ln, function(l) grepl(sn, l, fixed = TRUE) ||
                        (grepl(l, sn, fixed = TRUE) && nchar(l) >= 0.7 * nchar(sn)), logical(1))
  broad <- c("T_cell", "Macrophage", "DC", "NK_ILC", "Colonocyte", "Stem_TA", "Endothelial", "Fibroblast")
  generic <- !is.na(sig_class) && norm_name(ct) == norm_name(sig_class)   # the T_cell signature itself: CD4_T counts
  by_class <- !is.na(sig_class) & (!sig_class %in% broad | generic) & uc_classes(names(labels_n)) %in% sig_class
  match_lab <- names(labels_n)[by_name | by_class]
  match_lab <- match_lab[labels_n[match_lab] >= 20]
  lab_pos <- table(m$label[pos])
  pct_match <- if (n_pos) 100 * sum(lab_pos[intersect(names(lab_pos), match_lab)]) / n_pos else NA
  top_other <- sort(lab_pos[setdiff(names(lab_pos), match_lab)], decreasing = TRUE)
  status <- if (length(match_lab)) "found (label exists)" else
    if (n_pos >= max(50, 5e-4 * nrow(m))) "MISSING: co-expressing cells, but no label" else "not detectable (too few co-expressing cells)"
  miss_tab[[s]] <- data.frame(signature = s, cell_type = ct, genes_on_panel = paste(sigs[[s]], collapse = ";"),
                              n_coexpressing_cells = n_pos, pct_of_all_cells = round(100 * n_pos / nrow(m), 3),
                              matching_labels = paste(match_lab, collapse = ";"),
                              pct_coexpressing_in_matching_label = round(pct_match, 1),
                              top_other_labels = paste(sprintf("%s (%.0f%%)", names(utils::head(top_other, 3)),
                                                               100 * utils::head(top_other, 3) / max(n_pos, 1)), collapse = "; "),
                              status = status)
  if (n_pos) miss_hm[[s]] <- data.frame(signature = s, label = names(lab_pos), pct = 100 * as.numeric(lab_pos) / n_pos)
}
miss_tab <- dplyr::bind_rows(miss_tab)
save_table(miss_tab, "05uc_aim1_missing_cell_types")
message("Aim 1.3 -- signatures by status:")
print(table(miss_tab$status))
if (any(grepl("^MISSING", miss_tab$status)))
  message("Possibly missed cell types: ", paste(miss_tab$cell_type[grepl("^MISSING", miss_tab$status)], collapse = ", "))
hm <- dplyr::bind_rows(miss_hm)
if (nrow(hm)) {
  top_labels <- hm |> group_by(label) |> summarise(v = max(pct)) |> filter(v >= 10) |> pull(label)
  hm <- hm[hm$label %in% top_labels, ]
  hm$signature <- factor(hm$signature, levels = rev(miss_tab$signature))
  hm$label <- factor(hm$label, levels = intersect(names(lab_pal), top_labels))
  st <- stats::setNames(miss_tab$status, miss_tab$signature)
  ylab <- function(x) paste0(x, ifelse(grepl("^MISSING", st[x]), "  [MISSING]", ifelse(grepl("^not", st[x]), "  [rare]", "")))
  p <- ggplot(hm, aes(label, signature, fill = pct)) + geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct >= 15, sprintf("%.0f", pct), "")), size = 2.3) +
    scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = c(0, 100), name = "% of the\nco-expressing\ncells") +
    scale_y_discrete(labels = ylab) + theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank()) +
    labs(x = "label given in 04", y = "marker signature (cells with >= 2 of its genes)",
         title = "Aim 1.3: where do the cells co-expressing each signature end up?",
         subtitle = "a row without a red tile under its own label = cell type merged into other clusters; [MISSING] = no label for it at all")
  save_fig(p, "05uc_aim1_missing_cell_types", width = max(10, 0.32 * length(top_labels) + 6),
           height = max(7, 0.2 * length(sigs) + 3))
}
unclear <- grepl("unknown|unclear|mixed|doublet|^cluster|unlabelled|^other|low_quality", m$label, ignore.case = TRUE)
uc_tab <- m |> mutate(unclear = unclear) |> group_by(core_id, nancy) |>
  summarise(n_cells = n(), pct_unclear = 100 * mean(unclear), .groups = "drop") |> arrange(desc(pct_unclear))
save_table(uc_tab, "05uc_aim1_unclear_labels_per_core")
message("Cells with unclear / unlabelled labels: ", round(100 * mean(unclear), 1), "%")
} ## end AIM 1

## =============================================================================
## Per-core features (Aims 2 and 3)
## =============================================================================
if (any(c("2", "3") %in% AIMS)) {
message("Per-core spatial features ...")
cls_tab <- m |> count(class = ifelse(is.na(class), "(no class)", class), label, name = "n_cells") |> arrange(class, desc(n_cells))
save_table(cls_tab, "05uc_label_to_class")
message("Labels -> feature classes (check tables/05uc_label_to_class.csv; patterns: UC_CLASS_PATTERNS in R/uc_helpers.R):")
print(as.data.frame(cls_tab |> group_by(class) |> summarise(labels = paste(label, collapse = ", "), n_cells = sum(n_cells))),
      row.names = FALSE, right = FALSE)
feat <- uc_core_features(m)

## 2.2 niches: K nearest neighbours' label composition, k-means into N_NICHES
types <- sort(unique(m$label))
nbhd <- matrix(0, nrow(m), length(types), dimnames = list(rownames(m), types))
for (cid in unique(m$core_id)) {
  idx <- which(m$core_id == cid)
  if (length(idx) <= NICHE_K_NEIGHBOURS) next
  nn <- RANN::nn2(m[idx, c("x_um", "y_um")], k = NICHE_K_NEIGHBOURS + 1)$nn.idx[, -1]
  lab <- matrix(match(m$label[idx][nn], types), nrow = length(idx))
  for (t in seq_along(types)) nbhd[idx, t] <- rowMeans(lab == t)
}
km <- stats::kmeans(nbhd, centers = N_NICHES, nstart = 5, iter.max = 50)
cen <- km$centers
niche_names <- apply(cen, 1, function(r) paste(names(sort(r, decreasing = TRUE))[1:2], collapse = " + "))
niche_names <- paste0("N", seq_len(N_NICHES), ": ", niche_names)
m$niche <- factor(niche_names[km$cluster], levels = niche_names)
niche_core <- as.data.frame.matrix(prop.table(table(m$core_id, m$niche), 1) * 100)
colnames(niche_core) <- paste0("Niche_N", seq_len(N_NICHES), "_pct")
feat <- cbind(feat, niche_core[feat$core_id, , drop = FALSE])
save_table(data.frame(niche = niche_names, round(cen, 3), check.names = FALSE), "05uc_aim2_niche_composition")

feat <- cbind(core_tab[feat$core_id, c("slide_id", "patient_id", "nancy", "nancy_group",
                                         intersect(c("Nancy", "Nancy_TMA"), colnames(core_tab)))], feat)
feat_cols <- setdiff(colnames(feat), c("slide_id", "patient_id", "nancy", "nancy_group", "Nancy", "Nancy_TMA",
                                       "core_id", "n_cells", "area_mm2"))
feat_cols <- feat_cols[vapply(feat[feat_cols], function(v) mean(!is.na(v)) >= 0.5 && stats::sd(v, na.rm = TRUE) > 0, logical(1))]
if (!is.null(clin)) feat <- cbind(feat, clin[feat$core_id, setdiff(colnames(clin), "core_id"), drop = FALSE])
rownames(feat) <- feat$core_id
save_table(feat, "05uc_core_features")
saveRDS(list(features = feat, feature_cols = feat_cols, niche = stats::setNames(as.character(m$niche), rownames(m))),
        file.path(PROCESSED_DIR, "05uc_core_features.rds"))
message(length(feat_cols), " features per core -> tables/05uc_core_features.csv")
key_feats <- intersect(c("Neutrophil_per_mm2", "Neutrophil_intraepithelial_pct", "Neutrophil_aggregated_pct",
                         "Eosinophil_per_mm2", "Plasma_per_mm2", "Plasma_near_crypt_base_pct", "Immune_pct",
                         "Inflamed_epi_pct_epithelium", "Regenerative_pct_epithelium", "Goblet_ratio",
                         "Epithelium_pct", "Inflam_fibroblast_pct_stroma"), feat_cols)
grade_lab <- function(n) factor(ifelse(is.na(n), "NA", paste("Nancy", n)), levels = c(paste("Nancy", 0:4), "NA"))
feat$grade <- grade_lab(feat$nancy)
grade_cols <- c(stats::setNames(NANCY_COLOURS, paste("Nancy", 0:4)), "NA" = "grey60")
}

## =============================================================================
## AIM 2
## =============================================================================
if ("2" %in% AIMS) {
## 2.1 composition per core (cell types), cores grouped by Nancy grade
comp <- m |> count(core_id, label) |> group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
comp$grade <- feat[comp$core_id, "grade"]
core_ord <- feat |> arrange(grade, desc(Epithelium_pct)) |> pull(core_id)
comp$core_id <- factor(comp$core_id, levels = core_ord)
comp$label <- factor(comp$label, levels = names(lab_pal))
save_table(comp, "05uc_aim2_composition_per_core")
p <- ggplot(comp, aes(core_id, frac, fill = label)) + geom_col(width = 1) +
  facet_grid(~ grade, scales = "free_x", space = "free_x") + scale_fill_manual(values = lab_pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) + theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), panel.grid = element_blank()) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "core", y = "cells", title = "Aim 2: cell type composition per core, by Nancy grade")
save_fig(p, "05uc_aim2_composition_per_core", width = max(14, 0.12 * nrow(feat) + 6), height = 8)

## 2.2 niches
nc <- m |> count(core_id, niche) |> group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
nc$grade <- feat[nc$core_id, "grade"]; nc$core_id <- factor(nc$core_id, levels = core_ord)
niche_pal <- stats::setNames(grDevices::hcl.colors(N_NICHES, "Dark 3"), niche_names)
p1 <- ggplot(nc, aes(core_id, frac, fill = niche)) + geom_col(width = 1) +
  facet_grid(~ grade, scales = "free_x", space = "free_x") + scale_fill_manual(values = niche_pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) + theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), panel.grid = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2)) + labs(x = "core", y = "cells", title = "Niches per core, by Nancy grade")
cd <- as.data.frame(as.table(cen)); colnames(cd) <- c("niche", "label", "frac")
cd$niche <- factor(niche_names[as.integer(cd$niche)], levels = rev(niche_names))
cd <- cd[cd$label %in% names(which(apply(cen, 2, max) >= 0.05)), ]
p2 <- ggplot(cd, aes(label, niche, fill = frac)) + geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#2D004B", labels = scales::percent, name = "mean share of\nneighbours") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = sprintf("Niche = composition of a cell's %d nearest neighbours", NICHE_K_NEIGHBOURS))
save_fig(p2 / p1 + plot_layout(heights = c(1, 1.3)), "05uc_aim2_niches", width = max(14, 0.12 * nrow(feat) + 6), height = 12)

## 2.3 feature heatmap: z-scores, cores ordered by grade, clustered within grade
Z <- scale(as.matrix(feat[, feat_cols]))
Z[is.na(Z)] <- 0; Z <- pmax(pmin(Z, 3), -3)
ord <- unlist(lapply(split(feat$core_id, feat$grade), function(cc) {
  if (length(cc) < 3) return(cc)
  cc[stats::hclust(stats::dist(Z[cc, , drop = FALSE]))$order]
}))
fo <- if (ncol(Z) > 2) colnames(Z)[stats::hclust(stats::dist(t(Z)))$order] else colnames(Z)
zd <- data.frame(core_id = rep(rownames(Z), ncol(Z)), feature = rep(colnames(Z), each = nrow(Z)), z = as.vector(Z))
zd$core_id <- factor(zd$core_id, levels = ord); zd$feature <- factor(zd$feature, levels = fo)
zd$grade <- feat[as.character(zd$core_id), "grade"]
p <- ggplot(zd, aes(core_id, feature, fill = z)) + geom_tile() +
  facet_grid(~ grade, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3, 3), name = "z-score") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), panel.grid = element_blank()) +
  labs(x = "core", y = NULL, title = "Aim 2: per-core features, cores ordered by Nancy grade (clustered within grade)",
       subtitle = "columns that look different inside one grade = heterogeneity within the same Nancy score")
save_fig(p, "05uc_aim2_feature_heatmap", width = max(14, 0.12 * nrow(feat) + 7), height = max(8, 0.17 * length(feat_cols) + 3))

## 2.4 spread within each grade
if (length(key_feats)) {
  kd <- tidyr::pivot_longer(feat[, c("core_id", "patient_id", "grade", key_feats)], dplyr::all_of(key_feats),
                            names_to = "feature", values_to = "value")
  kd$feature <- factor(kd$feature, levels = key_feats)
  p <- ggplot(kd[!is.na(kd$value), ], aes(grade, value)) +
    geom_boxplot(aes(fill = grade), outlier.shape = NA, alpha = 0.35, colour = "grey40") +
    geom_point(aes(colour = grade), position = position_jitter(width = 0.15, seed = 1), size = 1.3) +
    scale_fill_manual(values = grade_cols, guide = "none") + scale_colour_manual(values = grade_cols, guide = "none") +
    facet_wrap(~ feature, scales = "free_y", ncol = 4) + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = NULL, title = "Aim 2: key UC features per core, by Nancy grade (one dot = one core)",
         subtitle = "spread inside a box = heterogeneity within the same Nancy score")
  save_fig(p, "05uc_aim2_key_features_by_nancy", width = 15, height = 3.2 * ceiling(length(key_feats) / 4) + 1.5)
}
Zf <- scale(as.matrix(feat[, feat_cols])); Zf[is.na(Zf)] <- 0
pc <- stats::prcomp(Zf)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pcd <- data.frame(core_id = rownames(Zf), PC1 = pc$x[, 1], PC2 = pc$x[, 2], grade = feat$grade, patient_id = feat$patient_id)
cent <- t(sapply(split(seq_len(nrow(Zf)), feat$grade), function(i) colMeans(Zf[i, , drop = FALSE])))
dist_c <- sqrt(rowSums((Zf - cent[as.character(feat$grade), , drop = FALSE])^2))
dev_feats <- apply(Zf - cent[as.character(feat$grade), , drop = FALSE], 1, function(r)
  paste(sprintf("%s %+.1f", names(r)[order(-abs(r))][1:3], r[order(-abs(r))][1:3]), collapse = "; "))
het <- data.frame(core_id = rownames(Zf), patient_id = feat$patient_id, nancy = feat$nancy,
                  distance_to_grade_centre = round(dist_c, 2), most_deviating_features_z = dev_feats)
het <- het[order(het$nancy, -het$distance_to_grade_centre), ]
save_table(het, "05uc_aim2_within_grade_heterogeneity")
lab_out <- pcd[rank(-dist_c) <= 10, ]
p_pca <- ggplot(pcd, aes(PC1, PC2, colour = grade)) + geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = grade_cols, name = NULL) + theme_bw(base_size = 12) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
       title = "Cores in feature space", subtitle = "labelled: 10 cores furthest from their grade's centre")
if (requireNamespace("ggrepel", quietly = TRUE))
  p_pca <- p_pca + ggrepel::geom_text_repel(data = lab_out, aes(label = core_id), size = 3, show.legend = FALSE)
het$grade <- grade_lab(het$nancy)
p_d <- ggplot(het, aes(grade, distance_to_grade_centre, colour = grade)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") + geom_point(position = position_jitter(width = 0.15, seed = 1)) +
  scale_colour_manual(values = grade_cols, guide = "none") + theme_bw(base_size = 12) +
  labs(x = NULL, y = "distance to the grade centre (z units)", title = "Heterogeneity within each grade")
pp <- p_pca + p_d + plot_layout(widths = c(1.3, 1))
## same patient vs different patients, same grade
pairs <- which(upper.tri(diag(nrow(Zf))), arr.ind = TRUE)
pairs <- pairs[feat$grade[pairs[, 1]] == feat$grade[pairs[, 2]], , drop = FALSE]
if (nrow(pairs) && any(feat$patient_id[pairs[, 1]] == feat$patient_id[pairs[, 2]])) {
  pdist <- data.frame(grade = feat$grade[pairs[, 1]],
                      relation = ifelse(feat$patient_id[pairs[, 1]] == feat$patient_id[pairs[, 2]], "same patient", "different patients"),
                      distance = sqrt(rowSums((Zf[pairs[, 1], , drop = FALSE] - Zf[pairs[, 2], , drop = FALSE])^2)))
  save_table(pdist |> group_by(grade, relation) |> summarise(n_pairs = n(), median_distance = stats::median(distance), .groups = "drop"),
             "05uc_aim2_same_vs_different_patient")
  pp <- pp + (ggplot(pdist, aes(grade, distance, fill = relation)) + geom_boxplot(outlier.size = 0.4) +
                scale_fill_manual(values = c("same patient" = "#F5C26B", "different patients" = "grey70"), name = NULL) +
                theme_bw(base_size = 12) + theme(legend.position = "top") +
                labs(x = NULL, y = "distance between two cores of the same grade", title = "Same patient vs different patients")) +
    plot_layout(widths = c(1.3, 1, 1.1))
}
save_fig(pp, "05uc_aim2_within_grade_spread", width = 19, height = 6.5)
} ## end AIM 2

## =============================================================================
## AIM 3
## =============================================================================
if ("3" %in% AIMS) {
## 3.1 agreement between Nancy scores
scores <- data.frame(core_id = core_tab$core_id, row.names = core_tab$core_id)
for (c in intersect(c("Nancy_TMA", "Nancy"), colnames(core_tab))) scores[[c]] <- as_num(core_tab[[c]])
if (!is.null(clin)) for (c in intersect(UC_CLIN_NANCY_COLS, colnames(clin))) scores[[c]] <- as_num(clin[core_tab$core_id, c])
sc <- setdiff(colnames(scores), "core_id")
sc <- sc[vapply(scores[sc], function(v) sum(!is.na(v)) >= 3, logical(1))]
if (length(sc) >= 2) {
  prs <- utils::combn(sc, 2)
  agree <- dplyr::bind_rows(lapply(seq_len(ncol(prs)), function(k) {
    a <- scores[[prs[1, k]]]; b <- scores[[prs[2, k]]]; ok <- !is.na(a) & !is.na(b)
    data.frame(score_a = prs[1, k], score_b = prs[2, k], n_cores = sum(ok), exact_agreement_pct = 100 * mean(a[ok] == b[ok]),
               within_1_pct = 100 * mean(abs(a[ok] - b[ok]) <= 1), weighted_kappa = wkappa(a, b))
  }))
  save_table(agree, "05uc_aim3_nancy_agreement")
  print(agree, row.names = FALSE)
  cf <- dplyr::bind_rows(lapply(seq_len(ncol(prs)), function(k) {
    d <- data.frame(a = scores[[prs[1, k]]], b = scores[[prs[2, k]]]); d <- d[stats::complete.cases(d), ]
    if (!nrow(d)) return(NULL)
    d |> count(a, b) |> mutate(pair = sprintf("%s vs %s\n(kappa_w = %.2f)", prs[1, k], prs[2, k], agree$weighted_kappa[k]))
  }))
  p <- ggplot(cf, aes(factor(b), factor(a), fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = n)) +
    scale_fill_gradient(low = "#F2F2F2", high = "#3B6FB6") + facet_wrap(~ pair, scales = "free") +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank()) +
    labs(x = "second score", y = "first score", title = "Aim 3.1: agreement between Nancy scores (cores)")
  save_fig(p, "05uc_aim3_nancy_agreement", width = 4.5 * min(3, ncol(prs)) + 1, height = 4.2 * ceiling(ncol(prs) / 3) + 1)
}

## 3.2 each feature vs Nancy
if (any(!is.na(feat$nancy))) {
  an <- assoc_nancy(feat, feat_cols)
  save_table(an, "05uc_aim3_features_vs_nancy")
  if (nrow(an)) {
    an$sig <- ifelse(an$padj_spearman < 0.05, "FDR < 5%", "n.s.")
    an$feature <- factor(an$feature, levels = an$feature[order(an$spearman_rho)])
    p1 <- ggplot(an, aes(spearman_rho, feature, colour = sig)) + geom_vline(xintercept = 0, colour = "grey60") +
      geom_segment(aes(x = 0, xend = spearman_rho, yend = feature), linewidth = 0.4) + geom_point(size = 2) +
      scale_colour_manual(values = c("FDR < 5%" = "#B2182B", "n.s." = "grey60"), name = NULL) +
      theme_bw(base_size = 9) + theme(legend.position = "top") +
      labs(x = "Spearman rho with Nancy (cores)", y = NULL, title = "Aim 3.2: which features follow the Nancy score?")
    top <- utils::head(as.character(an$feature[order(an$p_spearman)]), 8)
    td <- tidyr::pivot_longer(feat[!is.na(feat$nancy), c("core_id", "grade", top)], dplyr::all_of(top),
                              names_to = "feature", values_to = "value")
    td$feature <- factor(td$feature, levels = top)
    p2 <- ggplot(td[!is.na(td$value), ], aes(grade, value, colour = grade)) +
      geom_boxplot(outlier.shape = NA, colour = "grey50") + geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1) +
      scale_colour_manual(values = grade_cols, guide = "none") + facet_wrap(~ feature, scales = "free_y", ncol = 4) +
      theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = NULL, y = NULL, title = "Strongest 8")
    save_fig(p1 + p2 + plot_layout(widths = c(1, 2)), "05uc_aim3_features_vs_nancy", width = 20,
             height = max(8, 0.16 * nrow(an) + 2))
  }

  ## 3.3 ordinal model, leave one patient out
  mf <- intersect(UC_MODEL_FEATURES, feat_cols)
  d <- feat[!is.na(feat$nancy), c("core_id", "patient_id", "nancy", mf)]
  if (length(mf) >= 2 && length(unique(d$patient_id)) >= 6 && length(unique(d$nancy)) >= 2) {
    for (v in mf) if (grepl("_per_mm2$", v)) d[[v]] <- log1p(d[[v]])
    pred <- rep(NA_real_, nrow(d)); used <- character(nrow(d))
    for (p_out in unique(d$patient_id)) {
      tr <- d$patient_id != p_out; te <- !tr
      X <- d[, mf, drop = FALSE]
      med <- vapply(X[tr, , drop = FALSE], stats::median, numeric(1), na.rm = TRUE)
      for (v in mf) X[[v]][is.na(X[[v]])] <- med[[v]]
      mu <- colMeans(X[tr, , drop = FALSE]); s <- apply(X[tr, , drop = FALSE], 2, stats::sd); s[s == 0 | is.na(s)] <- 1
      X <- as.data.frame(scale(X, mu, s))
      dtr <- cbind(y = factor(d$nancy[tr], ordered = TRUE), X[tr, , drop = FALSE])
      fit <- tryCatch(suppressWarnings(MASS::polr(y ~ ., data = dtr, Hess = FALSE)), error = function(e) NULL)
      if (!is.null(fit)) {
        pr <- stats::predict(fit, newdata = X[te, , drop = FALSE], type = "probs")
        pr <- matrix(pr, nrow = sum(te))
        pred[te] <- as.vector(pr %*% as.numeric(fit$lev)); used[te] <- "ordinal"
      } else {
        dtr$y <- d$nancy[tr]
        fit <- stats::lm(y ~ ., data = dtr)
        pred[te] <- suppressWarnings(stats::predict(fit, newdata = X[te, , drop = FALSE])); used[te] <- "linear"
      }
    }
    lv <- sort(unique(d$nancy))
    d$spatial_nancy <- round(pred, 2)
    d$predicted_grade <- lv[pmax(1, pmin(length(lv), round(stats::approx(lv, seq_along(lv), pred, rule = 2)$y)))]
    d$model <- used
    d$difference <- d$spatial_nancy - d$nancy
    save_table(d[order(-abs(d$difference)), ], "05uc_aim3_spatial_nancy_per_core")
    metr <- data.frame(n_cores = nrow(d), n_patients = length(unique(d$patient_id)), features = paste(mf, collapse = ";"),
                       exact_pct = 100 * mean(d$predicted_grade == d$nancy),
                       within_1_pct = 100 * mean(abs(d$predicted_grade - d$nancy) <= 1),
                       weighted_kappa = wkappa(d$predicted_grade, d$nancy),
                       spearman = stats::cor(d$spatial_nancy, d$nancy, method = "spearman"))
    save_table(metr, "05uc_aim3_spatial_nancy_performance")
    message("Aim 3.3 spatial Nancy (leave one patient out): exact ", round(metr$exact_pct), "%, within 1 grade ",
            round(metr$within_1_pct), "%, weighted kappa ", round(metr$weighted_kappa, 2))
    d$grade <- grade_lab(d$nancy)
    disc <- d[abs(d$difference) >= 1.5, ]
    p1 <- ggplot(d, aes(grade, spatial_nancy, colour = grade)) +
      geom_boxplot(outlier.shape = NA, colour = "grey50") + geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1.6) +
      geom_abline(intercept = -1, slope = 1, linetype = 2, colour = "grey60") +
      scale_colour_manual(values = grade_cols, guide = "none") + theme_bw(base_size = 12) +
      labs(x = paste("pathologist:", nancy_col), y = "spatial Nancy (model, patient left out)",
           title = "Aim 3.3: spatial Nancy vs the pathologist's grade",
           subtitle = sprintf("exact %.0f%%, within 1 grade %.0f%%, weighted kappa %.2f; labelled = |difference| >= 1.5",
                              metr$exact_pct, metr$within_1_pct, metr$weighted_kappa))
    if (nrow(disc) && requireNamespace("ggrepel", quietly = TRUE))
      p1 <- p1 + ggrepel::geom_text_repel(data = disc, aes(label = core_id), size = 2.8, colour = "black", max.overlaps = 30)
    cm <- d |> count(nancy, predicted_grade)
    p2 <- ggplot(cm, aes(factor(predicted_grade), factor(nancy), fill = n)) + geom_tile(colour = "white") +
      geom_text(aes(label = n)) + scale_fill_gradient(low = "#F2F2F2", high = "#3B6FB6", guide = "none") +
      theme_bw(base_size = 12) + labs(x = "predicted grade", y = "pathologist's grade", title = "Confusion (cores)")
    save_fig(p1 + p2 + plot_layout(widths = c(1.4, 1)), "05uc_aim3_spatial_nancy", width = 14, height = 6)
    feat$spatial_nancy <- d$spatial_nancy[match(feat$core_id, d$core_id)]
  } else message("Aim 3.3 skipped: needs >= 2 model features, >= 6 patients and >= 2 Nancy grades.")
}

## 3.4 patient level: outcomes (Disease.Control, remission) vs features
oc <- if (is.null(clin)) character(0) else intersect(UC_OUTCOME_COLS, colnames(feat))
if (length(oc)) {
  pcols <- c(feat_cols, intersect("spatial_nancy", colnames(feat)))
  pat <- feat |> group_by(patient_id) |>
    summarise(n_cores = n(), nancy_max = suppressWarnings(max(nancy, na.rm = TRUE)), nancy_mean = mean(nancy, na.rm = TRUE),
              dplyr::across(dplyr::all_of(pcols), ~ mean(.x, na.rm = TRUE)),
              dplyr::across(dplyr::all_of(oc), ~ dplyr::first(stats::na.omit(.x))), .groups = "drop") |> as.data.frame()
  pat$nancy_max[!is.finite(pat$nancy_max)] <- NA
  for (o in oc) pat[[o]] <- as_outcome(pat[[o]])
  save_table(pat, "05uc_aim3_patient_table")
  res <- dplyr::bind_rows(lapply(oc, function(o) assoc_outcome(pat, c("nancy_max", pcols), o)))
  save_table(res, "05uc_aim3_features_vs_outcomes")
  ## does a feature add to Nancy? leave-one-patient-out AUC of outcome ~ Nancy (+ feature)
  cv_auc <- function(y, X) {
    pr <- vapply(seq_along(y), function(i) {
      fit <- suppressWarnings(stats::glm(y[-i] ~ ., data = X[-i, , drop = FALSE], family = stats::binomial()))
      suppressWarnings(stats::predict(fit, newdata = X[i, , drop = FALSE], type = "response"))
    }, numeric(1))
    auc_mw(pr[y == 1], pr[y == 0])
  }
  add_tab <- list()
  for (o in oc) {
    y0 <- pat[[o]]
    if (nlevels(droplevels(y0)) != 2) next
    ok <- !is.na(y0) & !is.na(pat$nancy_max)
    y <- as.integer(droplevels(y0[ok])) - 1
    if (min(table(y)) < 4) next
    base <- cv_auc(y, data.frame(nancy = pat$nancy_max[ok]))
    cand <- res$feature[res$outcome == o & res$feature != "nancy_max"]
    cand <- utils::head(cand[order(res$p[res$outcome == o & res$feature != "nancy_max"])], 10)
    for (f in cand) {
      x <- pat[[f]][ok]; x[is.na(x)] <- stats::median(x, na.rm = TRUE)
      add_tab[[paste(o, f)]] <- data.frame(outcome = o, feature = f, n_patients = length(y),
                                           cv_auc_nancy = base, cv_auc_nancy_plus_feature = cv_auc(y, data.frame(nancy = pat$nancy_max[ok], f = x)))
    }
  }
  add_tab <- dplyr::bind_rows(add_tab)
  if (nrow(add_tab)) {
    add_tab$gain <- add_tab$cv_auc_nancy_plus_feature - add_tab$cv_auc_nancy
    save_table(add_tab, "05uc_aim3_added_value_over_nancy")
  }
  main_o <- if ("Disease.Control" %in% oc) "Disease.Control" else oc[1]
  r <- res[res$outcome == main_o, ]
  if (nrow(r)) {
    top <- utils::head(r$feature[order(r$p)], 8)
    td <- tidyr::pivot_longer(pat[!is.na(pat[[main_o]]), c("patient_id", "nancy_max", main_o, top)], dplyr::all_of(top),
                              names_to = "feature", values_to = "value")
    td$feature <- factor(td$feature, levels = top)
    td$nancy_group <- nancy_group(td$nancy_max)
    lab <- r[match(top, r$feature), ]
    lab$feature <- factor(lab$feature, levels = top)
    p1 <- ggplot(td[!is.na(td$value), ], aes(.data[[main_o]], value)) +
      geom_boxplot(outlier.shape = NA, colour = "grey45") +
      geom_point(aes(colour = nancy_group), position = position_jitter(width = 0.15, seed = 1), size = 1.8) +
      geom_text(data = lab, aes(x = -Inf, y = Inf, label = sprintf("AUC %.2f, p %.3g", auc, p)), hjust = -0.05,
                vjust = 1.3, size = 3, inherit.aes = FALSE) +
      scale_colour_manual(values = GROUP_COLOURS(), name = "patient's highest\nNancy", na.value = "grey60") +
      facet_wrap(~ feature, scales = "free_y", ncol = 4) + theme_bw(base_size = 11) +
      labs(x = main_o, y = "patient mean over cores",
           title = sprintf("Aim 3.4: features vs %s (one dot = one patient)", main_o),
           subtitle = "AUC = chance that a random 'yes' patient has a higher value than a random 'no' patient (0.5 = no difference)")
    save_fig(p1, paste0("05uc_aim3_features_vs_", make.names(main_o)), width = 15, height = 7.5)
  }
  if (nrow(add_tab)) {
    ad <- add_tab[add_tab$outcome == main_o, ]
    if (nrow(ad)) {
      ad$feature <- factor(ad$feature, levels = ad$feature[order(ad$cv_auc_nancy_plus_feature)])
      p2 <- ggplot(ad, aes(cv_auc_nancy_plus_feature, feature)) +
        geom_vline(aes(xintercept = cv_auc_nancy), linetype = 2, colour = "#3B6FB6") +
        geom_vline(xintercept = 0.5, colour = "grey70") +
        geom_segment(aes(x = cv_auc_nancy, xend = cv_auc_nancy_plus_feature, yend = feature), colour = "grey60") +
        geom_point(size = 2.5, colour = "#C0392B") + theme_bw(base_size = 11) +
        labs(x = "leave-one-patient-out AUC", y = NULL,
             title = sprintf("Aim 3.4: does a spatial feature add to Nancy for %s?", main_o),
             subtitle = "blue dashed = Nancy alone (patient's highest); red = Nancy + feature")
      save_fig(p2, paste0("05uc_aim3_added_value_", make.names(main_o)), width = 9, height = 5.5)
    }
  }
} else message("Aim 3.4 skipped: no clinical outcomes joined (UC_CLINICAL_CSV, columns ",
               paste(UC_OUTCOME_COLS, collapse = ", "), ").")
} ## end AIM 3
message("Done.")
