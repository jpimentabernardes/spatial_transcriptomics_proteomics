## =============================================================================
## 05_tma_first_insights.R
##
## First insights from the multi-organ / multi-tumour TMAs, organised by aim.
##
## AIM 1 -- probe specificity: are the colon + custom probes gut-specific, or
##          do the immune targets generalise to other organs?
##   1.1 Signal over background per gene x core (vs negative-control probes).
##       "Detected" = clearly above the non-binding-probe rate in that core.
##   1.2 Organ specificity (tau) -- computed twice: on all cells AND within
##       immune cells only. An immune probe with high tau on all cells but low
##       tau within immune cells is NOT organ-specific; organs just differ in
##       how many immune cells they contain. Separating these two is the key
##       step for answering Aim 1 honestly.
##   1.3 Lineage specificity: does an immune probe's signal land in immune
##       cells in every organ, or does it appear in epithelium/stroma in some
##       organs (non-specific binding, or segmentation spill-over -> compare
##       with MECR in 03c)?
##   1.4 Generalisability: is the immune-target profile of e.g. T cells in
##       lung the same as in colon (correlation between organs per cell type)?
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
## INPUT : xenium_<seg>_annotated.rds (04), config/probe_expectations.csv
## OUTPUT: tables/05_*.csv, figures/05_*.png, xenium_<seg>_insights.rds
## =============================================================================

source("config/config.R")
source("R/helpers.R")

xen <- readRDS(obj_path("annotated"))
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")
meta <- xen[[]]
genes <- rownames(counts)

probes <- read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE)
probes$gene <- as_seurat_features(probes$gene)
missing_probes <- setdiff(probes$gene, genes)
if (length(missing_probes) > 0) message("In probe_expectations.csv but not on panel: ", paste(missing_probes, collapse = ", "))
probes <- probes[probes$gene %in% genes, ]
immune_targets <- probes$gene[probes$category == "immune_target"]
custom_genes <- probes$gene[probes$panel == "custom"]

core_info <- meta |>
  distinct(core_id, slide_id, patient_id, organ, tissue_type, diagnosis, location) |>
  as.data.frame()
rownames(core_info) <- core_info$core_id

## =============================================================================
## AIM 1
## =============================================================================

## -----------------------------------------------------------------------
## 1.1 Signal-to-background per gene x core
##     lambda = counts per negative-control probe in this core = what a
##     non-binding probe collects here (autofluorescence, off-target
##     hybridisation, mis-decoding). Background differs a lot between organs
##     (e.g. autofluorescent lung/liver), so it must be per core, not global.
##     CAVEAT: gene probe sets and control probes are not identical in probe
##     number/design, so treat SBR as a ranking, not an absolute FDR.
## -----------------------------------------------------------------------
ctrl_col <- if ("nCount_ControlProbe" %in% colnames(meta)) "ControlProbe" else "ControlCodeword"
message("Background model uses: ", ctrl_col)
pb_core <- pseudobulk(counts, meta$core_id)                     # genes x cores
bg_core <- meta |>
  group_by(core_id) |>
  summarise(ctrl = sum(.data[[paste0("nCount_", ctrl_col)]]),
            n_ctrl = first(.data[[paste0("nFeatTotal_", ctrl_col)]]),
            n_cells = n(), .groups = "drop") |>
  mutate(lambda = ctrl / n_ctrl)
lambda <- stats::setNames(bg_core$lambda, bg_core$core_id)[colnames(pb_core)]

sbr <- as.data.frame(as.table(as.matrix(pb_core)), stringsAsFactors = FALSE)
colnames(sbr) <- c("gene", "core_id", "counts")
sbr$lambda <- lambda[sbr$core_id]
sbr$log2_sbr <- log2((sbr$counts + 0.5) / (sbr$lambda + 0.5))
sbr$p <- stats::ppois(sbr$counts - 1, sbr$lambda, lower.tail = FALSE)
sbr$padj <- stats::p.adjust(sbr$p, method = "BH")
sbr$detected <- sbr$log2_sbr >= SBR_LOG2_MIN & sbr$padj < SBR_PADJ_MAX
sbr <- cbind(sbr, core_info[sbr$core_id, c("slide_id", "patient_id", "organ", "tissue_type")])
sbr$panel <- ifelse(sbr$gene %in% custom_genes, "custom", "colon")
save_table(sbr, "05_aim1_sbr_gene_by_core")

sbr_organ <- sbr |>
  group_by(gene, panel, organ) |>
  summarise(frac_cores_detected = mean(detected), median_log2_sbr = median(log2_sbr),
            n_cores = n(), .groups = "drop")
save_table(sbr_organ, "05_aim1_sbr_gene_by_organ")

## Panel-wide view: how many genes are detected above background in each organ?
p_det <- sbr_organ |>
  group_by(organ, panel) |>
  summarise(pct_genes_detected = 100 * mean(frac_cores_detected >= 0.5), .groups = "drop") |>
  ggplot(aes(organ, pct_genes_detected, fill = panel)) +
  geom_col(position = "dodge") + theme_bw() +
  labs(y = "% of panel genes above background in >= half the cores",
       title = "Aim 1: fraction of the panel that 'works' per organ")
save_fig(p_det, "05_aim1_pct_genes_detected_by_organ", width = 10, height = 5)

focus <- unique(c(immune_targets, custom_genes))
if (length(focus) > 0) {
  hm <- sbr_organ |> filter(gene %in% focus)
  p_hm <- ggplot(hm, aes(organ, gene, fill = median_log2_sbr)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.0f%%", 100 * frac_cores_detected)), size = 2.3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = SBR_LOG2_MIN) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") + theme_bw() +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(title = "Aim 1: immune + custom probes, median log2 signal/background per organ",
         subtitle = "label = % of cores in which the probe is detected above background")
  save_fig(p_hm, "05_aim1_sbr_heatmap_focus_genes", width = 10, height = max(5, 0.25 * length(focus) + 2))
}

## -----------------------------------------------------------------------
## 1.2 Organ specificity (tau): all cells vs within immune cells.
##     Each core weighted equally: CPM per core, then mean over the organ's
##     cores (otherwise one huge core dominates its organ).
## -----------------------------------------------------------------------
organ_mean_cpm <- function(cells) {
  pb <- pseudobulk(counts[, cells, drop = FALSE], meta[cells, "core_id"])
  cpm <- sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e6
  org <- core_info[colnames(cpm), "organ"]
  sapply(split(colnames(cpm), org), function(cc) rowMeans(cpm[, cc, drop = FALSE]))
}
cpm_all <- organ_mean_cpm(colnames(xen))
imm_cells <- colnames(xen)[meta$lineage == "Immune"]
cpm_imm <- if (length(imm_cells) > 100) organ_mean_cpm(imm_cells) else NULL

tau_tab <- data.frame(gene = genes, tau_all = apply(cpm_all, 1, tau_index))
tau_tab$top_organ_all <- colnames(cpm_all)[apply(cpm_all, 1, which.max)]
if (!is.null(cpm_imm)) {
  tau_tab$tau_immune <- apply(cpm_imm, 1, tau_index)
  tau_tab$top_organ_immune <- colnames(cpm_imm)[apply(cpm_imm, 1, which.max)]
}
tau_tab$panel <- ifelse(tau_tab$gene %in% custom_genes, "custom", "colon")
tau_tab$immune_target <- tau_tab$gene %in% immune_targets
save_table(tau_tab, "05_aim1_organ_specificity_tau")

if (!is.null(cpm_imm)) {
  lab <-tau_tab[tau_tab$immune_target | tau_tab$panel == "custom", ]
  p_tau <- ggplot(tau_tab, aes(tau_all, tau_immune)) +
    geom_abline(linetype = 2, color = "grey50") +
    geom_point(aes(color = panel, shape = immune_target), alpha = 0.7) +
    geom_text(data = lab, aes(label = gene), size = 2.5, vjust = -0.6) +
    theme_bw() +
    labs(x = "tau across organs, all cells", y = "tau across organs, immune cells only",
         title = "Aim 1: organ specificity -- below the diagonal = driven by immune cell abundance, not the probe")
  save_fig(p_tau, "05_aim1_tau_all_vs_immune", width = 10, height = 8)
}

## -----------------------------------------------------------------------
## 1.3 Lineage specificity per probe x organ.
##     enrichment = (share of the gene's counts in its expected lineage) /
##                  (share of ALL transcripts in that lineage). 1 = no
##                  specificity; >>1 = on target.
##     offtarget_vs_bg = mean counts per NON-expected-lineage cell / what a
##                  negative-control probe collects per cell. ~1 = the
##                  off-target signal is just background; >>1 = real signal
##                  (or spill-over) outside the expected cells.
## -----------------------------------------------------------------------
lin_probes <- probes[!is.na(probes$expected_lineage) & nzchar(probes$expected_lineage), ]
if (nrow(lin_probes) > 0) {
  tot_by_cell <- meta$nCount_Xenium
  ctrl_per_cell <- meta[[paste0("nCount_", ctrl_col)]] / meta[[paste0("nFeatTotal_", ctrl_col)]]
  lin_spec <- dplyr::bind_rows(lapply(split(seq_len(nrow(meta)), meta$organ), function(idx) {
    org <- meta$organ[idx[1]]
    dplyr::bind_rows(lapply(seq_len(nrow(lin_probes)), function(j) {
      g <- lin_probes$gene[j]
      on <- meta$lineage[idx] == lin_probes$expected_lineage[j]
      gc <- counts[g, idx]
      data.frame(
        organ = org, gene = g, panel = lin_probes$panel[j], expected_lineage = lin_probes$expected_lineage[j],
        n_expected_cells = sum(on),
        enrichment = (sum(gc[on]) / max(sum(gc), 1)) / (sum(tot_by_cell[idx][on]) / sum(tot_by_cell[idx])),
        mean_expected = mean(gc[on]),
        mean_offtarget = mean(gc[!on]),
        offtarget_vs_bg = mean(gc[!on]) / max(mean(ctrl_per_cell[idx][!on]), 1e-6)
      )
    }))
  }))
  save_table(lin_spec, "05_aim1_lineage_specificity")
  p_lin <- ggplot(lin_spec, aes(organ, gene, fill = log2(enrichment))) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f", offtarget_vs_bg)), size = 2.2) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
    facet_grid(expected_lineage ~ ., scales = "free_y", space = "free_y") + theme_bw() +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(title = "Aim 1: on-target enrichment of each probe in its expected lineage, per organ",
         subtitle = "label = off-target signal / negative-control background (per cell)")
  save_fig(p_lin, "05_aim1_lineage_specificity", width = 10, height = max(5, 0.25 * nrow(lin_probes) + 2))
}

## -----------------------------------------------------------------------
## 1.4 Generalisability of immune-target profiles across organs.
##     For each immune cell type, pseudobulk per organ over the immune
##     targets, then correlate organs. High correlation = the probes read out
##     the same cell type the same way everywhere.
## -----------------------------------------------------------------------
gen_genes <- if (length(immune_targets) >= 3) immune_targets else genes
imm_types <- setdiff(unique(meta$cell_type[meta$lineage == "Immune"]), NA)
gen_tab <- dplyr::bind_rows(lapply(imm_types, function(ct) {
  cells <- colnames(xen)[meta$cell_type == ct]
  org <- meta[cells, "organ"]
  ok_org <- names(which(table(org) >= MIN_CELLS_PSEUDOBULK))
  if (length(ok_org) < 2) return(NULL)
  cells <- cells[org %in% ok_org]
  pb <- pseudobulk(counts[gen_genes, cells, drop = FALSE], meta[cells, "organ"])
  lcpm <- log1p(sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e4)
  cm <- stats::cor(lcpm, method = "spearman")
  pr <- which(upper.tri(cm), arr.ind = TRUE)
  data.frame(cell_type = ct, organ_a = rownames(cm)[pr[, 1]], organ_b = colnames(cm)[pr[, 2]],
             spearman = cm[pr])
}))
if (nrow(gen_tab) > 0) {
  save_table(gen_tab, "05_aim1_generalisability_immune_profiles")
  p_gen <- ggplot(gen_tab, aes(cell_type, spearman)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Spearman r between organs (immune-target pseudobulk)",
         title = "Aim 1: do immune cell types look the same across organs through this panel?")
  save_fig(p_gen, "05_aim1_generalisability", width = 10, height = 6)
}

if (length(immune_targets) > 0 && length(imm_cells) > 0) {
  sub_imm <- CreateSeuratObject(counts = counts[, imm_cells], assay = "Xenium", meta.data = meta[imm_cells, ])
  sub_imm <- NormalizeData(sub_imm, scale.factor = median(xen$nCount_Xenium), verbose = FALSE)
  sub_imm$type_organ <- paste(sub_imm$cell_type, sub_imm$organ, sep = " | ")
  p_dot_imm <- DotPlot(sub_imm, features = immune_targets, group.by = "type_organ") +
    RotatedAxis() + ggtitle("Aim 1: immune targets within immune cell types, split by organ")
  save_fig(p_dot_imm, "05_aim1_immune_targets_dotplot_by_organ", width = 14,
           height = max(6, 0.18 * length(unique(sub_imm$type_organ))))
  rm(sub_imm)
}

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
xen$niche <- paste0("N", km$cluster)
meta$niche <- xen$niche
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
