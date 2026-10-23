## =============================================================================
## 04d_uc_eosinophils_from_granulocytes.R
##
## UC biopsies, after 04: eosinophils that 04 left inside the neutrophil cluster
## are called by their own markers and relabelled "Eosinophil".
##
## Why a marker rule and not a re-clustering: on this panel only SIGLEC8 and RNASE3
## are eosinophil-specific (ALOX15 also in epithelium), each with 0-2 transcripts per
## cell. A small population defined by 2-3 genes rarely forms its own cluster, and
## requiring BOTH genes misses most eosinophils (dropouts). So, only inside the
## granulocyte cluster (no mast cells -- they also express SIGLEC8 -- no epithelium):
##   strict    SIGLEC8 >= 1 AND RNASE3 >= 1
##   default   SIGLEC8 + RNASE3 >= 2 transcripts, OR one of them + ALOX15        <- used
##   lenient   SIGLEC8 or RNASE3 >= 1
## and cells with mast-cell genes (CPA3 / MS4A2 / HDC / KIT >= 2) are left out.
## Choose with UC_EOS_RULE in config.R ("strict", "default", "lenient").
##
## Checks (figure 04d_eosinophil_call.png):
##   - false-positive rate: the same rule on T and B cells (should be ~0); the
##     eosinophil call is credible when it is far above that
##   - marker dot plot: called eosinophils vs remaining neutrophils vs mast cells
##   - sub-clustering of the granulocyte cluster on granulocyte genes only: does
##     one sub-cluster carry the called eosinophils?
##   - per core: eosinophils called vs the pathology report (05_uc_genes_eosinophils.R)
##
## Re-runnable: the 04 labels are kept in cell_type_04 and every run starts from them.
## Re-run 04d after every run of 04 (04 rewrites cell_type).
##
## INPUT : xenium_<seg>_annotated.rds (04)
## OUTPUT: cell_type (eosinophils relabelled), cell_type_04, eos_call_<rule> columns in
##         xenium_<seg>_annotated.rds; tables/04d_*.csv; figures/04d_eosinophil_call.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_EOS_PARENT  <- get0("UC_EOS_PARENT", ifnotfound = "neutro|granulo")   # labels searched (regex, lower case)
UC_EOS_RULE    <- get0("UC_EOS_RULE", ifnotfound = "default")
UC_EOS_LABEL   <- get0("UC_EOS_LABEL", ifnotfound = "Eosinophil")
EOS_SPECIFIC   <- c("SIGLEC8", "RNASE3", "RNASE2", "CLC", "CCR3", "IL5RA", "PRG2", "EPX")   # those on the panel are used
EOS_SUPPORT    <- "ALOX15"
NEUT_GENES     <- c("CSF3R", "S100A9", "S100A12", "TREM1", "CXCL8", "IL1B", "FCGR3B", "CXCR2")
MAST_GENES     <- c("CPA3", "MS4A2", "HDC", "KIT", "TPSAB1")

xen <- readRDS(obj_path("annotated"))
if (!"cell_type_04" %in% colnames(xen[[]])) xen$cell_type_04 <- xen$cell_type
lab04 <- as.character(xen$cell_type_04)
panel <- rownames(xen[["Xenium"]])
on_panel <- function(g) { g <- match_case(g, panel); g[g %in% panel] }
eg <- on_panel(EOS_SPECIFIC); sg <- on_panel(EOS_SUPPORT); ng <- on_panel(NEUT_GENES); mg <- on_panel(MAST_GENES)
message("Eosinophil-specific genes on the panel: ", paste(eg, collapse = ", "),
        " | support: ", paste(sg, collapse = ", "), " | neutrophil: ", paste(ng, collapse = ", "),
        " | mast: ", paste(mg, collapse = ", "))
if (!length(eg)) stop("No eosinophil-specific gene (", paste(EOS_SPECIFIC, collapse = ", "), ") on the panel.")

parent <- grepl(UC_EOS_PARENT, tolower(lab04)) & xen$lineage %in% "Immune"
if (!any(parent)) stop("No immune cells with a label matching '", UC_EOS_PARENT, "' -- set UC_EOS_PARENT in config.R ",
                       "to the name of the cluster that holds the eosinophils (labels: ",
                       paste(sort(unique(lab04[xen$lineage %in% "Immune"])), collapse = ", "), ").")
message("Searching in: ", paste(sort(unique(lab04[parent])), collapse = ", "), " (", format(sum(parent), big.mark = ","), " cells)")

cnt <- GetAssayData(xen, assay = "Xenium", layer = "counts")
gsum <- function(g) if (length(g)) Matrix::colSums(cnt[g, , drop = FALSE]) else rep(0, ncol(cnt))
gpos <- function(g) if (length(g)) Matrix::colSums(cnt[g, , drop = FALSE] > 0) else rep(0, ncol(cnt))
eos_n <- gsum(eg); eos_genes_pos <- gpos(eg); sup_n <- gsum(sg); mast_n <- gsum(mg)
rules <- list(
  strict  = eos_genes_pos >= min(2, length(eg)),
  default = eos_n >= 2 | (eos_n >= 1 & sup_n >= 1),
  lenient = eos_n >= 1)
if (!UC_EOS_RULE %in% names(rules)) stop("UC_EOS_RULE must be one of ", paste(names(rules), collapse = ", "))
not_mast <- mast_n < 2

## -----------------------------------------------------------------------
## 1. calls per rule, and the same rule on T / B cells (false-positive rate)
## -----------------------------------------------------------------------
ref <- grepl("t_cell|cd4|cd8|treg|tfh|b_cell|naive", tolower(lab04)) & xen$lineage %in% "Immune"
rt <- dplyr::bind_rows(lapply(names(rules), function(r) data.frame(
  rule = r, n_called = sum(rules[[r]] & parent & not_mast),
  pct_of_granulocyte_cluster = 100 * mean((rules[[r]] & not_mast)[parent]),
  pct_of_T_B_cells_background = if (any(ref)) 100 * mean(rules[[r]][ref]) else NA)))
rt$fold_over_background <- rt$pct_of_granulocyte_cluster / pmax(rt$pct_of_T_B_cells_background, 0.01)
save_table(rt, "04d_eosinophil_rules")
print(rt, row.names = FALSE)
for (r in names(rules)) xen[[paste0("eos_call_", r)]] <- rules[[r]] & parent & not_mast
call <- rules[[UC_EOS_RULE]] & parent & not_mast
if (rt$fold_over_background[rt$rule == UC_EOS_RULE] < 3)
  message("!! The '", UC_EOS_RULE, "' rule fires in T/B cells almost as often as in the granulocyte cluster -- ",
          "the call is close to background (spill-over / noise); try UC_EOS_RULE <- \"strict\".")

## -----------------------------------------------------------------------
## 2. sub-clusters of the granulocyte cluster on granulocyte genes only
## -----------------------------------------------------------------------
feats <- unique(c(eg, sg, ng, mg))
sub <- CreateSeuratObject(counts = cnt[feats, parent, drop = FALSE], assay = "Xenium")
sub <- NormalizeData(sub, scale.factor = stats::median(xen$nCount_Xenium), verbose = FALSE)
sub <- ScaleData(sub, features = feats, verbose = FALSE)
npc <- max(2, min(8, length(feats) - 1))
sub <- RunPCA(sub, features = feats, npcs = npc, verbose = FALSE, approx = FALSE)
sub <- FindNeighbors(sub, dims = seq_len(npc), verbose = FALSE)
sub <- FindClusters(sub, resolution = get0("UC_EOS_SUBCLUSTER_RES", ifnotfound = 0.1), verbose = FALSE)
sc <- data.frame(subcluster = paste0("g", sub$seurat_clusters), called = call[parent])
st <- sc |> group_by(subcluster) |> summarise(n_cells = n(), pct_called_eosinophil = 100 * mean(called), .groups = "drop") |>
  arrange(desc(pct_called_eosinophil))
save_table(st, "04d_granulocyte_subclusters")
message("Granulocyte sub-clusters (on ", length(feats), " granulocyte / mast genes) and % of their cells called eosinophil:")
print(as.data.frame(st), row.names = FALSE)

## -----------------------------------------------------------------------
## 3. relabel and save
## -----------------------------------------------------------------------
new <- lab04
new[call] <- UC_EOS_LABEL
xen$cell_type <- new
tab <- data.frame(label_04 = lab04, label_04d = new, slide_id = xen$slide_id) |>
  count(slide_id, label_04, label_04d, name = "n_cells") |> filter(label_04 != label_04d | grepl(UC_EOS_PARENT, tolower(label_04)))
save_table(tab, "04d_relabelled")
message(sprintf("Relabelled %s cells as '%s' (rule '%s'): %.1f%% of the granulocyte cluster.",
                format(sum(call), big.mark = ","), UC_EOS_LABEL, UC_EOS_RULE, 100 * mean(call[parent])))
if (any(!is.na(xen$core_id))) {
  pc <- data.frame(core_id = xen$core_id, call = call, parent = parent) |> filter(!is.na(core_id)) |>
    group_by(core_id) |> summarise(n_granulocyte_cluster = sum(parent), n_eosinophils = sum(call), .groups = "drop") |>
    arrange(desc(n_eosinophils))
  save_table(pc, "04d_eosinophils_per_core")
}

## -----------------------------------------------------------------------
## 4. figure
## -----------------------------------------------------------------------
grp <- dplyr::case_when(call ~ "called eosinophils", parent ~ "rest of granulocyte cluster",
                        grepl("mast", tolower(lab04)) ~ "mast cells", ref ~ "T / B cells (background)", TRUE ~ NA_character_)
keep <- which(!is.na(grp))
if (length(keep) > 2e5) keep <- sort(sample(keep, 2e5))
genes_show <- unique(c(eg, sg, ng, mg))
dd <- dplyr::bind_rows(lapply(genes_show, function(g) {
  x <- cnt[g, keep]
  data.frame(group = grp[keep], x = x) |> group_by(group) |>
    summarise(pct = 100 * mean(x > 0), mean_counts = mean(x), n = n(), .groups = "drop") |> mutate(gene = g)
}))
dd$gene <- factor(dd$gene, levels = genes_show)
dd$set <- factor(ifelse(dd$gene %in% c(eg, sg), "eosinophil", ifelse(dd$gene %in% ng, "neutrophil", "mast")),
                 levels = c("eosinophil", "neutrophil", "mast"))
dd$group <- factor(dd$group, levels = rev(c("called eosinophils", "rest of granulocyte cluster", "mast cells",
                                            "T / B cells (background)")))
p1 <- ggplot(dd, aes(gene, group, size = pct, colour = mean_counts)) + geom_point() +
  scale_size_area(max_size = 10, name = "% cells\npositive") +
  scale_colour_gradient(low = "#D9D9D9", high = "#B2182B", name = "mean\ncounts") +
  facet_grid(~ set, scales = "free_x", space = "free_x") + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = NULL, title = sprintf("Called eosinophils (rule '%s') vs the rest", UC_EOS_RULE),
       subtitle = "called cells: high SIGLEC8 / RNASE3, low mast genes\nT / B cells = background")
rt$rule <- factor(rt$rule, levels = names(rules))
rl <- tidyr::pivot_longer(rt, c(pct_of_granulocyte_cluster, pct_of_T_B_cells_background), names_to = "where", values_to = "pct")
rl$where <- ifelse(rl$where == "pct_of_granulocyte_cluster", "granulocyte cluster", "T / B cells (background)")
p2 <- ggplot(rl, aes(rule, pct, fill = where)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(0.8), vjust = -0.3, size = 3.3) +
  scale_fill_manual(values = c("granulocyte cluster" = "#C0392B", "T / B cells (background)" = "grey60"), name = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top") +
  labs(x = NULL, y = "% of cells meeting the rule", title = "Rules: signal vs background")
sd2 <- dplyr::bind_rows(lapply(feats, function(g) {
  x <- GetAssayData(sub, layer = "data")[g, ]
  data.frame(subcluster = paste0("g", sub$seurat_clusters), x = x) |> group_by(subcluster) |>
    summarise(pct = 100 * mean(x > 0), mean = mean(x), .groups = "drop") |> mutate(gene = g)
}))
sd2$gene <- factor(sd2$gene, levels = genes_show)
sd2$subcluster <- factor(sd2$subcluster, levels = rev(st$subcluster))
p3 <- ggplot(sd2, aes(gene, subcluster, size = pct, colour = mean)) + geom_point() +
  scale_size_area(max_size = 8, name = "% cells\npositive") +
  scale_colour_gradient(low = "#D9D9D9", high = "#3B6FB6", name = "mean\nnormalised") +
  scale_y_discrete(labels = function(x) sprintf("%s (%.0f%% called)", x, st$pct_called_eosinophil[match(x, st$subcluster)])) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Granulocyte cluster re-clustered on granulocyte genes only",
       subtitle = "a sub-cluster with SIGLEC8 / RNASE3 and most of the called cells = the clustering agrees with the rule")
save_fig((p1 | p2) / p3 + plot_layout(heights = c(1, 1.1), widths = c(2, 1)), "04d_eosinophil_call",
         width = 17, height = 11)

saveRDS(xen, obj_path("annotated"))
message("Saved -> ", obj_path("annotated"), "  (cell_type updated; 04 labels in cell_type_04)")
