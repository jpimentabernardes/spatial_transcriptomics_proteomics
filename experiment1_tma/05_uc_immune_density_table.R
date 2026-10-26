## =============================================================================
## 05_uc_immune_density_table.R
##
## UC biopsies, after 04 / 04d: one row per core with the density (cells per mm^2
## of tissue) of each immune population, in the layout of the CellScape tables:
##
##   slide  tma_core  TMA_Sector_join  Koordinate_join  B_cells_density_per_mm2 ...
##   FL0035SectorA  A_1-c  A  1c  17.7 ...
##
## Populations come from the 04 / 04d cell type labels, mapped in
## config/uc_density_categories.csv (category, regular expression on the label --
## edit it to change what counts as e.g. M1 / M2 macrophages). On top of the labels:
##   - CD4 is not on the panel: T cells WITHOUT CD8A / CD8B count as CD4; a generic
##     "T_cell" label is split by CD8A / CD8B (>= 1 transcript = CD8)
##   - FOXP3: CD4 / CD8 T cells with >= UC_DENSITY_FOXP3_MIN FOXP3 transcripts go to
##     T_cells_CD4FOXP3 / T_cells_CD8FOXP3 (as a FOXP3 stain would); the Treg label is CD4FOXP3
##   - a category with an empty regex (Special_Macrophages by default) is written as NA
## Each cell is counted once (first matching category, in the order of the csv).
##
## Area = tissue area of the core from the cells (occupied 20 um bins), so holes and
## torn edges are not counted; CellScape may use the annotated core area instead --
## compare ratios between cores rather than absolute numbers if they differ.
##
## slide = the CellScape slide of the Xenium slide + "Sector" + sector (UC_DENSITY_SLIDE_IDS:
## TMA_UC_A = FL0033, TMA_UC_CD = FL0034 -> FL0033SectorA, FL0034SectorC, FL0034SectorD);
## tma_core = <sector>_<number>-<letter> (A_1-c), from the core_id.
##
## INPUT : xenium_<seg>_annotated.rds (04 / 04d), config/uc_density_categories.csv
## OUTPUT: data/export/05uc_immune_densities_per_core_<seg>.csv   (the table)
##         tables/05uc_immune_density_label_mapping.csv            (which label went where -- check it)
##         tables/05uc_immune_density_counts_per_core.csv          (cell counts and area behind the densities)
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/uc_helpers.R")

UC_DENSITY_CATEGORIES_CSV <- get0("UC_DENSITY_CATEGORIES_CSV", ifnotfound = file.path(CONFIG_DIR, "uc_density_categories.csv"))
UC_DENSITY_SLIDE_IDS      <- get0("UC_DENSITY_SLIDE_IDS", ifnotfound = c(TMA_UC_A = "FL0033", TMA_UC_CD = "FL0034"))
UC_DENSITY_FOXP3_MIN      <- get0("UC_DENSITY_FOXP3_MIN", ifnotfound = 1)
UC_DENSITY_SPLIT_CD8      <- get0("UC_DENSITY_SPLIT_CD8", ifnotfound = TRUE)    # generic T cells -> CD4 / CD8 by CD8A / CD8B
EXPORT_DIR <- file.path(PROJECT_DIR, "data", "export")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

cats <- read_config_csv(UC_DENSITY_CATEGORIES_CSV, stringsAsFactors = FALSE)
cats$label_regex[is.na(cats$label_regex)] <- ""
xen <- readRDS(obj_path("annotated"))
m <- uc_cell_table(xen)
counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")[, rownames(m)]
rm(xen); invisible(gc())
gene_n <- function(g) { g <- match_case(g, rownames(counts)); g <- g[g %in% rownames(counts)]
  if (length(g)) Matrix::colSums(counts[g, , drop = FALSE]) else rep(0, ncol(counts)) }
cd8 <- gene_n(c("CD8A", "CD8B")); foxp3 <- gene_n("FOXP3")

## -----------------------------------------------------------------------
## 1. category per cell (immune cells only; first matching row of the csv wins)
## -----------------------------------------------------------------------
lab <- tolower(m$label)
imm <- m$compartment == "Immune"
m$category <- NA_character_
for (k in seq_len(nrow(cats))) {
  rx <- cats$label_regex[k]
  if (!nzchar(rx)) next
  hit <- imm & is.na(m$category) & grepl(rx, lab, perl = TRUE)
  m$category[hit] <- cats$category[k]
}
is_cd4 <- m$category %in% "T_cells_CD4"; is_cd8 <- m$category %in% "T_cells_CD8"
if (UC_DENSITY_SPLIT_CD8) {
  generic <- is_cd4 & grepl("^t_cell|^t cell|^t$", lab)          # generic T cells: CD8A / CD8B decides
  m$category[generic & cd8 >= 1] <- "T_cells_CD8"
}
m$category[m$category %in% "T_cells_CD4" & foxp3 >= UC_DENSITY_FOXP3_MIN] <- "T_cells_CD4FOXP3"
m$category[m$category %in% "T_cells_CD8" & foxp3 >= UC_DENSITY_FOXP3_MIN] <- "T_cells_CD8FOXP3"

map <- m[imm, ] |> count(label, category = ifelse(is.na(category), "(not counted)", category), name = "n_cells") |>
  arrange(category, desc(n_cells))
save_table(map, "05uc_immune_density_label_mapping")
message("Immune labels -> categories (check tables/05uc_immune_density_label_mapping.csv; edit ",
        basename(UC_DENSITY_CATEGORIES_CSV), " to change):")
print(as.data.frame(map |> group_by(category) |> summarise(n_cells = sum(n_cells), labels = paste(unique(label), collapse = ", "))),
      row.names = FALSE, right = FALSE)
empty <- cats$category[!nzchar(cats$label_regex) & !cats$category %in% c("T_cells_CD8FOXP3")]
if (length(empty)) message("No labels defined for: ", paste(empty, collapse = ", "), " -> NA in the table.")

## -----------------------------------------------------------------------
## 2. densities per core
## -----------------------------------------------------------------------
cores <- sort(unique(m$core_id))
area <- vapply(cores, function(cid) { i <- m$core_id == cid; core_area_mm2(m$x_um[i], m$y_um[i]) }, numeric(1))
cnt <- as.data.frame.matrix(table(factor(m$core_id, levels = cores), factor(m$category, levels = cats$category)))
dens <- sweep(as.matrix(cnt), 1, area, "/")
dens[, cats$category[!nzchar(cats$label_regex) & !cats$category %in% "T_cells_CD8FOXP3"]] <- NA
key <- strsplit(core_key_from_id(cores), "|", fixed = TRUE)
sector <- toupper(vapply(key, `[`, "", 1)); num <- vapply(key, `[`, "", 2); let <- tolower(vapply(key, `[`, "", 3))
xslide <- m$slide_id[match(cores, m$core_id)]
fl <- ifelse(xslide %in% names(UC_DENSITY_SLIDE_IDS), UC_DENSITY_SLIDE_IDS[xslide], xslide)
if (any(!xslide %in% names(UC_DENSITY_SLIDE_IDS)))
  message("No CellScape slide name for ", paste(unique(xslide[!xslide %in% names(UC_DENSITY_SLIDE_IDS)]), collapse = ", "),
          " -- the Xenium slide_id is used (set UC_DENSITY_SLIDE_IDS in config.R).")
out <- data.frame(slide = paste0(fl, "Sector", sector), tma_core = sprintf("%s_%s-%s", sector, num, let),
                  TMA_Sector_join = sector, Koordinate_join = paste0(num, let), check.names = FALSE)
for (c in cats$category) out[[paste0(c, "_density_per_mm2")]] <- dens[, c]
o <- order(sector, as.integer(num), let)
out <- out[o, ]
f <- file.path(EXPORT_DIR, sprintf("05uc_immune_densities_per_core_%s.csv", SEGMENTATION))
utils::write.csv(out, f, row.names = FALSE, na = "NA")
message("Density table (", nrow(out), " cores) -> ", f)
qc <- data.frame(core_id = cores, xenium_slide = xslide, tma_core = sprintf("%s_%s-%s", sector, num, let), area_mm2 = round(area, 4),
                 n_cells = as.integer(table(factor(m$core_id, levels = cores))), cnt, check.names = FALSE)[o, ]
save_table(qc, "05uc_immune_density_counts_per_core")
print(utils::head(out, 3), row.names = FALSE)
message("Done.")
