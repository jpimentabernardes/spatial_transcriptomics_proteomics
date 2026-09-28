## =============================================================================
## 001_spatial_metabolomics_prep.R
##
## Combines the two-stage metabolomics preparation pipeline into one file:
##   PART 1: full, un-compromised
##     raw-data preprocessing + annotation. RUN ON YOUR CLUSTER, NOT A LAPTOP
##     OR BINDER -- 46,447 features x 193,200 pixels is genuinely large.
##   PART 2: builds the class-sized object
##     students actually use, from Part 1's output.
##
## Both parts still run on the cluster in one sitting, but Part 1's full
## annotated object is still saved to disk between parts -- so if you want to
## re-tune Part 2's subsetting parameters (OFF_TISSUE_CLUSTERS,
## KNOWN_FEATURES_ONLY, N_FEATURES_KEEP, ...) later, you can skip straight to
## Part 2 by loading that saved file, without repeating Part 1's expensive
## raw-data read and annotation.
##
## OUTPUT: data/processed/spamtp_metabolomics_FULL.rds (Part 1, full object)
##         data/processed/spamtp_metabolomics.rds      (Part 2, class-ready)
## =============================================================================

suppressPackageStartupMessages({
  library(Cardinal)
  library(SpaMTP)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)
setwd('/your/own/path/')

FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot_obj, name, width = 9, height = 7, dpi = 150) {
  path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
  message("Saved figure -> ", path)
  invisible(path)
}

## Same ggplot2-based pixel-map helper used elsewhere in this course
## (Cardinal's own formula-based image() dispatch has proven fragile across
## Cardinal versions/builds).
plot_pixel_map <- function(obj, column, title) {
  pd <- as.data.frame(Cardinal::pixelData(obj))
  stopifnot(column %in% colnames(pd))
  ggplot(pd, aes(x = x, y = y, fill = .data[[column]])) +
    geom_raster() +
    viridis::scale_fill_viridis(option = "magma") +
    coord_fixed() +
    scale_y_reverse() +
    theme_minimal() +
    ggtitle(title)
}

get_feature_meta <- function(obj) {
  tryCatch(
    obj[["Spatial"]][[]],          # Seurat v5 Assay5 feature-metadata accessor
    error = function(e) obj[["Spatial"]]@meta.data  # Seurat v4-style Assay fallback
  )
}
get_cell_meta <- function(obj) {
  tryCatch(obj[[]], error = function(e) obj@meta.data)
}

OUT_DIR <- "data/processed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
#full_path  <- file.path(OUT_DIR, "spamtp_metabolomics_FULL.rds")
class_path <- file.path(OUT_DIR, "spamtp_metabolomics.rds")


## =============================================================================
## PART 1: full raw-data preprocessing + annotation
## (was 00b_cluster_preprocess_metabolomics.R)
## =============================================================================

DATA_DIR   <- "data/metabolomics"
RUN_NAME   <- "run"
MASS_RANGE <- c(200, 1000)
POLARITY   <- "negative"

## A lot of RAM is needed, multi-core parallelism -- ADAPT to however many cores your job gets.
N_CORES <- min(16, parallel::detectCores())
setCardinalBPPARAM(BiocParallel::MulticoreParam(workers = N_CORES))

imzml_path <- file.path(DATA_DIR, paste0(RUN_NAME, ".imzML"))
message("Reading raw MSI data ...")
msi_raw <- Cardinal::readMSIData(imzml_path, mass.range = MASS_RANGE)
print(msi_raw)

msi_raw <- Cardinal::summarizePixels(msi_raw, c(TIC = "sum"))
save_fig(plot_pixel_map(msi_raw, "TIC", "Raw total ion image (pre-QC)"), "01_tic_raw_preQC")

## Filter peaks present in a real fraction of pixels
message("Normalising + frequency-filtering (full scan, affordable here) ...")
msi_peaks <- Cardinal::normalize(msi_raw, method = "tic")
msi_peaks <- Cardinal::peakFilter(msi_peaks, freq.min = 0.01)  # keep peaks in >=1% of pixels
msi_peaks <- Cardinal::process(msi_peaks)
print(msi_peaks)

sm <- SpaMTP::CardinalToSeurat(msi_peaks)
print(sm)

p_qc <- Seurat::ImageFeaturePlot(sm, features = "nFeature_Spatial", size = 1.5, dark.background = FALSE) +
  ggtitle("m/z features detected per pixel (post-filtering)")
save_fig(p_qc, "02_nfeature_per_pixel")

sm_annotated <- SpaMTP::AnnotateSM(
  sm,
  db        = rbind(SpaMTP::HMDB_db, SpaMTP::Lipidmaps_db, SpaMTP::Chebi_db),
  ppm_error = 3,     
  polarity  = POLARITY,
  adducts   = if (POLARITY == "positive") c("M+H", "M+K", "M+Na") else c("M-H")
)
sm_annotated[["Spatial"]]@meta.data <- SpaMTP::RefineLipids(
  sm_annotated[["Spatial"]]@meta.data,
  annotation.column = "all_IsomerNames",
  lipid_info = "simple"
)

## Annotation coverage -- how many of the retained m/z features actually got
## a database hit. 
meta <- tryCatch(
  sm_annotated[["Spatial"]][[]],
  error = function(e) sm_annotated[["Spatial"]]@meta.data
)
is_annotated <- !is.na(meta$all_IsomerNames) & meta$all_IsomerNames != ""
coverage_df <- data.frame(
  status = c("Annotated", "Unannotated"),
  n = c(sum(is_annotated), sum(!is_annotated))
)
p_coverage <- ggplot(coverage_df, aes(x = status, y = n, fill = status)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Annotation coverage", x = NULL, y = "Number of features")
save_fig(p_coverage, "03_annotation_coverage")
sm_annotated <- SpaMTP::NormalizeSMData(sm_annotated)

saveRDS(sm_annotated, full_path)


## =============================================================================
## PART 2: build the class-sized object
##
## FEATURE SELECTION STRATEGY (informed by actually inspecting the real
## annotation output, not guessed upfront):
##   - "Annotated" as a yes/no flag turned out to be useless here -- merging
##     HMDB+LipidMaps+ChEBI at 3 ppm finds SOME candidate for ~100% of
##     features by chance (ChEBI especially is a general chemical ontology,
##     not curated endogenous metabolites). What actually varies is
##     SPECIFICITY: number of candidate isomer names per feature ranges from
##     1 (unambiguous) to 1197 (meaningless). Fewer candidates = more
##     trustworthy call for teaching purposes.
##   - There are genuinely two panels here, visible directly in the data:
##     LIPID hits (Lipid.Maps.Category/Main.Class populated by RefineLipids())
##     vs SMALL-MOLECULE hits (everything else, from HMDB/ChEBI). The class
##     subset should represent both, not let one dominate.
##   - On top of panel + specificity, we prioritise features that are
##     actually useful for the exercises: spatially variable (so there's a
##     visible pattern to look at) and cluster markers (so the DE/heatmap
##     section has something real to show).
##
## Deliberately writes to the SAME filename script 02 (transcriptomics) does
## not depend on, and script 03 (integration) reads directly -- this lets
## students/Binder skip Part 1/2's raw-data handling entirely and go straight
## to integration.
## =============================================================================

sm_full <- sm_annotated

## -----------------------------------------------------------------------
## CONFIG -- tune these until the printed file size at the bottom is
## representative for the capabilities of you computer.
## if running the HPC no need to subset
## -----------------------------------------------------------------------
N_FEATURES_KEEP <- 1000    # total m/z features to keep in the class object
LIPID_FRACTION  <- 0.5     # target share of N_FEATURES_KEEP from the lipid panel 
MAX_CANDIDATES  <- 20      # soft filter: exclude features with more than this many candidate isomer
KNOWN_FEATURES_ONLY   <- TRUE   # confirmed: restrict to confidently-annotated features
MAX_CANDIDATES_KNOWN  <- 3      # threshold used only when KNOWN_FEATURES_ONLY = TRUE

OFF_TISSUE_CLUSTERS <- c("12", "17", "4", "5", "14", "10")  # confirmed via visual + TIC-based check

PIXEL_BIN_SIZE <- 1  # 1 = no binning; 2 = merge 2x2 pixel blocks (see final section)

## -----------------------------------------------------------------------
## 1. Optional hard pre-filter: known features only
##    Applied before clustering so a noisy tail of ambiguous features
##    doesn't influence the PCA/clustering used for everything downstream,
##    and so DE testing runs on fewer features too.
## -----------------------------------------------------------------------
if (KNOWN_FEATURES_ONLY) {
  meta0 <- get_feature_meta(sm_full)
  n_candidates0 <- lengths(strsplit(meta0$all_IsomerNames, ";"))
  known_features <- rownames(meta0)[n_candidates0 <= MAX_CANDIDATES_KNOWN]
  message("KNOWN_FEATURES_ONLY: keeping ", length(known_features), " of ", nrow(meta0),
          " features with <=", MAX_CANDIDATES_KNOWN, " candidate isomer names.")
  sm_full <- subset(sm_full, features = known_features)
}

## -----------------------------------------------------------------------
## 2. Cluster the full (or known-features-only) object -- Part 1 saves an
##    annotated but unclustered object; clustering happens here because this
##    part needs cluster identity for both off-tissue removal and
##    marker-based selection below.
## -----------------------------------------------------------------------
sm_full <- SpaMTP::NormalizeSMData(sm_full)
sm_full <- FindVariableFeatures(sm_full, verbose = FALSE)
sm_full <- ScaleData(sm_full, verbose = FALSE)
sm_full <- RunPCA(sm_full, npcs = 30, verbose = FALSE)
sm_full <- FindNeighbors(sm_full, dims = 1:30, verbose = FALSE)
sm_full <- FindClusters(sm_full, resolution = 0.5, cluster.name = "metabolite_clusters", verbose = FALSE)
sm_full$metabolite_clusters <- factor(sm_full$metabolite_clusters)

p_clusters_spatial <- Seurat::ImageDimPlot(
  sm_full, group.by = "metabolite_clusters", size = 1.5, dark.background = FALSE
) + ggtitle("Metabolite clusters (inspect this to pick OFF_TISSUE_CLUSTERS)")
save_fig(p_clusters_spatial, "00c_metabolite_clusters_spatial")

## Best-effort quantitative hint: off-tissue/matrix background tends to have
## HIGH mean TIC (matrix crystals often out-ionise real tissue signal). TIC
## as a per-pixel metadata column doesn't reliably survive CardinalToSeurat()
## -- if it's missing, compute it directly from the Spatial assay's counts
## instead of just skipping this.
cell_meta <- get_cell_meta(sm_full)
if (!"TIC" %in% colnames(cell_meta)) {
  counts <- Seurat::GetAssayData(sm_full, assay = "Spatial", slot = "counts")
  cell_meta$TIC <- Matrix::colSums(counts)
}
tic_by_cluster <- cell_meta |>
  dplyr::group_by(metabolite_clusters) |>
  dplyr::summarise(mean_TIC = mean(TIC, na.rm = TRUE), n_pixels = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_TIC))
print(tic_by_cluster, n = Inf)

## -----------------------------------------------------------------------
## 3. Optional: drop off-tissue background pixels, identified from the plot
##    (and/or TIC hint) above. Removing these HERE -- before DE/spatial-
##    variability testing -- means background noise can't influence which
##    features get flagged as markers or spatially variable, and makes those
##    steps themselves faster on a re-run, not just the final file smaller.
## -----------------------------------------------------------------------
if (length(OFF_TISSUE_CLUSTERS) > 0) {
  n_before <- ncol(sm_full)
  sm_full <- subset(sm_full, subset = metabolite_clusters %in% OFF_TISSUE_CLUSTERS, invert = TRUE)
  message("Removed off-tissue clusters ", paste(OFF_TISSUE_CLUSTERS, collapse = ", "), ": ",
          n_before, " -> ", ncol(sm_full), " pixels.")
} else {
  message("OFF_TISSUE_CLUSTERS is empty -- keeping all pixels this run. If ",
          "00c_metabolite_clusters_spatial.png shows an obvious background region, ",
          "set OFF_TISSUE_CLUSTERS and re-run for a smaller, cleaner object.")
}

## -----------------------------------------------------------------------
## 4. Feature selection
## -----------------------------------------------------------------------
meta <- get_feature_meta(sm_full)
meta$feature <- rownames(meta)
meta$n_candidates <- lengths(strsplit(meta$all_IsomerNames, ";"))
meta$panel <- ifelse(
  !is.na(meta$Lipid.Maps.Category) & meta$Lipid.Maps.Category != "",
  "lipid", "small_molecule"
)
message("Panel breakdown (all retained features): ",
        paste(names(table(meta$panel)), table(meta$panel), sep = "=", collapse = ", "))

## Diagnostic: is the small-molecule panel genuinely thin, or is something
## wrong with the panel split? DAN matrix in negative mode is a classic
## lipidomics-oriented MALDI method, so a real, heavy skew toward lipid hits
## here would not be surprising -- these plots let you confirm which it is
## rather than assuming either way.
lipid_meta <- meta[meta$panel == "lipid" & !is.na(meta$Lipid.Maps.Category) & meta$Lipid.Maps.Category != "", ]

p_lipid_category <- ggplot(
  as.data.frame(sort(table(lipid_meta$Lipid.Maps.Category), decreasing = TRUE)),
  aes(x = reorder(Var1, Freq), y = Freq)
) +
  geom_col(fill = "#4C72B0") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lipid panel: features per Lipid Maps Category (retained features)",
       x = NULL, y = "Number of features")
save_fig(p_lipid_category, "00c_lipid_category_breakdown")

## Main.Class is much finer-grained (many possible values) -- cap to the top
## 25 by count so the plot stays readable.
class_counts <- sort(table(lipid_meta$Lipid.Maps.Main.Class), decreasing = TRUE)
class_counts <- head(class_counts, 25)
p_lipid_class <- ggplot(
  as.data.frame(class_counts),
  aes(x = reorder(Var1, Freq), y = Freq)
) +
  geom_col(fill = "#55A868") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lipid panel: features per Main Class (top 25, retained features)",
       x = NULL, y = "Number of features")
save_fig(p_lipid_class, "00c_lipid_class_breakdown", height = 9)

## >>> EXERCISE 1 <<<
## Open 00c_lipid_category_breakdown.png. Which 1-2 Lipid Maps Categories
## dominate this dataset? DAN matrix in negative mode is a classic
## lipidomics-oriented MALDI method -- does the dominant category match what
## you'd expect it to favor (e.g. glycerophospholipids, sphingolipids)? If
## you have time, look up the top Main.Class within your dominant Category
## in 00c_lipid_class_breakdown.png and see if you recognise the lipid
## subclass from your own background.

## 4a. Cluster markers -- must-keep, regardless of annotation specificity
##     (a strong statistical DE signal is meaningful on its own).
##
##     IMPORTANT: FindAllDEMs()'s `n` argument controls pseudo-replicate
##     pooling for the underlying limma test ("Pooling one sample into N
##     replicates") -- it is NOT a top-N-markers cutoff. With
##     return.individual = FALSE, it returns a list whose $DEMs element is
##     the FULL feature x cluster DE table (every tested feature against
##     every cluster), not a pre-filtered marker list. We do the significance
##     filtering and top-N-per-cluster selection ourselves below.
TOP_MARKERS_PER_CLUSTER <- 15
DEM_FDR_MAX <- 0.05

message("Finding cluster markers (FindAllDEMs) ...")
dem_results <- tryCatch(
  SpaMTP::FindAllDEMs(
    data = sm_full, slot = "data", ident = "metabolite_clusters",
    n = 15, logFC_threshold = 1, return.individual = FALSE,
    run_name = "class_subset_markers", annotation.column = "all_IsomerNames"
  ),
  error = function(e) {
    message("FindAllDEMs() failed (", conditionMessage(e), ") -- continuing without marker features.")
    NULL
  }
)

marker_features <- tryCatch({
  dem_table <- if (is.list(dem_results) && "DEMs" %in% names(dem_results)) {
    dem_results$DEMs
  } else if (is.data.frame(dem_results)) {
    dem_results  # in case a future SpaMTP version returns the table directly
  } else {
    stop("unrecognised FindAllDEMs() return structure")
  }
  gene_col <- intersect(c("gene", "feature"), colnames(dem_table))[1]
  stopifnot(!is.na(gene_col), "cluster" %in% colnames(dem_table), "FDR" %in% colnames(dem_table))

  message("Full DE table: ", nrow(dem_table), " feature-cluster tests.")
  sig <- dem_table[!is.na(dem_table$FDR) & dem_table$FDR < DEM_FDR_MAX & abs(dem_table$logFC) >= 1, ]
  message("Significant (FDR<", DEM_FDR_MAX, ", |logFC|>=1): ", nrow(sig), " feature-cluster hits.")

  top_sig <- sig |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = abs(logFC), n = TOP_MARKERS_PER_CLUSTER) |>
    dplyr::ungroup()
  unique(top_sig[[gene_col]])
}, error = function(e) {
  message("Marker extraction from FindAllDEMs() output failed (", conditionMessage(e),
          ") -- continuing without marker features.")
  character(0)
})
message("Cluster-marker features (must-keep): ", length(marker_features))

## 4b. Spatially variable features -- ADAPT/VERIFY: FindSpatiallyVariableMetabolites()'s
##     exact argument names weren't confirmed against a live SpaMTP session while writing
##     this script -- if it errors below, check ?SpaMTP::FindSpatiallyVariableMetabolites
##     and adjust, or rely on the Seurat-variance fallback (still reasonable).
message("Finding spatially variable features ...")
spatial_var_features <- tryCatch({
  svg <- SpaMTP::FindSpatiallyVariableMetabolites(sm_full, assay = "Spatial", selection.method = "moransi")
  svg_meta <- get_feature_meta(svg)
  rank_col <- intersect(c("MoransI_observed", "moransi.spatial.cor", "SVM.rank"), colnames(svg_meta))
  if (length(rank_col) == 0) stop("couldn't find a Moran's I / ranking column in the returned object")
  rownames(svg_meta)[order(svg_meta[[rank_col[1]]], decreasing = TRUE)]
}, error = function(e) {
  message("FindSpatiallyVariableMetabolites() didn't run as expected (", conditionMessage(e),
          ") -- falling back to Seurat's FindVariableFeatures() ranking instead.")
  VariableFeatures(sm_full)
})

## 4c. Build the final feature list: must-keep markers first, then fill the
##     remaining budget from each panel (ranked by spatial variability,
##     restricted to reasonably specific annotations) at roughly
##     LIPID_FRACTION : (1 - LIPID_FRACTION).
fill_budget <- max(0, N_FEATURES_KEEP - length(marker_features))
n_lipid_target <- round(fill_budget * LIPID_FRACTION)
n_small_target <- fill_budget - n_lipid_target

specific_enough <- meta$feature[meta$n_candidates <= MAX_CANDIDATES]
rank_within <- function(panel_name, target_n) {
  candidates <- meta$feature[meta$panel == panel_name & meta$feature %in% specific_enough]
  candidates <- candidates[!candidates %in% marker_features]  # don't double count
  ranked <- intersect(spatial_var_features, candidates)       # rank by spatial variability
  ranked <- unique(c(ranked, candidates))                     # append any not covered by the SVM ranking
  head(ranked, target_n)
}
lipid_fill <- rank_within("lipid", n_lipid_target)
small_fill <- rank_within("small_molecule", n_small_target)

keep_features <- unique(c(marker_features, lipid_fill, small_fill))
keep_features <- intersect(keep_features, rownames(sm_full))  # safety: only real feature names
if (length(keep_features) > N_FEATURES_KEEP) keep_features <- keep_features[seq_len(N_FEATURES_KEEP)]

final_panel_counts <- table(meta$panel[meta$feature %in% keep_features])
message("Final selection: ", length(keep_features), " features (",
        length(marker_features), " cluster markers, ",
        length(lipid_fill), " lipid fill, ", length(small_fill), " small-molecule fill). ",
        "Panel breakdown: ",
        paste(names(final_panel_counts), final_panel_counts, sep = "=", collapse = ", "))

save_fig(
  ggplot(as.data.frame(final_panel_counts), aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() + geom_text(aes(label = Freq), vjust = -0.3) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = "Class subset: panel composition", x = NULL, y = "Number of features"),
  "00c_panel_composition"
)

## Same Category breakdown as before, but restricted to the final class
## subset -- confirms selection kept real lipid-class diversity rather than
## collapsing onto just one or two dominant classes.
final_lipid_meta <- lipid_meta[lipid_meta$feature %in% keep_features, ]
if (nrow(final_lipid_meta) > 0) {
  p_final_lipid_category <- ggplot(
    as.data.frame(sort(table(final_lipid_meta$Lipid.Maps.Category), decreasing = TRUE)),
    aes(x = reorder(Var1, Freq), y = Freq)
  ) +
    geom_col(fill = "#C44E52") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Lipid panel: features per Category (final class subset)",
         x = NULL, y = "Number of features")
  save_fig(p_final_lipid_category, "00c_final_lipid_category_breakdown")
}

sm_class <- subset(sm_full, features = keep_features)

## >>> EXERCISE 2 <<<
## KNOWN_FEATURES_ONLY, MAX_CANDIDATES_KNOWN, and OFF_TISSUE_CLUSTERS above
## were all tuned by actually inspecting this dataset (annotation
## specificity distribution, TIC-per-cluster table, spatial cluster plot) --
## not chosen blindly. Try setting MAX_CANDIDATES_KNOWN <- 1 instead of 3
## (only features with a SINGLE unambiguous candidate name) and re-run from
## section 1 onward. How many features survive? Look at
## 00c_panel_composition.png before and after -- does the lipid/small-
## molecule balance hold up, or does one panel nearly disappear?

## -----------------------------------------------------------------------
## 5. Optional: spatially bin pixels to reduce resolution further (NOT
##    random subsampling -- that would break the spatial image structure
##    students need for plots. Binning merges NxN blocks of adjacent pixels
##    into one, summing intensities, which keeps the image coherent at lower
##    resolution -- similar in spirit to how microscopy images get binned.)
##    Complements, rather than replaces, the off-tissue removal above --
##    off-tissue removal cuts irrelevant pixels; binning reduces resolution
##    of the remaining (relevant) ones.
## -----------------------------------------------------------------------
if (PIXEL_BIN_SIZE > 1) {
  message("Spatially binning pixels ", PIXEL_BIN_SIZE, "x", PIXEL_BIN_SIZE, " ...")
  coords <- SpaMTP::GetTissueCoordinates(sm_class)
  bin_x <- floor((coords$x - 1) / PIXEL_BIN_SIZE)
  bin_y <- floor((coords$y - 1) / PIXEL_BIN_SIZE)
  bin_id <- paste0(bin_x, "_", bin_y)

  expr <- Seurat::GetAssayData(sm_class, assay = "Spatial", slot = "counts")
  binned <- t(rowsum(t(as.matrix(expr)), group = bin_id))  # sum within each bin

  bin_coords <- data.frame(bin_id = bin_id, x = bin_x, y = bin_y) |>
    dplyr::distinct(bin_id, .keep_all = TRUE)
  rownames(bin_coords) <- bin_coords$bin_id
  bin_coords <- bin_coords[colnames(binned), c("x", "y")]

  sm_class <- Seurat::CreateSeuratObject(counts = binned, assay = "Spatial")
  sm_class[["fov"]] <- Seurat::CreateFOV(
    coords = bin_coords, type = "centroids", assay = "Spatial"
  )
  message("Binned to ", ncol(sm_class), " pixels (from ",
          nrow(coords), " original).")
}

print(sm_class)
saveRDS(sm_class, class_path)
