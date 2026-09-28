## =============================================================================
## 003b_manual_alignment.R
##
## Required before 003_integration_spamtp.R. SpaMTP::AlignSpatialOmics() does
## not work on Xenium/FOV-based ST objects (Visium-only @scale.factors slot)
## -- confirmed structural incompatibility, not fixable by arguments. This
## does the same underlying job via manually-identified landmark points and a
## standard least-squares affine registration instead of AlignSpatialOmics()'s
## interactive Shiny app -- a legitimate, standard registration approach, and
## arguably more appropriate here anyway since the MALDI and Xenium slides
## are ADJACENT sections, not the same slide.
##
## THIS SCRIPT NEEDS TWO PASSES the first time:
##   Pass 1: run through STEP 1, inspect figures/003b_landmark_picking_grids.png,
##           identify 3-6 matching tissue landmarks visible in both panels.
##   Pass 2: fill in landmarks_msi/landmarks_xenium with real coordinates and
##           re-run the rest of the script.
##
## OUTPUT: data/processed/sm_transformed_manual_new.rds
## =============================================================================

suppressPackageStartupMessages({
  library(SpaMTP)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

PROCESSED_DIR <- "data/processed"
SM_PATH <- file.path(PROCESSED_DIR, "spamtp_metabolomics.rds")
ST_PATH <- file.path(PROCESSED_DIR, "seurat_xenium.rds")
OUT_PATH <- file.path(PROCESSED_DIR, "sm_transformed_manual_new.rds")

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

sm <- readRDS(SM_PATH)
st <- readRDS(ST_PATH)

df_sm <- Seurat::GetTissueCoordinates(sm)  # sm is a plain Seurat object -- no SpaMTP:: prefix needed
df_st <- Seurat::GetTissueCoordinates(st)

## -----------------------------------------------------------------------
## STEP 1: plot both tissue outlines WITH COORDINATE GRIDLINES, so you can
## read off approximate (x,y) values for matching landmarks by eye -- there's
## no click-to-pick interactivity here, just visual coordinate reading. Look
## for the same recognisable feature in both (a tissue corner, a fold, the
## central lumen visible in the metabolomics TIC image, an obvious
## anatomical boundary) and note its approximate coordinates in each plot.
## -----------------------------------------------------------------------
p_grid_st <- ggplot(df_st, aes(x, y)) +
  geom_point(size = 0.15, color = "lightblue4", alpha = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major = element_line(color = "grey70")) +
  ggtitle("Xenium coordinates (pick landmark x,y from here)")

p_grid_sm <- ggplot(df_sm, aes(x, y)) +
  geom_point(size = 0.15, color = "indianred", alpha = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major = element_line(color = "grey70")) +
  ggtitle("MALDI coordinates (pick landmark x,y from here)")

save_fig(p_grid_st | p_grid_sm, "003b_landmark_picking_grids", width = 14)

## -----------------------------------------------------------------------
## STEP 2: matched landmark coordinates, same order (row 1 in MSI = row 1 in
## Xenium = the same physical tissue feature). Need at least 3 non-collinear
## pairs; 4-6 well-spread points give a more robust fit than the bare minimum.
## -----------------------------------------------------------------------
landmarks_msi <- data.frame(
  x = c(450, 410, 90, 120, 340),
  y = c(240, 250, 40, 190, 125)
)
landmarks_xenium <- data.frame(
  x = c(6600, 6080, 1950, 1800, 5450),
  y = c(4700, 4600, 500, 2950, 2500)
)

## -----------------------------------------------------------------------
## STEP 3: fit a 2D affine transform (translation + rotation + scale + shear)
## mapping MSI coordinates -> Xenium coordinates via least squares.
##
## IMPORTANT: build an explicit, cleanly-named data.frame for lm() rather
## than writing landmarks_msi$x directly in the formula -- the latter creates
## an internal variable name that predict() cannot reliably match against
## newdata later, and fails SILENTLY (no error, just recycles training data
## instead of transforming real per-pixel coordinates).
## -----------------------------------------------------------------------
landmark_df <- data.frame(
  sx = landmarks_msi$x, sy = landmarks_msi$y,
  tx = landmarks_xenium$x, ty = landmarks_xenium$y
)
fit_x <- lm(tx ~ sx + sy, data = landmark_df)
fit_y <- lm(ty ~ sx + sy, data = landmark_df)

print(data.frame(
  landmark = seq_len(nrow(landmark_df)),
  resid_x = round(residuals(fit_x), 1),
  resid_y = round(residuals(fit_y), 1)
))

apply_affine <- function(xy, fit_x, fit_y) {
  newdata <- data.frame(sx = xy$x, sy = xy$y)
  data.frame(
    x = predict(fit_x, newdata = newdata),
    y = predict(fit_y, newdata = newdata)
  )
}

## -----------------------------------------------------------------------
## STEP 4: apply the fitted transform to ALL metabolomics pixel coordinates,
## write them back into the SpaMTP object, and save.
## -----------------------------------------------------------------------
df_sm_transformed <- apply_affine(df_sm, fit_x, fit_y)

sm_transformed <- sm
new_coords <- df_sm_transformed
rownames(new_coords) <- df_sm$cell

## Sanity check BEFORE trusting this output.
stopifnot(
  "Row count mismatch -- transform did not run on all pixels" =
    nrow(new_coords) == nrow(df_sm),
  "Transformed x has suspiciously little variation" =
    length(unique(round(new_coords$x, 1))) > 10,
  "Transformed y has suspiciously little variation" =
    length(unique(round(new_coords$y, 1))) > 10
)

centroids_obj <- sm[["fov"]]@boundaries$centroids
orig_coords <- centroids_obj@coords
cell_order <- centroids_obj@cells
stopifnot(all(cell_order %in% rownames(new_coords)))

new_coords_matrix <- as.matrix(new_coords[cell_order, c("x", "y")])
dimnames(new_coords_matrix) <- dimnames(orig_coords)

sm_transformed[["fov"]]@boundaries$centroids@coords <- new_coords_matrix

old_head <- head(Seurat::GetTissueCoordinates(sm))
new_head <- head(Seurat::GetTissueCoordinates(sm_transformed))
stopifnot(!isTRUE(all.equal(old_head$x, new_head$x)))
message("Coordinate write confirmed.")

## -----------------------------------------------------------------------
## STEP 5: visual check -- overlay both tissues in the SAME coordinate space
## now. They won't overlap perfectly (adjacent sections, not the same
## slide), but the overall tissue OUTLINE and orientation should now roughly
## correspond -- similar overall shape/size, same rotation, no mirroring.
## -----------------------------------------------------------------------
df_check <- rbind(
  data.frame(x = df_st$x, y = df_st$y, sample = "Xenium"),
  data.frame(x = df_sm_transformed$x, y = df_sm_transformed$y, sample = "MALDI (transformed)")
)
p_after <- ggplot(df_check, aes(x, y, color = sample)) +
  geom_point(size = 0.2, alpha = 0.4) +
  scale_color_manual(values = c("Xenium" = "lightblue4", "MALDI (transformed)" = "indianred")) +
  coord_fixed() + theme_bw() +
  ggtitle("Overlay after manual alignment -- shapes should roughly correspond")
print(p_after)
save_fig(p_after, "003b_post_alignment_overlay")

saveRDS(sm_transformed, OUT_PATH)
message("Saved manually-aligned object -> ", OUT_PATH)
