## =============================================================================
## 05d_plot_gene_core.R
##
## One gene in one core, fluorescence style: the real DAPI image (optionally + a
## protein stain) on black, cells carrying the gene in yellow -> white (brighter =
## more transcripts). EXP1_05D_STAIN=dapi | dapi+boundary | dapi+protein | none.
##
## Run:   Rscript 05d_plot_gene_core.R NCAM1 "A 4F" TMA_TUMOR
##        Rscript 05d_plot_gene_core.R NCAM1 A_4-f            (slide optional)
## The core is matched ignoring case, spaces, "_" and "-" against core_id and, where
## 02 kept them, core_name / cellscape_core (A 4F = A_4f = A_4-f = a4f). If it is not
## found, the cores of the slide are listed.
##
## INPUT : xenium_<seg>_annotated.rds (or _cores.rds if there is no annotated object)
## OUTPUT: figures/05d_<gene>_<core>.png  (and _zoom.png: the ZOOM_UM window with most positive cells)
## Dot size / zoom: EXP1_05D_DOT=0.4 (smaller dots), EXP1_05D_ZOOM_UM=800 (wider zoom)
## =============================================================================

source("config/config.R")
source("R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 05d_plot_gene_core.R <GENE> "<core>" [slide_id]')
GENE <- args[1]; CORE <- args[2]; SLIDE <- if (length(args) >= 3) args[3] else NA
DAPI_BLUE <- "#2F6BFF"
DOT_SCALE <- as.numeric(Sys.getenv("EXP1_05D_DOT", "0.55"))     # dot width as a fraction of the cell (1 = full cell)
ZOOM_UM   <- as.numeric(Sys.getenv("EXP1_05D_ZOOM_UM", "600"))  # side of the zoom window

f <- obj_path("annotated"); if (!file.exists(f)) f <- obj_path("cores")
message("Object: ", f)
xen <- readRDS(f)
m <- xen[[]]
m <- m[!is.na(m$core_id), ]
if (!is.na(SLIDE)) m <- m[m$slide_id == SLIDE, ]
norm <- function(x) tolower(gsub("[^[:alnum:]]", "", as.character(x)))
## core names from 02 (core_name / cellscape_core), if the object does not carry them
alias <- data.frame(core_id = unique(m$core_id))
cores_file <- file.path(TABLES_DIR, "02_cores_all.csv")
if (file.exists(cores_file)) {
  ca <- utils::read.csv(cores_file, stringsAsFactors = FALSE)
  for (col in intersect(c("core_name", "cellscape_core"), colnames(ca))) alias[[col]] <- ca[[col]][match(alias$core_id, ca$core_id)]
}
hit <- apply(alias, 1, function(r) any(norm(r) == norm(CORE), na.rm = TRUE))
if (!any(hit)) {
  ## also allow the slide's sector letter in front of the core id (A + 4F)
  hit <- norm(alias$core_id) == norm(CORE) | endsWith(norm(alias$core_id), norm(CORE))
}
if (sum(hit) != 1) {
  show <- utils::head(alias[order(alias$core_id), , drop = FALSE], 60)
  stop(if (!any(hit)) "Core '" else "More than one core matches '", CORE, "'", if (!is.na(SLIDE)) paste0(" on ", SLIDE),
       ". Cores", if (!is.na(SLIDE)) paste0(" of ", SLIDE), " (first 60):\n",
       paste(utils::capture.output(print(show, row.names = FALSE)), collapse = "\n"), call. = FALSE)
}
cid <- alias$core_id[hit]
g <- match_case(GENE, rownames(xen[["Xenium"]]))
if (!g %in% rownames(xen[["Xenium"]])) stop(GENE, " is not on the panel.")
d <- m[m$core_id == cid, ]
d$expr <- as.numeric(GetAssayData(xen, assay = "Xenium", layer = "counts")[g, rownames(d)])
rm(xen); invisible(gc())
message(sprintf("%s in core %s (%s): %d of %d cells positive, %d transcripts", g, cid, d$slide_id[1],
                sum(d$expr > 0), nrow(d), sum(d$expr)))

## ---------------------------------------------------------------------------
## background: the Xenium morphology image (EXP1_05D_STAIN)
##   dapi            DAPI in blue (default)
##   dapi+boundary   + cell-boundary stain (ATP1A1 / CD45 / E-cadherin) in grey
##   dapi+protein    + interior protein stain (aSMA / vimentin) in grey
##   none            no image: every cell drawn as a blue dot (as before)
## The image is cut with py/crop_region.py (PYTHON_BIN, needs tifffile + zarr: the
## exp1-imaging environment). Xenium folder from config/slides.csv, or EXP1_05D_XENIUM_DIR.
## ---------------------------------------------------------------------------
STAIN <- tolower(Sys.getenv("EXP1_05D_STAIN", "dapi"))
xdir <- Sys.getenv("EXP1_05D_XENIUM_DIR", "")
if (!nzchar(xdir) && STAIN != "none") {
  sl <- read_slides()
  xdir <- if (d$slide_id[1] %in% sl$slide_id) resolve_path(sl$xenium_dir[sl$slide_id == d$slide_id[1]]) else ""
}
if (STAIN != "none" && (!nzchar(xdir) || !dir.exists(xdir))) {
  message("!! Xenium output folder of ", d$slide_id[1], " not found ('", xdir, "'): no stain, cells drawn as dots. ",
          "Set EXP1_05D_XENIUM_DIR=/path/to/output-XETG... (the folder with morphology_focus/).")
  STAIN <- "none"
}
chan_files <- if (STAIN != "none") sort(Sys.glob(file.path(xdir, "morphology_focus", "*.ome.tif"))) else character(0)
if (STAIN != "none" && !length(chan_files)) chan_files <- Sys.glob(file.path(xdir, "morphology_focus.ome.tif"))
stretch <- function(v) {                                   # 1st-99.8th percentile, mild gamma
  q <- stats::quantile(v[v > 0], c(0.01, 0.998), na.rm = TRUE)
  if (!all(is.finite(q)) || q[2] <= q[1]) return(v * 0)
  pmin(pmax((v - q[1]) / (q[2] - q[1]), 0), 1)^0.8
}
get_bg <- function(x0, y0, x1, y1, tag) {
  if (STAIN == "none") return(NULL)
  crop <- function(file, kind) {
    r <- crop_morphology(xdir, x0, y0, x1, y1, paste0(tag, "_", kind), image = file)
    list(img = r$image, meta = r$meta)
  }
  dapi_file <- chan_files[grepl("dapi", basename(chan_files), ignore.case = TRUE)][1]
  if (is.na(dapi_file)) dapi_file <- chan_files[1]
  dp <- crop(dapi_file, "dapi")
  extra <- NULL
  idx <- c("dapi+boundary" = 2, "dapi+protein" = 4)[STAIN]
  if (!is.na(idx)) {
    if (length(chan_files) >= idx) extra <- crop(chan_files[idx], sub("dapi\\+", "", STAIN))$img else
      message("!! No ", sub("dapi\\+", "", STAIN), " stain in ", xdir, " (only ", length(chan_files),
              " morphology channel(s): multimodal segmentation kit not used?) -- DAPI only.")
  }
  k <- max(1, ceiling(max(dim(dp$img)) / 2200))           # keep ~2000 px across
  sub_px <- function(a) a[seq(1, nrow(a), by = k), seq(1, ncol(a), by = k), drop = FALSE]
  dv <- stretch(sub_px(dp$img)); ev <- if (is.null(extra)) dv * 0 else stretch(sub_px(extra))
  rgb <- array(0, c(dim(dv), 3))
  rgb[, , 1] <- pmin(1, 0.10 * dv + 0.55 * ev)
  rgb[, , 2] <- pmin(1, 0.30 * dv + 0.55 * ev)
  rgb[, , 3] <- pmin(1, 1.00 * dv + 0.55 * ev)
  list(rgb = rgb, x0 = dp$meta$x0_um, x1 = dp$meta$x1_um, y0 = dp$meta$y0_um, y1 = dp$meta$y1_um)
}

## y is plotted as -y_um (image convention: y grows downwards)
plot_core <- function(d, panel_in, title, subtitle, bg = NULL) {
  ext <- max(diff(range(d$x_um)), diff(range(d$y_um)), 50)
  area <- if ("cell_area" %in% colnames(d)) d$cell_area else NA
  diam <- ifelse(is.na(area) | area <= 0, 10, 2 * sqrt(area / pi))
  d$size <- pmax(0.05, DOT_SCALE * pmin(diam, 25) * panel_in * 25.4 / ext / 0.75)
  pos <- d[d$expr > 0, ]; pos <- pos[order(pos$expr), ]
  bar <- if (ext > 600) 200 else if (ext > 200) 100 else 50
  p <- ggplot(mapping = aes(x_um, -y_um))
  if (!is.null(bg)) {
    p <- p + annotation_raster(bg$rgb, xmin = bg$x0, xmax = bg$x1, ymin = -bg$y1, ymax = -bg$y0) +
      geom_blank(data = data.frame(x_um = c(bg$x0, bg$x1), y_um = c(bg$y0, bg$y1)))
  } else {
    p <- p + geom_point(data = d, aes(size = size), colour = DAPI_BLUE, alpha = 0.75, stroke = 0)
  }
  p + geom_point(data = pos, aes(size = size * 1.4, fill = pmin(expr, 5)), shape = 21, colour = "black", stroke = 0.15) +
    scale_fill_gradientn(colours = c("#FFD000", "#FFF080", "#FFFFFF"), limits = c(1, 5), breaks = 1:5,
                         labels = c("1", "2", "3", "4", "5+"), name = paste(g, "\ntranscripts\nper cell")) +
    scale_size_identity() + coord_fixed(expand = FALSE) + theme_void(base_size = 12) +
    annotate("segment", x = min(d$x_um), xend = min(d$x_um) + bar, y = -max(d$y_um) + 0.05 * ext,
             yend = -max(d$y_um) + 0.05 * ext, colour = "white", linewidth = 0.9) +
    annotate("text", x = min(d$x_um) + bar / 2, y = -max(d$y_um) + 0.025 * ext, label = paste(bar, "um"),
             colour = "white", size = 3.5) +
    theme(plot.background = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          legend.background = element_rect(fill = "black", colour = NA),
          text = element_text(colour = "white"), legend.text = element_text(colour = "white"),
          plot.title = element_text(colour = "white", face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = element_text(colour = "grey80", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10)) +
    labs(title = title, subtitle = subtitle)
}
stain_txt <- c(none = "blue = all cells", dapi = "blue = DAPI", "dapi+boundary" = "blue = DAPI, grey = cell boundary stain",
               "dapi+protein" = "blue = DAPI, grey = interior protein stain")[STAIN]
sub <- sprintf("%s  |  yellow -> white = %s+ cells (%d of %d, %.1f%%)", stain_txt, g,
               sum(d$expr > 0), nrow(d), 100 * mean(d$expr > 0))
tag <- paste0("05d_", make.names(g), "_", make.names(cid))
pad <- 20
bg <- get_bg(min(d$x_um) - pad, min(d$y_um) - pad, max(d$x_um) + pad, max(d$y_um) + pad, tag)
save_fig(plot_core(d, 9, sprintf("%s  -  core %s  (%s)", g, cid, d$slide_id[1]), sub, bg), tag,
         width = 11, height = 10, dpi = 300)

## zoom: the ZOOM_UM window with the most positive cells (or the densest one)
w <- ZOOM_UM
if (diff(range(d$x_um)) > 1.5 * w) {
  b <- paste(floor(d$x_um / (w / 2)), floor(d$y_um / (w / 2)))
  score <- tapply(d$expr > 0, b, sum) + 1e-3 * tapply(rep(1, nrow(d)), b, sum)
  best <- strsplit(names(which.max(score)), " ")[[1]]
  cx <- (as.numeric(best[1]) + 0.5) * w / 2; cy <- (as.numeric(best[2]) + 0.5) * w / 2
  z <- d[abs(d$x_um - cx) <= w / 2 & abs(d$y_um - cy) <= w / 2, ]
  bz <- get_bg(cx - w / 2, cy - w / 2, cx + w / 2, cy + w / 2, paste0(tag, "_zoom"))
  save_fig(plot_core(z, 9, sprintf("%s  -  core %s, %d um zoom", g, cid, w),
                     sprintf("%s  |  %d of %d cells positive in this window", stain_txt, sum(z$expr > 0), nrow(z)), bz),
           paste0(tag, "_zoom"), width = 11, height = 10, dpi = 300)
}
message("Done.")
