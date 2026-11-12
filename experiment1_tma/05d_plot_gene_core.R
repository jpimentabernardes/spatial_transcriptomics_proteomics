## =============================================================================
## 05d_plot_gene_core.R
##
## One gene in one core, fluorescence style: black background, every cell a
## DAPI-blue dot (about the size of the cell), cells carrying the gene in yellow ->
## white (brighter = more transcripts).
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

plot_core <- function(d, panel_in, title, subtitle) {
  ext <- max(diff(range(d$x_um)), diff(range(d$y_um)), 50)
  area <- if ("cell_area" %in% colnames(d)) d$cell_area else NA
  diam <- ifelse(is.na(area) | area <= 0, 10, 2 * sqrt(area / pi))
  d$size <- pmax(0.05, DOT_SCALE * pmin(diam, 25) * panel_in * 25.4 / ext / 0.75)
  pos <- d[d$expr > 0, ]; pos <- pos[order(pos$expr), ]
  bar <- if (ext > 600) 200 else if (ext > 200) 100 else 50
  ggplot(mapping = aes(x_um, y_um)) +
    geom_point(data = d, aes(size = size), colour = DAPI_BLUE, alpha = 0.75, stroke = 0) +
    geom_point(data = pos, aes(size = size * 1.4, colour = pmin(expr, 5)), stroke = 0) +
    scale_colour_gradientn(colours = c("#FFD000", "#FFF080", "#FFFFFF"), limits = c(1, 5), breaks = 1:5,
                           labels = c("1", "2", "3", "4", "5+"), name = paste(g, "\ntranscripts\nper cell")) +
    scale_size_identity() + scale_y_reverse() + coord_fixed() + theme_void(base_size = 12) +
    annotate("segment", x = min(d$x_um), xend = min(d$x_um) + bar, y = max(d$y_um) + 0.03 * ext,
             yend = max(d$y_um) + 0.03 * ext, colour = "white", linewidth = 0.9) +
    annotate("text", x = min(d$x_um) + bar / 2, y = max(d$y_um) + 0.07 * ext, label = paste(bar, "um"),
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
sub <- sprintf("blue = all cells (%s)  |  yellow -> white = %s+ cells (%d of %d, %.1f%%)",
               format(nrow(d), big.mark = ","), g, sum(d$expr > 0), nrow(d), 100 * mean(d$expr > 0))
tag <- paste0("05d_", make.names(g), "_", make.names(cid))
save_fig(plot_core(d, 9, sprintf("%s  -  core %s  (%s)", g, cid, d$slide_id[1]), sub), tag, width = 11, height = 10, dpi = 300)

## zoom: the ZOOM_UM window with the most positive cells (or the densest one)
w <- ZOOM_UM
if (diff(range(d$x_um)) > 1.5 * w) {
  b <- paste(floor(d$x_um / (w / 2)), floor(d$y_um / (w / 2)))
  score <- tapply(d$expr > 0, b, sum) + 1e-3 * tapply(rep(1, nrow(d)), b, sum)
  best <- strsplit(names(which.max(score)), " ")[[1]]
  cx <- (as.numeric(best[1]) + 0.5) * w / 2; cy <- (as.numeric(best[2]) + 0.5) * w / 2
  z <- d[abs(d$x_um - cx) <= w / 2 & abs(d$y_um - cy) <= w / 2, ]
  save_fig(plot_core(z, 9, sprintf("%s  -  core %s, %d um zoom", g, cid, w),
                     sprintf("%d of %d cells positive in this window", sum(z$expr > 0), nrow(z))),
           paste0(tag, "_zoom"), width = 11, height = 10, dpi = 300)
}
message("Done.")
