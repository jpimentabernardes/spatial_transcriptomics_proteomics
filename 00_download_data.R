## 00_download_data.R
##
## Fallback / complement to postBuild: run this INSIDE a live Binder session
## (or locally) to fetch or refresh data without rebuilding the Docker image.
## Useful for:
##   - swapping in Monday's real Xenium run without touching postBuild/install.R
##     (which would force a full image rebuild for every student)
##   - testing new data URLs quickly during setup
##
## Safe to re-run: skips anything already downloaded.

## R's download.file() defaults to a 60-second timeout regardless of file
## size -- this is what killed your 14.65 GB download partway through, not a
## network problem. Raise it generously for any large downloads in this script.
options(timeout = max(3600, getOption("timeout")))

dir.create("data/metabolomics", showWarnings = FALSE, recursive = TRUE)
dir.create("data/xenium",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed",    showWarnings = FALSE, recursive = TRUE)

download_if_missing <- function(url, dest) {
  if (file.exists(dest)) {
    message("Already present, skipping: ", dest)
    return(invisible(TRUE))
  }
  message("Downloading ", url, " -> ", dest)
  ok <- tryCatch({
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE, method = "libcurl")
    TRUE
  }, error = function(e) {
    warning("Download failed for ", url, ": ", conditionMessage(e))
    FALSE
  })
  invisible(ok)
}

## ---------------------------------------------------------------------------
## 1. Metabolomics (MALDI imzML/ibd) -- fill these in with your hosting URLs
## ---------------------------------------------------------------------------
maldi_imzml_url <- Sys.getenv("MALDI_IMZML_URL", "https://REPLACE-ME/your_run.imzML")
maldi_ibd_url   <- Sys.getenv("MALDI_IBD_URL",   "https://REPLACE-ME/your_run.ibd")

if (!grepl("MALDI_data_Swiss/", maldi_imzml_url)) {
  download_if_missing(maldi_imzml_url, "data/metabolomics/20260820_MST_MLI_MS222_Xeniumslide_mousecolon_uppersection_DAN_200-1000_15um_A37_neg_460x420.imzml")
  download_if_missing(maldi_ibd_url,   "data/metabolomics/20260820_MST_MLI_MS222_Xeniumslide_mousecolon_uppersection_DAN_200-1000_15um_A37_neg_460x420.IBD")
} else {
  message("MALDI_IMZML_URL/MALDI_IBD_URL not set (or you already copied the",
          " files into data/metabolomics/ manually) -- skipping.")
}

## ---------------------------------------------------------------------------
## 2. Transcriptomics (Xenium) -- swap this URL/env var on Monday for the real
##    sequential-section run without touching postBuild.
##
##    DEFAULT CHANGED: the full mouse colon "outs" bundle is ~14.65 GB
##    (full-resolution morphology_focus/ images for a 220k-cell whole-colon
##    run) -- not workable to download live, or to bake into a shared Binder
##    image. Defaulting instead to 10x/Seurat's own small "tiny subset" Xenium
##    tutorial dataset (same file structure, tens of MB), confirmed from
##    Seurat's official vignette. This is what "today's placeholder" should
##    actually be. If you specifically want the mouse colon dataset, see the
##    README section "If you really want the full mouse colon dataset" for a
##    resumable download command and expect it to take a while.
## ---------------------------------------------------------------------------
xenium_url <- Sys.getenv(
  "XENIUM_URL",
  "https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip"
)
xenium_zip <- "data/xenium/xenium_outs.zip"
xenium_dir <- "data/xenium/outs"

if (!dir.exists(xenium_dir)) {
  if (download_if_missing(xenium_url, xenium_zip)) {
    utils::unzip(xenium_zip, exdir = xenium_dir)
    file.remove(xenium_zip)
  }
} else {
  message("Xenium data already present at ", xenium_dir, ", skipping download.")
}

message("00_download_data.R: done. Check data/metabolomics/ and data/xenium/outs/.")