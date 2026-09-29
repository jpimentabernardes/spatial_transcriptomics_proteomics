# Spatial Multi-Omics Class (MALDI + Xenium + SpaMTP)

> **New:** `experiment1_tma/` holds the research pipeline for Experiment 1:
> multi-organ and multi-tumour TMAs with Xenium, CellScape spatial proteomics on the consecutive
> section, and H&E. It covers QC, segmentation benchmark (10x / cellpose / segger), probe-based
> annotation, Aim 1/2 analyses, registration (GEASO / Spateo) and RNA↔protein integration.
> See `experiment1_tma/README.md`. The class material below is unchanged.
>
> **New:** `tma_toolkit/` is a Python package for large TMAs (e.g. 200 cores). It covers the
> core-metadata sheet, automatic dearraying, selecting and saving sets of cores, and
> core/group comparisons on AnnData. See `tma_toolkit/README.md`.

Teaching repo for a 3-part practical:

1. `scripts/01_spatial_metabolomics_cardinal.R` — MALDI imaging MS with **Cardinal** (+ SpaMTP wrappers)
2. `scripts/02_spatial_transcriptomics_seurat_xenium.R` — 10x Xenium with **Seurat**
3. `scripts/03_integration_spamtp.R` — integrating both with **SpaMTP**

Designed to run unmodified inside an RStudio session launched from **mybinder.org**, following the
same pattern as https://github.com/binder-examples/r (`install.R` + `runtime.txt`, RStudio opens
automatically in the Binder tab).

---

## 0. Read this before you do anything else

### (a) Your Xenium link needs a second click
The dataset page you have
(https://www.10xgenomics.com/datasets/fresh-frozen-mouse-colon-with-xenium-multimodal-cell-segmentation-1-standard)
shows a **Zarr bundle** (`cells.zarr.zip`, `cell_feature_matrix.zarr.zip`, ...) meant for the
Xenium Explorer desktop app. `Seurat::LoadXenium()` cannot read that format. On the same page,
under **"Output and supplemental files"**, look for the standard **Xenium Onboard Analysis "outs"
bundle** (a `*_outs.zip` containing `cell_feature_matrix.h5`, `cells.parquet`, `transcripts.parquet`,
`experiment.xenium`, `morphology_focus/`, etc.). That's the one the script expects. The likely
direct URL, following 10x's normal naming convention, is:

```
https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_mouse_Colon_FF/Xenium_V1_mouse_Colon_FF_outs.zip
```

**Verify this on the dataset page yourself** before class — I could not confirm the exact filename
programmatically (10x's file list is rendered by JavaScript). If it 404s, grab the real link from
the page and drop it into `postBuild` / `scripts/00_download_data.R`.

### (b) Data does not belong in the git repo
GitHub hard-blocks files >100 MB, and a 2.54 GB `.ibd` will make cloning miserable even with Git
LFS. Host the raw metabolomics files (and, ideally, a trimmed Xenium bundle) somewhere Binder can
`wget`/`curl` at build time — Zenodo, OSF, an institutional file server, or a Dropbox/Drive direct
link — and put those URLs in `postBuild` (see below). Only code, `install.R`, `apt.txt`,
`runtime.txt` and `postBuild` should be in the git repo itself.

### (c) Binder's free tier is small and slow to cold-start
mybinder.org gives each user session roughly 1-2 GB RAM and no GPU, and a fresh image build for a
stack this heavy (Cardinal + EBImage + Seurat + SpaMTP + their compiled dependencies: fftw, cairo,
gdal, hdf5...) can take 30-90 minutes the *first* time. After that first build, BinderHub caches
the image and subsequent launches (yours and every student's) are fast — **as long as nothing in
the repo changes**. Practical plan:

- Push this repo and trigger one Binder build **today**, not Thursday morning.
- Don't touch `install.R`/`apt.txt`/`postBuild` again once it builds successfully — any edit
  invalidates the cache and forces a full rebuild.
- Consider seriously downsampling the MALDI data for the *teaching* copy (crop to a smaller ROI,
  or reduce mass resolution/bin width) so peak-picking finishes in a few minutes on ~1-2 GB RAM.
  Keep your full-resolution file for your own research use outside class.
- Test the actual sequence (metabolomics → transcriptomics → integration) end-to-end on Binder,
  not just locally — local RStudio has no memory ceiling; Binder does.

### (d) Serial sections are not the same section
Your MALDI slide and the Xenium slide (arriving Monday) are **adjacent, not identical**, tissue
sections. Even with perfect coordinate alignment, you're matching two physically different slices
of tissue — cell boundaries and exact microanatomy will not correspond 1:1 the way they would on
the same section split for CyAssist-style co-registration. Expect the alignment step
(`AlignSpatialOmics()`, see script 3) to need real manual landmark matching, and expect the
integration results to look noisier than a same-section case study. That's worth saying out loud
to students, not just something to debug around.

### (e) One real unknown in the integration script — flagged, not hidden
SpaMTP's documented `MapSpatialOmics()` examples are all for spot-based Visium data ("lowres
mode", pooling many metabolomics pixels per spot). SpaMTP's own paper reports single-cell
resolution Xenium+MALDI integration is supported (their glioma case study), but I could not pull
the exact current parameter names/defaults for the single-cell ("highres") mode from the public
docs I have access to. Script 3 is written to the documented pattern and flags the one line you
must confirm with `?MapSpatialOmics` (or the "Handling Big Datasets" / "Single Cell Multi-Omics"
vignettes on the SpaMTP site) once you have Seurat ≥5.2 loaded — do this well before Thursday,
ideally the moment Monday's data lands.

---

## Repo layout

```
.
├── README.md
├── install.R              # R/Bioconductor/GitHub package installs (Binder R buildpack)
├── apt.txt                # system libraries needed to compile Cardinal/Seurat/sf deps
├── runtime.txt            # pins the R version / CRAN snapshot date
├── postBuild              # (optional) fetch + unpack data at image build time
├── data/                  # empty in git; populated by postBuild or scripts/00_download_data.R
└── scripts/
    ├── 00_download_data.R                       # fallback if you don't want to bake data into the image
    ├── 01_spatial_metabolomics_cardinal.R
    ├── 02_spatial_transcriptomics_seurat_xenium.R
    └── 03_integration_spamtp.R
```

## Running order

Each script writes its main result object to `data/processed/*.rds` so the next script can just
`readRDS()` it — nobody has to keep an 8 GB R session alive across three scripts.

1. `01_spatial_metabolomics_cardinal.R` → `data/processed/spamtp_metabolomics.rds`
2. `02_spatial_transcriptomics_seurat_xenium.R` → `data/processed/seurat_xenium.rds`
   (today: 10x public colon dataset; Monday: point `XENIUM_DIR` at the real sequential-section run)
3. `03_integration_spamtp.R` reads both `.rds` files and produces the integrated SpaMTP object.

## Before Thursday, in order of urgency

1. Get one successful Binder build (structure above) with placeholder/small data.
2. Confirm the real Xenium `outs.zip` URL from the 10x dataset page.
3. Decide on data hosting for the MALDI files and wire the URL into `postBuild`.
4. Once Monday's Xenium run lands, swap `XENIUM_DIR` in script 2 and re-run scripts 2 and 3 end to
   end on Binder (not just your laptop) to confirm memory/time budget before class.
5. Resolve the flagged `MapSpatialOmics()` mode question in script 3.
