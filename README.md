# Spatial transcriptomics + proteomics on TMAs

Analysis code for tissue microarrays measured with **10x Xenium** (spatial transcriptomics),
**CellScape** (spatial proteomics, consecutive section) and **H&E**.

| Folder | What it is |
|---|---|
| [`experiment1_tma/`](experiment1_tma/README.md) | The Experiment 1 pipeline, R + Python, run step by step: `01` Xenium QC and objects, `02` TMA dearraying (or virtual cores on whole sections, `02a`), `03a`-`03e` segmentation (10x / cellpose / segger) and its benchmark, `04` cell-type annotation, `05` first TMA insights, `06`-`09` CellScape, registration, RNA-protein integration and H&E regions. Settings live in `experiment1_tma/config/`. |
| [`tma_toolkit/`](tma_toolkit/README.md) | Python package for large TMAs (e.g. 200 cores): core-metadata sheet, automatic dearraying, selecting and saving sets of cores, core/group comparisons on AnnData. |
| `tools/make_notebooks.py` | Builds the Jupyter notebooks in `*/notebooks/` from the Python scripts (same code, one cell per step). |

Start with [`experiment1_tma/README.md`](experiment1_tma/README.md): it lists every step, its
inputs and outputs, the order to run them and how to set up the R and conda environments.

The earlier class scripts (MALDI + Xenium + SpaMTP practical) were removed from this branch; they
remain in the git history.
