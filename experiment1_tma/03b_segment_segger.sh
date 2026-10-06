#!/usr/bin/env bash
## =============================================================================
## 03b_segment_segger.sh -- segger (dpeerlab) transcript-based segmentation.
##
## segger treats segmentation as transcript-to-cell link prediction with a
## graph neural network, seeded by the 10x nuclei/cells. It needs a CUDA GPU
## and its own pixi environment (Python 3.11):
##
##   curl -fsSL https://pixi.sh/install.sh | sh
##   git clone https://github.com/dpeerlab/segger.git && cd segger
##   pixi install -e cuda121
##
## CLI verified against the segger README/docs (2026): `segger segment -i <xenium outs> -o <dir>`
## writes segger_segmentation.parquet (one row per transcript: row_index,
## segger_cell_id, x, y, feature_name, filtered, ...). `segger export` turns
## it into boundaries / anndata / spatialdata. Run `segger segment --help`
## once to see the options of the version you installed -- e.g. a scRNA-seq
## reference for gene-gene correlations (--gene-corr-reference-path) can
## help for the custom immune probes.
##
## Python / Jupyter version of this script: 03b_segment_segger.py (notebooks/03b_segment_segger.ipynb).
##
## Usage:  bash 03b_segment_segger.sh <slide_id> <xenium_outs_dir> <project_dir> [segger_repo_dir]
## Then :  python 03b_segger_to_common.py ... (run automatically at the end)
## =============================================================================
set -euo pipefail

SLIDE_ID="$1"
XENIUM_DIR="$2"
PROJECT_DIR="$3"
SEGGER_REPO="${4:-$HOME/segger}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RAW_OUT="${PROJECT_DIR}/data/segmentation/segger_raw/${SLIDE_ID}"
COMMON_OUT="${PROJECT_DIR}/data/segmentation/segger/${SLIDE_ID}"
CORES_CSV="${PROJECT_DIR}/tables/02_cores_all.csv"
mkdir -p "${RAW_OUT}" "${COMMON_OUT}"

echo "== segger segment: ${SLIDE_ID}"
( cd "${SEGGER_REPO}" && pixi run -e cuda121 segger segment -i "${XENIUM_DIR}" -o "${RAW_OUT}" )

echo "== segger export boundaries (cell polygons -> cell areas)"
( cd "${SEGGER_REPO}" && pixi run -e cuda121 segger export boundaries \
    -s "${RAW_OUT}/segger_segmentation.parquet" -o "${RAW_OUT}/export" ) \
  || echo "boundary export failed -- continuing without cell areas"

echo "== convert to the common format read by R/helpers.R"
( cd "${SEGGER_REPO}" && pixi run -e cuda121 python "${CODE_DIR}/03b_segger_to_common.py" \
    --slide-id "${SLIDE_ID}" \
    --xenium-dir "${XENIUM_DIR}" \
    --segger-parquet "${RAW_OUT}/segger_segmentation.parquet" \
    --boundaries "${RAW_OUT}/export/cell_boundaries.parquet" \
    --cores-csv "${CORES_CSV}" \
    --out-dir "${COMMON_OUT}" )
