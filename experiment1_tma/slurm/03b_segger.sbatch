#!/bin/bash
## 03b segger as a SLURM batch job: runs without a browser or an open terminal.
## Submit from Scripts_tma/:   sbatch slurm/03b_segger.sbatch TMA_ORGAN /path/to/output-XETG...
## Follow it:                   squeue -u $USER ;  tail -f segger_<jobid>.log
##                              (segger's own output: <project>/data/segmentation/segger_raw/<slide>/segger_segment.log)
#SBATCH --job-name=segger
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=segger_%j.log

SLIDE_ID=${1:?usage: sbatch slurm/03b_segger.sbatch <slide_id> <xenium_dir> [extra 03b options]}
XENIUM_DIR=${2:?usage: sbatch slurm/03b_segger.sbatch <slide_id> <xenium_dir> [extra 03b options]}
shift 2
PROJECT_DIR=${EXP1_PROJECT_DIR:-/work_ikmb/sukmb430/Spatial_TMA/work}     # ADAPT
SEGGER_REPO=${SEGGER_REPO:-/work_ikmb/sukmb430/segger}                    # ADAPT

source ~/.bashrc
conda activate exp1-imaging
nvidia-smi
python 03b_segment_segger.py --slide-id "$SLIDE_ID" --xenium-dir "$XENIUM_DIR" \
    --project-dir "$PROJECT_DIR" --segger-repo "$SEGGER_REPO" "$@"
