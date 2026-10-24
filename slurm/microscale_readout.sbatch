#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --job-name=microscale
#SBATCH --output=slurm-%A_%a.out
# Add your site's --partition / --account lines as needed.
#
# One array task per run. Each RUN is a directory under RUNS_DIR whose output/ holds nodelist.csv and edgelist.csv
# (and, optionally, config.ini next to output/, which lets the checks find the seed files and same_year_citations).
#
#   sbatch --array=0-2 --export=ALL,RUNS="r0 r1 r2",RUNS_DIR=/path/to/runs,OUT_DIR=/path/to/results,PYTHON=python3 \
#       slurm/microscale_readout.sbatch
#
# Needs numpy, pandas, scipy, pyarrow (requirements.txt). Add `module load ...` / venv activation below if needed.
set -euo pipefail
: "${RUNS:?set RUNS to a space-separated list of run directory names}"
: "${RUNS_DIR:?set RUNS_DIR}"; : "${OUT_DIR:?set OUT_DIR}"; PYTHON=${PYTHON:-python3}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
read -ra R <<< "$RUNS"; RUN=${R[${SLURM_ARRAY_TASK_ID:-0}]}
mkdir -p "$OUT_DIR"
$PYTHON "$HERE/microscale/microscale_readout.py" run "$RUNS_DIR/$RUN/output" "$OUT_DIR/$RUN" ${EXTRA_ARGS:-}
