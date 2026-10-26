#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One stage of the multiscale readout for one run. Submitted by slurm/run_pipeline.sh, which sets the job name, log file,
# resources and the environment below; do not submit this file by hand unless you set the same variables.
#   STAGE  micro | ikc | leiden | macro | merge
#   RUN    run directory name under RUNS_DIR (its output/ holds nodelist.csv and edgelist.csv)
#   REPO_DIR RUNS_DIR OUT_DIR PYTHON  IKC_PY IKC_K  LEIDEN_RES  TRIANGLES_BIN  MODULES (optional: modules to load)
# REPO_DIR is passed in because Slurm runs a spooled copy of this script, so its own path cannot be used to find the repo.
set -euo pipefail
: "${STAGE:?}" "${RUN:?}" "${REPO_DIR:?}" "${RUNS_DIR:?}" "${OUT_DIR:?}"
if [ -n "${MODULES:-}" ]; then source /etc/profile.d/modules.sh; module load ${MODULES}; fi
PYTHON=${PYTHON:-python3}; IN="$RUNS_DIR/$RUN/output"; OUT="$OUT_DIR/$RUN"
echo "stage=$STAGE run=$RUN host=$(hostname) start=$(date -Is)"
case "$STAGE" in
  micro)  "$PYTHON" "$REPO_DIR/microscale/microscale_readout.py" run "$IN" "$OUT/micro" ${MICRO_ARGS:-} ;;
  ikc)    "$PYTHON" "$REPO_DIR/mesoscale/mesoscale_readout.py" ikc "$IN" "$OUT/meso" --ikc-py "${IKC_PY:?set IKC_PY}" --k "${IKC_K:-20}" --python "$PYTHON" ;;
  leiden) "$PYTHON" "$REPO_DIR/mesoscale/mesoscale_readout.py" leiden "$IN" "$OUT/meso" --resolution "${LEIDEN_RES:-0.1}" ${LEIDEN_ARGS:-} ;;
  macro)  "$PYTHON" "$REPO_DIR/macroscale/macroscale_readout.py" run "$IN" "$OUT/macro" --triangles-bin "${TRIANGLES_BIN:-$REPO_DIR/macroscale/triangles}" ;;
  merge)  "$PYTHON" "$REPO_DIR/multiscale.py" merge "$OUT" ;;
  *) echo "unknown STAGE=$STAGE" >&2; exit 2 ;;
esac
echo "stage=$STAGE run=$RUN finished $(date -Is)"
