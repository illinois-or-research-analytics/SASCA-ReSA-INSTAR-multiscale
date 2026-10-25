#!/bin/bash
# Submit the multiscale readout (micro, meso, macro) for one or more simulation runs as a chain of Slurm jobs.
#
#   slurm/run_pipeline.sh [-e pipeline.env] RUN [RUN ...]
#
# For each RUN it submits, in parallel, the stages micro, ikc, leiden and macro (choose a subset with STAGES="micro macro"),
# then a merge job that waits for all of them (afterany, so a failed stage still lets the others merge) and writes
#   OUT_DIR/RUN/summary_all.csv       every measure as "<scale>.<measure>"
#   OUT_DIR/RUN/stages.txt            which stage outputs were found
# Stage logs are in OUT_DIR/RUN/logs/. Settings come from the environment or from the file given with -e (see
# slurm/pipeline.env.example): RUNS_DIR, OUT_DIR, PYTHON, IKC_PY, IKC_K, LEIDEN_RES, MODULES, SBATCH_EXTRA, and per-stage
# resource overrides RES_micro, RES_ikc, RES_leiden, RES_macro, RES_merge (sbatch options as one string).
# Compare runs afterwards:  python multiscale.py compare baseline=OUT_DIR/r0,OUT_DIR/r1  treatment=OUT_DIR/t0
set -euo pipefail
ENV_FILE=""
while getopts "e:" o; do case $o in e) ENV_FILE=$OPTARG ;; *) echo "usage: $0 [-e env_file] RUN [RUN ...]" >&2; exit 2 ;; esac; done
shift $((OPTIND - 1)); [ $# -ge 1 ] || { echo "usage: $0 [-e env_file] RUN [RUN ...]" >&2; exit 2; }
if [ -n "$ENV_FILE" ]; then set -a; source "$ENV_FILE"; set +a; fi
REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)        # correct here: this runs on the login node, not under Slurm
: "${RUNS_DIR:?set RUNS_DIR}" "${OUT_DIR:?set OUT_DIR}"
STAGES=${STAGES:-"micro ikc leiden macro"}
TRIANGLES_BIN=${TRIANGLES_BIN:-$REPO_DIR/macroscale/triangles}
case " $STAGES " in *" ikc "*) : "${IKC_PY:?set IKC_PY (path to IKC.py) or drop ikc from STAGES}"; [ -f "$IKC_PY" ] || { echo "IKC_PY not found: $IKC_PY" >&2; exit 1; } ;; esac
case " $STAGES " in *" macro "*) [ -x "$TRIANGLES_BIN" ] || { echo "building the triangles program (make -C macroscale)"; make -C "$REPO_DIR/macroscale"; } ;; esac

RES_micro=${RES_micro:-"--cpus-per-task=4 --mem=64G --time=01:00:00"}
RES_ikc=${RES_ikc:-"--cpus-per-task=4 --mem=32G --time=02:00:00"}
RES_leiden=${RES_leiden:-"--cpus-per-task=2 --mem=64G --time=04:00:00"}
RES_macro=${RES_macro:-"--cpus-per-task=16 --mem=32G --time=01:00:00"}
RES_merge=${RES_merge:-"--cpus-per-task=1 --mem=4G --time=00:10:00"}

submit() {  # stage run [dependency]
  local stage=$1 run=$2 dep=${3:-} res_var="RES_$1"
  sbatch --parsable --job-name="ms_${stage}_${run}" --output="$OUT_DIR/$run/logs/%x_%j.out" ${dep:+--dependency="$dep"} \
    ${!res_var} ${SBATCH_EXTRA:-} \
    --export="ALL,STAGE=$stage,RUN=$run,REPO_DIR=$REPO_DIR,RUNS_DIR=$RUNS_DIR,OUT_DIR=$OUT_DIR,PYTHON=${PYTHON:-python3},IKC_PY=${IKC_PY:-},IKC_K=${IKC_K:-20},LEIDEN_RES=${LEIDEN_RES:-0.1},TRIANGLES_BIN=$TRIANGLES_BIN,MODULES=${MODULES:-}" \
    "$REPO_DIR/slurm/multiscale_stage.sbatch"
}

for run in "$@"; do
  [ -f "$RUNS_DIR/$run/output/edgelist.csv" ] && [ -f "$RUNS_DIR/$run/output/nodelist.csv" ] || { echo "missing $RUNS_DIR/$run/output/{nodelist,edgelist}.csv" >&2; exit 1; }
  mkdir -p "$OUT_DIR/$run/logs"; ids=""
  for stage in $STAGES; do id=$(submit "$stage" "$run"); ids="$ids:$id"; echo "$run  $stage -> job $id"; done
  mid=$(submit merge "$run" "afterany$ids"); echo "$run  merge -> job $mid (after$ids)"
done
