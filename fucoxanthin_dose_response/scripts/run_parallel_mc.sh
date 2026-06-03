#!/usr/bin/env bash
# Parallel Monte Carlo runner for the bounded NF-kB simulator.
#
# This script parallelizes independent Monte Carlo chunks across Fucoxanthin doses.
# Each chunk writes its own *_mean.csv, *_sd.csv, and *_summary.csv files under
# results/chunks/.  After all jobs finish, scripts/merge_mc_chunks.py combines
# the chunks into final results/fucox{dose}_long_*.csv files.
#
# Example:
#   MC_TOTAL=1000 N_CHUNKS=8 MAX_JOBS=8 STEPS=500000 DT=0.001 SAVE_EVERY=500 \
#       bash scripts/run_parallel_mc.sh
#
# Optional environment variables:
#   DOSES="0 25 50 100"          doses to simulate
#   SCENARIO="combined"          simulator scenario
#   SIGMA=0.01                   stochastic noise strength
#   MAKE_FIGURES=1               call make figures after merge; set 0 to skip
#   SKIP_EXISTING=0              set 1 to skip chunk if its summary already exists

set -euo pipefail

mkdir -p results/chunks logs

MC_TOTAL=${MC_TOTAL:-1000}
N_CHUNKS=${N_CHUNKS:-8}
STEPS=${STEPS:-500000}
DT=${DT:-0.001}
SAVE_EVERY=${SAVE_EVERY:-500}
PROGRESS_EVERY=${PROGRESS_EVERY:-25}
MAX_JOBS=${MAX_JOBS:-8}
SIGMA=${SIGMA:-0.01}
DOSES=${DOSES:-"0 25 50 100"}
SCENARIO=${SCENARIO:-combined}
MAKE_FIGURES=${MAKE_FIGURES:-1}
SKIP_EXISTING=${SKIP_EXISTING:-0}
BASE_SEED=${BASE_SEED:-10000}

if [ ! -x ./nfkb_sim ]; then
    echo "[error] ./nfkb_sim not found or not executable. Run: make" >&2
    exit 1
fi

if [ "${N_CHUNKS}" -le 0 ]; then
    echo "[error] N_CHUNKS must be positive." >&2
    exit 1
fi

if [ "${MC_TOTAL}" -lt "${N_CHUNKS}" ]; then
    echo "[error] MC_TOTAL must be >= N_CHUNKS." >&2
    exit 1
fi

# Unique run id used for chunk manifest and logs.  Human-readable and sortable.
RUN_ID=${RUN_ID:-$(date +%Y%m%d_%H%M%S)}
MANIFEST="results/chunks/parallel_manifest_${RUN_ID}.csv"
LATEST_MANIFEST="results/chunks/parallel_manifest_latest.csv"

cat > "${MANIFEST}" <<EOF_MANIFEST
dose,chunk,mc,seed,prefix
EOF_MANIFEST

BASE_MC=$((MC_TOTAL / N_CHUNKS))
REMAINDER=$((MC_TOTAL % N_CHUNKS))

cat <<EOF_INFO
Parallel NF-kB Monte Carlo run
==============================
Run ID         = ${RUN_ID}
Scenario       = ${SCENARIO}
Doses          = ${DOSES}
Total MC/dose  = ${MC_TOTAL}
Chunks/dose    = ${N_CHUNKS}
Base MC/chunk  = ${BASE_MC}
Remainder      = ${REMAINDER}
Steps          = ${STEPS}
dt             = ${DT}
Save every     = ${SAVE_EVERY}
Progress every = ${PROGRESS_EVERY}
Max jobs       = ${MAX_JOBS}
Sigma          = ${SIGMA}
Make figures   = ${MAKE_FIGURES}
Skip existing  = ${SKIP_EXISTING}
Manifest       = ${MANIFEST}
EOF_INFO

wait_for_slot() {
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
}

launch_job() {
    local dose="$1"
    local chunk="$2"
    local mc="$3"
    local seed="$4"
    local prefix="$5"
    local logfile="$6"

    if [ "${SKIP_EXISTING}" = "1" ] && [ -f "${prefix}_summary.csv" ]; then
        echo "[skip existing] dose=${dose}, chunk=${chunk}, prefix=${prefix}"
        return 0
    fi

    echo "[launch] dose=${dose}, chunk=${chunk}/${N_CHUNKS}, mc=${mc}, seed=${seed}"
    ./nfkb_sim \
        --scenario "${SCENARIO}" \
        --sigma "${SIGMA}" \
        --mc "${mc}" \
        --weights random \
        --fucox "${dose}" \
        --steps "${STEPS}" \
        --dt "${DT}" \
        --save-every "${SAVE_EVERY}" \
        --progress-every "${PROGRESS_EVERY}" \
        --seed "${seed}" \
        --out-prefix "${prefix}" \
        > "${logfile}" 2>&1 &
}

# Launch jobs.
for dose in ${DOSES}; do
    for chunk in $(seq 1 "${N_CHUNKS}"); do
        wait_for_slot

        if [ "${chunk}" -le "${REMAINDER}" ]; then
            mc=$((BASE_MC + 1))
        else
            mc=${BASE_MC}
        fi

        # Avoid decimal arithmetic in seed if dose contains a decimal point.
        dose_seed=$(echo "${dose}" | tr -cd '0-9')
        if [ -z "${dose_seed}" ]; then dose_seed=0; fi
        seed=$((BASE_SEED + dose_seed * 1000 + chunk))
        safe_dose=$(echo "${dose}" | sed 's/[^A-Za-z0-9_.-]/_/g')
        prefix="results/chunks/fucox${safe_dose}_chunk${chunk}"
        logfile="logs/fucox${safe_dose}_chunk${chunk}.log"

        echo "${dose},${chunk},${mc},${seed},${prefix}" >> "${MANIFEST}"
        launch_job "${dose}" "${chunk}" "${mc}" "${seed}" "${prefix}" "${logfile}"
    done
done

wait

echo "All chunk jobs finished."
cp "${MANIFEST}" "${LATEST_MANIFEST}"

# Merge chunks into final files.
python3 scripts/merge_mc_chunks.py \
    --manifest "${MANIFEST}" \
    --out-dir results

# Generate figures if requested.  Do not fail the whole run if matplotlib is absent.
if [ "${MAKE_FIGURES}" = "1" ]; then
    echo "Generating figures..."
    if make figures; then
        echo "Figures generated."
    else
        echo "[warning] make figures failed. CSV results were still generated." >&2
    fi
fi

echo "Done. Final merged results are in results/."
