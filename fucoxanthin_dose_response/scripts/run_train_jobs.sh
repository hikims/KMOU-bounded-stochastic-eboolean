#!/usr/bin/env bash
set -euo pipefail

# Run this from the project root:
#   bash scripts/run_train_jobs.sh quick
#   bash scripts/run_train_jobs.sh moderate
#   bash scripts/run_train_jobs.sh full
#
# On macOS, to prevent sleep during a train ride:
#   caffeinate -i bash scripts/run_train_jobs.sh moderate

MODE="${1:-moderate}"

mkdir -p results logs figure
make

case "$MODE" in
  quick)
    echo "[train] Quick long-run test"
    make run_long_test
    make figures
    ;;
  moderate)
    echo "[train] Moderate long-run dose response"
    make run_long_mc LONG_MC=200 LONG_STEPS=100000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=20
    make run_long_fucox LONG_MC=200 LONG_STEPS=100000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=20
    make figures
    ;;
  full)
    echo "[train] Full long-run dose response"
    make run_long_mc LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
    make run_long_fucox LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
    make figures
    ;;
  custom)
    echo "[train] Custom mode: use make targets directly with LONG_MC/LONG_STEPS/etc."
    ;;
  *)
    echo "Unknown mode: $MODE"
    echo "Use: quick, moderate, full"
    exit 1
    ;;
esac

echo "[train] Finished. Figures are in ./figure and logs are in ./logs"
