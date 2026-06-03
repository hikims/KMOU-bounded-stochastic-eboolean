# NF-kB bounded simulator v4: long-run version

This version adds long-run support for train/overnight simulations.

## New options

- `--save-every N`: integrate every step but write CSV rows only every `N` steps.
- `--progress-every N`: print progress every `N` Monte Carlo replicates.

The AUC and final summary are computed over the full integration horizon, not only saved rows.

## Recommended train workflow

Compile:

```bash
make clean
make
```

Quick long-run smoke test:

```bash
make run_long_test
make figures
```

Moderate train run, good first production run:

```bash
caffeinate -i bash scripts/run_train_jobs.sh moderate
```

Full heavier train run:

```bash
caffeinate -i bash scripts/run_train_jobs.sh full
```

If you do not want to use `caffeinate`, simply run:

```bash
bash scripts/run_train_jobs.sh moderate
```

## Makefile variables

You can override long-run settings:

```bash
make run_long_mc LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
```

```bash
make run_long_fucox LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
```

## Output

Long-run CSV files are saved in `results/`, logs in `logs/`, and figures in `figure/`.

Important long-run outputs:

- `results/combined_long_mc_summary.csv`
- `results/fucox0_long_summary.csv`
- `results/fucox25_long_summary.csv`
- `results/fucox50_long_summary.csv`
- `results/fucox100_long_summary.csv`
- `figure/dose_response/fucoxanthin_dose_response_AUC.png`
- `figure/dose_response/fucoxanthin_dose_response_final_activity.png`
