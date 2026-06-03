# Parallel Monte Carlo / Fucoxanthin dose-response runs

This version can use multiple CPU cores by splitting independent Monte Carlo
replicates into chunks and running the chunks in parallel.  It also runs the
Fucoxanthin dose conditions in parallel.

## Quick test

```bash
make clean
make
make run_parallel_test
```

This should generate small chunk files under `results/chunks/`, merge them into
`results/fucox0_long_*`, `results/fucox25_long_*`, etc., and generate figures.

## Recommended long run on MacStudio

```bash
make run_parallel_mc PAR_MC=1000 PAR_CHUNKS=8 PAR_JOBS=8 PAR_STEPS=500000 PAR_DT=0.001 PAR_SAVE=500 PAR_PROGRESS=25
```

Meaning:

- total MC per dose: 1000
- chunks per dose: 8
- MC per chunk: about 125
- max simultaneously running jobs: 8
- steps per trajectory: 500000
- final time: 500
- save every 500 steps, so about 1000 saved time points per output file

## tmux workflow

```bash
tmux new -s nfkb_parallel
cd /path/to/nfkb_bounded_simulator_v5_parallel
make clean
make
make run_parallel_mc PAR_MC=1000 PAR_CHUNKS=8 PAR_JOBS=8 PAR_STEPS=500000 PAR_DT=0.001 PAR_SAVE=500 PAR_PROGRESS=25
```

Detach with `Ctrl-b` then `d`.

Reattach:

```bash
tmux attach -t nfkb_parallel
```

## Check progress

List running simulator jobs:

```bash
ps aux | grep nfkb_sim | grep -v grep
```

Watch all chunk logs:

```bash
tail -5 logs/fucox*_chunk*.log
```

Watch one log:

```bash
tail -f logs/fucox100_chunk1.log
```

## Final merged outputs

After all chunks finish, the script automatically creates:

```text
results/fucox0_long_mean.csv
results/fucox0_long_sd.csv
results/fucox0_long_summary.csv
results/fucox25_long_mean.csv
...
results/fucox100_long_summary.csv
```

Figures are generated under `figure/` if Python dependencies are available.

## If figure generation fails

The simulation and merge are still valid.  Later run:

```bash
make figures
```

or install plotting packages:

```bash
python3 -m pip install matplotlib pandas numpy
```
