# NF-kB bounded SDE simulator v3: branch-isolation version

This version adds branch-isolation options for topology validation.

## Compile

```bash
make clean
make
```

## Run biological crosstalk deterministic simulations

```bash
make run_det
```

This generates:

- `results/canonical_det_*.csv`
- `results/noncanonical_det_*.csv`
- `results/combined_det_*.csv`

## Run branch-isolation deterministic simulations

```bash
make run_iso_det
```

This generates:

- `results/canonical_iso_det_*.csv`
- `results/noncanonical_iso_det_*.csv`
- `results/combined_iso_det_*.csv`

In this mode:

- canonical-only suppresses canonical-to-noncanonical leakage.
- noncanonical-only suppresses noncanonical-to-canonical leakage.
- combined keeps biological crosstalk intact.

## Generate figures

```bash
make figures
```

Figures are saved to `figure/`.

Important folders:

- `figure/comparison/` for biological crosstalk comparison after `make run_det`
- `figure/comparison_iso/` for branch-isolation comparison after `make run_iso_det`

## Monte Carlo and Fucoxanthin

```bash
make run_mc
make run_fucox
make figures
```

## New command-line options

```bash
--isolate-branches 0|1
--cross-branch-scale X
```

Example:

```bash
./nfkb_sim --scenario noncanonical --sigma 0 --isolate-branches 1 --cross-branch-scale 0 --out-prefix results/noncanonical_iso_det
```
