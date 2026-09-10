# Parameter sensitivity and revision analyses

This directory contains the numerical analyses that were added during peer
review of

> **A Bounded Stochastic Extended Boolean Model for Signal Amplification, Noise
> Propagation, and Feedback Oscillation in Cell Signaling Pathways**

The code reuses the NF-κB network, edge list, Fucoxanthin-like perturbation,
Waller–Kraft operator, and projected Euler–Maruyama update from
[`../fucoxanthin_dose_response/`](../fucoxanthin_dose_response/). It adds
explicit parameter controls, replicate-level outputs, paired Monte Carlo
comparisons, and projection diagnostics, together with two validation
experiments that share the same build.

| Component | Manuscript location | Main outputs |
|---|---|---|
| `nfkb_revision_sim` + `scripts/run_parameter_sensitivity.py` + `scripts/plot_parameter_sensitivity.py` | Section 4.3.4; Figures 11–12; Tables 6–7 | one-at-a-time scans of the activation rates, deactivation rates, noise amplitude, and Waller–Kraft parameter `r`: matched-control AUC and final-activity suppression, paired-bootstrap 95% CIs, ensemble CVs, projection diagnostics |
| `low_activation_validation` + `scripts/plot_low_activation_validation.py` | Section 4.2.1; Figure 4 | full bounded two-node motif versus its affine low-activation SDE and closed moment ODEs: mean/CV trajectories and approximation errors as activity increases |
| `dt_convergence` + `scripts/plot_dt_convergence.py` | Section 4.3.5; Figure 13 | coupled time-step refinement of the full NF-κB model: differences from the finest step, paired pathwise RMSE, and pre-projection overshoot rates |

The analyses are kept separate from the original figure-generating directories
so that the history of the originally submitted code remains easy to audit.
The regression test in `tests/` verifies that, with every new multiplier set to
one, `nfkb_revision_sim` reproduces the original `nfkb_sim` output.

## Directory layout

```text
parameter_sensitivity/
├── Makefile
├── README.md
├── requirements.txt
├── config/
│   └── sensitivity_grid.csv           # one-at-a-time parameter grid (17 settings)
├── src/
│   ├── nfkb_model.hpp / nfkb_model.cpp   # NF-κB network with parameter controls
│   ├── main_nfkb.cpp                  # nfkb_revision_sim
│   ├── low_activation_validation.cpp  # low_activation_validation
│   └── dt_convergence.cpp             # dt_convergence
├── scripts/
│   ├── run_parameter_sensitivity.py   # launches, merges, and bootstraps the scan
│   ├── plot_parameter_sensitivity.py  # Figures 11–12 and the Table 7 CSV/LaTeX
│   ├── plot_low_activation_validation.py   # Figure 4
│   └── plot_dt_convergence.py         # Figure 13
├── tests/
│   └── check_original_equivalence.py  # regression check against ../fucoxanthin_dose_response
├── reference_outputs/                 # summary tables of the manuscript runs (see below)
├── figures/main_text/                 # PNG versions of Figures 4, 11, 12, and 13
├── results/                           # generated; not tracked (except .gitkeep)
└── logs/                              # generated; not tracked (except .gitkeep)
```

## Requirements

- C++17 compiler (`g++`, `clang++`, or Apple Clang)
- GNU/BSD `make`
- Python 3.9 or later with the packages in `requirements.txt`
  (`numpy`, `pandas`, `matplotlib`)

```bash
python3 -m venv .venv
source .venv/bin/activate                 # Windows: .venv\Scripts\activate
python3 -m pip install -r requirements.txt
```

## Build and test

```bash
make                  # builds nfkb_revision_sim, low_activation_validation, dt_convergence
make check            # compiles the Python scripts and checks --help of each executable
make regression-check # builds ../fucoxanthin_dose_response/nfkb_sim and compares outputs
make smoke            # short end-to-end run of all three workflows (minutes)
```

The smoke profile confirms that compilation, simulation, merging, bootstrap
calculation, and plotting work. Its sample size and time horizon are far too
small for scientific interpretation.

By default the plotting scripts write only the main-text figures listed below.
Pass `--include-extra-figures` to any plotting script to also produce the
diagnostic panels that were used during the analysis.

| Script output | Manuscript item |
|---|---|
| `Figure_A_low_activation_T20.{png,eps}` | Figure 4 |
| `Figure_B_NFkB_sensitivity_T500.{png,eps}` | Figure 11 |
| `Figure_B2_NFkB_noise_scan_T500.{png,eps}` | Figure 12 |
| `Figure_C_time_step_refinement_T20.{png,eps}` | Figure 13 |
| `noise_scan_table.{csv,tex}` | Table 7 (its `sigma = 0.01` rows are Table 6) |

---

## 1. NF-κB parameter sensitivity

### Parameter grid

`config/sensitivity_grid.csv` varies one family at a time around the baseline
(`alpha_scale = beta_scale = sigma_scale = 1`, `r = 0.75`):

- global activation-rate multiplier `alpha_scale`: `0.5, 0.75, 1.0, 1.25, 1.5`;
- global deactivation-rate multiplier `beta_scale`: `0.5, 0.75, 1.0, 1.25, 1.5`;
- noise multiplier `sigma_scale`: `0, 1, 10, 30, 70`, applied to the baseline
  amplitude `sigma_base = 0.01`, so the effective amplitudes are
  `sigma = 0, 0.01, 0.1, 0.3, 0.7`;
- Waller–Kraft parameter `r`: `0.5, 0.625, 0.75, 0.875, 1.0`.

This gives 17 unique settings (the baseline is shared by all four families).
The alpha and beta scans are **global family-level perturbations**: all node
activation rates or all node deactivation rates are multiplied together. They
are not a ranking of individual-node sensitivities. `sigma = 0` is a
deterministic reference in which only the sampled interaction weights vary
between replicates.

### Paired simulation design

Each setting is run at the untreated control (`C = 0`) and at the high
perturbation level (`C = 100`). The random seed depends only on the Monte Carlo
chunk index, not on the setting or the perturbation level, so matching
replicates across settings and levels use the same sampled interaction weights
and the same Brownian streams. Suppression is computed within each setting,

```text
100 × [1 - metric(C = 100; theta) / metric(C = 0; theta)]
```

which prevents a change in baseline activity from being misclassified as a
perturbation response. Confidence intervals are pointwise 95% percentile
intervals from a paired bootstrap over replicate indices (2,000 resamples in
the production profile). Ensemble CVs are reported as point estimates.

The scan uses random edge weights (`--weights random`), matching the ensemble
interpretation of the original implementation; the sampled realization is held
fixed during a single trajectory.

### Waller–Kraft caveat

In the original NF-κB code each individual edge stores a one-element weight
vector. For a one-element vector `min = max`, so `r` cannot change that edge's
Waller–Kraft value. In this network implementation `r` still affects the
Waller–Kraft pooling of **multiple active outgoing edge strengths** used in the
downstream-associated inactivation term, provided the weights are
heterogeneous. For that reason the `r` scan uses random weights; with identical
fixed scalar weights the scan is numerically degenerate, and the executable
prints a warning in that case. The original gate construction is retained
rather than silently redefined.

### Profiles and the manuscript run

```bash
make sensitivity-smoke                          # workflow test only
make sensitivity-pilot                          # T = 50, 100 replicates
make sensitivity-production RUN_ID=production_v2   # T = 500, 200 paired replicates
make sensitivity-figures RUN_ID=production_v2      # Figures 11–12 and Table 7
```

The production profile uses the manuscript time step (`dt = 0.001`) and
horizon (`T = 500`) with 200 paired replicates per setting and perturbation
level (`--chunks 8`, `--jobs 8`, `--bootstrap 2000`, `--base-seed 41000`).
The manuscript results are the run `production_v2`. Because every seed is fixed
(`--base-seed 41000` for the chunks, `--bootstrap-seed 99173` for the
bootstrap), a fresh `make sensitivity-production` with the current grid
reproduces the same numbers on the same platform (the C++ standard library's
`normal_distribution` is implementation-defined, so other platforms may differ
in the last digits); `production_v2` itself was obtained by extending an
earlier run with the original narrower noise grid,

```bash
make sensitivity-extend RUN_ID=production_v2 SOURCE_RUN_ID=production_v1
```

which copies completed chunk outputs of `SOURCE_RUN_ID` into `RUN_ID`, simulates
only the new settings, and then merges and bootstraps all settings.

Custom runs are possible, for example

```bash
python3 scripts/run_parameter_sensitivity.py \
  --profile production --run-id production_v3 --exe ./nfkb_revision_sim \
  --doses 0,25,50,100 --mc 500 --chunks 20 --jobs 10
```

and selected settings can be run alone (the baseline is added automatically):

```bash
python3 scripts/run_parameter_sensitivity.py \
  --profile pilot --exe ./nfkb_revision_sim \
  --conditions alpha_0p50,alpha_1p50,sigma_0p00,sigma_70p00,r_0p500,r_1p000
```

The runner is restartable: completed chunk outputs are skipped unless
`--no-skip-existing` is given, `--run-only` launches jobs without merging, and
`--merge-only` rebuilds the tables from an existing job manifest.

### Output files of a completed run

- `sensitivity_node_metrics.csv`: final-activity and AUC mean, SD, SE, and CV
  per setting, perturbation level, and node;
- `sensitivity_replicates.csv`: replicate-level values and pairing keys;
- `sensitivity_effects.csv`: matched-control reductions with paired 95% CIs;
- `sensitivity_relative_to_baseline.csv`: raw-output changes from the baseline
  setting;
- `sensitivity_diagnostics.csv`: pre-projection predictor violations and
  overshoot magnitudes;
- `job_manifest.csv`, `run_settings.csv`, `conditions_used.csv`, and the
  per-setting merged files.

`plot_parameter_sensitivity.py` reads `sensitivity_effects.csv`,
`sensitivity_node_metrics.csv`, `sensitivity_diagnostics.csv`, and
`run_settings.csv` from `--run-dir` and writes the heatmap figure, the
noise-scan figure, and `noise_scan_table.{csv,tex}` to `--out-dir`.

---

## 2. Low-activation approximation validation

The bounded two-node motif is

```text
dX1 = [u(1 - X1) - d1 X1] dt + sigma X1(1 - X1) dW1,
dX2 = [g X1(1 - X2) - d2 X2] dt + sigma X2(1 - X2) dW2,
```

and its first-order expansion at `(X1, X2) = (0, 0)` is

```text
dY1 = [u - (u + d1) Y1] dt + sigma Y1 dW1,
dY2 = [g Y1 - d2 Y2] dt + sigma Y2 dW2.
```

The program compares the full bounded SDE, the affine linearized SDE, and the
closed first- and second-moment ODEs. The two stochastic systems receive
identical Brownian increments. Increasing `u` moves the system from a
low-activity regime toward saturation, so that the approximation error is
reported rather than assumed to be small at all activity levels.

Manuscript run (`T = 20`, `dt = 0.001`, 5,000 paired paths per input,
`sigma = 0.7`, `g = d1 = d2 = 1`, inputs `0.01, 0.03, 0.1, 0.3, 1, 3`):

```bash
make low-activation-production
```

or manually

```bash
./low_activation_validation \
  --inputs 0.01,0.03,0.1,0.3,1,3 \
  --dt 0.001 --time 20 --mc 5000 --save-every 100 \
  --sigma 0.7 --gain 1 --d1 1 --d2 1 \
  --out-dir results/low_activation

python3 scripts/plot_low_activation_validation.py \
  --data-dir results/low_activation \
  --out-dir figures/low_activation
```

Outputs include time courses, final summaries, full-versus-affine mean errors,
moment-ODE-versus-affine Monte Carlo errors, bounded-predictor overshoots, and
the fraction of linearized trajectories leaving `[0, 1]`.

---

## 3. Coupled time-step refinement

`dt_convergence` generates Brownian increments at the finest requested time
step and sums them to obtain every coarser increment. All step sizes therefore
use the same interaction weights and the same stochastic paths, which is
stronger than reusing an integer seed with different step sizes.

The manuscript figure (Figure 13) uses the `T = 20` check with 100 paired
replicates, steps `0.002, 0.001, 0.0005`, and both `C = 0` and `C = 100`:

```bash
make dt-pilot
```

Full-horizon (`T = 500`) checks are provided as separate targets because they
are expensive; run them individually or on a compute server and then plot both
prefixes:

```bash
make dt-production-control
make dt-production-high

python3 scripts/plot_dt_convergence.py \
  --prefix results/dt_convergence/control \
  --prefix results/dt_convergence/high_dose \
  --out-dir figures/dt_convergence
```

The convergence tables report final-activity and AUC differences from the
finest step, paired pathwise RMSE against the finest-step trajectory, final-CV
differences, and the frequency and magnitude of predictor values outside
`[0, 1]` before projection. The last diagnostic matters because post-projection
values are always bounded; checking only the projected trajectory would hide
how often the numerical safeguard was active.

---

## Reference outputs

`reference_outputs/` holds the summary tables of the runs reported in the
manuscript, so that the tables can be checked and the figures regenerated
without repeating the simulations:

- `parameter_sensitivity/production_v2/`: `run_settings.csv`,
  `conditions_used.csv`, `job_manifest.csv`, `sensitivity_effects.csv`,
  `sensitivity_node_metrics.csv`, `sensitivity_diagnostics.csv`,
  `sensitivity_relative_to_baseline.csv`, and the generated
  `noise_scan_table.{csv,tex}`. The replicate-level file and the chunk outputs
  are not included because of their size; they are regenerated by the
  production run.
- `low_activation_T20/`: the four CSV files written by
  `low_activation_validation` for the manuscript settings.
- `dt_convergence_T20/`: the `control_*` and `high_dose_*` files written by
  `dt_convergence` for the `T = 20` check.

Regenerate the main-text figures and Table 7 from these files:

```bash
python3 scripts/plot_parameter_sensitivity.py \
  --run-dir reference_outputs/parameter_sensitivity/production_v2 \
  --out-dir figures/main_text
python3 scripts/plot_low_activation_validation.py \
  --data-dir reference_outputs/low_activation_T20 \
  --out-dir figures/main_text
python3 scripts/plot_dt_convergence.py \
  --prefix reference_outputs/dt_convergence_T20/control \
  --prefix reference_outputs/dt_convergence_T20/high_dose \
  --out-dir figures/main_text
```

`figures/main_text/` contains the PNG versions of Figures 4, 11, 12, and 13;
the EPS files submitted with the manuscript are produced by the same commands.

## Reproducibility notes

1. Keep `job_manifest.csv`, `run_settings.csv`, and the exact Git commit hash
   for every result reported in the manuscript.
2. Do not compare a perturbed condition at altered parameters against the
   baseline parameter set. Use the matched control generated for the same
   setting.
3. Report the scanned range. A one-at-a-time scan does not establish global
   robustness to arbitrary joint parameter changes.
4. The production profile is a reproducible starting point, not a guarantee
   that 200 replicates are sufficient for every output. Use the paired
   confidence intervals to justify any larger confirmation run.
5. Generated `results/`, `figures/` (other than `figures/main_text/`), and
   `logs/` are excluded from version control.
