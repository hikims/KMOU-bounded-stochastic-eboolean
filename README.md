# Bounded Stochastic Extended Boolean Model for Cell Signaling Pathways

Code accompanying the manuscript:

> **A Bounded Stochastic Extended Boolean Model for Signal Amplification, Noise
> Propagation, and Feedback Oscillation in Cell Signaling Pathways**
> Minsoo Kim, Department of Data Science, Korea Maritime and Ocean University.
> *(under review at Mathematical Biosciences, 2026; manuscript MBS-D-26-00706)*

This repository contains the simulation and analysis code that reproduces all
computational figures and tables in the paper, including the analyses added
during peer review. The model is a bounded stochastic extended Boolean
formulation in which activation and deactivation/inhibition inputs are separated
and the diffusion is degenerate at both boundaries, so that the biological state
space `[0,1]^N` is positively invariant.

## Repository structure

Each directory is self-contained and builds with its own `Makefile`. See the
`README.md` inside each directory for detailed run instructions.

| Directory | Purpose | Figures in the paper |
|-----------|---------|----------------------|
| [`theory_benchmarks/`](theory_benchmarks/) | Motif-level numerical benchmarks: boundary preservation, moment/CV amplification, and the Hopf-type loop-gain threshold | Figures 1–3 |
| [`branch_validation/`](branch_validation/) | Branch-specific verification of the NF-κB topology under canonical, noncanonical, and combined stimulation | Figures 6–7 |
| [`fucoxanthin_dose_response/`](fucoxanthin_dose_response/) | NF-κB pathway-level Fucoxanthin dose-response simulations and representative output dynamics | Figures 8–10 |
| [`parameter_sensitivity/`](parameter_sensitivity/) | Analyses added in revision: one-at-a-time parameter sensitivity of the NF-κB model with paired-bootstrap intervals and an extended noise scan, validation of the low-activation approximation, and a coupled time-step refinement check | Figures 4, 11–13; Tables 6–7 |

## Requirements

- A C/C++ toolchain with `make` (e.g. `gcc`/`clang`; C++17 for `parameter_sensitivity/`).
- Plotting scripts in each `scripts/` directory (see the per-directory README for
  any Python/Julia dependencies; `parameter_sensitivity/requirements.txt` lists
  the Python packages used there).

## Building and running

From within any of the directories:

```bash
make            # build the simulation executable
make run        # run the simulation (see the local README for options)
make clean      # remove build artifacts
```

Refer to each directory's `README.md` for the exact targets, parameters, and the
scripts that regenerate the figures. `parameter_sensitivity/` uses its own
targets (`make check`, `make regression-check`, `make smoke`, and the
analysis-specific targets listed in its README).

## Reproducing the figures

The compiled executables and raw outputs are **not** tracked in this repository;
they are regenerated from source. Build each module with `make`, run it, and use
the scripts in the corresponding `scripts/` directory to produce the figures in
`figure/` (or `figures/` for `parameter_sensitivity/`). The summary tables
behind Figures 4 and 11–13 and Tables 6–7 are kept in
`parameter_sensitivity/reference_outputs/`, so those figures and tables can be
regenerated without repeating the long simulations.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Kim2026BoundedSEB,
  author  = {Kim, Minsoo},
  title   = {A Bounded Stochastic Extended Boolean Model for Signal Amplification,
             Noise Propagation, and Feedback Oscillation in Cell Signaling Pathways},
  journal = {Mathematical Biosciences},
  year    = {2026},
  note    = {Under review}
}
```

<!-- After archiving a release on Zenodo, add the DOI badge and update the citation:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).
