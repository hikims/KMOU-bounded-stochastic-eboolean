# Bounded Stochastic Extended Boolean Model for Cell Signaling Pathways

Code accompanying the manuscript:

> **A Bounded Stochastic Extended Boolean Model for Signal Amplification, Noise
> Propagation, and Feedback Oscillation in Cell Signaling Pathways**
> Minsoo Kim, Department of Data Science, Korea Maritime and Ocean University.
> *(submitted to the Journal of Theoretical Biology, 2026)*

This repository contains the simulation and analysis code that reproduces all
figures in the paper. The model is a bounded stochastic extended Boolean
formulation in which activation and deactivation/inhibition inputs are separated
and the diffusion is degenerate at both boundaries, so that the biological state
space `[0,1]^N` is positively invariant.

## Repository structure

Each directory is self-contained and builds with its own `Makefile`. See the
`README.md` inside each directory for detailed run instructions.

| Directory | Purpose | Figures in the paper |
|-----------|---------|----------------------|
| [`theory_benchmarks/`](theory_benchmarks/) | Motif-level numerical benchmarks: boundary preservation, moment/CV amplification, and the Hopf-type loop-gain threshold | Figures 1–3 |
| [`branch_validation/`](branch_validation/) | Branch-specific verification of the NF-κB topology under canonical, noncanonical, and combined stimulation | Figures 5–6 |
| [`fucoxanthin_dose_response/`](fucoxanthin_dose_response/) | NF-κB pathway-level Fucoxanthin dose-response simulations and representative output dynamics | Figures 7–9 |

## Requirements

- A C/C++ toolchain with `make` (e.g. `gcc`/`clang`).
- Plotting scripts in each `scripts/` directory (see the per-directory README for
  any Python/Julia dependencies).

## Building and running

From within any of the three directories:

```bash
make            # build the simulation executable
make run        # run the simulation (see the local README for options)
make clean      # remove build artifacts
```

Refer to each directory's `README.md` for the exact targets, parameters, and the
scripts that regenerate the figures.

## Reproducing the figures

The compiled executables and raw outputs are **not** tracked in this repository;
they are regenerated from source. Build each module with `make`, run it, and use
the scripts in the corresponding `scripts/` directory to produce the figures in
`figure/`.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Kim2026BoundedSEB,
  author  = {Kim, Minsoo},
  title   = {A Bounded Stochastic Extended Boolean Model for Signal Amplification,
             Noise Propagation, and Feedback Oscillation in Cell Signaling Pathways},
  journal = {Journal of Theoretical Biology},
  year    = {2026},
  note    = {Submitted}
}
```

<!-- After archiving a release on Zenodo, add the DOI badge and update the citation:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).
