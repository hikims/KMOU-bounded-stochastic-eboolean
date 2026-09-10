#!/usr/bin/env python3
"""Plot the bounded-versus-linear low-activation validation experiment.

By default this script generates only the main-text figure used in the
revision (Figure A). Pass --include-extra-figures to additionally create the
older diagnostic panels and representative time-course figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Embed text as TrueType (Type 42) so that EPS/PDF text remains editable and
# searchable, matching plot_parameter_sensitivity.py.
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["pdf.fonttype"] = 42


MODEL_LABELS = {
    "bounded_mc": "Full bounded SDE (Monte Carlo)",
    "linear_mc": "Affine linearized SDE (Monte Carlo)",
    "moment_ode": "Linear moment ODE",
}
MODEL_MARKERS = {"bounded_mc": "o", "linear_mc": "s", "moment_ode": "^"}
MODEL_LINEWIDTHS = {"bounded_mc": 1.8, "linear_mc": 1.2, "moment_ode": 1.2}
ERROR_SERIES = (
    ("linear_vs_bounded_rel_error", "bounded_mean", "o", "Linearized SDE vs bounded SDE"),
    ("moment_vs_linear_rel_error", "linear_mc_mean", "s", "Moment ODE vs linearized SDE"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("results/low_activation"))
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--formats", default="png,eps")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--representative-inputs",
        default=None,
        help="comma-separated input levels; default: minimum, middle, maximum",
    )
    parser.add_argument(
        "--include-extra-figures",
        action="store_true",
        help="also render legacy standalone diagnostic panels and time-course figures",
    )
    return parser.parse_args()


def save_figure(fig: plt.Figure, stem: Path, formats: list[str], dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = dpi
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)


def nearest_available(requested: list[float], available: np.ndarray) -> list[float]:
    selected: list[float] = []
    for value in requested:
        nearest = float(available[np.argmin(np.abs(available - value))])
        if nearest not in selected:
            selected.append(nearest)
    return selected


def format_time_horizon(value: float) -> str:
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def final_comparison(summary: pd.DataFrame, value_col: str, ylabel: str, stem: Path, formats: list[str], dpi: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    for axis, node in zip(axes, (1, 2)):
        sub = summary[summary["node"] == node]
        for model in ("bounded_mc", "linear_mc", "moment_ode"):
            line = sub[sub["model"] == model].sort_values("input")
            axis.plot(
                line["input"],
                line[value_col],
                marker=MODEL_MARKERS[model],
                linewidth=MODEL_LINEWIDTHS[model],
                label=MODEL_LABELS[model],
            )
        axis.set_xscale("log")
        axis.set_xlabel("Input strength")
        axis.set_ylabel(ylabel)
        axis.set_title(f"Node {node}")
        axis.grid(True, linewidth=0.4)
    axes[0].legend(loc="best", fontsize=8, framealpha=1.0)
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def error_figure(errors: pd.DataFrame, out_dir: Path, formats: list[str], dpi: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    for axis, node in zip(axes, (1, 2)):
        line = errors[errors["node"] == node].sort_values("bounded_mean")
        for error_col, x_col, marker, label in ERROR_SERIES:
            axis.scatter(
                line[x_col],
                100.0 * line[error_col],
                marker=marker,
                label=label,
            )
        axis.axhline(0.0, linestyle="--", linewidth=0.8)
        axis.set_xlabel("Corresponding final mean activity")
        axis.set_ylabel("Relative mean error (%)")
        axis.set_title(f"Node {node}")
        axis.grid(True, linewidth=0.4)
    axes[0].legend(loc="best", fontsize=8, framealpha=1.0)
    fig.suptitle("Breakdown of the low-activation approximation")
    save_figure(fig, out_dir / "low_activation_relative_error", formats, dpi)
    plt.close(fig)


def timecourse_figure(
    timecourse: pd.DataFrame,
    representative_inputs: list[float],
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    for input_value in representative_inputs:
        sub_input = timecourse[np.isclose(timecourse["input"], input_value)]
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
        for axis, node in zip(axes, (1, 2)):
            sub = sub_input[sub_input["node"] == node]
            for model in ("bounded_mc", "linear_mc", "moment_ode"):
                line = sub[sub["model"] == model].sort_values("time")
                axis.plot(
                    line["time"],
                    line["mean"],
                    linewidth=MODEL_LINEWIDTHS[model],
                    label=MODEL_LABELS[model],
                )
            axis.set_xlabel("Time")
            axis.set_ylabel("Mean activity")
            axis.set_title(f"Node {node}")
            axis.grid(True, linewidth=0.4)
        axes[0].legend(loc="best", fontsize=8, framealpha=1.0)
        fig.suptitle(f"Low-activation comparison at input = {input_value:g}")
        safe = str(input_value).replace("-", "m").replace(".", "p")
        save_figure(fig, out_dir / f"low_activation_timecourse_input_{safe}", formats, dpi)
        plt.close(fig)


def plot_main_figure(
    summary: pd.DataFrame,
    errors: pd.DataFrame,
    final_time: float,
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 11.5), constrained_layout=True)

    # Top row: final mean versus input strength.
    for axis, node in zip(axes[0], (1, 2)):
        sub = summary[summary["node"] == node]
        for model in ("bounded_mc", "linear_mc", "moment_ode"):
            line = sub[sub["model"] == model].sort_values("input")
            axis.plot(
                line["input"],
                line["final_mean"],
                marker=MODEL_MARKERS[model],
                linewidth=MODEL_LINEWIDTHS[model],
                label=MODEL_LABELS[model],
            )
        axis.set_xscale("log")
        axis.set_xlabel("Input strength")
        axis.set_ylabel("Final mean activity")
        axis.set_title(f"Node {node}")
        axis.grid(True, linewidth=0.4)

    # Middle row: final CV versus input strength.
    for axis, node in zip(axes[1], (1, 2)):
        sub = summary[summary["node"] == node]
        for model in ("bounded_mc", "linear_mc", "moment_ode"):
            line = sub[sub["model"] == model].sort_values("input")
            axis.plot(
                line["input"],
                line["final_cv"],
                marker=MODEL_MARKERS[model],
                linewidth=MODEL_LINEWIDTHS[model],
            )
        axis.set_xscale("log")
        axis.set_xlabel("Input strength")
        axis.set_ylabel("Final coefficient of variation")
        axis.grid(True, linewidth=0.4)

    # Bottom row: approximation errors shown as points only.
    for axis, node in zip(axes[2], (1, 2)):
        line = errors[errors["node"] == node].sort_values("bounded_mean")
        for error_col, x_col, marker, label in ERROR_SERIES:
            axis.scatter(
                line[x_col],
                100.0 * line[error_col],
                marker=marker,
                label=label,
            )
        axis.axhline(0.0, linestyle="--", linewidth=0.8)
        axis.set_xlabel("Corresponding final mean activity")
        axis.set_ylabel("Relative mean error (%)")
        axis.grid(True, linewidth=0.4)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        axes[0, 0].legend(handles, labels, loc="best", fontsize=8, framealpha=1.0)
    bottom_handles, bottom_labels = axes[2, 0].get_legend_handles_labels()
    if bottom_handles:
        axes[2, 0].legend(bottom_handles, bottom_labels, loc="best", fontsize=8, framealpha=1.0)

    fig.suptitle(
        f"Low-activation validation (T = {format_time_horizon(final_time)})",
        fontsize=16,
    )
    stem = out_dir / f"Figure_A_low_activation_T{format_time_horizon(final_time)}"
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir or (args.data_dir / "figures")
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    summary = pd.read_csv(args.data_dir / "low_activation_summary.csv")
    errors = pd.read_csv(args.data_dir / "low_activation_errors.csv")
    timecourse = pd.read_csv(args.data_dir / "low_activation_timecourse.csv")
    final_time = float(timecourse["time"].max())

    plot_main_figure(summary, errors, final_time, out_dir, formats, args.dpi)

    if args.include_extra_figures:
        final_comparison(
            summary,
            "final_mean",
            "Final mean activity",
            out_dir / "low_activation_final_mean",
            formats,
            args.dpi,
        )
        final_comparison(
            summary,
            "final_cv",
            "Final coefficient of variation",
            out_dir / "low_activation_final_cv",
            formats,
            args.dpi,
        )
        error_figure(errors, out_dir, formats, args.dpi)

        available = np.sort(timecourse["input"].astype(float).unique())
        if args.representative_inputs:
            requested = [float(item.strip()) for item in args.representative_inputs.split(",") if item.strip()]
        else:
            requested = [float(available[0]), float(available[len(available) // 2]), float(available[-1])]
        timecourse_figure(timecourse, nearest_available(requested, available), out_dir, formats, args.dpi)

    print(f"Low-activation figures saved to {out_dir}")


if __name__ == "__main__":
    main()
