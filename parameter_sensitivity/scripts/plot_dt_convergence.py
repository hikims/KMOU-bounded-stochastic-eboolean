#!/usr/bin/env python3
"""Create the time-step refinement figure.

By default this script generates only the main-text figure used in the
revision (Figure C). Pass --include-extra-figures to additionally create the
legacy convergence and diagnostic plots.
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


DISPLAY = {
    "p50_p65": "p50/p65",
    "TNFa": r"TNF$\alpha$",
    "p52": "p52",
    "BAFF": "BAFF",
}
DEFAULT_NODES = tuple(DISPLAY.keys())
MARKERS = {"p50_p65": "o", "TNFa": "s", "p52": "^", "BAFF": "D"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix",
        action="append",
        type=Path,
        required=True,
        help="prefix of a dt_convergence run, e.g. results/dt_convergence/control",
    )
    parser.add_argument(
        "--label",
        action="append",
        default=None,
        help="optional human-readable label per prefix",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("figures/dt_convergence"))
    parser.add_argument("--nodes", default=",".join(DEFAULT_NODES))
    parser.add_argument("--formats", default="png,eps")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--include-extra-figures",
        action="store_true",
        help="also render legacy line-based convergence diagnostics",
    )
    return parser.parse_args()


def save_figure(fig: plt.Figure, stem: Path, formats: list[str], dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = dpi
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)


def load_run(prefix: Path, label: str | None) -> tuple[str, float, pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(Path(str(prefix) + "_summary.csv"))
    diagnostics = pd.read_csv(Path(str(prefix) + "_diagnostics.csv"))
    final_time = np.nan
    if label is None:
        settings_path = Path(str(prefix) + "_settings.csv")
        if settings_path.exists():
            settings = pd.read_csv(settings_path).iloc[0]
            label = f"{settings['scenario']}, C={settings['fucoxanthin']:g}"
            final_time = float(settings.get("final_time", np.nan))
        else:
            label = prefix.name
    else:
        settings_path = Path(str(prefix) + "_settings.csv")
        if settings_path.exists():
            settings = pd.read_csv(settings_path).iloc[0]
            final_time = float(settings.get("final_time", np.nan))
    return label, final_time, summary, diagnostics


def format_horizon(value: float) -> str:
    value = float(value)
    if np.isnan(value):
        return "unknown"
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def plot_grid(
    runs: list[tuple[str, float, pd.DataFrame, pd.DataFrame]],
    nodes: list[str],
    columns: tuple[str, str],
    ylabels: tuple[str, str],
    titles: tuple[str, str],
    suptitle: str,
    stem: Path,
    formats: list[str],
    dpi: int,
) -> None:
    fig, axes = plt.subplots(
        len(runs), 2, figsize=(12.0, 4.2 * len(runs)), squeeze=False, constrained_layout=True
    )
    for row_index, (run_label, _, summary, _) in enumerate(runs):
        for col_index, (column, ylabel, title) in enumerate(zip(columns, ylabels, titles)):
            axis = axes[row_index, col_index]
            for node in nodes:
                line = summary[summary["node"] == node].sort_values("dt")
                if line.empty:
                    continue
                values = line[column].to_numpy(dtype=float)
                plot_values = np.where(values > 0.0, values, np.nan)
                axis.plot(line["dt"], plot_values, marker="o", label=DISPLAY.get(node, node))
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_xlabel("Time step Δt")
            axis.set_ylabel(ylabel)
            axis.set_title(f"{run_label}: {title}")
            axis.grid(True, linewidth=0.4)
        axes[row_index, 0].legend(loc="best", fontsize=8, framealpha=1.0)
    fig.suptitle(suptitle)
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def plot_main_figure(
    runs: list[tuple[str, float, pd.DataFrame, pd.DataFrame]],
    nodes: list[str],
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    fig, axes = plt.subplots(len(runs), 2, figsize=(12.0, 4.0 * len(runs)), squeeze=False, constrained_layout=True)
    metrics = (
        ("final_mean_relative_difference_vs_finest", "Final activity"),
        ("auc_mean_relative_difference_vs_finest", "AUC"),
    )
    for row_index, (run_label, _, summary, _) in enumerate(runs):
        for col_index, (column, title) in enumerate(metrics):
            axis = axes[row_index, col_index]
            for node in nodes:
                line = summary[summary["node"] == node].sort_values("dt")
                if line.empty:
                    continue
                axis.scatter(
                    line["dt"],
                    100.0 * line[column].to_numpy(dtype=float),
                    marker=MARKERS.get(node, "o"),
                    s=36,
                    label=DISPLAY.get(node, node),
                )
            axis.axhline(0.0, linestyle="--", linewidth=0.8)
            axis.set_xscale("log")
            axis.set_xlabel("Time step Δt")
            axis.set_ylabel("Relative difference from finest step (%)")
            axis.set_title(f"{run_label}: {title}")
            axis.grid(True, linewidth=0.4)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=min(4, len(labels)), framealpha=1.0)
    horizons = [final_time for _, final_time, _, _ in runs if np.isfinite(final_time)]
    final_time = horizons[0] if horizons else np.nan
    fig.suptitle(
        f"Time-step refinement against the finest step (T = {format_horizon(final_time)})",
        fontsize=16,
    )
    stem = out_dir / f"Figure_C_time_step_refinement_T{format_horizon(final_time)}"
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def plot_projection_diagnostic(
    runs: list[tuple[str, float, pd.DataFrame, pd.DataFrame]],
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    fig, axis = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    for label, _, _, diagnostics in runs:
        line = diagnostics.sort_values("dt")
        axis.plot(line["dt"], line["violation_rate"], marker="o", label=label)
    axis.set_xscale("log")
    if any((diagnostics["violation_rate"] > 0.0).any() for _, _, _, diagnostics in runs):
        axis.set_yscale("log")
    axis.set_xlabel("Time step Δt")
    axis.set_ylabel("Pre-projection predictor violation rate")
    axis.set_title("Projection diagnostic")
    axis.grid(True, linewidth=0.4)
    axis.legend(framealpha=1.0)
    save_figure(fig, out_dir / "dt_convergence_projection_rate", formats, dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    labels = args.label or []
    if labels and len(labels) != len(args.prefix):
        raise SystemExit("Provide either no --label values or one label per --prefix")
    runs = [load_run(prefix, labels[i] if labels else None) for i, prefix in enumerate(args.prefix)]
    nodes = [item.strip() for item in args.nodes.split(",") if item.strip()]
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]

    plot_main_figure(runs, nodes, args.out_dir, formats, args.dpi)

    if args.include_extra_figures:
        plot_grid(
            runs,
            nodes,
            (
                "final_mean_relative_difference_vs_finest",
                "auc_mean_relative_difference_vs_finest",
            ),
            (
                "Relative difference from finest step (fraction)",
                "Relative difference from finest step (fraction)",
            ),
            ("Final activity", "AUC"),
            "Coupled Brownian-path time-step convergence",
            args.out_dir / "dt_convergence_relative_difference",
            formats,
            args.dpi,
        )

        plot_grid(
            runs,
            nodes,
            ("final_pathwise_rmse_vs_finest", "auc_pathwise_rmse_vs_finest"),
            (
                "Paired pathwise RMSE from finest step",
                "Paired pathwise RMSE from finest step",
            ),
            ("Final activity", "AUC"),
            "Pathwise convergence under aggregated Brownian increments",
            args.out_dir / "dt_convergence_pathwise_rmse",
            formats,
            args.dpi,
        )

        plot_projection_diagnostic(runs, args.out_dir, formats, args.dpi)

    print(f"Time-step convergence figures saved to {args.out_dir}")


if __name__ == "__main__":
    main()
