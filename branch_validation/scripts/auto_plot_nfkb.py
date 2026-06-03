#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent if HERE.parent.name == "scripts" else Path.cwd().resolve()

RESULTS = ROOT / "results"
FIGURE = ROOT / "figure"
OUTDIR = FIGURE / "branch_validation"


# ============================================================
# Settings
# ============================================================

RUNS = {
    "canonical_iso_det": "Canonical only",
    "noncanonical_iso_det": "Noncanonical only",
    "combined_iso_det": "Combined",
}

PHENOTYPE_LABELS = {
    "Phenotype_canonical": "Canonical phenotype",
    "Phenotype_noncanonical": "Noncanonical phenotype",
}

PHENOTYPE_COLORS = {
    "Phenotype_canonical": "#1f77b4",
    "Phenotype_noncanonical": "#ff7f0e",
}

plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 20,
        "axes.labelsize": 20,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "figure.titlesize": 22,
        "lines.linewidth": 2.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


# ============================================================
# Utilities
# ============================================================

def read_mean_csv(run_name: str) -> Optional[pd.DataFrame]:
    path = RESULTS / f"{run_name}_mean.csv"

    if not path.exists():
        print(f"[skip] missing: {path}")
        return None

    df = pd.read_csv(path)

    if "time" not in df.columns:
        print(f"[skip] no time column: {path}")
        return None

    return df


def save_png_eps(fig: plt.Figure, outpath_without_ext: Path) -> None:
    outpath_without_ext.parent.mkdir(parents=True, exist_ok=True)

    png_path = outpath_without_ext.with_suffix(".png")
    eps_path = outpath_without_ext.with_suffix(".eps")

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", bbox_inches="tight")
    plt.close(fig)

    print("saved", png_path)
    print("saved", eps_path)


def collect_branch_validation_data() -> Dict[str, pd.DataFrame]:
    data: Dict[str, pd.DataFrame] = {}

    for run_name in RUNS:
        df = read_mean_csv(run_name)
        if df is not None:
            data[run_name] = df

    if len(data) < 2:
        raise RuntimeError(
            "Not enough branch-validation CSV files. "
            "Expected files like:\n"
            "  results/canonical_iso_det_mean.csv\n"
            "  results/noncanonical_iso_det_mean.csv\n"
            "  results/combined_iso_det_mean.csv\n"
            "Run the branch-isolation simulations first, e.g. make run_iso_det."
        )

    return data


# ============================================================
# Figure 1: final phenotype bar plot
# ============================================================

def save_final_bar_plot(data: Dict[str, pd.DataFrame]) -> None:
    """
    Main branch-validation figure.

    Shows final canonical/noncanonical phenotype activity under:
        canonical only
        noncanonical only
        combined
    """

    run_names = [r for r in RUNS if r in data]
    x = np.arange(len(run_names))
    width = 0.34

    canonical_vals = []
    noncanonical_vals = []

    for run in run_names:
        df = data[run]

        canonical_vals.append(
            float(df["Phenotype_canonical"].iloc[-1])
            if "Phenotype_canonical" in df.columns else np.nan
        )

        noncanonical_vals.append(
            float(df["Phenotype_noncanonical"].iloc[-1])
            if "Phenotype_noncanonical" in df.columns else np.nan
        )

    fig, ax = plt.subplots(figsize=(10.8, 6.6))

    bars1 = ax.bar(
        x - width / 2,
        canonical_vals,
        width,
        label="Canonical phenotype",
        color=PHENOTYPE_COLORS["Phenotype_canonical"],
        alpha=0.95,
    )

    bars2 = ax.bar(
        x + width / 2,
        noncanonical_vals,
        width,
        label="Noncanonical phenotype",
        color=PHENOTYPE_COLORS["Phenotype_noncanonical"],
        alpha=0.95,
    )

    ax.set_title(
        "Branch-specific phenotype activation",
        fontweight="normal",
        pad=16,
    )

    ax.set_ylabel("Final normalized activity", fontweight="normal")

    ax.set_xticks(x)
    ax.set_xticklabels([RUNS[r] for r in run_names], rotation=0)

    ax.set_ylim(0, 1.05)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)

    ax.tick_params(axis="both", which="major", labelsize=17, length=6, width=1.2)

    # Legend inside the plotting area, around y = 0.8~0.9
    ax.legend(
        loc="center",
        bbox_to_anchor=(0.5, 0.88),
        bbox_transform=ax.transAxes,
        ncol=2,
        frameon=True,
        framealpha=0.85,
        facecolor="white",
        edgecolor="none",
        fontsize=16,
        handlelength=2.4,
        columnspacing=1.2,
        borderpad=0.4,
    )

    # Value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if np.isfinite(h):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + 0.025,
                    f"{h:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=15,
                    fontweight="bold",
                )

    fig.tight_layout()
    save_png_eps(fig, OUTDIR / "branch_validation_final_bar")


# ============================================================
# Figure 2: phenotype time-course validation
# ============================================================

def save_phenotype_timecourse(data: Dict[str, pd.DataFrame]) -> None:
    """
    Supplementary branch-validation figure.

    Shows phenotype trajectories for the branch-isolation runs.
    Legends are placed inside each subplot.
    """

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), sharey=True)

    panels = [
        ("Phenotype_canonical", "Canonical phenotype"),
        ("Phenotype_noncanonical", "Noncanonical phenotype"),
    ]

    for ax, (node, panel_title) in zip(axes, panels):
        line_handles = []
        line_labels = []

        for run_name, run_label in RUNS.items():
            if run_name not in data:
                continue

            df = data[run_name]
            if node not in df.columns:
                continue

            line, = ax.plot(
                df["time"],
                df[node],
                label=run_label,
                linewidth=2.8,
            )

            line_handles.append(line)
            line_labels.append(run_label)

        ax.set_title(panel_title, pad=10, fontweight="normal")
        ax.set_xlabel("time", fontweight="normal")

        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, which="major", linestyle="--", alpha=0.35)
        ax.minorticks_on()
        ax.grid(True, which="minor", linestyle=":", alpha=0.18)

        ax.tick_params(axis="both", which="major", labelsize=17, length=6, width=1.2)
        ax.tick_params(axis="both", which="minor", length=3, width=0.8)

        # Legend inside each subplot, around y = 0.8~0.9
        ax.legend(
            handles=line_handles,
            labels=line_labels,
            loc="center",
            bbox_to_anchor=(0.55, 0.86),
            bbox_transform=ax.transAxes,
            ncol=1,
            frameon=True,
            framealpha=0.85,
            facecolor="white",
            edgecolor="none",
            fontsize=15,
            handlelength=2.4,
            borderpad=0.4,
            labelspacing=0.35,
        )

    axes[0].set_ylabel("normalized activity", fontweight="normal")

    # Main title with reduced gap to subplot titles
    fig.suptitle(
        "Branch validation under isolated stimulation",
        y=0.98,
        fontweight="normal",
    )

    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.14,
        top=0.82,
        wspace=0.25,
    )

    save_png_eps(fig, OUTDIR / "branch_validation_phenotype_timecourse")


# ============================================================
# Index
# ============================================================

def write_index() -> None:
    files = sorted(list(OUTDIR.rglob("*.png")) + list(OUTDIR.rglob("*.eps")))
    index_path = OUTDIR / "branch_validation_index.txt"

    with open(index_path, "w", encoding="utf-8") as f:
        for p in files:
            f.write(str(p.relative_to(OUTDIR)) + "\n")

    print("saved", index_path)


# ============================================================
# Main
# ============================================================

def main() -> None:
    if not RESULTS.exists():
        raise SystemExit("results/ does not exist.")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    data = collect_branch_validation_data()

    save_final_bar_plot(data)
    save_phenotype_timecourse(data)
    write_index()

    print("Done. Branch-validation figures saved in", OUTDIR)


if __name__ == "__main__":
    main()