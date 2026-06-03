#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# Node groups
# ============================================================

OUTPUT_NODES = [
    "p50_p65",
    "TNFa",
    "Phenotype_canonical",
    "p100_RelB",
    "p52",
    "BAFF",
    "Phenotype_noncanonical",
]

CANONICAL_OUTPUTS = [
    "p50_p65",
    "TNFa",
    "Phenotype_canonical",
]

NONCANONICAL_OUTPUTS = [
    "p100_RelB",
    "p52",
    "BAFF",
    "Phenotype_noncanonical",
]

NODE_LABELS = {
    "p50_p65": "p50/p65",
    "TNFa": "TNFα",
    "Phenotype_canonical": "Can. phenotype",
    "p100_RelB": "p100:RelB",
    "p52": "p52",
    "BAFF": "BAFF",
    "Phenotype_noncanonical": "Noncan. phenotype",
}


# ============================================================
# Plot style
# ============================================================

plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 22,
        "axes.labelsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 17,
        "figure.titlesize": 24,
        "lines.linewidth": 2.8,

        # EPS/PDF font handling
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


# ============================================================
# Paths
# ============================================================

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent if HERE.parent.name == "scripts" else Path.cwd().resolve()
RESULTS = ROOT / "results"
FIGURE = ROOT / "figure"


# ============================================================
# Basic utilities
# ============================================================

def read_csv(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None or not path.exists():
        return None
    return pd.read_csv(path)


def read_timecourse(path: Optional[Path]) -> Optional[pd.DataFrame]:
    df = read_csv(path)
    if df is None or "time" not in df.columns:
        return None
    return df


def stem_from_mean(path: Path) -> str:
    return path.name.replace("_mean.csv", "")


def stem_from_summary(path: Path) -> str:
    return path.name.replace("_summary.csv", "")


def sd_path(mean_path: Path) -> Optional[Path]:
    p = mean_path.with_name(mean_path.name.replace("_mean.csv", "_sd.csv"))
    return p if p.exists() else None


def available(df: pd.DataFrame, nodes: List[str]) -> List[str]:
    return [n for n in nodes if n in df.columns]


def nice_node(node: str) -> str:
    return NODE_LABELS.get(node, node)


def format_dose(dose: float) -> str:
    return str(int(dose)) if float(dose).is_integer() else f"{dose:g}"


def dose_from_name(name: str) -> Optional[float]:
    m = re.search(r"fucox(?:anthin)?_?(\d+(?:\.\d+)?)", name.lower())
    return float(m.group(1)) if m else None


def lighten_color(color, amount: float = 0.18):
    """
    Return a lightened version of a matplotlib color.
    This avoids EPS transparency warnings from fill_between alpha.
    """
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple((1.0 - amount) * white + amount * rgb)


def save_png_eps(fig: plt.Figure, outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)

    png_path = outpath.with_suffix(".png")
    eps_path = outpath.with_suffix(".eps")

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", bbox_inches="tight")
    plt.close(fig)

    print("saved", png_path)
    print("saved", eps_path)


# ============================================================
# Split output time-course figures
# ============================================================

def save_split_output_timecourse(mean_path: Path) -> None:
    """
    Create a two-panel output time-course figure:
        left  = canonical outputs
        right = noncanonical outputs

    Layout:
        - y-axis upper limit is expanded to 1.30
        - legends are placed inside each axes around y = 1.16
        - canonical legend is one row
        - noncanonical legend is two rows
        - no global title
        - each panel title includes the run name, e.g. fucox25_long
    """

    df = read_timecourse(mean_path)
    if df is None:
        return

    canonical_nodes = available(df, CANONICAL_OUTPUTS)
    noncanonical_nodes = available(df, NONCANONICAL_OUTPUTS)

    if not canonical_nodes and not noncanonical_nodes:
        return

    sd_df = read_timecourse(sd_path(mean_path)) if sd_path(mean_path) else None
    stem = stem_from_mean(mean_path)

    outdir = FIGURE / stem
    outpath = outdir / f"{stem}_split_outputs.png"

    fig, axes = plt.subplots(1, 2, figsize=(16.5, 6.9), sharey=True)
    t = df["time"]
    x_mid = 0.5 * (float(t.min()) + float(t.max()))

    panel_specs = [
        (axes[0], canonical_nodes, f"{stem}: Canonical outputs"),
        (axes[1], noncanonical_nodes, f"{stem}: Noncanonical outputs"),
    ]

    for ax, nodes, title in panel_specs:
        for node in nodes:
            line, = ax.plot(
                t,
                df[node],
                label=nice_node(node),
                linewidth=2.8,
            )

            if sd_df is not None and node in sd_df.columns:
                light_color = lighten_color(line.get_color(), amount=0.16)
                ax.fill_between(
                    t,
                    df[node] - sd_df[node],
                    df[node] + sd_df[node],
                    color=light_color,
                    linewidth=0,
                )

        ax.set_title(
            title,
            fontsize=22,
            fontweight="normal",
            pad=12,
        )

        ax.set_xlabel("time", fontsize=20, fontweight="normal")

        if ax is axes[0]:
            ax.set_ylabel("normalized activity", fontsize=20, fontweight="normal")
        else:
            ax.set_ylabel("")

        # Reserve a smaller blank legend region.
        # Legend will sit around y = 1.16, within 1.1–1.2.
        ax.set_ylim(-0.02, 1.30)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        ax.tick_params(axis="both", which="major", labelsize=18, length=6, width=1.2)
        ax.tick_params(axis="both", which="minor", length=3, width=0.8)

        ax.grid(True, which="major", linestyle="--", color="#d9d9d9", linewidth=1.0)
        ax.grid(True, which="minor", linestyle=":", color="#eeeeee", linewidth=0.8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        # canonical outputs: 3개라서 한 줄
        # noncanonical outputs: 4개라서 2줄
        if len(nodes) <= 3:
            legend_ncol = len(nodes)
        else:
            legend_ncol = 2

        ax.legend(
            loc="center",
            bbox_to_anchor=(x_mid, 1.16),
            bbox_transform=ax.transData,
            ncol=legend_ncol,
            frameon=False,
            fontsize=18,
            handlelength=2.2,
            columnspacing=1.2,
            labelspacing=0.55,
            borderaxespad=0.0,
        )

    fig.subplots_adjust(
        left=0.07,
        right=0.985,
        bottom=0.13,
        top=0.88,
        wspace=0.08,
    )

    save_png_eps(fig, outpath)


# ============================================================
# Fucoxanthin summary table utilities
# ============================================================

def collect_fucox_summary_files() -> List[Tuple[float, Path]]:
    pairs: List[Tuple[float, Path]] = []

    for control in [
        RESULTS / "fucox0_long_summary.csv",
        RESULTS / "fucox0_summary.csv",
        RESULTS / "combined_long_mc_summary.csv",
        RESULTS / "combined_mc_summary.csv",
        RESULTS / "combined_det_summary.csv",
    ]:
        if control.exists():
            pairs.append((0.0, control))
            break

    for path in sorted(RESULTS.glob("*_summary.csv")):
        d = dose_from_name(stem_from_summary(path))
        if d is not None:
            pairs.append((d, path))

    unique: Dict[float, Path] = {}

    for dose, path in pairs:
        if dose not in unique:
            unique[dose] = path
        else:
            old = unique[dose]
            if "long" in path.name and "long" not in old.name:
                unique[dose] = path

    return sorted(unique.items(), key=lambda x: x[0])


def build_fucox_summary_tables() -> Optional[Dict[str, pd.DataFrame]]:
    dose_files = collect_fucox_summary_files()

    if len(dose_files) < 2:
        print("skip fucox publication figures: not enough fucox summary files")
        return None

    doses = [dose for dose, _ in dose_files]

    tables = {
        "AUC_mean": pd.DataFrame(index=OUTPUT_NODES, columns=doses, dtype=float),
        "AUC_sd": pd.DataFrame(index=OUTPUT_NODES, columns=doses, dtype=float),
        "final_mean": pd.DataFrame(index=OUTPUT_NODES, columns=doses, dtype=float),
        "final_sd": pd.DataFrame(index=OUTPUT_NODES, columns=doses, dtype=float),
    }

    for dose, path in dose_files:
        df = read_csv(path)

        if df is None or "node" not in df.columns:
            continue

        df = df.set_index("node")

        for node in OUTPUT_NODES:
            if node not in df.index:
                continue

            if "auc_mean" in df.columns:
                tables["AUC_mean"].loc[node, dose] = float(df.loc[node, "auc_mean"])

            if "auc_sd" in df.columns:
                tables["AUC_sd"].loc[node, dose] = float(df.loc[node, "auc_sd"])

            if "final_mean" in df.columns:
                tables["final_mean"].loc[node, dose] = float(df.loc[node, "final_mean"])

            if "final_sd" in df.columns:
                tables["final_sd"].loc[node, dose] = float(df.loc[node, "final_sd"])

    return tables


def relative_change_percent(df: pd.DataFrame) -> pd.DataFrame:
    if 0.0 not in df.columns:
        raise ValueError("A 0-dose control is required for relative heatmaps.")

    control = df[0.0].replace(0.0, np.nan)
    rel = df.copy()

    for col in df.columns:
        rel[col] = 100.0 * (df[col] - control) / control

    return rel


# ============================================================
# Heatmap
# ============================================================

def draw_heatmap(
    ax,
    table: pd.DataFrame,
    title: str,
    *,
    vmax: float,
    show_ylabel: bool = True,
):
    plot_table = table.T
    arr = plot_table.to_numpy(dtype=float)

    im = ax.imshow(
        arr,
        aspect="auto",
        cmap="Blues",
        vmin=0.0,
        vmax=vmax,
    )

    ax.set_title(title, pad=14, fontsize=22, fontweight="normal")
    ax.set_xlabel("")

    if show_ylabel:
        ax.set_ylabel(
            "Fucoxanthin concentration",
            fontsize=18,
            fontweight="normal",
        )
    else:
        ax.set_ylabel("")

    ax.set_xticks(np.arange(len(plot_table.columns)))
    ax.set_xticklabels(
        [nice_node(col) for col in plot_table.columns],
        rotation=35,
        ha="right",
        fontsize=15,
        fontweight="normal",
    )

    ax.set_yticks(np.arange(len(plot_table.index)))
    ax.set_yticklabels(
        [format_dose(d) for d in plot_table.index],
        fontsize=15,
        fontweight="normal",
    )

    ax.set_xticks(np.arange(-0.5, len(plot_table.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(plot_table.index), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.4)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if np.isfinite(val):
                color = "white" if val >= 0.55 * vmax else "black"
                ax.text(
                    j,
                    i,
                    f"{val:.1f}",
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=15,
                    fontweight="bold",
                )

    return im


# ============================================================
# Ranked bar plot
# ============================================================

def plot_ranked_reduction_bar(
    auc_reduction: pd.DataFrame,
    final_reduction: pd.DataFrame,
    outpath: Path,
    dose: float = 100.0,
) -> None:
    if dose not in auc_reduction.columns:
        print(f"skip ranked bar: dose {dose} not found")
        return

    auc_vals = auc_reduction[dose].copy()
    final_vals = final_reduction[dose].copy()

    data = pd.DataFrame(
        {
            "node": auc_vals.index,
            "AUC reduction": auc_vals.values,
            "Final activity reduction": final_vals.values,
        }
    )

    data["label"] = data["node"].map(lambda x: nice_node(x))
    data = data.sort_values("AUC reduction", ascending=True)

    y = np.arange(len(data))

    fig, ax = plt.subplots(figsize=(11.8, 6.8))

    ax.barh(
        y - 0.18,
        data["AUC reduction"],
        height=0.34,
        label="AUC reduction",
        alpha=1.0,
    )

    ax.barh(
        y + 0.18,
        data["Final activity reduction"],
        height=0.34,
        label="Final activity reduction",
        alpha=1.0,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(data["label"], fontsize=18)

    ax.set_xlabel(
        f"Reduction at Fucoxanthin {format_dose(dose)} (%)",
        fontsize=20,
        fontweight="normal",
    )

    ax.set_title(
        "Ranked output suppression at high Fucoxanthin dose",
        fontsize=22,
        fontweight="normal",
        pad=14,
    )

    ax.tick_params(axis="both", which="major", labelsize=18, length=6, width=1.2)
    ax.grid(True, axis="x", linestyle="--", color="#d9d9d9", linewidth=1.0)
    ax.legend(loc="lower right", frameon=False, fontsize=17)

    xmax = max(data["AUC reduction"].max(), data["Final activity reduction"].max())
    xoffset = max(1.2, xmax * 0.025)

    for i, (_, row) in enumerate(data.iterrows()):
        ax.text(
            row["AUC reduction"] + xoffset,
            i - 0.18,
            f"{row['AUC reduction']:.1f}",
            va="center",
            ha="left",
            fontsize=16,
            fontweight="bold",
            color="black",
            clip_on=False,
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=1.0,
                pad=1.5,
            ),
        )

        ax.text(
            row["Final activity reduction"] + xoffset,
            i + 0.18,
            f"{row['Final activity reduction']:.1f}",
            va="center",
            ha="left",
            fontsize=16,
            fontweight="bold",
            color="black",
            clip_on=False,
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=1.0,
                pad=1.5,
            ),
        )

    ax.set_xlim(0, xmax * 1.30)

    save_png_eps(fig, outpath)


# ============================================================
# Publication figures
# ============================================================

def save_fucox_publication_figures() -> None:
    tables = build_fucox_summary_tables()

    if tables is None:
        return

    dose_dir = FIGURE / "dose_response"
    dose_dir.mkdir(parents=True, exist_ok=True)

    auc_rel = relative_change_percent(tables["AUC_mean"])
    final_rel = relative_change_percent(tables["final_mean"])

    auc_reduction = -auc_rel
    final_reduction = -final_rel

    auc_reduction = auc_reduction.clip(lower=0.0)
    final_reduction = final_reduction.clip(lower=0.0)

    common_vmax = max(
        float(np.nanmax(auc_reduction.to_numpy(dtype=float))),
        float(np.nanmax(final_reduction.to_numpy(dtype=float))),
        1.0,
    )

    # ------------------------------------------------------------
    # 1) Heatmap
    # ------------------------------------------------------------
    fig = plt.figure(figsize=(21.0, 7.2))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.0, 1.0, 0.045],
        wspace=0.18,
    )

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    im1 = draw_heatmap(
        ax0,
        auc_reduction,
        "AUC reduction (%)",
        vmax=common_vmax,
        show_ylabel=True,
    )

    im2 = draw_heatmap(
        ax1,
        final_reduction,
        "Final activity reduction (%)",
        vmax=common_vmax,
        show_ylabel=False,
    )

    cbar = fig.colorbar(im2, cax=cax)
    cbar.set_label(
        "Reduction vs 0-dose (%)",
        rotation=90,
        fontsize=17,
        fontweight="normal",
    )
    cbar.ax.tick_params(labelsize=14)

    fig.suptitle(
        "Fucoxanthin selectively suppresses canonical NF-κB outputs",
        fontsize=24,
        fontweight="normal",
        y=0.98,
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.94,
        bottom=0.22,
        top=0.84,
        wspace=0.18,
    )

    save_png_eps(fig, dose_dir / "fucoxanthin_heatmap_reduction.png")

    # ------------------------------------------------------------
    # 2) Ranked bar plot
    # ------------------------------------------------------------
    plot_ranked_reduction_bar(
        auc_reduction=auc_reduction,
        final_reduction=final_reduction,
        outpath=dose_dir / "fucoxanthin_ranked_reduction_dose100.png",
        dose=100.0,
    )


# ============================================================
# Index
# ============================================================

def write_index() -> None:
    pngs = sorted(FIGURE.rglob("*.png"))
    epss = sorted(FIGURE.rglob("*.eps"))
    path = FIGURE / "figure_index.txt"

    with open(path, "w", encoding="utf-8") as f:
        f.write("PNG files\n")
        f.write("=========\n")
        for p in pngs:
            f.write(str(p.relative_to(FIGURE)) + "\n")

        f.write("\nEPS files\n")
        f.write("=========\n")
        for p in epss:
            f.write(str(p.relative_to(FIGURE)) + "\n")

    print("saved", path)


# ============================================================
# Main
# ============================================================

def main() -> None:
    if not RESULTS.exists():
        raise SystemExit("results/ does not exist. Run simulations first.")

    FIGURE.mkdir(parents=True, exist_ok=True)

    # Only save publication-relevant time-course figures:
    # fucox0_long, fucox25_long, fucox50_long, fucox100_long
    for dose in [0, 25, 50, 100]:
        mean_path = RESULTS / f"fucox{dose}_long_mean.csv"
        if mean_path.exists():
            save_split_output_timecourse(mean_path)
        else:
            print(f"skip split output: {mean_path} not found")

    # Main publication figures
    save_fucox_publication_figures()

    write_index()

    print("Done. Publication figures saved in", FIGURE)


if __name__ == "__main__":
    main()