#!/usr/bin/env python3
"""Create publication-ready figures from a parameter-sensitivity run.

By default this script generates the main-text heatmap (Figure B), a
noise-intensity scan figure (Figure B2) that shows matched-control suppression
with paired 95% CIs and the final-time ensemble CV across the sigma rows, and
a CSV/LaTeX table of the sigma rows (suppression, CIs, CVs, and pre-projection
overshoot diagnostics). Pass --include-extra-figures to additionally create the
legacy line panels, change-from-baseline heatmap, and sigma-CV diagnostic.

Sigma rows are labelled with the effective noise amplitude sigma = q_sigma *
sigma_base (multiplier in parentheses); sigma_base is read from
run_settings.csv.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


# Embed text as TrueType (Type 42) so that labels in the EPS/PDF exports remain
# searchable and extractable instead of being converted to Type 3 outlines.
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["pdf.fonttype"] = 42


FAMILIES = ("alpha", "beta", "sigma", "r")
DEFAULT_SIGMA_BASE = 0.01
BASELINE_X = {"alpha": 1.0, "beta": 1.0, "sigma": 1.0, "r": 0.75}
X_LABEL = {
    "alpha": r"activation-rate multiplier $q_\alpha$",
    "beta": r"deactivation-rate multiplier $q_\beta$",
    "sigma": r"noise-intensity multiplier $q_\sigma$",
    "r": r"Waller--Kraft parameter $r$",
}
DISPLAY = {
    "p50_p65": "p50/p65",
    "TNFa": r"TNF$\alpha$",
    "Phenotype_canonical": "Canonical phenotype",
    "p100_RelB": "p100/RelB",
    "p52": "p52",
    "BAFF": "BAFF",
    "Phenotype_noncanonical": "Noncanonical phenotype",
}
DEFAULT_MAIN_NODES = ("p50_p65", "TNFa", "p52", "BAFF")
DEFAULT_LINE_NODES = ("p50_p65", "TNFa", "p52", "BAFF")
DEFAULT_HEATMAP_NODES = (
    "p50_p65",
    "TNFa",
    "Phenotype_canonical",
    "p100_RelB",
    "p52",
    "BAFF",
    "Phenotype_noncanonical",
)
FAMILY_SYMBOL = {"alpha": "α", "beta": "β", "sigma": "σ", "r": "r"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument(
        "--results-root", type=Path, default=Path("results/parameter_sensitivity")
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--dose", type=float, default=None, help="default: largest nonzero dose")
    parser.add_argument(
        "--nodes",
        default=",".join(DEFAULT_LINE_NODES),
        help="comma-separated node names for optional legacy line figures",
    )
    parser.add_argument(
        "--main-nodes",
        default=",".join(DEFAULT_MAIN_NODES),
        help="comma-separated node names to include in Figure B heatmaps",
    )
    parser.add_argument("--formats", default="png,eps", help="comma-separated formats")
    parser.add_argument("--dpi", type=int, default=400)
    parser.add_argument(
        "--annotation-fontsize",
        type=float,
        default=9.0,
        help="font size of the numbers printed inside heatmap cells",
    )
    parser.add_argument(
        "--annotation-decimals",
        type=int,
        default=1,
        help="number of decimal places shown inside heatmap cells",
    )
    parser.add_argument(
        "--annotation-stroke-width",
        type=float,
        default=1.15,
        help="width of the contrasting outline around heatmap numbers",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="colormap for the main sensitivity heatmap",
    )
    parser.add_argument(
        "--change-cmap",
        default="coolwarm",
        help="diverging colormap for the optional change-from-baseline heatmap",
    )
    parser.add_argument(
        "--include-extra-figures",
        action="store_true",
        help="also render legacy line panels and supplementary diagnostics",
    )
    return parser.parse_args()


def latest_run(root: Path) -> Path:
    candidates = [
        path
        for path in root.iterdir()
        if path.is_dir() and (path / "sensitivity_effects.csv").exists()
    ] if root.exists() else []
    if not candidates:
        raise SystemExit(f"No completed sensitivity run found under {root}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def save_figure(fig: plt.Figure, stem: Path, formats: list[str], dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt.lower() == "png":
            kwargs["dpi"] = dpi
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)


def isclose(series: pd.Series, value: float) -> pd.Series:
    return np.isclose(series.to_numpy(dtype=float), value, rtol=0.0, atol=1e-12)


def family_rows(frame: pd.DataFrame, family: str) -> pd.DataFrame:
    selected = frame[(frame["family"] == family) | (frame["condition_id"] == "baseline")].copy()
    selected["x_value"] = selected["level"].astype(float)
    selected.loc[selected["condition_id"] == "baseline", "x_value"] = BASELINE_X[family]
    return selected.sort_values("x_value")


def format_level(value: float, family: str) -> str:
    if family == "r":
        return f"{value:.3f}"
    return f"{value:.2f}"


def format_sigma(value: float) -> str:
    """Compact representation of an effective noise amplitude (0, 0.01, 0.1, ...)."""

    return f"{float(value):.3g}"


def format_multiplier(value: float) -> str:
    return f"{float(value):g}"


def row_label(family: str, level: float, sigma_base: float) -> str:
    """Heatmap row label.

    Rate families show the multiplier as before.  Sigma rows show the effective
    noise amplitude sigma = level * sigma_base with the multiplier in
    parentheses, because the amplitude, not the multiplier, is what a reader
    compares with the motif benchmarks.
    """

    if family == "sigma":
        return f"σ = {format_sigma(level * sigma_base)} (×{format_multiplier(level)})"
    return f"{FAMILY_SYMBOL[family]} = {format_level(level, family)}"


def format_horizon(value: float) -> str:
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def read_run_setting(settings_path: Path, key: str, default: float) -> float:
    if not settings_path.exists():
        return default
    settings = pd.read_csv(settings_path).iloc[0]
    if key not in settings or pd.isna(settings[key]):
        return default
    return float(settings[key])


def build_row_spec(
    data: pd.DataFrame, sigma_base: float = DEFAULT_SIGMA_BASE
) -> list[tuple[str, float, str]]:
    row_spec: list[tuple[str, float, str]] = []
    for family in FAMILIES:
        family_levels = sorted(
            set(data.loc[data["family"] == family, "level"].astype(float).tolist())
            | {BASELINE_X[family]}
        )
        for level in family_levels:
            row_spec.append((family, float(level), row_label(family, level, sigma_base)))
    return row_spec


def matrix_for_metric(
    effects: pd.DataFrame,
    metric: str,
    dose: float,
    nodes: list[str],
    row_spec: list[tuple[str, float, str]],
) -> np.ndarray:
    data = effects[(effects["metric"] == metric) & isclose(effects["dose"], dose)].copy()
    baseline_rows = data[data["condition_id"] == "baseline"].set_index("node")
    matrix = np.full((len(row_spec), len(nodes)), np.nan, dtype=float)

    for row_index, (family, level, _) in enumerate(row_spec):
        if np.isclose(level, BASELINE_X[family]):
            row_frame = baseline_rows
        else:
            row_frame = data[(data["family"] == family) & isclose(data["level"], level)].set_index("node")
        for col_index, node in enumerate(nodes):
            if node in row_frame.index:
                matrix[row_index, col_index] = float(row_frame.loc[node, "reduction_percent"])
    return matrix


def draw_group_separators(axis: plt.Axes, row_spec: list[tuple[str, float, str]]) -> None:
    boundary_indices: list[int] = []
    prev_family = row_spec[0][0]
    for idx, (family, _, _) in enumerate(row_spec[1:], start=1):
        if family != prev_family:
            boundary_indices.append(idx)
            prev_family = family
    for boundary in boundary_indices:
        axis.axhline(boundary - 0.5, color="white", linewidth=2.0)
        axis.axhline(boundary - 0.5, color="black", linewidth=0.8)


def relative_luminance(rgba: tuple[float, ...] | np.ndarray) -> float:
    """Return WCAG relative luminance for an RGB(A) color."""

    rgb = np.asarray(rgba[:3], dtype=float)
    linear = np.where(
        rgb <= 0.04045,
        rgb / 12.92,
        ((rgb + 0.055) / 1.055) ** 2.4,
    )
    return float(0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2])


def contrasting_text_colors(image: object, value: float) -> tuple[str, str]:
    """Choose text and outline colors with the stronger background contrast."""

    rgba = image.cmap(image.norm(value))
    luminance = relative_luminance(rgba)
    contrast_with_black = (luminance + 0.05) / 0.05
    contrast_with_white = 1.05 / (luminance + 0.05)
    if contrast_with_black >= contrast_with_white:
        return "black", "white"
    return "white", "black"


def annotate_heatmap(
    axis: plt.Axes,
    image: object,
    values: np.ndarray,
    fontsize: float,
    decimals: int,
    stroke_width: float,
) -> None:
    """Print high-contrast values inside heatmap cells."""

    number_format = f"{{:.{max(0, decimals)}f}}"
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            value = values[row, col]
            if not np.isfinite(value):
                continue
            text_color, outline_color = contrasting_text_colors(image, float(value))
            label = axis.text(
                col,
                row,
                number_format.format(value),
                ha="center",
                va="center",
                color=text_color,
                fontsize=fontsize,
                fontweight="bold",
                clip_on=True,
            )
            label.set_path_effects(
                [
                    path_effects.withStroke(
                        linewidth=stroke_width,
                        foreground=outline_color,
                    )
                ]
            )


def draw_cell_grid(axis: plt.Axes, n_rows: int, n_cols: int) -> None:
    """Add subtle cell boundaries without overpowering the group separators."""

    axis.set_xticks(np.arange(-0.5, n_cols, 1.0), minor=True)
    axis.set_yticks(np.arange(-0.5, n_rows, 1.0), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.35, alpha=0.45)
    axis.tick_params(which="minor", bottom=False, left=False)


def plot_main_heatmap(
    effects: pd.DataFrame,
    dose: float,
    nodes: list[str],
    final_time: float,
    out_dir: Path,
    formats: list[str],
    dpi: int,
    annotation_fontsize: float,
    annotation_decimals: int,
    annotation_stroke_width: float,
    cmap: str,
    sigma_base: float = DEFAULT_SIGMA_BASE,
) -> None:
    row_spec = build_row_spec(effects, sigma_base)
    row_labels = [label for _, _, label in row_spec]
    metric_titles = {"auc": "AUC reduction (%)", "final": "Final activity reduction (%)"}
    matrices = {
        metric: matrix_for_metric(effects, metric, dose, nodes, row_spec)
        for metric in ("auc", "final")
    }

    all_values = np.concatenate([matrix.ravel() for matrix in matrices.values()])
    finite = all_values[np.isfinite(all_values)]
    vmin = float(np.nanmin(finite)) if finite.size else 0.0
    vmax = float(np.nanmax(finite)) if finite.size else 1.0
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0

    fig_height = max(10.0, 0.40 * len(row_labels) + 3.0)
    fig, axes = plt.subplots(1, 2, figsize=(16.0, fig_height), constrained_layout=True)
    axes = np.atleast_1d(axes)
    image = None
    for axis, metric in zip(axes, ("auc", "final")):
        matrix = matrices[metric]
        image = axis.imshow(
            matrix,
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            interpolation="nearest",
        )
        axis.set_xticks(
            np.arange(len(nodes)),
            [DISPLAY.get(node, node) for node in nodes],
            rotation=25,
            ha="right",
        )
        axis.set_yticks(np.arange(len(row_labels)), row_labels)
        axis.tick_params(axis="both", labelsize=10)
        axis.set_title(metric_titles[metric], fontsize=13, pad=10)
        axis.set_xlabel("NF-κB output", fontsize=11)
        draw_cell_grid(axis, matrix.shape[0], matrix.shape[1])
        draw_group_separators(axis, row_spec)
        annotate_heatmap(
            axis,
            image,
            matrix,
            fontsize=annotation_fontsize,
            decimals=annotation_decimals,
            stroke_width=annotation_stroke_width,
        )
    axes[0].set_ylabel("Sensitivity condition", fontsize=11)
    if image is not None:
        cbar = fig.colorbar(image, ax=axes, shrink=0.92, pad=0.025)
        cbar.set_label("Reduction relative to matched control (%)", fontsize=11)
        cbar.ax.tick_params(labelsize=10)

    fig.suptitle(
        f"NF-κB sensitivity summary at perturbation level {dose:g} (T = {format_horizon(final_time)})",
        fontsize=16,
    )
    stem = out_dir / f"Figure_B_NFkB_sensitivity_T{format_horizon(final_time)}"
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def plot_reduction_panels(
    effects: pd.DataFrame,
    metric: str,
    dose: float,
    nodes: list[str],
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    data = effects[(effects["metric"] == metric) & isclose(effects["dose"], dose)].copy()
    if data.empty:
        raise SystemExit(f"No {metric} effects found at dose {dose:g}")

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    axes = axes.ravel()
    for axis, family in zip(axes, FAMILIES):
        sub = family_rows(data, family)
        for node in nodes:
            line = sub[sub["node"] == node].sort_values("x_value")
            if line.empty:
                continue
            y = line["reduction_percent"].to_numpy(dtype=float)
            low = line["ci95_low"].to_numpy(dtype=float)
            high = line["ci95_high"].to_numpy(dtype=float)
            lower_error = np.maximum(0.0, y - low)
            upper_error = np.maximum(0.0, high - y)
            yerr = np.vstack([lower_error, upper_error])
            if not np.isfinite(yerr).any():
                yerr = None
            axis.errorbar(
                line["x_value"],
                y,
                yerr=yerr,
                marker="o",
                capsize=3,
                linewidth=1.4,
                label=DISPLAY.get(node, node),
            )
        axis.axvline(BASELINE_X[family], linestyle=":", linewidth=1.0)
        axis.axhline(0.0, linestyle="--", linewidth=0.8)
        axis.set_xlabel(X_LABEL[family])
        axis.set_ylabel("Reduction relative to matched control (%)")
        axis.set_title(f"{family} sensitivity")
        axis.grid(True, linewidth=0.4)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.03), ncol=min(4, len(labels)), framealpha=1.0)
    label = "AUC" if metric == "auc" else "final activity"
    fig.suptitle(f"NF-κB {label} reduction at perturbation level {dose:g}")
    save_figure(fig, out_dir / f"sensitivity_{metric}_reduction", formats, dpi)
    plt.close(fig)


def plot_reduction_change_heatmap(
    effects: pd.DataFrame,
    dose: float,
    out_dir: Path,
    formats: list[str],
    dpi: int,
    annotation_fontsize: float,
    annotation_decimals: int,
    annotation_stroke_width: float,
    cmap: str,
) -> None:
    data = effects[(effects["metric"] == "auc") & isclose(effects["dose"], dose)].copy()
    baseline = data[data["condition_id"] == "baseline"].set_index("node")["reduction_percent"]
    data = data[data["condition_id"] != "baseline"].copy()
    data["change_from_baseline"] = data.apply(
        lambda row: row["reduction_percent"] - baseline.get(row["node"], np.nan), axis=1
    )

    condition_order: list[str] = []
    for family in FAMILIES:
        rows = data[data["family"] == family][["condition_id", "level"]].drop_duplicates()
        rows = rows.sort_values("level")
        condition_order.extend(rows["condition_id"].tolist())
    nodes = [node for node in DEFAULT_HEATMAP_NODES if node in set(data["node"])]
    pivot = data.pivot_table(
        index="condition_id", columns="node", values="change_from_baseline", aggfunc="first"
    ).reindex(index=condition_order, columns=nodes)

    values = pivot.to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    max_abs = float(np.max(np.abs(finite_values))) if finite_values.size else 1.0
    max_abs = max(max_abs, np.finfo(float).eps)
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    fig_height = max(6.5, 0.43 * max(1, len(pivot.index)) + 2.0)
    fig, axis = plt.subplots(figsize=(12.5, fig_height), constrained_layout=True)
    image = axis.imshow(
        values,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )
    axis.set_xticks(
        np.arange(len(nodes)),
        [DISPLAY.get(node, node) for node in nodes],
        rotation=30,
        ha="right",
    )
    axis.set_yticks(np.arange(len(pivot.index)), pivot.index)
    axis.tick_params(axis="both", labelsize=10)
    axis.set_xlabel("NF-κB output", fontsize=11)
    axis.set_ylabel("Sensitivity condition", fontsize=11)
    axis.set_title(
        f"Change in AUC reduction from baseline at perturbation level {dose:g}\n"
        "(percentage points; each condition uses its own matched untreated control)"
    )
    cbar = fig.colorbar(image, ax=axis, pad=0.025)
    cbar.set_label("Change from baseline (percentage points)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    draw_cell_grid(axis, values.shape[0], values.shape[1])
    if values.size <= 180:
        annotate_heatmap(
            axis,
            image,
            values,
            fontsize=annotation_fontsize,
            decimals=annotation_decimals,
            stroke_width=annotation_stroke_width,
        )
    save_figure(fig, out_dir / "sensitivity_auc_reduction_change_heatmap", formats, dpi)
    plt.close(fig)


def plot_noise_cv(
    node_metrics: pd.DataFrame,
    dose: float,
    nodes: list[str],
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    data = node_metrics[(node_metrics["family"] == "sigma") | (node_metrics["condition_id"] == "baseline")].copy()
    data["x_value"] = data["level"].astype(float)
    data.loc[data["condition_id"] == "baseline", "x_value"] = 1.0
    panel_doses = [0.0]
    if dose not in panel_doses:
        panel_doses.append(dose)

    fig, axes = plt.subplots(1, len(panel_doses), figsize=(6.0 * len(panel_doses), 4.5), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, panel_dose in zip(axes, panel_doses):
        sub = data[isclose(data["dose"], panel_dose)]
        for node in nodes:
            line = sub[sub["node"] == node].sort_values("x_value")
            if line.empty:
                continue
            axis.plot(
                line["x_value"],
                line["final_cv"],
                marker="o",
                linewidth=1.4,
                label=DISPLAY.get(node, node),
            )
        axis.axvline(1.0, linestyle=":", linewidth=1.0)
        axis.set_xlabel(X_LABEL["sigma"])
        axis.set_ylabel("Final-time ensemble CV")
        axis.set_title(f"Perturbation level {panel_dose:g}")
        axis.grid(True, linewidth=0.4)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.03), ncol=min(4, len(labels)), framealpha=1.0)
    fig.suptitle("Output variability under noise-intensity scaling")
    save_figure(fig, out_dir / "sensitivity_noise_final_cv", formats, dpi)
    plt.close(fig)


def sigma_rows(frame: pd.DataFrame, sigma_base: float) -> pd.DataFrame:
    """Rows of the sigma family plus the shared baseline, ordered by multiplier."""

    data = frame[(frame["family"] == "sigma") | (frame["condition_id"] == "baseline")].copy()
    data["multiplier"] = data["level"].astype(float)
    data.loc[data["condition_id"] == "baseline", "multiplier"] = BASELINE_X["sigma"]
    data["sigma"] = data["multiplier"] * sigma_base
    return data.sort_values("multiplier")


def sigma_levels(frame: pd.DataFrame, sigma_base: float) -> list[tuple[float, float]]:
    rows = sigma_rows(frame, sigma_base)
    return sorted({(float(m), float(s)) for m, s in zip(rows["multiplier"], rows["sigma"])})


def plot_noise_scan(
    effects: pd.DataFrame,
    node_metrics: pd.DataFrame,
    dose: float,
    nodes: list[str],
    sigma_base: float,
    final_time: float,
    out_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    """Figure B2: matched-control suppression and ensemble CV across the sigma rows.

    The horizontal axis is categorical (one position per tested amplitude) so
    that the deterministic reference sigma = 0 can be shown next to the
    logarithmically spaced amplitudes.
    """

    levels = sigma_levels(effects, sigma_base)
    if len(levels) < 2:
        print("Noise-scan figure skipped: fewer than two sigma levels available")
        return
    positions = np.arange(len(levels))
    tick_labels = [format_sigma(s) for _, s in levels]
    multipliers = [m for m, _ in levels]
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    effect_rows = sigma_rows(effects[isclose(effects["dose"], dose)], sigma_base)
    metric_rows = sigma_rows(node_metrics, sigma_base)
    control_dose = 0.0

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.9), constrained_layout=True)
    metric_titles = {"auc": "AUC reduction (%)", "final": "Final activity reduction (%)"}
    panel_titles = {"auc": "AUC suppression", "final": "Final-activity suppression"}
    for axis, metric in zip(axes[:2], ("auc", "final")):
        sub = effect_rows[effect_rows["metric"] == metric]
        for color, node in zip(colors, nodes):
            line = sub[sub["node"] == node].set_index("multiplier")
            if line.empty:
                continue
            y = np.array([line["reduction_percent"].get(m, np.nan) for m in multipliers], dtype=float)
            low = np.array([line["ci95_low"].get(m, np.nan) for m in multipliers], dtype=float)
            high = np.array([line["ci95_high"].get(m, np.nan) for m in multipliers], dtype=float)
            yerr = np.vstack([np.maximum(0.0, y - low), np.maximum(0.0, high - y)])
            axis.errorbar(
                positions,
                y,
                yerr=None if not np.isfinite(yerr).any() else yerr,
                marker="o",
                capsize=3,
                linewidth=1.4,
                color=color,
                label=DISPLAY.get(node, node),
            )
        if BASELINE_X["sigma"] in multipliers:
            axis.axvline(multipliers.index(BASELINE_X["sigma"]), linestyle=":", linewidth=1.0, color="grey")
        axis.set_xticks(positions, tick_labels)
        axis.set_xlabel("effective noise amplitude σ")
        axis.set_ylabel(metric_titles[metric])
        axis.set_title(f"{panel_titles[metric]} at perturbation level {dose:g}")
        axis.grid(True, linewidth=0.4)

    axis = axes[2]
    for color, node in zip(colors, nodes):
        for panel_dose, style, tag in ((control_dose, "-", "C = 0"), (dose, "--", f"C = {dose:g}")):
            line = metric_rows[isclose(metric_rows["dose"], panel_dose) & (metric_rows["node"] == node)]
            line = line.set_index("multiplier")
            if line.empty:
                continue
            y = np.array([line["final_cv"].get(m, np.nan) for m in multipliers], dtype=float)
            axis.plot(
                positions,
                y,
                linestyle=style,
                marker="o" if style == "-" else "s",
                linewidth=1.4,
                color=color,
                label=f"{DISPLAY.get(node, node)}, {tag}",
            )
    axis.set_xticks(positions, tick_labels)
    axis.set_xlabel("effective noise amplitude σ")
    axis.set_ylabel("Final-time ensemble CV")
    axis.set_title("Replicate variability at T = " + format_horizon(final_time))
    axis.grid(True, linewidth=0.4)
    axis.legend(fontsize=8, ncol=2, framealpha=1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=min(4, len(labels)), framealpha=1.0)
    fig.suptitle(
        f"Noise-intensity scan (σ = q_σ × {format_sigma(sigma_base)}) at perturbation level {dose:g}, "
        f"T = {format_horizon(final_time)}",
        fontsize=14,
    )
    stem = out_dir / f"Figure_B2_NFkB_noise_scan_T{format_horizon(final_time)}"
    save_figure(fig, stem, formats, dpi)
    plt.close(fig)


def write_noise_scan_table(
    effects: pd.DataFrame,
    node_metrics: pd.DataFrame,
    diagnostics: pd.DataFrame | None,
    dose: float,
    nodes: list[str],
    sigma_base: float,
    out_dir: Path,
) -> None:
    """Write the sigma rows as a CSV and a booktabs LaTeX fragment.

    The CSV also carries the pre-projection predictor diagnostics of each
    sigma condition (violation rate and maximum overshoot for the control and
    treated ensembles), which are needed to state whether projection was ever
    active at the larger amplitudes.
    """

    effect_rows = sigma_rows(effects[isclose(effects["dose"], dose)], sigma_base)
    metric_rows = sigma_rows(node_metrics, sigma_base)
    diag_rows = sigma_rows(diagnostics, sigma_base) if diagnostics is not None else None

    records: list[dict[str, object]] = []
    for multiplier, sigma in sigma_levels(effects, sigma_base):
        e_sub = effect_rows[isclose(effect_rows["multiplier"], multiplier)]
        m_sub = metric_rows[isclose(metric_rows["multiplier"], multiplier)]
        diag: dict[str, float] = {}
        if diag_rows is not None:
            d_sub = diag_rows[isclose(diag_rows["multiplier"], multiplier)]
            for panel_dose, tag in ((0.0, "control"), (dose, "treated")):
                row = d_sub[isclose(d_sub["dose"], panel_dose)]
                if row.empty:
                    continue
                row = row.iloc[0]
                diag[f"violation_rate_{tag}"] = float(row.get("violation_rate", np.nan))
                diag[f"max_lower_overshoot_{tag}"] = float(row.get("max_lower_overshoot", np.nan))
                diag[f"max_upper_overshoot_{tag}"] = float(row.get("max_upper_overshoot", np.nan))
        for node in nodes:
            record: dict[str, object] = {
                "sigma_multiplier": multiplier,
                "sigma": sigma,
                "node": node,
                "node_label": DISPLAY.get(node, node),
            }
            for metric in ("auc", "final"):
                row = e_sub[(e_sub["metric"] == metric) & (e_sub["node"] == node)]
                if row.empty:
                    continue
                row = row.iloc[0]
                record[f"{metric}_reduction_percent"] = float(row["reduction_percent"])
                record[f"{metric}_ci95_low"] = float(row["ci95_low"])
                record[f"{metric}_ci95_high"] = float(row["ci95_high"])
            for panel_dose, tag in ((0.0, "control"), (dose, "treated")):
                row = m_sub[isclose(m_sub["dose"], panel_dose) & (m_sub["node"] == node)]
                if row.empty:
                    continue
                row = row.iloc[0]
                record[f"final_mean_{tag}"] = float(row["final_mean"])
                record[f"final_cv_{tag}"] = float(row["final_cv"])
                record[f"auc_mean_{tag}"] = float(row["auc_mean"])
                record[f"auc_cv_{tag}"] = float(row["auc_cv"])
            record.update(diag)
            records.append(record)
    table = pd.DataFrame.from_records(records)
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / "noise_scan_table.csv", index=False)

    def ci(point: object, low: object, high: object) -> str:
        if not all(np.isfinite([float(point), float(low), float(high)])):
            return "--"
        return f"{float(point):.2f} ({float(low):.2f}--{float(high):.2f})"

    def cv(value: object) -> str:
        return "--" if not np.isfinite(float(value)) else f"{float(value):.3f}"

    lines = [
        "% Generated by plot_parameter_sensitivity.py; sigma rows of the sensitivity ensemble.",
        "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrr@{}}",
        "\\toprule",
        "$\\sigma$ & Output & AUC reduction, \\% (95\\% CI) & Final reduction, \\% (95\\% CI) "
        "& CV ($C=0$) & CV ($C=" + f"{dose:g}" + "$) \\\\",
        "\\midrule",
    ]
    for (multiplier, sigma), group in table.groupby(["sigma_multiplier", "sigma"], sort=True):
        first = True
        for row in group.itertuples(index=False):
            label = f"{format_sigma(sigma)} ($\\times{format_multiplier(multiplier)}$)" if first else ""
            first = False
            lines.append(
                f"{label} & {row.node_label} & "
                f"{ci(getattr(row, 'auc_reduction_percent', np.nan), getattr(row, 'auc_ci95_low', np.nan), getattr(row, 'auc_ci95_high', np.nan))} & "
                f"{ci(getattr(row, 'final_reduction_percent', np.nan), getattr(row, 'final_ci95_low', np.nan), getattr(row, 'final_ci95_high', np.nan))} & "
                f"{cv(getattr(row, 'final_cv_control', np.nan))} & {cv(getattr(row, 'final_cv_treated', np.nan))} \\\\"
            )
    lines += ["\\bottomrule", "\\end{tabular*}"]
    (out_dir / "noise_scan_table.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir or latest_run(args.results_root)
    out_dir = args.out_dir or (run_dir / "figures")
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    nodes = [item.strip() for item in args.nodes.split(",") if item.strip()]
    main_nodes = [item.strip() for item in args.main_nodes.split(",") if item.strip()]

    effects_path = run_dir / "sensitivity_effects.csv"
    metrics_path = run_dir / "sensitivity_node_metrics.csv"
    settings_path = run_dir / "run_settings.csv"
    if not effects_path.exists() or not metrics_path.exists():
        raise SystemExit(f"Run is incomplete: expected {effects_path} and {metrics_path}")
    effects = pd.read_csv(effects_path)
    node_metrics = pd.read_csv(metrics_path)
    diagnostics_path = run_dir / "sensitivity_diagnostics.csv"
    diagnostics = pd.read_csv(diagnostics_path) if diagnostics_path.exists() else None
    final_time = read_run_setting(settings_path, "final_time", np.nan)
    sigma_base = read_run_setting(settings_path, "sigma_base", DEFAULT_SIGMA_BASE)

    nonzero_doses = sorted(value for value in effects["dose"].astype(float).unique() if value != 0.0)
    if not nonzero_doses:
        raise SystemExit("No nonzero perturbation dose available in sensitivity_effects.csv")
    dose = args.dose if args.dose is not None else nonzero_doses[-1]

    plot_main_heatmap(
        effects,
        dose,
        main_nodes,
        final_time,
        out_dir,
        formats,
        args.dpi,
        annotation_fontsize=args.annotation_fontsize,
        annotation_decimals=args.annotation_decimals,
        annotation_stroke_width=args.annotation_stroke_width,
        cmap=args.cmap,
        sigma_base=sigma_base,
    )
    plot_noise_scan(
        effects, node_metrics, dose, main_nodes, sigma_base, final_time, out_dir, formats, args.dpi
    )
    write_noise_scan_table(effects, node_metrics, diagnostics, dose, main_nodes, sigma_base, out_dir)

    if args.include_extra_figures:
        plot_reduction_panels(effects, "auc", dose, nodes, out_dir, formats, args.dpi)
        plot_reduction_panels(effects, "final", dose, nodes, out_dir, formats, args.dpi)
        plot_reduction_change_heatmap(
            effects,
            dose,
            out_dir,
            formats,
            args.dpi,
            annotation_fontsize=args.annotation_fontsize,
            annotation_decimals=args.annotation_decimals,
            annotation_stroke_width=args.annotation_stroke_width,
            cmap=args.change_cmap,
        )
        plot_noise_cv(node_metrics, dose, nodes, out_dir, formats, args.dpi)

    print(f"Sensitivity figures saved to {out_dir}")


if __name__ == "__main__":
    main()
