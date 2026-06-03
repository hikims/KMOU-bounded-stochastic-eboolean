#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# Paths
# ============================================================

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent if HERE.parent.name == "scripts" else Path.cwd().resolve()

RESULTS = ROOT / "results"
FIGURE = ROOT / "figure"
OUTDIR = FIGURE / "theory_benchmarks"
PREFIX = RESULTS / "theory"


# ============================================================
# Plot style
# ============================================================

plt.rcParams.update(
    {
        "font.size": 15,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "legend.fontsize": 13,
        "figure.titlesize": 22,
        "lines.linewidth": 2.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


# ============================================================
# Utilities
# ============================================================

def save_png_eps(fig: plt.Figure, path_no_ext: Path) -> None:
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)

    png = path_no_ext.with_suffix(".png")
    eps = path_no_ext.with_suffix(".eps")

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(eps, format="eps", bbox_inches="tight")

    plt.close(fig)

    print("saved", png)
    print("saved", eps)


def style_axis(ax, *, grid_axis: str = "both") -> None:
    if grid_axis in ("both", "x"):
        ax.grid(
            True,
            axis="x",
            which="major",
            linestyle="--",
            color="#d9d9d9",
            linewidth=1.0,
        )

    if grid_axis in ("both", "y"):
        ax.grid(
            True,
            axis="y",
            which="major",
            linestyle="--",
            color="#d9d9d9",
            linewidth=1.0,
        )

    ax.minorticks_on()

    if grid_axis == "both":
        ax.grid(
            True,
            which="minor",
            linestyle=":",
            color="#eeeeee",
            linewidth=0.8,
        )

    ax.tick_params(axis="both", which="major", labelsize=15, length=5, width=1.1)


# ============================================================
# Benchmark I
# ============================================================

def plot_boundary() -> None:
    ts_path = PREFIX.with_name(PREFIX.name + "_boundary_timeseries.csv")
    sm_path = PREFIX.with_name(PREFIX.name + "_boundary_summary.csv")

    if not ts_path.exists() or not sm_path.exists():
        print("skip boundary: missing CSV")
        return

    ts = pd.read_csv(ts_path)
    sm = pd.read_csv(sm_path)

    model_order = [
        "Original EM",
        "Tanh EM",
        "Revised projected EM",
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.9))

    # A. maximum activity over time
    ax = axes[0]

    for model in model_order:
        sub = ts[ts["model"] == model]
        if sub.empty:
            continue
        ax.plot(sub["time"], sub["max"], label=model)

    ax.axhline(1.0, linestyle="--", linewidth=1.5, color="black")
    ax.set_title("Maximum activity over ensemble", pad=10)
    ax.set_xlabel("time")
    ax.set_ylabel("max activity")

    style_axis(ax)
    ax.legend(frameon=False, fontsize=12, loc="upper left")

    # B. maximum boundary violation
    ax = axes[1]

    sm = sm.set_index("model").loc[model_order].reset_index()
    x = np.arange(len(sm))
    vals = sm["Vmax"].to_numpy(float)

    ax.bar(x, vals)
    ax.set_xticks(x)
    ax.set_xticklabels(model_order, rotation=25, ha="right")
    ax.set_title("Maximum boundary violation", pad=10)
    ax.set_ylabel(r"$V_{\max}$")

    ymax = max(float(np.nanmax(vals)) * 1.22, 0.1)
    ax.set_ylim(0, ymax)

    style_axis(ax, grid_axis="y")

    for i, v in enumerate(vals):
        ax.text(
            i,
            v + ymax * 0.025,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
        )

    # C. exit fraction
    ax = axes[2]

    vals = sm["exit_fraction"].to_numpy(float)

    ax.bar(x, vals)
    ax.set_xticks(x)
    ax.set_xticklabels(model_order, rotation=25, ha="right")
    ax.set_title("Trajectories leaving [0,1]", pad=10)
    ax.set_ylabel("fraction")
    ax.set_ylim(0, 1.15)

    style_axis(ax, grid_axis="y")

    for i, v in enumerate(vals):
        ax.text(
            i,
            v + 0.035,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
        )

    fig.suptitle(
        "Benchmark I: boundary preservation",
        y=1.04,
        fontweight="normal",
    )

    fig.tight_layout()

    save_png_eps(fig, OUTDIR / "boundary_preservation_benchmark")


# ============================================================
# Benchmark II
# ============================================================

def plot_moment() -> None:
    ts_path = PREFIX.with_name(PREFIX.name + "_moment_timeseries.csv")
    sm_path = PREFIX.with_name(PREFIX.name + "_moment_summary.csv")

    if not ts_path.exists() or not sm_path.exists():
        print("skip moment: missing CSV")
        return

    df = pd.read_csv(ts_path)
    sm = pd.read_csv(sm_path).iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.2))

    # A. mean dynamics
    ax = axes[0]

    ax.plot(df["time"], df["mc_m1"], label="MC mean node 1")
    ax.plot(df["time"], df["mc_m2"], label="MC mean node 2")
    ax.plot(df["time"], df["ode_m1"], "--", label="moment ODE node 1")
    ax.plot(df["time"], df["ode_m2"], "--", label="moment ODE node 2")

    ax.set_title("Mean dynamics", pad=10)
    ax.set_xlabel("time")
    ax.set_ylabel("mean activity")

    style_axis(ax)
    ax.legend(frameon=False, fontsize=11, loc="lower right")

    # B. CV amplification
    ax = axes[1]

    ax.plot(df["time"], df["mc_cv1"], label="MC CV node 1")
    ax.plot(df["time"], df["mc_cv2"], label="MC CV node 2")
    ax.plot(df["time"], df["ode_cv1"], "--", label="moment ODE CV node 1")
    ax.plot(df["time"], df["ode_cv2"], "--", label="moment ODE CV node 2")

    ax.axhline(
        float(sm["theory_CV1"]),
        color="black",
        linestyle=":",
        linewidth=1.5,
    )

    ax.axhline(
        float(sm["theory_CV2"]),
        color="black",
        linestyle=":",
        linewidth=1.5,
    )

    ax.text(
        0.03,
        0.96,
        rf"Theory $\mathrm{{CV}}_2/\mathrm{{CV}}_1={float(sm['theory_CV2_over_CV1']):.2f}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        bbox=dict(facecolor="white", edgecolor="none", alpha=1.0),
    )

    ax.set_title("Downstream CV amplification", pad=10)
    ax.set_xlabel("time")
    ax.set_ylabel("coefficient of variation")

    style_axis(ax)
    ax.legend(frameon=False, fontsize=10, loc="lower right")

    fig.suptitle(
        "Benchmark II: moment equations and downstream CV amplification",
        y=1.04,
        fontweight="normal",
    )

    fig.tight_layout()

    save_png_eps(fig, OUTDIR / "moment_cv_amplification")


# ============================================================
# Benchmark III
# ============================================================

def plot_hopf() -> None:
    amp_path = PREFIX.with_name(PREFIX.name + "_hopf_amplitude.csv")
    ts_path = PREFIX.with_name(PREFIX.name + "_hopf_timeseries.csv")

    if not amp_path.exists() or not ts_path.exists():
        print("skip hopf: missing CSV")
        return

    amp = pd.read_csv(amp_path)
    ts = pd.read_csv(ts_path)

    Lc = 8.0

    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.4))

    # ------------------------------------------------------------
    # A. representative time series
    # ------------------------------------------------------------
    ax = axes[0]

    for L in sorted(ts["L"].unique()):
        sub = ts[ts["L"] == L]
        ax.plot(
            sub["time_after_transient"],
            sub["X1"],
            linewidth=2.8,
            label=rf"$L = {L:.1f}$",
        )

    ax.axhline(
        0.5,
        linestyle=":",
        color="black",
        linewidth=1.3,
    )

    ax.set_title("Representative motif activity", pad=12)
    ax.set_xlabel("time after transient")
    ax.set_ylabel(r"$X_1 = 0.5 + z_1$")

    # y축 상한 0.7
    ax.set_ylim(0.38, 0.70)

    # y축 tick은 0.4, 0.5, 0.6만 표시
    ax.set_yticks([0.4, 0.5, 0.6])

    style_axis(ax)

    # 레전드를 그림 내부 상단에 배치
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=3,
        frameon=True,
        framealpha=0.95,
        facecolor="white",
        edgecolor="none",
        fontsize=13,
        handlelength=2.6,
        columnspacing=1.5,
    )

    # ------------------------------------------------------------
    # B. amplitude scan
    # ------------------------------------------------------------
    ax = axes[1]

    L = amp["L"].to_numpy(float)
    A = amp["amplitude_z1"].to_numpy(float)

    ax.plot(
        L,
        A,
        marker="o",
        linewidth=2.5,
        label="measured amplitude",
    )

    lc_line = ax.axvline(
        Lc,
        linestyle="--",
        color="black",
        linewidth=1.5,
        label=rf"theory $L_c={Lc:.1f}$",
    )

    post = L > Lc
    xfit = np.sqrt(np.maximum(L[post] - Lc, 0.0))
    yfit = A[post]

    sqrt_line = None

    if len(xfit) >= 2 and np.dot(xfit, xfit) > 0:
        coef = float(np.dot(xfit, yfit) / np.dot(xfit, xfit))

        sqrt_line, = ax.plot(
            L[post],
            coef * np.sqrt(L[post] - Lc),
            "--",
            linewidth=2.5,
            label=r"$\propto \sqrt{L-L_c}$ guide",
        )

    ax.set_title("Oscillation amplitude scan", pad=12)
    ax.set_xlabel("loop gain L")
    ax.set_ylabel(r"amplitude of $z_1$")

    style_axis(ax)

    right_handles = [lc_line]
    right_labels = [rf"theory $L_c={Lc:.1f}$"]

    if sqrt_line is not None:
        right_handles.append(sqrt_line)
        right_labels.append(r"$\propto \sqrt{L-L_c}$ guide")

    ax.legend(
        right_handles,
        right_labels,
        loc="upper left",
        frameon=True,
        framealpha=0.90,
        facecolor="white",
        edgecolor="none",
        fontsize=13,
        handlelength=2.6,
    )

    # 메인 타이틀을 더 아래로
    fig.suptitle(
        "Benchmark III: Hopf-type threshold in a negative-feedback motif",
        y=0.955,
        fontweight="normal",
    )

    # suptitle과 subplot 사이 간격 조정
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.90])

    save_png_eps(fig, OUTDIR / "hopf_threshold_scan")


# ============================================================
# Index file
# ============================================================

def write_index() -> None:
    files = sorted(
        list(OUTDIR.rglob("*.png")) +
        list(OUTDIR.rglob("*.eps"))
    )

    path = OUTDIR / "theory_benchmarks_index.txt"

    with open(path, "w", encoding="utf-8") as f:
        for p in files:
            f.write(str(p.relative_to(OUTDIR)) + "\n")

    print("saved", path)


# ============================================================
# Main
# ============================================================

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    plot_boundary()
    plot_moment()
    plot_hopf()
    write_index()

    print("Done. Theory benchmark figures saved in", OUTDIR)


if __name__ == "__main__":
    main()