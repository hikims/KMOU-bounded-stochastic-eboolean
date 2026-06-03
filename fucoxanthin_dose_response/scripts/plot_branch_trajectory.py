#!/usr/bin/env python3
"""Plot branch-specific trajectories from nfkb_sim CSV output.

Examples from the project root:
  python3 scripts/plot_branch_trajectory.py results/canonical_det_mean.csv --branch canonical --show
  python3 scripts/plot_branch_trajectory.py results/noncanonical_det_mean.csv --branch noncanonical --show
  python3 scripts/plot_branch_trajectory.py results/combined_det_mean.csv --branch all --outdir results/plots

Examples from inside scripts/:
  python3 plot_branch_trajectory.py ../results/canonical_det_mean.csv --branch canonical --show
  python3 plot_branch_trajectory.py ../results/combined_det_mean.csv --branch all --outdir ../results/plots

If a matching *_sd.csv file exists, the script can draw mean +/- SD bands.
For deterministic runs, SD is usually zero.
"""

import argparse
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt


# Branch presets. Receptors are excluded by default because they often stay at 1.
# Use --include-receptors if you want to show them.
RECEPTORS = ["TNFR", "IL1R", "CD40", "LTB", "BAFFR"]

BRANCH_NODES: Dict[str, List[str]] = {
    "canonical": [
        "TRAF2_5",
        "TRAF6",
        "Btx",
        "TAK1_TAB",
        "NEMO_IKKa_IKKb",
        "IkBa_release_signal",
        "p50_p65",
        "TNFa",
        "p100",
        "Phenotype_canonical",
    ],
    "noncanonical": [
        "TRAF2_3",
        "NIK",
        "IKKa_noncanonical",
        "p100_RelB",
        "p52",
        "BAFF",
        "Phenotype_noncanonical",
    ],
    "outputs": [
        "p50_p65",
        "TNFa",
        "p100",
        "Phenotype_canonical",
        "p100_RelB",
        "p52",
        "BAFF",
        "Phenotype_noncanonical",
    ],
    "early": [
        "TRAF2_5",
        "TRAF6",
        "TRAF2_3",
        "Btx",
        "TAK1_TAB",
        "NEMO_IKKa_IKKb",
        "NIK",
        "IKKa_noncanonical",
    ],
    "receptors": RECEPTORS,
}

DEFAULT_ALL_BRANCHES = ["canonical", "noncanonical", "outputs"]


def infer_sd_path(mean_csv: Path) -> Optional[Path]:
    """Infer matching *_sd.csv path from *_mean.csv."""
    s = str(mean_csv)
    if s.endswith("_mean.csv"):
        candidate = Path(s[:-9] + "_sd.csv")
        if candidate.exists():
            return candidate
    return None


def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"File not found: {path}")
    df = pd.read_csv(path)
    if "time" not in df.columns:
        raise SystemExit(f"CSV must contain a 'time' column: {path}")
    return df


def select_nodes(branch: str, df: pd.DataFrame, include_receptors: bool) -> List[str]:
    if branch == "all_nodes":
        nodes = [c for c in df.columns if c != "time"]
    else:
        nodes = list(BRANCH_NODES[branch])
        if include_receptors and branch != "receptors":
            if branch == "canonical":
                nodes = ["TNFR", "IL1R"] + nodes
            elif branch == "noncanonical":
                nodes = ["CD40", "LTB", "BAFFR"] + nodes
            else:
                nodes = RECEPTORS + nodes

    available = [n for n in nodes if n in df.columns]
    missing = [n for n in nodes if n not in df.columns]
    if missing:
        print(f"Warning: missing columns for branch '{branch}': {missing}")
    if not available:
        raise SystemExit(f"No requested nodes are available for branch '{branch}'.")
    return available


def plot_one_branch(
    mean_df: pd.DataFrame,
    sd_df: Optional[pd.DataFrame],
    branch: str,
    out_path: Optional[Path],
    include_receptors: bool,
    title: Optional[str],
    show: bool,
    ylim: Tuple[float, float],
) -> None:
    nodes = select_nodes(branch, mean_df, include_receptors)

    plt.figure(figsize=(10, 5.6))
    t = mean_df["time"]

    for node in nodes:
        line, = plt.plot(t, mean_df[node], label=node, linewidth=1.8)
        if sd_df is not None and node in sd_df.columns:
            lower = mean_df[node] - sd_df[node]
            upper = mean_df[node] + sd_df[node]
            plt.fill_between(t, lower, upper, alpha=0.15, color=line.get_color(), linewidth=0)

    display_title = title if title else f"NF-kB {branch} branch"
    plt.title(display_title)
    plt.xlabel("time")
    plt.ylabel("normalized activity")
    plt.ylim(*ylim)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=250)
        print(f"Saved: {out_path}")

    if show:
        plt.show()
    else:
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot NF-kB branch-specific trajectories.")
    parser.add_argument("csv", help="Mean CSV file produced by nfkb_sim, e.g. results/combined_det_mean.csv")
    parser.add_argument(
        "--branch",
        default="canonical",
        choices=["canonical", "noncanonical", "outputs", "early", "receptors", "all_nodes", "all"],
        help="Branch to plot. Use 'all' to save canonical, noncanonical, and outputs plots.",
    )
    parser.add_argument("--sd", default="auto", help="SD CSV path, 'auto', or 'none'. Default: auto")
    parser.add_argument("--include-receptors", action="store_true", help="Include receptor nodes in branch plots.")
    parser.add_argument("--out", default=None, help="Output PNG path for a single branch.")
    parser.add_argument("--outdir", default=None, help="Output directory when --branch all is used.")
    parser.add_argument("--title", default=None, help="Custom title for a single branch plot.")
    parser.add_argument("--show", action="store_true", help="Show interactive plot window.")
    parser.add_argument("--ylim", nargs=2, type=float, default=[-0.02, 1.02], help="Y-axis limits. Default: -0.02 1.02")
    args = parser.parse_args()

    mean_csv = Path(args.csv)
    mean_df = load_csv(mean_csv)

    sd_df = None
    if args.sd != "none":
        if args.sd == "auto":
            sd_path = infer_sd_path(mean_csv)
        else:
            sd_path = Path(args.sd)
        if sd_path is not None and sd_path.exists():
            sd_df = load_csv(sd_path)
            print(f"Using SD file: {sd_path}")

    ylim = (float(args.ylim[0]), float(args.ylim[1]))

    if args.branch == "all":
        outdir = Path(args.outdir) if args.outdir else mean_csv.parent / "plots"
        stem = mean_csv.name.replace("_mean.csv", "").replace(".csv", "")
        for branch in DEFAULT_ALL_BRANCHES:
            out_path = outdir / f"{stem}_{branch}.png"
            plot_one_branch(
                mean_df=mean_df,
                sd_df=sd_df,
                branch=branch,
                out_path=out_path,
                include_receptors=args.include_receptors,
                title=f"{stem}: {branch}",
                show=False,
                ylim=ylim,
            )
        if args.show:
            print("Note: --branch all saves plots but does not open multiple interactive windows.")
    else:
        out_path = Path(args.out) if args.out else None
        plot_one_branch(
            mean_df=mean_df,
            sd_df=sd_df,
            branch=args.branch,
            out_path=out_path,
            include_receptors=args.include_receptors,
            title=args.title,
            show=args.show or out_path is None,
            ylim=ylim,
        )


if __name__ == "__main__":
    main()
