#!/usr/bin/env python3
"""Plot selected node trajectories from a simulator CSV file.

Example:
  python3 scripts/plot_trajectory.py results/canonical_det_mean.csv \
      --nodes TNFR TRAF2_5 TAK1_TAB NEMO_IKKa_IKKb p50_p65 Phenotype_canonical \
      --out results/canonical_det_mean.png
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("csv", help="CSV produced by nfkb_sim, e.g. *_mean.csv")
parser.add_argument("--nodes", nargs="*", default=None, help="Node names to plot. Default: all nodes except receptors.")
parser.add_argument("--out", default=None, help="Output image path. If omitted, show interactively.")
args = parser.parse_args()

df = pd.read_csv(args.csv)
if "time" not in df.columns:
    raise SystemExit("CSV must contain a time column")

nodes = args.nodes
if not nodes:
    skip = {"TNFR", "IL1R", "CD40", "LTB", "BAFFR"}
    nodes = [c for c in df.columns if c != "time" and c not in skip]

missing = [n for n in nodes if n not in df.columns]
if missing:
    raise SystemExit(f"Missing node columns: {missing}")

plt.figure(figsize=(9, 5))
for n in nodes:
    plt.plot(df["time"], df[n], label=n)
plt.xlabel("time")
plt.ylabel("normalized activity")
plt.ylim(-0.02, 1.02)
plt.legend(loc="best", fontsize=8)
plt.tight_layout()
if args.out:
    plt.savefig(args.out, dpi=200)
else:
    plt.show()
