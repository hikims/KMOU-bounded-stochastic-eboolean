#!/usr/bin/env python3
"""Merge parallel Monte Carlo chunks produced by run_parallel_mc.sh.

The simulator writes chunk-level mean, SD, and summary files.  This script merges
chunks exactly at the level of first and second moments:

    E[X]    = weighted average of chunk means
    E[X^2]  = weighted average of (chunk variance + chunk mean^2)
    SD      = sqrt(E[X^2] - E[X]^2)

The same rule is applied to final values and AUC values in the summary files.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def read_manifest(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Manifest not found: {path}")
    df = pd.read_csv(path)
    required = {"dose", "chunk", "mc", "seed", "prefix"}
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"Manifest missing columns: {sorted(missing)}")
    if df.empty:
        raise SystemExit("Manifest is empty.")
    df["dose"] = df["dose"].astype(str)
    df["mc"] = df["mc"].astype(int)
    return df


def safe_dose_label(dose: str) -> str:
    return "".join(ch for ch in str(dose) if ch.isalnum() or ch in "_.-")


def read_timecourse(prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mean_path = Path(prefix + "_mean.csv")
    sd_path = Path(prefix + "_sd.csv")
    if not mean_path.exists():
        raise FileNotFoundError(mean_path)
    if not sd_path.exists():
        raise FileNotFoundError(sd_path)
    mean = pd.read_csv(mean_path)
    sd = pd.read_csv(sd_path)
    if "time" not in mean.columns or "time" not in sd.columns:
        raise ValueError(f"Missing time column in {mean_path} or {sd_path}")
    return mean, sd


def read_summary(prefix: str) -> pd.DataFrame:
    summary_path = Path(prefix + "_summary.csv")
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)
    summary = pd.read_csv(summary_path)
    required = {"node", "final_mean", "final_sd", "auc_mean", "auc_sd"}
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"Summary {summary_path} missing columns: {sorted(missing)}")
    return summary


def merge_timecourses(rows: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    means: List[pd.DataFrame] = []
    sds: List[pd.DataFrame] = []
    ns: List[int] = []

    for _, row in rows.iterrows():
        mean, sd = read_timecourse(str(row["prefix"]))
        means.append(mean)
        sds.append(sd)
        ns.append(int(row["mc"]))

    time = means[0]["time"].to_numpy()
    cols = [c for c in means[0].columns if c != "time"]

    for idx, mean in enumerate(means[1:], start=2):
        if not np.allclose(mean["time"].to_numpy(), time):
            raise ValueError(f"Time grid mismatch in chunk {idx}")
        if [c for c in mean.columns if c != "time"] != cols:
            raise ValueError(f"Node columns mismatch in chunk {idx}")

    total_n = float(sum(ns))
    out_mean = pd.DataFrame({"time": time})
    out_sd = pd.DataFrame({"time": time})

    for col in cols:
        combined_mean = np.zeros_like(time, dtype=float)
        combined_second = np.zeros_like(time, dtype=float)

        for n, mean, sd in zip(ns, means, sds):
            mu = mean[col].to_numpy(dtype=float)
            sig = sd[col].to_numpy(dtype=float)
            combined_mean += n * mu
            combined_second += n * (sig * sig + mu * mu)

        combined_mean /= total_n
        combined_second /= total_n
        var = np.maximum(combined_second - combined_mean * combined_mean, 0.0)

        out_mean[col] = combined_mean
        out_sd[col] = np.sqrt(var)

    return out_mean, out_sd


def merge_summary(rows: pd.DataFrame) -> pd.DataFrame:
    summaries: List[pd.DataFrame] = []
    ns: List[int] = []

    for _, row in rows.iterrows():
        summaries.append(read_summary(str(row["prefix"])))
        ns.append(int(row["mc"]))

    nodes = summaries[0]["node"].tolist()
    for idx, summary in enumerate(summaries[1:], start=2):
        if summary["node"].tolist() != nodes:
            raise ValueError(f"Summary node mismatch in chunk {idx}")

    total_n = float(sum(ns))
    records: List[Dict[str, float | str]] = []

    for node_idx, node in enumerate(nodes):
        rec: Dict[str, float | str] = {"node": node}

        for metric in ["final", "auc"]:
            mean_key = f"{metric}_mean"
            sd_key = f"{metric}_sd"
            combined_mean = 0.0
            combined_second = 0.0

            for n, summary in zip(ns, summaries):
                mu = float(summary.loc[node_idx, mean_key])
                sig = float(summary.loc[node_idx, sd_key])
                combined_mean += n * mu
                combined_second += n * (sig * sig + mu * mu)

            combined_mean /= total_n
            combined_second /= total_n
            var = max(combined_second - combined_mean * combined_mean, 0.0)
            rec[mean_key] = combined_mean
            rec[sd_key] = float(np.sqrt(var))

        final_mean = float(rec["final_mean"])
        final_sd = float(rec["final_sd"])
        rec["final_cv"] = final_sd / final_mean if abs(final_mean) > 1e-12 else np.nan
        records.append(rec)

    # Keep the original column order used by the simulator.
    return pd.DataFrame(records)[["node", "final_mean", "final_sd", "final_cv", "auc_mean", "auc_sd"]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge parallel NF-kB Monte Carlo chunks.")
    parser.add_argument("--manifest", required=True, help="CSV manifest generated by run_parallel_mc.sh")
    parser.add_argument("--out-dir", default="results", help="Output directory for merged CSV files")
    parser.add_argument(
        "--prefix-suffix",
        default="_long",
        help="Suffix for merged output prefix, e.g. fucox50_long. Default: _long",
    )
    args = parser.parse_args()

    manifest = read_manifest(Path(args.manifest))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading manifest: {args.manifest}")
    print(f"Total chunk jobs in manifest: {len(manifest)}")

    for dose, rows in manifest.groupby("dose", sort=False):
        rows = rows.sort_values("chunk")
        total_mc = int(rows["mc"].sum())
        print(f"\nMerging dose={dose} ({len(rows)} chunks, total MC={total_mc})")

        mean_df, sd_df = merge_timecourses(rows)
        summary_df = merge_summary(rows)

        dose_label = safe_dose_label(str(dose))
        out_prefix = out_dir / f"fucox{dose_label}{args.prefix_suffix}"

        mean_path = Path(str(out_prefix) + "_mean.csv")
        sd_path = Path(str(out_prefix) + "_sd.csv")
        summary_path = Path(str(out_prefix) + "_summary.csv")

        mean_df.to_csv(mean_path, index=False)
        sd_df.to_csv(sd_path, index=False)
        summary_df.to_csv(summary_path, index=False)

        print(f"  saved {mean_path}")
        print(f"  saved {sd_path}")
        print(f"  saved {summary_path}")

    print("\nMerge complete.")


if __name__ == "__main__":
    main()
