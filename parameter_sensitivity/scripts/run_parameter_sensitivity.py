#!/usr/bin/env python3
"""Run and merge the NF-kB parameter-sensitivity analysis.

The script varies one parameter family at a time while reusing the same Monte
Carlo stream seeds across parameter settings and Fucoxanthin doses.  This common-
random-number design makes the comparisons paired and reduces Monte Carlo noise.

Because the seed depends only on the chunk index, a grid can be extended
without recomputing unchanged conditions: pass ``--reuse-chunks-from RUN_ID``
to copy the completed chunk outputs of a previous run (with identical
numerical settings) into the new run directory before simulating whatever is
new.  The paired design is preserved because the copied chunks used the same
weight and Brownian streams.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROFILES = {
    "smoke": {"steps": 200, "dt": 0.01, "mc": 4, "chunks": 1, "jobs": 2, "bootstrap": 100},
    "pilot": {"steps": 50_000, "dt": 0.001, "mc": 100, "chunks": 4, "jobs": 4, "bootstrap": 1000},
    # The production profile has the same T=500 horizon as the manuscript but
    # fewer replicates per condition than the original dose-response run.  Use
    # the paired bootstrap intervals to decide whether selected cases need a
    # larger confirmation run.
    "production": {
        "steps": 500_000,
        "dt": 0.001,
        "mc": 200,
        "chunks": 8,
        "jobs": 8,
        "bootstrap": 2000,
    },
}


@dataclass(frozen=True)
class Job:
    condition_id: str
    family: str
    level: float
    alpha_scale: float
    beta_scale: float
    sigma_scale: float
    r: float
    dose: float
    chunk: int
    mc: int
    seed: int
    prefix: Path
    log: Path
    command: tuple[str, ...]


def safe_number(value: float) -> str:
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def parse_float_list(text: str) -> list[float]:
    values = [float(token.strip()) for token in text.split(",") if token.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one numeric value")
    return values


def split_mc(total: int, chunks: int) -> list[int]:
    if total <= 0 or chunks <= 0 or total < chunks:
        raise ValueError("mc must be positive and at least as large as chunks")
    base, remainder = divmod(total, chunks)
    return [base + (1 if i < remainder else 0) for i in range(chunks)]


def read_grid(path: Path, selected: set[str] | None) -> pd.DataFrame:
    grid = pd.read_csv(path)
    required = {
        "condition_id",
        "family",
        "level",
        "alpha_scale",
        "beta_scale",
        "sigma_scale",
        "r",
    }
    missing = required.difference(grid.columns)
    if missing:
        raise SystemExit(f"Grid {path} is missing columns: {sorted(missing)}")
    if grid["condition_id"].duplicated().any():
        duplicates = grid.loc[grid["condition_id"].duplicated(), "condition_id"].tolist()
        raise SystemExit(f"Duplicate condition_id values: {duplicates}")
    if "baseline" not in set(grid["condition_id"]):
        raise SystemExit("The grid must contain a condition with condition_id=baseline")
    if selected:
        selected = set(selected) | {"baseline"}
        unknown = selected.difference(set(grid["condition_id"]))
        if unknown:
            raise SystemExit(f"Unknown selected conditions: {sorted(unknown)}")
        grid = grid[grid["condition_id"].isin(selected)].copy()
    numeric = ["level", "alpha_scale", "beta_scale", "sigma_scale", "r"]
    for col in numeric:
        grid[col] = pd.to_numeric(grid[col], errors="raise")
    if not grid["r"].between(0.5, 1.0).all():
        raise SystemExit("Every r value must lie in [0.5, 1]")
    return grid.reset_index(drop=True)


def create_jobs(
    grid: pd.DataFrame,
    doses: list[float],
    exe: Path,
    run_dir: Path,
    steps: int,
    dt: float,
    total_mc: int,
    chunks: int,
    base_seed: int,
    sigma_base: float,
    scenario: str,
    progress_every: int,
) -> list[Job]:
    jobs: list[Job] = []
    chunk_sizes = split_mc(total_mc, chunks)
    for row in grid.itertuples(index=False):
        for dose in doses:
            dose_label = safe_number(dose)
            for chunk_index, chunk_mc in enumerate(chunk_sizes, start=1):
                # The seed depends only on chunk.  Consequently condition and dose
                # runs with the same chunk use identical weight and Brownian streams.
                seed = base_seed + chunk_index - 1
                prefix = (
                    run_dir
                    / "chunks"
                    / str(row.condition_id)
                    / f"dose_{dose_label}"
                    / f"chunk_{chunk_index:03d}"
                )
                log = (
                    run_dir
                    / "logs"
                    / str(row.condition_id)
                    / f"dose_{dose_label}_chunk_{chunk_index:03d}.log"
                )
                command = (
                    str(exe),
                    "--scenario",
                    scenario,
                    "--weights",
                    "random",
                    "--steps",
                    str(steps),
                    "--dt",
                    f"{dt:.17g}",
                    "--save-every",
                    str(steps),
                    "--mc",
                    str(chunk_mc),
                    "--progress-every",
                    str(progress_every),
                    "--seed",
                    str(seed),
                    "--sigma",
                    f"{sigma_base:.17g}",
                    "--alpha-scale",
                    f"{row.alpha_scale:.17g}",
                    "--beta-scale",
                    f"{row.beta_scale:.17g}",
                    "--sigma-scale",
                    f"{row.sigma_scale:.17g}",
                    "--rb",
                    f"{row.r:.17g}",
                    "--fucox",
                    f"{dose:.17g}",
                    "--write-timecourse",
                    "0",
                    "--write-replicates",
                    "1",
                    "--out-prefix",
                    str(prefix),
                )
                jobs.append(
                    Job(
                        condition_id=str(row.condition_id),
                        family=str(row.family),
                        level=float(row.level),
                        alpha_scale=float(row.alpha_scale),
                        beta_scale=float(row.beta_scale),
                        sigma_scale=float(row.sigma_scale),
                        r=float(row.r),
                        dose=float(dose),
                        chunk=chunk_index,
                        mc=chunk_mc,
                        seed=seed,
                        prefix=prefix,
                        log=log,
                        command=command,
                    )
                )
    return jobs


CHUNK_SUFFIXES = ("_summary.csv", "_replicates.csv", "_diagnostics.csv")
OPTIONAL_CHUNK_SUFFIXES = ("_network.csv",)

# Settings that must agree between two runs before chunk outputs from one run
# can be reused in the other.  They fix the numerical scheme, the horizon, the
# replicate layout, and the common random numbers (seed depends on chunk only).
REUSE_KEYS = (
    "steps",
    "dt",
    "mc_per_condition_dose",
    "chunks",
    "base_seed",
    "sigma_base",
    "scenario",
    "doses",
)


def output_complete(job: Job) -> bool:
    return all(Path(str(job.prefix) + suffix).exists() for suffix in CHUNK_SUFFIXES)


def check_reuse_compatibility(source_dir: Path, settings: dict[str, object]) -> None:
    """Abort unless the source run used the same numerical settings."""

    settings_path = source_dir / "run_settings.csv"
    if not settings_path.exists():
        raise SystemExit(f"Cannot reuse chunks: {settings_path} not found")
    source = pd.read_csv(settings_path).iloc[0].to_dict()
    mismatches = []
    for key in REUSE_KEYS:
        if key not in source:
            mismatches.append(f"{key}: missing in source run")
            continue
        current = settings[key]
        previous = source[key]
        if isinstance(current, str):
            same = str(previous) == current
        else:
            same = math.isclose(float(previous), float(current), rel_tol=0.0, abs_tol=1e-12)
        if not same:
            mismatches.append(f"{key}: source={previous!r}, current={current!r}")
    if mismatches:
        raise SystemExit(
            "Cannot reuse chunks from "
            f"{source_dir}: numerical settings differ.\n  " + "\n  ".join(mismatches)
        )


def reuse_chunks(jobs: Iterable[Job], source_dir: Path, run_dir: Path) -> int:
    """Copy completed chunk outputs of matching jobs from a previous run.

    A job matches when the same condition_id, dose, and chunk exist under
    ``source_dir/chunks``.  Jobs whose outputs are already present in the new
    run directory are left untouched, and jobs with no source outputs are
    simply simulated afresh by the normal job runner.
    """

    copied = 0
    for job in jobs:
        if output_complete(job):
            continue
        relative = job.prefix.relative_to(run_dir)
        source_prefix = source_dir / relative
        if not all(Path(str(source_prefix) + suffix).exists() for suffix in CHUNK_SUFFIXES):
            continue
        job.prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in CHUNK_SUFFIXES + OPTIONAL_CHUNK_SUFFIXES:
            source_file = Path(str(source_prefix) + suffix)
            if source_file.exists():
                shutil.copy2(source_file, Path(str(job.prefix) + suffix))
        copied += 1
    return copied


def run_job(job: Job, skip_existing: bool) -> tuple[Job, str]:
    if skip_existing and output_complete(job):
        return job, "skipped"
    job.prefix.parent.mkdir(parents=True, exist_ok=True)
    job.log.parent.mkdir(parents=True, exist_ok=True)
    with job.log.open("w", encoding="utf-8") as handle:
        handle.write("command: " + shlex.join(job.command) + "\n\n")
        handle.flush()
        completed = subprocess.run(
            job.command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        tail = job.log.read_text(encoding="utf-8", errors="replace")[-4000:]
        raise RuntimeError(
            f"Job failed ({job.condition_id}, dose={job.dose}, chunk={job.chunk}).\n"
            f"Log: {job.log}\n{tail}"
        )
    return job, "completed"


def write_manifest(jobs: Iterable[Job], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "condition_id",
                "family",
                "level",
                "alpha_scale",
                "beta_scale",
                "sigma_scale",
                "r",
                "dose",
                "chunk",
                "mc",
                "seed",
                "prefix",
                "log",
                "command",
            ]
        )
        for job in jobs:
            writer.writerow(
                [
                    job.condition_id,
                    job.family,
                    job.level,
                    job.alpha_scale,
                    job.beta_scale,
                    job.sigma_scale,
                    job.r,
                    job.dose,
                    job.chunk,
                    job.mc,
                    job.seed,
                    job.prefix,
                    job.log,
                    shlex.join(job.command),
                ]
            )


def merge_summary(group: pd.DataFrame) -> pd.DataFrame:
    summaries: list[pd.DataFrame] = []
    counts: list[int] = []
    for row in group.itertuples(index=False):
        summaries.append(pd.read_csv(str(row.prefix) + "_summary.csv"))
        counts.append(int(row.mc))
    nodes = summaries[0]["node"].tolist()
    if any(df["node"].tolist() != nodes for df in summaries[1:]):
        raise RuntimeError("Node order differs between summary chunks")
    total = float(sum(counts))
    records: list[dict[str, float | str]] = []
    for node_idx, node in enumerate(nodes):
        record: dict[str, float | str] = {"node": node}
        for prefix in ("final", "auc"):
            mean_col = f"{prefix}_mean"
            sd_col = f"{prefix}_sd"
            mean = sum(n * float(df.loc[node_idx, mean_col]) for n, df in zip(counts, summaries)) / total
            second = sum(
                n
                * (
                    float(df.loc[node_idx, sd_col]) ** 2
                    + float(df.loc[node_idx, mean_col]) ** 2
                )
                for n, df in zip(counts, summaries)
            ) / total
            sd = math.sqrt(max(second - mean * mean, 0.0))
            record[mean_col] = mean
            record[sd_col] = sd
            record[f"{prefix}_se"] = sd / math.sqrt(total)
            record[f"{prefix}_cv"] = sd / abs(mean) if abs(mean) > 1e-12 else np.nan
        records.append(record)
    return pd.DataFrame.from_records(records)


def merge_replicates(group: pd.DataFrame) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for row in group.itertuples(index=False):
        frame = pd.read_csv(str(row.prefix) + "_replicates.csv")
        frame["chunk"] = int(row.chunk)
        frame["pair_key"] = (
            "c"
            + frame["chunk"].astype(str).str.zfill(3)
            + "_s"
            + frame["stream_seed"].astype(str)
        )
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    if result.duplicated(["pair_key", "node"]).any():
        raise RuntimeError("Duplicate pair_key/node rows after merging replicate chunks")
    return result


def merge_diagnostics(group: pd.DataFrame) -> pd.DataFrame:
    frames = [pd.read_csv(str(row.prefix) + "_diagnostics.csv") for row in group.itertuples(index=False)]
    data = pd.concat(frames, ignore_index=True)
    values = int(data["predictor_values"].sum())
    violations = int(data["predictor_violations"].sum())
    total_overshoot = float((data["mean_overshoot_per_value"] * data["predictor_values"]).sum())
    first = data.iloc[0].to_dict()
    first.update(
        {
            "replicates": int(data["replicates"].sum()),
            "predictor_violations": violations,
            "predictor_values": values,
            "violation_rate": violations / values if values else 0.0,
            "max_lower_overshoot": float(data["max_lower_overshoot"].max()),
            "max_upper_overshoot": float(data["max_upper_overshoot"].max()),
            "mean_overshoot_per_value": total_overshoot / values if values else 0.0,
        }
    )
    return pd.DataFrame([first])


def paired_bootstrap_reduction(
    control: np.ndarray,
    treated: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    if control.size == 0 or treated.size != control.size:
        return np.nan, np.nan, np.nan
    control_mean = float(control.mean())
    if not np.isfinite(control_mean) or abs(control_mean) <= 1e-12:
        # A percentage reduction is undefined when the matched control is zero.
        return np.nan, np.nan, np.nan
    point = 100.0 * (1.0 - float(treated.mean()) / control_mean)
    if n_bootstrap <= 0:
        return point, np.nan, np.nan
    samples = np.full(n_bootstrap, np.nan, dtype=float)
    n = control.size
    block = 256
    for start in range(0, n_bootstrap, block):
        stop = min(start + block, n_bootstrap)
        indices = rng.integers(0, n, size=(stop - start, n))
        cmean = control[indices].mean(axis=1)
        tmean = treated[indices].mean(axis=1)
        valid = np.isfinite(cmean) & np.isfinite(tmean) & (np.abs(cmean) > 1e-12)
        block_samples = np.full(stop - start, np.nan, dtype=float)
        block_samples[valid] = 100.0 * (1.0 - tmean[valid] / cmean[valid])
        samples[start:stop] = block_samples
    finite_samples = samples[np.isfinite(samples)]
    if finite_samples.size == 0:
        return point, np.nan, np.nan
    low, high = np.percentile(finite_samples, [2.5, 97.5])
    return point, float(low), float(high)


def merge_all(
    manifest: pd.DataFrame,
    grid: pd.DataFrame,
    run_dir: Path,
    n_bootstrap: int,
    bootstrap_seed: int,
) -> None:
    merged_root = run_dir / "merged"
    merged_root.mkdir(parents=True, exist_ok=True)
    summary_frames: list[pd.DataFrame] = []
    replicate_frames: list[pd.DataFrame] = []
    diagnostic_frames: list[pd.DataFrame] = []

    metadata = grid.set_index("condition_id")
    for (condition_id, dose), group in manifest.groupby(["condition_id", "dose"], sort=False):
        summary = merge_summary(group)
        replicates = merge_replicates(group)
        diagnostics = merge_diagnostics(group)
        dose_label = safe_number(float(dose))
        out_dir = merged_root / condition_id
        out_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_dir / f"dose_{dose_label}_summary.csv", index=False)
        replicates.to_csv(out_dir / f"dose_{dose_label}_replicates.csv", index=False)
        diagnostics.to_csv(out_dir / f"dose_{dose_label}_diagnostics.csv", index=False)

        row = metadata.loc[condition_id]
        common_metadata = {
            "condition_id": condition_id,
            "family": str(row["family"]),
            "level": float(row["level"]),
            "alpha_scale": float(row["alpha_scale"]),
            "beta_scale": float(row["beta_scale"]),
            "sigma_scale": float(row["sigma_scale"]),
            "r": float(row["r"]),
            "dose": float(dose),
        }
        # Diagnostics already contain some parameter columns. Assignment, rather
        # than DataFrame.insert(), avoids duplicate-column failures.
        def attach_metadata(frame: pd.DataFrame) -> pd.DataFrame:
            frame = frame.copy()
            for key, value in common_metadata.items():
                frame[key] = value
            metadata_columns = list(common_metadata)
            remaining_columns = [col for col in frame.columns if col not in metadata_columns]
            return frame[metadata_columns + remaining_columns]

        summary = attach_metadata(summary)
        replicates = attach_metadata(replicates)
        diagnostics = attach_metadata(diagnostics)
        summary_frames.append(summary)
        replicate_frames.append(replicates)
        diagnostic_frames.append(diagnostics)

    node_metrics = pd.concat(summary_frames, ignore_index=True)
    replicate_metrics = pd.concat(replicate_frames, ignore_index=True)
    diagnostics = pd.concat(diagnostic_frames, ignore_index=True)
    node_metrics.to_csv(run_dir / "sensitivity_node_metrics.csv", index=False)
    replicate_metrics.to_csv(run_dir / "sensitivity_replicates.csv", index=False)
    diagnostics.to_csv(run_dir / "sensitivity_diagnostics.csv", index=False)

    doses = sorted(float(value) for value in replicate_metrics["dose"].unique())
    if 0.0 not in doses:
        raise RuntimeError("A zero-dose control is required to calculate suppression")

    rng = np.random.default_rng(bootstrap_seed)
    effect_records: list[dict[str, float | int | str]] = []
    for condition_id, condition_data in replicate_metrics.groupby("condition_id", sort=False):
        control = condition_data[condition_data["dose"] == 0.0]
        row = metadata.loc[condition_id]
        for dose in [value for value in doses if value != 0.0]:
            treated = condition_data[condition_data["dose"] == dose]
            for node in sorted(set(control["node"]).intersection(treated["node"])):
                cnode = control[control["node"] == node]
                tnode = treated[treated["node"] == node]
                paired = cnode.merge(
                    tnode,
                    on=["pair_key", "node"],
                    suffixes=("_control", "_treated"),
                    validate="one_to_one",
                )
                for metric in ("auc", "final_value"):
                    point, low, high = paired_bootstrap_reduction(
                        paired[f"{metric}_control"].to_numpy(dtype=float),
                        paired[f"{metric}_treated"].to_numpy(dtype=float),
                        n_bootstrap,
                        rng,
                    )
                    effect_records.append(
                        {
                            "condition_id": condition_id,
                            "family": row["family"],
                            "level": row["level"],
                            "alpha_scale": row["alpha_scale"],
                            "beta_scale": row["beta_scale"],
                            "sigma_scale": row["sigma_scale"],
                            "r": row["r"],
                            "dose": dose,
                            "node": node,
                            "metric": "auc" if metric == "auc" else "final",
                            "reduction_percent": point,
                            "ci95_low": low,
                            "ci95_high": high,
                            "paired_replicates": len(paired),
                        }
                    )
    effects = pd.DataFrame.from_records(effect_records)
    effects.to_csv(run_dir / "sensitivity_effects.csv", index=False)

    # Relative change of each raw metric from the baseline parameter condition.
    baseline = node_metrics[node_metrics["condition_id"] == "baseline"].set_index(["dose", "node"])
    deviation_records: list[dict[str, float | str]] = []
    for row in node_metrics.itertuples(index=False):
        base = baseline.loc[(row.dose, row.node)]
        for metric in ("final_mean", "final_cv", "auc_mean", "auc_cv"):
            base_value = float(base[metric])
            value = float(getattr(row, metric))
            relative = 100.0 * (value / base_value - 1.0) if abs(base_value) > 1e-12 else np.nan
            deviation_records.append(
                {
                    "condition_id": row.condition_id,
                    "family": row.family,
                    "level": row.level,
                    "dose": row.dose,
                    "node": row.node,
                    "metric": metric,
                    "value": value,
                    "baseline_value": base_value,
                    "relative_change_percent": relative,
                }
            )
    pd.DataFrame.from_records(deviation_records).to_csv(
        run_dir / "sensitivity_relative_to_baseline.csv", index=False
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=PROFILES, default="smoke")
    parser.add_argument("--grid", type=Path, default=Path("config/sensitivity_grid.csv"))
    parser.add_argument("--exe", type=Path, default=Path("./nfkb_revision_sim"))
    parser.add_argument("--results-root", type=Path, default=Path("results/parameter_sensitivity"))
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--doses", type=parse_float_list, default=parse_float_list("0,100"))
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--mc", type=int, default=None, help="total replicates per condition and dose")
    parser.add_argument("--chunks", type=int, default=None)
    parser.add_argument("--jobs", type=int, default=None)
    parser.add_argument("--bootstrap", type=int, default=None)
    parser.add_argument("--base-seed", type=int, default=41000)
    parser.add_argument("--bootstrap-seed", type=int, default=99173)
    parser.add_argument("--sigma-base", type=float, default=0.01)
    parser.add_argument("--scenario", default="combined")
    parser.add_argument("--progress-every", type=int, default=0)
    parser.add_argument("--conditions", default=None, help="comma-separated condition_id subset")
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--reuse-chunks-from",
        default=None,
        metavar="RUN_ID",
        help=(
            "copy completed chunk outputs from a previous run under --results-root "
            "(same numerical settings required); only conditions absent from that "
            "run are simulated, then everything is merged"
        ),
    )
    parser.add_argument("--run-only", action="store_true", help="run jobs but do not merge")
    parser.add_argument("--merge-only", action="store_true", help="merge an existing run without launching jobs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    profile = PROFILES[args.profile]
    steps = args.steps if args.steps is not None else profile["steps"]
    dt = args.dt if args.dt is not None else profile["dt"]
    total_mc = args.mc if args.mc is not None else profile["mc"]
    chunks = args.chunks if args.chunks is not None else profile["chunks"]
    max_jobs = args.jobs if args.jobs is not None else profile["jobs"]
    n_bootstrap = args.bootstrap if args.bootstrap is not None else profile["bootstrap"]

    if steps <= 0 or dt <= 0.0 or total_mc <= 0 or chunks <= 0 or max_jobs <= 0:
        raise SystemExit("steps, dt, mc, chunks, and jobs must be positive")
    if total_mc < chunks:
        raise SystemExit("mc must be at least as large as chunks")
    if 0.0 not in args.doses:
        raise SystemExit("--doses must include 0 for the matched control")

    selected = None
    if args.conditions:
        selected = {item.strip() for item in args.conditions.split(",") if item.strip()}
    grid = read_grid(args.grid, selected)

    run_id = args.run_id or f"{args.profile}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_dir = args.results_root / run_id
    manifest_path = run_dir / "job_manifest.csv"

    exe = args.exe.resolve()
    if not args.merge_only and not exe.exists():
        raise SystemExit(f"Executable not found: {exe}. Run `make` first.")

    jobs = create_jobs(
        grid=grid,
        doses=args.doses,
        exe=exe,
        run_dir=run_dir,
        steps=steps,
        dt=dt,
        total_mc=total_mc,
        chunks=chunks,
        base_seed=args.base_seed,
        sigma_base=args.sigma_base,
        scenario=args.scenario,
        progress_every=args.progress_every,
    )
    write_manifest(jobs, manifest_path)
    grid.to_csv(run_dir / "conditions_used.csv", index=False)

    settings = {
        "profile": args.profile,
        "steps": steps,
        "dt": dt,
        "final_time": steps * dt,
        "mc_per_condition_dose": total_mc,
        "chunks": chunks,
        "jobs": max_jobs,
        "bootstrap": n_bootstrap,
        "base_seed": args.base_seed,
        "sigma_base": args.sigma_base,
        "scenario": args.scenario,
        "doses": ";".join(f"{value:g}" for value in args.doses),
    }
    if args.reuse_chunks_from:
        source_dir = args.results_root / args.reuse_chunks_from
        if source_dir.resolve() == run_dir.resolve():
            raise SystemExit("--reuse-chunks-from must name a different run than --run-id")
        check_reuse_compatibility(source_dir, settings)
        settings["reused_chunks_from"] = args.reuse_chunks_from
    pd.DataFrame([settings]).to_csv(run_dir / "run_settings.csv", index=False)

    if args.reuse_chunks_from and not args.merge_only:
        copied = reuse_chunks(jobs, args.results_root / args.reuse_chunks_from, run_dir)
        print(f"Reused {copied} completed chunk outputs from {args.reuse_chunks_from}")

    if not args.merge_only:
        print(f"Run directory: {run_dir}")
        print(f"Jobs: {len(jobs)}; parallel processes: {max_jobs}")
        completed = skipped = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_jobs) as executor:
            futures = {executor.submit(run_job, job, args.skip_existing): job for job in jobs}
            for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                job, status = future.result()
                completed += status == "completed"
                skipped += status == "skipped"
                print(
                    f"[{index}/{len(jobs)}] {status}: {job.condition_id}, "
                    f"dose={job.dose:g}, chunk={job.chunk}"
                )
        print(f"Completed jobs: {completed}; skipped jobs: {skipped}")

    if args.run_only:
        return

    manifest = pd.read_csv(manifest_path)
    missing = [
        prefix
        for prefix in manifest["prefix"]
        if not Path(str(prefix) + "_summary.csv").exists()
    ]
    if missing:
        raise SystemExit(f"Cannot merge: {len(missing)} job outputs are missing. First: {missing[0]}")

    print("Merging chunks and calculating paired bootstrap intervals...")
    merge_all(manifest, grid, run_dir, n_bootstrap, args.bootstrap_seed)
    print(f"Sensitivity analysis complete: {run_dir}")


if __name__ == "__main__":
    main()
