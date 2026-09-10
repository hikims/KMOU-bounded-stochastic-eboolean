#!/usr/bin/env python3
"""Regression check: the baseline revision simulator matches the paper code.

The test uses random edge weights and stochastic forcing.  With all new
sensitivity multipliers set to one, the refactored executable must reproduce the
original time courses and summary columns for the same seed.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--original",
        type=Path,
        default=Path("../fucoxanthin_dose_response/nfkb_sim"),
    )
    parser.add_argument("--revision", type=Path, default=Path("./nfkb_revision_sim"))
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("--atol", type=float, default=1e-12)
    return parser.parse_args()


def run(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(command)
            + "\nstdout:\n"
            + completed.stdout
            + "\nstderr:\n"
            + completed.stderr
        )


def compare_csv(
    original: Path,
    revision: Path,
    columns: list[str] | None,
    rtol: float,
    atol: float,
) -> None:
    left = pd.read_csv(original)
    right = pd.read_csv(revision)
    if columns is None:
        columns = [column for column in left.columns if column in right.columns]
    if not columns:
        raise AssertionError(f"No common columns in {original} and {revision}")

    for column in columns:
        numeric = pd.api.types.is_numeric_dtype(left[column]) and pd.api.types.is_numeric_dtype(
            right[column]
        )
        if not numeric:
            # Text columns (node names, edge labels) must match exactly.  The
            # dtype test is used instead of ``dtype == object`` because recent
            # pandas versions store text in a dedicated string dtype.
            if left[column].astype(str).tolist() != right[column].astype(str).tolist():
                raise AssertionError(f"Text column differs: {column} ({original.name})")
        else:
            a = left[column].to_numpy(dtype=float)
            b = right[column].to_numpy(dtype=float)
            if a.shape != b.shape or not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
                difference = np.nanmax(np.abs(a - b)) if a.shape == b.shape else np.nan
                raise AssertionError(
                    f"Numeric column differs: {column} ({original.name}); max abs diff={difference}"
                )


def main() -> None:
    args = parse_args()
    original_exe = args.original.resolve()
    revision_exe = args.revision.resolve()
    if not original_exe.exists():
        raise SystemExit(f"Original executable not found: {original_exe}")
    if not revision_exe.exists():
        raise SystemExit(f"Revision executable not found: {revision_exe}")

    common = [
        "--scenario", "combined",
        "--sigma", "0.01",
        "--steps", "1000",
        "--dt", "0.001",
        "--save-every", "250",
        "--mc", "3",
        "--weights", "random",
        "--fucox", "100",
        "--seed", "12345",
    ]

    with tempfile.TemporaryDirectory(prefix="nfkb_regression_") as tmp:
        tmpdir = Path(tmp)
        original_prefix = tmpdir / "original" / "run"
        revision_prefix = tmpdir / "revision" / "run"
        run([str(original_exe), *common, "--out-prefix", str(original_prefix)])
        run(
            [
                str(revision_exe),
                *common,
                "--alpha-scale", "1",
                "--beta-scale", "1",
                "--sigma-scale", "1",
                "--rb", "0.75",
                "--write-timecourse", "1",
                "--write-replicates", "1",
                "--out-prefix", str(revision_prefix),
            ]
        )

        summary_columns = [
            "node",
            "final_mean",
            "final_sd",
            "final_cv",
            "auc_mean",
            "auc_sd",
        ]
        compare_csv(
            Path(str(original_prefix) + "_summary.csv"),
            Path(str(revision_prefix) + "_summary.csv"),
            summary_columns,
            args.rtol,
            args.atol,
        )
        for suffix in ("_mean.csv", "_sd.csv", "_rep0.csv", "_network.csv"):
            compare_csv(
                Path(str(original_prefix) + suffix),
                Path(str(revision_prefix) + suffix),
                None,
                args.rtol,
                args.atol,
            )

    print("PASS: baseline revision simulator reproduces the original stochastic run.")


if __name__ == "__main__":
    main()
