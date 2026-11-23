#!/usr/bin/env python3
"""Replay saved Z-syndrome failures with a different decoder, frame by frame.

The selected frames are sampled with their original seeds and trial indices.
Only a mode that is identical to conditioned_robust until its Z-syndrome
failure fallback may be extrapolated to the original 300-frame aggregate.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one TSV result row: {path}")
    return rows[0]


def source_failures(root: Path, basis: str) -> list[tuple[int, int]]:
    failures: list[tuple[int, int]] = []
    for log in sorted((root / basis).glob("conditioned_robust_5_*_0.07125.log")):
        match = re.fullmatch(
            r"conditioned_robust_5_(\d+)_0\.07125\.log", log.name)
        if match is None:
            continue
        seed = int(match.group(1))
        tsv = log.with_suffix(".tsv")
        row = read_row(tsv)
        text = log.read_text(encoding="utf-8")
        trials = [int(trial) for trial in re.findall(
            r"\[Z-syndrome-failure\] trial=(\d+)", text)]
        if int(row["ZSyndromeFailures"]) != len(trials) or int(row["Total"]) != 5:
            raise ValueError(f"Baseline log/TSV mismatch: {log}")
        failures.extend((seed, trial) for trial in trials)
    return sorted(failures)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", choices=("star", "path"), required=True)
    parser.add_argument("--limit", type=int, default=0,
                        help="Replay only the first N failures; 0 means all")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mode", choices=("conditioned_robust_joint_projection",
                                           "conditioned_robust_adaptive"),
                        default="conditioned_robust_joint_projection")
    parser.add_argument("--baseline-dir", type=Path, default=Path(
        "experiments/gf256_cascade/results/g2_basis_path_n12_300frames"))
    parser.add_argument("--simulator", type=Path,
                        default=Path("/tmp/gf256_cascade_build/gf256_cascade_simulate"))
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.limit < 0 or args.workers < 1 or args.workers > (os.cpu_count() or 1):
        parser.error("Invalid limit or worker count")
    failures = source_failures(args.baseline_dir, args.basis)
    expected = 38 if args.basis == "star" else 43
    if len(failures) != expected:
        parser.error(f"Expected {expected} baseline Z failures, found {len(failures)}")
    selected = failures[:args.limit] if args.limit else failures
    matrix_dir = args.baseline_dir / args.basis / "build"
    config = {"basis": args.basis, "mode": args.mode,
              "baseline_dir": str(args.baseline_dir.resolve()),
              "simulator_sha256": sha256(args.simulator),
              "hx_sha256": sha256(matrix_dir / "HX_Cascade_GF256.txt"),
              "hz_sha256": sha256(matrix_dir / "HZ_Cascade_GF256.txt"),
              "frames": [list(frame) for frame in selected],
              "p_D": 0.07125, "inner_length": 12}
    args.run_dir.mkdir(parents=True, exist_ok=True)
    config_path = args.run_dir / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        parser.error(f"Run directory has different inputs: {config_path}")
    config_path.write_text(json.dumps(config, indent=2) + "\n")

    def run_one(frame: tuple[int, int]) -> dict:
        seed, trial = frame
        prefix = args.run_dir / f"seed_{seed}_trial_{trial}"
        tsv, log = prefix.with_suffix(".tsv"), prefix.with_suffix(".log")
        if not (tsv.exists() and log.exists()):
            command = [str(args.simulator), "5", str(seed), "0.07125",
                       args.mode, str(trial), "8", str(matrix_dir), "12"]
            env = os.environ.copy()
            env["GF256_CASCADE_SKIP_TARGET_DIAG"] = "1"
            temporary_tsv = tsv.with_suffix(".tsv.tmp")
            temporary_log = log.with_suffix(".log.tmp")
            with temporary_tsv.open("w", encoding="utf-8") as output, \
                    temporary_log.open("w", encoding="utf-8") as errors:
                subprocess.run(command, stdout=output, stderr=errors,
                               env=env, check=True)
            temporary_tsv.replace(tsv)
            temporary_log.replace(log)
        row = read_row(tsv)
        text = log.read_text(encoding="utf-8")
        if (int(row["Total"]) != 1 or float(row["p_D"]) != 0.07125 or
                f"seed={seed} construction={args.mode}" not in text):
            raise ValueError(f"Invalid targeted replay: {prefix}")
        logical_fail = int(row["Failures"])
        z_fail = int(row["ZSyndromeFailures"])
        if logical_fail not in (0, 1) or z_fail not in (0, 1):
            raise ValueError(f"Invalid targeted outcome: {prefix}")
        return {"seed": seed, "trial": trial, "logical_failure": logical_fail,
                "z_syndrome_failure": z_fail,
                "undetected_logical_failure": int(row["UndetectedLogicalFailures"])}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run_one, selected))
    rescued = sum(not item["logical_failure"] for item in results)
    summary = {"basis": args.basis, "mode": args.mode,
               "baseline_z_failure_frames": len(failures),
               "replayed": len(results), "rescued": rescued,
               "still_logical_failures": len(results) - rescued,
               "still_z_syndrome_failures": sum(
                   item["z_syndrome_failure"] for item in results),
               "undetected_logical_failures_in_replay": sum(
                   item["undetected_logical_failure"] for item in results),
               "frames": results}
    if len(results) == len(failures):
        baseline = json.loads((args.baseline_dir / "aggregate_summary.json").read_text())
        baseline_failures = baseline["results"][args.basis]["logical_failures"]
        summary["derived_full_300_frame_failures"] = baseline_failures - rescued
        summary["note"] = ("Valid only because this mode exactly matches "
                           "conditioned_robust until Z-syndrome failure.")
    path = args.run_dir / "summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(path)
    print(f"Rescued {rescued}/{len(results)} selected baseline Z failures")
    if "derived_full_300_frame_failures" in summary:
        print("Full 300-frame logical failures:",
              summary["derived_full_300_frame_failures"])


if __name__ == "__main__":
    main()
