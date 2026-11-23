#!/usr/bin/env python3
"""Replay baseline Z-syndrome failures with the experimental RB-style fallback."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FAILURES = [(20283000, 10), (20283001, 1), (20283002, 6),
            (20283003, 12), (20283003, 24)]
DEFAULT_MATRIX = ROOT / "experiments/gf256_cascade/results/" \
    "g2_basis_path_n12_300frames/path/build"


def replay(seed: int, trial: int, args: argparse.Namespace) -> dict:
    stem = f"seed{seed}_trial{trial}"
    tsv = args.run_dir / f"{stem}.tsv"
    log = args.run_dir / f"{stem}.log"
    env = os.environ.copy()
    env["GF256_CASCADE_CONSTRAINT_PENALTY_EXPONENT"] = str(args.penalty)
    env["GF256_CASCADE_RB_BRANCHES"] = str(args.branches)
    env["GF256_CASCADE_RB_MAX_STEPS"] = str(args.steps)
    command = [str(args.executable), str(trial + 1), str(seed),
               str(args.p_depolar), "conditioned_robust_restart_belief",
               str(trial), "8", str(args.matrix_dir), "12"]
    with tsv.open("w", encoding="utf-8") as out, \
            log.open("w", encoding="utf-8") as err:
        subprocess.run(command, cwd=ROOT, env=env, stdout=out, stderr=err,
                       check=True)
    with tsv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1 or int(rows[0]["Total"]) != 1:
        raise RuntimeError(f"Unexpected simulator output in {tsv}")
    logs = log.read_text(encoding="utf-8")
    rb_line = next((line for line in logs.splitlines()
                    if line.startswith("[restart-belief-fallback]")), "")
    return {"seed": seed, "trial": trial,
            "logical_failures": int(rows[0]["Failures"]),
            "z_syndrome_failures": int(rows[0]["ZSyndromeFailures"]),
            "exact_recovery_frames": int(rows[0]["ExactRecoveryFrames"]),
            "exact_x_frames": int(rows[0]["ExactXFrames"]),
            "rb_log": rb_line, "tsv": str(tsv), "log": str(log)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "experiments/"
                        "gf256_cascade/results/restart_belief_known_failures")
    parser.add_argument("--executable", type=Path,
                        default=Path("/tmp/ldpc-h2-build/gf256_cascade_simulate"))
    parser.add_argument("--matrix-dir", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--p-depolar", type=float, default=0.07125)
    parser.add_argument("--penalty", type=float, default=80.0)
    parser.add_argument("--branches", type=int, default=4)
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.run_dir.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(replay, seed, trial, args)
                   for seed, trial in FAILURES if not (seed == 20283000 and trial == 10)]
        rows = [future.result() for future in concurrent.futures.as_completed(futures)]
    # Include the separately completed first replay from this same audit set.
    first_tsv = ROOT / "experiments/gf256_cascade/results/"
    first_tsv /= "restart_belief_known_fail_seed20283000_trial10.tsv"
    with first_tsv.open(newline="", encoding="utf-8") as stream:
        first = next(csv.DictReader(stream, delimiter="\t"))
    first_log = first_tsv.with_suffix(".log")
    first_text = first_log.read_text(encoding="utf-8")
    rows.append({"seed": 20283000, "trial": 10,
        "logical_failures": int(first["Failures"]),
        "z_syndrome_failures": int(first["ZSyndromeFailures"]),
        "exact_recovery_frames": int(first["ExactRecoveryFrames"]),
        "exact_x_frames": int(first["ExactXFrames"]),
        "rb_log": next((line for line in first_text.splitlines()
                        if line.startswith("[restart-belief-fallback]")), ""),
        "tsv": str(first_tsv), "log": str(first_log)})
    rows.sort(key=lambda row: (row["seed"], row["trial"]))
    summary = {"p_D": args.p_depolar, "matrix": str(args.matrix_dir),
        "branches": args.branches, "max_steps": args.steps,
        "constraint_penalty_exponent": args.penalty,
        "development_only_preselected_baseline_failures": True,
        "frames": len(rows),
        "logical_failures": sum(row["logical_failures"] for row in rows),
        "z_syndrome_failures": sum(row["z_syndrome_failures"] for row in rows),
        "exact_recovery_frames": sum(row["exact_recovery_frames"] for row in rows),
        "frames_detail": rows}
    path = args.run_dir / "summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
