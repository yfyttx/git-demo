#!/usr/bin/env python3
"""Measure whether cascade projection disagreement predicts logical failures."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MATRIX_DIR = ROOT / "experiments/gf256_cascade/results/"
MATRIX_DIR /= "g2_basis_path_n12_300frames/path/build"
EXECUTABLE = Path("/tmp/ldpc-h2-build/gf256_cascade_simulate")
RISK_PATTERN = re.compile(r"\[projection-risk-audit\] trial=(\d+) (.*)")


def parse_fields(text: str) -> dict[str, int]:
    return {key: int(value) for key, value in
            re.findall(r"(\w+)=(-?\d+)", text)}


def run_batch(args: tuple[int, int, Path]) -> dict:
    seed, frames, run_dir = args
    log_path = run_dir / f"seed{seed}.log"
    tsv_path = run_dir / f"seed{seed}.tsv"
    env = os.environ.copy()
    env["GF256_CASCADE_PROJECTION_RISK_AUDIT"] = "1"
    start = time.perf_counter()
    command = [str(EXECUTABLE), str(frames), str(seed), "0.07125",
               "conditioned_robust_adaptive", "-1", "8", str(MATRIX_DIR), "12"]
    with tsv_path.open("w", encoding="utf-8") as stdout, \
            log_path.open("w", encoding="utf-8") as stderr:
        subprocess.run(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr,
                       check=True)
    elapsed = time.perf_counter() - start
    records = []
    for line in log_path.read_text(encoding="utf-8").splitlines():
        match = RISK_PATTERN.search(line)
        if match:
            fields = parse_fields(match.group(2))
            records.append({"seed": seed, "trial": int(match.group(1)), **fields})
    if len(records) != frames:
        raise RuntimeError(f"Expected {frames} risk rows in {log_path}, got {len(records)}")
    return {"seed": seed, "frames": frames, "seconds": elapsed,
            "log": str(log_path), "tsv": str(tsv_path), "records": records}


def auc_positive_higher(records: list[dict], score_name: str) -> float | None:
    positives = [r[score_name] for r in records if r["logical_failure"]]
    negatives = [r[score_name] for r in records if not r["logical_failure"]]
    if not positives or not negatives:
        return None
    wins = sum(1.0 if p > n else 0.5 if p == n else 0.0
               for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=100)
    parser.add_argument("--seed-start", type=int, default=20283000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--run-dir", type=Path, default=ROOT /
                        "experiments/gf256_cascade/results/projection_risk_audit/holdout100")
    args = parser.parse_args()
    if args.frames < 1 or args.workers < 1:
        parser.error("frames and workers must be positive")
    if not EXECUTABLE.is_file():
        parser.error(f"simulator not found: {EXECUTABLE}")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    batches = min(args.workers, args.frames)
    base, remainder = divmod(args.frames, batches)
    jobs = []
    first_seed = args.seed_start
    for index in range(batches):
        size = base + (index < remainder)
        jobs.append((first_seed + index, size, args.run_dir))

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=batches) as pool:
        futures = [pool.submit(run_batch, job) for job in jobs]
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            results.append(result)
            failures = sum(r["logical_failure"] for r in result["records"])
            print(f"[{done}/{batches}] seed={result['seed']} "
                  f"fail={failures}/{result['frames']} "
                  f"elapsed={result['seconds']:.1f}s", flush=True)

    records = [r for batch in results for r in batch["records"]]
    records.sort(key=lambda r: (r["seed"], r["trial"]))
    csv_path = args.run_dir / "projection_risk_frames.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    total = len(records)
    failures = sum(r["logical_failure"] for r in records)
    risk_key = "risk_hard_disagreements"
    ranked = sorted(records, key=lambda r: r[risk_key], reverse=True)
    top_count = max(1, (total + 4) // 5)
    top = ranked[:top_count]
    lower = ranked[top_count:]
    summary = {
        "p_D": 0.07125,
        "frames": total,
        "logical_failures": failures,
        "logical_fer": failures / total,
        "risk_metric": risk_key,
        "risk_mean_success": (sum(r[risk_key] for r in records
                                   if not r["logical_failure"]) /
                              max(1, total - failures)),
        "risk_mean_failure": (sum(r[risk_key] for r in records
                                   if r["logical_failure"]) / max(1, failures)),
        "risk_auc": auc_positive_higher(records, risk_key),
        "top_20_percent": {"frames": len(top), "failures": sum(
            r["logical_failure"] for r in top)},
        "remaining_80_percent": {"frames": len(lower), "failures": sum(
            r["logical_failure"] for r in lower)},
        "batches": [{k: v for k, v in batch.items() if k != "records"}
                    for batch in results],
        "frame_csv": str(csv_path),
    }
    summary_path = args.run_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
