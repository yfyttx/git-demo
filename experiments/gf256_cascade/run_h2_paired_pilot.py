#!/usr/bin/env python3
"""Paired pilot for two precomputed GF(256), length-12 inner codes."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODE = "conditioned_robust_adaptive"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_one(label: str, matrix_dir: Path, seed: int, frames: int,
            p_depolar: float, executable: Path, run_dir: Path,
            mode: str) -> dict:
    base = run_dir / f"{label}_{seed}"
    tsv = base.with_suffix(".tsv")
    log = base.with_suffix(".log")
    if not tsv.exists() or not log.exists():
        command = [str(executable), str(frames), str(seed), str(p_depolar),
                   mode, "-1", "8", str(matrix_dir), "12"]
        with tsv.open("w", encoding="utf-8") as out, \
                log.open("w", encoding="utf-8") as err:
            subprocess.run(command, cwd=ROOT, stdout=out, stderr=err, check=True)
    with tsv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1 or int(rows[0]["Total"]) != frames:
        raise ValueError(f"Incomplete result in {tsv}")
    failures = {int(index) for index in re.findall(
        r"\[failure\] trial=(\d+)", log.read_text(encoding="utf-8"))}
    if len(failures) != int(rows[0]["Failures"]):
        raise ValueError(f"Failure count disagrees with {log}")
    return {"label": label, "seed": seed, "frames": frames,
            "failures": len(failures), "failure_trials": sorted(failures),
            "z_syndrome_failures": int(rows[0]["ZSyndromeFailures"]),
            "x_syndrome_failures": int(rows[0]["XSyndromeFailures"]),
            "mean_physical_errors": float(rows[0]["MeanPhysicalErrors"]),
            "tsv": str(tsv), "log": str(log)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames-per-seed", type=int, default=5)
    parser.add_argument("--seeds", type=int, default=4)
    parser.add_argument("--seed-start", type=int, default=20272000)
    parser.add_argument("--p-depolar", type=float, default=0.07125)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mode", default=DEFAULT_MODE,
                        help="decoder mode used identically for both matrices")
    parser.add_argument("--executable", type=Path,
                        default=Path("/tmp/ldpc-h2-build/gf256_cascade_simulate"))
    parser.add_argument("--old-matrices", type=Path,
                        default=ROOT / "experiments/gf256_cascade/results/"
                        "g2_basis_path_n12_300frames/path/build")
    parser.add_argument("--new-matrices", type=Path,
                        default=ROOT / "experiments/gf256_cascade/results/"
                        "h2_path2_n12/build")
    parser.add_argument("--run-dir", type=Path,
                        default=ROOT / "experiments/gf256_cascade/results/"
                        "h2_path2_n12/paired_20")
    args = parser.parse_args()
    if args.frames_per_seed < 1 or args.seeds < 1 or args.workers < 1:
        parser.error("frames, seeds, and workers must be positive")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    inputs = {"mode": args.mode, "frames_per_seed": args.frames_per_seed,
              "seeds": args.seeds, "seed_start": args.seed_start,
              "p_depolar": args.p_depolar,
              "executable_sha256": sha256(args.executable),
              "matrices": {label: {"directory": str(directory),
                "hx_sha256": sha256(directory / "HX_Cascade_GF256.txt"),
                "hz_sha256": sha256(directory / "HZ_Cascade_GF256.txt")}
                for label, directory in (("old", args.old_matrices),
                                         ("new", args.new_matrices))}}
    config_path = args.run_dir / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != inputs:
        parser.error("Existing run directory has different inputs")
    config_path.write_text(json.dumps(inputs, indent=2) + "\n")
    tasks = [(label, directory, args.seed_start + offset,
              args.frames_per_seed, args.p_depolar, args.executable, args.run_dir,
              args.mode)
             for offset in range(args.seeds)
             for label, directory in (("old", args.old_matrices),
                                      ("new", args.new_matrices))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_one, *task) for task in tasks]
        results = []
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            results.append(result)
            print(f"[{index}/{len(tasks)}] {result['label']} seed={result['seed']} "
                  f"fail={result['failures']}/{result['frames']}", flush=True)
    by_key = {(item["label"], item["seed"]): item for item in results}
    comparison = {label: {"frames": args.seeds * args.frames_per_seed,
        "logical_failures": sum(by_key[label, args.seed_start + i]["failures"]
                                for i in range(args.seeds)),
        "z_syndrome_failures": sum(
            by_key[label, args.seed_start + i]["z_syndrome_failures"]
            for i in range(args.seeds))}
        for label in ("old", "new")}
    paired = {"old_failed_new_succeeded": 0,
              "old_succeeded_new_failed": 0, "both_failed": 0}
    for i in range(args.seeds):
        seed = args.seed_start + i
        old = set(by_key["old", seed]["failure_trials"])
        new = set(by_key["new", seed]["failure_trials"])
        paired["old_failed_new_succeeded"] += len(old - new)
        paired["old_succeeded_new_failed"] += len(new - old)
        paired["both_failed"] += len(old & new)
    summary = {"comparison": comparison, "paired": paired,
               "batches": sorted(results, key=lambda item: (item["seed"], item["label"]))}
    (args.run_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"comparison": comparison, "paired": paired}, indent=2))


if __name__ == "__main__":
    main()
