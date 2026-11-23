#!/usr/bin/env python3
"""Resumable, paired GF(256) cascade decoder comparison.

Each decoder receives the same matrix, channel parameter, seed, and frame
indices.  Completed batches are reused only if their recorded inputs match.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("GF256_CASCADE_BUILD", "/tmp/gf256_cascade_build"))
MATRICES = (ROOT / "experiments/gf256_cascade/results/gf256_cascade_curve"
            / "run_20261025_001136/matrices")
MODES = ("conditioned_robust", "conditioned_robust_refresh",
         "conditioned_robust_joint_projection",
         "conditioned_robust_joint_list",
         "conditioned_robust_adaptive", "conditioned_robust_guess",
         "conditioned_robust_osd2_weak",
         "conditioned_robust_bpgd", "conditioned_robust_restart_belief",
         "conditioned_robust_restart_belief_global",
         "conditioned_robust_tuned_restart_belief_global",
         "conditioned_robust_tuned_syndrome_flip",
         "conditioned_robust_tuned_syndrome_flip_layered",
         "conditioned_robust_long", "conditioned_robust_tuned",
         "conditioned_robust_softx", "conditioned_robust_structured_softx",
         "conditioned_robust_augmented",
         "conditioned_robust_layered",
         "conditioned_robust_overcomplete",
         "conditioned_robust_si",
         "conditioned_robust_joint_projection_safe",
         "conditioned_fullhx_bp", "conditioned_projection_independent")
FAILURE_RE = re.compile(r"\[failure\] trial=(\d+)")
FAILURE_TYPE_RE = re.compile(r"\[failure\] trial=(\d+).*?Xsyn=(\d+) Zsyn=(\d+)")
TIME_RE = re.compile(r"user=(\S+) system=(\S+) elapsed=(\S+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_result(path: Path, log_path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 2:
        raise ValueError(f"Unexpected TSV format: {path}")
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    log = log_path.read_text(encoding="utf-8")
    failures = {int(value) for value in FAILURE_RE.findall(log)}
    if len(failures) != int(row["Failures"]):
        raise ValueError(f"Failure log and TSV disagree: {path}")
    failure_types = {}
    for trial, x_syndrome, z_syndrome in FAILURE_TYPE_RE.findall(log):
        failure_types[int(trial)] = ("x_syndrome" if x_syndrome == "0" else
                                     "z_syndrome" if z_syndrome == "0" else
                                     "undetected_logical")
    if set(failure_types) != failures:
        raise ValueError(f"Failure type log and TSV disagree: {path}")
    match = TIME_RE.search(log)
    if match is None:
        raise ValueError(f"Missing CPU timing: {log_path}")
    row["failure_trials"] = failures
    row["failure_types"] = failure_types
    row["cpu_seconds"] = float(match.group(1)) + float(match.group(2))
    row["wall_seconds"] = float(match.group(3))
    return row


def wilson(failures: int, frames: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = failures / frames
    denom = 1 + z * z / frames
    center = (p + z * z / (2 * frames)) / denom
    radius = z * math.sqrt(p * (1 - p) / frames + z * z / (4 * frames**2)) / denom
    lower = 0.0 if failures == 0 else max(0.0, center - radius)
    upper = 1.0 if failures == frames else min(1.0, center + radius)
    return lower, upper


def mcnemar_exact(rescued: int, regressed: int) -> float:
    n = rescued + regressed
    if n == 0:
        return 1.0
    # Exact, two-sided conditional binomial test for paired binary outcomes.
    lower = sum(math.comb(n, k) for k in range(min(rescued, regressed) + 1))
    return min(1.0, 2 * lower / 2**n)


def task_paths(run_dir: Path, p: float, seed: int, mode: str) -> tuple[Path, Path, Path]:
    stem = f"pD_{p:.5f}_seed_{seed}_{mode}"
    return (run_dir / f"{stem}.tsv", run_dir / f"{stem}.log",
            run_dir / f"{stem}.meta.json")


def run_one(task: tuple, run_dir: Path, executable: Path, matrices: Path,
            executable_hash: str, matrix_hashes: dict[str, str],
            inner_length: int, independent_outer_iterations: int = 20) -> tuple:
    p, frames, seed, mode = task
    tsv, log, meta = task_paths(run_dir, p, seed, mode)
    expected = {"p_D": p, "frames": frames, "seed": seed, "mode": mode,
                "executable_sha256": executable_hash, "matrix_sha256": matrix_hashes}
    if inner_length != 6:
        expected["inner_length"] = inner_length
    if mode == "conditioned_projection_independent":
        expected["independent_outer_iterations"] = independent_outer_iterations
    if tsv.exists() and log.exists() and meta.exists():
        if json.loads(meta.read_text(encoding="utf-8")) != expected:
            raise ValueError(f"Existing batch has different inputs: {tsv}")
        row = read_result(tsv, log)
        if int(row["Total"]) != frames or float(row["p_D"]) != p:
            raise ValueError(f"Existing batch has different output: {tsv}")
        return (task, row, True)
    tmp_tsv = tsv.with_suffix(".tsv.partial")
    tmp_log = log.with_suffix(".log.partial")
    command = ["/usr/bin/time", "-f", "user=%U system=%S elapsed=%e",
               str(executable), str(frames), str(seed), str(p), mode]
    if inner_length != 6:
        command.extend(("-1", "8", "build", str(inner_length)))
    with tmp_tsv.open("w", encoding="utf-8") as output, \
         tmp_log.open("w", encoding="utf-8") as errors:
        env = os.environ.copy()
        if mode == "conditioned_projection_independent":
            env["GF256_INDEPENDENT_OUTER_ITERS"] = str(independent_outer_iterations)
        subprocess.run(command, cwd=matrices, env=env,
                       stdout=output, stderr=errors, check=True)
    tmp_tsv.replace(tsv)
    tmp_log.replace(log)
    meta.write_text(json.dumps(expected, indent=2) + "\n", encoding="utf-8")
    return (task, read_result(tsv, log), False)


def summarize(run_dir: Path, points: list[float], seeds: int,
              seed_start: int, modes: list[str]) -> None:
    rows = []
    pairwise_rows = []
    transition_rows = []
    for index, p in enumerate(points):
        paired = {}
        for batch in range(seeds):
            seed = seed_start + index * seeds + batch
            paired[seed] = {}
            for mode in modes:
                tsv, log, _ = task_paths(run_dir, p, seed, mode)
                paired[seed][mode] = read_result(tsv, log)
            reference = paired[seed][modes[0]]
            for mode in modes[1:]:
                other = paired[seed][mode]
                # A decoder variant may intentionally alter its X estimate.
                # Pairing requires identical sampled frames, not identical
                # decoding outcomes. Physical-error totals are a direct check
                # that the seeded channel draws match; X performance is an
                # outcome metric and is summarized separately.
                if (reference["Total"] != other["Total"] or
                    reference["MeanPhysicalErrors"] != other["MeanPhysicalErrors"]):
                    raise ValueError(f"Paired channel check failed for seed {seed}")
        reference_fails = {(seed, trial) for seed, records in paired.items()
                           for trial in records[modes[0]]["failure_trials"]}
        failure_sets = {mode: {(seed, trial) for seed, records in paired.items()
                               for trial in records[mode]["failure_trials"]}
                        for mode in modes}
        for index_a, mode_a in enumerate(modes):
            for mode_b in modes[index_a + 1:]:
                a_only = len(failure_sets[mode_a] - failure_sets[mode_b])
                b_only = len(failure_sets[mode_b] - failure_sets[mode_a])
                pairwise_rows.append({"p_D": p, "mode_a": mode_a, "mode_b": mode_b,
                                      "a_fail_b_success": a_only,
                                      "a_success_b_fail": b_only,
                                      "mcnemar_two_sided_p": mcnemar_exact(a_only, b_only)})
                transitions = {}
                for seed, records in paired.items():
                    for trial in range(int(records[mode_a]["Total"])):
                        category_a = records[mode_a]["failure_types"].get(trial, "success")
                        category_b = records[mode_b]["failure_types"].get(trial, "success")
                        key = (category_a, category_b)
                        transitions[key] = transitions.get(key, 0) + 1
                for (category_a, category_b), count in sorted(transitions.items()):
                    transition_rows.append({"p_D": p, "mode_a": mode_a,
                                            "mode_b": mode_b,
                                            "outcome_a": category_a,
                                            "outcome_b": category_b,
                                            "frames": count})
        for mode in modes:
            records = [paired[seed][mode] for seed in paired]
            failures = failure_sets[mode]
            total = sum(int(record["Total"]) for record in records)
            lo, hi = wilson(len(failures), total)
            rescued = len(reference_fails - failures)
            regressed = len(failures - reference_fails)
            rows.append({"p_D": p, "f_m": 2*p/3, "mode": mode, "frames": total,
                         "logical_failures": len(failures), "FER": len(failures)/total,
                         "wilson_95_lower": lo, "wilson_95_upper": hi,
                         "z_syndrome_failures": sum(int(r["ZSyndromeFailures"]) for r in records),
                         "undetected_logical_failures": sum(int(r["UndetectedLogicalFailures"]) for r in records),
                         "cpu_seconds": sum(r["cpu_seconds"] for r in records),
                         "cpu_seconds_per_frame": sum(r["cpu_seconds"] for r in records)/total,
                         "rescued_vs_baseline": rescued, "regressed_vs_baseline": regressed,
                         "mcnemar_two_sided_p": mcnemar_exact(rescued, regressed)})
    with (run_dir / "summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (run_dir / "pairwise.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(pairwise_rows[0]))
        writer.writeheader()
        writer.writerows(pairwise_rows)
    with (run_dir / "failure_transitions.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(transition_rows[0]))
        writer.writeheader()
        writer.writerows(transition_rows)
    print(f"Summary: {run_dir / 'summary.csv'}")
    for row in rows:
        print(f"p_D={row['p_D']:.5f} {row['mode']} "
              f"FER={row['logical_failures']}/{row['frames']}={row['FER']:.5g} "
              f"CPU={row['cpu_seconds_per_frame']:.2f}s/frame "
              f"rescued={row['rescued_vs_baseline']} regressed={row['regressed_vs_baseline']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-depolar", nargs="+", type=float,
                        default=[0.0675, 0.07125, 0.075])
    parser.add_argument("--frames-per-point", type=int, default=500)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--seed-start", type=int, default=20263000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--modes", nargs="+", choices=MODES,
                        default=["conditioned_robust", "conditioned_robust_adaptive",
                                 "conditioned_robust_long"])
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--build-dir", type=Path, default=BUILD)
    parser.add_argument("--matrix-dir", type=Path, default=MATRICES)
    parser.add_argument("--inner-length", type=int, default=6)
    parser.add_argument("--independent-outer-iters", type=int, default=20,
                        help="BP iterations for the direct G2-to-outer mode (1..200)")
    parser.add_argument("--allow-custom-matrices", action="store_true",
                        help="Use matrix hashes from --matrix-dir instead of the fixed reference set")
    args = parser.parse_args()
    args.build_dir = args.build_dir.resolve()
    args.matrix_dir = args.matrix_dir.resolve()
    if not 1 <= args.independent_outer_iters <= 200:
        parser.error("--independent-outer-iters must be in 1..200")
    if args.frames_per_point < 1 or args.batch_size < 1 or args.workers < 1 or args.inner_length < 2:
        parser.error("frame count, batch size, and workers must be positive")
    if args.frames_per_point % args.batch_size:
        parser.error("--frames-per-point must be divisible by --batch-size")
    if not 1 <= args.workers <= os.cpu_count():
        parser.error("--workers exceeds available CPU count")
    if len(args.p_depolar) != len(set(args.p_depolar)) or any(
            not 0 < p < 2/3 for p in args.p_depolar):
        parser.error("p_D values must be distinct and in (0, 2/3)")
    modes = list(dict.fromkeys(args.modes))
    if modes[0] not in MODES:
        parser.error("first mode must be a supported decoder mode")
    if len(modes) < 2:
        parser.error("at least two distinct decoder modes are required")
    executable = args.build_dir / "gf256_cascade_simulate"
    if not executable.is_file():
        parser.error(f"missing executable: {executable}")
    matrix_hashes = {name: sha256(args.matrix_dir / "build" / name) for name in (
        "HX_Cascade_GF256.txt", "HZ_Cascade_GF256.txt",
        "HX_TIT_Outer_For_HZ_GF256.txt", "HZ_TIT_Outer_For_HX_GF256.txt")}
    if not args.allow_custom_matrices:
        subprocess.run(["sha256sum", "-c", str(ROOT / "experiments/gf256_cascade/matrix_hashes.sha256")],
                       cwd=args.matrix_dir, check=True, stdout=subprocess.DEVNULL)
    run_dir = args.run_dir or (ROOT / "experiments/gf256_cascade/results/paired_decoder_study"
                               / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    run_dir.mkdir(parents=True, exist_ok=True)
    config = {"p_depolar": args.p_depolar, "frames_per_point": args.frames_per_point,
              "batch_size": args.batch_size, "seed_start": args.seed_start,
              "modes": modes, "executable_sha256": sha256(executable),
              "matrix_sha256": matrix_hashes, "matrix_dir": str(args.matrix_dir)}
    if args.inner_length != 6 or args.allow_custom_matrices:
        config.update({"inner_length": args.inner_length,
                       "allow_custom_matrices": args.allow_custom_matrices})
    if "conditioned_projection_independent" in modes:
        config["independent_outer_iterations"] = args.independent_outer_iters
    config_path = run_dir / "config.json"
    if config_path.exists() and json.loads(config_path.read_text(encoding="utf-8")) != config:
        parser.error(f"run directory has different configuration: {config_path}")
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    seeds = args.frames_per_point // args.batch_size
    tasks = [(p, args.batch_size, args.seed_start + index*seeds + batch, mode)
             for index, p in enumerate(args.p_depolar)
             for batch in range(seeds) for mode in modes]
    print(f"Running {len(tasks)} batches in {run_dir}", flush=True)
    started = time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, task, run_dir, executable, args.matrix_dir,
                               config["executable_sha256"], matrix_hashes,
                               args.inner_length, args.independent_outer_iters): task
                   for task in tasks}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            task, row, cached = future.result()
            print(f"[{done}/{len(tasks)}] p_D={task[0]:.5f} seed={task[2]} "
                  f"{task[3]} fail={row['Failures']}/{row['Total']} "
                  f"CPU={row['cpu_seconds']:.1f}s cached={cached}", flush=True)
    summarize(run_dir, args.p_depolar, seeds, args.seed_start, modes)
    print(f"Elapsed wall time: {(time.monotonic()-started)/60:.1f} min")


if __name__ == "__main__":
    main()
