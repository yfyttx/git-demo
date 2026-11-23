#!/usr/bin/env python3
"""Resumable paired star/path G2 experiment for hundreds of frames.

Each batch uses one seed for both G2 representations. Completed paired
batches are reused only when the inputs, binaries, and matrix hashes match.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_SCRIPT = HERE / "run_g2_basis_comparison.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wilson(failures: int, frames: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = failures / frames
    denominator = 1 + z*z/frames
    center = (p + z*z/(2*frames))/denominator
    radius = z*math.sqrt(p*(1-p)/frames + z*z/(4*frames*frames))/denominator
    return max(0.0, center-radius), min(1.0, center+radius)


def summary_path(args: argparse.Namespace, seed: int) -> Path:
    return args.run_dir / (
        f"summary_{args.mode}_{args.batch_size}_{seed}_{args.p_depolar}.json")


def read_cached(args: argparse.Namespace, seed: int) -> dict | None:
    path = summary_path(args, seed)
    if not path.exists():
        return None
    report = json.loads(path.read_text(encoding="utf-8"))
    expected = {"matrix_seed": args.matrix_seed, "channel_seed": seed,
                "inner_length": args.inner_length, "p_D": args.p_depolar,
                "mode": args.mode}
    if any(report.get(key) != value for key, value in expected.items()):
        raise ValueError(f"Cached batch has different inputs: {path}")
    for basis in ("star", "path"):
        result = report["results"][basis]
        matrix = args.run_dir / basis / "build" / "HX_Cascade_GF256.txt"
        if result["frames"] != args.batch_size or sha256(matrix) != result["hx_sha256"]:
            raise ValueError(f"Cached batch has different frame count or matrix: {path}")
    return report


def run_batch(args: argparse.Namespace, seed: int) -> dict:
    command = [sys.executable, str(PAIR_SCRIPT), "--inner-length",
               str(args.inner_length), "--frames", str(args.batch_size),
               "--seed", str(seed), "--p-depolar", str(args.p_depolar),
               "--matrix-seed", str(args.matrix_seed), "--mode", args.mode,
               "--build-dir", str(args.build_dir), "--run-dir", str(args.run_dir)]
    log = args.run_dir / f"batch_{seed}.log"
    error = args.run_dir / f"batch_{seed}.err"
    with log.open("w", encoding="utf-8") as out, \
            error.open("w", encoding="utf-8") as err:
        subprocess.run(command, stdout=out, stderr=err, check=True)
    result = read_cached(args, seed)
    if result is None:
        raise ValueError(f"Batch finished without a summary: {seed}")
    return result


def aggregate(args: argparse.Namespace, seeds: list[int]) -> dict:
    batches = [read_cached(args, seed) for seed in seeds]
    if any(batch is None for batch in batches):
        raise ValueError("Cannot aggregate incomplete batches")
    results = {}
    failed = {}
    for basis in ("star", "path"):
        entries = [batch["results"][basis] for batch in batches]
        frames = sum(entry["frames"] for entry in entries)
        failures = sum(entry["logical_failures"] for entry in entries)
        failed[basis] = {(seed, trial) for seed, entry in zip(seeds, entries)
                         for trial in entry["failure_trials"]}
        results[basis] = {
            "frames": frames,
            "logical_failures": failures,
            "logical_fer": failures / frames,
            "wilson_95": wilson(failures, frames),
            "x_syndrome_failures": sum(e["x_syndrome_failures"] for e in entries),
            "z_syndrome_failures": sum(e["z_syndrome_failures"] for e in entries),
            "undetected_logical_failures": sum(
                e["undetected_logical_failures"] for e in entries),
            "exact_x_frames": sum(e["exact_x_frames"] for e in entries),
            "mean_physical_errors": sum(
                e["mean_physical_errors"]*e["frames"] for e in entries)/frames,
            "four_cycles_symbol_graph": entries[0]["four_cycles_symbol_graph"],
            "hx_sha256": entries[0]["hx_sha256"],
        }
        if len(failed[basis]) != failures:
            raise ValueError(f"Failure count mismatch for {basis}")
    if results["star"]["mean_physical_errors"] != results["path"]["mean_physical_errors"]:
        raise ValueError("Channel frames were not paired")
    rescued = len(failed["star"] - failed["path"])
    regressed = len(failed["path"] - failed["star"])
    discordant = rescued + regressed
    p_value = (min(1.0, 2*sum(math.comb(discordant, k)
               for k in range(min(rescued, regressed)+1))/2**discordant)
               if discordant else 1.0)
    return {"p_D": args.p_depolar, "f_m": 2*args.p_depolar/3,
            "inner_length": args.inner_length, "mode": args.mode,
            "matrix_seed": args.matrix_seed, "seed_start": args.seed_start,
            "batch_size": args.batch_size, "batches": len(seeds),
            "results": results, "star_failed_path_succeeded": rescued,
            "star_succeeded_path_failed": regressed,
            "both_failed": len(failed["star"] & failed["path"]),
            "mcnemar_exact_two_sided_p": p_value}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--seed-start", type=int, default=20268200)
    parser.add_argument("--p-depolar", type=float, default=0.07125)
    parser.add_argument("--matrix-seed", type=int, default=20261018)
    parser.add_argument("--inner-length", type=int, default=12)
    parser.add_argument("--mode", default="conditioned_robust")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--build-dir", type=Path, default=Path("/tmp/gf256_cascade_build"))
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    if (args.frames < 1 or args.batch_size < 1 or args.frames % args.batch_size or
            args.inner_length < 2 or not 0 < args.p_depolar < 2/3 or
            args.workers < 1 or args.workers > (os.cpu_count() or 1)):
        parser.error("invalid frame, batch, inner-length, channel, or worker count")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    config = {"frames": args.frames, "batch_size": args.batch_size,
              "seed_start": args.seed_start, "p_D": args.p_depolar,
              "matrix_seed": args.matrix_seed, "inner_length": args.inner_length,
              "mode": args.mode, "simulator_sha256": sha256(
                  args.build_dir / "gf256_cascade_simulate"),
              "constructor_sha256": sha256(
                  args.build_dir / "gf256_cascade_construct"),
              "batch_script_sha256": sha256(PAIR_SCRIPT)}
    config_path = args.run_dir / "batched_config.json"
    if config_path.exists():
        if json.loads(config_path.read_text(encoding="utf-8")) != config:
            parser.error(f"Run directory belongs to different inputs: {config_path}")
    else:
        config_path.write_text(json.dumps(config, indent=2)+"\n", encoding="utf-8")
    seeds = [args.seed_start + index for index in range(args.frames//args.batch_size)]
    started = time.monotonic()

    # Build matrices once before concurrent workers read them.
    first = read_cached(args, seeds[0])
    if first is None:
        first = run_batch(args, seeds[0])
    print(f"[1/{len(seeds)}] seed={seeds[0]} "
          f"star={first['results']['star']['logical_failures']}/{args.batch_size} "
          f"path={first['results']['path']['logical_failures']}/{args.batch_size}",
          flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {}
        done_count = 1
        for seed in seeds[1:]:
            cached = read_cached(args, seed)
            if cached is None:
                futures[pool.submit(run_batch, args, seed)] = seed
            else:
                done_count += 1
                print(f"[{done_count}/{len(seeds)}] seed={seed} cached "
                      f"star={cached['results']['star']['logical_failures']}/{args.batch_size} "
                      f"path={cached['results']['path']['logical_failures']}/{args.batch_size}",
                      flush=True)
        for future in concurrent.futures.as_completed(futures):
            seed = futures[future]
            batch = future.result()
            done_count += 1
            print(f"[{done_count}/{len(seeds)}] seed={seed} "
                  f"star={batch['results']['star']['logical_failures']}/{args.batch_size} "
                  f"path={batch['results']['path']['logical_failures']}/{args.batch_size}",
                  flush=True)
    summary = aggregate(args, seeds)
    output = args.run_dir / "aggregate_summary.json"
    output.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    with (args.run_dir / "aggregate_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("basis", "frames", "logical_failures", "FER",
                         "x_syndrome_failures", "z_syndrome_failures",
                         "undetected_logical_failures", "wilson_95_lower",
                         "wilson_95_upper"))
        for basis in ("star", "path"):
            item = summary["results"][basis]
            writer.writerow((basis, item["frames"], item["logical_failures"],
                             item["logical_fer"], item["x_syndrome_failures"],
                             item["z_syndrome_failures"], item["undetected_logical_failures"],
                             *item["wilson_95"]))
    print(f"Finished in {(time.monotonic()-started)/60:.1f} min: {output}")


if __name__ == "__main__":
    main()
