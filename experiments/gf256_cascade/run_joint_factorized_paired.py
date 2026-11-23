#!/usr/bin/env python3
"""Compare staged X->Z decoding with simultaneous factorized X/Z BP on paired frames."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path


FRAME_PATTERN = re.compile(
    r"\[frame-result\] trial=(\d+) logical_failure=([01]) "
    r"x_syndrome_ok=([01]) z_syndrome_ok=([01]) "
    r"exact_x=([01]) physical_errors=(\d+)"
)
MODES = ("conditioned_robust", "joint_factorized")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_mode(args: argparse.Namespace, mode: str) -> tuple[dict[int, dict], dict]:
    env = os.environ.copy()
    env["GF256_FRAME_RESULTS"] = "1"
    env["GF256_JOINT_ITERS"] = str(args.joint_iterations)
    command = [
        str(args.binary), str(args.frames), str(args.seed), str(args.p_depolar),
        mode, "-1", "8", str(args.matrix_dir), "6",
    ]
    started = time.perf_counter()
    completed = subprocess.run(command, env=env, text=True, capture_output=True)
    wall_seconds = time.perf_counter() - started
    (args.run_dir / f"{mode}.stdout.tsv").write_text(completed.stdout, encoding="utf-8")
    (args.run_dir / f"{mode}.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{mode} exited {completed.returncode}; see stderr log")
    frames: dict[int, dict] = {}
    for match in FRAME_PATTERN.finditer(completed.stderr):
        trial, fail, x_ok, z_ok, exact_x, errors = map(int, match.groups())
        if trial in frames:
            raise RuntimeError(f"Duplicate frame {trial} in {mode}")
        frames[trial] = {
            "logical_failure": fail,
            "x_syndrome_ok": x_ok,
            "z_syndrome_ok": z_ok,
            "exact_x": exact_x,
            "physical_errors": errors,
        }
    if set(frames) != set(range(args.frames)):
        raise RuntimeError(f"{mode} produced {len(frames)} of {args.frames} frame records")
    summary = {
        "mode": mode,
        "wall_seconds": wall_seconds,
        "seconds_per_frame": wall_seconds / args.frames,
        "logical_failures": sum(row["logical_failure"] for row in frames.values()),
        "x_syndrome_failures": sum(1 - row["x_syndrome_ok"] for row in frames.values()),
        "z_syndrome_failures": sum(1 - row["z_syndrome_ok"] for row in frames.values()),
        "exact_x_frames": sum(row["exact_x"] for row in frames.values()),
    }
    return frames, summary


def exact_mcnemar(rescued: int, regressed: int) -> float:
    n = rescued + regressed
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, k) for k in range(min(rescued, regressed) + 1))
    return min(1.0, 2.0 * tail / 2**n)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20265000)
    parser.add_argument("--p-depolar", type=float, default=0.07125)
    parser.add_argument("--joint-iterations", type=int, default=12)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.matrix_dir = args.matrix_dir.resolve()
    args.run_dir.mkdir(parents=True, exist_ok=False)
    if args.frames <= 0 or args.joint_iterations <= 0:
        parser.error("frames and joint-iterations must be positive")

    results = {}
    summaries = []
    for mode in MODES:
        print(f"Running {mode} on {args.frames} frames...", flush=True)
        results[mode], summary = run_mode(args, mode)
        summaries.append(summary)
        print(f"  failures={summary['logical_failures']}/{args.frames}, "
              f"Z syndrome={summary['z_syndrome_failures']}, "
              f"wall={summary['wall_seconds']:.1f}s", flush=True)

    baseline, joint = (results[mode] for mode in MODES)
    rows = []
    for trial in range(args.frames):
        a, b = baseline[trial], joint[trial]
        if a["physical_errors"] != b["physical_errors"]:
            raise RuntimeError(f"Noise mismatch at frame {trial}")
        rows.append({
            "trial": trial,
            "physical_errors": a["physical_errors"],
            **{f"baseline_{key}": value for key, value in a.items()
               if key != "physical_errors"},
            **{f"joint_{key}": value for key, value in b.items()
               if key != "physical_errors"},
        })
    with (args.run_dir / "paired_frames.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    rescued = sum(row["baseline_logical_failure"] and not row["joint_logical_failure"]
                  for row in rows)
    regressed = sum(not row["baseline_logical_failure"] and row["joint_logical_failure"]
                    for row in rows)
    metadata = {
        "p_depolar": args.p_depolar,
        "f_m": 2 * args.p_depolar / 3,
        "frames": args.frames,
        "seed": args.seed,
        "joint_iterations": args.joint_iterations,
        "binary_sha256": sha256(args.binary),
        "matrix_sha256": {
            path.name: sha256(path)
            for path in sorted(args.matrix_dir.glob("*.txt"))
        },
        "summary": summaries,
        "joint_rescued_vs_baseline": rescued,
        "joint_regressed_vs_baseline": regressed,
        "mcnemar_two_sided_p": exact_mcnemar(rescued, regressed),
    }
    (args.run_dir / "summary.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Paired frames verified; rescued={rescued}, regressed={regressed}, "
          f"McNemar p={metadata['mcnemar_two_sided_p']:.4g}")
    print(args.run_dir)


if __name__ == "__main__":
    main()
