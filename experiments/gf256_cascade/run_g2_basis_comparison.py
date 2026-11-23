#!/usr/bin/env python3
"""Compare star and adjacent-pair G2 on paired GF(256) cascade frames.

The two G2 bases generate the same inner code. The script verifies that the
outer matrices and HZ are byte identical, then reports Tanner four-cycles,
syndrome failures, logical FER, and paired failure transitions.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import re
import subprocess
import time
from pathlib import Path


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def four_cycles(matrix_path: Path) -> int:
    with matrix_path.open(encoding="utf-8") as stream:
        cols, rows, _ = map(int, stream.readline().split())
        stream.readline()  # variable degrees
        stream.readline()  # check degrees
        shared = collections.Counter()
        for _ in range(rows):
            columns = [int(pair.split()[0]) for pair in re.findall(
                r"\d+\s+\d+", stream.readline())]
            if any(not 1 <= col <= cols for col in columns):
                raise ValueError(f"Invalid matrix entry in {matrix_path}")
            for i, first in enumerate(columns):
                for second in columns[i + 1:]:
                    shared[(first, second)] += 1
    return sum(count * (count - 1) // 2 for count in shared.values())


def matrix_rows(matrix_path: Path) -> list[dict[int, int]]:
    with matrix_path.open(encoding="utf-8") as stream:
        _, count, _ = map(int, stream.readline().split())
        stream.readline()
        stream.readline()
        return [{int(col): int(value) for col, value in re.findall(
            r"(\d+)\s+(\d+)", stream.readline())} for _ in range(count)]


def check_row_equivalence(star_path: Path, path_path: Path, inner_length: int) -> None:
    star = matrix_rows(star_path)
    path = matrix_rows(path_path)
    inner_dimension = inner_length - 1
    if len(star) != len(path) or len(star) % inner_dimension:
        raise ValueError("Unexpected HX row dimensions")
    for row, actual in enumerate(path):
        expected = star[row].copy()
        if row % inner_dimension:
            for col, coefficient in star[row - 1].items():
                if col in expected:
                    if expected[col] != coefficient:
                        raise ValueError("GF coefficients disagree in row transform")
                    del expected[col]  # equal GF(2^8) values sum to zero
                else:
                    expected[col] = coefficient
        if expected != actual:
            raise ValueError(f"HX row-equivalence failed at row {row}")


def read_result(tsv: Path, log: Path) -> tuple[dict[str, str], set[int]]:
    with tsv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one result row: {tsv}")
    failures = {int(index) for index in re.findall(
        r"\[failure\] trial=(\d+)", log.read_text(encoding="utf-8"))}
    if len(failures) != int(rows[0]["Failures"]):
        raise ValueError(f"Failure log disagrees with TSV: {log}")
    return rows[0], failures


def wilson(failures: int, frames: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = failures / frames
    denom = 1 + z * z / frames
    midpoint = (p + z * z / (2 * frames)) / denom
    halfwidth = z * math.sqrt(p * (1 - p) / frames + z*z / (4*frames*frames)) / denom
    return max(0.0, midpoint - halfwidth), min(1.0, midpoint + halfwidth)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20268000)
    parser.add_argument("--p-depolar", type=float, default=0.07125)
    parser.add_argument("--matrix-seed", type=int, default=20261018)
    parser.add_argument("--inner-length", type=int, default=6)
    parser.add_argument("--mode", default="conditioned_robust")
    parser.add_argument("--build-dir", type=Path, default=Path("/tmp/gf256_cascade_build"))
    parser.add_argument("--run-dir", type=Path, default=Path("/tmp/g2_path_comparison"))
    args = parser.parse_args()
    if args.frames <= 0 or args.inner_length < 2 or not 0 < args.p_depolar < 2/3:
        parser.error("frames must be positive, inner length >= 2, and p_D in (0, 2/3)")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    constructor = args.build_dir / "gf256_cascade_construct"
    simulator = args.build_dir / "gf256_cascade_simulate"
    for basis in ("star", "path"):
        root = args.run_dir / basis
        root.mkdir(exist_ok=True)
        matrix = root / "build" / "HX_Cascade_GF256.txt"
        if not matrix.exists():
            command = [str(constructor), str(args.matrix_seed), "2", "6", "97",
                       "35", "2", "8", str(root / "build"),
                       str(args.inner_length), basis]
            with (root / "construct.log").open("w", encoding="utf-8") as output:
                subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=True)
        construction_log = (root / "construct.log").read_text(encoding="utf-8")
        if (f"seed={args.matrix_seed}\n" not in construction_log or
                f"inner generator: {args.inner_length-1}x{args.inner_length} ({basis})"
                not in construction_log):
            raise ValueError(f"Existing {basis} matrices have a different construction")
    same_files = ("HX_TIT_Outer_For_HZ_GF256.txt",
                  "HZ_TIT_Outer_For_HX_GF256.txt", "HZ_Cascade_GF256.txt")
    for name in same_files:
        if file_hash(args.run_dir / "star" / "build" / name) != file_hash(
                args.run_dir / "path" / "build" / name):
            raise ValueError(f"Matrix is not paired: {name}")
    check_row_equivalence(
        args.run_dir / "star" / "build" / "HX_Cascade_GF256.txt",
        args.run_dir / "path" / "build" / "HX_Cascade_GF256.txt",
        args.inner_length)
    results = {}
    for basis in ("star", "path"):
        root = args.run_dir / basis
        result_file = root / f"{args.mode}_{args.frames}_{args.seed}_{args.p_depolar}.tsv"
        log_file = result_file.with_suffix(".log")
        command = [str(simulator), str(args.frames), str(args.seed),
                   str(args.p_depolar), args.mode, "-1", "8",
                   str(root / "build"), str(args.inner_length)]
        start = time.monotonic()
        with result_file.open("w", encoding="utf-8") as out, \
                log_file.open("w", encoding="utf-8") as err:
            subprocess.run(command, stdout=out, stderr=err, check=True)
        (root / f"{result_file.stem}.seconds").write_text(
            f"{time.monotonic() - start:.3f}\n", encoding="utf-8")
        row, failures = read_result(result_file, log_file)
        if int(row["Total"]) != args.frames or float(row["p_D"]) != args.p_depolar:
            raise ValueError(f"Result inputs do not match request: {result_file}")
        count = int(row["Failures"])
        results[basis] = {
            "frames": args.frames,
            "logical_failures": count,
            "logical_fer": count / args.frames,
            "wilson_95": wilson(count, args.frames),
            "x_syndrome_failures": int(row["XSyndromeFailures"]),
            "z_syndrome_failures": int(row["ZSyndromeFailures"]),
            "undetected_logical_failures": int(row["UndetectedLogicalFailures"]),
            "mean_physical_errors": float(row["MeanPhysicalErrors"]),
            "exact_x_frames": int(row["ExactXFrames"]),
            "four_cycles_symbol_graph": four_cycles(root / "build" / "HX_Cascade_GF256.txt"),
            "hx_sha256": file_hash(root / "build" / "HX_Cascade_GF256.txt"),
            "failure_trials": sorted(failures),
        }
    if (results["star"]["mean_physical_errors"] !=
            results["path"]["mean_physical_errors"] or
            results["star"]["exact_x_frames"] != results["path"]["exact_x_frames"]):
        raise ValueError("Paired channel/X outcomes do not match")
    star_fail = set(results["star"]["failure_trials"])
    path_fail = set(results["path"]["failure_trials"])
    rescued = len(star_fail - path_fail)
    regressed = len(path_fail - star_fail)
    discordant = rescued + regressed
    mcnemar_p = (min(1.0, 2 * sum(math.comb(discordant, k) for k in
                 range(min(rescued, regressed) + 1)) / 2**discordant)
                 if discordant else 1.0)
    report = {"matrix_seed": args.matrix_seed, "channel_seed": args.seed,
              "inner_length": args.inner_length,
              "p_D": args.p_depolar, "f_m": 2*args.p_depolar/3,
              "mode": args.mode, "results": results,
              "star_failed_path_succeeded": rescued,
              "star_succeeded_path_failed": regressed,
              "both_failed": len(star_fail & path_fail),
              "mcnemar_exact_two_sided_p": mcnemar_p}
    output = args.run_dir / f"summary_{args.mode}_{args.frames}_{args.seed}_{args.p_depolar}.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
