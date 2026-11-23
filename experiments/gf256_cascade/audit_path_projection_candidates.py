"""Offline projection-list coverage on fixed GF(256) failure frames.

The decoder and OSD list see only channel priors, matrix, and syndrome. The
simulator compares with the sampled error after decoding, for diagnostics.
This script is not a FER estimator or a decoder improvement.
"""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/gf256_cascade/results/z_mechanism_audit_20261116"
STUDY = ROOT / ("experiments/gf256_cascade/results/"
                "matrix20261019_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20353000")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path,
                        default=Path("/tmp/ldpc-h2-build/gf256_cascade_simulate"))
    parser.add_argument("--limit", type=int, default=3)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--free", type=int, default=32)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--alternatives", type=int, default=4)
    parser.add_argument("--mode", choices=("conditioned_robust_tuned",
                                            "conditioned_robust_layered"),
                        default="conditioned_robust_tuned")
    parser.add_argument("--source", choices=("outer", "joint", "both"),
                        default="both")
    parser.add_argument("--controls", action="store_true",
                        help="Audit the first rescued and successful frame")
    args = parser.parse_args()
    if args.limit < 1 or args.workers < 1:
        parser.error("limit and workers must be positive")
    config = json.loads((STUDY / "config.json").read_text())
    matrix = ROOT / config["matrix_dir"] / "build"
    for name, digest in config["matrix_sha256"].items():
        actual = hashlib.sha256((matrix / name).read_bytes()).hexdigest()
        if actual != digest:
            raise RuntimeError(f"matrix checksum mismatch: {name}")
    all_selected = json.loads((BASE / "mechanism_summary.json").read_text())["selected"]
    if args.controls:
        selected = [x for x in all_selected if x["group"] in ("rescued", "success")]
    else:
        selected = [x for x in all_selected if x["group"] == "persistent_z"][:args.limit]
    out = BASE / "path_projection_candidate_audit"
    if args.mode != "conditioned_robust_tuned":
        out = out / args.mode
    if args.controls:
        out = out / "controls"
    out.mkdir(exist_ok=True)
    sources = ("outer", "joint") if args.source == "both" else (args.source,)
    run_config = {
        "frames": selected, "sources": sources,
        "free": args.free, "order": args.order,
        "alternatives": args.alternatives,
        "mode": args.mode, "p_D": 0.07125,
        "executable_sha256": hashlib.sha256(args.executable.read_bytes()).hexdigest(),
        "matrix_sha256": config["matrix_sha256"],
    }
    (out / "config.json").write_text(json.dumps(run_config, indent=2) + "\n")

    def run(job):
        item, source = job
        stem = out / f"{source}_seed{item['seed']}_trial{item['trial']}"
        log_path, tsv_path = stem.with_suffix(".log"), stem.with_suffix(".tsv")
        if not (log_path.exists() and log_path.stat().st_size and
                tsv_path.exists() and tsv_path.stat().st_size):
            env = os.environ.copy()
            env.update({
                "GF256_CASCADE_SKIP_TARGET_DIAG": "1",
                "GF256_CASCADE_PROJECTION_RISK_AUDIT": "1",
                "GF256_CASCADE_PATH_CANDIDATE_AUDIT": "1",
                "GF256_CASCADE_PATH_CANDIDATE_SOURCE": source,
                "GF256_CASCADE_PATH_CANDIDATE_FREE": str(args.free),
                "GF256_CASCADE_PATH_CANDIDATE_ORDER": str(args.order),
                "GF256_CASCADE_PATH_CANDIDATE_ALTERNATIVES": str(args.alternatives),
            })
            cmd = [str(args.executable), "5", str(item["seed"]), ".07125",
                   args.mode, str(item["trial"]), "8",
                   str(matrix), "12"]
            with tsv_path.open("w") as stdout, log_path.open("w") as stderr:
                subprocess.run(cmd, cwd=ROOT, env=env, stdout=stdout,
                               stderr=stderr, check=True)
        log = log_path.read_text()
        matches = re.findall(
            r"\[outer-candidate-audit\] trial=\d+ stage=(\w+) "
            r"coordinate=(\d+) candidates=(\d+) "
            r"true_projection_covered=(\d+) min_symbol_distance=(-?\d+)", log)
        if not matches:
            raise RuntimeError(f"candidate audit missing: {log_path}")
        return [{"seed": item["seed"], "trial": item["trial"],
                 "group": item["group"],
                 "source": source, "stage": stage, "coordinate": int(coordinate),
                 "candidates": int(candidates), "covered": int(covered),
                 "min_symbol_distance": int(distance)}
                for stage, coordinate, candidates, covered, distance in matches]

    jobs = [(item, source) for item in selected for source in sources]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = [row for group in pool.map(run, jobs) for row in group]
    fields = ("seed", "trial", "group", "source", "stage", "coordinate",
              "candidates", "covered", "min_symbol_distance")
    with (out / "projection_coverage.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for item in selected:
        for source in sources:
            group = [r for r in rows if r["seed"] == item["seed"]
                     and r["source"] == source and r["stage"] == "initial"]
            print(f"seed={item['seed']} trial={item['trial']} source={source} "
                  f"covered={sum(r['covered'] for r in group)}/{len(group)} "
                  f"nearest={min(r['min_symbol_distance'] for r in group)}")


if __name__ == "__main__":
    main()
