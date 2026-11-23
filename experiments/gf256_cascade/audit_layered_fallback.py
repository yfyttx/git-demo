"""Paired offline audit of layered BP as a fallback after tuned Z failure.

Only baseline Z-syndrome-failure frames are replayed. A layered result would
replace the baseline only if it passes the full Z syndrome. The criterion is
available to an online decoder; actual errors are consulted only afterward.
"""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/gf256_cascade/results/z_mechanism_audit_20261116"
STUDY = ROOT / ("experiments/gf256_cascade/results/"
                "matrix20261019_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20353000")


def classification(result):
    if int(result["ZSyndromeFailures"]):
        return "z_syndrome"
    if int(result["XSyndromeFailures"]):
        return "x_syndrome"
    if int(result["UndetectedLogicalFailures"]):
        return "undetected_logical"
    if int(result["Failures"]):
        raise RuntimeError("failure not explained by syndrome/logic fields")
    return "success"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path,
                        default=Path("/tmp/ldpc-h2-build/gf256_cascade_simulate"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=0,
                        help="First N Z-failure frames; 0 means all")
    args = parser.parse_args()
    if args.workers < 1 or args.limit < 0:
        parser.error("invalid workers or limit")
    original = list(csv.DictReader((BASE / "all_1000_frame_mechanisms.csv").open()))
    selected = [row for row in original if row["baseline"] == "z_syndrome"]
    if args.limit:
        selected = selected[:args.limit]
    config = json.loads((STUDY / "config.json").read_text())
    matrix = ROOT / config["matrix_dir"] / "build"
    for name, digest in config["matrix_sha256"].items():
        if hashlib.sha256((matrix / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"matrix checksum mismatch: {name}")
    out = BASE / "layered_fallback"
    out.mkdir(exist_ok=True)
    (out / "config.json").write_text(json.dumps({
        "selection": "baseline Z-syndrome failure only",
        "selected_frames": selected,
        "mode": "conditioned_robust_layered",
        "replacement_rule": "choose layered only if its full Z syndrome passes",
        "p_D": 0.07125, "matrix_sha256": config["matrix_sha256"],
        "executable_sha256": hashlib.sha256(args.executable.read_bytes()).hexdigest(),
    }, indent=2) + "\n")

    def replay(row):
        seed, trial = int(row["seed"]), int(row["trial"])
        stem = out / f"seed{seed}_trial{trial}"
        log, tsv = stem.with_suffix(".log"), stem.with_suffix(".tsv")
        if not (log.exists() and log.stat().st_size and
                tsv.exists() and tsv.stat().st_size):
            env = os.environ.copy()
            env.update({"GF256_CASCADE_SKIP_TARGET_DIAG": "1",
                        "GF256_FRAME_RESULTS": "1"})
            cmd = [str(args.executable), "5", str(seed), ".07125",
                   "conditioned_robust_layered", str(trial), "8",
                   str(matrix), "12"]
            with tsv.open("w") as stdout, log.open("w") as stderr:
                subprocess.run(cmd, cwd=ROOT, env=env, stdout=stdout,
                               stderr=stderr, check=True)
        with tsv.open() as stream:
            data = list(csv.DictReader(stream, delimiter="\t"))
        if len(data) != 1 or int(data[0]["Total"]) != 1:
            raise RuntimeError(f"expected one frame: {tsv}")
        layered = classification(data[0])
        choose_layered = layered not in ("z_syndrome", "x_syndrome")
        choose_after_sf = row["syndrome_flip"] == "z_syndrome" and choose_layered
        return {"seed": seed, "trial": trial, "baseline": row["baseline"],
                "syndrome_flip": row["syndrome_flip"], "layered": layered,
                "choose_layered": int(choose_layered),
                "hybrid": layered if choose_layered else row["baseline"],
                "syndrome_flip_then_layered": layered if choose_after_sf
                    else row["syndrome_flip"]}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        replays = list(pool.map(replay, selected))
    with (out / "selected_transitions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=replays[0].keys())
        writer.writeheader()
        writer.writerows(replays)
    full = args.limit == 0 or len(selected) == sum(
        row["baseline"] == "z_syndrome" for row in original)
    report = {
        "selected": len(replays), "complete_baseline_z_failures": full,
        "layered_syndrome_valid": sum(row["choose_layered"] for row in replays),
        "rescued": sum(row["hybrid"] == "success" for row in replays),
        "converted_to_undetected": sum(row["hybrid"] == "undetected_logical"
                                       for row in replays),
        "still_z_syndrome": sum(row["hybrid"] == "z_syndrome" for row in replays),
    }
    if full:
        transitions = {(row["seed"], row["trial"]): row for row in replays}
        hybrid = [transitions.get((int(row["seed"]), int(row["trial"])),
                                  {"hybrid": row["baseline"]})["hybrid"]
                  for row in original]
        combined = [transitions.get((int(row["seed"]), int(row["trial"])),
                                    {"syndrome_flip_then_layered": row["syndrome_flip"]})
                    ["syndrome_flip_then_layered"] for row in original]
        report.update({
            "all_frames": len(original),
            "baseline_logical_failures": sum(row["baseline"] != "success"
                                             for row in original),
            "syndrome_flip_logical_failures": sum(row["syndrome_flip"] != "success"
                                                  for row in original),
            "hybrid_logical_failures": sum(value != "success" for value in hybrid),
            "hybrid_undetected": sum(value == "undetected_logical" for value in hybrid),
            "syndrome_flip_then_layered_failures": sum(value != "success"
                                                        for value in combined),
            "syndrome_flip_then_layered_undetected": sum(
                value == "undetected_logical" for value in combined),
        })
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
