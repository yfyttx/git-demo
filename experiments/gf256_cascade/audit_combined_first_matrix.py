"""Independent-matrix paired audit of SF followed by layered fallback.

Replays only frames on which the saved SF run failed the Z syndrome. All other
frames are unchanged because the combined mode invokes layered BP only then.
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
STUDY = ROOT / ("experiments/gf256_cascade/results/"
                "path_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20352000")
OUT = ROOT / ("experiments/gf256_cascade/results/"
              "z_mechanism_audit_20261116/combined_first_matrix_layered_only")
MODE = "conditioned_robust_layered"


def frame_categories(log, count=5):
    categories = ["success"] * count
    for trial, x_syn, z_syn in re.findall(
            r"\[failure\] trial=(\d+).*? Xsyn=(\d+) Zsyn=(\d+)", log):
        trial = int(trial)
        if not int(x_syn):
            categories[trial] = "x_syndrome"
        elif not int(z_syn):
            categories[trial] = "z_syndrome"
        else:
            categories[trial] = "undetected_logical"
    return categories


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--executable", type=Path,
                        default=Path("/tmp/ldpc-h2-build/gf256_cascade_simulate"))
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    config = json.loads((STUDY / "config.json").read_text())
    matrix = ROOT / config["matrix_dir"] / "build"
    for name, digest in config["matrix_sha256"].items():
        if hashlib.sha256((matrix / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"matrix checksum mismatch: {name}")
    original = []
    for seed in range(config["seed_start"],
                      config["seed_start"] + config["frames_per_point"] // config["batch_size"]):
        filename = STUDY / (f"pD_0.07125_seed_{seed}_"
                            "conditioned_robust_tuned_syndrome_flip.log")
        categories = frame_categories(filename.read_text(), config["batch_size"])
        original.extend({"seed": seed, "trial": trial, "syndrome_flip": category}
                        for trial, category in enumerate(categories))
    baseline_z = sum(row["syndrome_flip"] == "z_syndrome" for row in original)
    baseline_fail = sum(row["syndrome_flip"] != "success" for row in original)
    if len(original) != 1000 or baseline_z != 32 or baseline_fail != 35:
        raise RuntimeError("saved 1000-frame SF run does not match published summary")
    selected = [row for row in original if row["syndrome_flip"] == "z_syndrome"]
    OUT.mkdir(exist_ok=True)
    (OUT / "config.json").write_text(json.dumps({
        "study": str(STUDY.relative_to(ROOT)), "mode": MODE,
        "selected": selected, "p_D": 0.07125,
        "reconstruction_rule": "for saved SF Z-syndrome failures only, choose layered if full syndrome passes; otherwise retain SF",
        "matrix_sha256": config["matrix_sha256"],
        "executable_sha256": hashlib.sha256(args.executable.read_bytes()).hexdigest(),
    }, indent=2) + "\n")

    def replay(row):
        seed, trial = row["seed"], row["trial"]
        stem = OUT / f"seed{seed}_trial{trial}"
        tsv, log = stem.with_suffix(".tsv"), stem.with_suffix(".log")
        if not (tsv.exists() and tsv.stat().st_size and
                log.exists() and log.stat().st_size):
            env = os.environ.copy()
            env.update({"GF256_CASCADE_SKIP_TARGET_DIAG": "1",
                        "GF256_FRAME_RESULTS": "1"})
            cmd = [str(args.executable), "5", str(seed), ".07125", MODE,
                   str(trial), "8", str(matrix), "12"]
            with tsv.open("w") as stdout, log.open("w") as stderr:
                subprocess.run(cmd, env=env, cwd=ROOT, stdout=stdout,
                               stderr=stderr, check=True)
        with tsv.open() as stream:
            data = list(csv.DictReader(stream, delimiter="\t"))
        if len(data) != 1 or int(data[0]["Total"]) != 1:
            raise RuntimeError(f"expected single-frame TSV: {tsv}")
        categories = frame_categories(log.read_text(),5)
        layered = categories[trial]
        if (layered == "success") != (int(data[0]["Failures"]) == 0):
            raise RuntimeError(f"log/TSV mismatch: {stem}")
        chosen = layered not in ("z_syndrome", "x_syndrome")
        return {**row, "layered": layered,
                "layered_chosen": int(chosen),
                "combined": layered if chosen else row["syndrome_flip"]}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        replays = list(pool.map(replay, selected))
    with (OUT / "selected_transitions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=replays[0].keys())
        writer.writeheader()
        writer.writerows(replays)
    improved = sum(row["combined"] == "success" for row in replays)
    undetected = sum(row["combined"] == "undetected_logical" for row in replays)
    report = {"frames": 1000, "syndrome_flip_failures": baseline_fail,
              "syndrome_flip_z_syndrome_failures": baseline_z,
              "selected_replays": len(replays), "rescued": improved,
              "converted_to_undetected": undetected,
              "still_z_syndrome": sum(row["combined"] == "z_syndrome"
                                      for row in replays),
              "combined_failures": baseline_fail - improved,
              "combined_undetected": 3 + undetected}
    (OUT / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
