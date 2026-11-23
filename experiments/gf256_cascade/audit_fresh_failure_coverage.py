"""Replay all remaining failures of the frozen three-mode study.

Candidate lists are diagnostic OSD-1 probes, not online decoder candidates.
Exact sampled projections and their outer-row-space cosets are scored only
after decoding. Per-projection coverage is not a chosen full-frame recovery.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
MODE = "conditioned_robust_tuned_syndrome_flip_layered"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=ROOT / "experiments/gf256_cascade/results/paired_three_modes_fresh1000_matrix20261019")
    parser.add_argument("--out", type=Path, default=ROOT / "experiments/gf256_cascade/results/fresh34_logical_coverage_audit_20261117")
    parser.add_argument("--executable", type=Path, default=Path("/tmp/gf256_cascade_build/gf256_cascade_simulate"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--neighbor", action="store_true", help="Diagnostic shared-symbol list feedback")
    parser.add_argument("--pair", action="store_true", help="Diagnostic joint physical-window OSD on one neighboring pair")
    parser.add_argument("--pair-all", action="store_true", help="Audit all neighboring physical windows; implies --pair")
    parser.add_argument("--limit", type=int, default=0, help="First N persistent failures; zero selects all plus controls")
    parser.add_argument("--start", type=int, default=0, help="Skip this many persistent failures before applying limit")
    args = parser.parse_args()
    if args.pair_all:
        args.pair = True
    cfg = json.loads((args.study / "config.json").read_text())
    matrix = Path(cfg["matrix_dir"]) / "build"
    for name, expected in cfg["matrix_sha256"].items():
        if digest(matrix / name) != expected:
            raise ValueError(f"Changed matrix: {name}")
    records = []
    for seed in range(cfg["seed_start"], cfg["seed_start"] + cfg["frames_per_point"] // cfg["batch_size"]):
        logs = {}
        for mode in cfg["modes"]:
            path = args.study / f"pD_0.07125_seed_{seed}_{mode}.log"
            logs[mode] = path.read_text()
        failures = {int(t): (int(x), int(z), int(exact)) for t, exact, x, z in re.findall(
            r"\[failure\] trial=(\d+).*?Xexact=(\d+) Xsyn=(\d+) Zsyn=(\d+)", logs[MODE])}
        base_failures = {int(t) for t in re.findall(r"\[failure\] trial=(\d+)", logs[cfg["modes"][0]])}
        sf = {int(t): (int(n), int(v)) for t, n, v in re.findall(
            r"\[syndrome-flip\] trial=(\d+) attempts=(\d+).*?valid_candidates=(\d+)", logs[MODE])}
        for trial in range(cfg["batch_size"]):
            x, z, exact = failures.get(trial, (1, 1, -1))
            group = ("persistent_z" if not z else "undetected" if trial in failures else
                     "rescued" if trial in base_failures else "success")
            attempts, valid = sf.get(trial, (0, 0))
            records.append(dict(seed=seed, trial=trial, group=group,
                                logical_failure=int(trial in failures), x_syndrome_ok=x,
                                z_syndrome_ok=z, exact_x=exact,
                                sf_attempts=attempts, sf_valid_candidates=valid))
    selected = [r for r in records if r["group"] == "persistent_z"]
    if args.limit:
        selected = selected[args.start:args.start+args.limit]
    for group in (() if args.limit else ("undetected", "rescued", "success")):
        selected.extend([r for r in records if r["group"] == group][:1])
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = dict(study=str(args.study), selected=selected, source="joint",
                    mode=MODE, p_D=0.07125, free=32, order=1, alternatives=4,
                    sf_bits=8, sf_max_weight=1, executable_sha256=digest(args.executable),
                    matrix_sha256=cfg["matrix_sha256"],
                    limitation="Per-projection outer-coset coverage of diagnostic lists; not a selected full-frame recovery or FER estimate")
    if args.neighbor or args.pair or args.limit or args.start:
        manifest.update(neighbor=args.neighbor, pair=args.pair, limit=args.limit,
                        start=args.start)
    if args.pair_all:
        manifest["pair_all"] = True
    config_path = args.out / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != manifest:
        raise ValueError("Existing audit configuration differs")
    config_path.write_text(json.dumps(manifest, indent=2) + "\n")
    executable = args.out / "gf256_cascade_simulate"
    shutil.copy2(args.executable, executable)
    shutil.copy2(Path(__file__), args.out / Path(__file__).name)
    shutil.copy2(ROOT / "experiments/gf256_cascade/simulate.cpp", args.out / "simulate.cpp")
    write_csv(args.out / "all_frames.csv", records)

    def run(item):
        stem = args.out / f"seed{item['seed']}_trial{item['trial']}"
        log, tsv = stem.with_suffix(".log"), stem.with_suffix(".tsv")
        if not tsv.exists():
            env = {k: v for k, v in os.environ.items() if not k.startswith("GF256_")}
            env.update(GF256_CASCADE_SKIP_TARGET_DIAG="1", GF256_FRAME_RESULTS="1",
                       GF256_CASCADE_PROJECTION_RISK_AUDIT="1", GF256_CASCADE_PATH_CANDIDATE_AUDIT="1",
                       GF256_CASCADE_PATH_CANDIDATE_SOURCE="joint", GF256_CASCADE_SF_BITS="8",
                       GF256_CASCADE_SF_MAX_WEIGHT="1")
            if args.neighbor:
                env["GF256_CASCADE_NEIGHBOR_CANDIDATE_AUDIT"] = "1"
            if args.pair:
                env["GF256_CASCADE_PAIR_CANDIDATE_AUDIT"] = "1"
            if args.pair_all:
                env["GF256_CASCADE_PAIR_ALL_AUDIT"] = "1"
            partial = tsv.with_suffix(".partial")
            command = [str(executable), str(cfg["batch_size"]), str(item["seed"]), ".07125",
                       MODE, str(item["trial"]), "8", str(matrix), "12"]
            with partial.open("w") as out, log.open("w") as err:
                subprocess.run(command, cwd=ROOT, env=env, stdout=out, stderr=err, check=True)
            partial.replace(tsv)
        data = log.read_text()
        match = re.search(r"\[frame-result\] trial=\d+ logical_failure=(\d+) x_syndrome_ok=(\d+) z_syndrome_ok=(\d+) exact_x=(\d+)", data)
        if not match:
            raise ValueError(f"Missing frame outcome: {log}")
        actual = tuple(map(int, match.groups()))
        expected = tuple(item[k] for k in ("logical_failure", "x_syndrome_ok", "z_syndrome_ok"))
        if actual[:3] != expected or (item["exact_x"] >= 0 and actual[3] != item["exact_x"]):
            raise ValueError(f"Replay differs from frozen study: {log}")
        coverage = [dict(seed=item["seed"], trial=item["trial"], group=item["group"],
                         stage=s, coordinate=int(c), candidates=int(n), covered=int(v), min_symbol_distance=int(d), coset_covered=int(k))
                    for s, c, n, v, d, k in re.findall(
                        r"\[outer-candidate-audit\] trial=\d+ stage=(\w+) coordinate=(\d+) candidates=(\d+) true_projection_covered=(\d+) min_symbol_distance=(-?\d+) true_outer_coset_covered=(\d+)", data)]
        projections = [dict(seed=item["seed"], trial=item["trial"], group=item["group"],
                            coordinate=int(c), row_weight=int(w), truth_syndrome=int(t),
                            decoded_syndrome=int(s), symbol_differences=int(d))
                       for c, w, t, s, d in re.findall(
                           r"\[outer-projection\].*?coordinate=(\d+) g2_row_weight=(\d+) true_syndrome=(\d+) decoded_syndrome=(\d+) symbol_differences=(\d+)", data)]
        if item["group"] == "persistent_z" and not coverage:
            raise ValueError(f"Missing failure coverage: {log}")
        return item, actual, coverage, projections

    coverage, projections, outcomes = [], [], []
    start = time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run, item) for item in selected]
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            item, actual, c, p = future.result()
            coverage.extend(c)
            projections.extend(p)
            outcomes.append({**item, "replay_exact_x": actual[3], "replay_matches": True})
            write_csv(args.out / "projection_coverage.csv", coverage)
            write_csv(args.out / "final_projections.csv", projections)
            write_csv(args.out / "replay_outcomes.csv", outcomes)
            print(f"[{i}/{len(selected)}] seed={item['seed']} trial={item['trial']} {item['group']} verified", flush=True)
    summary = dict(selected=len(selected), persistent_z=sum(r["group"] == "persistent_z" for r in selected),
                   replay_matches=len(outcomes), elapsed_seconds=time.monotonic()-start, stages={})
    for stage in sorted({r["stage"] for r in coverage}):
        rows = [r for r in coverage if r["stage"] == stage and r["group"] == "persistent_z"]
        keys = {(r["seed"], r["trial"]) for r in rows}
        summary["stages"][stage] = dict(frames=len(keys), projections=len(rows),
            truth_covered=sum(r["covered"] for r in rows),
            outer_coset_covered=sum(r["coset_covered"] for r in rows),
            frames_all_truth_covered=sum(all(r["covered"] for r in rows if (r["seed"],r["trial"])==key) for key in keys),
            frames_all_outer_cosets_covered=sum(all(r["coset_covered"] for r in rows if (r["seed"],r["trial"])==key) for key in keys))
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
