#!/usr/bin/env python3
"""Run reproducible GF(256) cascade simulations at one or several p_D values."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILD_DEFAULT = Path(os.environ.get("GF256_CASCADE_BUILD", "/tmp/gf256_cascade_build"))
RESULTS = ROOT / "experiments/gf256_cascade/results/gf256_cascade_curve"
MATRIX_SEED = 20261018
# The p_D=0.06 baseline is reused from the verified 3050-frame run.
# Increase frame counts here for tighter uncertainty, especially at low FER.
SWEEP = [(0.06375, 500), (0.06750, 500), (0.07125, 300), (0.07500, 150)]
BATCH_SIZE_DEFAULT = 50
WORKERS_DEFAULT = min(os.cpu_count() or 1, 8)


def run_logged(command: list[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                       check=True, text=True)


def prior_point_is_available() -> int:
    paths = sorted((ROOT / "experiments/gf256_cascade/results").glob(
        "robust_202610*_50.tsv")) + sorted((ROOT / "experiments/gf256_cascade/results").glob(
        "robust_202610*_350.tsv"))
    frames = failures = 0
    for path in paths:
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) != 2:
            continue
        fields = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        if abs(float(fields["p_D"]) - 0.06) < 1e-12:
            frames += int(fields["Total"])
            failures += int(fields["Failures"])
    if frames != 3050 or failures != 0:
        raise RuntimeError(
            f"Expected the verified p_D=0.06 point (0/3050); found {failures}/{frames}. "
            "Restore the saved robust results before making the curve."
        )
    return frames


def make_tasks(run_dir: Path, batch_size: int, build: Path,
               sweep: list[tuple[float, int]]):
    tasks = []
    seed = 20262000
    for p_depolar, total_frames in sweep:
        point_dir = run_dir / "runs" / f"pD_{p_depolar:.5f}"
        point_dir.mkdir(parents=True, exist_ok=True)
        remaining, batch_id = total_frames, 0
        while remaining:
            count = min(batch_size, remaining)
            stem = f"batch_{batch_id:03d}_seed_{seed}"
            tsv_path = point_dir / f"{stem}.tsv"
            log_path = point_dir / f"{stem}.log"
            command = [str(build / "gf256_cascade_simulate"), str(count),
                       str(seed), f"{p_depolar:.8f}", "conditioned_robust"]
            tasks.append((command, run_dir / "matrices", log_path, tsv_path))
            remaining -= count
            batch_id += 1
            seed += 1
    return tasks


def execute_tasks(tasks, workers: int) -> None:
    stop_requested = threading.Event()
    children: set[subprocess.Popen] = set()
    children_lock = threading.Lock()

    def run_one(task):
        command, cwd, log_path, tsv_path = task
        tmp_path = tsv_path.with_suffix(".tsv.partial")
        started = time.monotonic()
        with tmp_path.open("w", encoding="utf-8") as output, \
                log_path.open("w", encoding="utf-8") as log:
            with children_lock:
                if stop_requested.is_set():
                    return tsv_path, 130, 0.0
                child = subprocess.Popen(command, cwd=cwd, stdout=output,
                                         stderr=log, text=True)
                children.add(child)
            try:
                return_code = child.wait()
            finally:
                with children_lock:
                    children.discard(child)
        elapsed = time.monotonic() - started
        if return_code == 0:
            tmp_path.replace(tsv_path)
        else:
            tmp_path.unlink(missing_ok=True)
        return tsv_path, return_code, elapsed

    # On SIGINT, terminate active child processes and preserve completed batches.
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_one, task): task[3] for task in tasks}
        try:
            for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
                tsv_path, rc, elapsed = future.result()
                print(f"[{index}/{len(tasks)}] {tsv_path.name}: exit={rc}, "
                      f"elapsed={elapsed / 60:.1f} min", flush=True)
                if rc:
                    raise RuntimeError(f"Simulation failed; inspect {tsv_path.with_suffix('.log')}")
        except KeyboardInterrupt:
            stop_requested.set()
            with children_lock:
                for child in children:
                    child.terminate()
            for future in futures:
                future.cancel()
            print("Interrupted. Completed batch TSV files are preserved; partial batches "
                  "are marked .tsv.partial. Active decoder processes were stopped.",
                  file=sys.stderr, flush=True)
            raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=WORKERS_DEFAULT,
                        help=f"parallel CPU processes (default: {WORKERS_DEFAULT})")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE_DEFAULT,
                        help=f"frames per subprocess (default: {BATCH_SIZE_DEFAULT})")
    parser.add_argument("--p-depolar", type=float,
                        help="run one p_D point instead of the default sweep")
    parser.add_argument("--frames", type=int,
                        help="frame count for --p-depolar")
    parser.add_argument("--build-dir", type=Path, default=BUILD_DEFAULT)
    args = parser.parse_args()
    if args.workers < 1 or args.batch_size < 1:
        parser.error("--workers and --batch-size must be positive")
    if (args.p_depolar is None) != (args.frames is None):
        parser.error("--p-depolar and --frames must be supplied together")
    if args.frames is not None and args.frames < 1:
        parser.error("--frames must be positive")
    if args.p_depolar is not None and not 0.0 < args.p_depolar < 2.0 / 3.0:
        parser.error("--p-depolar must be between 0 and 2/3")
    sweep = ([(args.p_depolar, args.frames)] if args.p_depolar is not None
             else SWEEP)

    baseline_frames = prior_point_is_available()
    subprocess.run(["cmake", "-S", str(ROOT), "-B", str(args.build_dir),
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run(["cmake", "--build", str(args.build_dir), "--target",
                    "gf256_cascade_construct", "gf256_cascade_simulate", "-j", "4"],
                   check=True)

    run_dir = RESULTS / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    matrix_dir = run_dir / "matrices"
    matrix_dir.mkdir(parents=True)
    run_logged([str(args.build_dir / "gf256_cascade_construct"), str(MATRIX_SEED)],
               matrix_dir, run_dir / "construct.log")
    subprocess.run(["sha256sum", "-c", str(
        ROOT / "experiments/gf256_cascade/matrix_hashes.sha256")],
        cwd=matrix_dir, check=True)

    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              check=True, capture_output=True, text=True).stdout.strip()
    config = {
        "matrix_seed": MATRIX_SEED,
        "decoder": "conditioned_robust",
        "workers": args.workers,
        "batch_size": args.batch_size,
        "baseline_p_D_0.06_frames": baseline_frames,
        "sweep_p_D_and_frames": sweep,
        "x_axis": "f_m = 2 p_D / 3",
        "git_revision": revision,
        "source_status_note": "Working-tree changes, if any, are not represented by git_revision.",
    }
    (run_dir / "run_config.json").write_text(json.dumps(config, indent=2) + "\n",
                                               encoding="utf-8")
    tasks = make_tasks(run_dir, args.batch_size, args.build_dir, sweep)
    print(f"Run directory: {run_dir}")
    print(f"Baseline: p_D=0.06, 0/{baseline_frames} failures (saved run)")
    print(f"Starting {sum(count for _, count in sweep)} new frames in {len(tasks)} batches "
          f"using {args.workers} CPU processes. Ctrl+C stops at batch boundaries.", flush=True)
    execute_tasks(tasks, args.workers)
    subprocess.run([sys.executable,
                    str(ROOT / "experiments/gf256_cascade/plot_gf256_cascade_sweep.py"),
                    "--run-dir", str(run_dir)], cwd=ROOT, check=True)
    print(f"Finished. Figure and source data are in {run_dir}")


if __name__ == "__main__":
    main()
