#!/usr/bin/env python3
"""Run a fixed GF(256) cascade code over Z-biased Pauli noise and plot FER.

For each (p, eta), p_X = p_Y = p/(eta+2), p_Z = eta*p/(eta+2).
eta is p_Z/p_X; eta=1 is the depolarizing channel. This convention differs
from the Y-bias convention used in Tuckett et al. (2019).
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
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MATRICES = (ROOT / "experiments/gf256_cascade/results/"
                    "g2_basis_path_n12_300frames/path")
DEFAULT_RUN = (ROOT / "experiments/gf256_cascade/results/"
               "z_bias_fer_curve_n12_path")
REQUIRED_MATRICES = (
    "HX_Cascade_GF256.txt", "HZ_Cascade_GF256.txt",
    "HX_Cascade_Physical_Binary.txt", "HZ_Cascade_Physical_Binary.txt",
    "HX_TIT_Outer_For_HZ_GF256.txt", "HZ_TIT_Outer_For_HX_GF256.txt",
)
MODE = "conditioned_robust_adaptive"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one TSV result row: {path}")
    return rows[0]


def row_p(row: dict[str, str]) -> float:
    return float(row["p_total"] if "p_total" in row else row["p_D"])


def wilson95(failures: int, frames: int) -> tuple[float, float]:
    z = 1.959963984540054
    phat = failures / frames
    denom = 1 + z*z/frames
    center = (phat + z*z/(2*frames))/denom
    radius = z*math.sqrt(phat*(1-phat)/frames + z*z/(4*frames*frames))/denom
    return max(0.0, center-radius), min(1.0, center+radius)


def run_batch(task: tuple[float, float, int, int], args: argparse.Namespace,
              executable: Path, run_dir: Path) -> tuple[float, float, int, int]:
    bias, p, seed, frames = task
    stem = f"eta_{bias:g}_p_{p:.5f}_seed_{seed}"
    tsv = run_dir / "batches" / f"{stem}.tsv"
    log = run_dir / "batches" / f"{stem}.log"
    if tsv.exists():
        row = read_row(tsv)
        if int(row["Total"]) != frames or abs(row_p(row)-p) > 1e-10:
            raise ValueError(f"Saved batch has inconsistent inputs: {tsv}")
        return task
    temporary = tsv.with_suffix(".tsv.partial")
    env = os.environ.copy()
    env["GF256_Z_BIAS"] = str(bias)
    command = [str(executable), str(frames), str(seed), str(p), MODE,
               "-1", "8", "build", "12"]
    with temporary.open("w", encoding="utf-8") as out, log.open("w", encoding="utf-8") as err:
        subprocess.run(command, cwd=args.matrix_dir, env=env,
                       stdout=out, stderr=err, check=True)
    row = read_row(temporary)
    if int(row["Total"]) != frames or abs(row_p(row)-p) > 1e-10:
        raise ValueError(f"Simulator returned inconsistent batch: {temporary}")
    temporary.replace(tsv)
    return task


def summarize(tasks: list[tuple[float, float, int, int]], run_dir: Path) -> list[dict]:
    totals: dict[tuple[float, float], dict[str, int]] = {}
    for bias, p, seed, frames in tasks:
        row = read_row(run_dir / "batches" / f"eta_{bias:g}_p_{p:.5f}_seed_{seed}.tsv")
        item = totals.setdefault((bias, p), {"frames": 0, "logical_failures": 0,
                                             "x_syndrome_failures": 0,
                                             "z_syndrome_failures": 0,
                                             "undetected_logical_failures": 0})
        for target, source in (("frames", "Total"), ("logical_failures", "Failures"),
                               ("x_syndrome_failures", "XSyndromeFailures"),
                               ("z_syndrome_failures", "ZSyndromeFailures"),
                               ("undetected_logical_failures", "UndetectedLogicalFailures")):
            item[target] += int(row[source])
    result = []
    for (bias, p), counts in sorted(totals.items()):
        n = counts["frames"]
        k = counts["logical_failures"]
        low, high = wilson95(k, n)
        result.append({"eta_z_over_x": bias, "p_total": p,
                       "p_x": p/(bias+2), "p_y": p/(bias+2),
                       "p_z": bias*p/(bias+2), **counts,
                       "logical_fer": k/n,
                       "wilson_95_lower": low, "wilson_95_upper": high,
                       "zero_fail_one_sided_95_upper":
                           1 - 0.05**(1/n) if k == 0 else ""})
    csv_path = run_dir / "source_data.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    return result


def plot(result: list[dict], run_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharex=True)
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
    markers = ("o", "s", "D", "^", "P")
    biases = sorted({r["eta_z_over_x"] for r in result})
    for index, bias in enumerate(biases):
        subset = [r for r in result if r["eta_z_over_x"] == bias]
        c, m = colors[index % len(colors)], markers[index % len(markers)]
        positive = [r for r in subset if r["logical_failures"]]
        if positive:
            xs = [r["p_total"] for r in positive]
            ys = [r["logical_fer"] for r in positive]
            axes[0].errorbar(xs, ys,
                yerr=[[r["logical_fer"]-r["wilson_95_lower"] for r in positive],
                      [r["wilson_95_upper"]-r["logical_fer"] for r in positive]],
                color=c, marker=m, linestyle="-", capsize=2.5, label=fr"$\eta={bias:g}$")
        else:
            axes[0].plot([], [], color=c, marker=m, label=fr"$\eta={bias:g}$")
        zeros = [r for r in subset if not r["logical_failures"]]
        if zeros:
            axes[0].scatter([r["p_total"] for r in zeros],
                            [r["zero_fail_one_sided_95_upper"] for r in zeros],
                            marker="v", facecolors="none", edgecolors=c, s=55)
        axes[1].plot([r["p_total"] for r in subset],
                     [r["z_syndrome_failures"]/r["frames"] for r in subset],
                     color=c, marker=m, label=fr"$\eta={bias:g}$")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Quantum logical FER")
    axes[1].set_ylabel("Z-syndrome failure fraction")
    for ax in axes:
        ax.set_xlabel("Total physical Pauli error probability $p$")
        ax.grid(True, alpha=0.22, which="both")
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    for extension in ("png", "pdf", "svg"):
        fig.savefig(run_dir / f"fer_vs_p_z_bias.{extension}", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--biases", nargs="+", type=float, default=[1, 3, 10])
    parser.add_argument("--p-values", nargs="+", type=float,
                        default=[0.060, 0.065, 0.070, 0.075])
    parser.add_argument("--frames-per-point", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed-start", type=int, default=20271000)
    parser.add_argument("--matrix-dir", type=Path, default=DEFAULT_MATRICES)
    parser.add_argument("--build-dir", type=Path, default=Path("/tmp/gf256_cascade_build"))
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    args.matrix_dir = args.matrix_dir.resolve()
    args.run_dir = args.run_dir.resolve()
    args.build_dir = args.build_dir.resolve()
    if (args.frames_per_point < 1 or args.batch_size < 1 or
        args.frames_per_point % args.batch_size or args.workers < 1 or
        args.workers > (os.cpu_count() or 1)):
        parser.error("Require positive frames divisible by batch size and valid worker count")
    if (any(not math.isfinite(x) or x <= 0 for x in args.biases) or
        any(not math.isfinite(x) or not 0 < x < 1 for x in args.p_values) or
        len(set(args.biases)) != len(args.biases) or
        len(set(args.p_values)) != len(args.p_values)):
        parser.error("Biases must be distinct positive finite values; p must be distinct in (0,1)")
    if any(not (args.matrix_dir / "build" / name).is_file() for name in REQUIRED_MATRICES):
        parser.error(f"Missing fixed path/n2=12 GF(256) matrix files under {args.matrix_dir / 'build'}")
    subprocess.run(["cmake", "-S", str(ROOT), "-B", str(args.build_dir),
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run(["cmake", "--build", str(args.build_dir), "--target",
                    "gf256_cascade_simulate", "-j", "4"], check=True)
    executable = args.build_dir / "gf256_cascade_simulate"
    config = {"matrix_dir": str(args.matrix_dir),
              "matrix_sha256": {name: sha256(args.matrix_dir / "build" / name)
                                for name in REQUIRED_MATRICES},
              "simulator_sha256": sha256(executable), "mode": MODE,
              "bias_definition": "eta=p_Z/p_X, p_X=p_Y=p/(eta+2)",
              "biases": args.biases, "p_values": args.p_values,
              "frames_per_point": args.frames_per_point,
              "batch_size": args.batch_size, "seed_start": args.seed_start}
    args.run_dir.mkdir(parents=True, exist_ok=True)
    config_path = args.run_dir / "run_config.json"
    if config_path.exists() and json.loads(config_path.read_text(encoding="utf-8")) != config:
        parser.error(f"Existing run has different code or parameters; choose another --run-dir: {config_path}")
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (args.run_dir / "batches").mkdir(exist_ok=True)
    batches = args.frames_per_point // args.batch_size
    tasks = [(bias, p, args.seed_start + i*batches + batch, args.batch_size)
             for i, (bias, p) in enumerate((b, q) for b in args.biases
                                           for q in args.p_values)
             for batch in range(batches)]
    print(f"{len(args.biases)} biases x {len(args.p_values)} p values x "
          f"{args.frames_per_point} frames = {sum(t[3] for t in tasks)} frames", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_batch, task, args, executable, args.run_dir): task
                   for task in tasks}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            bias, p, seed, frames = future.result()
            print(f"[{done}/{len(tasks)}] eta={bias:g} p={p:.5f} "
                  f"seed={seed} frames={frames}", flush=True)
    result = summarize(tasks, args.run_dir)
    plot(result, args.run_dir)
    (args.run_dir / "figure_notes.md").write_text(
        "# GF(256) n2=12 path Z-biased Pauli FER\n\n"
        "Source: `source_data.csv`; raw batches in `batches/`; fixed inputs and "
        "hashes in `run_config.json`.\n\n"
        "Pauli noise: p_X=p_Y=p/(eta+2), p_Z=eta*p/(eta+2); eta=p_Z/p_X. "
        "Thus eta=1 is depolarizing. This is not Tuckett et al.'s Y-bias eta.\n\n"
        "Panel A: quantum logical FER (Wilson 95% intervals for positive counts). "
        "Open downward triangles mark the exact one-sided 95% upper bound "
        "for zero failures, not measured FER. Panel B: Z-syndrome failure fraction. "
        "This plot uses one fixed matrix and is not a threshold estimate.\n",
        encoding="utf-8")
    print(f"Wrote {args.run_dir / 'source_data.csv'} and "
          f"{args.run_dir / 'fer_vs_p_z_bias.pdf'}", flush=True)


if __name__ == "__main__":
    main()
