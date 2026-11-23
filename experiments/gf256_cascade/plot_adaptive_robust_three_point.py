#!/usr/bin/env python3
"""Aggregate paired adaptive-vs-robust runs and plot the three-point comparison.

The two non-overlapping p_D=0.07125 runs are pooled because their matrix,
executable, decoder modes, frame count and channel protocol match. All figures
are regenerated from each run's summary.csv/config.json files.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "experiments/gf256_cascade/results"
INPUT_RUNS = [
    RESULTS / "paired_adaptive_pD_0p0705_300",
    RESULTS / "paired_adaptive_vs_robust_fresh300",
    RESULTS / "paired_adaptive_vs_robust_holdout300",
    RESULTS / "paired_adaptive_pD_0p072_300",
]
OUTPUT = RESULTS / "adaptive_robust_three_point_curve"
MODES = ("conditioned_robust", "conditioned_robust_adaptive")
MODE_LABELS = {
    "conditioned_robust": "Robust",
    "conditioned_robust_adaptive": "Adaptive",
}
Z95 = 1.959963984540054


def wilson(k: int, n: int) -> tuple[float, float]:
    p = k / n
    z2 = Z95 * Z95
    den = 1 + z2 / n
    center = (p + z2 / (2 * n)) / den
    half = Z95 * np.sqrt((p * (1 - p) + z2 / (4 * n)) / n) / den
    return float(center - half), float(center + half)


def load_and_validate() -> tuple[dict, dict]:
    pooled: dict[tuple[float, str], dict] = defaultdict(lambda: {
        "frames": 0, "logical_failures": 0, "z_syndrome_failures": 0,
        "undetected_logical_failures": 0, "cpu_seconds": 0.0,
        "seed_ranges": [], "run_dirs": [], "matrix_hash": None,
        "executable_hash": None,
    })
    paired: dict[float, dict] = defaultdict(lambda: {"rescued": 0, "regressed": 0})
    expected_matrix = None
    expected_executable = None

    for run_dir in INPUT_RUNS:
        config_path = run_dir / "config.json"
        summary_path = run_dir / "summary.csv"
        pairwise_path = run_dir / "pairwise.csv"
        if not all(p.exists() for p in (config_path, summary_path, pairwise_path)):
            raise FileNotFoundError(f"Incomplete run directory: {run_dir}")
        config = json.loads(config_path.read_text())
        if config["modes"] != list(MODES):
            raise ValueError(f"Unexpected decoder modes in {config_path}")
        matrix_hash = config["matrix_sha256"]
        executable_hash = config["executable_sha256"]
        if expected_matrix is None:
            expected_matrix = matrix_hash
            expected_executable = executable_hash
        if matrix_hash != expected_matrix or executable_hash != expected_executable:
            raise ValueError(f"Incompatible matrix/executable in {config_path}")
        if config["inner_length"] != 12:
            raise ValueError(f"Unexpected inner length in {config_path}")

        with summary_path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        if {row["mode"] for row in rows} != set(MODES):
            raise ValueError(f"Missing decoder mode in {summary_path}")
        p_depolar = float(rows[0]["p_D"])
        seed_start = int(config["seed_start"])
        seed_end = seed_start + int(config["frames_per_point"]) // int(config["batch_size"]) - 1
        for row in rows:
            key = (float(row["p_D"]), row["mode"])
            item = pooled[key]
            item["frames"] += int(row["frames"])
            item["logical_failures"] += int(row["logical_failures"])
            item["z_syndrome_failures"] += int(row["z_syndrome_failures"])
            item["undetected_logical_failures"] += int(row["undetected_logical_failures"])
            item["cpu_seconds"] += float(row["cpu_seconds"])
            item["seed_ranges"].append(f"{seed_start}-{seed_end}")
            item["run_dirs"].append(str(run_dir.relative_to(ROOT)))
            item["matrix_hash"] = matrix_hash["HX_Cascade_GF256.txt"]
            item["executable_hash"] = executable_hash

        with pairwise_path.open(newline="") as stream:
            pair_rows = list(csv.DictReader(stream))
        if len(pair_rows) != 1:
            raise ValueError(f"Expected one paired comparison in {pairwise_path}")
        pair = pair_rows[0]
        paired[p_depolar]["rescued"] += int(pair["a_fail_b_success"])
        paired[p_depolar]["regressed"] += int(pair["a_success_b_fail"])

    source_rows = []
    for (p_depolar, mode), item in sorted(pooled.items()):
        n = item["frames"]
        failures = item["logical_failures"]
        lower, upper = wilson(failures, n)
        source_rows.append({
            "p_D": p_depolar,
            "f_m": 2 * p_depolar / 3,
            "mode": mode,
            "frames": n,
            "logical_failures": failures,
            "FER": failures / n,
            "wilson_95_lower": lower,
            "wilson_95_upper": upper,
            "z_syndrome_failures": item["z_syndrome_failures"],
            "z_syndrome_failure_rate": item["z_syndrome_failures"] / n,
            "undetected_logical_failures": item["undetected_logical_failures"],
            "undetected_logical_failure_rate": item["undetected_logical_failures"] / n,
            "cpu_seconds_per_frame": item["cpu_seconds"] / n,
            "seed_ranges": ";".join(item["seed_ranges"]),
            "run_dirs": ";".join(item["run_dirs"]),
            "matrix_HX_sha256": item["matrix_hash"],
            "executable_sha256": item["executable_hash"],
        })
    return {"rows": source_rows, "paired": paired}, {
        "matrix_hashes": expected_matrix,
        "executable_hash": expected_executable,
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data, provenance = load_and_validate()
    source_csv = OUTPUT / "source_data.csv"
    with source_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])

    colors = {"conditioned_robust": "#5B6573", "conditioned_robust_adaptive": "#007C78"}
    markers = {"conditioned_robust": "o", "conditioned_robust_adaptive": "s"}
    series_by_mode = {
        mode: sorted((r for r in data["rows"] if r["mode"] == mode), key=lambda r: r["p_D"])
        for mode in MODES
    }
    fig, axes = plt.subplots(1, 4, figsize=(8.0, 2.65), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.24, top=0.76, wspace=0.48)
    metrics = [
        ("FER", "Logical FER", True),
        ("z_syndrome_failure_rate", "Z syndrome\nfailure rate", False),
        ("undetected_logical_failure_rate", "Undetected logical\nfailure rate", False),
        ("cpu_seconds_per_frame", "CPU time\n(s/frame)", False),
    ]
    for ax, (metric, ylabel, with_ci) in zip(axes, metrics):
        for mode in MODES:
            rows = series_by_mode[mode]
            x = [r["f_m"] for r in rows]
            y = [r[metric] for r in rows]
            kwargs = dict(color=colors[mode], marker=markers[mode], linewidth=1.6,
                          markersize=4, label=MODE_LABELS[mode])
            if with_ci:
                low = [r[metric] - r["wilson_95_lower"] for r in rows]
                high = [r["wilson_95_upper"] - r[metric] for r in rows]
                ax.errorbar(x, y, yerr=[low, high], capsize=2, **kwargs)
            else:
                ax.plot(x, y, **kwargs)
        ax.set_xlabel(r"$f_m$")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_xticks([0.047, 0.0475, 0.048], labels=["0.0470", "0.0475", "0.0480"])
        ax.tick_params(labelsize=7, direction="in")
        ax.grid(True, which="major", color="#D8DDE3", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, loc="upper left")
    fig.supxlabel(r"Marginal flip probability $f_m$", y=0.035, fontsize=8)
    for idx, ax in enumerate(axes):
        ax.text(-0.19, 1.08, f"({chr(97 + idx)})", transform=ax.transAxes,
                fontsize=8, fontweight="bold", va="bottom")
    fig.savefig(OUTPUT / "figure.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT / "figure.svg", bbox_inches="tight")
    fig.savefig(OUTPUT / "figure.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    comparison = []
    for p_depolar in sorted(data["paired"]):
        pair = data["paired"][p_depolar]
        n_discordant = pair["rescued"] + pair["regressed"]
        p_value = (binomtest(pair["rescued"], n_discordant, 0.5).pvalue
                   if n_discordant else 1.0)
        comparison.append({"p_D": p_depolar, "f_m": 2*p_depolar/3,
                           **pair, "mcnemar_exact_two_sided_p": p_value})
    with (OUTPUT / "paired_tests.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)

    notes = f"""# Adaptive vs robust: three-point experiment figure

- Source status: `ACTUAL_RUN`.
- Claim supported: on the fixed GF(256), inner-length-12 path cascade matrix, the adaptive decoder has lower observed logical FER than the raw robust decoder at the three tested channel points; compute cost is also shown.
- Matrix: one fixed matrix. HX hash `{provenance['matrix_hashes']['HX_Cascade_GF256.txt']}`; HZ hash `{provenance['matrix_hashes']['HZ_Cascade_GF256.txt']}`.
- Executable SHA256: `{provenance['executable_hash']}`.
- Points: `p_D=0.0705, 0.07125, 0.072`, with `f_m=2p_D/3`; each point has 300 paired frames except `p_D=0.07125`, which pools two non-overlapping 300-frame runs (600 frames total).
- Aggregation: pooled failure/frame counts by mode; Wilson 95% binomial intervals for plotted FER; arithmetic total CPU seconds divided by frames; paired exact two-sided McNemar/binomial test from the paired rescue/regression counts. No run was excluded.
- Source files: `source_data.csv`, `paired_tests.csv`; generated from each run's `summary.csv`, `pairwise.csv`, and `config.json`.
- Generator: `experiments/gf256_cascade/plot_adaptive_robust_three_point.py`.
- Output: `figure.pdf`, `figure.svg`, `figure.png`.
- Placement plan: two-column-spanning exploratory figure, four meaningful panels in one row: (a) logical FER, (b) Z syndrome failure rate, (c) undetected logical failure rate, (d) CPU seconds per frame. Canvas width 8.0 in; inspect at the intended manuscript placement before treating as final.
- Caption draft: Logical frame error rate (FER), Z-syndrome failure rate, undetected logical failure rate, and CPU time for `conditioned_robust` and `conditioned_robust_adaptive` on a fixed GF(256), inner-length-12 path cascade matrix. The x-axis is the marginal flip probability `f_m`; each point uses paired noise frames (`n=300`, except `n=600` at `p_D=0.07125`). FER bars show Wilson 95% confidence intervals. CPU time is total measured decoder CPU seconds per frame. Results describe this matrix and tested channel points only.
- QA: generated by the source-driven script. Visual inspection pending.
"""
    (OUTPUT / "figure_notes.md").write_text(notes, encoding="utf-8")
    print(f"Wrote {OUTPUT / 'figure.pdf'}, source_data.csv, paired_tests.csv")
    for row in comparison:
        print(row)


if __name__ == "__main__":
    main()
