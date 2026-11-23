#!/usr/bin/env python3
"""Aggregate saved robust decoder TSVs and draw logical FER versus f_m."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import beta
except ImportError as exc:
    raise SystemExit("Plotting requires matplotlib and scipy. Install with: "
                     "python3 -m pip install matplotlib scipy") from exc


ROOT = Path(__file__).resolve().parents[2]
HISTORICAL = ROOT / "experiments/gf256_cascade/results"


def read_tsv(path: Path) -> dict[str, str] | None:
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(rows) != 1 or not {"p_D", "Total", "Failures"} <= rows[0].keys():
            return None
        return rows[0]
    except (OSError, csv.Error):
        return None


def collect(run_dir: Path) -> dict[float, dict[str, int]]:
    paths = list(HISTORICAL.glob("robust_202610*_50.tsv"))
    paths += list(HISTORICAL.glob("robust_202610*_350.tsv"))
    paths += list((run_dir / "runs").glob("pD_*/*.tsv"))
    totals: dict[float, dict[str, int]] = {}
    for path in sorted(set(paths)):
        row = read_tsv(path)
        if row is None:
            continue
        try:
            p = round(float(row["p_D"]), 8)
            item = totals.setdefault(p, {"frames": 0, "failures": 0,
                                         "syndrome_failures": 0,
                                         "exact_recovery_frames": 0})
            item["frames"] += int(row["Total"])
            item["failures"] += int(row["Failures"])
            item["syndrome_failures"] += int(row.get("SyndromeFailures", 0))
            item["exact_recovery_frames"] += int(row.get("ExactRecoveryFrames", 0))
        except (TypeError, ValueError):
            continue
    return totals


def interval(failures: int, frames: int) -> tuple[float, float, float]:
    estimate = failures / frames
    if failures == 0:
        # Exact one-sided 95% Clopper–Pearson upper bound, not a zero FER claim.
        return estimate, 0.0, float(beta.ppf(0.95, 1, frames))
    lower = 0.0 if failures == 0 else float(beta.ppf(0.025, failures, frames - failures + 1))
    upper = 1.0 if failures == frames else float(beta.ppf(0.975, failures + 1, frames - failures))
    return estimate, lower, upper


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True,
                        help="run folder created by run_gf256_cascade_sweep.py")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        parser.error(f"run directory does not exist: {run_dir}")
    data = collect(run_dir)
    if not data:
        raise SystemExit(f"No complete TSV results found under {run_dir}")

    rows = []
    for p, counts in sorted(data.items()):
        n, k = counts["frames"], counts["failures"]
        estimate, low, high = interval(k, n)
        rows.append({
            "p_D": p,
            "f_m": 2 * p / 3,
            "frames": n,
            "logical_failures": k,
            "logical_FER": estimate,
            "CI95_lower": low,
            "CI95_upper": high,
            "interval_note": "one-sided exact 95% upper bound; zero observed failures"
                            if k == 0 else "two-sided exact 95% Clopper-Pearson interval",
            "syndrome_failures": counts["syndrome_failures"],
            "exact_recovery_frames": counts["exact_recovery_frames"],
        })

    csv_path = run_dir / "source_data.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    xs = [row["f_m"] for row in rows]
    # Plot zero-count points at their exact 95% upper limit and mark them with
    # downward triangles; never place zero on a logarithmic axis.
    ys = [row["CI95_upper"] if row["logical_failures"] == 0
          else row["logical_FER"] for row in rows]
    ax.plot(xs, ys, color="#a33b32", linewidth=1.5, zorder=1)
    for row, x, y in zip(rows, xs, ys):
        k, n = row["logical_failures"], row["frames"]
        if k == 0:
            ax.scatter([x], [y], marker="v", facecolors="white", edgecolors="#a33b32",
                       linewidths=1.6, s=58, zorder=3)
            ax.annotate("0 failures", (x, y), xytext=(4, 6),
                        textcoords="offset points", fontsize=8, color="#555555")
        else:
            ax.errorbar(x, y, yerr=[[y - row["CI95_lower"]],
                                    [row["CI95_upper"] - y]],
                        fmt="o", color="#a33b32", ecolor="#a33b32",
                        capsize=3, markersize=5, zorder=3)

    ax.set_yscale("log")
    ax.set_xlabel(r"Marginal flip probability $f_m=2p_D/3$")
    ax.set_ylabel("Quantum logical frame error rate")
    ax.set_title("GF(256) concatenated CSS code — conditioned robust BP\n"
                 "P=97; 3,492 GF(256) symbols (27,936 physical qubits); "
                 "p_D=0.06 baseline: 0/3,050",
                 fontsize=13)
    ax.grid(True, which="both", linestyle=":", linewidth=0.7, alpha=0.65)
    min_x, max_x = min(xs), max(xs)
    pad = max((max_x - min_x) * 0.08, 0.0005)
    ax.set_xlim(min_x - pad, max_x + pad)
    min_y = min(row["CI95_upper"] for row in rows if row["CI95_upper"] > 0)
    ax.set_ylim(max(1e-6, min_y * 0.55), 1.0)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(run_dir / f"fer_vs_fm.{suffix}", dpi=220)
    plt.close(fig)

    (run_dir / "figure_notes.md").write_text(
        "# Figure notes\n\n"
        "The plotted outcome is quantum logical frame error rate under the "
        "`conditioned_robust` decoder. The horizontal coordinate uses "
        "`f_m = 2 p_D / 3` for the depolarizing channel. Positive failure counts "
        "show exact two-sided 95% Clopper–Pearson intervals. Downward triangles "
        "mark zero observed failures and are positioned at the exact one-sided "
        "95% upper bound; they do not mean the true FER is zero. The saved "
        "`source_data.csv` contains the aggregated counts and intervals.\n",
        encoding="utf-8")
    print(f"Wrote {csv_path}, fer_vs_fm.png/.pdf/.svg, and figure_notes.md")


if __name__ == "__main__":
    main()
