"""Aggregate compatible simulation TSV files and report an exact one-sided bound."""

import csv
import math
import sys
from pathlib import Path


def binomial_cdf_at_most(k: int, n: int, p: float) -> float:
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0 if k < n else 1.0
    term = math.exp(n * math.log1p(-p))
    total = term
    for j in range(1, k + 1):
        term *= (n - j + 1) / j * p / (1 - p)
        total += term
    return total


def upper_confidence_bound(k: int, n: int, alpha: float = 0.05) -> float:
    if k == n:
        return 1.0
    low, high = 0.0, 1.0
    for _ in range(80):
        middle = (low + high) / 2
        if binomial_cdf_at_most(k, n, middle) > alpha:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def main(paths: list[str]) -> None:
    if not paths:
        raise SystemExit("Usage: summarize_results.py file1.tsv [file2.tsv ...]")
    totals = {key: 0 for key in (
        "Total", "Failures", "SyndromeFailures", "XSyndromeFailures",
        "ZSyndromeFailures", "UndetectedLogicalFailures", "ExactXFrames",
        "ExactRecoveryFrames", "XFallbackAttempts", "XFallbackSuccesses",
        "ZSoftFallbackAttempts", "ZSoftFallbackSuccesses",
        "RobustAttempts", "RobustSelections",
    )}
    depolarization_rate = None
    for name in paths:
        with Path(name).open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(rows) != 1:
            raise SystemExit(f"Expected one data row in {name}")
        row = rows[0]
        if depolarization_rate is None:
            depolarization_rate = row["p_D"]
        elif row["p_D"] != depolarization_rate:
            raise SystemExit("All rows must use the same p_D")
        for key in totals:
            totals[key] += int(row.get(key) or 0)
    n, k = totals["Total"], totals["Failures"]
    print(f"p_D={depolarization_rate} files={len(paths)}")
    print(f"FER={k}/{n}={k/n:.9g}")
    print(f"exact_one_sided_95_percent_upper={upper_confidence_bound(k, n):.9g}")
    for key, value in totals.items():
        if key not in ("Total", "Failures"):
            print(f"{key}={value}")


if __name__ == "__main__":
    main(sys.argv[1:])
