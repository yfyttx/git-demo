import argparse
import os
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MATRIX_DIR = os.path.join(SCRIPT_DIR, "matrices")
DECODER = os.path.join(SCRIPT_DIR, "build-win", "simulate_qec.exe")

# Quick verification config: use the default J=2, L=6, P=7 matrix set.
CONFIGS = [
    ("N=8400 (J=2, L=6, P=7, SPC length=25)", 2, 6, 7, 25),
]

ERROR_RATES = [0.02, 0.015, 0.012, 0.01, 0.008, 0.006, 0.004]


def legacy_suffix(suffix):
    parts = suffix.split("_")
    p_part = next((part for part in parts if part.startswith("P")), None)
    n_part = next((part for part in parts if part.startswith("n")), None)
    if p_part and n_part:
        return f"_{p_part}_{n_part}"
    return suffix


def matrix_exists(suffix):
    current = os.path.join(MATRIX_DIR, f"HX_Phys{suffix}.txt")
    if os.path.exists(current):
        return True
    return os.path.exists(os.path.join(MATRIX_DIR, f"HX_Phys{legacy_suffix(suffix)}.txt"))


def parse_args():
    parser = argparse.ArgumentParser(description="Run a quick QLDPC performance simulation.")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--max-iters", type=int, default=50)
    parser.add_argument("--target-failures", type=int, default=100)
    parser.add_argument(
        "--fm",
        default=",".join(str(rate) for rate in ERROR_RATES),
        help="Comma-separated f_m values.",
    )
    parser.add_argument(
        "--decoder-mode",
        choices=["split"],
        default="split",
        help="Decoder mode to pass to simulate_qec.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    error_rates = [float(item) for item in args.fm.split(",") if item.strip()]
    plt.figure(figsize=(10, 7))

    for label, j_val, l_val, p_val, spc_length in CONFIGS:
        print(f"\n>>> Running quick QEC simulation for {label}")

        suffix = f"_J{j_val}_L{l_val}_P{p_val}_spc{spc_length}_GF256"
        if not matrix_exists(suffix):
            print(f"  [Skip] Missing matrix file for suffix: {suffix}")
            continue

        cmd = [
            DECODER,
            str(args.trials),
            str(args.max_iters),
            str(args.target_failures),
            suffix,
        ]
        cmd.extend(str(rate) for rate in error_rates)
        print(f"  Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, cwd=SCRIPT_DIR, capture_output=True, text=True, check=False)
            output = result.stdout
            print(output)

            fer_x = []
            fer_y = []
            zero_failure_x = []
            zero_failure_upper = []
            for line in output.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) >= 8:
                    try:
                        fm = float(parts[0])
                        fer = float(parts[1])
                        upper95 = float(parts[7])
                    except ValueError:
                        continue
                    if fer == 0.0 and upper95 > 0.0:
                        zero_failure_x.append(fm)
                        zero_failure_upper.append(upper95)
                    else:
                        fer_x.append(fm)
                        fer_y.append(fer)

            if fer_y or zero_failure_x:
                if fer_y:
                    plt.plot(fer_x, fer_y, marker="s", label=label)
                if zero_failure_x:
                    plt.scatter(
                        zero_failure_x,
                        zero_failure_upper,
                        marker="v",
                        facecolors="none",
                        edgecolors="tab:blue",
                        label="0 failures: 95% upper bound",
                    )
            else:
                print("  [Warning] Could not parse FER data")

        except Exception as exc:
            print(f"  [Error] Simulation failed: {exc}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Physical Error Rate (p)")
    plt.ylabel("Total Quantum FER")
    plt.title("Quick QLDPC Performance Test")
    plt.grid(True, which="both", ls="-", alpha=0.3)
    plt.legend()
    plt.gca().invert_xaxis()

    output_png = os.path.join(SCRIPT_DIR, "quantum_final_comparison.png")
    plt.savefig(output_png)
    print(f"\n[Done] Saved comparison plot to {output_png}")


if __name__ == "__main__":
    main()
