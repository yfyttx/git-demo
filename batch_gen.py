import os
import subprocess


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MATRIX_DIR = os.path.join(SCRIPT_DIR, "matrices")
EXECUTABLE = os.path.join(SCRIPT_DIR, "build-win", "ldpc_constructor.exe")

# Parameter sets: (J, L, P, m2, n2)
CONFIGS = [
    (2, 6, 7, 5, 30),
]


def main():
    os.makedirs(MATRIX_DIR, exist_ok=True)

    if not os.path.exists(EXECUTABLE):
        print(f"Error: executable not found: {EXECUTABLE}")
        print("Please build ldpc_constructor.exe in build-win first.")
        return 1

    for J, L, P, m2, n2 in CONFIGS:
        cmd = [EXECUTABLE, str(J), str(L), str(P), str(m2), str(n2)]
        print(f"\n[Running] {' '.join(cmd)}")
        subprocess.run(cmd, cwd=SCRIPT_DIR, check=True)

    print("\n--- All matrix sets generated. ---")
    print(f"Matrix files are in: {MATRIX_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
