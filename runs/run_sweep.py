import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path


def run_command(command, cwd, capture_to=None):
    with subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as proc:
        lines = []
        for line in proc.stdout:
            print(line, end="")
            lines.append(line)
        code = proc.wait()

    if code != 0:
        raise subprocess.CalledProcessError(code, command)

    if capture_to is not None:
        capture_to.write_text("".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Run QEC simulation sweeps.")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--max-iters", type=int, default=50)
    parser.add_argument("--target-failures", type=int, default=100)
    parser.add_argument(
        "--fm",
        default="0.010,0.012,0.014,0.016,0.018,0.020",
        help="Comma-separated f_m values.",
    )
    parser.add_argument("--construct", action="store_true", help="Regenerate matrices first.")
    parser.add_argument("--no-build", action="store_true", help="Skip cmake build step.")
    parser.add_argument(
        "--suffix",
        default="_J2_L6_P7_spc25_GF256",
        help="Matrix filename suffix to simulate.",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    build_dir = project_root / "build-win"
    simulator = build_dir / "simulate_qec.exe"
    constructor = build_dir / "ldpc_constructor.exe"
    run_dir = project_root / "runs"
    run_dir.mkdir(exist_ok=True)

    if not args.no_build:
        print("Building simulate_qec...")
        subprocess.run(
            ["cmake", "--build", str(build_dir), "--target", "simulate_qec"],
            cwd=project_root,
            check=True,
        )
        if args.construct:
            print("Building ldpc_constructor...")
            subprocess.run(
                ["cmake", "--build", str(build_dir), "--target", "ldpc_constructor"],
                cwd=project_root,
                check=True,
            )

    if args.construct:
        if not constructor.exists():
            raise FileNotFoundError(f"Constructor not found: {constructor}")
        print("Regenerating matrices...")
        subprocess.run([str(constructor)], cwd=project_root, check=True)

    if not simulator.exists():
        raise FileNotFoundError(f"Simulator not found: {simulator}")

    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = run_dir / f"sim_{timestamp}.tsv"

    print()
    print("Running simulation")
    print(f"  Trials         = {args.trials}")
    print(f"  MaxIters       = {args.max_iters}")
    print(f"  TargetFailures = {args.target_failures}")
    print(f"  Matrix suffix   = {args.suffix}")
    print(f"  f_m list       = {args.fm}")
    print(f"  Output         = {output_file}")
    print()

    start = time.perf_counter()
    command = [
        str(simulator),
        str(args.trials),
        str(args.max_iters),
        str(args.target_failures),
        args.suffix,
        args.fm,
    ]
    run_command(command, cwd=project_root, capture_to=output_file)
    elapsed = time.perf_counter() - start

    print()
    print(f"Done in {elapsed:.1f} seconds")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
