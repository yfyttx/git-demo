#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="${KASAI_BUILD_DIR:-/tmp/git-demo-kasai-direct-build}"
run_dir="${KASAI_RUN_DIR:-$repo_root/experiments/gf256_cascade/results/kasai_direct_20261117_full}"
frames="${KASAI_FRAMES:-1000}"
seed="${KASAI_SEED:-20261019}"
mkdir -p "$run_dir"
cmake -S "$repo_root" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release > "$run_dir/configure.log"
cmake --build "$build_dir" --target kasai_paper -j4 > "$run_dir/build.log" 2>&1
for fm in 0.0400 0.0450 0.0475; do
  echo "Running f_m=$fm, frames=$frames, seed=$seed"
  "$build_dir/kasai_paper" "$frames" "$seed" "$fm" 8 97 35 2 6 \
    "$run_dir/fm_${fm}" > "$run_dir/fm_${fm}.log" 2>&1
  cat "$run_dir/fm_${fm}.log"
done
python3 "$repo_root/experiments/gf256_cascade/summarize_kasai_direct.py" "$run_dir"
