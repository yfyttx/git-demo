#!/usr/bin/env bash
# Fixed baseline versus bounded physical-bit guessing. New run directory per build.
set -euo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$project_root"
study_build="$project_root/experiments/gf256_cascade/.build-guided"
study_run="$project_root/experiments/gf256_cascade/results/guided_study_$(date +%Y%m%d_%H%M%S)"
cmake -S . -B "$study_build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$study_build" --target gf256_cascade_simulate -j4
python3 experiments/gf256_cascade/run_paired_decoder_study.py \
  --p-depolar 0.07125 --frames-per-point 100 --batch-size 5 \
  --seed-start 20281102 --workers 4 \
  --modes conditioned_robust_adaptive conditioned_robust_guess \
  --matrix-dir experiments/gf256_cascade/results/g2_basis_path_n12_300frames/path \
  --inner-length 12 --allow-custom-matrices \
  --build-dir "$study_build" --run-dir "$study_run" "$@"
