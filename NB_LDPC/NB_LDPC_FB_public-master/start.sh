#!/bin/bash
binary="./essai"
matrix="../../build/HX_Safe_GF64.txt"
common_args="50 $matrix" # 迭代次数 50, 矩阵

for SNR in 5.0 5.1 5.2 5.3; do
    echo "Running LOW SNR $SNR (100 frames)..."
    $binary 100 $common_args $SNR 20 20 0.5 64
done

echo "Running MID SNR 5.4 (1000 frames)..."
$binary 1000 $common_args 5.4 20 20 0.5 64

for SNR in 5.5 5.6 5.7 5.8; do
    echo "Running HIGH SNR $SNR (5000 frames)..."
    $binary 5000 $common_args $SNR 20 20 0.5 64
done
