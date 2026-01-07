#!/bin/bash

# 仿真参数
NB_MONTE=10000 # 增大仿真帧数，提高精度（建议10000帧）
NB_ITER=100    # 增加迭代次数（GF(64)需更多迭代）
MATRIX_FILE="../../build/H_gamma.txt"
OFFSET=0.5
NB_OPER=64
N_VC=2
N_CV=4

# 待测试的Eb/No列表
SNR_LIST=(8 10 12 14 16 18 20)

# 批量运行仿真
for snr in "${SNR_LIST[@]}"; do
    echo "=== 开始仿真 Eb/No = $snr dB ==="
    ./essai $NB_MONTE $NB_ITER $MATRIX_FILE $snr $N_VC $N_CV $OFFSET $NB_OPER
    # 等待仿真完成，避免资源竞争
    sleep 2
done

echo "=== 所有SNR点仿真完成 ==="
