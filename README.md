# 量子 LDPC 码构造与译码

本仓库包含 Kasai 等人的非二进制 CSS 码复现实验，以及基于同类外码的级联构造和软判决译码。仿真在理想综合征条件下运行。分量块错误率、物理 WER 和量子逻辑 FER 分别记录。

## 构建与运行

在仓库根目录运行。C++ 程序需要 CMake、支持 C++17 的编译器和 Eigen3。主要 Python 实验需要 NumPy、SciPy，绘图需要 Matplotlib。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
```

| 入口 | 用途 |
| --- | --- |
| `kasai_paper` | 单独生成并译码 `H_Gamma`、`H_Delta`，统计论文图 4 口径的分量块错误率 |
| `gf256_cascade_construct`、`gf256_cascade_simulate` | 同域 Kronecker 积级联，支持多种内码和译码方案 |
| `experiments/gf256_cascade/construct_companion_cascade.py`、`run_companion_staged_xz.py` | 用 GF(256) 外码系数对应的二进制乘法块与 `[9,8,2]` 内码构造校验矩阵，分阶段估计 X/Z 错误 |
| `construct`、`simulate_qec_gf64` | GF(64) 构造与仿真 |
| `ldpc_constructor`、`simulate_qec` | 原有 GF(256) 矩阵格式与批量扫描入口，供 `generate_matrices.ps1` 和 `runs/run_sweep.py` 调用 |

直接外码、条件先验和同外码级联对照的快速检查：

```bash
./build/kasai_paper 1 20261019 0.04 8 97 35 2 6

python3 experiments/gf256_cascade/run_companion_staged_xz.py \
  --frames 1 --p-depolar 0.07125 --seed 20261122 --max-iters 100 \
  --z-prior conditioned

python3 experiments/gf256_cascade/run_inner_code_comparison.py \
  --p-total 0.12 --alpha 100 --frames-per-matrix 1 --batch-size 1 \
  --seed-start 2026112300 --workers 2 --output-dir /tmp/ldpc-cascade-smoke
```

已保存的 5,000 帧对照使用每套矩阵 5,000 个配对噪声帧。要重跑该实验，在第三条命令中指定 `--frames-per-matrix 5000 --batch-size 100` 和新的输出目录。

## 实验汇总

- **直接外码**。GF(256)、`P=97`、物理码长 4,656、码率 1/3。在 `f_m=0.04` 运行 1,000 帧，Gamma 和 Delta 的分量块错误分别为 0 次和 1 次。[复现说明](KASAI_PAPER_REPRO.md)
- **Z 侧条件先验**。在 `p_D=0.07125` 的单矩阵 100 帧配对试验中，用实际的 X 侧译码估计构造 Z 先验，逻辑失败由 70 次降至 4 次。[实验设置与结果](experiments/gf256_cascade/COMPANION_CONDITIONED_PRIOR_20261118.md)
- **九位内码与同外码基线**。在总 Pauli 错误概率 `p=0.12`、偏置 `alpha=100`、相同局部 Hadamard 和相同码率下，两个外码实例的逻辑失败分别从 53/5000 降到 37/5000、从 49/5000 降到 34/5000。该对照覆盖这一固定信道点。[配对汇总数据](experiments/gf256_cascade/results/inner_code_comparison_alpha100_p0.12_confirm5000/summary.json)

更多构造和译码选项见[级联实验说明](experiments/gf256_cascade/README.md)，当前进展与未完成问题见[研究任务记录](TASKS.md)。仓库跟踪复现所需的输入矩阵和选定的实验汇总。运行生成的逐帧数据、日志和图像保存在各实验的输出目录。
