# Kasai 等人 2012 TIT 论文复现模式

论文：K. Kasai, M. Hagiwara, H. Imai, K. Sakaniwa, "Quantum Error Correction beyond the Bounded Distance Decoding Limit," IEEE Transactions on Information Theory 58(2), 1223–1230 (2012). [arXiv 全文](https://arxiv.org/html/1007.1778v2)。

`kasai_paper` 是独立运行入口。矩阵构造对应论文第 II 节：按 Theorem 1 生成两个 `(J=2,L,P)` 的 QC 二进制支撑；按式 (16) 在模 `2^p-1` 上求 `H_Gamma` 的非零系数；按式 (13) 求 `H_Delta`，并验证 `H_Gamma H_Delta^T=0`。物理二进制校验块分别是 `A(Gamma)` 和 `A(Delta)^T`。

译码对应论文第 III 节：对两个分量分别使用 `GF(2^p)` sum-product，按式 (18)–(19) 对行节点做 XOR 卷积，用 Walsh–Hadamard 变换计算；列节点乘信道先验和其余校验消息；每轮做最大后验符号判决，并检查综合征。没有使用项目的联合 X/Z 消息或级联后处理。每个物理位独立经过去极化信道，`p_D=3 f_m/2`，每个 X/Z 分量的边缘翻转率是 `f_m`。

```bash
cmake -S . -B /tmp/git-demo-kasai-paper-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/git-demo-kasai-paper-build --target kasai_paper -j2
/tmp/git-demo-kasai-paper-build/kasai_paper 1000 20261018 0.04 8 31 5 2 6
```

参数依次为 `帧数 随机种子 f_m 有限域指数p P sigma tau L [输出文件前缀]`。支持 `p=6,8,9,10`；默认 `20 20261018 0.02 8 31 5 2 6`。指定输出前缀后，会保存 `Gamma`、`Delta` 及二进制展开 `HC`、`HD` 四个 Matrix Market 文件。`Gamma/Delta` 文件中的整数是有限域元素的多项式表示，矩阵乘法仍应按 GF(2^p) 计算。程序还会验证二进制展开满足 `HC HD^T=0`（模 2）。

输出中的 `Gamma_block_error` / `Delta_block_error` 定义为该分量估计错误向量与真实错误向量不相同的帧比例；`*_syndrome_fail` 表示 50 轮后估计仍不满足综合征的比例。`actual_RQ` 使用生成矩阵的实际秩，`nominal_RQ` 使用论文的 `1-2J/L`。

论文公开了构造族和译码公式，没有给出每条仿真曲线的完整矩阵系数与随机种子。这里可以复现方法并比较同参数趋势；单次生成的矩阵和论文图中的具体矩阵不保证逐元素相同。尤其 `P=31` 和 50 轮是本项目可复现实验的设定，不能未经核对便把输出点当作论文图中同一个点。

## GF(256)、码长 4656 的直接译码实验（2026-11-17）

只使用 `H_Gamma`、`H_Delta`；`G2`、`H2` 和级联译码器均不参与。配置为 `J=2, L=6, P=97, sigma=35, tau=2, q=256`，矩阵种子 `20261019`，两个分量分别使用非二进制 sum-product BP，最多 50 轮。退极化信道满足 `p_D=1.5 f_m`。两个 GF(256) 矩阵分别是 `194×582`；二进制展开后的量子码物理码长为 4656，实际码率为 1/3。

| `f_m` | `p_D` | 帧数 | Gamma 分量块错误 | Delta 分量块错误 | Gamma / Delta 综合征失败 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0400 | 0.06000 | 1000 | 0/1000 | 1/1000 | 0 / 1 |
| 0.0450 | 0.06750 | 1000 | 8/1000 | 10/1000 | 8 / 9 |
| 0.0475 | 0.07125 | 1000 | 63/1000 | 66/1000 | 63 / 64 |

三个点的 Gamma 和 Delta Matrix Market 文件各自 SHA-256 完全相同，确认为同一矩阵扫描。结果及矩阵在 `experiments/gf256_cascade/results/kasai_direct_20261117_full/`。一键复现：

```bash
experiments/gf256_cascade/run_kasai_direct_20261117.sh
```

块错误按估计错误向量是否逐符号等于真实错误向量统计，匹配论文图 4 的**分量**口径；它不是稳定子等价的完整量子逻辑 FER。`0/1000` 只表示本次没有观察到错误。论文没有公开该曲线的完整矩阵系数和随机种子，且当前 50 轮上限尚未核实与图 4 完全相同；本实验复现了方法和同类参数，不声称逐点复现原图。
