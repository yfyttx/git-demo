# 生成 GF(256) 级联码 FER 曲线

在仓库根目录运行：

```bash
python3 experiments/gf256_cascade/run_gf256_cascade_sweep.py
```

脚本会用 Release 模式编译构造器和译码器，按固定种子重新构造矩阵、核对矩阵哈希，然后运行 `conditioned_robust`。已有的 `p_D=0.06`、0/3050 结果作为曲线起点；新仿真默认参数写在 `run_gf256_cascade_sweep.py` 顶部的 `SWEEP` 中（合计 1450 帧）。完成后自动汇总并生成图和源数据。

可以调整并行进程数和每个批次的帧数：

```bash
python3 experiments/gf256_cascade/run_gf256_cascade_sweep.py --workers 4 --batch-size 25
```

默认最多使用 8 个并行 CPU 进程。每个进程本身单线程。机器比较忙时可降低 `--workers`；CPU 核心较多、内存充足时可调高，但当前译码程序没有 GPU/CUDA 实现。

如果某次仿真已完成、只想重新绘图：

```bash
python3 experiments/gf256_cascade/plot_gf256_cascade_sweep.py \
  --run-dir experiments/gf256_cascade/results/gf256_cascade_curve/run_时间戳
```

绘图依赖 `matplotlib` 和 `scipy`。缺少时安装：

```bash
python3 -m pip install matplotlib scipy
```

每次运行会新建 `experiments/gf256_cascade/results/gf256_cascade_curve/run_时间戳/`，里面包括：

- `runs/`：各个 (p_D) 点的 TSV 原始结果和逐批日志；
- `matrices/`、`construct.log`：本次使用的矩阵；
- `run_config.json`：矩阵种子、译码方式、批量大小和帧数配置；
- `source_data.csv`：画图实际使用的汇总计数、FER 和置信区间；
- `fer_vs_fm.png`、`fer_vs_fm.pdf`、`fer_vs_fm.svg`：图像文件；
- `figure_notes.md`：图中指标和置信区间的说明。

横轴按去极化模型换算为 (f_m=2p_D/3)。纵轴是量子逻辑帧错误率。零失败点不会画成 FER=0，而是用倒三角标记其单侧 95% 精确置信上限；非零点显示双侧 95% Clopper–Pearson 区间。低错误率点要想把 FER 上限压得更低，需要增加该点的帧数，可直接修改脚本顶部 `SWEEP`，例如把 `(0.06375, 500)` 改为 `(0.06375, 5000)`。

如果按 Ctrl+C，已完成批次的 TSV 会保留在本次目录中，未完成批次会带 `.partial` 后缀。绘图脚本只汇总完整的 `.tsv` 文件；再次运行主脚本会创建一个新的 run 目录。
