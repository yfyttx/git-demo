# Z 侧综合征失败：消息刷新与联合投影诊断（2026-10-26）

## 固定实例与结果

矩阵沿用 `../gf256_cascade_curve/run_20261025_001136/matrices`，即 GF(256)、
`H_X^cas = H_Z^TIT ⊗ [1|I_5]`、3492 个 GF(256) 符号（27,936 物理比特）。
退极化参数 `p_D=0.07125`，对应边缘翻转率 `f_m=0.0475`。
三个模式使用相同的随机种子 `20263200` 和前 10 帧。

| 模式 | 逻辑失败 | Z 综合征失败 | 综合征通过但逻辑失败 |
| --- | ---: | ---: | ---: |
| 原 `conditioned_robust` | 4/10 | 4/10 | 0/10 |
| 末轮外层消息刷新 `conditioned_robust_refresh` | 4/10 | 4/10 | 0/10 |
| 失败后五投影联合候选 `conditioned_robust_joint_projection` | 4/10 | 0/10 | 4/10 |
| 失败后五投影有限候选列表 `conditioned_robust_joint_list` | 4/10 | 0/10 | 4/10 |

刷新模式的 CPU 时间为 245.82 秒，原模式为 232.92 秒，联合模式为
276.42 秒，有限候选列表模式为 259.05 秒（均为 10 帧总和）。这些小批次
只用于定位机制，不能估计总体 FER 或声称门限改善。

诊断日志用 `bp_valid` 记录五个外层 BP 输出中通过校验的个数，
`raw_valid` 记录组合后硬判决通过的个数，`final_valid` 记录最终物理
Z 估计是否通过完整 `H_X^cas` 校验，`outer_refreshed` 表示是否刷新了末轮消息。
原模式失败的 trial 1、2、5、9 在初次、重试、备选路径中均为
`bp_valid=0, raw_valid=0, final_valid=0`；刷新后四帧仍失败。
因此末轮消息不同步是真实的实现差异，但不是这四帧失败的主因。

新联合候选模式只在原译码器 Z 综合征失败后运行：先用外层 OSD-0
生成五个综合征有效的投影，再按共享内层符号的精确局部似然评分，
尝试条件 OSD 候选更新：

\[
L(u_0,\ldots,u_4)=\sum_j\log\!\left[
 \sum_{t\in GF(256)}P_j(z_0=t)
 \prod_{a=0}^{4}P_j(z_{a+1}=t+u_a)
\right].
\]

四个失败帧最终都得到满足综合征的候选，但仍属于错误的逻辑类；
联合坐标更新次数为零。随后扩展为从有限 OSD 候选列表逐一用联合似然
评分，同批 10 帧的逻辑失败仍为 4/10，联合坐标更新仍为零。
在 trial 1 的重新运行中，日志确认列表模式实际评估了 85 个联合评分候选，
仍没有接受更新。
也就是说，这两个小候选集都没有找到比 OSD-0 更好的投影组合，
不能将其作为有效的 FER 改进。
单帧 `osd_first2.tsv` 也单独复现了“综合征修复而逻辑未修复”。

## 信息量核对与后续方向

在独立退极化噪声模型下、假定物理 X 错误完全已知，逐位枚举
6 个 X 位和 5 个投影位的所有状态，在 `p_D=0.07125` 得到：

- 单个投影的条件熵约 1617 bit，高于其 194 行 GF(256)
  综合征最多提供的 1552 bit。
- 五个投影的**联合**条件熵约 5803 bit，低于完整
  `H_X^cas` 的 970 行 GF(256) 综合征最多提供的 7760 bit。

这说明逐投影单独译码的信息预算紧张，而共享符号提供重要相关性。
熵计算是信道模型下的平均量，不是对有限码长 FER 的保证。

下一步应在失败帧上构造**多候选联合搜索**，避免每个投影只保留一个
OSD-0 解；候选选择只能使用综合征与信道先验，模拟真实错误只能用于
事后评价。同时以原译码器和已验证的自适应/加长迭代模式为配对基线，
报告逻辑 FER、Z 综合征失败、未检出逻辑失败和 CPU 时间。
如多候选联合搜索仍无法降低逻辑 FER，再比较增加外层冗余或改变内码
`H_2,G_2` 的独立矩阵构造，不应单纯调整 OSD 使综合征通过。

## 复现

从仓库根目录编译：

```bash
cmake -S . -B /tmp/gf256_cascade_refresh_build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/gf256_cascade_refresh_build --target gf256_cascade_simulate -j 4
```

配对刷新实验使用 `run_paired_decoder_study.py --modes conditioned_robust
conditioned_robust_refresh`，其中 `--seed-start 20263200 --frames-per-point 10
--batch-size 10 --p-depolar 0.07125`。设置环境变量
`GF256_CASCADE_PROJECTION_DIAG=1` 可输出投影诊断。完整参数及矩阵哈希
记录在 `config.json`；原始 TSV、日志和逐帧结果保存在本目录。

`osd_frame1.log` 是一次带 `target_trial` 的诊断运行留下的中断日志；
完成的 OSD 小批次是 `osd_first2.tsv` 和 `osd_first2.log`。
