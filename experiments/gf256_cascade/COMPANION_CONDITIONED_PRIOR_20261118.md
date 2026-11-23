# 新伴随块构造：实际 X 估计条件先验配对实验

## 问题与冻结设置

验证在保持矩阵和 BP 不变时，由实际 X 侧估计产生的条件符号先验能否改善 Z 外层译码。固定100帧，未按中途结果调整停止条件。

- GF(256)、外层 P=97、两外层矩阵均194×582，矩阵种子20261019。
- 二进制内码 G2=[I8|1]，8×9；H2为1×9全一行。
- HX=1552×5238，HZ=2134×5238；物理码长5238，逻辑比特1552，码率8/27。
- 退极化 p_D=0.07125，f_m=0.0475；噪声种子20261122；两侧BP最多100轮。
- X侧保持原分阶段译码；没有完整HZ回退。Z侧保持同一个外层BP，无OSD、syndrome-flip等后处理。
- 两组各100帧，随机种子、噪声抽样流程完全相同。比较脚本核对参数、源文件哈希、矩阵清单及每帧X阶段标签。

## 方法

固定先验组使用每帧相同的256项概率表。条件先验组仅在完整X综合征通过时使用估计 xhat；否则沿用固定先验。通过综合征不保证估计正确，故此门控并非正确性证明。

对估计位 xhat_i，令

\[
p_i=\begin{cases}p_D/(3-2p_D),&\hat x_i=0,\\1/2,&\hat x_i=1.\end{cases}
\]

设 u 为8位输出模式、z9=b，则

\[
P_U(u\mid\hat x)=\sum_{b=0}^1 p_9^b(1-p_9)^{1-b}
\prod_{i=1}^8 p_i^{u_i\oplus b}(1-p_i)^{1-(u_i\oplus b)}.
\]

每个9位块生成自己的256项符号先验，共582张表。仅首次构造信道先验，BP每轮按标准变量更新使用该先验，不重复累计先验权重。

这是把X估计视为给定真值的硬条件模型，不是给定全部测量综合征的精确贝叶斯后验。不能由一次成功实验断言硬条件化始终优于软信息方法。

在线输入只有测量综合征、矩阵、信道参数、X译码估计。真实X/Z错误用于综合征生成、映射核对和离线评估，不用于先验生成或分支选择。

## 实现检查

`check_companion_conditioned_prior.py` 枚举512种X位型，并针对每种枚举512种Z错误，核对G2输出分布；在本信道点最大绝对差为0。额外核对共享先验与逐变量复制先验的BP输出、停止轮数和综合征状态一致。

每帧验证二进制/域综合征映射；输出综合征状态与完整物理矩阵核对。完整逻辑失败由两侧残差的二进制稳定子行空间成员关系判断，不以逐位精确恢复代替逻辑成功。

## 复现

在仓库根目录运行：

```bash
python3 experiments/gf256_cascade/check_companion_conditioned_prior.py
python3 experiments/gf256_cascade/run_companion_staged_xz.py \
  --frames 100 --p-depolar 0.07125 --seed 20261122 --max-iters 100 \
  --z-prior marginal \
  --output experiments/gf256_cascade/results/companion_conditioned_repeat/marginal.json
python3 experiments/gf256_cascade/run_companion_staged_xz.py \
  --frames 100 --p-depolar 0.07125 --seed 20261122 --max-iters 100 \
  --z-prior conditioned \
  --output experiments/gf256_cascade/results/companion_conditioned_repeat/conditioned.json
python3 experiments/gf256_cascade/summarize_companion_conditioned_pair.py \
  experiments/gf256_cascade/results/companion_conditioned_repeat
```

本次实际结果目录：`results/companion_conditioned_pair100_20261118/`。两组并发运行，保存的耗时是墙钟时间，不能视为干净的CPU开销对照。

## 完成结果

| 指标 | 固定边缘先验 | 实际X估计条件先验 |
|---|---:|---:|
| 完整逻辑失败 |70/100（70%）|4/100（4%）|
| X综合征失败 |4|4|
| Z综合征失败 |70|4|
| 漏检逻辑失败 |0|0|
| X精确恢复 |96/100|96/100|
| 使用条件先验 |0|96|

配对救回66帧，退步0帧，共同失败4帧，精确McNemar双侧p=2.710505431213761e-20。两组X阶段标签一致。条件组的96个X综合征通过帧在离线评估中均精确恢复X，且最终逻辑成功；剩余4帧X综合征失败，按预定规则沿用固定先验，Z也失败。

这说明在本固定矩阵、信道点及样本上，X估计提供的相关信息有明显帮助。尚不证明更低FER区间、多矩阵或X综合征通过但估计错误时同样有效。下一步优先冻结此方案做独立扩大样本；另可单独加入X阶段失败后的完整HZ回退，配对评估其是否救回剩余失败，不能直接假设能救回。

本次固定组70/100不能与此前边缘伯努利审计的80/100逐帧相减：随机种子和噪声生成不同，且此前主要是外层块错误指标。本次表格两列才是同帧完整逻辑FER对照。
