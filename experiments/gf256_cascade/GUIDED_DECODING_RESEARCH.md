# GF(256) 级联码：文献驱动的 Z 侧译码改进

日期：2026-11-02。状态：研究与先导实验，尚未证明达到可发表水平。

### 2026-11-03 文献补充：relay 调度与级联层级消息传递

1. **Relay-BP**：Müller et al., *Improved belief propagation is sufficient for real-time decoding of quantum memory*, arXiv:2506.01779 (2025)。作者先定义带每变量记忆强度 \(\gamma_j\) 的 DMem-BP，再把多个 DMem-BP leg 串接；下一 leg 从上一 leg 的软边缘继续，并在候选满足 syndrome 时按信道权重比较。其经验结果来自二进制 min-sum、CSS 的 XZ/XYZ 解码和 circuit-level noise；论文明确指出问题相关的 \(\gamma_j\) 分布（包括负值）很关键。它不是本项目 q-ary sum-product/Turbo 解码的直接性能基线。
2. **级联量子码层级消息传递**：Wu et al., *Learning to Decode Concatenated Quantum Codes with Hierarchical Message Passing*, arXiv:2608.28571 (2026-08-28 预印本)。它在量子 Hamming、many-hypercube 和非 CSS \([[8,3,3]]\) 递归级联码上双向传递向量软信念，并学习消息聚合；声称在其模型和码族上优于相应 hard-decision 基线。它确实支持“把高层约束反馈到底层、保留软信息”这一一般方向，但其递归编码图、量子码族、训练式组合器与本项目的 classical-GF(256) 外码 × path 内码矩阵/Turbo 译码并不相同，不能把其阈值或性能迁移为本项目预期。
3. **与当前 fixed-tuned 结果的区别**：代码中 `conditioned_robust_tuned` 的 `0.2` 是内外 Turbo 概率向量交换的混合系数（`decodeHxCascadeTurbo`），不是 Relay-BP 的 per-variable memory strength \(\gamma_j\)，也没有 Relay-BP 的逐 leg 软边缘传递。当前 tuned 优于 raw robust 只能先视为实现/调度参数改善，不能单独表述成 Relay-BP 或新的译码器贡献。
4. **Rigby 等人的相关性恢复与增强 BP**：Rigby、Olivier、Jarvis, *Modified belief propagation decoders for quantum low-density parity-check codes*, Phys. Rev. A 100, 012330 (2019), [arXiv:1903.07404](https://arxiv.org/abs/1903.07404)。该文在 CSS 码上用 GF(2) 分开译 X/Z 后，通过重新调整错误概率来恢复退极化信道的 X-Z 相关性；另一种 augmented decoder 在 BP 重启时重复一部分校验行，以提高这些校验消息的影响。作者指出，尝试次数带来的 FER 降幅依赖码族，而且仅降低检测到的 BP 失败不等于降低逻辑 FER。它没有为本项目 GF(256) 级联码提供现成参数。
5. **CAMEL 码设计与 ensemble BP**：Miao et al., *A Joint Code and Belief Propagation Decoder Design for Quantum LDPC Codes*, [arXiv:2401.06874](https://arxiv.org/abs/2401.06874)。其 ensemble decoder 对四元变量作不同硬猜测，并行运行多条 BP，再从满足 syndrome 的候选中按权重选择；码构造还将短环集中到被固定节点。可借鉴有限预算的多起点候选机制，但不能直接把其性能迁移到本项目，也不能简单遍历 GF(256) 符号的 256 种取值。

### 2026-11-03 Rigby-inspired augmented-check pilot

实现了一个仅供实验的 `conditioned_robust_augmented`：原始 tuned Z-side Turbo 未满足完整 syndrome 时，随机选择约 10% 的外层校验行并重复，做 3 次额外 BP/Turbo 尝试；重复行仅改变 BP 消息的相对影响，最终候选仍须满足未修改的完整 `H_X^{cas}` syndrome，并按已知 X 条件信道先验排序。该实现是机制移植，不是 Rigby 2019 二进制 decoder 的复现。

GF(8) 精确后验单测、编译和正常帧 smoke 均通过。预选困难帧 `20272001/trial=2` 上三次重试都没有找到完整 syndrome-valid 候选。20 帧配对 pilot (`p_D=0.07125`，矩阵种子 `20261018`，噪声种子 `20298000–20298003`) 中 tuned 与 augmented 都为 `1/20`，同一帧 Z syndrome failure，配对 rescue/regression `0/0`、McNemar `p=1.0`；CPU `9.49` vs `10.77 s/frame`。因此这个简单的随机重复外层行策略没有显示收益，且增加约 13.5% CPU；不扩大该实现的样本量。结果位于 `results/augmented_check_paired_pilot20_matrix20261018_pD0p07125/`。

### 2026-11-03 固定低阻尼调度的跨矩阵配对复验

第一矩阵（构造种子 `20261018`）在 `p_D=0.07125` 的新配对 300 帧中，raw robust `41/300`、fixed tuned `16/300`、adaptive `15/300`；raw→tuned 配对救回 26、退步 1（精确 McNemar 双侧 `p=4.17×10^-7`），raw→adaptive 救回 26、退步 0（`p=2.98×10^-8`）。tuned/adaptive 仅有 1 个 discordant pair（`p=1.0`），adaptive CPU `10.04 s/frame`，tuned `8.38 s/frame`。此前第二矩阵（种子 `20261019`）为 raw `46/300`、tuned `15/300`、adaptive `17/300`，tuned 相对 raw 救回 32、退步 1（`p=7.92×10^-9`）；adaptive 也未显著优于 tuned，且更慢。两组独立矩阵结果一致支持 fixed tuned 在这个特定信道点优于 raw robust，但目前只验证一个 `p_D`，不能外推成阈值或完整曲线，也不能当作新的理论译码算法贡献。逐批日志、同帧转移和 hash 分别在 `results/paired_three_way_matrix20261018_300/`、`results/paired_adaptive_vs_tuned_matrix20261019_300/`。

### 2026-11-03 文献补充：用内层 syndrome 改写外层先验

Berthusen et al., *Adaptive Syndrome Extraction*, **PRX Quantum 6, 030307 (2025)**, [DOI: 10.1103/ps3r-wf84](https://doi.org/10.1103/ps3r-wf84)。该文研究的是内层 `[[4,2,2]]` Iceberg error-detecting code 与外层 HGP qLDPC 的连接，不是本项目的 GF(256) 级联矩阵。与本项目最相关的是其软信息处理：内层非零 syndrome 表示相应小块有错误；应用固定内层修正后，作者据此上调受影响外层逻辑变量的错误先验（文中将两个可能受影响逻辑位的概率设为 `0.5`，其余位取 `p_2`），再交给外层解码器。核心启发是把内层观测保留为外层概率信息，而不是只把内层硬判决当成最终答案。该数值规则依赖 Iceberg 的特定错误映射，不能照搬到本项目。

**对本项目的可检验含义**：当前 `decodeHxCascadeTurbo` 在 Z 侧内部交换 GF(256) 概率向量，但旧的 X→Z 跨分量接口仅使用硬判决 `corrected_x` 的基坐标，经 `P(Z_b=1|X_b=\hat x_b)` 构造条件先验。已实现的 `conditioned_robust_softx` 不再依赖这一步硬位作为主要先验：它在完整 `H_Z` GF(256) 图上另跑 sum-product，若该软译码的硬候选通过完整 X 综合征检查，就将 APP 符号边缘通过 Pauli 转移映射为 Z 先验；否则回退到原硬 X 条件先验。它检验“保留 X 不确定性是否帮助 Z”，但 full-H_Z BP 可能受短环影响，也没有向外层投影迭代传播 H2 软消息。更贴合级联结构的后续候选仍是构造结构化 X-side Turbo posterior：外层各投影和 H2 因子双向交换 extrinsic GF(256) 消息，随后将物理 X 符号边缘映射成 Z 先验；需先做小规模精确后验核对，并与现有 softx/tuned 同预算比较。

本项目当前的基准 X 译码器对 `H_Z` 上半部逐坐标译外层，再用每个内块下半部 `H_2` 校验恢复 selector 未覆盖的符号；该分阶段步骤产出硬 `corrected_x`。新 softx 实验另对完整 `H_Z` 做 BP APP，但不是分阶段结构的软化版本，也没有把 `H_2` 消息反馈给前五个外层投影。更完整的候选是做结构化 X-side Turbo：外层每个投影输出 GF(256) extrinsic/posterior，内层按实际 `H_2` 行系数和 syndrome 用 GF 卷积更新所有连接符号，再将 cavity 消息送回外层迭代。收敛后得到每个物理 X 符号的边缘 `P(X_v=x|s_X)`；逐符号通过实际 `FieldBasis` 得到

\[
P(Z_v=z\mid s_X)=\sum_{x\in\mathrm{GF}(256)}P(X_v=x\mid s_X)
\prod_{b=0}^{7}P(Z_{v,b}=z_b\mid X_{v,b}=x_b),
\]

其中每个 `P(Z_b|X_b)` 由当前 `PauliChannel` 的 `p_i,p_x,p_y,p_z` 精确给出。该映射在符号内保留 X 位相关性；独立地把各物理位边缘相乘会额外丢失相关性，不应作为首选。实现时可用 8 层二元转移矩阵做 tensor transform，避免对每个物理符号显式做 `256×256` 求和。然后将得到的每符号 Z 先验送进现有 Z-side path Turbo，而不是仅传硬 X 向量。

这仍是近似解码：外层 BP 可能有 loopy bias，按每个 GF 符号边缘独立送往 Z 侧会丢掉不同物理符号间的联合相关；先在 GF(8) 小码上逐帧与包含上、下层 syndrome 的精确枚举后验核对，再做 GF(256) 配对 FER。必须保证 extrinsic cavity 语义和总计算预算公平。这样 `H_2` syndrome 经 Pauli X-Z 通道相关性影响 Z 译码，但不把 X/Z 的 CSS 综合征误当成同一约束。

**实现状态（2026-11-03）**：新增实验模式 `conditioned_robust_softx`，在完整 `H_Z` GF(256) Tanner 图上运行 sum-product，只有其硬候选通过 X 综合征时才将 APP 边缘转为符号级 Z 先验；若未通过，则回退到现有硬 X 条件先验。Z 侧保持 tuned Turbo schedule。边缘的 Pauli 通道变换用 8 层二元转移实现，避免 `q^2` 暴力枚举；新增 `test_soft_x_prior.cpp`，用 GF(8) 树因子图穷举校验 X 后验和 X→Z 符号先验。隔离 build 编译通过；精确枚举最大 X 后验差为 `3.33e-16`，X→Z 变换最大差为 `1.67e-16`。10 帧配对 pilot 中 tuned 和 soft-X 均为 `0/10`，配对救回/退步 `0/0`、McNemar `p=1`；软后验有效 `10/10`，CPU `7.955` vs `8.612 s/frame`（约慢 `8.3%`）。样本没有失败帧，不能判断 FER 增益；结果仅验证接线与先导耗时。后续需在有足够失败事件的新配对数据上检查逻辑 FER，额外全图 X BP 成本计入总耗时。

以上均是待验证假设。要确保 cavity/extrinsic 消息没有重复计入输入先验；并与固定 `conditioned_robust_tuned` 同帧、同计算预算比较量子逻辑 FER。不能仅凭软后验或 syndrome 通过宣称收益。

**100 帧配对开发结果（2026-11-03）**：矩阵种子 `20261018`、`p_D=0.07125`、噪声种子 `20296000–20296019`，tuned 与 softx 各 100 帧。逻辑失败同为 `4/100`（各 Wilson 95% CI `[0.0157,0.0984]`），配对救回/退步 `0/0`，McNemar `p=1.0`；两者具体是同一批 3 次 Z syndrome failure 和 1 次 syndrome-pass logical failure。softx 的软 X 后验有效 `100/100`，但平均 CPU `9.43` 对 tuned `8.77 s/frame`，约慢 `7.6%`。所以当前完整 `H_Z` BP-APP 版本在这组配对帧上没有改善逻辑 FER，只增加耗时；不继续做该版本的独立 FER 扩大验证。结果 `results/softx_paired_dev100_matrix20261018_pD0p07125/`。

**结构化 X 软译码原型（2026-11-03）**：完整 `H_Z` BP-APP 负结果后，新增 `runCascadeXSoftBP`：每个 selector 坐标使用外层 GF(256) Tanner 图，每个物理块使用实际重复的 `H_2` 子矩阵；两层仅交换各自校验消息的 extrinsic product，信道先验只在变量节点计入一次。每轮生成每个 GF(256) X 符号的 APP，之后通过精确 Pauli X→Z 符号转移提供给 tuned Z-side Turbo。代码会验证 `H_Z` 上块确实为 `H_outer⊗selector`、下块为逐块重复的 block-diagonal H2，否则拒绝运行。GF(8) 的无环小级联例与全枚举 APP 最大误差 `1.11e-16`；本地 GF(256)、内长 12 矩阵单帧 smoke 的 X syndrome 通过。配对 100 帧开发实验正在运行：矩阵种子 `20261018`、`p_D=0.07125`、噪声种子 `20297000–20297019`，比较 fixed tuned。尚无 FER 结果；即使 syndrome 通过，也必须看稳定子等价逻辑 FER。实现不声称与完整 loopy Tanner graph BP 完全等价，其消息交换近似按投影图与内块图分组的 Turbo-BP。



**本项目可检验的移植顺序**：fixed-tuned 已在两个矩阵种子的中心点配对实验中重复优于 raw robust；softx APP 的 100 帧开发比较为负。下一步若实现新候选，应优先评估 Rigby-inspired 的有限 augmented-check 重启，或与 cascade 结构一致的外层软后验反馈，并且只在固定 tuned 失败时增量尝试；设计一个不读取真值的 relay-inspired 原型：多个有预算上限的 Turbo leg，明确记录各 leg 的初始化信息、交换/记忆更新和候选选择规则；候选必须满足完整 syndrome，选择只能用译码器可得的先验/软信息。第一阶段只做定向帧和小配对开发集，固定总 BP/Turbo 调用预算并与 `conditioned_robust_tuned` 同预算比较。若无法在 GF(256) 概率向量层写出一致的 extrinsic 更新，或仅得到 syndrome 通过而逻辑 FER 不降，就停止该移植并记录负结果。

## 1. 问题与评价目标

固定 12 符号 path 内码、GF(256) 外码及矩阵种子 20261018。当前主要瓶颈是 `H_X^cas` 对 Z 分量的综合征不收敛。内码换成两个短校验组或对 path 系数加权的 20 帧实验均未显示净改善；保留这些负结果，不将它们当作已验证的优化。

研究假设：部分失败来自耦合 BP 的不良固定点。根据综合征和软信息约束少量物理比特，重新运行内层精确链消息及外层 GF-BP，可能恢复这些帧。

主要指标是稳定子等价意义下的量子逻辑 FER，同时报告 X/Z 综合征失败、漏检逻辑失败、CPU 秒/帧及配对转变。综合征通过不能替代逻辑成功。

## 2. 检索到的相关方法

| 方法 | 论文与来源 | 可借鉴部分 | 本项目的移植限制 |
|---|---|---|---|
| BPGD | Yao 等，2024，[Belief Propagation Decoding of Quantum LDPC Codes with Guided Decimation](https://arxiv.org/abs/2312.10950) | 按 BP 可靠性依次固定变量以改善收敛 | 原文二进制/四元 Pauli 译码，不是本项目的 GF(256) 级联译码 |
| GDG | Gong 等，2024，[Toward Low-latency Iterative Decoding of QLDPC Codes Under Circuit-Level Noise](https://arxiv.org/abs/2403.18901) | 保留若干变量的两种猜测，组成候选集合 | 原文电路噪声及滑动窗口；本项目先做无测量噪声实验 |
| LSD | Hillmann 等，[Localized statistics decoding](https://arxiv.org/abs/2406.18655)，2025 年期刊版 | 用可靠性引导局部区域的线性求解，减少全局 OSD 成本 | 需检查二进制展开与 GF 符号映射，并控制局部区域扩张 |
| MBP | Kuo、Lai，[Exploiting degeneracy in belief propagation decoding of quantum codes](https://doi.org/10.1038/s41534-022-00623-2)，2022 | 记忆项及抑制消息，缓解短环和简并带来的问题 | GF 外层迁移需单独推导，不能直接把二进制公式替换系数 |
| Relay-BP | Müller 等，2025，[Improved belief propagation is sufficient for real-time decoding of quantum memory](https://arxiv.org/abs/2506.01779) | 用变化的记忆强度接续迭代，打破振荡和对称性 | 原论文性能不等于本级联码性能，需独立基线实现 |
| Stabilizer Inactivation (SI) | Du Crest, Mhalla, Savin，[Stabilizer Inactivation for Message-Passing Decoding of Quantum LDPC Codes](https://arxiv.org/abs/2205.06125) | 在 MP 失败后，按低可靠度排序，逐一失活某个对偶校验支持；对不接触失活变量的子图重跑 MP，再解一个规模受校验权重限制的局部线性系统。作者报告其多种 qLDPC 构造上优于 MP-OSD 且复杂度较低 | 原文为二进制 CSS 码；本项目以 GF(256) 符号后验构造可靠度、以对偶外层行作支持，并对局部多解按后验评分，属于未验证的 q-ary 适配。源文见 [arXiv HTML 第 3 节](https://arxiv.org/html/2205.06125v3) |
| 内层 syndrome 到外层先验 | Berthusen 等，2025，[Adaptive Syndrome Extraction](https://doi.org/10.1103/ps3r-wf84) | 内层检测结果转成受影响外层变量的软错误先验 | 内码为 `[[4,2,2]]` 检错码，先验规则需按本项目 GF(256) 结构重新推导 |
| 逻辑空间 OSD | Shang，[Logical-Space Theory of Ordered-Statistics Decoding for Quantum LDPC Codes](https://arxiv.org/abs/2610.24503)，2026-10-17 预印本 | 将有限列表失败拆成逻辑陪集未覆盖与候选代表选择错误；用逻辑签名指导列表覆盖 | 新预印本，当前结构是 GF(256) 级联码；逻辑映射需从本项目 CSS/stabilizer 关系重新构造，不能直接套二进制公式 |

检索不构成穷尽性综述。上述论文的阈值、码长和噪声模型不可与本项目直接比较。

**后续实验优先级（2026-11-03 更新）**：fixed-tuned 是目前主要译码基线；adaptive 在两个矩阵种子的中心点没有显示优于 tuned。完整图和结构化 `softx` 已有 100 帧配对开发结果，均未改善逻辑 FER且更慢；随机增广校验在 20 帧 pilot 无增益，overcomplete 方案在 5 帧试验中约慢 12 倍。SI 的 GF(8) 局部精确测试通过，但在 `p_D=0.072` 的 100 帧 tuned 配对开发比较中逻辑 FER 同为 `13/100`、失败集合完全相同，且 SI 稍慢；故停止直接逐投影 SI 的预算调优和该版本 holdout。若继续 SI，应先设计共享物理符号下的多投影协调并验证；否则转向其他 Z 侧机制。任何开发增益都必须用冻结实现、新种子、独立矩阵/信道点复验；若只修 syndrome 而不改善逻辑 FER，应记为负结果。CAMEL 式多硬猜测仍暂缓，因为 q=256 候选空间大，当前结构不符合其短环集中前提。

### 2026-11-03：Stabilizer Inactivation 的本项目原型

SI 原文讨论二进制 CSS Tanner 图。译码某一误差分量时，选择对应对偶稳定子矩阵中的低可靠度校验行，暂时移除其支持变量参与 MP；若余下子图能满足保留 syndrome，再解回涉及失活变量的小线性系统。文中 SI 是 MP 的后处理，逐个增加失活校验，局部系统大小由被选 check 的支撑权重限制。其成功标准仍是 syndrome-valid 候选，量子逻辑成功还需再检查稳定子等价。

本项目外层投影为 GF(256) 校验。对 Z 侧外层矩阵 `HZ_TIT_Outer_For_HX`，候选支持来自其对偶 `HX_TIT_Outer_For_HZ`；当前程序在 GF(256) 上运行时验证该支持向量确实满足 `H r^T=0`。每个外层符号的可靠度暂用 `log(P_1/P_2)`（最大和次大符号概率之比）求和近似二进制 SI 的 `sum |LLR|`；局部 affine 解用 BP 符号 APP 评分。原文没有给出这套 q-ary 排序或评分，因此需单独精确核对。

代码模式 `conditioned_robust_si` 已接入并编译。20 帧配对 pilot（每个失败投影至多尝试 4 个对偶支持）在 `p_D=0.07125`、内长 12 path 矩阵上 tuned/SI 均 `0/20`，配对 rescue/regression `0/0`，CPU `9.3145/9.3115 s/frame`；样本无失败，不能估计收益。困难帧 `20272001/trial=2` 的预算 4 版本首轮尝试 9 次、修复 1 个投影；整帧仍有 40 个未满足 Z checks；重试尝试 8 次且无局部修复。之后将预算增至 12，首轮/重试分别尝试 25/24 次，仍只修复首轮同一个投影，整帧结果不变。增加尝试数没有解决该帧。

GF(8) 受控例中，以 `H=[[1,1,0],[0,0,1]]`、对偶支持 `[1,1,0]` 构造局部失活问题；SI 恢复的候选满足完整 syndrome，且在 8 个可能局部解中与精确后验 MAP 解一致。该核对验证小型局部消元/评分实现，不证明 GF(256) loopy cascade 中的可靠度排序或整体 FER。

随后在独立噪声种子 `20312000–20312099`、`p_D=0.072` 上完成 100 帧配对开发比较（当前矩阵、内长 12、每个失败投影最多失活 12 个 check）。tuned 与 SI 均 `13/100`，95% Wilson 区间均 `[0.0776,0.2098]`；失败集合完全相同，含 12 次 Z syndrome failure 和 1 次 undetected logical failure，配对救回/退步 `0/0`、McNemar `p=1.0`。SI 平均 `12.7608`、tuned `12.4831 s/frame`（+2.2%）。SI 共尝试 747 个首轮局部子问题并修复 24 个投影；重试阶段尝试 432 次且没有修复，整帧逻辑结果未变。这个版本没有逻辑 FER 收益，不安排独立 holdout，也不再盲目增加失活预算。结果 `results/si_paired_dev100_matrix20261018_pD0p072/`。

主要疑点是各外层投影共用物理符号。单独修复一个投影并不保证其他投影的边缘一致，更不保证重组后的完整级联 syndrome 通过。GF(8) 局部求解已通过穷举，但 100 帧逻辑 FER 比较为负，因此若继续此方向应研究多投影协调，而非单投影预算。结果目录 `results/si_paired_pilot20_matrix20261018_pD0p07125/`、`results/si_paired_dev100_matrix20261018_pD0p072/`；诊断日志 `/tmp/si_known_failure_replay_20261103.log` 和 `/tmp/si_known_failure_replay_budget12_20261103.log`。

## 3. 本次实现：有界物理比特分支

新增模式 `conditioned_robust_guess`：

1. 完整运行既有 `conditioned_robust_adaptive`；若综合征已通过，直接保留其结果。
2. 失败时，通过 GF 外层与内层链消息计算物理符号边缘概率。链中使用外层 extrinsic 消息，避免把内层先验重复计入。
3. 找到未满足校验所邻接的物理符号，将其 GF 概率分布边缘化成各物理 Z 比特的概率。
4. 按接近 0.5 的程度排序，选前四个比特；每个分别尝试 0、1，先试较可能的值。
5. 与猜测冲突的物理符号先验乘 `exp(-25)`，再运行 10 轮内外交换，每轮外层至多 20 次，阻尼 0.2。
6. 第一个通过完整综合征的候选被接受；最多八次分支。全部失败则保留原失败结果。

这是受 BPGD/GDG 启发的局部分支实验，并非原论文算法的完整复现。选择最不确定的局部比特也不同于 BPGD 的最可靠变量固定规则。

真实错误不参与分支选择、先验构造或停止条件。真实错误仅在最后用于稳定子等价与 exact recovery 统计。X 条件先验来自译码估计，仍可能含错，因此需要同时报告 X 结果。

## 4. 正确性核对与实验设计

`test_guided_cavity.cpp` 使用 GF(8)、单外层校验、三个物理符号，枚举全部 512 个错误向量，对照链边缘概率。覆盖普通/加权 path、无固定/固定 0/固定 1，以及有/无 X 条件先验，共十二组；最大误差约 1.18e-14，重建综合征均通过。这验证局部计算，不证明整个 GF(256) loopy BP 的统计正确性。

### 开发集

- `p_D=0.07125`，`f_m=0.0475`。
- 种子 20272000–20272007，每种子 5 帧，共 40 帧。
- 相同矩阵、相同物理噪声，比较 adaptive 与 guess。
- 这批帧之前被用于诊断，必须标为开发集，不用于宣称泛化。
- 先重放的两帧：20272003/2 被新分支完全恢复；20272001/2 八次分支后仍失败。不能用这两帧估计 FER。

### 独立验证与消融

冻结代码、参数、矩阵哈希后，先用新种子 20281102 起的独立 100 帧筛查；若无净改善，记录负结果并修改方法，不反复在同一验证集挑参数。每次调参后使用新的保留种子。

若有改善，扩大到每点 1000–5000 帧，并使用至少三个矩阵种子，覆盖瀑布区多个噪声点。除 adaptive 外，应加入更长迭代的计算预算基线、移除分支/移除内层交换的消融，以及可行的 OSD/LSD 或 Relay 类强基线。报告 CPU 与尾延迟，不仅比较迭代次数。

配对统计采用 McNemar 检验及 FER 置信区间；零失败报告上限。固定样本量的验证不根据中途好坏提前结束。若后续要按错误数停止，需要明确采用相应统计方案。

## 5. 研究结论的门槛

没有统一的“可发论文 FER”。需要：可重复的独立测试增益、适当的强基线、计算代价说明、清晰的贡献定位，以及无真实错误泄漏的实现。只移植已有猜测方法并得到少量改善，还不足以证明创新。当前先验证此结构上的困难与有效机制，再判断是否形成可发表贡献。

## 6. 条件信道中的歧义来源

若退极化模型为 P(X)=P(Y)=P(Z)=p_D/3，写物理 Pauli 错误为 (x,z)，则

- P(z=1|x=1)=P(Y)/(P(X)+P(Y))=1/2；
- P(z=1|x=0)=(p_D/3)/(1-2p_D/3)。

在 p_D=0.07125 时，后一概率约为 0.0249344。path 投影是相邻物理符号的和：在二进制同一坐标上，若两个 x 都为 0，投影错误概率为 2r(1-r)；只要其中一个 x 为 1，该投影的局部先验就是 1/2。这只是局部先验，校验仍可提供信息。共享符号造成投影间的相关性，不能将这些边缘概率当作互相独立的联合信道。

这些公式以正确的 X 信息为条件；实际译码器使用估计的 X，因此应核对其失败传播。后续可考察：对条件概率为 1/2 的物理比特做局部可靠性引导线性求解，并与当前分支搜索比较。但不能用仿真中的真实 X 替代估计 X 来生成主要结果。

## 7. 已完成的 40 帧开发集结果

| 模式 | 逻辑失败 | Z 综合征失败 | 漏检逻辑失败 | CPU 秒/帧 |
|---|---:|---:|---:|---:|
| adaptive | 4/40（10%） | 3 | 1 | 14.40 |
| adaptive + guided guess | 3/40（7.5%） | 2 | 1 | 18.43 |

配对救回 1 帧、退步 0 帧，CPU 增加约 28%。McNemar 双侧精确 p=1.0，尚无统计显著性。新旧 FER 的 Wilson 95% 区间分别为 [2.58%,19.86%] 与 [3.96%,23.05%]。不可宣称稳定提升 25% 或已经具备发表水平。

原始结果：`results/guided_guess_paired_development40/summary.csv`、`pairwise.csv`、`failure_transitions.csv`。两个残留 Z 综合征失败帧为 20272001/2、20272006/2；综合征通过但逻辑失败的 20272000/3 被有意保留，未用真值选择候选。

## 8. 独立 100 帧验证

种子 20281102–20281121，矩阵/程序哈希与配置见 `results/guided_guess_paired_holdout100/config.json`。adaptive 与 guess 均有 6/100 逻辑失败（FER 6%），其中 5 次 Z 综合征失败、1 次漏检逻辑失败。两者失败集合完全相同：配对救回 0、退步 0，McNemar p=1.0。Wilson 95% 区间均为 [2.78%,12.48%]。guess 耗时 15.65 秒/帧，adaptive 为 12.31 秒/帧，慢约 27%。因此本轮有界分支猜测无统计或实用收益，暂不建议将其当作译码器改进。

下一步已实现顺序可靠度固定的 BPGD 式试验模式 `conditioned_robust_bpgd`。先重放开发集里的两帧残余 Z 综合征失败（20272001/2、20272006/2）；若没有恢复，再检查条件局部线性修复/LSD 思路。任何方法若要证明 FER 改进，须在新的冻结 holdout 上复验，并同时比较耗时、漏检逻辑失败及等预算基线。

## 9. BPGD 式顺序固定的反例重放

新模式 `conditioned_robust_bpgd` 对两帧残余 Z 综合征失败样本各最多固定 8 个比特，再次运行耦合 BP。种子 `20272001/2` 的未满足校验数保持 135；`20272006/2` 从 344 降至 101，仍不满足完整综合征。两帧真实 X 均被译码精确恢复，但真实错误不用于选择固定变量。8 步选择的比特置信度均恰为 0.5，揭示该失败邻域存在大量条件先验本征不确定的符号位。由于样本只含两个预先选定的困难帧，不能据此估 FER；当前只说明此固定策略未修复它们。日志位于 `results/guided_bpgd_dev_seed20272001_trial2.log` 和 `...seed20272006_trial2.log`。

## 10. 投影列表与 syndrome 修复的逻辑类诊断

### 代表帧 oracle 诊断

对开发帧 `p_D=0.07125, seed=20272001, trial=2` 运行完整诊断（此处 true-X 仅用于 oracle 对照，不进入候选选择）：

- 当前 Turbo 候选有 135 个未满足 `H_X` 校验。
- 使用真实 X 先验仍有 102 个未满足校验；因此这帧的失败不能只归因于 X 估计错误。
- 直接对完整 `H_X` 做 GF(256) BP，200 次迭代后有 1922 个未满足校验。
- syndrome beam 从 135 降至 130；localized OSD 后仍有 133 个。
- 完整矩阵 OSD-1（65,281 个候选）和局部 OSD-2 都找到有效 syndrome 候选，但两者仍属于错误逻辑类。
- 当前 Turbo Z 候选相对仿真真实 Z 的物理信道 log-likelihood 低约 2826.68。真实错误仅用于诊断，不能用来选择算法候选。

这说明问题同时包含级联投影搜索的收敛困难和逻辑类选择困难。单纯提高 syndrome 通过率会把许多失败变成漏检逻辑失败。

### 三种候选路径重放

同一矩阵、同一信道种子下，预选两帧 `20272001/2`、`20272006/2`：

| 方法 | Syndrome 结果 | 逻辑结果 | 观察 |
|---|---|---|---|
| `conditioned_robust_joint_list`，原 OSD-1 列表 | 两帧均修到有效 | 两帧均逻辑失败 | 旧小列表在两个 syndrome 失败帧均无联合更新；对原本 syndrome 通过的帧不触发 |
| `conditioned_robust_joint_list`，有界 OSD-2 列表 | 两帧均修到有效 | 两帧均逻辑失败 | `20272001/2` 发生 7 次联合更新、联合 score 增加 334.506、残余符号从 687 降至 642；`20272006/2` 无更新。只对两个预选失败帧测试 |
| `conditioned_robust_outer_osd` | 两帧均修到有效 | 两帧均逻辑失败 | BP 投影有效数分别 6/11、8/11；OSD 补齐其余投影，但选到错误逻辑类 |
| `conditioned_robust_joint_osd` | 最终 syndrome 有效 | 仍逻辑失败 | 联合全矩阵 OSD 候选 score 比当前 robust 候选低 16340.5，未被接受 |

这些定向重放只诊断已知困难帧，不是随机样本 FER 评估。完整日志：`results/joint_list_known_fails_20272001_trial2.*`、`joint_list_known_fails_20272006_trial2.*`、`joint_list_undetected_20272000_trial3.*`、`outer_osd_replay_seed20272001.log`、`outer_osd_replay_seed20272006.log`、`joint_osd_replay_seed20272001_trial2.log`、`joint_list_osd2_replay_seed20272001_trial2.log`、`joint_list_osd2_replay_seed20272006_trial2.log`。

## 11. 扩大 OSD 列表后的陪集覆盖审计

逻辑空间 OSD 预印本将列表失败分成“真实逻辑签名未覆盖”和“虽覆盖但选择了错误代表”，并通过自由变量到逻辑商空间的线性映射来设计候选列表。该论文推导和实验针对二进制 HGP/BB 码；本项目只借用覆盖诊断思路，GF(256) 的陪集签名由本项目自己的 CSS 稳定子商映射计算，不能把该论文的二进制结论直接当作本级联码定理。参见 [Shang, 2026, arXiv:2610.24503](https://arxiv.org/abs/2610.24503)。

对 `p_D=0.07125, seed=20272001, trial=2`，完整 `H_X^{cas}` 的 OSD 消元秩为 2,134，符号变量数为 6,984，即有 4,850 个自由符号。离线真值审计显示真实错误在这组自由坐标中有 566 个符号偏离 BP 硬判决；低可靠度前 256/512 个自由符号中分别已有 120/213 个偏离。该统计只解释有限阶列表为什么可能漏掉目标类，真实错误没有进入在线候选生成或评分。

| 候选列表 | syndrome-valid 词数 | 不同逻辑陪集数 | 覆盖真实陪集 |
|---|---:|---:|---:|
| 弱枢轴 OSD-1，256 个自由符号、每符号 8 个替代值 | 2,049 | 2,008 | 否 |
| 弱枢轴 OSD-1，512 个自由符号、每符号 8 个替代值 | 4,097 | 4,014 | 否 |
| 弱枢轴 OSD-2，32 个自由符号、每符号 2 个替代值 | 2,049 | 1,921 | 否 |
| 弱枢轴 OSD-2，64 个自由符号、每符号 2 个替代值 | 8,193 | 7,687 | 否 |
| 弱枢轴 OSD-2，128 个自由符号、每符号 1 个替代值 | 8,257 | 8,002 | 否 |

所以该帧的主要障碍不是当前列表中的代表排序，而是这些列表都没有产生真实逻辑陪集的代表。把 OSD-2 支撑增至 128 个自由符号、并将候选数保持在约 8.2k，仍未覆盖真实类。继续无差别扩大物理候选会增加时间和内存；后续优先研究逻辑签名引导的候选组合/流式列表，或改善 BP/Turbo 软信息，使候选中心接近正确逻辑类。以上仍是单个预选困难帧的机制诊断，不代表 FER 改善。

审计日志：`results/osd2_logical_class_audit/seed20272001_osd1_weak_free256_trial2.log`、`seed20272001_osd1_weak_free512_trial2.log`、`seed20272001_osd2_weak_free32_alt2_trial2.log`、`seed20272001_osd2_weak_free64_alt2_trial2.log`、`seed20272001_osd2_weak_free128_alt1_trial2.log`。

## 12. Syndrome 仿射空间 Gibbs 采样原型

### 12.1 级联码的逻辑商空间规模

重新按当前 `construct.log` 与块矩阵结构核对了维数。外层 TIT 两个 GF(256) 矩阵尺寸均为 `194×582`、行秩均为 194；路径内码长 12，`H2` 为 `1×12` 全一校验，`G2` 为 `11×12` 满行秩生成矩阵，selector 为 `11×12` 满行秩。于是

\[
n=582\cdot12=6984,
\qquad \operatorname{rank}H_X=194\cdot11=2134,
\]
\[
\operatorname{rank}H_Z=194\cdot11+582\cdot1=2716.
\]

第二个秩式中，HZ 上块在每个内块的第 12 列为零，下块 `I_{582}\otimes H2` 在对应第 12 列有唯一单位支撑，因此下块的 582 行与上块行空间直和。由 CSS 正交关系，编码维数为

\[
k_{GF}=n-\operatorname{rank}H_X-\operatorname{rank}H_Z=2134
\]

个 GF(256) 符号，即物理二进制展开后 55,872 个物理量子比特编码 17,072 个逻辑量子比特。对固定 Z syndrome，`ker(H_X)` 有 `6984-2134=4850` 个 GF 自由度；模去 `row(H_Z)` 后仍有 2,134 维 GF(256) 逻辑商空间。这意味着可区分逻辑类数量为 \(256^{2134}\)，通用 MCMC 即便保持 syndrome，也必须在极高维逻辑商中找到正确类。

这是根据当前构造日志中外层秩、内码/selector 秩及块支撑分离得到的结构秩；不是对任意改过的矩阵都成立。之后若矩阵文件改变，应重新验证两外层秩和 CSS 正交性，再重算这些参数。

2026 年逻辑空间 OSD 预印本把二进制 OSD 自由坐标映到逻辑签名，并将候选失败拆成“正确类未覆盖”和“候选类选择错误”；其分析针对 HGP/BB 二进制码，不能直接套成本项目的 GF(256) 定理，但其**先构造商空间签名再设计列表**的方向与本项目当前诊断相符。[原文](https://arxiv.org/abs/2610.24503)

为避免继续构造大量低阶 OSD 物理候选，新增了定向研究原型：先用 GF(256) 消元参数化所有满足 `H_X^{cas} z=s_Z` 的 correction，再对一个自由符号的一个多项式基比特做热浴更新，并同步修改主元符号。每次更新都沿 `ker(H_X^{cas})` 的方向移动，因此保持完整 syndrome。更新概率按 `P(z_i|x_i)` 的物理信道先验计算；本目标帧的译码 X 恰好与仿真真值一致，但算法只接收译码器给出的 X 估计。

保留的样本按精确 `cascadeZCosetKey` 分组；采样频率近似逻辑类后验质量，输出最高频类中似然最高的已见代表。此方法的数学目标合理，但有限链的类频率只有在充分混合后才能解释为后验类质量；单链、短链或高自相关样本可能严重偏差。

`seed=20272001/trial=2` 的定向结果：

| sweep 数 | 保留样本 | 采到的逻辑类数 | 最高频类样本数 | 采到真实类 | 输出 syndrome / 逻辑 |
|---:|---:|---:|---:|---|---|
| 30 | 20 | 8 | 6 | 否 | 通过 / 失败 |
| 200 | 38 | 25 | 5 | 否 | 通过 / 失败 |
| 200（独立 seed 20261111） | 38 | 17 | 10 | 否 | 通过 / 失败 |
| 200，整符号热浴每 sweep 更新 32 个变量 | 38 | 33 | 2 | 否 | 通过 / 失败 |

这只能说明链可以产生多个 syndrome-valid 逻辑类，尚不能说明它能采到或选中正确类。两条 200-sweep 独立链的类别数为 25 和 17，最高频分别为 5/38 和 10/38，显示当前样本量下分布依赖链种子。未计算自相关/有效样本量、未证明多链收敛、未运行随机帧 FER。故目前不应接入默认译码或宣称提升。

增加整 GF 符号 heat-bath 后，链探索到 33 个不同类，但最高频仅 2/38；6,400 次符号更新只有 175 次改变。当前两种更新方式均未采到正确逻辑类，且没有稳定的类后验估计。下一步如继续 MCMC，应测试能跨越高能垒的块提议或温度交换，并加入链间分布、自相关和有效样本量诊断。

#### 并行温度和稳定子更新审计（2026-11-02）

已实现可选 replica-exchange 版本：每条链仍在精确 syndrome 仿射空间内按逆温度 `β` 采样，交替尝试相邻温度状态交换；只用 `β=1` 冷链的样本统计逻辑类。冷链记录类连续保持比例及后验分数 lag-1 相关，用来暴露明显的链黏滞，但这些单一诊断量不能证明收敛或给出可靠 ESS。

在同一困难帧 `seed=20272001/trial=2` 上：

| 逆温度 β | sweeps / 冷链保留样本 | 冷链不同类 / 众数 | 交换接受 | 类持续率 / score lag-1 | 真类采到 |
|---|---:|---:|---:|---:|---|
| `1,.5,.2` | 200 / 38 | 29 / 4 | 3/200 (1.5%) | 未记录 / 未记录 | 否 |
| `1,.8,.6,.4` | 200 / 38 | 31 / 3 | 23/300 (7.7%) | 未记录 / 未记录 | 否 |
| `1,.9,.8,.7,.6,.5` | 400 / 75 | 30 / 13 | 167/1000 (16.7%) | 0.608 / 0.946 | 否 |
| `1,.8,.6,.4`，再加每 sweep 32 次 HZ 稳定子更新 | 200 / 75 | 21 / 18 | 14/300 (4.7%) | 0.730 / 0.962 | 否 |

更密温度梯让交换更频繁，但 score lag-1 仍为 0.946；仅 75 个冷链样本不足以支撑类后验估计。新增稳定子更新通过了逐行验证：每个 HZ 行都保持完整 H_X syndrome，且保持当前级联逻辑类签名。它能沿稳定子方向改善同一类内的代表元，但这次实验未能改善不同逻辑类之间的混合。四组均输出 syndrome-valid 但逻辑错误的候选；这不是 FER 评估，也没有逻辑 FER 收益。

因此，当前不应把普通并行温度链扩大为 FER 译码器。若继续 MCMC，下一步需要设计能在逻辑商空间之间有效移动的提议，再做多链独立性、类频数稳定性和长链诊断；若没有结构化的低能逻辑移动，应转向其它有明确验证路径的译码器，而非继续微调温度间距。

#### 稀疏逻辑坐标分解与定向提议

利用 RREF 的 4,850 个自由坐标，将完整 `H_Z` 的每行限制到这些坐标上。由于每个 `H_Z` 行在 `ker(H_X)` 中，且自由坐标对该核空间是一一参数化，这个限制矩阵与 `row(H_Z)` 同秩。对它做稀疏 GF(256) 消元得到秩 2,716；其非主元自由坐标组成 2,134 个独立逻辑补空间坐标。于是 MCMC 可以分别做 HZ 行热浴（稳定子、类内移动）和这些逻辑坐标热浴（跨逻辑类移动）。模拟器输出并验证了 `kernel_dimension=4850, stabilizer_rank=2716, logical_dimension=2134`。

在困难帧 `seed=20272001/trial=2`、温度 `1,.8,.6,.4`、120 sweeps、每 sweep 32 个一般自由符号提议及 32 个稳定子更新下，增添 32 个定向逻辑坐标更新得到 90 个样本、26 个类、众数 9/90；关闭定向逻辑更新的同 seed 消融得到 90 个样本、24 个类、众数 11/90。两边都未采到真实类，score lag-1 分别为 0.963 和 0.968。这个差别不足以证明逻辑覆盖增加；加逻辑更新也增加了一半提议数，所以需要等预算对照和更多独立链。它仍然是一次单帧链混合审计，不是 FER 结果。

这次坐标分解提供了明确的逻辑维数和可用的逻辑方向，但还没有构造低代价/高后验的逻辑候选列表。下一步应利用方向的物理支撑和信道代价设计候选排序/组合，在定向失败帧上先核对真实类覆盖，再决定是否做随机帧 FER。

采样器位于 `simulate.cpp` 的 `sampleGfSyndromeCosets`。定向启用环境变量 `GF256_CASCADE_ENABLE_GIBBS_AUDIT=1`；可设置 `GF256_CASCADE_GIBBS_AUDIT_SWEEPS`、`GF256_CASCADE_GIBBS_AUDIT_SAMPLE_EVERY`、`GF256_CASCADE_GIBBS_AUDIT_SEED`、`GF256_CASCADE_GIBBS_AUDIT_BLOCK_UPDATES`（0 表示逐位更新，正数表示每 sweep 的整符号更新数）、`GF256_CASCADE_GIBBS_LOGICAL_UPDATES`、`GF256_CASCADE_GIBBS_BETAS`（逗号分隔的逆温度，首项须为 1）以及 `GF256_CASCADE_GIBBS_STABILIZER_UPDATES`。`GF256_CASCADE_GIBBS_AUDIT_ONLY=1` 在输出采样诊断后直接结束，跳过后续 OSD/调度扫描。对应日志位于 `results/osd2_logical_class_audit/`，文件名前缀 `seed20272001_gibbs_`。

### OSD-2 联合列表配对比较

OSD-2 列表在固定的 20 帧开发集（`p_D=0.07125`, `f_m=0.0475`, seeds 20272000–20272003）与相同基础模式 `conditioned_robust` 配对：

| 模式 | 逻辑 FER | Z syndrome failures | CPU 秒/帧 |
|---|---:|---:|---:|
| `conditioned_robust` | 8/20 (0.40) | 7 | 11.66 |
| `conditioned_robust_joint_list` | 6/20 (0.30) | 0 | 13.32 |

配对救回 2、退步 0，McNemar 精确双侧 `p=0.5`。两次救回都来自 syndrome 失败转成逻辑成功；另 5 次 syndrome 失败变成漏检逻辑失败。该开发结果有正向信号，但样本过小。

在冻结源代码和矩阵哈希后，用全新种子 20282000–20282003 做 20 帧独立配对：两种模式均 4/20 逻辑失败；OSD-2 把 4 次 Z syndrome failures 全部变成漏检逻辑失败，救回 0、退步 0，McNemar `p=1.0`。CPU 从 8.94 增至 9.87 秒/帧（约 10.3%）。因此开发集提升没有复现，OSD-2 目前只改善 syndrome 成功率，并未改善量子逻辑 FER。配置与逐批日志在 `results/joint_list_osd2_holdout20/`，开发集在 `results/joint_list_osd2_vs_base20/`。不能将该方法报告为有效译码改进，也不能将其提升到默认模式。

### 双投影搜索负结果

为越过单投影贪心局部最优，曾短暂加入每次同时替换两个外层投影的有界搜索。10 帧开发配对集（seeds 20272000–20272001）基础模式与双投影模式均为 3/10 逻辑失败，救回 0、退步 0，McNemar `p=1.0`；CPU 从 11.07 增至 12.68 秒/帧。独立样本 `20282000/2` 未作任何双投影更新；`20282000/4` 做 6 次更新，分数仅比单投影 OSD-2 多约 0.7，仍是逻辑失败。未发现 FER 收益，因此已撤回该额外搜索代码。结果位于 `results/joint_list_pairwise_dev10/` 与 `results/joint_list_pairwise_replay_*`。

### 当前解释与下一步

当前证据支持 Z 侧级联译码器存在失败；还不足以断定矩阵构造本身错误。特别是 outer OSD 与全矩阵 OSD 都能恢复 syndrome、却仍选错逻辑类，说明必须对稳定子陪集进行候选聚合或可靠性引导局部解，而不能仅选择最可能的单一物理错误向量。`joint_list` 已加一个有界 OSD-2 列表：20 帧开发集的 FER 有小幅改善，但独立 20 帧没有任何配对救回，且 syndrome failures 转成漏检逻辑失败。双投影搜索也未改善 10 帧开发 FER，已撤回。下一步应先解决跨投影逻辑类的候选聚合/可靠性引导局部解，再开展更大规模统计验证。现阶段没有达到可发表改进。

## 11. 物理二进制 BP+LSD 定向先导

### 映射检查

级联矩阵使用 `HX_Cascade_Physical_Binary.txt`，形状为 `17072 × 55872`；对应 GF(256) 矩阵为 `2134 × 6984`，每个 GF 符号展开成 8 个物理比特。模拟器的 `validateBinaryExpansion` 在启动时逐行核对 GF 系数乘法块；Z syndrome 由该二进制矩阵直接计算。导出文件将 syndrome、译码得到的 X 估计与仿真真值分开，BP+LSD 只接收 syndrome 和 X 估计导出的条件 Z 信道；真值只在译码后用于逻辑等价检查。

### 单帧结果

使用 `ldpc==2.1.0` 官方实现，在开发失败帧 `p_D=0.07125, seed=20272001, trial=2` 上对完整物理二进制 (H_X) 做 Z 侧 BP+LSD：30 次 min-sum BP 迭代，`LSD_0`，并行更新，缩放系数 0.625。该帧本地译码得到的 X 估计恰好等于仿真 X 真值（真值未提供给译码器），所以这里 Z 条件信道的 X 侧信息准确。

| 译码器 | Z syndrome | Z 稳定子等价 | 帧逻辑结果 | 单帧耗时 |
|---|---|---|---|---:|
| 本地 `conditioned_robust_adaptive` | 失败，135 条未满足校验 | 不适用 | 失败 | 约 12 秒量级（既有基线） |
| 官方 BP+LSD `LSD_0` | 通过 | 失败，落在错误逻辑陪集 | 失败 | 约 3.5 分钟 |

此为预选困难帧，不是 FER 样本。BP+LSD 证明完整二进制矩阵可以找到 syndrome 解，但未解决逻辑类选择，也没有实用的当前点耗时优势。由于该参数点的条件 Z 比特先验在 `0.0249344` 至 `0.5` 之间，且出现大量高不确定比特，LSD 所利用的局部小簇假设可能不适配此操作区间；需要以簇大小统计验证，不能仅据单帧推断。

同帧尝试 `LSD_CS` 一阶局部重处理超过 5 分钟，内存约 1.5 GB 后人工终止；无解码结果可报告。这个测量只说明该设置/帧计算代价过高，不等价于算法理论上无法完成。

### 新的算法设计指向

2026 年 10 月的逻辑空间 OSD 预印本将有限列表的逻辑失败分解为 **logical coverage miss** 与 **wrong representative selection**。本项目已在 `simulate.cpp` 中实现 GF(256) 级联稳定子商空间的规范余数签名，并对完整矩阵 OSD-2 列表做了事后覆盖审计。两个预选失败帧（`seed=20272001` 和 `20272006`，均 `trial=2`）各生成 7,285 个 syndrome-valid 候选、7,285 个不同逻辑陪集；真值所在陪集在两帧都未覆盖。第一帧另测并列项随机排序，以及 8、32、64、256、1024 项可靠度窗口扰动。4 个窗口列表各有 7,285 个候选，合并后有 29,140 个不同陪集，仍未覆盖真实陪集；窗口 1024 的单候选后验分数更高，但逻辑类仍错误。前 64 个候选的规范标签与直接稳定子等价测试一致（0 个不匹配）。列表中每个陪集仅有一个候选，因此 class-mass 聚合无法修复这个覆盖缺口。审计只在候选生成后使用真值计算覆盖标签，不参与译码。

所以当前下一步应先扩大候选**覆盖**：对同一 syndrome 更换可靠度排序/消元基，合并多个 OSD 候选列表并统计逻辑陪集覆盖、运行时间和内存；只有某个列表覆盖了正确逻辑类，才值得研究类质量聚合或代表元选择。两帧结果是定向困难样本，不能代表总体 FER，也不能据此宣称新译码器性能改善。

先导数据与候选输出位于 `results/binary_lsd_pilot/`；运行脚本为 `run_binary_lsd_pilot.py`，其输入/真值文件分离。Python 依赖 `ldpc==2.1.0` 位于隔离环境 `/tmp/gf256-lsd-venv`，未增加项目运行依赖。

### Masked-Turbo beam 约束修复与叶节点修复先导（2026-11-02）

在困难帧 `p_D=0.07125, seed=20272001, trial=2` 上复核 masked-Turbo beam 时发现一个实现错误：内长 12 的两组 path 被 decoder 递归拆成两个内长 6 子问题，但递归调用没有传入各组的物理 bit constraints。因此第一轮宽度 4、深度 4/8 的 beam 日志不能用于评价正确实现，已在 `simulate.cpp` 中映射并传递子问题约束。

修复后，分支输出一致性按返回的物理 GF 符号逐位检查；不满足已分支条件的节点会被淘汰。约束惩罚仍通过信道符号先验中的 `exp(-λ)` 实现，默认 `λ=25`。另外增加 `GF256_CASCADE_CONSTRAINT_PENALTY_EXPONENT`、beam 内部 Turbo 预算环境参数，以及只允许修改未固定物理位的叶节点 syndrome beam 修复。

| beam 设置 | 调用数 / 满足分支约束 | 保留节点 | 叶修复成功 | syndrome-valid 候选 | 完整定向诊断耗时 / 峰值内存 |
|---|---:|---:|---:|---:|---:|
| width 4, depth 4, 3×10 Turbo, λ=25 | 21/20 | 13 | 未启用 | 0 | 98.66 s / 312 MB |
| width 4, depth 4, 10×20 Turbo, λ=25 | 21/16 | 13 | 未启用 | 0 | 186.18 s / 312 MB |
| width 4, depth 4, 3×10 Turbo, λ=80 | 23/23 | 14 | 未启用 | 0 | 101.28 s / 312 MB |
| width 4, depth 4, 3×10 Turbo, λ=80，受约束 leaf repair | 23/23 | 14 | 0/4 | 0 | 100.48 s / 312 MB |
| width 8, depth 6, 3×10 Turbo, λ=80，受约束 leaf repair | 63/63 | 38 | 0/8 | 0 | 148.76 s / 465 MB |

表中耗时包含目标困难帧原有的 BP、OSD 与 Gibbs 诊断，不等于生产译码器每帧时延；所有行都只是同一预选失败帧的机制诊断，不是 FER 比较。提高惩罚能让分支约束得到满足，但在当前 Turbo 硬判决周围增加受约束 syndrome 修复，仍没有产生可行词。高迭代预算甚至降低了约束一致性，说明“多跑 BP”本身不是解决方案。当前暂停继续扩大 beam；下一步应设计更贴合 path 级联因子分解的投影变量更新/精确 syndrome 求解，或另行构造 PEG 内码并同步开发对应译码器。真实 Z 仅用于 syndrome-valid 后的逻辑类覆盖审计，未参与在线候选选择。

实验日志：`results/logical_candidate_audit/seed20272001_masked_turbo_beam_constraint_fix_w4_d4_audit_only_trial2.log`、`seed20272001_masked_turbo_beam_constraint_fix_w4_d4_turbo10_outer20_trial2.log`、`seed20272001_masked_turbo_beam_penalty80_w4_d4_trial2.log`、`seed20272001_masked_turbo_leaf_repair_w4_d4_trial2.log`、`seed20272001_masked_turbo_leaf_repair_w8_d6_trial2.log`。

## 投影分歧风险指标的 100 帧审计

在固定矩阵、`p_D=0.07125` 下，用新种子 `20283000–20283003` 评估每帧初始外层投影硬判决与 BP 硬判决之间的分歧数。100 帧有 5 次逻辑失败，AUC 为 0.953，风险最高的 20 帧包含 5 次失败；Wilson 95% FER 区间约为 `[0.0215, 0.1118]`。失败帧的分数为 40、43、62、99、209，而分数为 0 的帧有 54 个。

该 AUC 不能解释为逻辑失败预警成功：5 次失败全部是已被 syndrome check 检出的 Z 综合征失败，没有 syndrome-valid 的逻辑失败。因此分歧数只显示与当前已知 syndrome 失败高度相关，不能证明它提供额外译码决策信息，也不能证明降低 FER。结果保存在 `results/projection_risk_audit/holdout100/summary.json`、`projection_risk_frames.csv` 和各批次日志。下一步先测试更长 Turbo 迭代方案在全新配对种子上的逻辑 FER 与 CPU 代价；风险指标只有在足量 syndrome-pass 逻辑错误上显示条件区分能力后，才考虑用于自适应触发。

## 新矩阵上的长迭代配对试验

在当前内长 12 path 矩阵（校验矩阵哈希见 `results/paired_adaptive_vs_long_fresh100/config.json`）上，使用 `p_D=0.07125`、新种子 `20284000–20284003`，对 `conditioned_robust_adaptive` 与 `conditioned_robust_long` 做相同噪声帧配对比较，各 100 帧。adaptive 为 `3/100`（Wilson 95% CI `[0.0103,0.0845]`），long 为 `5/100`（`[0.0215,0.1118]`）；相对 adaptive，long 救回 1 帧、退步 3 帧，精确 McNemar `p=0.625`。CPU 分别约 `10.23` 与 `8.63` 秒/帧。该样本不支持加长 Turbo 迭代改善 FER；long 作为负结果记录，不推荐替换 adaptive。完整汇总见 `results/paired_adaptive_vs_long_fresh100/summary.csv`、`pairwise.csv` 和 `failure_transitions.csv`。

## Restart Belief 的适配方向

[Restart Belief](https://arxiv.org/abs/2511.13281) 的核心不是单纯增加 BP 迭代：根 BP 失败后，按输出 LLR 排序，逐步固定候选错误位、更新 syndrome 并从头运行 BP；不同初始猜测形成独立分支，候选按合法 syndrome 与错误权重/早停条件选择。论文报告其在多种 binary/GF(4) QLDPC 码上的仿真表现，并提供 [QBelief 参考实现](https://github.com/DiegoForlivesi/QBelief)。

Restart Belief 原文描述为：根 BP 输出 LLR 全局升序排序；每个外层分支从对应根排序位置开始，固定该错误位为 1；分支内 BP 未收敛时继续固定当前输出 LLR 最低的**全体变量位**，更新 syndrome 并重启 BP。参考 [arXiv:2511.13281，算法 1 与 III-A](https://arxiv.org/html/2511.13281v1)。文中还依赖码距纠错半径做早停，本项目码距未知，不能照搬该保证。

旧实验模式 `conditioned_robust_restart_belief` 与原文不同：其根分支及分支内新位都限制为当前未满足校验所接触的位。新增 `conditioned_robust_restart_belief_global`，改为对 Turbo 产生的物理比特后验 LLR 全局排序。两者都在级联 Turbo 上使用有限惩罚近似固定约束、以 Z syndrome 为可行性判据，并按条件信道 log-prior 选候选；因 Turbo 软后验与原论文 binary BP 不同，仍是 Restart-Belief-inspired 适配，而非论文算法的逐项复现。

在固定 path 内长 12 矩阵、`p_D=0.07125` 上，对 5 个预选基线 Z syndrome 失败帧运行 4 分支、每分支最多 4 步的重放（惩罚指数 80）。2/5 帧找到并选中正确候选（exact recovery），另 3 帧没有产生 syndrome-valid 候选；所有帧的 X 分量均 exact。真值只在译码后用于 exact-recovery 审计，没有参与候选生成或评分。逐帧结果见 `results/restart_belief_known_failures/summary.json` 和对应日志。由于样本是按基线失败预选的，这只能证明该回退能救回部分已知困难帧，不能估计总体 FER。

在未参与开发的新种子 `20285000–20285009` 上完成 100 帧同噪声配对比较（`p_D=0.07125`、4 workers、局限于未满足校验的 4×4 版本）。adaptive 为 `5/100`（Wilson 95% `[0.0215,0.1118]`），旧 Restart-Belief 版本为 `4/100`（`[0.0157,0.0984]`）；配对救回 1 帧、退步 0 帧，精确 McNemar 双侧 `p=1.0`。Z syndrome failures 从 4 降至 3，undetected logical failure 两边均为 1。平均 CPU 从 `10.50` 增至 `13.93 s/frame`，慢约 33%。因此该样本没有证明 FER 增益；结果在 `results/paired_restart_belief_fresh100/`。

全局排序版在另一组新种子 `20286000–20286009` 完成同样 100 帧配对。adaptive 为 `4/100`（Wilson 95% `[0.0157,0.0984]`），全局排序版为 `2/100`（`[0.0055,0.0700]`）；救回 2、退步 0，精确 McNemar 双侧 `p=0.5`；Z syndrome failures `4→2`，漏检逻辑失败均为 0。CPU 从 `11.26` 增至 `14.53 s/frame`，慢约 29%。点估计有改善但只有两个配对差异，证据仍弱。结果 `results/paired_restart_belief_global_fresh100/`。

下一步将 4×4 全局排序版本冻结，在新的 300 帧同噪声配对样本上复验，记录 FER、配对转移、综合征失败、漏检逻辑失败、CPU 和回退触发率。未获得可重复的逻辑 FER 增益前不扩大常规仿真部署。

对 3 个在 4×4 预算下未救回的预选帧，将预算扩到 8 分支×8 步（各最多 64 次 BP）后，`20283001/trial=1` 与 `20283003/trial=24` 均仍为 0 个 syndrome-valid 候选；`20283002/trial=6` 找到 1 个候选并 exact recovery。`20283003/trial=24` 的目标序列共运行 25 帧，耗时约 14 分钟；这显示该搜索预算代价很高。结果在 `results/restart_belief_8x8_known_failures/`。3 帧定向开发样本只说明更大预算有时能覆盖正确候选，不能说明总体 FER 改善。下一步停止无差别加深搜索，优先研究对失败帧更有区分力的分支变量排序或利用 H2 投影的软信息更新；新方法必须用独立配对种子验证。

## 2026-11-03 Augmented-check 与 layered outer BP 进展

Rigby 等人 2019 年的 augmented decoder 通过随机复制检查节点形成增强因子图，按原综合征复制相应 syndrome 约束，并多次重启。当前代码移植了“重复外层检查节点”机制；均匀随机选行和按残差数加权选行两个版本在同一 20 帧配对 pilot 都未优于 tuned，各为 `1/20`，同一帧失败，且分别慢约 13.5% 和 14.1%。代表性困难帧的三次重试都没生成完整 syndrome-valid 候选。此结果说明当前简单行重复策略对该实例无效，不否定论文原算法在论文码族中的结果；本项目不再扩大这两个版本。

另外新增 `conditioned_robust_layered`：在外层 GF(256) BP 中，每处理一行校验便立即刷新相邻变量的外发消息，使后续校验能使用本轮最新消息；其余级联 Turbo 参数沿用 fixed-tuned (`10` rounds、`20` outer iterations、mixing `0.2`)。该算法已编译并通过普通 GF(256) 单帧 smoke，20 帧 pilot 与 flooding tuned FER 相同（均 `1/20`，同一帧 Z syndrome failure），CPU 为 `9.44` vs `10.06 s/frame`。现用 100 个不重叠新噪声帧验证 FER 非劣性及速度差异；不得将当前小样本速度下降表述为稳定加速或性能提升。

## 2026 年新文献：利用 BP 与 OSD-0 分歧自适应分配搜索

Wang 与 Coveney, *Adaptive decoding of quantum LDPC codes through decoder disagreement*, [arXiv:2609.37629](https://arxiv.org/abs/2609.37629)（2026-10-01 预印本）。文章用 BP 硬输出和相同 BP 后验可靠度构造的 OSD-0 syndrome-consistent 输出之间的 Hamming 距离作风险信号，只把高风险样本升级到截断的一阶自由变量搜索；并与 syndrome weight、BP residual 和随机路由作同预算比较。其结果来自二进制检测模型、特定 BB/radial 码与 circuit-level noise，不能直接转移到 GF(256) 级联码。

对本项目的启发是预测“额外搜索是否值得”，而非仅预测 syndrome failure。已有 projection-disagreement 风险审计的 AUC 高，但 5 次失败全是已被 syndrome check 检出的失败，未证明对漏检逻辑错误有额外信息。移植前需用真实配对开发数据检验 BP 与 syndrome-consistent 级联/OSD-0 候选的物理位分歧，是否能区分“更强搜索可救回”的帧；全矩阵 OSD-0 当前依赖稠密消元，必须先解决每帧预处理成本，才能进行合理的路由复杂度对照。


## 2026-11-03 文献补充：逻辑陪集覆盖与退化概率估计

1. Shang, *Logical-Space Theory of Ordered-Statistics Decoding for Quantum LDPC Codes*, [arXiv:2610.24503](https://arxiv.org/abs/2610.24503)（2026 预印本）。它把 OSD 的逻辑失败严格拆成两部分：候选列表没有覆盖真实逻辑签名（coverage miss），以及覆盖了真实签名但成本准则选错类（selection failure）。论文用固定 RREF 的自由坐标到逻辑签名线性映射，构造 coverage-aware 的权重二候选表，并在其二进制 BB/HGP 实例上分别验证候选覆盖和类代表选择的影响。作者明确限定了实验模型和代码族；这不是 GF(256) 级联码上的现成算法或 FER 结论。
2. 与本项目的直接对应：`cascadeZCosetSignature` / `cascadeZCosetKey` 提供级联码 Z 稳定子商空间的精确标签；已有 OSD 审计曾发现困难帧候选列表不含真实类，且许多 syndrome-valid 词属于不同类。因此下一步应先建立“GF(256) 自由符号扰动 → 物理二进制 Pauli 差分 → 级联逻辑签名”的可验证映射，统计固定候选预算下不同列表的真实类覆盖率，再设计不使用真值的列表与类选择。二进制论文的公式不能直接照搬：GF(256) 一个自由符号有 8 个二进制自由位，场基和物理缩放会影响扰动映射。
3. Krishnamoorthy et al., *Certified decoding of quantum LDPC codes*, [arXiv:2608.25545](https://arxiv.org/abs/2608.25545)（2026 预印本）。其关键点是退化 ML 比较的是逻辑类的配分函数（同一逻辑类下全部 stabilizer-equivalent 物理错误概率之和），不是单个最可能错误词；作者用 coset MRF、Bethe 近似与 AIS/common-random-numbers 估计类概率，并在部分 BB 实例上报告验证。此框架解释了为什么本项目“找到任意 syndrome-valid 候选”不够。它仍需改造成 GF(256) 级联矩阵与当前条件 Pauli 信道，并评估 17k 左右逻辑维数下的候选类数及运行成本；先做小型精确实例和稳定子变量局部因子验证，不能直接声称能达到其结果。
4. 当前实验决策：`joint_projection_safe` 将 10/100 个 Z syndrome failures 全变成 10/100 undetected logical failures、总 FER 不变。这正是 syndrome feasibility 与 logical-class selection 的区别。先停止继续放宽 syndrome-valid 候选接受条件；将下一项开发任务改为逻辑陪集覆盖/类概率诊断，在线选择始终只依赖 syndrome、矩阵、信道模型和 decoder 输出，真实错误仅用于离线覆盖审计。

## 2026-11-03：Restart-Belief 与 syndrome-flip 配对进展

### Tuned + global Restart-Belief

新增的 `conditioned_robust_tuned_restart_belief_global` 已接入与 tuned 相同的 `10` Turbo rounds、`20` outer iterations、mixing `0.2` 初始译码，并只在 Z syndrome 不通过时执行全局低 LLR bit 重启。编译与单帧 smoke 检查通过。

在固定内长 12 path 矩阵、`p_D=0.07125`、新噪声 seeds `20320000–20320003` 的 20 帧配对 pilot 上，tuned 与新模式均 `1/20`，救回/退步 `0/0`、精确 McNemar `p=1.0`；平均 CPU `10.06→14.13 s/frame`。预选难帧重放中 4 分支×4 步、16 次额外 BP 没找到 syndrome-valid 候选。该版本没有表现出逻辑 FER 收益，不扩大。

结果目录：`results/paired_tuned_restart_global_pilot20_seed20320000/`。

### Reliability-guided GF(256) syndrome-flip

参考 Wang, Li, Mueller, *Fully Parallelized BP Decoding for Quantum LDPC Codes Can Outperform BP-OSD*, [arXiv:2507.00254](https://arxiv.org/abs/2507.00254)。论文的 BP-SF 用 BP 迭代中的位翻转频率识别候选位，枚举低重量位型 `t`，对修改综合征 `s' = s + Ht` 做试探 BP，再将 `t` 加回候选。原文主要在二进制 BB/QCode 实例验证，因此这里实现的是 GF(256) 级联适配原型，不是原文算法复现：目前按 Turbo BP 的物理位后验不确定度排序，枚举最不确定 4 个物理位中重量 1 和 2 的模式。

映射使用物理位到域符号的 `FieldBasis::fromDual` 及每列 `inner_scales`。令 `t` 是相应 GF(256) 坐标差分，采用特征 2 的恒等式：

\[
H_X(e'+t)=(s+H_Xt)+H_Xt=s.
\]

每个 shifted-Turbo 候选加回 `t` 后仍需通过原始 GF(256) Z syndrome；在线候选打分用当前 X 估计条件下的物理 Pauli 信道 log-likelihood，不使用真实错误。模式名为 `conditioned_robust_tuned_syndrome_flip`，参数由 `GF256_CASCADE_SF_BITS` 与 `GF256_CASCADE_SF_MAX_WEIGHT` 控制。

在预选基线失败批次 seed `20268300` 的 5 帧上，重量≤1 的前 4 位没有有效候选；改为重量≤2 后有 1 个 syndrome-valid 候选，结果从 tuned `1/5` 变为 `0/5`。这是按失败挑选的开发数据，只验证映射/分支，不计作 FER 证据。

固定 4 位、最大重量 2 后，新 pilot seeds `20322000–20322003` 的 20 帧两侧均 `0/20`，没有失败可比较。另在独立 seeds `20323000–20323019` 的 100 帧同帧配对验证中：

| 解码器 | 逻辑失败 | Wilson 95% 区间 | Z syndrome failure | undetected logical | CPU / frame |
|---|---:|---:|---:|---:|---:|
| `conditioned_robust_tuned` | 7/100 | `[0.0343, 0.1375]` | 7 | 0 | 9.76 s |
| tuned + syndrome-flip | 5/100 | `[0.0215, 0.1118]` | 5 | 0 | 14.93 s |

配对结果救回 2、退步 0，精确 McNemar 双侧 `p=0.5`；均为 Z syndrome failure 的救回，没有漏检逻辑失败。点估计 FER 从 0.07 到 0.05 是探索性正向信号，但事件太少，置信区间重叠，且平均 CPU 增加约 53%，不能称为可重复或统计显著的性能改善。

结果目录：`results/paired_tuned_syndrome_flip_pilot20_seed20322000/`、`results/paired_tuned_syndrome_flip_holdout100_seed20323000/`。

第一组独立 300 帧复制在固定代码、矩阵及预算下完成，种子 `20324000–20324059`：

| 解码器 | 逻辑失败 | Wilson 95% 区间 | Z syndrome failure | undetected logical | CPU / frame |
|---|---:|---:|---:|---:|---:|
| `conditioned_robust_tuned` | 15/300 | `[0.0305, 0.0808]` | 15 | 0 | 9.84 s |
| tuned + syndrome-flip | 13/300 | `[0.0255, 0.0727]` | 13 | 0 | 13.67 s |

配对救回 2、退步 0，精确 McNemar 双侧 `p=0.5`。300 帧中 syndrome-flip fallback 触发 15 次，共评估 150 个扰动模式，其中 10 个产生原 syndrome-valid 候选，2 个最终消除逻辑失败，另 13 帧仍失败。平均 CPU 增加约 39%；失败批次中观测到的每个 5 帧 batch CPU 最高约 236 秒。结合此前独立 100 帧结果，两组都表现为少量救回，但每组的配对证据均不显著。由于方法只在 syndrome-failure 帧触发，它不能改变基线 syndrome-valid 帧结果；重点是它救回失败帧的频率是否足以抵消显著时延。

将这两组互不重叠的独立样本作描述性合并：tuned 为 `22/400=5.5%`，syndrome-flip 为 `18/400=4.5%`；配对救回 4、退步 0，精确 McNemar 双侧 `p=0.125`。合并分析仍未达到常用 0.05 显著性阈值，不可据此声称已有统计显著改进。

第二组独立 300 帧复制已完成，种子 `20325000–20325059`：tuned `14/300`（Wilson 95% `[0.0280,0.0768]`），syndrome-flip `13/300`（`[0.0255,0.0727]`）；配对救回 1、退步 0，精确 McNemar 双侧 `p=1.0`。Z syndrome failure 从 14 降到 13，漏检逻辑失败两边均为 0；CPU `10.01→13.65 s/frame`，墙钟约 31.1 分钟。

三组互不重叠的独立样本合并为描述性统计：tuned `36/700=5.14%`，syndrome-flip `31/700=4.43%`；救回 5、退步 0，精确 McNemar 双侧 `p=0.0625`。改善方向一致，但当前事件数仍未达到 `p<0.05`，且 decoder 平均耗时有明显增加，不能宣称 FER 已显著改善或满足论文结论。方法仅救回少数基线 syndrome failure，当前保留为可选 fallback，后续优先投入更合适的级联投影/内码因子译码。

完整结果：`results/paired_tuned_syndrome_flip_holdout100_seed20323000/`、`results/paired_tuned_syndrome_flip_validation300_seed20324000/`、`results/paired_tuned_syndrome_flip_replication300_seed20325000/`。


### 级联 GF(256) 下的逻辑类配分函数（待验证的推导）

设 Z 错误估计为 `r∈GF(256)^N`，满足 `H_X r=s`，Z 稳定子生成矩阵为 `H_Z`，且 `H_X H_Z^T=0`。它的逻辑等价类为 `r+row(H_Z)`。若先取 `H_Z` 的满行秩生成基，令稳定子系数 `u∈GF(256)^{rank(H_Z)}`，每个物理 GF 符号位置的条件先验为 `π_i(·)`（由当前 X 估计、物理位序/基变换和信道模型得到），则候选逻辑类的非归一化概率可写成

\[
Z_{[r]}=\sum_{u\in GF(256)^{\operatorname{rank}(H_Z)}}
\prod_{i=1}^{N}\pi_i\!\left(r_i+\sum_a u_a(H_Z)_{ai}\right).
\]

行秩基保证没有稳定子参数冗余；即使保留冗余生成行，各逻辑类会乘上相同常数，不影响比较。该式把类内所有稳定子等价代表的概率合并起来，对应二进制文献的 coset partition function。它可表示为以稳定子系数 `u_a` 为 GF(256) 变量、以物理符号 `i` 为局部因子的因子图；每个因子只依赖 `H_Z` 第 `i` 列中的非零稳定子系数。CSS 正交性保证任意 `u` 都不改变 `H_X` syndrome。

这只是代数上的候选路线。GF(8) toy 核对已完成：`test_coset_partition_identity.py` 对 `H_X=H_Z=[1,1,0]` 的 64 个 syndrome/logical 类逐一比较直接 512 个错误向量的分组和与 8 个稳定子参数求和，最大质量差 `1.388e-17`；并对三种不同 field scaling、固定 X 位型和退极化信道 `P(Z|X)` 构造物理位先验。另外，C++ `test_soft_x_prior.cpp` 使用项目 `FieldBasis`、Pauli `P(Z|X)` 和 3 个非零尺度逐项对照 Turbo 内部 Z 变量与物理位先验，最大误差为 0；既有 GF(8) 后验测试最大误差 `1.11e-16`。随后新增 `test_cascade_coset_invariants.cpp`，已在两套实际 GF(256) path 矩阵上逐行验证 HZ 生成元与 HX 正交、稳定子成员判定及陪集键不变，并用一个非平凡逻辑见证检查分类；细节见 `results/coset_invariant_validation/validation.log`。这些结果验证了矩阵级结构不变量，但没有验证 GF(256) 类概率近似或译码 FER。下一步应对少量困难帧的候选逻辑类实现低成本 Bethe/局部枚举近似，与可精确得到的类分组核验误差，再评估复杂度。若近似只拟合单个词概率或产生无校准排序，就不能用作在线选择。

## 2026-11-15 PEG H2 与通用 G2 Turbo 译码诊断

构造器新增 `peg3` 内码：GF(256) 上 `H2=3×12`、`G2=9×12`，内码 Tanner girth 为 6，构造校验 `H2 G2^T=0`、CSS 正交及二进制展开正交均通过。矩阵为 `HX=1746×6984`、`HZ=3492×6984`，推得量子率 `1/4`。既有单行 path 内码是 `H2=1×12`，`HX=2134×6984`、`HZ=2716×6984`，率 `11/36`，所以这是不同码率的结构试验，不能当同码率 H2 消融。

更正早期 selector 解释：定理 3.1 的内校验矩阵写成 `H_C1=(P1,P2)`，并要求 `P2` 可逆；这是分阶段 Z 译码在外层恢复前 `k1` 个坐标后，用内综合征求出遗漏的 `r1` 个坐标所需的条件。它不等价于“任意零空间生成基 `G2` 必须逐项为 `[I_k|A]`”。若只对同一 `H2` 的零空间基作可逆行变换，内码零空间相同；Kronecker 生成矩阵的行空间也相同。不同基仍可能改变 BP 图的消息路径和外层扩域符号的坐标约定，所以实现要核实坐标映射，但不能将换基说成修复 CSS 构造错误。[论文定理与译码段落](https://opus.lib.uts.edu.au/bitstream/10453/167215/2/2012.00226v2.pdf)给出 `H_C1=(P1,P2)`、`P2` 可逆及恢复公式。

首次生成时复用了 `computeGeneratorMatrix` 的稀疏主元策略，得到的是同一 `H2` 零空间的另一组基。将 `G2` 换成系统形基会改变 Tanner 图中的行表示，但若行空间相同，不能据此称旧内码或 CSS 矩阵构造错误。论文式分阶段译码需要的是 `H2=[P1|P2]` 且右侧方块 `P2` 可逆，以便从已译出的前 `k` 个坐标和内综合征恢复剩余坐标；这和强制 `G2` 本身逐项等于 `[I_k|A]` 是不同陈述。需要另行核对的是：`G2` 的行基/符号坐标如何与外层扩域符号映射配套，以及当前 generic Turbo 是否按实际系数传递消息。

另构造了显式系统形零空间基：对 `H2=[A|I_3]` 直接生成 `G2=[I_9|A^T]`，并断言 `H2G2^T=0`。矩阵保存在 `results/h2_peg3_n12_systematic/build/`，没有覆盖旧基。它的 GF CSS 正交和二进制展开正交均通过，H2 girth=6；尺寸为 `HX=1746×6984`、`HZ=3492×6984`，推得率 `1/4`。该操作是方便坐标解释和构造特定 factor graph 的换基，不是已证实的码构造修复。旧 path 率 `11/36`，不能与 PEG 结果直接作同码率比较。

为支持 `HX=H_outer⊗G2` 的一般 `G2`，`decodeHxCascadeTurbo` 新增通用内因子图：每个物理块使用局部约束 `G2 e + y=0`，外层投影与内层物理符号通过 GF(256) extrinsic 概率迭代。通用分支 `GF256_CASCADE_GENERIC_INNER_ITERS` 默认 5，取值范围 `[1,50]`。代码能编译。

在 `p_D=0.07125`、仿真主种子 `20272000`、20 帧、`conditioned_robust_tuned` 下做了同序列小批诊断：

| 内码/硬判决方式 | 逻辑失败 | syndrome failure | 漏检逻辑失败 | 说明 |
|---|---:|---:|---:|---|
| 旧 path `H2=1×12` | 2/20 | Z: 2/20 | 0 | 17 帧 exact recovery；旧内码基线 |
| 原始零空间基 PEG，Turbo APP 各符号独立硬判决 | 20/20 | Z: 20/20 | 0 | 通用内因子图未能收敛到 syndrome-valid 输出；不能由该结果单独归因于坐标基 |
| 原始非系统形 PEG，外层 syndrome 投影 + 内层精确线性恢复 | 20/20 | Z: 0/20 | 20 | 全部满足 syndrome，却都在错误逻辑类；仅诊断了投影硬判决策略 |

结果文件在 `results/h2_peg3_n12/projected_recovery_pilot20/`：`path_baseline_sameframes.tsv/log`、`peg3_original_default.tsv/log`、`peg3.tsv/log`。三组均报告均值物理错误数 `3959.9`，支持输入噪声序列一致。旧基与系统形基来自同一个 `H2` 零空间，故它们不能作为不同内码/不同距离的比较；旧 path 又有不同 `H2`、码率，亦不是控制变量充分的 PEG 消融。投影硬恢复只作为机制诊断，不作为默认行为；由 `GF256_CASCADE_GENERIC_EXACT_RECOVERY=1` 显式启用，默认保留原 Turbo APP 硬判决。

换成系统形 G2 后，以同信道点和同主种子对前三帧做了机制 smoke：`results/h2_peg3_n12_systematic/pilot20/peg3_systematic.log` 中 3/3 帧仍有 Z syndrome failure（分别 1053、512、274 个未满足校验）。原计划 20 帧因早期连续失败而中断；没有产生完整 TSV，因此不能报告这个表示下的 FER。显式开启精确投影跑同一帧后，`exact_recovery_smoke.tsv/log` 显示 syndrome 通过但逻辑仍失败。该结果说明单纯换成系统形基不足以修复当前译码器，也未解决逻辑类判决；它不证明前一基在构造上错误。

之后用同一仿真器、`p_D=0.07125`、主种子 `20272000` 配对重跑前三帧：旧基/系统形基两侧均 `3/3` syndrome 和逻辑失败，均值物理错误数同为 `3955`；每帧未满足校验数旧基为 `1016,1076,222`，系统形为 `1259,1205,441`。这进一步表明换基会改变 loopy BP 的残差轨迹，但此样本不足以比较逻辑 FER。可复核的 TSV、逐帧日志和矩阵/二进制哈希在 `results/peg_g2_basis_same3_seed20272000/`；不将其解释成系统形或非系统形基的性能胜负。

结论：目前没有证据表明多行 PEG H2 改善译码。多行内码的 generic Turbo 仍未得到有效 syndrome 输出；“系统形基”前三帧也全部 syndrome failure，故应优先核对一般 `G2` factor-message、外层符号映射和坐标顺序，而不是继续把基对齐当作修复。GF(8) 树形 toy 的 24 个精确先验/syndrome 组合现已通过；接下来应检查 GF(256) 上真实失败帧的变量顺序、局部消息残差和候选逻辑类覆盖，再设计同码率、同码空间的配对 FER 比较。旧候选上的精确投影可消除 syndrome failure，却把失败全部变成漏检逻辑失败；所有候选仍须按稳定子等价逻辑 FER 评价。

### Generic G2 因子图的精确 toy 对照（2026-11-15）

新增 `test_generic_g2_turbo_exact.cpp`：构造 GF(8) 的树形级联因子图，`G2=[1, α, α²]`，外层只有一个双变量校验；用 3 组先验（无 X 条件先验及两组不同 X 估计）与全部 8 个 syndrome 形成 24 个条件。每例独立枚举 `8^6` 物理错误，直接计算每个物理符号的后验边缘，并与 generic Turbo APP 对照。最大绝对边缘误差 `7.30638e-13`，显式 X→Z 先验与从 `x_estimate` 内部生成的先验差为 `0`。测试退出条件同时约束这两项误差，构建与执行通过。它只验证这个无环 toy 图上的因子消息和先验坐标映射；不能证明真实 GF(256) 多环图会收敛，也不能解释其大量失败。

### 论文构造边界及 OSD syndrome-pass 负结果（2026-11-15）

参考 Fan 等人的构造时需区分域与 Kronecker 记号：论文 II-C/定理证明把内码 `C1` 的矩阵放在 `GF(q)`，外码 `C2/C3` 放在 `GF(q^{k1})`，并将扩域系数用伴随矩阵形式展开；式 (41) 给出级联校验的分阶段结构。参见 [arXiv HTML，II-C 与定理证明](https://arxiv.org/html/2012.00226)。当前 `construct.cpp` 则将 TIT 外矩阵和 `G2/H2` 直接作为 GF(256) 元素，用普通标量 GF(256) Kronecker 积生成矩阵。当前矩阵已经通过本项目的 CSS 正交核验，但它不是上述论文构造的逐项复现；论文定理的参数、距离和分阶段译码性能不可直接转移到当前模型。若要复现论文，应单独实现域嵌入/伴随矩阵展开、坐标映射和相应分阶段译码，并重新验证正交性与逻辑 FER。

对系统形 PEG 候选单帧运行 `conditioned_robust_osd`：`p_D=0.07125`、矩阵 seed `20261018`、噪声 seed `20272000`。OSD 检查 65,281 个候选并得到完整 syndrome-valid 候选，但量子逻辑仍失败；该帧被记为 undetected logical failure，残余差错有 2,879 个非零符号，运行时间约 45 秒。配置、结果和日志保存在 `results/h2_peg3_n12_systematic/osd_pilot1_seed20272000/`。这只是一个机制反例，不估计 FER；说明将 OSD 从 projection failure 中恢复 syndrome 仍不能替代逻辑陪集判决，而且该候选成本较高。下一步应优先做离线候选类覆盖率审计，在线选择不得读取仿真真值；若无法在合理候选预算下覆盖正确逻辑类，应先修正构造/坐标或软信息，而不是继续堆叠 syndrome 修复。

### generic-G2 内因子迭代次数扫描（2026-11-15）

为确认多行 PEG 内码失败是否只是局部内因子迭代次数不足，在相同 systematic PEG 矩阵、相同噪声 seed `20272000` 的第 0 帧上固定 `p_D=0.07125`、generic damping `0.9`，扫描 `GF256_CASCADE_GENERIC_INNER_ITERS=1,5,15,30`。未满足 Z 校验数分别为 `1549,912,910,910`，报告残余符号数分别为 `2000,1306,1313,1313`；四种设置均为 syndrome failure 和 logical failure。5 次比 1 次显著减少该帧残差，15 至 30 次则停在相同平台，未解决失败。逐次日志和命令参数摘要见 `results/peg3_generic_inner_iters_single_frame_seed20272000/`。这是同一帧上的机制扫描，不是 FER 实验，也不应据此把 damping 或迭代数作为已优化参数。下一步应记录每轮 inner/outer extrinsic 变化、熵及 syndrome 残差，以区分未收敛与错误固定点；若残差稳定而 syndrome 仍失败，再考虑协调多个外层投影或局部精确更新。

### PEG 多行内码 layered 外层 schedule 配对 pilot（2026-11-15）

为检验更换外层更新顺序能否缓解系统形 PEG `G2` 的 syndrome 残差，在完全相同矩阵和噪声帧上比较 `conditioned_robust_tuned` 与 `conditioned_robust_layered`。矩阵目录 `results/h2_peg3_n12_systematic/`，其 `H2=3×12`、`G2=9×12`；`p_D=0.07125`（`f_m=0.0475`），局部内迭代 5、generic damping 0.9，噪声 seeds `20330000–20330003`，共 20 帧。实际物理错误均值为 3998/帧，确认配对信道输入相同。

| 译码器 | 逻辑失败 | Z syndrome failure | 漏检逻辑失败 | 平均未满足 Z 校验 | 平均残余符号 | CPU 秒/帧 |
|---|---:|---:|---:|---:|---:|---:|
| tuned | 20/20 | 20/20 | 0 | 490.2 | 724.6 | 29.67 |
| layered | 20/20 | 20/20 | 0 | 472.2 | 685.2 | 29.89 |

所有帧两种译码器都发生相同类别的失败，没有配对救回或退步（McNemar `p=1`）。layered 的平均残差略低，但没有一个 syndrome 通过，因而这不是逻辑 FER 改善。样本只说明在此码/信道点上 schedule 替换无效；20/20 失败对应 Wilson 95% 区间 `[0.839,1]`，该高失败率不能外推到更低 `p_D`。结果及逐帧日志保存在 `results/peg3_layered_vs_tuned_pilot20_seed20330000/`。

文献核对也表明，Fan et al. 的 AQCTPC 译码不是把完整短环图交给 BP：论文 Algorithm 2 先按内校验矩阵逐子块取得 syndrome 并选取内码代表，再把修正后的外层 syndrome 映射到扩域并对 `C_3` 做 syndrome decoding，最后用 `H_{c_1}=(P_1,P_2)` 中可逆的 `P_2` 解出遗漏坐标；见 [论文原文的 III 节及 Algorithm 2](https://arxiv.org/abs/2012.00226)。但这不能直接替换当前 generic Turbo：该论文的级联码使用 `GF(q)` 内码与 `GF(q^{k_1})` 外码及其扩域展开，而项目当前的 `construct.cpp` 是同一个 `GF(256)` 上普通 Kronecker 积。要移植其分阶段算法，必须先证明项目矩阵具备对应的内码/扩域坐标分解，不能只照搬步骤或把下半部 HZ 校验当作 Z 错误的额外 syndrome。

下一步把实验移到 PEG 候选上两种译码器都不饱和的较低 `p_D`，同时记录 outer projection 与真实投影之间的符号错误数；先看是否有可恢复的逻辑 FER 差异，再决定是否投入论文式跨域构造。不可用 `p=1` 失败平台上的低校验残差代替成功率改进。

### 多行 PEG 内码的稳定子判定修复与同帧复测（2026-11-15）

在比较多行 PEG `H2=3×12` 时发现逻辑评估器存在结构限制：旧 `isCascadeZStabilizer` 对一般内码默认只有一个 parity 坐标（`n-k=1`），而 `cascadeZCosetSignature` 也只把 `word[j*n2+k2]` 当成唯一尾符号。PEG 候选有三个 parity 符号，导致 syndrome 已通过、但差错其实是 `H2` 行空间内的帧被误报为 logical failure；基于旧 signature 的候选分组也不可信。之前的 PEG syndrome-failure 数仍然有效，但所有 syndrome-pass logical 标签及 FER 都需按修正实现重新计算。

当前矩阵由构造器明确保证系统形 `H2=[A|I_3]`、`G2=[I_9|B]`。Z 稳定子等价判定对每个内块用尾部三坐标消去 `H2` 行空间，随后检查九个信息坐标跨外块向量是否属于 `HX_TIT_Outer_For_HZ` 的行空间；Z coset signature 使用同一消元后的坐标再作外层约化。X 稳定子检查验证每块完整向量属于 `G2` 行空间并检查九个外层系数向量。当前一般多行分支遇到不符合这两个系统形条件的矩阵会报错，不会静默套用单校验近似。

实现位置：`simulate.cpp` 的 `FieldRowSpace` 现在可持有内层 `H2/G2`；generic systematic 分支用于 `isCascadeXStabilizer`、`isCascadeZStabilizer` 和 `cascadeZCosetSignature`。实际 PEG 矩阵的 `test_gf256_cascade_coset_invariants` 结果：`H_Z` 生成行与 `H_X` 的正交错误 `0/3492`，被 Z 稳定子判定拒绝 `0/3492`，改变 Z 陪集键 `0/3492`；一个非平凡逻辑见证的 syndrome、稳定子判定和陪集键三者一致。原 1×12 path 矩阵也重新执行了该测试，`0/2716` 异常，旧路径未受影响。

用修正后二进制同帧实验复核 `p_D=0.04`、seed `20330200–20330207`，相同 8 个物理错误帧下：旧 path tuned `0/8` 逻辑失败；PEG 3×12 tuned `3/8`，全部是 Z syndrome failure，漏检逻辑失败 `0`。PEG 旧结果曾报 `6/8`；重新核验后其中 3 帧其实是稳定子等价恢复。两矩阵的平均物理错误数均为 `2221.125`。配对结果是旧 path-only 成功 3 帧、反向 0 帧，McNemar `p=0.25`；这个小样本且码率不同的构造比较只作诊断，不能证明 path 普遍优于 PEG。修正后 PEG 同一批 8 帧的 tuned 与 layered 也都是 `3/8`，失败集合完全相同、McNemar `p=1`，CPU 分别 `9.42` 与 `9.77 s/frame`。因此 schedule 改动没有给出可见收益。

可复核文件：`results/path_vs_peg_same_frames_pD0p04_seed20330200_corrected_eval/`、`results/peg3_layered_vs_tuned_corrected_eval_pilot8_pD0p04_seed20330200/`、旧 path 的 `results/g2_basis_path_n12_300frames/path/coset_invariants_corrected.log`、PEG 的 `results/h2_peg3_n12_systematic/coset_invariants_corrected.log`，以及对应的 `test_gf256_cascade_coset_invariants` 输出。旧多行 PEG 实验中的 syndrome-pass 逻辑标签已作废，待重跑的还有 OSD/Gibbs/逻辑列表等用到旧陪集 signature 的历史 PEG 结果。

修复评估器后，进一步在 `p_D=0.04`（`f_m=0.02667`）对两套矩阵跑完全相同的 20 个通道帧（seed `20330400–20330419`），均用 `conditioned_robust_tuned`：

| 矩阵 | 量子率 | 逻辑失败 | Z syndrome failure | 漏检逻辑失败 | HX 符号图 4 环 |
|---|---:|---:|---:|---:|---:|
| path `H2=1×12` | `11/36` | 0/20 | 0 | 0 | 35,502 |
| PEG `H2=3×12` | `1/4` | 18/20 | 18 | 0 | 61,110 |

两边平均物理错误数都是 `2237.15`。seed `20330400` 的 `.truth.tsv` 导出 SHA256 完全一致（`85dd2575…b7ef61b5`），逐位确认 X/Z 通道噪声相同。配对 discordant frames 为 path 成功而 PEG 失败 18 帧、反向 0 帧；精确 McNemar 双侧 `p=7.63×10^{-6}`。这是对当前固定码和噪声点很强的负面证据，但码率、内码行数和 HX 四环数都变了，只支持“这个 PEG3×12 候选及其当前译码器组合不合适”，不证明 PEG 构造普遍较差。结果：`results/path_same_frames_as_peg20_pD0p04_seed20330400/`（path）和 `results/peg3_syndrome_flip_corrected_eval_pilot20_pD0p04_seed20330400/`（PEG tuned 数据）。

在 PEG 同一矩阵、相同 20 帧上再比较 tuned 与 tuned + syndrome-flip：两者均 `18/20`，失败集完全相同、全为 Z syndrome failure；配对救回/退步 `0/0`、McNemar `p=1`。平均 CPU 从 `12.77` 增到 `39.92 s/frame`，因此这次后处理没有帮助该失效模式。PEG layered schedule 在另一组修正评估的 8 帧 pilot 也是 `3/8` 对 `3/8`，未见改善。该组完整路径对照、同帧噪声哈希、码率及短环数同时列于上表；这些结果现已同步到根目录 `TASKS.md`。

下一步不应只按 H2 Tanner girth 挑候选。先在同一码率、固定外层矩阵和正确的多行稳定子评估下生成/筛选多个 H2；同步报告 `G2` 的结构、全局 `H_X` 四环和小重量逻辑算符，再对有希望的候选做固定参数的配对 FER。只有 syndrome failure 降低且逻辑类正确率改善，才扩大到独立数百/数千帧。

### path 矩阵上的独立 tuned/adaptive 复验与失败帧回放（2026-11-15）

为检查此前 50 帧自适应译码器改善能否复现，在固定 GF(256)、内长 12 path 矩阵上使用新种子 `20350000–20350019`，`p_D=0.07125`（`f_m=0.0475`），逐帧配对比较 `conditioned_robust_tuned` 与 `conditioned_robust_adaptive`，每模式 100 帧。矩阵哈希为 `HX=cc8fe822…c2b858c`、`HZ=c626eecd…e5dcbb8`；可执行文件 SHA256 `d924600a…28d0986`。

| 方法 | 逻辑失败 | Z syndrome failure | 漏检逻辑失败 | Wilson 95% FER 区间 | CPU 秒/帧 |
|---|---:|---:|---:|---:|---:|
| tuned | 6/100 | 6 | 0 | `[0.0278, 0.1248]` | 9.49 |
| adaptive | 6/100 | 6 | 0 | `[0.0278, 0.1248]` | 11.55 |

失败帧完全相同，配对救回/退步均为 0，McNemar `p=1`。因此此前 50 帧上的 adaptive 改善未在这个独立样本复现，且平均 CPU 时间增加约 22%。结果目录：`results/path_tuned_vs_adaptive_holdout100_pD0p07125_seed20350000/`。

进一步对该组 tuned 的 6 个 Z syndrome 失败帧（`20350001/t3`, `20350003/t1`, `20350007/t0`, `20350011/t4`, `20350012/t4`, `20350019/t2`）做目标帧回放。该回放是按失败筛选的机制诊断，不能直接作为 FER 样本：

| 回放方法 | 逻辑失败 | Z syndrome failure | 漏检逻辑失败 | 六帧总 CPU |
|---|---:|---:|---:|---:|
| `conditioned_robust_long` | 6/6 | 4/6 | 2/6 | 171.7 s |
| `conditioned_robust_tuned_syndrome_flip` | 4/6 | 4/6 | 0/6 | 465.6 s |

long 在两帧上把 syndrome 修到通过，但仍选择了错误逻辑陪集；这再次说明 syndrome-valid 不等于量子逻辑成功。syndrome-flip 在 6 个预选失败帧中救回 2 帧且没有退步，其余 4 帧仍 Z syndrome failure。代码检查确认 syndrome-flip 与 tuned 使用相同初始 Turbo 参数，且只在 Z syndrome 不通过时进入后处理；所以将其应用于这组 100 帧可精确推出 `4/100` 对 `6/100`，救回 2、退步 0，McNemar `p=0.5`。这是一个独立样本上的正向点估计，但统计证据不足，且失败帧回退成本高，不能据此宣称已改善 FER。

逐帧输出和参数保存在 `results/path_failed6_long_vs_syndrome_flip_diagnostic/`。下一步应冻结该 syndrome-flip 配置，在另一个独立种子区间做确认性配对实验，并按完整帧数报告平均 CPU；不应继续在这 6 个预选失败帧上挑参数。

### syndrome-flip 独立 300 帧确认（2026-11-15）

固定参数，不在本批数据上调参；使用相同 path 矩阵、`p_D=0.07125`（`f_m=0.0475`），新种子 `20351000–20351059`，tuned 对 tuned+syndrome-flip，各 300 帧，`batch-size=5`、4 workers。矩阵和可执行文件哈希与上面的 100 帧实验一致。

| 译码器 | 逻辑 FER | Wilson 95% 区间 | Z syndrome failures | 漏检逻辑失败 | CPU 秒/帧 |
|---|---:|---:|---:|---:|---:|
| tuned | `9/300=3.00%` | `[1.59%,5.60%]` | 9 | 0 | 8.90 |
| tuned + syndrome-flip | `7/300=2.33%` | `[1.13%,4.74%]` | 7 | 0 | 10.72 |

逐帧配对结果：基线失败、syndrome-flip 成功 2 帧；基线成功、syndrome-flip 失败 0 帧；精确 McNemar 双侧 `p=0.5`。平均 CPU 增加约 20.5%。点估计方向与前一组 100 帧同向，但这批单独未达到显著性。将两组不重叠种子作描述性合并为 tuned `15/400`、syndrome-flip `11/400`，总救回 4、退步 0、McNemar `p=0.125`，仍不够支持稳定 FER 改进。

结论：syndrome-flip 有初步、方向一致的救援信号，但目前统计证据弱，且计算成本增加。下一步冻结设置，在新独立种子上继续配对验证；结果目录 `results/path_tuned_vs_syndrome_flip_holdout300_pD0p07125_seed20351000/`，包括 `summary.csv`、`pairwise.csv`、`failure_transitions.csv`、逐批日志和输入哈希。

### syndrome-flip 新独立 1000 帧配对确认（2026-11-15）

在不改变 decoder 参数的条件下，用新的噪声种子 `20352000–20352199` 做 1000 帧配对实验；矩阵、信道点、可执行文件哈希和上一节一致，固定 `p_D=0.07125`、`f_m=0.0475`、8 个最不确定物理位、只枚举重量 1 的 syndrome-flip 候选。

| 译码器 | 逻辑失败/FER | Wilson 95% 区间 | Z syndrome failure | 漏检逻辑失败 | CPU 秒/帧 |
|---|---:|---:|---:|---:|---:|
| tuned | `41/1000=4.10%` | `[3.04%,5.51%]` | 38 | 3 | 13.37 |
| tuned + syndrome-flip | `35/1000=3.50%` | `[2.53%,4.83%]` | 32 | 3 | 16.73 |

逐帧转移为：959 帧两者成功，3 帧两者均漏检逻辑失败，32 帧两者均 Z syndrome failure，另有 6 帧从基线 Z syndrome failure 变成 syndrome-flip 成功；没有成功退化成失败的帧。精确 McNemar 双侧 `p=0.03125`。这构成固定矩阵、固定参数、固定信道点上的正向配对证据，点估计 FER 降低 0.6 个百分点（相对降低约 14.6%），但 CPU 增加约 25.2%。

结合此前独立 300 帧，方向相同（该组单独为 `9/300→7/300`、救回 2、退步 0、`p=0.5`）；100 帧筛选失败回放亦有 2 帧救援。应把 1000 帧段作为主要确认结果，其他批次视为补充；不把所有连续观察批次的 p 值合并成替代确认性检验。当前结论仍只覆盖一个 path 矩阵和一个 `p_D`，并且每次 fallback 的失败帧耗时可能很高。下一步应在独立矩阵种子上复验，随后在预定的邻近信道点检验稳定性和复杂度；不能据此宣称已达到发表标准。

结果目录：`results/path_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20352000/`；关键汇总文件为 `summary.csv`、`pairwise.csv`、`failure_transitions.csv` 与 `config.json`，逐批 TSV/log 保留 1000 帧的完整运行记录。


### 第二矩阵1000帧补完及评估器修复（2026-11-16）

矩阵seed=20261019，GF(256)，内长12，p_D=0.07125（f_m=0.0475），噪声seeds=20353000–20353199。原有193个完整配对批次共965帧；另7对批次在终端逻辑评估时报错，因此中间965帧不是完整随机样本，不能用于提前结论。

恢复原可执行文件后SHA256与配置完全一致（d924600a3cfd070ce9214537321d222781375de0304d635c78ef58d2a28d0986），重跑复现异常：`Generic cascade X-stabilizer check requires G2=[I|A]`。原因是加载了显式path G2的单校验内码也进入了为多行H2添加的系统形评估分支。现将该分支限定于k2<n2-1，单校验情形继续使用由H2缩放确定、与G2行基无关的原有检查。译码、噪声和候选选择均未改动。

实际第二矩阵核验：2134条HX生成行及零向量均被X稳定子检查接受；2716条HZ生成行正交、Z稳定子成员和陪集键均无异常；非平凡逻辑见证分类一致。此前完成的386个批次逐一满足ExactRecoveryFrames+SyndromeFailures=Total，故从未进入修复的终端分支，可复用。只重跑缺失的14个批次（每方法35帧）；保留原metadata，不伪造新的程序哈希。`completion_provenance.json`记录旧/新程序、复用依据及补跑任务，`reproducibility/`保存程序和修复后源码。

| 方法 | 逻辑失败/FER | Wilson 95%区间 | Z综合征失败 | 漏检逻辑失败 | CPU秒/帧 |
|---|---:|---:|---:|---:|---:|
| tuned | 48/1000=4.8% | 3.64%–6.31% | 47 | 1 | 13.26 |
| tuned+syndrome-flip | 36/1000=3.6% | 2.61%–4.94% | 35 | 1 | 17.26 |

逐帧配对：952帧均成功，35帧均Z综合征失败，1帧均漏检逻辑失败，12帧由Z综合征失败变为成功，没有退步。精确McNemar双侧p=0.00048828125；FER绝对降低1.2个百分点，相对降低25%，平均CPU增加约30.2%。CPU为原批次与补跑批次累计的进程CPU时间，非统一时刻的延迟基准；中断和修复调试开销未计入算法CPU。

第一矩阵同信道点独立1000帧为41→35次逻辑失败，第二矩阵为48→36。两个矩阵均出现正向增益，增强了跨矩阵重复证据，但尚未覆盖多信道点、第三矩阵、同计算预算强基线或尾延迟；不能声称解决全部Z错误或已满足发表标准。两个实验的哈希、评估适用分支和统计边界必须保留，不能将这次评估器修复计作译码算法改进。

结果目录：`results/matrix20261019_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20353000/`，包括summary.csv、pairwise.csv、failure_transitions.csv及补跑记录。
