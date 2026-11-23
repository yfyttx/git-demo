# 研究任务与进度

更新时间：2026-11-20

本文件记录 GF(256) 级联量子 LDPC 项目的完成事项、下一步任务和当前问题。判断译码改进以逻辑 FER 为准；综合征通过、局部投影修复或小样本现象不能单独作为成功证据。

近期工作汇总：[近期尝试、结果与原因分析](experiments/gf256_cascade/RECENT_WORK_REVIEW_20261116.md)。
Z 侧机制审计：[精确小例与实际困难帧](experiments/gf256_cascade/Z_FAILURE_MECHANISM_20261116.md)。

## 当前状态

- **9 位内块的同外码 5000 帧对照已复核**。在 `alpha=100, p=0.120` 使用新噪声种子 `2026112300–2026112349`。每套外码 50 批，同一帧分别运行独立辅助位基线和级联版。矩阵 1 的逻辑失败为 `53→37`，矩阵 2 为 `49→34`。级联版各独自救回 30、31 帧，基线各独自救回 14、16 帧。级联减基线的 FER 差为 −0.32 和 −0.30 个百分点，配对 bootstrap 95% 区间分别为 [−0.58, −0.06] 和 [−0.58, −0.04] 个百分点。精确 McNemar 双侧检验原始 `p=0.0226/0.0400`，两项 Holm 校正后均为 `p=0.0453`。级联版消除了本批 X 综合征失败，同时 Z 综合征失败由 27→37、19→34。整批 CPU 时间每帧约低 22%。主张范围为理想综合征下这一信道点和同外码、同码率、相同局部 Hadamard 的内部基线。[配对结果](experiments/gf256_cascade/results/inner_code_comparison_alpha100_p0.12_confirm5000/summary.json)

- **同外码、同码率的 9 位内块先导对照已完成**。固定 `alpha=100, p=0.120`、理想综合征、Hadamard 位置 1、3、5、7、9 和每侧 100 轮 BP。两套矩阵各 300 个配对噪声帧。独立辅助位基线与级联方案的逻辑失败分别为矩阵 1 的 `5→5` 和矩阵 2 的 `2→1`，配对双侧 McNemar 检验均为 `p=1`。级联版 X 综合征失败 `3→0`、`2→0`，Z 综合征失败 `2→5`、`0→1`。本次整批 CPU 时间每帧约低 22% 和 19%，所有计时包含仿真和评分。当前样本尚不能确认逻辑 FER 改善。主指标为逻辑 FER，独立辅助位上的稳定子等价残差会记入物理 WER。[配对数据](experiments/gf256_cascade/results/inner_code_comparison_alpha100_p0.12_pilot300/summary.json) · [运行脚本](experiments/gf256_cascade/run_inner_code_comparison.py)

- **同数量局部 Hadamard 对照已完成**。固定 alpha=100、两套 5238 位矩阵、100 轮 BP。第 1、3、5、7 位固定选中。p=.1000 时第五个选第 8 或第 9 位，两套各 400 帧均无逻辑或物理残差。p=.1300 的独立 200 帧，第 2/8/9 位的逻辑失败分别为矩阵 1 的 34/30/30 和矩阵 2 的 33/25/32。未发现第 9 位的稳定优势。此结果修正原先四位加至五位的解释。按同数量掩码，块内分阶段与投影的输入熵在位置置换下相同。[配对报告](experiments/gf256_cascade/PRA_EQUAL_COUNT_HADAMARD_20261119.md)

- **级联码局部 Hadamard 对照和码长复验完成**。固定每内块第 1、3、5、7、9 位共轭、相同物理信道。Fuentes 图 10 的同点 p=.0700, alpha=100 下，两套 5238 位矩阵的 1000 帧同噪声逻辑失败从 356 和 346 降为 0，另各 1000 帧独立验证也是 0 物理错误帧。两套在 p=.0800 均 0/1000。第 9 位消融在 p=.1000 使逻辑失败 161/400 和 152/400 降为 0。两套 10422 位同码率矩阵在 p=.0800/.1000 各 200 帧均零失败。比较采用理想综合征，物理稳定子经共轭成为混合 Pauli 校验。因此数据支持偏置定制编码的贡献和同点 WER 基准，电路噪声与同参数外部基线仍是更广泛主张的待验证项。[报告及复现](experiments/gf256_cascade/PRA_ALPHA100_LOCAL_H_20261119.md)

- **强 Z 偏置下的 GF(256) OSD-1 配对完成**。交换方向、内层奇偶先验、`alpha=100, p=.0700`。两套矩阵各 300 帧同噪声独立配对。Z 综合征失败从 `106/300` 和 `111/300` 降为 0，逻辑失败仍为 `106/300` 和 `111/300`，对应全部成为漏检逻辑失败。各 20 个已知失败帧的 129 个候选列表均未覆盖正确 Z 逻辑类。OSD 中位开销约 8.5 ms/触发帧。保留实验开关，默认关闭。详见 [OSD 对照与覆盖审计](experiments/gf256_cascade/PRA_ALPHA100_Z_OSD_20261119.md)。

- **完成 Fuentes 图 10 的同点偏置仿真**。交换 CSS 方向并使用内块奇偶先验，`alpha=100, p=.0700` 两套矩阵各 1000 帧，逻辑失败和物理 WER 分别为 `356/1000` 和 `346/1000`，物理 QBER 为 1.790% 和 1.728%。论文图 10 同点 WER 约 3.46% 和 0.405%。另补 `p=.0600` 各 200 帧零失败及 `p=.0800` 的 198/200 和 197/200。在这些点 X 综合征全部通过，剩余瓶颈在 Z。详见 [同点对照报告](experiments/gf256_cascade/PRA_ALPHA100_MATCHED_POINTS_20261119.md)。

- **新版 X 侧 GF(256) 奇偶先验回退已验证**。同一 9 位内码、退极化 `p_D=.07125`、`f_m=.0475`，一级 BP 失败后，以每块测得的奇偶值更新 256 个首部模式的概率，再运行最多 100 轮 GF(256) BP。通过完整 X 综合征的输出才能进入 Z 条件先验。穷举 512 个 9 位模式的最大概率差为 `1.06e-22`。矩阵 1 的旧样本回放中，X 综合征失败 `68→0`，Z 综合征失败 `68→1`，逻辑失败 `69→1/1000`；此前二进制回退对同一批 68 帧救回为 0。独立 300 帧配对中，矩阵 1 逻辑失败 `16→0`，矩阵 2 为 `17→0`，均无退步帧。共触发 101 次，全部逐位恢复真实 X。每套新样本的零失败单侧 95% 上限约为 0.994%。此结果针对单个退极化点。复现入口为 `experiments/gf256_cascade/run_companion_staged_xz.py --p-depolar .07125 --z-prior conditioned --x-fallback parity_gf256`，汇总 `experiments/gf256_cascade/results/companion_x_parity_gf256_20261119/paired_summary.json`。

- **偏置信道下一轮已排定**。先汇总 `p=.0705, alpha=100` 下已启动的同帧对照，两套矩阵分别比较原方向、交换方向且仅按奇偶校验恢复尾位、交换方向且用奇偶校验更新 Gamma 先验。每套前 200 帧用于三组配对，交换方向加先验另有 800 帧用于固定方案的精度估计。已在 Fuentes 图 10 的原始横坐标 `p=.0700` 启动交换方向加先验的首批各 100 帧，预定各 1000 帧。使用 `plot_pra_alpha100.py` 从完成的 JSON 批次重建汇总、McNemar 配对表和与论文的物理指标对比图。若两套仍主要失于 Z 综合征，则在开发帧上比较 Z 侧候选回退，验证时用新种子和固定预算。论文点位比较注明码长和码率差异。

- **Fuentes 偏置信道实验已完成原始方向基线**。在 `p=.0705, alpha=100`（Pauli 概率比 X:Y:Z=1:1:100），两套 5238 位矩阵各 1000 帧，X 综合征失败均为 0，Z 综合征失败和逻辑失败均为 1000。Z 代表按似然恢复后，物理 QBER 为 7.013% / 7.015%，物理 WER 都是 100%。200 帧退极化回放在同一外层输出上由 WER 200/200 降至 11/200，逻辑失败仍为 11/200。原构造有 1552 个独立的 Z 错误综合征比特，偏置点的物理 QBER 和 WER 下界分别约为 0.805% 和 99.99975%，适用于逐位恢复。交换 CSS 的 X/Z 方向后，Z 错误有 2134 个独立综合征比特。现已启动交换方向的同噪声对照，分别使用及省略内层奇偶校验对先验的更新。详见 [物理恢复界限](experiments/gf256_cascade/pra_physical_bounds.md) 和 [批次计划](experiments/gf256_cascade/results/pra_alpha100_20261118/plan.csv)。

- **与 Fuentes 图 10 对照前的代表选择核查**。9 位内码的 Z 译码当前返回 `(u,0)`。同一投影的另一个代表是其逐位取反，两者相差一个内层 Z 稳定子。在独立退极化噪声 `p_D=.0705` 下，582 个尾位至少有一个真实 Z 分量的概率为 `1-(1-.047)^582`，所以固定零尾位会使物理 WER 几乎为 1。下一轮在固定矩阵、相同噪声与相同外层输出上比较固定代表和按信道似然选代表，核对逐帧综合征和逻辑判定一致，记录物理 QBER/WER。论文对照采用 `pX=pY=p/102, pZ=100p/102`，同步更新两侧先验。此项为代码分析，尚无新偏置信道实测值。

- **按 PRA 103, 022617 的指标完成固定矩阵六点扫描并扩充转折区样本**：逻辑 FER 随 `p_D` 从 `.0600` 的 `0/200` 上升到 `.07125` 的 `69/1000`；转折区 `.0675/.0685/.0705` 分别为 `9/1000`、`85/4200`、`116/2700`。同时保留物理 residual QBER/WER。中高端转折点已得到约85和116个逻辑失败事件，95% Wilson区间约为`1.64–2.50%`和`3.59–5.13%`；低端 `.0600/.0650` 仍只有0–1个事件，尚不能定阈值或支持投稿结论。报告与图在 `experiments/gf256_cascade/PRA_103_022617_METRIC_BENCHMARK_20261118.md` 和 `results/pra_metric_reference_sweep_20261118/`；新增源数据为`companion_pra_sweep_pD0p0685_seed20261200.json`、`companion_pra_sweep_pD0p0705_seed20261201.json`。
- **第二组外层矩阵的四点复核**：同一`[9,8,2]`内码，码长5238、逻辑量子比特1552、`R_Q=8/27`，CSS正交及秩验证通过。矩阵2在`p_D=.0650/.0685/.0705/.07125`各跑1000帧，逻辑失败分别`3/11/41/60`，趋势与矩阵1相近；每点仍只有3–60次事件，误差区间较宽。共同点上区间重叠，不能判定矩阵优劣。详见`experiments/gf256_cascade/PRA_103_022617_METRIC_BENCHMARK_20261118.md`；矩阵和逐点结果在`results/companion_binary_n9_k8_outer_path2/`及`results/companion_pra_matrix2_*.json`。
- **新9位构造 X 侧完整矩阵 BP 回退已完成独立1000帧配对（2026-11-18）**：在分阶段译码 X 综合征失败后，运行完整二进制 HZcas sum-product BP，最多100轮。固定矩阵、pD=.07125、seed20261125/20261126两批共1000帧，X综合征失败68次，完整矩阵回退68次均未通过；两组完整逻辑失败均69/1000=6.9%，配对救回/退步均0。之前已知困难帧4次回放也0次救回。说明该完整矩阵硬判决BP回退没有解决当前失败；不证明不存在其他候选或最优译码也失败。详见[X侧回退审计](experiments/gf256_cascade/X_FULL_HZ_FALLBACK_20261118.md)，汇总`experiments/gf256_cascade/results/companion_xfallback_independent1000_20261118.json`。

- **新9位伴随块构造：实际X估计条件先验配对完成（2026-11-18）**：矩阵不变，pD=.07125、fm=.0475、100轮、seed20261122、固定100帧；完整逻辑FER从固定先验70/100降为条件先验4/100，配对救回66、退步0。Z综合征失败70→4；两组X综合征失败均4，漏检逻辑失败均0。96个X综合征通过帧使用估计生成条件先验，全部逻辑成功；其余4帧按门控沿用固定先验、两侧均失败。512×512局部穷举核对通过。结果仅限该矩阵/信道先导，待独立扩大样本；可随后独立评估X完整HZ回退。详见[条件先验实验](experiments/gf256_cascade/COMPANION_CONDITIONED_PRIOR_20261118.md)，结果`results/companion_conditioned_pair100_20261118/`（相对于实验目录）。

- **同输入 BP 与 G2 信道对照完成（2026-11-18）**：C++/Python 共100例（含Gamma、多项式基Delta、对偶基Delta）最终硬输出、停止轮数、综合征状态差异均0。固定同一C++内核、Delta矩阵、fm=.0475、100轮、100个共享底层噪声帧，原始错误块失败5/100，G2=[I8|1]变换后使用精确先验为80/100，故意使用原始失配先验为100/100；均为综合征失败。这是外层块错误率，不是完整量子逻辑FER。性能下降在同一内核下仍存在，尚不能据此证明最优译码不可行。详见[BP核对报告](experiments/gf256_cascade/BP_IMPLEMENTATION_AUDIT_20261118.md)。后续待分离迭代预算与实际X条件信息的作用。

- **9位内码完整两侧简单译码已接入（2026-11-17）**：`run_companion_staged_xz.py`，X上半部Gamma BP+H2补尾位，Z经G2符号先验交给Delta BP；不含回退和X条件反馈。每点20帧、100轮、seed20261120：pD=.06完整逻辑0/20失败，pD=.07125为17/20，均为Z综合征失败；两点X均20/20精确恢复。结果`experiments/gf256_cascade/results/companion_staged_xz_20261117/`。已核对两侧域/二进制综合征映射、CSS正交，并按真实残差的稳定子成员关系评估完整逻辑FER。小样本尚不能宣称性能稳定。


- **直接 Kasai 双分量实验已完成（2026-11-17）**：独立 `kasai_paper` 入口，GF(256)、J=2、L=6、P=97、n=4656、RQ=1/3，固定一套 Gamma/Delta 矩阵，每个 f_m 点1000帧、50轮BP。f_m=.04/.045/.0475 的 Gamma 块错误为0/8/63，Delta为1/10/66。三点矩阵哈希一致。无 G2/H2、无级联处理；分量块错误口径见 `KASAI_PAPER_REPRO.md`，结果在 `experiments/gf256_cascade/results/kasai_direct_20261117_full/`。


- **单位内码退化对照与Kasai直接复现已完成先导（2026-11-17）**：G2=I8，H2为空校验；HX/HZ均1552×4656，RQ=1/3。100个相同噪声帧比较二进制矩阵综合征入口与直接外层入口，矩阵、综合征、译码输出差异均0；两分量各6/100块错误。另Kasai专用入口每点200帧，fm=.04时Gamma/Delta为0/1，fm=.0475时18/10；直接生成的Gamma/Delta与当前外层系数完全一致。原论文具体矩阵种子和迭代上限未完全确认，不能声称复现了原图数值。详见 `experiments/gf256_cascade/IDENTITY_KASAI_COMPARISON_20261117.md`。


- **新伴随块构造的简单 Z 侧 BP 已实现并试跑（2026-11-17）**：`run_companion_z_bp.py` 采用9个物理位经G2映射成一个GF(256)符号的精确信道先验，读取真实测量模型生成的综合征，外层 HZ_TIT 上最多100轮概率域/FWHT BP。不使用X侧条件信息、不加回退、不向译码器传真实错误。20帧、seed20261117：pD=0.07125时 Z逻辑失败17/20，均为Z综合征失败，耗时约23.5秒；低噪声pD=0.03时0/20。每帧核对二进制HX综合征与外层GF(256)综合征的一致性，残差用二进制HZ行空间判断稳定子等价。这是Z分量先导实验，不是完整X/Z逻辑FER，也不可与旧构造的条件译码结果直接比较。


- **伴随块乘法构造已按用户要求更新（2026-11-17）**：外层 GF(256) 元素展开为 8×8 二进制乘法块；内码 H2 是 1×9 全一行，G2=[I8|1] 为 8×9，selector=[I8|0]。此前生成的 [12,8,3] 内码不再是当前方案。新构造脚本 `experiments/gf256_cascade/construct_companion_cascade.py` 输出到 `results/companion_binary_n9_k8_matrix20261019/`（相对于实验目录）。实算 HX=1552×5238、HZ=2134×5238，秩分别1552、2134；物理码长5238、逻辑比特1552、码率8/27。G2 H2^T 和 HX HZ^T 均为零。新格式为二进制校验矩阵；另设简单Z侧BP入口（见上条），旧完整仿真入口尚未适配，不能沿用旧构造的性能结果。Fan论文式(12)位于PDF第3页/期刊3973页右栏，式(13)位于PDF第4页/期刊3974页左栏。

- **旧构造的简单外层BP对照已运行，不能评价新构造**：目录`experiments/gf256_cascade/results/g2_then_outer_bp_pilot20_20261117/`和`g2_then_outer_bp_iters200_pilot20_20261117/`。这些运行沿用原12列GF(256) path矩阵，仅切换译码流程，发生在上述构造要求澄清前。后续不得将它们表述成老师建议的块乘法构造实验。

- **全相邻窗口审计（2026-11-17）**：前三个持续失败帧每帧检查全部10个相邻对、共970个局部物理候选。相对 initial 阶段正确外层陪集覆盖总数12/33→13/33；新增的是20355020/2的坐标5，其他阶段原已覆盖该坐标，因此没有补上跨阶段仍缺失的关键坐标。三帧在线标签均与冻结实验一致，未产生新FER收益。结果：`experiments/gf256_cascade/results/pair_all_pilot3_20261117/`。用户随后要求优先验证简单的G2变换后分别运行TIT外层BP流程，当前按这一方向做对照。

- **相邻投影联合物理候选（2026-11-17）**：已加入默认关闭的诊断分支，取一对外层硬判决残差最大的相邻投影，直接在共享三列物理符号的局部 `H_X` 上做 GF(256) OSD-1，同一次求解同时满足两组外层综合征。冻结版前三个持续Z失败帧各生成97个物理候选；三帧选中的6个投影原正确陪集覆盖0/6，加入联合候选后仍0/6。在线输出3/3与原实验一致，165条原候选审计记录完全相同。说明此有限预算及后验中心未补上缺失候选，不是新FER结果。详见[联合物理候选先导](experiments/gf256_cascade/PAIR_WINDOW_CANDIDATES_20261117.md)。

- **相邻投影候选生成先导（2026-11-17）**：已实现默认关闭的 `GF256_CASCADE_NEIGHBOR_CANDIDATE_AUDIT`，按译码后验选邻居列表前4项，经path共享物理符号链生成新先验，再做外层BP和OSD，两轮累计。仅离线诊断，未改变在线判决；固定排序的前3个Z失败帧已完成：跨阶段正确陪集覆盖分别7/11→7/11、7/11→7/11、6/11→6/11，无新增覆盖；3帧在线结果和165条原候选审计记录均保持一致。此版本不接入在线译码。列表反馈可能重复计数，属于未见收益的启发式，不能断言更强联合方法无效。结果目录 `experiments/gf256_cascade/results/neighbor_candidates_pilot3_20261117/`。

- **已完成：冻结版全新1000帧三模式配对**：第二套 path 矩阵、内长12、p_D=0.07125、seed 20355000–20355199。tuned 51/1000（5.1%）、SF 46/1000（4.6%）、SF+layered 35/1000（3.5%）；组合模式仍有34次Z综合征失败和1次漏检逻辑失败。CPU分别9.25、12.35、13.07秒/帧。来源：`experiments/gf256_cascade/results/paired_three_modes_fresh1000_matrix20261019/summary.csv`。
- **Z失败候选覆盖审计（2026-11-17）**：34个持续Z失败帧已全部回放，X均精确恢复。合并初始、重试及layered各阶段的离线候选后，每帧仍有2–6个投影缺少正确外层陪集；33/34帧的缺失坐标形成连续一段。这支持优先改进共享内层符号参与的联合候选生成，不能仅靠重排当前列表。结果仅适用于当前有限候选预算，未产生新的FER增益。34个失败帧及3个对照共37项已全部完成，回放结果标签37/37与冻结实验一致。诊断报告：[候选覆盖审计](experiments/gf256_cascade/Z_FAILURE_COVERAGE_20261117.md)。目录：`experiments/gf256_cascade/results/fresh34_logical_coverage_audit_20261117/`。

- **研究目标**：定位并改善级联结构 Z 侧译码失败，取得可重复的逻辑 FER 改进。
- **仿真基线**：GF(256)、内长 12 的 path 级联矩阵，矩阵种子 `20261018`；常用点 `p_D=0.07125`（`f_m=0.0475`），诊断点 `p_D=0.072`（`f_m=0.048`）。
- **基线译码器**：`conditioned_robust_tuned`。当前仍是对照基准；Z 侧存在 syndrome failure 和少量 syndrome-pass 的漏检逻辑失败。
- **最近结论**：fixed tuned 相对 raw robust 在两个矩阵的各300帧比较中有明确改善；syndrome-flip 相对 tuned 在第一矩阵独立1000帧上为41→35次逻辑失败，CPU约增加25%。第二矩阵1000帧也有增益（48→36），尚缺多信道点及同预算强基线确认。joint_projection_safe 的100帧负结果仍成立：只把10次综合征失败转为10次漏检逻辑失败。
- **最近完成**：PEG 多行 `H2` 的 3×12 候选及另一组系统形 `G2` 均已构造。现在更正此前的记录：`G2` 的某一零空间基不与 `[I_9|0]` 逐项一致，并不自动意味着 CSS 构造错误；同一 `H2` 的满秩零空间基仅换基时，内码和级联校验行空间不变，但 BP 图表示和符号坐标解释可能变化。系统形 `G2=[I_9|A]` 是便于按指定信息集解释/译码的表示约定，不是单凭它就能修复译码的证明。系统形候选通用 Turbo 的 3 帧机制检查仍为 `3/3` Z syndrome failure；精确投影的一帧 smoke 满足 syndrome 但逻辑失败。均非 FER 证据。
- **最近验证**：通用 `G2` Turbo 的 GF(8) 树形 toy 对照覆盖 3 组 X 条件先验 × 8 个 syndrome，共 24 例；每个物理符号后验边缘与独立穷举一致，最大绝对误差 `7.31e-13`，由 `x_estimate` 生成先验与显式先验的最大差为 `0`。测试退出条件现同时检查后验和先验映射误差。这验证局部因子和 syndrome 切片在无环小例上的消息计算，不证明 GF(256) loopy BP 的收敛或逻辑 FER。
- **修正评估器后的 20 帧同噪声码比较**：在 `p_D=0.04`、seed `20330400–20330419` 上，path `0/20`、PEG `18/20` 逻辑失败；配对为 path-only 成功 18 帧、无反向救回，精确 McNemar 双侧 `p=7.63e-6`。两边平均物理错误均为 `2237.15`，并用 seed `20330400` 的导出真值逐位确认 X/Z 噪声向量完全相同。PEG 的 18 次失败全是 Z syndrome failure，漏检逻辑失败 0。当前这只比较两个码率不同的单一矩阵（path `11/36`、PEG `1/4`），因此说明该 PEG 候选在此点很差，不足以断言 PEG 方法普遍差。HX 符号 Tanner 图 4 环数为 path `35,502`、PEG `61,110`，提示 PEG 候选图更易受短环影响，但尚不能单独证明因果。
- **PEG 上的后处理负结果**：相同 `p_D=0.04`、20 帧配对中，`conditioned_robust_tuned` 与 tuned + syndrome-flip 都是 `18/20`、相同 18 帧 Z syndrome failure，救回/退步 `0/0`，McNemar `p=1`；CPU 从 `12.77` 增至 `39.92 s/frame`。该 syndrome-flip 设置没有解决这个矩阵的失败模式。
- **path 矩阵独立 100 帧译码器复验**：在固定 GF(256) path `H2=1×12` 矩阵、`p_D=0.07125`、(f_m=0.0475)、新种子 `20350000–20350019` 上，tuned 与 adaptive 均 `6/100`（Wilson 95% `[0.0278,0.1248]`），失败集完全相同且都是 Z syndrome failure；adaptive CPU `11.55` 对 tuned `9.49 s/frame`。此前 50 帧的 adaptive 改善未在此独立样本复现。
- **独立样本中的失败帧机制回放**：对 tuned 的 6 个 Z syndrome 失败帧逐帧回放，long 模式仍有 `6/6` 逻辑失败，其中 2 帧 syndrome 通过但逻辑失败；syndrome-flip 模式为 `4/6` 逻辑失败（2 帧救回、0 退步），剩余 4 帧仍为 Z syndrome failure。按代码分支，syndrome-flip 在基线 syndrome-valid 帧不运行且与 tuned 共用初始译码，因此可在同一 100 帧样本上精确组合为 tuned `6/100` 对 tuned+syndrome-flip `4/100`；配对救回 2、退步 0，McNemar `p=0.5`，不能声称显著改善。6 个预选困难帧上的回放 CPU 合计分别为 `171.7 s`（long）和 `465.6 s`（syndrome-flip），属于机制诊断，非无条件平均复杂度。
- **独立 300 帧配对复验已完成**：同一 path 矩阵、`p_D=0.07125`、种子 `20351000–20351059`，tuned `9/300=3.0%`，tuned+syndrome-flip `7/300=2.33%`；Wilson 95% 区间分别 `[0.0159,0.0560]`、`[0.0113,0.0474]`。配对救回 2、退步 0，McNemar `p=0.5`；Z syndrome failure `9→7`，漏检逻辑失败均为 0。CPU `8.90→10.72 s/frame`，增加约 20.5%。方向与前一组 100 帧一致，但合并这两个互不重叠样本也只有 4 次救回、0 次退步，精确 McNemar `p=0.125`，仍不能声称显著收益。完整结果目录 `experiments/gf256_cascade/results/path_tuned_vs_syndrome_flip_holdout300_pD0p07125_seed20351000/`。
- **固定设置的 1000 帧独立确认已完成**：同一 GF(256) path 矩阵、`p_D=0.07125`（`f_m=0.0475`）、种子 `20352000–20352199`，tuned `41/1000=4.1%`（Wilson 95% `[3.04%,5.51%]`），tuned+syndrome-flip `35/1000=3.5%`（`[2.53%,4.83%]`）。两边均有 3 次漏检逻辑失败；Z syndrome failures `38→32`。配对救回 6、退步 0，精确 McNemar 双侧 `p=0.03125`。CPU `13.37→16.73 s/frame`（约 +25.2%）。这是单一矩阵、单一信道点上有统计信号的复验，不足以单独证明跨矩阵/信道的稳健提升。目录 `experiments/gf256_cascade/results/path_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20352000/`。
- **第二矩阵1000帧已补完（2026-11-16）**：tuned 48/1000（4.8%），syndrome-flip 36/1000（3.6%）；Z syndrome failures 47→35，漏检逻辑失败均为1；救回12、退步0，McNemar p=0.00048828125，CPU13.26→17.26秒/帧。中断原因已复现并修复：显式path G2误入系统形X稳定子检查。单校验分支恢复按H2缩放判定；实际2134条HX生成行及零向量全部通过，既有Z侧核验也通过。原965帧均未进入错误分支，可复用，缺35帧/方法已重跑；completion_provenance.json保留两版程序哈希和复用依据。输出目录 `experiments/gf256_cascade/results/matrix20261019_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20353000/`。
- **Z 侧失败机制审计（2026-11-16）**：第二矩阵独立1000帧中，两方法共同的35个 Z syndrome failure 各自运行8次 syndrome-flip 扰动，280次尝试均无 syndrome-valid 候选；12个救回帧各有有效候选；另1帧 syndrome-pass 但落入错误逻辑类。全部48个 tuned 逻辑失败帧 X 分量均精确恢复。选取的3个持续失败帧各有2–6个外层投影与真值不同且对应外层 syndrome 不通过；额外尝试每帧16个单比特及36个最多双比特扰动，仍均无有效候选。GF(8) 六符号 toy 的精确逻辑 MAP 对照显示信道歧义和现有 Turbo 内核算法损失并存，但 toy 数值不能外推到 GF(256) 长码。详见上述机制审计文档。
- **path 外层候选覆盖（2026-11-16）**：对最早3个持续 Z 失败帧，外层 BP 后验与共享符号融合后验各做129候选/投影的 OSD-1 审计，正确投影覆盖均为9/11、9/11、5/11；未覆盖投影的最近候选仍与真值相差38–141个 GF(256) 符号。现有列表重排不足以修复这些帧。已有 layered 调度在相同三帧中分别产生1次漏检逻辑失败、1次逻辑成功、1次 Z 综合征失败，不能据此估 FER。单投影信息量诊断：外层码率至少2/3；在 `p_D=0.07125`、正确 X 侧信息与独立单投影模型下，条件容量约0.652745，低于2/3。这不是整体级联码容量界，反而强调跨投影联合信息的重要性。详见机制审计文档及其可复算脚本。
- **分层回退配对结果（2026-11-16）**：新增 `conditioned_robust_tuned_syndrome_flip_layered` 模式；只在 tuned、syndrome-flip 均未通过完整 Z 综合征时尝试分层 BP。第二矩阵原1000帧离线合并为 tuned 48、SF 36、SF+layered 27 次逻辑失败；第一矩阵另一批1000帧为 SF 35→SF+layered 30。分层分别在两批中新增9/5次逻辑救回，也分别新增1次漏检逻辑失败。固定帧直接运行新模式与离线合并一致。第二矩阵是开发数据，第一矩阵历史结果已被观察；这支持进一步验证，不是论文级独立确认。初步 CPU 开销和新噪声100帧已测；更大独立样本与新信道点仍待测。
- **新种子100帧先导（2026-11-16）**：第二矩阵 `p_D=0.07125`、种子 `20354000–20354019`，SF `9/100` 逻辑失败，SF+layered `5/100`；配对救回4、退步0，McNemar双侧 `p=0.125`；Z综合征失败8→4，漏检逻辑失败均1；CPU15.04→15.98秒/帧（约+6.2%）。这是方向正向的新样本，但事件少，尚不能称为统计确认。完整数据在 `experiments/gf256_cascade/results/paired_sf_layered_fresh100_matrix20261019_seed20354000/`。
- **逻辑等价评估已修复**：此前稳定子检查只支持单行 `H2` 或一个固定双组特例；多行 PEG `H2=3×12` 被错误地按“一个尾部坐标”处理。现对系统形 `H2=[A|I_3]`、`G2=[I_9|B]` 使用完整的内码投影检查，Z 陪集键也先消去三个 parity 坐标，再按外层行空间约化；X 稳定子检查核验全部 `G2` 坐标。actual PEG 矩阵上 3492 条 `H_Z` 生成行全部通过，加入每行后陪集键不变，另一个非平凡逻辑见证被正确区分。旧 PEG 结果中 syndrome failure 计数仍有效，但 syndrome-pass 逻辑分类和使用旧陪集键的结果需重跑。
- **修正评估器后的同帧矩阵/译码器比较**：在 `p_D=0.04`（`f_m=0.02667`）、8 个相同噪声帧（seed `20330200–20330207`）上，旧 path 矩阵 `0/8` 逻辑失败；PEG 3×12 矩阵 `3/8`，三次都是 Z syndrome failure、undetected logical `0`。两侧平均物理错误完全相同（`2221.125`）；精确配对差异是 path-only 胜出 3 帧、无反向救回，McNemar `p=0.25`，样本仍很小且两矩阵码率不同。此前 PEG 的 `6/8` 是旧评估器把 3 个稳定子等价恢复误判为逻辑失败，已弃用。相同 PEG 矩阵与帧上，tuned 和 layered 均 `3/8`，失败集合相同；CPU `9.42` 与 `9.77 s/frame`，schedule 没有可见 FER 增益。结果分别见 `experiments/gf256_cascade/results/path_vs_peg_same_frames_pD0p04_seed20330200_corrected_eval/` 和 `experiments/gf256_cascade/results/peg3_layered_vs_tuned_corrected_eval_pilot8_pD0p04_seed20330200/`。
- **构造适用范围**：当前 `construct.cpp` 使用同一个 GF(256) 上的普通 Kronecker 积，将 TIT 外层矩阵与内层 `G2/H2` 组合。它满足已核验的 CSS 正交关系，但不能直接称为 Fan 等人 AQCTPC 定理构造的逐项复现：参考构造中外层码和内层码分别位于 GF(q) 与 GF(q^k) 等不同域，并要求内校验矩阵特定右块可逆以支持分阶段恢复。现有性能结论只适用于本项目这个同域构造；若要检验论文定理的预期，须另做域与坐标忠实的构造及其译码器，不能把换译码器当成等价替换。
- **OSD 机制负结果**：systematic PEG 候选上的单帧 OSD 后处理测试（`p_D=0.07125`、seed `20272000`）枚举 65,281 个候选并找到 syndrome-valid 候选，但该帧仍为逻辑失败，且属于漏检逻辑失败；残余符号数 2,879。它直接说明满足 syndrome 不足以判译码成功，也不提供 FER 估计。结果在 `experiments/gf256_cascade/results/h2_peg3_n12_systematic/osd_pilot1_seed20272000/`。
- **generic-G2 内因子迭代扫描**：systematic PEG 同一矩阵、同一噪声 seed `20272000` 的第 0 帧，固定 `p_D=0.07125`、generic damping `0.9`，将局部内因子 BP 迭代设为 1/5/15/30。未满足 Z 校验数依次为 `1549/912/910/910`，四种设置均 syndrome failure 和 logical failure。增加到 5 次有明显残差改善，但 15/30 进入平台；这是单帧机制证据，不能估 FER，也没有找到译码成功配置。逐次日志和汇总在 `experiments/gf256_cascade/results/peg3_generic_inner_iters_single_frame_seed20272000/`。
- **generic-G2 迭代跟踪与长重放**：新增默认关闭的 `GF256_CASCADE_GENERIC_TRACE`，记录每轮 outer→inner 与 inner→outer 外信息的平均/最大总变差和熵，不改变常规译码行为。相同困难帧上，20 轮后 outer 消息仍有明显变化；`conditioned_robust_long` 的 60 轮重试到 round 46–60 平均总变差约 `0.34–0.40`，最终仍有 995 个未满足 Z 校验，逻辑失败。开启跟踪后去掉诊断行，与关闭跟踪的同参数输出逐字相同；GF(8) 24 例精确测试也重新通过。详见同一结果目录的 `trace_iters5.log`、`long_turbo_trace.log`。
- **投影级故障拆分**：trace 的末轮 outer 投影未满足校验向量为 `[0,0,60,58,67,84,122,62,63]`；默认物理 APP 硬判决的矩阵导出 G2 投影中仅坐标 0、1 与真值一致且 syndrome 通过，坐标 2–8 均错且不通过。启用 exact projection recovery 后 9 个 projection 均 syndrome-valid，但坐标 2–8 仍与真值有 `343–454` 个符号差异，整帧 syndrome 通过仍逻辑失败。这分别暴露了内外硬判决一致性缺失和 syndrome-valid 投影的逻辑类选择问题。projection audit 已改为从实际 HX 提取 G2 并校验 `HX=H_outer⊗G2`，修复了旧 star/path 专用公式对一般 PEG G2 的误诊。完整结果在 `experiments/gf256_cascade/results/peg3_generic_inner_iters_single_frame_seed20272000/`。
- **OSD 投影候选覆盖审计**：仅用 decoder posterior/syndrome 构造候选，真实错误只在译码后做覆盖检查。OSD-1 的 129 个候选（32 free × 4 alternatives）和 1,025 个候选（128 × 8），以及 OSD-2 的 289 个候选（12 free × 2 alternatives），在初始/retry 阶段都只覆盖投影 0、1 的真实外层词；投影 2–8 均漏覆盖，最近候选仍距真实词约 205–391 个符号。说明对这些列表做类重排不能修复 coverage miss。仅为同一困难帧的离线机制证据，不是 FER 结果。逐坐标日志位于上述目录的 `outer_osd_candidate_coverage*.log`。
- **同码换基对照**：PEG 目录中的 `H2`、外层矩阵和 `HZ` 哈希一致，旧/系统形 `G2` 跨 3 个相同噪声种子配对运行；两边都是 `3/3` Z syndrome failure、`3/3` logical failure，均值物理错误数同为 `3955`。未满足校验数分别为旧基 `1016,1076,222`、系统形 `1259,1205,441`。小样本显示 BP 残差会随行基/图表示变化，但没有 FER 优劣证据；旧基不能因 syndrome residual 较小就判为更优。输入、矩阵哈希、TSV 和逐帧日志保存在 `experiments/gf256_cascade/results/peg_g2_basis_same3_seed20272000/`。
- **当前新路线**：受可靠度引导的 GF(256) syndrome-flip 在第三组独立 300 帧中为 `14/300→13/300`，救回 1、退步 0，McNemar `p=1.0`，CPU 慢约 36%。合并 700 帧探索样本为 tuned `36/700`、syndrome-flip `31/700`，配对救回 5、退步 0，精确 McNemar 双侧 `p=0.0625`；只能称为有待复验的正向趋势，尚非显著或论文级证据。下一阶段转向 H2/级联结构及投影译码的独立设计实验。
- **矩阵映射核查**：实际 GF(256) 的两套 path 矩阵通过了 HZ 行正交、稳定子成员、陪集键不变性检查；并找到非平凡逻辑见证验证逻辑类检查能区分它与稳定子。这排除了当前检查范围内的稳定子/陪集映射错误，但不构成译码性能证明。
- **当前研究焦点**：正确逻辑陪集没有进入候选列表（coverage miss），或候选代表被选到错误陪集（selection error）。PEG 试验进一步表明，强制满足完整综合征也可能把失败从可检测转成漏检逻辑失败；下一步要先核查投影变量映射、内码维数/码率及逻辑类选择，再尝试新算法。

## 已完成

- [x] 建立 GF(256) path 级联矩阵、可复现的仿真和同帧配对比较流程；保留多矩阵种子、信道点和译码器的历史结果。
- [x] 比较 `conditioned_robust_tuned`、SI、joint projection、safe joint projection、soft-X、layered outer BP、augmented-check、overcomplete 等译码原型；这些原型多数未展示可重复的逻辑 FER 改善；fixed tuned和syndrome-flip的独立证据及限制见当前状态，负结果均保留。
- [x] SI 局部精确求解通过 GF(8) 穷举核对；100 帧配对开发集 tuned 与 SI 均为 `13/100`，无配对救回，故不继续验证该版本。
- [x] 原始 joint projection 的 100 帧配对结果：tuned `10/100`、joint `26/100`；救回 1 帧、退步 17 帧，McNemar `p=0.000145`。降低 syndrome failure 的同时逻辑 FER 变差。
- [x] safe joint projection 的同帧复测：两者均 `10/100`；救回/退步 `0/0`，McNemar `p=1.0`。tuned 有 10 次 Z syndrome failure、0 次漏检；safe 有 0 次 syndrome failure、10 次漏检逻辑失败；CPU `11.6025` 对 `12.0405 s/frame`。
- [x] overcomplete 外层校验 5 帧 pilot：两者均 `0/5`，样本不足以比较 FER；新方法约慢 12 倍，不扩展当前实现。
- [x] 新增并运行 GF(8) 逻辑类配分函数恒等式精确核验：64 个 syndrome/逻辑类中，直接枚举与稳定子参数求和的最大质量差为 `1.388e-17`。
- [x] 用项目 C++ `FieldBasis`、条件 Pauli 先验和 3 个非零缩放核对内部 Z 符号先验与物理位先验的映射，最大误差为 `0`；既有后验检查误差不超过 `1.11e-16`。
- [x] 对两套实际 GF(256) path 矩阵逐行核验：每个 `H_Z` 行均与 `H_X` 正交、被 `isCascadeZStabilizer` 接受，并且将该行加到任意确定性测试向量后 `cascadeZCosetKey` 不变；两套矩阵均为 `0/2716` 异常。
- [x] 额外从外层 `H_X` 的零空间构造一个 `H_X` 综合征为零、但不属于 `H_Z` 行空间的逻辑见证；实际矩阵上 `isCascadeZStabilizer` 判为非稳定子，且 `cascadeZCosetKey` 非平凡。两套矩阵都通过。
- [x] 记录两篇近期方法作为候选思路： [Logical-Space Theory of OSD](https://arxiv.org/abs/2610.24503) 与 [Certified decoding of quantum LDPC codes](https://arxiv.org/abs/2608.25545)。它们是预印本；其算法和性能结论不能直接视为适用于本项目。
- [x] **定位新 9 位构造的 X 侧失败并做 GF(8) exact-MAP 对照**：独立1000帧中原有68个 X syndrome failures 全部表现为外层符号估计错误；中位每帧错112/582个外层符号，外层194行中位52行 syndrome residual。构造8比特 GF(8) 级联类比小码并穷举全部256个错误：BSC 0.0475 下 exact MAP block-error 2.96%，完整二进制 BP 5.22%，其中 BP syndrome failure mass 5.01%。该玩具是经典物理块错误指标，不能外推成量子逻辑 FER。报告与逐帧审计见 `experiments/gf256_cascade/X_FAILURE_MAP_DIAGNOSIS_20261118.md`。
- [x] **对照 PRA 103, 022617 的QBER/WER图并给当前译码器增加残余物理指标**：在与原先1000帧相同的两个500帧种子上重跑，新增 X/Z/Pauli-union residual QBER、physical WER，同时保留 syndrome failures 与 stabilizer-coset logical FER。汇总为 QBER 0.325%/4.849%/5.108%，physical WER 6.8%/100%/100%，logical FER 69/1000（95% Wilson `[5.49%,8.64%]`）。物理残余 WER 不等于逻辑 FER；当前单点还不是性能曲线。参考论文指标映射、图计划与缺口见 `experiments/gf256_cascade/PRA_103_022617_METRIC_BENCHMARK_20261118.md`。
- [x] **完成 PRA 风格的固定矩阵信道扫描**：矩阵、100轮、conditioned Z prior、无 X fallback 固定；`p_D=.0600/.0650/.0675/.0685/.0705/.07125` 的帧数分别 `200/200/1000/1200/1200/1000`，逻辑失败 `0/1/9/27/48/69`。已输出双面板 PDF/SVG/PNG、逐点 CSV 和绘图脚本。低端0–1次失败仍不足以估计低FER，转折区点目前为单矩阵数据，标记为先导曲线而非最终投稿证据。路径 `results/pra_metric_reference_sweep_20261118/`（相对实验目录），解释见 PRA 指标报告。

## 未完成任务

- [x] **核验实际 GF(256) 稳定子/陪集结构**：两套矩阵的所有 `H_Z` 生成行都保持 `H_X` syndrome，并且保持当前 `cascadeZCosetKey` 不变；`20261019` 矩阵使用保存的 `H2` 缩放，主实验矩阵按仿真代码中缺省单位缩放处理。
- [ ] 将上述核验扩展到随机稳定子线性组合，并继续核对物理二进制展开与 syndrome 映射；当前测试覆盖所有生成行，尚未单独抽查随机组合和所有展开位映射。
- [ ] 将 GF(8) 条件先验及逻辑陪集配分函数核验扩展到实际 GF(256) 映射，先通过小规模/逐项精确核对，再开发近似类概率算法。
- [ ] 拆解困难帧的候选列表覆盖率与类代表选择错误；候选生成和评分只能使用译码可用信息，真实错误仅用于离线审计。
- [ ] 对多行 `H2` 的 GF(256) loopy generic Turbo 做困难帧消息/残差审计；GF(8) 树形 exact 对照已通过，下一步比较固定矩阵上的分层 schedule、阻尼及低树宽消元/局部精确更新，并用 syndrome-valid 与逻辑陪集 FER 双指标筛选。
- [x] 将 generic-G2 Turbo 跟踪扩展为每轮投影 syndrome 残差及从实际 G2 计算的硬判决投影审计；已区分默认 APP 的内外一致性失败与精确投影后仍发生的错误逻辑类选择。
- [ ] 改善投影 2–8 的外层软后验或候选中心，并协调共享内层符号约束；先在多个独立失败帧上确认正确投影/逻辑类进入候选集，再做固定参数的配对开发 FER，成功后再做独立验证。
- [ ] 研究 coverage-aware 候选生成或逻辑类概率估计；先在可穷举小实例验证正确性和复杂度，再做固定矩阵、固定噪声种子的配对开发实验。
- [ ] 检查 Z 侧未满足校验子图、低可靠度物理比特簇，以及 X 估计误差对条件 Z 先验的传播。
- [ ] 只有当开发集出现可重复的逻辑 FER 收益，才冻结代码、矩阵、参数与哈希，并在独立噪声种子/矩阵上做配对验证；报告置信区间、两侧 syndrome failure、漏检逻辑失败和 CPU 时间。
- [x] **完成 restart-belief global 对照入口及 pilot**：新模式使用 tuned 主译码参数，只在 Z syndrome 不通过时进入 global restart-belief 回退；已加入配对脚本并完成编译和单帧检查。20 帧 pilot 与 tuned 均 `1/20`，救回/退步 `0/0`，McNemar `p=1.0`，CPU `10.06→14.13 s/frame`，不继续扩大该版本。
- [x] **实现 GF(256) syndrome-flip 原型并完成先导验证**：以 tuned 级联 Turbo 的物理位后验可靠度排序，对最不确定 4 位枚举重量 1、2 的 (t)，计算移位 syndrome (s+H_Xt)，重译后将 (t) 加回；固定使用条件信道对数似然选候选，并对原 syndrome 做最终检查。预选困难帧上找到有效候选；独立 100 帧配对为 `7/100→5/100`、救回 2/退步 0、McNemar `p=0.5`、平均 CPU `9.76→14.93 s/frame`。这是探索性正向信号，不足以证明改进。
- [x] **完成 syndrome-flip 第一组独立 300 帧复验**：固定矩阵、`p_D`、代码哈希及环境预算（4 个低可靠度物理位，max weight 2）；tuned `15/300`、syndrome-flip `13/300`，救回 2、退步 0，McNemar `p=0.5`，Z syndrome failures `15→13`，漏检逻辑失败均为 0，CPU `9.84→13.67 s/frame`。
- [x] **完成 syndrome-flip 第二组独立 300 帧复制**：固定矩阵、信道点与参数；tuned `14/300`、syndrome-flip `13/300`，救回 1、退步 0，McNemar `p=1.0`；Wilson 95% 区间分别 `[0.0280,0.0768]`、`[0.0255,0.0727]`；CPU `10.01→13.65 s/frame`。与前两组独立样本合并后仍未达 `p<0.05`，因此将其保留为低频回退原型，不再单独扩展帧数。
- [x] **构造 PEG 内码的系统形零空间基表示**：对 `H2=[A|I_3]` 增加 `G2=[I_9|A^T]` 直接构造，并核对 `H2G2^T=0`、GF/二进制 CSS 正交为 0、H2 girth=6。该表示让前 9 列成为明确的信息坐标；它是坐标/图表示选择，不能称作旧基有错或已修复 FER。新目录 `experiments/gf256_cascade/results/h2_peg3_n12_systematic/build/`，矩阵尺寸 `HX=1746×6984`、`HZ=3492×6984`，推得率 `1/4`。通用 Turbo 同参数前三帧均 Z syndrome failure；精确投影一帧 syndrome 通过但逻辑失败。日志 `results/h2_peg3_n12_systematic/pilot20/peg3_systematic.log`、`exact_recovery_smoke.tsv/log`。
- [x] **保留 PEG/基表示试验的负结果与适用边界**：原 PEG 生成基与系统形基来自同一个 H2 零空间；已有 20 帧和 3 帧译码数据分别是不同因子图/基表示上的小样本结果，不能解释成不同内码性能，也不能用来估计 PEG 相对 path 的 FER。精确投影分支 `GF256_CASCADE_GENERIC_EXACT_RECOVERY=1` 默认关闭；旧候选上 syndrome 通过但 `20/20` 漏检逻辑失败，说明该硬投影不能作为性能改进。
- [ ] **定位并修复级联译码的候选覆盖问题**：先用精确小码对照验证投影/内码因子消息，再把真实 GF(256) 困难帧作离线候选覆盖审计；只使用译码可用信息选在线候选，真实错误仅用于事后统计。之后在冻结矩阵、固定参数的独立配对种子上验证逻辑 FER、两侧 syndrome failure、漏检逻辑失败和运行时间。
- [x] **在独立 1000 帧上确认 syndrome-flip 救回信号**：参数固定，完整配对样本 tuned `41/1000`、syndrome-flip `35/1000`，救回 6、退步 0，McNemar `p=0.03125`，CPU 慢约 25%。
- [x] **分类第二矩阵1000帧失败并做小例精确 MAP 对照**：逐帧区分35个候选缺失的 Z syndrome failure、12个救回帧和1个漏检逻辑失败；导出6帧机制回放，检查真实外层投影；完成 GF(8) toy 全枚举与三帧扩预算定向诊断。原始数据、脚本和局限详见 `experiments/gf256_cascade/Z_FAILURE_MECHANISM_20261116.md`。
- [x] **审计 path 外层候选覆盖与单投影信息量**：新增默认关闭的 path 候选列表离线审计；固定3个失败帧与2个控制帧比较 outer/joint 先验下的正确投影覆盖；用独立脚本计算退极化信道下单条单位 path 投影的条件容量。未把诊断信息输入在线译码器。
- [x] **实现 SF 后分层回退并做两个历史样本的同帧诊断**：组合模式只用完整综合征选择候选，固定帧实现输出与离线合并一致；第二矩阵1000帧 SF 36→27，第一矩阵1000帧 SF 35→30。漏检逻辑失败分别1→2、3→4，须与逻辑 FER 一并报告。详见机制审计文档。
- [ ] **冻结新组合模式并做全新配对验证**：在未用于方法选择的噪声种子、至少两个矩阵和邻近多个 `p_D` 点上，比较 tuned、SF、SF+layered，报告逻辑 FER、Z syndrome failure、漏检逻辑失败与 CPU/帧；预先规定选择规则，不能按真值挑分支。
- [x] **完成组合模式的全新100帧先导**：在固定第二矩阵上直接配对运行 SF 与 SF+layered，9→5 次逻辑失败、4 次救回、0 次退步、`p=0.125`，CPU约+6.2%；继续扩大新样本的任务仍未完成。
- [ ] **检验 syndrome-flip 的跨矩阵泛化**：在不同 TIT 外层/矩阵种子上复现完整配对实验；同时保持内码、码率、译码设置、`p_D` 与帧数一致，并报告逻辑 FER、Z syndrome failure、漏检逻辑失败及 CPU。一个矩阵上的显著结果不能当作跨矩阵结论。
- [ ] **PEG 多行内码下诊断译码失效域**：layered schedule 的 20 帧配对 pilot 与 tuned 均 `20/20` Z syndrome failure；先在较低 `p_D` 找到两种译码均进入非饱和区的运行点，再比较投影错误与逻辑 FER。若 PEG 候选仍持续全失败，需把它作为不同码率/内码距离的构造问题处理，不能把 schedule 残差变化当作译码器改进。
- [x] **完成修复后 PEG/path 与 syndrome-flip 配对诊断**：同一 `p_D=0.04`、20 个相同噪声帧上，path `0/20`、PEG `18/20`；PEG tuned 与 tuned+syndrome-flip 均 `18/20`，救回/退步 `0/0`。这表明该 PEG 候选与当前译码器组合失配，但由于码率也不同，不能推广成 PEG 方法普遍较差。完整数据见研究记录。
- [ ] **优化多行 PEG 内码时先控制码率和图结构**：已完成一个 20 帧同噪声对照，单个 PEG3×12 候选明显落后 path，但内码维数/码率与短环数同时变化。下一轮需在同一码率/相同外层矩阵下生成多个 H2 候选，检查其真实 H2/G2、HX/HZ 短环和小重量逻辑算符，再选择候选做 FER；不要再仅按 H2 girth 选码。

## 当前问题与限制

- syndrome-valid 只表示校验方程满足，不保证恢复到正确逻辑陪集；safe projection 已明确展示了“检测失败减少、逻辑 FER 不变”的情况。
- 现有 OSD 类候选审计中，若干困难帧的正确逻辑陪集没有进入候选列表。只重排已有列表无法解决 coverage miss。
- GF(8) 精确测试只验证 toy 码的代数恒等式及先验坐标，**尚未证明实际 GF(256) 级联矩阵的陪集映射、类概率计算或 FER**。
- 当前同域 GF(256) Kronecker 构造与引用的异域 AQCTPC 定理条件不完全一致；不得将相关论文的距离/译码结论直接套用到当前码。
- 实际 GF(256) 检查验证了生成行层面的 CSS 正交、陪集键不变量，并抽查一个已知非平凡逻辑见证；尚未验证类概率近似、困难帧候选覆盖率或逻辑 FER 改善。
- 最近几种复杂译码原型尚未展示 FER 收益；在候选覆盖和实现映射没有验证前，不应扩展高成本大规模 FER 仿真。
- syndrome-flip 三组独立样本总计 tuned `36/700`、syndrome-flip `31/700`，救回 5、退步 0，精确 McNemar 双侧 `p=0.0625`；点估计改善仍未达到常用显著性门槛，且 CPU/frame 增加约 40%。其二进制参考论文也没有验证本项目的 GF(256) path 级联模型。
- 目前最强基线仍是 `conditioned_robust_tuned`；其性能是否跨矩阵种子和信道点稳定，仍需独立数据验证。

## 最近验证记录

### 新 9 位构造 X 侧失败诊断与 GF(8) MAP/BP 对照（2026-11-18）

- 报告：[`experiments/gf256_cascade/X_FAILURE_MAP_DIAGNOSIS_20261118.md`](experiments/gf256_cascade/X_FAILURE_MAP_DIAGNOSIS_20261118.md)
- 脚本：[`experiments/gf256_cascade/audit_companion_xfailures.py`](experiments/gf256_cascade/audit_companion_xfailures.py)；数据：`experiments/gf256_cascade/results/companion_xfailure_map_audit_20261118.json`。
- 在两组独立500帧原实验中重放68个 X syndrome failures；真实 X 错误仅用于事后审计。外层长度582，错误符号数中位数112（范围55–199），194个外层检查的残差数中位数52（范围27–92）。
- 8比特 GF(8) 级联类比小码对全部256个物理错误穷举；BSC `p=0.0475`，exact physical MAP block error `2.96%`，同一完整矩阵 binary BP `5.22%`，BP syndrome-failure mass `5.01%`。
- 解释边界：小码说明 BP 算法损失可以与信道歧义并存；它不是量子逻辑陪集译码，也不能把 FER 数值外推到实际码。下一步需围绕实际失败帧设计更强的 syndrome-valid 外层候选生成，并固定后用新种子验证。

### PRA 103, 022617 指标对照与残余物理错误统计（2026-11-18）

- 指标与图形参照：[`experiments/gf256_cascade/PRA_103_022617_METRIC_BENCHMARK_20261118.md`](experiments/gf256_cascade/PRA_103_022617_METRIC_BENCHMARK_20261118.md)。该文将 QBER 定义为译码后物理错误位比例，WER 为至少有一物理错误的块比例；本项目另将 stabilizer-equivalence 逻辑 FER 作为量子主指标。
- `run_companion_staged_xz.py` 现在逐帧记录 X/Z 残余权重及 Pauli 并集残余，汇总 QBER、physical WER、syndrome failures 和 logical FER。
- 在冻结的 `20261019` 新9位矩阵、depolarizing `p_D=0.07125`、`f_m=0.0475`、100轮、conditioned Z prior、无 fallback 下，对旧种子 20261125/20261126 各500帧重复输出新指标：logical FER `69/1000=6.9%`，Wilson 95% `[5.49%,8.64%]`；QBER X/Z/union `0.325%/4.849%/5.108%`；物理残余 WER X/Z/union `6.8%/100%/100%`。WER_Z=100% 不等于逻辑 FER=100%，因为 stabilizer-equivalent 残差可以无逻辑损伤。
- 新指标数据：`results/companion_pra_metrics_seed20261125.json`、`results/companion_pra_metrics_seed20261126.json`、聚合 `results/companion_pra_metrics_1000_20261118.json`。仅一个 channel point，不能画成可信的性能曲线；多点扫描与独立矩阵种子仍待完成。
- 多点扫描数据与图：`results/pra_metric_reference_sweep_20261118/source_data.csv`、`pra_metric_reference_pilot.pdf/.svg/.png`；生成脚本 `plot_pra_metric_reference_sweep.py`。曲线趋势已显现，但低FER点需要追加帧，且缺独立矩阵种子。

### GF(8) 配分函数与条件先验验证

- 脚本：[`experiments/gf256_cascade/test_coset_partition_identity.py`](experiments/gf256_cascade/test_coset_partition_identity.py)
- C++ 核验：[`experiments/gf256_cascade/test_soft_x_prior.cpp`](experiments/gf256_cascade/test_soft_x_prior.cpp)
- 输出：[`experiments/gf256_cascade/results/coset_partition_toy/validation.log`](experiments/gf256_cascade/results/coset_partition_toy/validation.log)
- 结果：GF(8) toy 码 64 个 syndrome/逻辑类全部匹配；最大质量误差 `1.388e-17`。条件先验缩放映射最大误差为 `0`。这是精确性单元验证，不是译码性能结果。

### 实际 GF(256) 稳定子/陪集键不变量

- 验证程序：[`experiments/gf256_cascade/test_cascade_coset_invariants.cpp`](experiments/gf256_cascade/test_cascade_coset_invariants.cpp)
- 输出：[`experiments/gf256_cascade/results/coset_invariant_validation/validation.log`](experiments/gf256_cascade/results/coset_invariant_validation/validation.log)
- 构建目标：`cmake --build <build-dir> --target test_gf256_cascade_coset_invariants`；运行时将矩阵构建目录作为唯一参数。
- 主实验矩阵目录 `g2_basis_path_n12_300frames/path/build`：`HX=2134×6984`、`HZ=2716×6984`；非正交 HZ 行 `0/2716`、稳定子判定拒绝 `0/2716`、改变陪集键 `0/2716`。该目录没有单独保存 H2 文件，测试按仿真缺省单位缩放。
- 两套矩阵均由尺寸推得内维数 `11`。第二矩阵种子 `20261019` 读取保存的 `H2=1×12` 并按仿真代码计算缩放；三项逐行异常均为 `0/2716`。另各自找到一个非平凡逻辑见证并正确区分。
- 结论范围：验证了当前有限域矩阵及级联 Z 陪集映射的生成行结构，并确认一个非平凡逻辑类能与稳定子区分；不证明候选列表覆盖率，也不代表 FER 已改善。

### Safe joint projection 配对复测

- 目录：[`experiments/gf256_cascade/results/joint_projection_safe_same_frames100_matrix20261018_pD0p072/`](experiments/gf256_cascade/results/joint_projection_safe_same_frames100_matrix20261018_pD0p072/)
- 参数：矩阵种子 `20261018`，`p_D=0.072`，相同噪声种子 `20313000–20313099`，100 帧。
- 结果：tuned/safe 逻辑 FER 均 `10/100`；救回/退步 `0/0`；safe 的 Z syndrome failure 从 10 降至 0，但漏检逻辑失败从 0 增至 10；CPU 增加约 3.8%。

### tuned + global restart-belief 配对对照（已完成 pilot）

- 计划比较：`conditioned_robust_tuned` 对 `conditioned_robust_tuned_restart_belief_global`。
- 目标：新方法保留 tuned 主译码设置，只在主译码 syndrome 不通过时增加 global restart-belief 搜索；使用同一矩阵、同一信道参数和同一组噪声帧配对比较。
- 结果：矩阵种子 `20261018` 对应的 12 内长 path 矩阵，`p_D=0.07125`，新种子 `20320000–20320003`，20 帧。tuned 与 global restart 两者均 `1/20`；救回/退步 `0/0`，McNemar `p=1.0`；CPU `10.06→14.13 s/frame`。同一批失败帧上 4×4 global 重启没有找到有效候选。该 pilot 不支持扩大该分支。
- 结果目录：`experiments/gf256_cascade/results/paired_tuned_restart_global_pilot20_seed20320000/`。

### Syndrome-flip 译码原型和配对实验

- 文献动机：[BP-SF 论文](https://arxiv.org/abs/2507.00254)通过不可靠位生成 syndrome-flip 试探；当前实现用每个 GF(256) 符号的物理位后验不确定度作排序，属于适配原型，不是逐项复现原文的 oscillation ranking。
- 算法模式：`conditioned_robust_tuned_syndrome_flip`；在 `run_paired_decoder_study.py` 注册。有限域特征为特征 2，所用移位关系为 (H_X(e'+t)=s+H_Xt+H_Xt=s)。只有校验原始 syndrome 通过的候选才参与 log-likelihood 选择。
- 小样本：预选失败批次 `20268300` 的 5 帧中，tuned `1/5`，4 位、重量≤2 的 syndrome-flip `0/5`；该失败帧找到 1 个原 syndrome-valid 候选。此数据用于机制调试，不用于 FER 估计。
- 独立 pilot：`results/paired_tuned_syndrome_flip_pilot20_seed20322000/`，两者均 `0/20`，没有失败事件可比较。
- 独立 100 帧配对：`results/paired_tuned_syndrome_flip_holdout100_seed20323000/`，tuned `7/100`（Wilson 95% `[0.0343, 0.1375]`），syndrome-flip `5/100`（`[0.0215, 0.1118]`）；配对救回 2、退步 0、McNemar `p=0.5`。Z syndrome failure `7→5`，漏检逻辑失败均为 0。CPU `9.76→14.93 s/frame`。
- 第一组独立 300 帧配对：`results/paired_tuned_syndrome_flip_validation300_seed20324000/`，tuned `15/300`（Wilson 95% `[0.0305,0.0808]`），syndrome-flip `13/300`（`[0.0255,0.0727]`）；救回 2、退步 0、McNemar `p=0.5`；CPU `9.84→13.67 s/frame`。
- 第二组独立 300 帧复制：`results/paired_tuned_syndrome_flip_replication300_seed20325000/`，tuned `14/300`（`[0.0280,0.0768]`），syndrome-flip `13/300`（`[0.0255,0.0727]`）；救回 1、退步 0、McNemar `p=1.0`；CPU `10.01→13.65 s/frame`。
- 三组独立样本合计：tuned `36/700=5.14%`，syndrome-flip `31/700=4.43%`；配对救回 5、退步 0，精确 McNemar `p=0.0625`。这是描述性合并，不能作为已证明显著增益的结论。

### 其他近期负结果索引

- 原始 joint projection：[`experiments/gf256_cascade/results/joint_projection_paired_dev100_matrix20261018_pD0p072/`](experiments/gf256_cascade/results/joint_projection_paired_dev100_matrix20261018_pD0p072/)
- SI 100 帧配对开发集：[`experiments/gf256_cascade/results/si_paired_dev100_matrix20261018_pD0p072/`](experiments/gf256_cascade/results/si_paired_dev100_matrix20261018_pD0p072/)
- overcomplete 外层校验 pilot：[`experiments/gf256_cascade/results/overcomplete_outer_paired_pilot5_matrix20261018_pD0p07125/`](experiments/gf256_cascade/results/overcomplete_outer_paired_pilot5_matrix20261018_pD0p07125/)
- 更完整的文献与译码路线记录：[`experiments/gf256_cascade/GUIDED_DECODING_RESEARCH.md`](experiments/gf256_cascade/GUIDED_DECODING_RESEARCH.md)

## 维护约定

- 每次重要实现或实验后更新本文件，保持“当前状态、已完成、未完成、当前问题、最近结果”一致。
- 只有实现及对应验证完成才勾选 `[x]`；区分精确性验证、机制诊断、开发集 FER 和独立验证。
- 新实验记录矩阵/信道/噪声种子、帧数、逻辑 FER、配对救回与退步、综合征失败、漏检逻辑失败、运行时间和结果目录。
- 性能改善需用同一批噪声帧配对比较；只有经过独立数据验证后，才能表述为具有泛化证据的提升。
