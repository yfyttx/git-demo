# 单位内码退化与 Kasai 原码复现

## 构造和约定

使用当前外层矩阵 seed=20261019，GF(256)、J=2、L=6、P=97、sigma=35、tau=2。新生成的直接Kasai程序Gamma/Delta与已保存外层矩阵逐元素比较，差异均为0。

G2=I8，H2为0×8空校验矩阵，选择矩阵S=I8。新伴随块构造退化成外层展开；HX、HZ均1552×4656，秩均1552，物理码长4656、逻辑维数1552、RQ=1/3。

必须区分约定：当前级联入口是Delta的列乘法展开与Gamma的转置展开；Kasai专用入口按Gamma列乘法展开、Delta转置展开。二者的具体二进制表示不相同，不把它们的不同样本结果当作同帧配对。

## 单位内码与直接外层的同帧核对

`run_identity_outer_comparison.py`固定当前级联约定，对同一物理退极化噪声分别从二进制完整矩阵和GF(256)外层生成综合征，两次调用同一确定性外层BP核对输出。X侧使用trace-dual基映射，Z侧使用多项式基；每侧先验按对应的物理位重量计算。

100帧已完成：矩阵差异0、综合征差异0、译码输出差异0。两种入口的Delta分量块错误均6/100，Gamma分量均6/100。

这验证矩阵和综合征映射的退化一致性，不是两套独立BP实现的交叉验证。

```bash
python3 experiments/gf256_cascade/construct_companion_cascade.py --identity-inner
OPENBLAS_NUM_THREADS=1 python3 experiments/gf256_cascade/run_identity_outer_comparison.py --frames 100 --fm 0.0475 --max-iters 50
```

## 直接Kasai算法先导结果

固定同一矩阵、seed20261019、每点200帧，两个分量分别非二进制sum-product译码，FFT/FWHT校验更新，不使用X/Z相关性或级联回退。最多50轮；该上限是当前实现设置，未确认与原论文仿真上限一致。

|fm|pD|Gamma块错误|Delta块错误|Gamma综合征失败|Delta综合征失败|
|---|---|---|---|---|---|
|0.04|0.06|0/200|1/200 (0.5%)|0|1|
|0.0475|0.07125|18/200 (9%)|10/200 (5%)|18|10|

块错误定义为估计向量不等于真实向量，匹配分量码块错误统计口径；不是稳定子等价意义下完整量子逻辑FER。0/200不能断言真实错误率为零。此参数与论文GF(256)、n=4656、RQ=1/3曲线相符，但论文未公开该曲线的完整系数/随机种子，200帧也不足以复现低错误率区间，故不能声称已完全复现原图。

编译运行：
```bash
cmake -S . -B /tmp/kasai_identity_build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/kasai_identity_build --target kasai_paper -j4
/tmp/kasai_identity_build/kasai_paper 200 20261019 0.0475 8 97 35 2 6
```

数据在results/kasai_direct_20261117/；同帧核对数据在results/identity_outer_comparison_fm0.0475_seed20261118/。
