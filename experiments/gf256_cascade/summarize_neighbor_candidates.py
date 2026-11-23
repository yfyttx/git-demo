"""Compare baseline and shared-symbol candidate coverage on selected frames."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    folder = args.directory
    status = json.loads((folder / 'summary.json').read_text())
    rows = list(csv.DictReader((folder / 'projection_coverage.csv').open()))
    keys = sorted({(r['seed'], r['trial']) for r in rows})
    frames = []
    for seed, trial in keys:
        rr = [r for r in rows if (r['seed'], r['trial']) == (seed, trial)]
        coverage = {}
        for kind in ('baseline', 'neighbor'):
            subset = [r for r in rr if r['stage'].endswith('_neighbor') == (kind == 'neighbor')]
            assert subset, (seed, trial, kind)
            coverage[kind] = {int(r['coordinate']) for r in subset if int(r['coset_covered'])}
        assert coverage['baseline'] <= coverage['neighbor']
        frames.append(dict(seed=int(seed), trial=int(trial),
                           baseline_covered=len(coverage['baseline']),
                           neighbor_covered=len(coverage['neighbor']),
                           gained=sorted(coverage['neighbor']-coverage['baseline'])))
    result = dict(replay_matches=status['replay_matches'], frames=frames,
                  limitation='Offline lists, selected failures, increased candidate/computation budget; no FER claim')
    (folder / 'neighbor_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    text = ['# 相邻投影共享符号候选生成：先导结果', '',
            '固定原矩阵、信道与噪声帧；真值仅用于候选生成结束后的覆盖评估。', '',
            '方法：每阶段原列表按原投影后验评分，取前4个候选，均匀形成符号边缘分布；以0.5权重与原外层外信息混合。利用path内层链前向/后向消息，为每个目标投影生成排除自身消息的先验，重新运行外层BP和OSD-1。执行两轮并保留原列表。', '',
            '此为候选反馈启发式：列表本身包含循环图信息，混合权重不保证消除重复计数，不能称为精确联合BP。每次OSD预算为32自由符号、每符号4替代值；两轮增加候选与计算量，若出现增益仍需同预算消融。', '',
            '| seed | trial | 原列表跨阶段覆盖/11 | 新列表跨阶段覆盖/11 | 新覆盖坐标 |',
            '|---|---:|---:|---:|---|']
    for r in frames:
        text.append(f"| {r['seed']} | {r['trial']} | {r['baseline_covered']} | {r['neighbor_covered']} | {r['gained']} |")
    text += ['', f"回放结果标签一致：{status['replay_matches']}/{status['selected']}。新增候选逐一检查外层综合征。", '',
             '这不是FER实验。即使11个坐标均覆盖，也仍需选择相互一致的整帧候选并检查稳定子等价。当前诊断未接入在线选择，因此未宣称修复任何帧。', '']
    (folder / 'report.md').write_text('\n'.join(text))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
