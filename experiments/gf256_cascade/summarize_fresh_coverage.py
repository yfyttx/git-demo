"""Summarize completed offline coverage diagnostics without estimating FER."""
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "experiments/gf256_cascade/results/fresh34_logical_coverage_audit_20261117"


def main():
    summary = json.loads((OUT / "summary.json").read_text())
    rows = list(csv.DictReader((OUT / "projection_coverage.csv").open()))
    projections = list(csv.DictReader((OUT / "final_projections.csv").open()))
    outcomes = list(csv.DictReader((OUT / "replay_outcomes.csv").open()))
    failed = [r for r in outcomes if r["group"] == "persistent_z"]
    assert len(failed) == 34 and summary["replay_matches"] == 37
    assert all(int(r["covered"]) <= int(r["coset_covered"]) for r in rows)
    assert all(r["truth_syndrome"] == "1" for r in projections)
    details = []
    for frame in failed:
        rr = [r for r in rows if (r["seed"], r["trial"]) == (frame["seed"], frame["trial"])]
        missing = [a for a in range(11) if not any(int(r["coset_covered"]) for r in rr if int(r["coordinate"]) == a)]
        exact_missing = [a for a in range(11) if not any(int(r["covered"]) for r in rr if int(r["coordinate"]) == a)]
        item = dict(seed=frame["seed"], trial=frame["trial"],
                    union_missing_outer_cosets=len(missing), union_missing_exact=len(exact_missing),
                    missing_coordinates=" ".join(map(str,missing)),
                    contiguous_missing_runs=sum(a-1 not in missing for a in missing))
        for stage in ("initial", "retry", "contender", "layered10", "layered20"):
            stage_rows = [r for r in rr if r["stage"] == stage]
            if stage_rows:
                assert len(stage_rows) == 11
            item[stage+"_covered"] = sum(int(r["coset_covered"]) for r in stage_rows) if stage_rows else "not_run"
        details.append(item)
    with (OUT / "frame_mechanisms.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(details[0]))
        writer.writeheader(); writer.writerows(sorted(details,key=lambda r:(r["seed"],r["trial"])))
    histogram=Counter(r["union_missing_outer_cosets"] for r in details)
    distance=[int(r["min_symbol_distance"]) for r in rows if r["group"]=="persistent_z" and r["coset_covered"]=="0"]
    final_bad=Counter()
    for frame in failed:
        pp=[r for r in projections if (r["seed"],r["trial"])==(frame["seed"],frame["trial"])]
        assert len(pp)==11
        final_bad[sum(r["decoded_syndrome"]=="0" for r in pp)]+=1
    aggregate=dict(union_missing_histogram=dict(sorted(histogram.items())),
                   equivalent_but_not_exact=sum(r["coset_covered"]=="1" and r["covered"]=="0" for r in rows),
                   failure_exact_x=sum(r["replay_exact_x"]=="1" for r in failed),
                   single_contiguous_missing_run=sum(r["contiguous_missing_runs"]==1 for r in details),
                   min_missing_distance=min(distance) if distance else None,
                   max_missing_distance=max(distance) if distance else None,
                   final_invalid_projection_histogram=dict(sorted(final_bad.items())))
    (OUT/"coverage_summary.json").write_text(json.dumps(aggregate,indent=2)+"\n")
    print(json.dumps(aggregate,indent=2))
    text=["# 新1000帧中的34个Z综合征失败：候选覆盖审计", "", "日期：2026-11-17。此处为预选困难帧机制回放，不是新FER实验。", "",
          "原始独立样本固定GF(256)、12列path、矩阵种子20261019、p_D=0.07125。组合译码器35/1000逻辑失败中，34帧Z综合征失败、1帧漏检。此次回放全部34个Z失败帧，并加入各1个漏检、救回和成功对照。37个回放的逻辑及综合征标签均与冻结实验一致。", "",
          "34个Z失败帧的X错误全部精确恢复；每帧8次syndrome-flip都未产生有效候选。故这批失败不能归因于X估计错误传播。", "",
          "## 各阶段候选覆盖", "", "离线OSD-1使用联合后验、32个自由符号、每符号4个替代值，最多129个候选/投影。候选不传回译码选择。真值只用于译码后的覆盖检查。", "",
          "| 阶段 | 有审计列表的失败帧 | 正确外层陪集覆盖/投影数 | 11个坐标全部覆盖的帧数 |", "|---|---:|---:|---:|"]
    for stage,value in summary["stages"].items():
        text.append(f"| {stage} | {value['frames']} | {value['outer_coset_covered']}/{value['projections']} | {value['frames_all_outer_cosets_covered']} |")
    text += ["", "## 合并各阶段列表后仍缺失的坐标", "", "| 缺失正确外层陪集的坐标数 | 帧数 |", "|---:|---:|"]
    text += [f"| {n} | {count} |" for n,count in sorted(histogram.items())]
    text += ["",f"共有{sum(r['union_missing_outer_cosets']>0 for r in details)}/34帧即使合并各阶段列表仍有缺失。{aggregate['single_contiguous_missing_run']}/34帧的缺失坐标构成一段连续链。连续性是观察结果，尚未与随机缺失模型比较，不能单独证明因果。",
             "",f"所有审计记录中，正确外层陪集已覆盖但真实投影未覆盖的记录数为{aggregate['equivalent_but_not_exact']}。未覆盖记录的最近候选仍与真实投影相差{aggregate['min_missing_distance']}–{aggregate['max_missing_distance']}个GF(256)符号。", "",
             "## 对后续工作的含义", "", "若某投影列表没有正确外层陪集，仅重排当前列表无法恢复稳定子等价的结果。下一步应优先改候选生成中心及相邻投影的共享符号消息，先用开发帧验证候选覆盖，再用全新噪声样本验证逻辑FER。不能用本次真值选择在线分支，也不能把本次回放当作译码增益。",
             "", "这些结果仅针对当前固定预算及阶段后验，不证明更大的列表或其他译码方法无效。候选覆盖检查的线性代数依据和完整局限见结果目录interpretation.md。", "",
             "## 数据与复现", "", "- 结果目录：`results/fresh34_logical_coverage_audit_20261117/`。", "- 逐坐标：`projection_coverage.csv`；逐帧：`frame_mechanisms.csv`；最终投影：`final_projections.csv`。", "- 原始基线：`results/paired_three_modes_fresh1000_matrix20261019/summary.csv`。", "- `config.json`记录输入哈希和选帧；保存实际运行可执行文件与源码。", "- 初始仅精确覆盖的探索回放保存在`results/fresh34_coverage_audit_20261117/`，已停止；它不是本报告的完整数据来源。", "",
             "```bash", "python3 experiments/gf256_cascade/audit_fresh_failure_coverage.py --workers 4", "python3 experiments/gf256_cascade/summarize_fresh_coverage.py", "```", ""]
    (ROOT/"experiments/gf256_cascade/Z_FAILURE_COVERAGE_20261117.md").write_text("\n".join(text))


if __name__=="__main__":
    main()
