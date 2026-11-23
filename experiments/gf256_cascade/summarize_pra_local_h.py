"""Aggregate and plot the biased-noise companion-code experiments."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest

from run_pra_local_h_batches import OUT, ROOT


def read_runs(directory):
    groups = defaultdict(list)
    for path in sorted(directory.glob('*_a*_p*_s*.json')):
        key, _ = path.stem.rsplit('_s', 1)
        data = json.loads(path.read_text())
        if len(data['records']) != data['frames']:
            raise ValueError(f'Incomplete records in {path}')
        groups[key].append(data)
    return groups


def summarize(groups):
    rows = []
    for key, runs in sorted(groups.items()):
        first = runs[0]
        settings = ('p_total','alpha','matrix_manifest','hadamards','max_iters')
        if any(any(r[k] != first[k] for k in settings) for r in runs):
            raise ValueError(f'Incompatible group {key}')
        if len({r['seed'] for r in runs}) != len(runs):
            raise ValueError(f'Duplicate seed in {key}')
        records = [rec for r in runs for rec in r['records']]
        n = len(records)
        k = sum(not r['logical_ok'] for r in records)
        w = np.array([r['physical_residual_weight'] for r in records])
        block_errors = int(np.count_nonzero(w))
        fer_ci = binomtest(k,n).proportion_ci()
        wer_ci = binomtest(block_errors,n).proportion_ci()
        length = first['matrix_manifest']['physical_qubits']
        rows.append(dict(group=key, p=first['p_total'], alpha=first['alpha'], n_phys=length,
                         rate=first['matrix_manifest']['quantum_rate'], frames=n, failures=k, fer=k/n,
                         fer_lo=fer_ci.low, fer_hi=fer_ci.high, wer_events=block_errors,
                         wer=block_errors/n, wer_lo=wer_ci.low, wer_hi=wer_ci.high,
                         qber=float(w.mean()/length),
                         staged_failures=sum(not r['staged_syndrome_ok'] for r in records),
                         projected_failures=sum(not r['projected_syndrome_ok'] for r in records),
                         undetected=sum(r['undetected'] for r in records),
                         mean_seconds=float(np.mean([r['seconds'] for r in records])),
                         zero_upper95=float(-np.expm1(np.log(.05)/n)) if not block_errors else '',
                         seeds=';'.join(str(r['seed']) for r in runs)))
    return rows


def paired(groups):
    comparisons = {}
    def accumulate(key, baseline, run):
        for name in ('noise_sha256','p_total','alpha','frames','matrix_manifest','max_iters'):
            if run[name] != baseline[name]:
                raise ValueError(f'Pair mismatch {name} for {key}')
        row = comparisons.setdefault(key, dict(frames=0, baseline_logical_failures=0,
                                               local_logical_failures=0, rescued=0, regressed=0))
        row['frames'] += run['frames']
        for old, new in zip(baseline['records'],run['records'],strict=True):
            if old['frame'] != new['frame']:
                raise ValueError(f'Frame order differs for {key}')
            row['baseline_logical_failures'] += not old['logical_ok']
            row['local_logical_failures'] += not new['logical_ok']
            row['rescued'] += not old['logical_ok'] and new['logical_ok']
            row['regressed'] += old['logical_ok'] and not new['logical_ok']

    for m in (1,2):
        key = f'paired_m{m}_a100_p0.07'
        for run in groups.get(key, []):
            baseline_path = ROOT/'results/pra_alpha100_20261118/swapped'/f'm{m}_p0.07_s{run["seed"]}.json'
            baseline = json.loads(baseline_path.read_text())
            accumulate(str(m),baseline,run)
        without_tail = {r['seed']:r for r in groups.get(f'tail_ablation_m{m}_a100_p0.1',[])}
        for run in groups.get(f'curve_m{m}_a100_p0.1',[]):
            baseline = without_tail.get(run['seed'])
            if baseline is not None:
                accumulate(f'tail_m{m}',baseline,run)
    return comparisons


def plot(directory, rows):
    with (ROOT/'pra_fig10_reference.csv').open() as file:
        reference = list(csv.DictReader(file))
    fig, axes = plt.subplots(1, 3, figsize=(9,3.25), layout='constrained', sharex=True)
    for ax, key, label in zip(axes, ('fer','wer','qber'), ('Logical FER','Physical WER','Physical QBER')):
        ax.set(xlabel=r'Physical Pauli error probability $p$', ylabel=label,
               yscale='log', xlim=(.055,.137), ylim=(1e-6,1.5))
        ax.grid(alpha=.2)
        if key != 'fer':
            for size,style in ((19014,'--'),(38028,':')):
                values = [r for r in reference if int(r['N_legend'])==size]
                ax.plot([float(r['p']) for r in values],
                        [float(r[key+'_digitized']) for r in values], style, color='.45',
                        lw=1.15, label=f'Fuentes {size:,}')
        for m,marker,color in ((1,'o','#146f9c'),(2,'s','#208375')):
            matching = [r for r in rows if r['alpha']==100 and f'_m{m}_' in r['group']
                        and r['group'].startswith(('curve_','independent_','extension_'))]
            for i,r in enumerate(sorted(matching,key=lambda r:r['p'])):
                y = r[key]
                legend = f'Local H, matrix {m}' if i==0 else None
                if y==0:
                    bound = r['zero_upper95'] if key!='fer' or r['zero_upper95'] else -np.expm1(np.log(.05)/r['frames'])
                    ax.plot(r['p'],bound,'v',color=color,mfc='none',ms=7 if m==1 else 4,
                            label=legend)
                else:
                    if key == 'qber':
                        ax.plot(r['p'],y,marker,color=color,ms=4,label=legend)
                    else:
                        ax.errorbar(r['p'],y,yerr=[[y-r[key+'_lo']],[r[key+'_hi']-y]],
                                    fmt=marker,color=color,capsize=2,ms=4,label=legend)
        # Compare the primary same-noise anchors against global H.
        for m,marker in ((1,'D'),(2,'P')):
            baseline = sorted((ROOT/'results/pra_alpha100_20261118/swapped').glob(f'm{m}_p0.07_*.json'))
            samples = [json.loads(path.read_text()) for path in baseline]
            fname={'fer':'logical_fer','wer':'wer_physical','qber':'qber_pauli'}[key]
            if samples:
                y=np.average([r[fname] for r in samples],weights=[r['frames'] for r in samples])
                ax.plot(.07,y,marker,color='#a23f33',ms=5,label=f'Global H, matrix {m}')
        ax.legend(loc='lower right',fontsize=6,frameon=False)
    fig.suptitle(r'Biased Pauli noise, $alpha=100$; triangles show 95% upper bounds', fontsize=10)
    fig.savefig(directory/'same_channel_metrics.pdf')
    fig.savefig(directory/'same_channel_metrics.png',dpi=200)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=OUT)
    args=parser.parse_args()
    groups=read_runs(args.directory)
    rows=summarize(groups)
    if not rows:
        parser.error('No completed batches')
    with (args.directory/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.directory/'paired.json').write_text(json.dumps(paired(groups),indent=2)+'\n')
    plot(args.directory,rows)
    for row in rows:
        print(row['group'],row['frames'],row['failures'],row['wer_events'],f"QBER={row['qber']:.6g}")


if __name__=='__main__':
    main()
