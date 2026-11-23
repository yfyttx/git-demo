#!/usr/bin/env python3
"""Plot the current companion-code pilot using PRA-style physical metrics plus logical FER.

This is an internal, exploratory two-panel plot, not yet a final manuscript figure.
Run from the repository root: python3 experiments/gf256_cascade/plot_pra_metric_reference_sweep.py
"""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'experiments/gf256_cascade/results/pra_metric_reference_sweep_20261118'
SOURCES = [
    ('experiments/gf256_cascade/results/companion_pra_scan_pD0p0600_seed20261160.json',),
    ('experiments/gf256_cascade/results/companion_pra_scan_pD0p0650_seed20261161.json',),
    ('experiments/gf256_cascade/results/companion_pra_scan_pD0p0685_seed20261162.json',
     'experiments/gf256_cascade/results/companion_pra_sweep_pD0p0685_seed20261166.json',
     'experiments/gf256_cascade/results/companion_pra_sweep_pD0p0685_seed20261200.json'),
    ('experiments/gf256_cascade/results/companion_pra_sweep_pD0p0675_seed20261165.json',),
    ('experiments/gf256_cascade/results/companion_pra_scan_pD0p0705_seed20261163.json',
     'experiments/gf256_cascade/results/companion_pra_sweep_pD0p0705_seed20261167.json',
     'experiments/gf256_cascade/results/companion_pra_sweep_pD0p0705_seed20261201.json'),
    ('experiments/gf256_cascade/results/companion_pra_metrics_seed20261125.json',
     'experiments/gf256_cascade/results/companion_pra_metrics_seed20261126.json'),
]
SECOND_MATRIX_SOURCES = [
    'experiments/gf256_cascade/results/companion_pra_matrix2_pD0650_seed20261301.json',
    'experiments/gf256_cascade/results/companion_pra_matrix2_pD0p0685_seed20261300.json',
    'experiments/gf256_cascade/results/companion_pra_matrix2_pD0705_seed20261302.json',
    'experiments/gf256_cascade/results/companion_pra_matrix2_pD07125_seed20261303.json',
]


def wilson(k, n, confidence=.95):
    z = 1.959963984540054
    p = k/n
    den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return max(0., center-half), min(1., center+half)


def exact_zero_upper(n, confidence=.95):
    return 1 - (1-confidence)**(1/n)


def read_rows(paths):
    return [json.loads(Path(ROOT / p).read_text()) for p in paths]


def build_rows():
    rows = []
    for paths in SOURCES:
        runs = read_rows(paths)
        n = sum(run['frames'] for run in runs)
        frames = [record for run in runs for record in run['records']]
        p = runs[0]['p_D']
        if any(run['p_D'] != p or run['matrix_manifest'] != runs[0]['matrix_manifest'] for run in runs):
            raise ValueError(f'Incompatible source runs at p_D={p}')
        counts = {
            key: sum(not f[field] for f in frames)
            for key, field in [('logical_fer','logical_ok'), ('x_logical_fer','x_logical_ok'),
                               ('z_logical_fer','z_logical_ok'), ('x_syndrome_failure','x_syndrome_ok'),
                               ('z_syndrome_failure','z_syndrome_ok')]
        }
        nphys = 9*runs[0]['matrix_manifest']['outer_columns']
        qber_frame = {
            'qber_x': np.array([f['x_residual_weight']/nphys for f in frames]),
            'qber_z': np.array([f['z_residual_weight']/nphys for f in frames]),
            'qber_pauli': np.array([f['pauli_residual_weight']/nphys for f in frames]),
        }
        rng = np.random.default_rng(771000 + round(p*1e6))
        boot = 5000
        for key, vals in qber_frame.items():
            samples = rng.integers(0, n, size=(boot,n))
            means = vals[samples].mean(axis=1)
            counts[key] = float(vals.mean())
            counts[key+'_ci_low'], counts[key+'_ci_high'] = map(float, np.quantile(means,[.025,.975]))
        for key in ('logical_fer','x_logical_fer','z_logical_fer','x_syndrome_failure','z_syndrome_failure'):
            count=counts[key]
            counts[key+'_events']=count
            counts[key]=count/n
            counts[key+'_ci_low'],counts[key+'_ci_high']=wilson(count,n)
            counts[key+'_zero_upper']=exact_zero_upper(n) if count==0 else ''
        counts.update(p_D=p, f_m=runs[0]['f_m'], frames=n,
            physical_wer_x=sum(f['x_physical_block_error'] for f in frames)/n,
            physical_wer_z=sum(f['z_physical_block_error'] for f in frames)/n,
            physical_wer_pauli=sum(f['physical_block_error'] for f in frames)/n,
            seeds=';'.join(str(run['seed']) for run in runs))
        rows.append(counts)
    return sorted(rows, key=lambda row: row['p_D'])


def plot(rows):
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path=OUT/'source_data.csv'
    with csv_path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)
    p=np.array([r['p_D'] for r in rows])
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.0),constrained_layout=True)
    ax=axes[0]
    colors={'logical_fer':'#153e75','x_logical_fer':'#118a8a','z_logical_fer':'#b5473a'}
    markers={'logical_fer':'o','x_logical_fer':'s','z_logical_fer':'^'}
    labels={'logical_fer':'Overall logical FER','x_logical_fer':'X logical failure','z_logical_fer':'Z logical failure'}
    for key in colors:
        y=np.array([r[key] for r in rows],float)
        lo=np.array([r[key+'_ci_low'] for r in rows],float)
        hi=np.array([r[key+'_ci_high'] for r in rows],float)
        pos=y>0
        if np.any(pos):
            ax.errorbar(p[pos],y[pos],yerr=[np.maximum(0,y[pos]-lo[pos]),np.maximum(0,hi[pos]-y[pos])],
                color=colors[key],marker=markers[key],lw=1.35,ms=4,capsize=2,label=labels[key])
        for i in np.flatnonzero(~pos):
            upper=rows[i][key+'_zero_upper']
            ax.scatter([p[i]],[upper],marker='v',facecolors='none',edgecolors=colors[key],s=34,zorder=4)
            ax.annotate('0 events', (p[i],upper), xytext=(3,5), textcoords='offset points',fontsize=6,color=colors[key])
    ax.set_yscale('log'); ax.set_ylim(5e-4,0.2)
    ax.set_xlabel(r'Gross depolarizing probability $p_D$')
    ax.set_ylabel('Frame failure probability')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=6,loc='lower right',frameon=False)
    ax.text(.02,.98,'(a)',transform=ax.transAxes,va='top',fontweight='bold')

    ax=axes[1]
    qcolors={'qber_x':'#118a8a','qber_z':'#b5473a','qber_pauli':'#153e75'}
    qmarkers={'qber_x':'s','qber_z':'^','qber_pauli':'o'}
    qlabels={'qber_x':'Residual X QBER','qber_z':'Residual Z QBER','qber_pauli':'Pauli-union QBER'}
    for key in qcolors:
        y=np.array([r[key] for r in rows],float)
        lo=np.array([r[key+'_ci_low'] for r in rows],float)
        hi=np.array([r[key+'_ci_high'] for r in rows],float)
        pos=y>0
        ax.errorbar(p[pos],y[pos],yerr=[np.maximum(0,y[pos]-lo[pos]),np.maximum(0,hi[pos]-y[pos])],
            color=qcolors[key],marker=qmarkers[key],lw=1.35,ms=4,capsize=2,label=qlabels[key])
        for i in np.flatnonzero(~pos):
            ax.scatter([p[i]],[max(1e-5,hi[i])],marker='v',facecolors='none',edgecolors=qcolors[key],s=34)
    ax.set_yscale('log'); ax.set_ylim(1e-5,.2)
    ax.set_xlabel(r'Gross depolarizing probability $p_D$')
    ax.set_ylabel('Residual physical-qubit fraction')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=6,loc='lower right',frameon=False)
    ax.text(.02,.98,'(b)',transform=ax.transAxes,va='top',fontweight='bold')
    fig.savefig(OUT/'pra_metric_reference_pilot.pdf')
    fig.savefig(OUT/'pra_metric_reference_pilot.svg')
    fig.savefig(OUT/'pra_metric_reference_pilot.png',dpi=240)
    plt.close(fig)
    note=("Internal pilot figure only; not a final manuscript figure. Fixed matrix 20261019, "
          "conditioned Z prior, 100 BP iterations, no X fallback. Panel (a) uses Wilson 95% intervals; "
          "zero-event points are shown at the one-sided exact 95% upper bound and labeled. Panel (b) "
          "uses frame-cluster bootstrap percentile 95% intervals (5,000 resamples); open markers indicate "
          "zero observed residuals. Source data are in source_data.csv. Physical QBER is not logical FER.")
    (OUT/'figure_notes.md').write_text(note+'\n')
    return csv_path


def plot_matrix_comparison(primary_rows):
    second_runs = read_rows(SECOND_MATRIX_SOURCES)
    z = 1.959963984540054
    secondary = []
    for run in second_runs:
        n = run['frames']
        k = run['logical_failures']
        p = k/n
        den = 1 + z*z/n
        center = (p + z*z/(2*n))/den
        half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
        secondary.append(dict(matrix='outer path2',p_D=run['p_D'],frames=n,failures=k,
                              fer=p,ci_low=max(0.,center-half),ci_high=min(1.,center+half),seed=run['seed']))
    primary = [dict(matrix='outer path1',p_D=r['p_D'],frames=r['frames'],failures=int(r['logical_fer_events']),
                    fer=float(r['logical_fer']),ci_low=float(r['logical_fer_ci_low']),
                    ci_high=float(r['logical_fer_ci_high']),seed=r['seeds'])
               for r in primary_rows if r['p_D'] in {0.065,0.0675,0.0685,0.0705,0.07125}]
    comparison = primary + secondary
    with (OUT/'two_matrix_fer_data.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(comparison[0].keys()))
        writer.writeheader(); writer.writerows(comparison)
    fig,ax=plt.subplots(figsize=(4.5,3.2),constrained_layout=True)
    styles=[('outer path1','#153e75','o',primary),('outer path2','#b5473a','s',secondary)]
    for label,color,marker,series in styles:
        series=sorted(series,key=lambda r:r['p_D'])
        x=np.array([r['p_D'] for r in series]); y=np.array([r['fer'] for r in series])
        lo=np.array([r['ci_low'] for r in series]); hi=np.array([r['ci_high'] for r in series])
        ax.errorbar(x,y,yerr=[y-lo,hi-y],color=color,marker=marker,lw=1.35,ms=4,capsize=2,label=label)
    ax.set_yscale('log'); ax.set_ylim(5e-4,.12)
    ax.set_xlabel(r'Gross depolarizing probability $p_D$')
    ax.set_ylabel('Overall logical FER')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=8,frameon=False)
    ax.set_title('Fixed decoder, two outer-matrix realizations',fontsize=9)
    fig.savefig(OUT/'pra_metric_two_matrix_fer.pdf')
    fig.savefig(OUT/'pra_metric_two_matrix_fer.svg')
    fig.savefig(OUT/'pra_metric_two_matrix_fer.png',dpi=240)
    plt.close(fig)
    (OUT/'two_matrix_figure_notes.md').write_text(
        'Exploratory comparison only. Same [9,8,2] inner code, decoder, channel model, and 100 iterations; '
        'the second matrix has 1000 frames per point. Error bars are Wilson 95% intervals. '
        'One matrix realization per outer construction is not enough to claim a general matrix-family effect. '
        'Data are in two_matrix_fer_data.csv.\n')


if __name__=='__main__':
    rows=build_rows(); path=plot(rows); plot_matrix_comparison(rows)
    print(f'Wrote {path.parent}: {len(rows)} channel points plus two-matrix comparison')
    for r in rows:
        print(f"pD={r['p_D']:.5f} n={r['frames']} logical={r['logical_fer']:.5f} "
              f"({r['logical_fer_events']} events) X/Z syndrome="
              f"{r['x_syndrome_failure_events']}/{r['z_syndrome_failure_events']} "
              f"QBER={r['qber_x']:.5f}/{r['qber_z']:.5f}/{r['qber_pauli']:.5f}")
