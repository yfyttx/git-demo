"""Aggregate completed alpha=100 shards, including paired, same-noise comparisons.

Usage: python3 experiments/gf256_cascade/plot_pra_alpha100.py [--plan plan.csv]
Optional plan columns: protocol,matrix,p,frames. Omitting the plan reports only
observed points. Matching log files are counted as pending until a full JSON exists.
"""
import argparse
from collections import defaultdict
import csv
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest, beta

HERE = Path(__file__).resolve().parent
DEFAULT_RUN = HERE / 'results/pra_alpha100_20261118'
PROTOCOLS = {'original': 'Original orientation',
             'swapped': 'Swapped + parity information',
             'swapped_plain': 'Swapped, plain'}
PATTERN = re.compile(r'm([12])_p(0\.\d+)_s(\d+)\.(json|log)$')
CI_COUNT_KEYS = ['logical', 'x_logical', 'z_logical', 'x_syndrome', 'z_syndrome', 'undetected', 'physical_wer']


def load_csv(path):
    with path.open(newline='') as fp:
        return list(csv.DictReader(fp))


def export(path, rows, fields=None):
    with path.open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=fields or list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read_shards(run_dir):
    frames = defaultdict(dict)
    shard_log = []
    signatures = {}
    for protocol in PROTOCOLS:
        folder = run_dir / protocol
        stems = sorted({path.stem for path in folder.glob('m*') if PATTERN.fullmatch(path.name)})
        for stem in stems:
            path = folder / (stem + '.json')
            matrix, probability, seed, _ = PATTERN.fullmatch(path.name).groups()
            point = protocol, int(matrix), float(probability)
            frame_set = frames[point]
            metadata = dict(protocol=protocol, matrix=matrix, p=probability, seed=seed,
                            file=str(path.relative_to(run_dir)), records=0, sha256='', source_hashes='')
            shard_log.append(metadata)
            if not path.exists():
                metadata['status'] = 'log only'
                continue
            raw = path.read_bytes()
            try:
                run = json.loads(raw)
            except json.JSONDecodeError as error:
                raise ValueError(f'Incomplete or malformed JSON: {path}; rerun after shard completion') from error
            nphys = run['matrix_manifest']['physical_qubits']
            check(nphys == 5238 and run['p_total'] == float(probability) and run['seed'] == int(seed)
                  and run['alpha'] == 100 and run['max_iters'] == 100, f'{path}: channel, length, seed or budget')
            if protocol == 'original':
                check(run['z_prior'] == 'conditioned' and run['x_fallback'] == 'none', f'{path}: protocol')
            else:
                check(run['orientation'] == 'global_hadamard' and run['parity_prior'] == (protocol == 'swapped'),
                      f'{path}: orientation or parity prior')
            hashes = run['matrix_manifest']['source_sha256']
            check(signatures.setdefault(int(matrix), hashes) == hashes, f'{path}: different outer matrix')
            records = run['records']
            check(run['frames'] > 0 and len(records) == run['frames'], f'{path}: frame count')
            weights = []
            for i, record in enumerate(records):
                frame_id = (int(seed), i)
                check(record['frame'] == i and frame_id not in frame_set, f'{path}: duplicate/nonsequential ID')
                weight = record['lift_pauli_weight' if protocol == 'original' else 'pauli_residual_weight']
                check(0 <= weight <= nphys, f'{path}: residual physical weight')
                x_ok, z_ok = record['x_syndrome_ok'], record['z_syndrome_ok']
                logical_ok = record['logical_ok']
                xlogical, zlogical = record['x_logical_ok'], record['z_logical_ok']
                check(logical_ok == (xlogical and zlogical) and (not logical_ok or (x_ok and z_ok)),
                      f'{path}: logical/syndrome flags')
                draws = tuple(record[k] for k in ('physical_x_only','physical_y','physical_z_only'))
                frame_set[frame_id] = (weight, draws, (not logical_ok, not xlogical, not zlogical,
                                           not x_ok, not z_ok, x_ok and z_ok and not logical_ok, weight > 0))
                weights.append(weight)
            summary = run['likelihood_lift'] if protocol == 'original' else run
            check(np.isclose(summary['qber_pauli'], np.mean(weights)/nphys, atol=1e-15, rtol=0)
                  and summary['wer_physical'] == np.count_nonzero(weights)/len(weights), f'{path}: summary mismatch')
            metadata.update(records=len(records), sha256=sha256(raw).hexdigest(), status='completed',
                            source_hashes=json.dumps(run['source_sha256'], sort_keys=True))
    if 1 in signatures and 2 in signatures:
        check(signatures[1] != signatures[2], 'Two matrix labels identify identical outer matrices')
    return frames, shard_log


def qber_interval(weights, nphys, seed):
    n = len(weights)
    if not np.any(weights):
        # QBER <= probability of a nonzero physical residual. Upper bound remains valid
        # for arbitrary errors correlated within a frame, even when bootstrap collapses.
        return 0., -np.expm1(np.log(.025)/n), 'bounded by exact physical WER interval'
    if n < 2 or np.all(weights == weights[0]):
        # Distribution-free two-sided Hoeffding interval for bounded frame means.
        mean = weights.mean()/nphys
        width = np.sqrt(np.log(40)/(2*n))
        return max(0, mean-width), min(1, mean+width), 'bounded-frame Hoeffding'
    rng = np.random.default_rng(seed)
    samples = np.empty(10000)
    chunk = max(1, min(200, 500000//n))
    for start in range(0, len(samples), chunk):
        batch = min(chunk, len(samples)-start)
        samples[start:start+batch] = weights[rng.integers(0,n,(batch,n))].mean(axis=1)/nphys
    lo, hi = np.quantile(samples,[.025,.975])
    return lo, hi, 'frame bootstrap percentile'


def summarize(frames, shard_log, plan):
    output = []
    for ordinal, key in enumerate(sorted(set(frames) | set(plan))):
        proto, matrix, p = key
        records = frames.get(key,{})
        n = len(records)
        expected = plan.get(key)
        logs = [s for s in shard_log if (s['protocol'],int(s['matrix']),float(s['p'])) == key]
        pending = sum(s['status'] != 'completed' for s in logs)
        row = dict(protocol=proto, matrix=matrix,p=p, observed_frames=n,target_frames=expected,
                   completed_shards=len(logs)-pending, pending_shards=pending,
                   status=('pending' if not n else 'partial' if pending or (expected and n < expected)
                           else 'over target' if expected and n > expected else 'target reached' if expected
                           else 'observed; target unspecified'))
        for metric in CI_COUNT_KEYS:
            for suffix in ('events','rate','lo','hi'):
                row[metric+'_'+suffix] = None
        row.update(residual_qubits=None, qubit_observations=None, qber=None,qber_lo=None,qber_hi=None,
                   qber_interval_method=None)
        if n:
            vals = list(records.values())
            counts = np.array([v[2] for v in vals]).sum(axis=0)
            for metric, k in zip(CI_COUNT_KEYS,counts):
                ci = binomtest(int(k),n).proportion_ci(confidence_level=.95)
                row.update({metric+'_events':int(k),metric+'_rate':k/n,
                            metric+'_lo':ci.low,metric+'_hi':ci.high})
            weights = np.array([v[0] for v in vals])
            lo,hi,method = qber_interval(weights,5238,120042+ordinal)
            row.update(residual_qubits=int(weights.sum()),qubit_observations=n*5238,
                       qber=weights.mean()/5238,qber_lo=lo,qber_hi=hi,qber_interval_method=method)
        output.append(row)
    return output


def paired_rows(frames):
    rows = []
    for key_a, key_b in combinations(sorted(frames),2):
        ap,am,p = key_a
        bp,bm,pb = key_b
        if p != pb or (ap != bp and am != bm):
            continue
        a,b = frames[key_a],frames[key_b]
        ids = a.keys() & b.keys()
        rescued = regressed = both_failed = 0
        for frame in ids:
            wa,draw_a,flags_a = a[frame]
            wb,draw_b,flags_b = b[frame]
            check(draw_a == draw_b, f'Noise composition mismatch at {key_a}/{key_b}, frame {frame}')
            fa,fb = flags_a[0], flags_b[0]
            rescued += int(fa and not fb)
            regressed += int(fb and not fa)
            both_failed += int(fa and fb)
        n = len(ids)
        discord = rescued+regressed
        rows.append(dict(protocol_a=ap,matrix_a=am,protocol_b=bp,matrix_b=bm,p=p,paired_frames=n,
                         unpaired_a=len(a)-n,unpaired_b=len(b)-n,
                         a_fail_b_ok=rescued if n else None, a_ok_b_fail=regressed if n else None,
                         both_fail=both_failed if n else None, both_ok=n-both_failed-discord if n else None,
                         exact_mcnemar_p=(binomtest(rescued,discord).pvalue if discord else 1.) if n else None))
    return rows


def plot(run_dir, rows, reference, bounds):
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.1), layout='constrained', sharex=True)
    colors = {'original': '#a62f24', 'swapped': '#126a9d', 'swapped_plain': '#44723b'}
    markers = {1: 'o', 2: 's'}
    metrics = ('qber', 'physical_wer', 'logical')
    for axis, metric, label in zip(axes, metrics, ('Physical QBER', 'Physical WER', 'Logical FER')):
        axis.set(xlabel=r'Pauli error probability $p$', ylabel=label, yscale='log', ylim=(1e-7, 1.6))
        axis.grid(alpha=.2)
        axis.set_axisbelow(True)
        if metric != 'logical':
            column = 'qber_digitized' if metric == 'qber' else 'wer_digitized'
            for length, style in [('19014','--'), ('38028',':')]:
                points = [r for r in reference if r['N_legend'] == length]
                axis.plot([float(r['p']) for r in points], [float(r[column]) for r in points],
                          style, color='.45', lw=1.1, label=f'Fuentes, N={length}')
            limit = 'physical_qber_lower_entropy' if metric == 'qber' else 'physical_wer_lower_holder'
            axis.plot([float(r['p']) for r in bounds],
                      [max(1e-8,float(r[limit])) for r in bounds],
                      color=colors['original'], lw=1, alpha=.5, label='Original orientation lower bound')
        for protocol in PROTOCOLS:
            for m, marker in markers.items():
                series = sorted((r for r in rows if r['protocol'] == protocol
                                 and r['matrix'] == m and r['observed_frames']), key=lambda r: r['p'])
                if not series:
                    continue
                value_key = metric if metric == 'qber' else metric+'_rate'
                xs = np.asarray([r['p'] for r in series])
                ys = np.asarray([r[value_key] for r in series])
                lows = np.asarray([r[metric+'_lo'] for r in series])
                highs = np.asarray([r[metric+'_hi'] for r in series])
                positive = ys > 0
                color = colors[protocol]
                if positive.any():
                    y = ys[positive]
                    axis.errorbar(xs[positive], y, [y-np.minimum(y,np.maximum(lows[positive],1e-8)),
                                                   np.maximum(0,highs[positive]-y)],
                                  color=color, marker=marker, mfc=color if m == 1 else 'white',
                                  ls='', capsize=2, ms=4, label=f'{PROTOCOLS[protocol]}, M{m}')
                for x, h in zip(xs[~positive], highs[~positive]):
                    axis.plot(x, h, 'v', color=color, mfc='white', ms=4)
    for axis in axes:
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(handles, labels, fontsize=6, loc='lower right')
    fig.savefig(run_dir/'pra_alpha100_comparison.pdf')
    fig.savefig(run_dir/'pra_alpha100_comparison.png', dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--plan', type=Path)
    args = parser.parse_args()
    plan_path = args.plan or args.run_dir/'plan.csv'
    plan = {}
    if plan_path.exists():
        plan = {(r['protocol'],int(r['matrix']),float(r['p'])):int(r['frames']) for r in load_csv(plan_path)}
    frames, shard_log = read_shards(args.run_dir)
    rows = summarize(frames, shard_log, plan)
    if not rows:
        parser.error('No completed or planned batches')
    export(args.run_dir/'summary.csv', rows)
    if shard_log:
        export(args.run_dir/'shards.csv', shard_log)
    pairs = paired_rows(frames)
    if pairs:
        export(args.run_dir/'paired.csv', pairs)
    plot(args.run_dir, rows, load_csv(HERE/'pra_fig10_reference.csv'),
         load_csv(HERE/'pra_physical_bounds.csv'))
    print(f'Wrote summary and comparison to {args.run_dir}')


if __name__ == '__main__':
    main()
