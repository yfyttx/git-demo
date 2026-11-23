"""Same-noise comparison of five-position Hadamard masks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from math import expm1, log
import os
from pathlib import Path
import subprocess
import sys

from scipy.stats import binomtest

from run_pra_batches import MATRICES


ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results/pra_equal_count_hadamard_20261119'
TAIL = 9
COMMON = [1, 3, 5, 7]
PLANS = {
    'main': (0.1, range(20261800, 20261804), (TAIL, 8)),
    'stress': (0.13, range(20261810, 20261812), (TAIL, 8, 2)),
}


def run_batch(job):
    phase, matrix, p, seed, position = job
    path = OUT / f'{phase}_m{matrix}_pos{position}_s{seed}.json'
    mask = sorted(COMMON + [position])
    matrix_dir = MATRICES[matrix - 1]
    if not path.exists():
        command = [sys.executable, str(ROOT / 'run_companion_local_hadamard.py'),
                   '--frames','100','--max-iters','100','--alpha','100','--p-total',str(p),
                   '--seed',str(seed),'--matrix-dir',str(matrix_dir),'--output',str(path),
                   '--hadamards',*map(str,mask)]
        with path.with_suffix('.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'})
        print(f'Finished {path.name}', flush=True)
    data = json.loads(path.read_text())
    if (data['hadamards'] != mask or data['p_total'] != p or data['alpha'] != 100
            or data['frames'] != 100 or data['max_iters'] != 100 or data['seed'] != seed
            or data['matrix_manifest'] != json.loads((matrix_dir / 'manifest.json').read_text())
            or len(data['records']) != data['frames']):
        raise ValueError(f'Unexpected batch settings in {path}')
    return (matrix, seed, position), data


def arm_summary(records):
    n = len(records)
    failures = sum(not r['logical_ok'] for r in records)
    wer_events = sum(r['physical_residual_weight'] > 0 for r in records)
    ci = binomtest(failures, n).proportion_ci()
    return dict(logical_failures=failures, logical_fer=failures/n, fer_exact_95=[ci.low,ci.high],
                wer_events=wer_events, wer_physical=wer_events/n,
                qber_pauli=sum(r['physical_residual_weight'] for r in records)/(5238*n),
                staged_failures=sum(not r['staged_syndrome_ok'] for r in records),
                projected_failures=sum(not r['projected_syndrome_ok'] for r in records),
                undetected=sum(r['undetected'] for r in records),
                zero_upper95=-expm1(log(.05)/n) if not failures else None)


def compare(results, seeds, positions):
    comparisons = []
    for matrix in (1,2):
        for position in positions:
            if position == TAIL:
                continue
            head, tail = [], []
            for seed in seeds:
                a, b = results[matrix,seed,position], results[matrix,seed,TAIL]
                for key in ('p_total','alpha','seed','frames','noise_sha256','matrix_manifest',
                            'max_iters','source_sha256'):
                    if a[key] != b[key]:
                        raise ValueError(f'Pair mismatch {key}, matrix {matrix}, seed {seed}')
                for x,y in zip(a['records'],b['records'],strict=True):
                    if x['frame'] != y['frame']:
                        raise ValueError('Paired frame indices differ')
                    head.append(x)
                    tail.append(y)
            rescued = sum(not a['logical_ok'] and b['logical_ok'] for a,b in zip(head,tail,strict=True))
            regressed = sum(a['logical_ok'] and not b['logical_ok'] for a,b in zip(head,tail,strict=True))
            test = binomtest(rescued,rescued+regressed).pvalue if rescued+regressed else 1.0
            comparisons.append(dict(matrix=matrix, head_position=position, frames=len(head),
                                    head=arm_summary(head), tail=arm_summary(tail),
                                    rescued=rescued, regressed=regressed, mcnemar_exact_p=test))
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=PLANS)
    parser.add_argument('--workers',type=int,default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    p,seeds,positions = PLANS[args.phase]
    OUT.mkdir(parents=True,exist_ok=True)
    jobs = [(args.phase,m,p,seed,position) for seed in seeds
            for position in reversed(positions) for m in (1,2)]
    with ThreadPoolExecutor(args.workers) as pool:
        results = dict(pool.map(run_batch,jobs))
    report = dict(p=p,alpha=100,seeds=list(seeds),comparisons=compare(results,seeds,positions))
    (OUT/f'{args.phase}_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
