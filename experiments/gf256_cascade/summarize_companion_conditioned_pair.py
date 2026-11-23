#!/usr/bin/env python3
"""Summarize matched marginal/estimated-X-conditioned companion simulations."""
import argparse
import json
from pathlib import Path
from scipy.stats import binomtest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir', type=Path)
    args = parser.parse_args()
    a = json.loads((args.run_dir/'marginal.json').read_text())
    b = json.loads((args.run_dir/'conditioned.json').read_text())
    for key in ['frames', 'p_D', 'seed', 'max_iters', 'matrix_manifest', 'source_sha256']:
        if a[key] != b[key]:
            raise ValueError(f'Unmatched {key}')
    if a['z_prior'] != 'marginal' or b['z_prior'] != 'conditioned':
        raise ValueError('Unexpected comparison modes')
    rescued = worsened = both_failed = 0
    pairs = []
    for x, y in zip(a['records'], b['records']):
        for key in ['frame', 'x_syndrome_ok', 'x_logical_ok', 'exact_x', 'x_iterations']:
            if x[key] != y[key]:
                raise ValueError(f'Unmatched X stage at {x["frame"]}: {key}')
        rescued += not x['logical_ok'] and y['logical_ok']
        worsened += x['logical_ok'] and not y['logical_ok']
        both_failed += not x['logical_ok'] and not y['logical_ok']
        pairs.append(dict(frame=x['frame'], marginal_logical_ok=x['logical_ok'],
                          conditioned_logical_ok=y['logical_ok']))
    fields = ['logical_failures', 'logical_fer', 'x_syndrome_failures', 'z_syndrome_failures',
              'x_logical_failures', 'z_logical_failures', 'undetected_logical_failures',
              'exact_x_frames', 'conditioned_frames', 'elapsed_seconds']
    summary = dict(frames=a['frames'], p_D=a['p_D'], f_m=a['f_m'], seed=a['seed'],
                   max_iters=a['max_iters'],
                   marginal={k:a[k] for k in fields}, conditioned={k:b[k] for k in fields},
                   rescued=rescued, worsened=worsened, both_failed=both_failed,
                   mcnemar_exact_two_sided_p=float(binomtest(rescued, rescued+worsened).pvalue)
                       if rescued+worsened else 1.,
                   timing_note='Wall time; the two runs were concurrent. Not a clean CPU-cost comparison.',
                   scope='Complete X/Z logical FER with binary stabilizer membership evaluation.',
                   records=pairs)
    (args.run_dir/'paired_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='records'}, indent=2))

if __name__ == '__main__':
    main()
