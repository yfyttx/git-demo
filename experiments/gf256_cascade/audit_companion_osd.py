#!/usr/bin/env python3
"""Use simulation truth offline to audit which OSD candidates were reachable."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from companion_gf256_osd import SyndromeOSD
from run_companion_swapped import (BIT_WEIGHTS, HEAD_BITS, HEAD_PARITY, ROOT,
                                   gamma_prior_table, mul, read_outer, trace)
from run_companion_z_bp import OuterBP, member, row_basis


def audit_matrix(matrix):
    folder = ('companion_binary_n9_k8_matrix20261019' if matrix == 1
              else 'companion_binary_n9_k8_outer_path2')
    matrix_dir = ROOT/'results'/folder
    baseline = ROOT/'results/pra_alpha100_20261118/swapped'/f'm{matrix}_p0.07_s20261500.json'
    saved = json.loads(baseline.read_text())
    p, alpha = saved['p_total'], saved['alpha']
    px = p/(alpha+2)
    fz = (alpha+1)*px
    hx = sparse.load_npz(matrix_dir/'HZ_binary.npz').astype(np.int64)
    zbasis = row_basis(sparse.load_npz(matrix_dir/'HX_binary.npz'))
    n, rows = read_outer(matrix_dir/'HX_TIT_Outer_For_HZ_GF256.txt')
    bp = OuterBP(rows, n)
    osd = SyndromeOSD(bp)
    dual = np.array([sum(trace(mul(v, 1 << i)) << i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)
    priors = gamma_prior_table(fz, dual, True)
    rng = np.random.default_rng(saved['seed'])
    results = []
    for record in saved['records']:
        draw = rng.random((n,9))
        if record['z_syndrome_ok']:
            continue
        z = ((draw >= px) & (draw < p)).astype(np.int64)
        syn = (hx @ z.ravel()) % 2
        lower = syn[-n:]
        upper = inverse_dual[syn[:-n].reshape(-1,8) @ BIT_WEIGHTS]
        hard, _, valid, posterior = bp.decode(upper, priors[lower], saved['max_iters'],
                                              return_log_posterior=True)
        if valid:
            raise ValueError('Historical failure failed to replay')
        winner, expected_count = osd.decode(upper, hard, posterior, priors[lower], order=1)
        # Independently reconstruct the candidate set. Only the checks below
        # use the error drawn by the simulator.
        choices = np.argsort(-posterior, axis=1, kind='stable')[:, :5]
        gap = posterior[np.arange(n), choices[:,0]] - posterior[np.arange(n), choices[:,1]]
        columns = np.argsort(gap, kind='stable')
        reduced, positions = osd.reduce(upper, columns)
        pivots = columns[positions]
        free_pos = np.setdiff1d(np.arange(n), positions)
        free = columns[free_pos]
        truth = inverse_dual[z[:,:8] @ BIT_WEIGHTS]
        base = hard.copy()
        base[pivots] = reduced[:, -1] ^ np.bitwise_xor.reduce(
            osd.product[reduced[:, free_pos], hard[free]], axis=1)
        candidates = [base]
        for position, symbol in zip(free_pos[:32], free[:32]):
            for value in choices[symbol,1:]:
                trial = base.copy()
                trial[symbol] = value
                trial[pivots] ^= osd.product[reduced[:, position], hard[symbol] ^ value]
                candidates.append(trial)
        if len(candidates) != expected_count or not any(np.array_equal(c, winner) for c in candidates):
            raise ValueError('Audit list disagrees with online OSD')
        covered = 0
        possible_xor = z.ravel()
        for c in candidates:
            head = HEAD_BITS[dual[c]]
            tail = lower ^ HEAD_PARITY[dual[c]]
            physical = np.column_stack((head,tail)).ravel()
            if not np.array_equal((hx @ physical) % 2, syn):
                raise ValueError('Candidate fails complete physical syndrome')
            covered += member(possible_xor ^ physical, zbasis)
        out = dict(matrix=matrix, frame=record['frame'], candidates=len(candidates),
                   free_symbols_wrong=int(np.count_nonzero(truth[free] != hard[free])),
                   correct_z_logical_candidates=covered)
        results.append(out)
        if len(results) == 20:
            break
    if len(results) != 20:
        raise ValueError('Expected 20 historical failures per matrix')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with ProcessPoolExecutor(2) as pool:
        records = sum(pool.map(audit_matrix, (1,2)), [])
    args.output.write_text(json.dumps(records, indent=2)+'\n')
    print('Frames with correct Z logical class:', sum(r['correct_z_logical_candidates'] > 0 for r in records))


if __name__ == '__main__':
    main()
