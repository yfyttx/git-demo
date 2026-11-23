#!/usr/bin/env python3
"""Offline audit of staged-X failures and exact MAP vs binary BP on a tiny GF(8) cascade.

True errors are reconstructed only for post-decoding diagnosis. They are never
passed to either decoder. The toy comparison measures classical block-error
probability, not quantum logical FER.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from construct_companion_cascade import ROOT as EXP_ROOT, read_outer, mul
from run_companion_staged_xz import decode_x
from run_companion_z_bp import OuterBP
from run_companion_staged_xz import decode_binary_bp_fallback


def trace(a):
    x = int(a)
    out = 0
    for _ in range(8):
        out ^= x & 1
        x >>= 1
    return out


def audit_failures():
    matrix_dir = EXP_ROOT / 'results/companion_binary_n9_k8_matrix20261019'
    n, _ = read_outer(matrix_dir / 'HZ_TIT_Outer_For_HX_GF256.txt')
    _, gamma = read_outer(matrix_dir / 'HX_TIT_Outer_For_HZ_GF256.txt')
    hz = sparse.load_npz(matrix_dir / 'HZ_binary.npz').astype(np.int64).tocsr()
    dual = np.array([sum(trace(mul(v, 1 << i)) << i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)
    bp = OuterBP(gamma, n)
    p_d = .07125
    f = 2*p_d/3
    w = np.array([v.bit_count() for v in range(256)])
    px = f**w * (1-f)**(8-w)
    batches = [('companion_xfallback_pair500_20261118', 20261125),
               ('companion_xfallback_pair500b_20261118', 20261126)]
    out_frames = []
    for folder, seed in batches:
        summary = json.load(open(EXP_ROOT / 'results' / folder / 'no_fallback.json'))
        fail_ids = [r['frame'] for r in summary['records'] if not r['x_syndrome_ok']]
        rng = np.random.default_rng(seed)
        for frame in range(max(fail_ids)+1 if fail_ids else 0):
            draw = rng.random((n, 9))
            if frame not in fail_ids:
                continue
            x = ((draw < p_d/3) | ((draw >= 2*p_d/3) & (draw < p_d))).astype(np.int64)
            sx = np.asarray(hz @ x.ravel()).ravel() % 2
            upper_size = 8*len(gamma)
            upper = sx[:upper_size].reshape(-1, 8) @ (1 << np.arange(8))
            target_outer_syndrome = inverse_dual[upper]
            estimate, iters, valid = decode_x(sx, bp, px, inverse_dual, dual, n, 100)
            estimate = estimate.reshape(n, 9)
            estimated_outer = inverse_dual[estimate[:, :8] @ (1 << np.arange(8))]
            residual = bp.syndrome(estimated_outer) ^ target_outer_syndrome
            out_frames.append(dict(batch=folder, seed=seed, frame=frame,
                true_outer_error_positions=np.flatnonzero(
                    inverse_dual[x[:, :8] @ (1 << np.arange(8))] != 0).tolist(),
                wrong_outer_symbols=int(np.count_nonzero(
                    inverse_dual[x[:, :8] @ (1 << np.arange(8))] != estimated_outer)),
                wrong_physical_head_bits=int(np.count_nonzero(x[:, :8] != estimate[:, :8])),
                violated_outer_checks=int(np.count_nonzero(residual)),
                violated_outer_check_positions=np.flatnonzero(residual).tolist(),
                iterations=int(iters), valid=bool(valid)))
    n_wrong = [r['wrong_outer_symbols'] for r in out_frames]
    n_checks = [r['violated_outer_checks'] for r in out_frames]
    return dict(frame_count=len(out_frames), p_D=p_d, f_m=f,
        wrong_outer_symbol_quantiles={str(q): float(np.quantile(n_wrong,q)) for q in (0,.25,.5,.75,1)},
        violated_outer_check_quantiles={str(q): float(np.quantile(n_checks,q)) for q in (0,.25,.5,.75,1)},
        frames_with_any_wrong_outer_symbol=int(sum(v>0 for v in n_wrong)),
        frames_with_any_unsatisfied_outer_check=int(sum(v>0 for v in n_checks)),
        frame_records=out_frames)


def gf8_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 0b1011  # x^3 + x + 1
    return out


def gf8_matrix(a):
    # Binary matrix for multiplication by a in polynomial basis.
    return np.array([[(gf8_mul(a, 1 << j) >> i) & 1 for j in range(3)]
                     for i in range(3)], dtype=np.int64)


def gf2_rank(a):
    rows = [sum(int(bit) << j for j, bit in enumerate(row) if bit)
            for row in np.asarray(a, dtype=np.int64)]
    rank = 0
    for col in range(a.shape[1]-1, -1, -1):
        pivot = next((i for i in range(rank, len(rows)) if (rows[i] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def tiny_map_vs_bp():
    # Inner [4,3] code: G2=[I3|1], H2=[1 1 1 1].
    g2 = np.hstack((np.eye(3, dtype=np.int64), np.ones((3,1), dtype=np.int64)))
    h2 = np.ones((1,4), dtype=np.int64)
    if np.any(h2 @ g2.T % 2):
        raise AssertionError('inner H2 and G2 are not orthogonal')
    # One GF(8) outer check u0 + alpha*u1, alpha=x (binary value 2).
    outer = np.hstack((g2, gf8_matrix(2) @ g2 % 2))
    lower = np.zeros((2,8), dtype=np.int64)
    lower[0,:4] = h2
    lower[1,4:] = h2
    h = np.vstack((outer, lower)) % 2
    # Exhaustively assign each error to its measured syndrome.
    p = .0475
    errors = np.arange(256, dtype=np.uint16)
    bits = ((errors[:,None] >> np.arange(8)) & 1).astype(np.int64)
    syndromes = (bits @ h.T) % 2
    weights = bits.sum(axis=1)
    probabilities = p**weights * (1-p)**(8-weights)
    unique = {tuple(x) for x in syndromes.tolist()}
    total_prob = float(probabilities.sum())
    map_fail = bp_fail = syndrome_fail_mass = 0.0
    syndrome_count = 0
    max_map_posterior = 0.0
    for s_tuple in unique:
        s = np.array(s_tuple, dtype=np.int64)
        mask = np.all(syndromes == s, axis=1)
        mass = float(probabilities[mask].sum())
        idx = np.flatnonzero(mask)
        best = idx[np.argmax(probabilities[idx])]
        map_fail += mass - float(probabilities[best])
        cand, _, ok = decode_binary_bp_fallback(sparse.csr_matrix(h), s, p, 100)
        syndrome_count += 1
        if not ok:
            syndrome_fail_mass += mass
            bp_fail += mass
        else:
            # A unique exact error is the block-success criterion for this classical toy.
            got = int(sum(int(cand[i]) << i for i in range(8)))
            bp_match = np.flatnonzero(np.all(bits[idx] == cand[None, :], axis=1))
            if len(bp_match):
                bp_fail += mass - float(probabilities[idx[bp_match[0]]])
            else:
                bp_fail += mass
        max_map_posterior = max(max_map_posterior, float(probabilities[best]/mass))
    return dict(field='GF(8), polynomial x^3+x+1', alpha=2,
        inner_G2=g2.tolist(), inner_H2=h2.tolist(), outer_check='u0 + alpha*u1',
        binary_check_matrix=h.tolist(), n_bits=8, check_rank=gf2_rank(h),
        bsc_p=p, syndromes=syndrome_count, map_block_error_probability=map_fail/total_prob,
        bp_block_error_probability=bp_fail/total_prob,
        bp_syndrome_failure_probability=syndrome_fail_mass/total_prob,
        map_is_lower_bound=True,
        note='Classical exact-error toy; no CSS stabilizer quotient and not directly extrapolatable to the GF(256) quantum code.')


def main():
    out = dict(actual_x_failure_audit=audit_failures(), tiny_gf8_map_vs_bp=tiny_map_vs_bp())
    dest = EXP_ROOT / 'results/companion_xfailure_map_audit_20261118.json'
    dest.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({
        'audit_frames': out['actual_x_failure_audit']['frame_count'],
        'actual_summary': {k:v for k,v in out['actual_x_failure_audit'].items() if k != 'frames'},
        'toy': out['tiny_gf8_map_vs_bp'], 'output': str(dest)}, indent=2))


if __name__ == '__main__':
    main()
