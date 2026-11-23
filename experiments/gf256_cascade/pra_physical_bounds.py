"""Finite-block physical-error converses for the original GF(256) [9,8,2] cascade.

Reproduce the reference table with python3 experiments/gf256_cascade/pra_physical_bounds.py.
All rows use an iid Pauli channel and noiseless CSS syndromes.
"""
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp, xlogy
from scipy.stats import binom

BLOCKS = 582
M = 1552
ALPHA = 100
OUTPUT = Path(__file__).with_suffix('.csv')


def binary_entropy(x):
    return -(xlogy(x, x) + xlogy(1-x, 1-x)) / np.log(2)


def projected_source(p):
    """Group X masks by weight; permuting 9 coordinates permutes complement pairs."""
    x_marginal = 2*p/(ALPHA+2)
    z_if_no_x = ALPHA*p/(ALPHA+2)/(1-x_marginal)
    mask_logp = binom.logpmf(np.arange(10), 9, x_marginal)
    symbols = np.arange(512)
    projection_logp = []
    expected_inner_error = 0.
    for k in range(10):
        zlogp = np.zeros(512)
        for j in range(9):
            bit_p = .5 if j < k else z_if_no_x
            zlogp += np.where((symbols >> j) & 1, np.log(bit_p), np.log1p(-bit_p))
        member0, member1 = zlogp[:256], zlogp[::-1][:256]
        projection_logp.append(np.logaddexp(member0, member1))
        expected_inner_error += np.exp(mask_logp[k] + logsumexp(np.minimum(member0, member1)))
    return mask_logp, np.asarray(projection_logp), expected_inner_error


def head_success_upper(n, bit_p):
    """Mass of the 2^M most likely iid Bernoulli words, including fractional last shell."""
    remaining = 1 << M
    if remaining >= (1 << n):
        return 1.
    weight, shell = 0, 1
    while remaining >= shell:
        remaining -= shell
        shell = shell * (n-weight) // (weight+1)
        weight += 1
    return binom.cdf(weight-1, n, bit_p) + (remaining/shell)*binom.pmf(weight, n, bit_p)


def counting_upper(x_marginal, z_if_no_x):
    heads = BLOCKS * 8
    # Reveal all tails and all Z values at X=1. Bound the numerical tail by success 1.
    cutoff = int(binom.isf(1e-15, heads, x_marginal))
    tail = binom.sf(cutoff, heads, x_marginal)
    return min(1., tail + math.fsum(
        binom.pmf(k, heads, x_marginal)*head_success_upper(heads-k, z_if_no_x)
        for k in range(cutoff+1)))


def summarize(p):
    x_marginal = 2*p/(ALPHA+2)
    z_if_no_x = ALPHA*p/(ALPHA+2)/(1-x_marginal)
    mask_logp, source_logp, inner_error = projected_source(p)
    entropy_u = -np.sum(np.exp(mask_logp[:,None]+source_logp)*source_logp) / np.log(2)

    def log_holder(s):
        log_norm = logsumexp(s*source_logp, axis=1)/s
        return M*(1-1/s)*np.log(2) + BLOCKS*logsumexp(mask_logp+log_norm)

    # Any fixed s>1 yields a valid upper bound; missed numerical optimum can only weaken it.
    s = minimize_scalar(log_holder, bounds=(1.00000001, 3.), method='bounded').x
    log_success = min(0., log_holder(s))
    head_success = counting_upper(x_marginal, z_if_no_x)
    entropy_z = x_marginal + (1-x_marginal)*binary_entropy(z_if_no_x)
    remaining = max(0., entropy_z-M/(9*BLOCKS))
    distortion = brentq(lambda q: binary_entropy(q)-remaining, 0., .5)
    return dict(p=p, alpha=ALPHA, n_blocks=BLOCKS, z_syndrome_bits=M,
                entropy_projected_bits_per_block=entropy_u, entropy_budget_per_block=M/BLOCKS,
                holder_order=s, projected_success_upper_holder=np.exp(log_success),
                projected_success_upper_counting=head_success, physical_wer_lower_holder=-np.expm1(log_success),
                physical_wer_lower_representative=-np.expm1(BLOCKS*np.log1p(-inner_error)),
                physical_qber_lower_representative=inner_error,
                physical_qber_lower_entropy=distortion)


if __name__ == '__main__':
    points = sorted({i/1000 for i in range(20,101,5)} | {.0705, .0675, .0685, .07125})
    rows = [summarize(p) for p in points]
    with OUTPUT.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f'Wrote {len(rows)} points to {OUTPUT}')
