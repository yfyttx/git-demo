#!/usr/bin/env python3
"""Exact GF(8) check of stabilizer-coset partition-function identity.

Toy CSS checks: H_X = H_Z = [1, 1, 0].  They are orthogonal in
characteristic two, with one nontrivial Z stabilizer and one logical symbol.
The script compares direct enumeration of each syndrome/logical class against
summing over the stabilizer coefficient u in GF(8).
"""
from itertools import product

# GF(8) in polynomial basis with primitive polynomial x^3 + x + 1.
Q = 8
REDUCTION = 0b1011


def gf_mul(a: int, b: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & Q:
            a ^= REDUCTION
    return out


def gf_inv(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError
    for b in range(1, Q):
        if gf_mul(a, b) == 1:
            return b
    raise AssertionError("nonzero GF(8) element has no inverse")


def physical_bits(symbol: int, inverse_scale: int) -> int:
    return gf_mul(symbol, inverse_scale)


# Distinct nonzero scale factors ensure the physical bit mapping is not simply
# the same symbol basis at every position.
scales = (1, 2, 4)
inverse_scales = tuple(gf_inv(scale) for scale in scales)
p_depolar = 0.18
p_z_given_x0 = (p_depolar / 3.0) / (1.0 - 2.0 * p_depolar / 3.0)
p_z_given_x1 = 0.5
# Fixed decoded X-symbol bit patterns; each bit conditions the corresponding
# Z probability, as in the depolarizing-channel prior used by the decoder.
x_bits = (0b001, 0b101, 0b010)


def prior(position: int, symbol: int) -> float:
    bits = physical_bits(symbol, inverse_scales[position])
    probability = 1.0
    for bit in range(3):
        p_z = p_z_given_x1 if ((x_bits[position] >> bit) & 1) else p_z_given_x0
        probability *= p_z if ((bits >> bit) & 1) else 1.0 - p_z
    return probability


# H_X H_Z^T = 1*1 + 1*1 = 0 in characteristic two.
assert (gf_mul(1, 1) ^ gf_mul(1, 1)) == 0

# Directly enumerate every error vector and group by its X-check syndrome and
# logical coordinate (the third symbol). The stabilizer shifts (u,u,0), so the
# third coordinate labels the eight logical classes within each syndrome.
direct = {}
for error in product(range(Q), repeat=3):
    syndrome = error[0] ^ error[1]
    logical = error[2]
    mass = prior(0, error[0]) * prior(1, error[1]) * prior(2, error[2])
    direct[(syndrome, logical)] = direct.get((syndrome, logical), 0.0) + mass

max_abs_error = 0.0
for syndrome in range(Q):
    for logical in range(Q):
        # r=(syndrome,0,logical) is a syndrome-compatible representative.
        # Add the Z stabilizer u*(1,1,0) and sum its class probability.
        partition = 0.0
        for u in range(Q):
            shifted = (syndrome ^ u, u, logical)
            partition += (prior(0, shifted[0]) * prior(1, shifted[1])
                          * prior(2, shifted[2]))
        max_abs_error = max(max_abs_error,
                            abs(partition - direct[(syndrome, logical)]))

        # All eight stabilizer shifts must preserve the syndrome and logical
        # class, independently of their physical probability.
        for u in range(Q):
            shifted = (syndrome ^ u, u, logical)
            assert (shifted[0] ^ shifted[1]) == syndrome
            assert shifted[2] == logical

assert max_abs_error < 1e-14, max_abs_error
print("GF(8) coset partition identity passed")
print("code: H_X=H_Z=[1,1,0], ranks=1/1, logical dimension=1")
print(f"classes checked: {Q * Q}; stabilizer terms per class: {Q}")
print(f"symbol scales: {scales}; depolarizing probability: {p_depolar}")
print(f"conditioned X bit patterns: {x_bits}")
print(f"maximum absolute mass difference: {max_abs_error:.3e}")
