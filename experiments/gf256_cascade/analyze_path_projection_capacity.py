"""Conditional additive-channel capacity of one unit-path projection.

This is an information diagnostic for independently decoded projections,
not a capacity calculation for the jointly decoded cascade code.
"""

import argparse
import math


def binary_entropy(probability):
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return -(probability * math.log2(probability)
             + (1.0 - probability) * math.log2(1.0 - probability))


def projection_capacity(p_depolar):
    # Depolarizing Pauli model: P(X)=P(Y)=P(Z)=p_D/3.
    px = 2.0 * p_depolar / 3.0
    z_given_x0 = (p_depolar / 3.0) / (1.0 - px)
    states = ((1.0 - px, z_given_x0), (px, 0.5))
    conditional_entropy = sum(
        prob_a * prob_b * binary_entropy(pa + pb - 2.0 * pa * pb)
        for prob_a, pa in states for prob_b, pb in states)
    return 1.0 - conditional_entropy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-depolar", type=float, nargs="*",
                        default=(0.06, 0.0675, 0.07125))
    parser.add_argument("--outer-cols", type=int, default=582)
    parser.add_argument("--outer-rows", type=int, default=194)
    args = parser.parse_args()
    if args.outer_cols <= 0 or not 0 <= args.outer_rows < args.outer_cols:
        parser.error("invalid outer matrix dimensions")
    rate_lower_bound = 1.0 - args.outer_rows / args.outer_cols
    print(f"outer_rate_lower_bound={rate_lower_bound:.9f}")
    for probability in args.p_depolar:
        if not 0.0 < probability < 1.0:
            parser.error("p_D must lie in (0,1)")
        capacity = projection_capacity(probability)
        print(f"p_D={probability:.8f} "
              f"conditioned_projection_capacity={capacity:.9f} "
              f"capacity_minus_rate_lower_bound={capacity-rate_lower_bound:+.9f}")
    low, high = 0.0, 0.5
    for _ in range(80):
        middle = (low + high) / 2.0
        if projection_capacity(middle) > rate_lower_bound:
            low = middle
        else:
            high = middle
    print(f"rate_bound_crossing_p_D={(low+high)/2.0:.9f}")


if __name__ == "__main__":
    main()
