"""Best possible exact physical recovery with true X and true G2*Z in every block.

Uses the nine-qubit SPC construction, independent Pauli noise, and the paper's
pX=pY, pZ=alpha*pX convention. Only the physical residual is scored.
"""
import argparse
from math import expm1, log1p
import numpy as np


def oracle_physical_rates(p, alpha, blocks):
    bit = p / (2 + alpha)
    pauli = np.array([[1-p, alpha*bit], [bit, bit]])
    bits = (np.arange(512)[:, None] >> np.arange(9)) & 1
    zero_tail = bits[:256]
    joint0 = pauli[bits[:, None, :], zero_tail[None, :, :]].prod(axis=2)
    joint1 = pauli[bits[:, None, :], 1-zero_tail[None, :, :]].prod(axis=2)
    if abs((joint0 + joint1).sum()-1) > 1e-10:
        raise ArithmeticError("Joint masses must sum to one")
    qubit_error = float(np.minimum(joint0, joint1).sum())
    frame_error = -expm1(blocks * log1p(-qubit_error))
    return qubit_error, frame_error


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p', type=float, nargs='+', default=[0.07, 0.0705])
    parser.add_argument('--alpha', type=float, default=100)
    parser.add_argument('--blocks', type=int, default=582)
    args = parser.parse_args()
    if args.blocks < 1 or not np.isfinite(args.alpha) or args.alpha < 0:
        parser.error('Supply a positive block count and a finite, nonnegative alpha.')
    print('p,alpha,blocks,oracle_physical_qber,oracle_physical_wer')
    for p in args.p:
        if not 0 <= p <= 1:
            parser.error('p must lie in [0,1].')
        qber, wer = oracle_physical_rates(p, args.alpha, args.blocks)
        print(f'{p},{args.alpha},{args.blocks},{qber:.12g},{wer:.12g}')
