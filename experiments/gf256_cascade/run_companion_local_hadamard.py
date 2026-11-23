#!/usr/bin/env python3
"""Decode a companion cascade conjugated by a repeated local Hadamard mask."""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import sparse

from construct_companion_cascade import ROOT, mul, read_outer
from run_companion_staged_xz import decode_z, likelihood_lift
from run_companion_swapped import (BLOCK_WEIGHTS, HEAD_BITS, HEAD_PARITY,
                                   decode_parity_side, trace)
from run_companion_z_bp import OuterBP, member, row_basis


def parity_prior(probabilities, dual):
    head = np.prod(np.where(HEAD_BITS, probabilities[:8], 1-probabilities[:8]), axis=1)
    tail = np.arange(2)[:, None] ^ HEAD_PARITY[dual]
    table = head[dual] * np.where(tail == 1, probabilities[8], 1-probabilities[8])
    return table / table.sum(axis=1, keepdims=True)


def project_prior(bit_prob):
    bit_prob = np.asarray(bit_prob).reshape(-1, 9)
    head = np.prod(np.where(HEAD_BITS, bit_prob[:, None, :8], 1-bit_prob[:, None, :8]), axis=2)
    tail = bit_prob[:, 8, None]
    return (1-tail)*head + tail*head[:, ::-1]


def simulate(args):
    x_check = sparse.load_npz(args.matrix_dir/'HX_binary.npz').astype(np.int64)
    z_check = sparse.load_npz(args.matrix_dir/'HZ_binary.npz').astype(np.int64)
    manifest = json.loads((args.matrix_dir/'manifest.json').read_text())
    n, gamma_rows = read_outer(args.matrix_dir/'HX_TIT_Outer_For_HZ_GF256.txt')
    n_delta, delta_rows = read_outer(args.matrix_dir/'HZ_TIT_Outer_For_HX_GF256.txt')
    if n_delta != n or x_check.shape != (8*len(delta_rows),9*n) or z_check.shape != (8*len(gamma_rows)+n,9*n):
        raise ValueError('Unexpected companion matrix dimensions')
    if np.any((x_check @ z_check.T).data % 2):
        raise ValueError('Matrix does not satisfy CSS commutation')
    mask = np.array([i in args.hadamards for i in range(1,10)])
    opposite = np.tile(mask,n).astype(np.int64)
    # Symplectic check of the physical mixed-Pauli generators.
    gx_x, gx_z = x_check.multiply(1-opposite), x_check.multiply(opposite)
    gz_x, gz_z = z_check.multiply(opposite), z_check.multiply(1-opposite)
    for first_x,first_z,second_x,second_z in [
        (gx_x,gx_z,gx_x,gx_z), (gz_x,gz_z,gz_x,gz_z), (gx_x,gx_z,gz_x,gz_z)
    ]:
        if np.any((first_x @ second_z.T + first_z @ second_x.T).data % 2):
            raise ValueError('Physical Pauli generators do not commute')
    x_basis, z_basis = row_basis(x_check), row_basis(z_check)
    gamma, delta = OuterBP(gamma_rows,n), OuterBP(delta_rows,n)
    dual = np.array([sum(trace(mul(v,1<<i))<<i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)
    px = py = args.p_total/(args.alpha+2)
    pz = args.alpha*px
    prob_staged = np.where(mask,pz+py,px+py)
    prob_projected = np.where(mask,px+py,pz+py)
    prob_z_only = np.where(mask,px,pz)
    prior_staged = parity_prior(prob_staged,dual)
    marginal = project_prior(prob_projected)[0]
    staged_patterns = (np.arange(512)[:, None] >> np.arange(9)) & 1
    conditional = np.where(staged_patterns, py/prob_staged, prob_z_only/(1-prob_staged))
    condition_prior = project_prior(conditional)
    rng = np.random.default_rng(args.seed)
    noise_hash = hashlib.sha256()
    records = []
    started = perf_counter()
    for frame in range(args.frames):
        draw = rng.random((n,9))
        physical_x = (draw < px) | ((draw >= pz+py) & (draw < args.p_total))
        physical_z = (draw >= px) & (draw < args.p_total)
        noise_hash.update((physical_x + 2*physical_z).astype(np.uint8).tobytes())
        x = np.where(mask, physical_z, physical_x).ravel().astype(np.int64)
        z = np.where(mask, physical_x, physical_z).ravel().astype(np.int64)
        sx = (z_check @ x) % 2
        sz = (x_check @ z) % 2
        if not (np.array_equal(sx,(gz_x@physical_z.ravel()+gz_z@physical_x.ravel())%2)
                and np.array_equal(sz,(gx_x@physical_z.ravel()+gx_z@physical_x.ravel())%2)):
            raise ValueError('Physical and virtual syndromes disagree')
        begin = perf_counter()
        xhat, xi, xv = decode_parity_side(sx, gamma, prior_staged, dual, inverse_dual, n, args.max_iters)
        if xv:
            lookup = xhat.reshape(n,9) @ BLOCK_WEIGHTS
            q, prior = conditional[lookup], condition_prior[lookup]
        else:
            q, prior = prob_projected, marginal
        zhat, zi, zv = decode_z(sz, delta, prior, n, args.max_iters)
        zhat = likelihood_lift(zhat,q)
        seconds = perf_counter() - begin
        if xv != np.array_equal((z_check@xhat)%2,sx) or zv != np.array_equal((x_check@zhat)%2,sz):
            raise ValueError('Invalid syndrome status')
        xr, zr = x ^ xhat, z ^ zhat
        xo, zo = bool(xv and member(xr,x_basis)), bool(zv and member(zr,z_basis))
        union = xr | zr  # Coordinate-wise X/Z swaps leave the union unchanged.
        records.append(dict(frame=frame, staged_syndrome_ok=bool(xv), projected_syndrome_ok=bool(zv),
                            logical_ok=xo and zo, undetected=bool(xv and zv and not (xo and zo)),
                            physical_residual_weight=int(union.sum()), seconds=seconds,
                            staged_iterations=xi, projected_iterations=zi))
    sources = [Path(__file__), ROOT/'run_companion_swapped.py', ROOT/'run_companion_z_bp.py',
               ROOT/'run_companion_staged_xz.py', ROOT/'construct_companion_cascade.py',
               ROOT/'run_identity_outer_comparison.py']
    nfail = sum(not r['logical_ok'] for r in records)
    return dict(frames=args.frames, seed=args.seed, p_total=args.p_total, alpha=args.alpha,
                hadamards=args.hadamards, max_iters=args.max_iters, matrix_manifest=manifest,
                noise_sha256=noise_hash.hexdigest(), logical_failures=nfail, logical_fer=nfail/args.frames,
                wer_physical=sum(r['physical_residual_weight']>0 for r in records)/args.frames,
                qber_pauli=sum(r['physical_residual_weight'] for r in records)/(9*n*args.frames),
                staged_syndrome_failures=sum(not r['staged_syndrome_ok'] for r in records),
                projected_syndrome_failures=sum(not r['projected_syndrome_ok'] for r in records),
                undetected_logical_failures=sum(r['undetected'] for r in records),
                source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
                elapsed_seconds=perf_counter()-started, records=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames',type=int,required=True)
    parser.add_argument('--p-total',type=float,default=.07)
    parser.add_argument('--alpha',type=float,default=100.)
    parser.add_argument('--seed',type=int,required=True)
    parser.add_argument('--max-iters',type=int,default=100)
    parser.add_argument('--matrix-dir',type=Path,default=ROOT/'results/companion_binary_n9_k8_matrix20261019')
    parser.add_argument('--hadamards',type=int,nargs='+',default=[1,3,5,7,9])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if (not np.isfinite(args.alpha) or args.alpha <= 0 or not 0 < args.p_total < 1
            or min(args.frames, args.max_iters)<1 or not set(args.hadamards)<=set(range(1,10))
            or len(args.hadamards)!=len(set(args.hadamards))):
        parser.error('Invalid channel, budget or Hadamard mask')
    result = simulate(args)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    tmp = args.output.with_suffix('.partial')
    tmp.write_text(json.dumps(result,indent=2)+'\n')
    tmp.replace(args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','matrix_manifest','source_sha256')},
                     indent=2))


if __name__ == '__main__':
    main()
