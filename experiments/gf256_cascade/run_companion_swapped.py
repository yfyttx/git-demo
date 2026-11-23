#!/usr/bin/env python3
"""Global Hadamard of the [9,8] companion cascade: staged Z and projected X."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from construct_companion_cascade import ROOT, mul, read_outer
from run_companion_staged_xz import decode_z, likelihood_lift
from run_companion_z_bp import OuterBP, member, row_basis
from run_identity_outer_comparison import trace


BIT_WEIGHTS = 1 << np.arange(8)
HEAD_BITS = (np.arange(256)[:, None] >> np.arange(8)) & 1
HEAD_PARITY = HEAD_BITS.sum(axis=1) & 1
BLOCK_WEIGHTS = 1 << np.arange(9)


def gamma_prior_table(fz, dual, parity_aware):
    """Table indexed by the observed block parity and then the Gamma field symbol."""
    head_prior = np.prod(np.where(HEAD_BITS, fz, 1-fz), axis=1)
    result = np.broadcast_to(head_prior[dual], (2, 256)).copy()
    if parity_aware:
        for syndrome in (0, 1):
            tail = syndrome ^ HEAD_PARITY[dual]
            result[syndrome] *= np.where(tail, fz, 1-fz)
            result[syndrome] /= result[syndrome].sum()
    return result


def x_conditional_table(px, py, pz):
    """Tables indexed by all 512 estimated Z patterns, independent of the true errors."""
    z_patterns = (np.arange(512)[:, None] >> np.arange(9)) & 1
    z_marginal = py + pz
    bit_prob = np.where(z_patterns, py/z_marginal, px/(1-z_marginal))
    head_prob = np.prod(np.where(HEAD_BITS[None,:,:], bit_prob[:,None,:8],
                                 1-bit_prob[:,None,:8]), axis=2)
    tail_prob = bit_prob[:,8,None]
    projection = (1-tail_prob)*head_prob + tail_prob*head_prob[:,::-1]
    return bit_prob, projection


def decode_parity_side(syndrome, gamma, table, dual, inverse_dual, n, max_iters, *,
                       osd=None, order=0, layered=None):
    lower = syndrome[-n:]
    upper = inverse_dual[syndrome[:-n].reshape(-1, 8) @ BIT_WEIGHTS]
    prior = table[lower]
    if osd is None:
        symbols, iters, valid = gamma.decode(upper, prior, max_iters)
    else:
        symbols, iters, valid, posterior = gamma.decode(upper, prior, max_iters,
                                                         return_log_posterior=True)
    candidates, osd_seconds = 0, 0.0
    retry_count, retry_iterations, retry_seconds = 0, 0, 0.0
    if not valid and layered is not None:
        original_symbols = symbols
        tick = time.perf_counter()
        for reverse, damping in ((False, 0.15), (True, 0.35)):
            symbols, used, valid, posterior = layered.decode(
                upper, prior, max_iters, damping=damping, reverse=reverse)
            retry_count += 1
            retry_iterations += used
            if valid:
                break
        retry_seconds = time.perf_counter() - tick
        if not valid and osd is None:
            symbols = original_symbols
    if not valid and osd is not None:
        tick = time.perf_counter()
        symbols, candidates = osd.decode(upper, symbols, posterior, prior, order=order)
        osd_seconds = time.perf_counter() - tick
        valid = True  # SyndromeOSD verifies the outer syndrome before returning.
    head = HEAD_BITS[dual[symbols]]
    tail = lower ^ HEAD_PARITY[dual[symbols]]
    result = np.column_stack((head, tail)).ravel(), iters, valid
    if osd is None and layered is None:
        return result
    return (*result, candidates, osd_seconds, retry_count, retry_iterations, retry_seconds)


def simulate(args):
    raw_hx = sparse.load_npz(args.matrix_dir/'HX_binary.npz').astype(np.int64)
    raw_hz = sparse.load_npz(args.matrix_dir/'HZ_binary.npz').astype(np.int64)
    manifest = json.loads((args.matrix_dir/'manifest.json').read_text())
    n, delta = read_outer(args.matrix_dir/'HZ_TIT_Outer_For_HX_GF256.txt')
    n_gamma, gamma = read_outer(args.matrix_dir/'HX_TIT_Outer_For_HZ_GF256.txt')
    if n_gamma != n or raw_hz.shape != (8*len(gamma)+n,9*n) or raw_hx.shape != (8*len(delta),9*n):
        raise ValueError('Unexpected companion matrix dimensions')
    g = np.loadtxt(args.matrix_dir/'G2_binary.txt',dtype=int,ndmin=2)
    h = np.loadtxt(args.matrix_dir/'H2_binary.txt',dtype=int,ndmin=2)
    if not np.array_equal(g,np.column_stack((np.eye(8,dtype=int),np.ones(8,dtype=int)))) or not np.array_equal(h,np.ones((1,9),dtype=int)):
        raise ValueError('Expected [I8|1] and all-one H2')
    # A global Hadamard interchanges stabilizer rows while channel axes stay physical.
    hx, hz = raw_hz, raw_hx
    if np.any((hx @ hz.T).data % 2):
        raise ValueError('Swapped checks violate CSS orthogonality')
    x_basis, z_basis = row_basis(hx), row_basis(hz)
    dual = np.array([sum(trace(mul(v, 1 << i)) << i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)
    kernel = OuterBP
    if args.native:
        from companion_fast_bp import NativeBP
        kernel = NativeBP
    staged_bp, projected_bp = kernel(gamma, n), kernel(delta, n)
    osd = None
    if args.z_osd != 'none':
        from companion_gf256_osd import SyndromeOSD
        osd = SyndromeOSD(staged_bp)
    layered = None
    if args.z_layered_retry:
        from companion_layered_bp import LayeredBP
        layered = LayeredBP(staged_bp)
    px = py = args.p_total/(args.alpha+2)
    pz = args.alpha*px
    fx, fz = px+py, pz+py
    staged_table = gamma_prior_table(fz, dual, args.parity_prior)
    bit_table, projection_table = x_conditional_table(px, py, pz)
    marginal_head = np.prod(np.where(HEAD_BITS,fx,1-fx),axis=1)
    marginal_projection = (1-fx)*marginal_head + fx*marginal_head[::-1]
    rng = np.random.default_rng(args.seed)
    fingerprints = hashlib.sha256()
    records = []
    start = time.perf_counter()
    for frame in range(args.frames):
        draw = rng.random((n,9))
        x = ((draw < px) | ((draw >= fz) & (draw < args.p_total))).astype(np.int64)
        z = ((draw >= px) & (draw < args.p_total)).astype(np.int64)
        fingerprints.update((x+2*z).astype(np.uint8).tobytes())
        x, z = x.ravel(), z.ravel()
        sx = (hz @ x)%2
        sz = (hx @ z)%2
        tick = time.perf_counter()
        if osd is None and layered is None:
            zhat, zi, zv = decode_parity_side(sz, staged_bp, staged_table, dual, inverse_dual, n, args.max_iters)
            osd_candidates, osd_seconds = 0, 0.0
            retries, retry_iterations, retry_seconds = 0, 0, 0.0
        else:
            zhat, zi, zv, osd_candidates, osd_seconds, retries, retry_iterations, retry_seconds = decode_parity_side(
                sz, staged_bp, staged_table, dual, inverse_dual, n, args.max_iters,
                osd=osd, order=int(args.z_osd) if osd is not None else 0, layered=layered)
        if zv:
            patterns = zhat.reshape(n,9) @ BLOCK_WEIGHTS
            bit_probs = bit_table[patterns]
            projected_prior = projection_table[patterns]
        else:
            bit_probs, projected_prior = fx, marginal_projection
        x0, xi, xv = decode_z(sx, projected_bp, projected_prior, n, args.max_iters)
        xhat = likelihood_lift(x0, bit_probs)
        decoder_seconds = time.perf_counter()-tick
        if xv != np.array_equal((hz @ xhat)%2,sx) or zv != np.array_equal((hx @ zhat)%2,sz):
            raise ValueError('Physical syndrome does not match decoder status')
        if np.any((xhat^x0).reshape(n,9) != (xhat^x0).reshape(n,9)[:,:1]):
            raise ValueError('Representative difference must consist of all-one inner stabilizers')
        xr, zr = x ^ xhat, z ^ zhat
        xo, zo = bool(xv and member(xr,x_basis)), bool(zv and member(zr,z_basis))
        union = xr | zr
        record = dict(frame=frame, x_syndrome_ok=bool(xv), z_syndrome_ok=bool(zv),
                      x_logical_ok=xo, z_logical_ok=zo, logical_ok=xo and zo,
                      conditioned_prior_used=bool(zv), exact_x=not xr.any(), exact_z=not zr.any(),
                      x_residual_weight=int(xr.sum()), z_residual_weight=int(zr.sum()),
                      pauli_residual_weight=int(union.sum()),
                      x_physical_block_error=bool(xr.any()), z_physical_block_error=bool(zr.any()),
                      physical_block_error=bool(union.any()), x_iterations=xi, z_iterations=zi,
                      seconds=decoder_seconds, physical_x_only=int((x & (1-z)).sum()),
                      physical_y=int((x & z).sum()), physical_z_only=int((z & (1-x)).sum()))
        record.update(z_osd_candidates=osd_candidates, z_osd_seconds=osd_seconds)
        record.update(z_retries=retries, z_retry_iterations=retry_iterations, z_retry_seconds=retry_seconds)
        records.append(record)
        print(f'{frame+1}/{args.frames}: X syndrome={xv}, Z syndrome={zv}, logical={xo and zo}', flush=True)
    source_files=[Path(__file__), ROOT/'run_companion_staged_xz.py', ROOT/'run_companion_z_bp.py',
                  ROOT/'run_identity_outer_comparison.py', ROOT/'construct_companion_cascade.py']
    if args.native:
        source_files.append(ROOT/'companion_fast_bp.py')
    if osd is not None:
        source_files.append(ROOT/'companion_gf256_osd.py')
    if layered is not None:
        source_files.append(ROOT/'companion_layered_bp.py')
    result = dict(frames=args.frames, p_total=args.p_total, alpha=args.alpha, f_x=fx, f_z=fz,
                  seed=args.seed, max_iters=args.max_iters, mode='staged_Z_plus_G2_outer_X',
                  orientation='global_hadamard', parity_prior=args.parity_prior, kernel=kernel.__name__,
                  z_osd=args.z_osd,
                  z_layered_retry=args.z_layered_retry,
                  x_prior='conditioned on valid Z estimate, else marginal', x_lift='likelihood',
                  matrix_dir=str(args.matrix_dir.resolve()), matrix_manifest=manifest,
                  physical_hx_rank=len(x_basis), physical_hz_rank=len(z_basis),
                  noise_sha256=fingerprints.hexdigest(), elapsed_seconds=time.perf_counter()-start,
                  source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},
                  records=records)
    for axis in ('x','z'):
        result[f'{axis}_syndrome_failures'] = sum(not r[f'{axis}_syndrome_ok'] for r in records)
        result[f'{axis}_logical_failures'] = sum(not r[f'{axis}_logical_ok'] for r in records)
        result[f'wer_{axis}'] = sum(r[f'{axis}_physical_block_error'] for r in records)/args.frames
    for axis in ('x','z','pauli'):
        result[f'qber_{axis}'] = sum(r[f'{axis}_residual_weight'] for r in records)/(9*n*args.frames)
    result.update(logical_failures=sum(not r['logical_ok'] for r in records),
                  wer_physical=sum(r['physical_block_error'] for r in records)/args.frames,
                  conditioned_frames=sum(r['conditioned_prior_used'] for r in records),
                  undetected_logical_failures=sum(r['x_syndrome_ok'] and r['z_syndrome_ok']
                                                   and not r['logical_ok'] for r in records))
    result['logical_fer'] = result['logical_failures']/args.frames
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames',type=int,default=20)
    parser.add_argument('--p-total',type=float,default=.0705)
    parser.add_argument('--alpha',type=float,default=100.)
    parser.add_argument('--seed',type=int,default=20261310)
    parser.add_argument('--max-iters',type=int,default=100)
    parser.add_argument('--matrix-dir',type=Path,default=ROOT/'results/companion_binary_n9_k8_matrix20261019')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--parity-prior',action='store_true')
    parser.add_argument('--z-osd', choices=['none','0','1'], default='none',
                        help='Apply GF(256) order-0 or order-1 OSD only after failed Z BP')
    parser.add_argument('--z-layered-retry', action='store_true',
                        help='Retry failed Z BP with two damped layered schedules')
    parser.add_argument('--native',action='store_true')
    args = parser.parse_args()
    if not np.isfinite(args.alpha) or args.alpha <= 0 or not 0 < args.p_total < 1 or min(args.frames,args.max_iters)<1:
        parser.error('Require positive alpha and frame/iteration counts and p in (0,1)')
    if (args.z_osd != 'none' or args.z_layered_retry) and (args.native or not args.parity_prior):
        parser.error('Z fallback requires parity-aware Python BP')
    result = simulate(args)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    partial = args.output.with_suffix(args.output.suffix+'.partial')
    partial.write_text(json.dumps(result,indent=2)+'\n')
    partial.replace(args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','source_sha256','matrix_manifest')},indent=2))


if __name__ == '__main__':
    main()
