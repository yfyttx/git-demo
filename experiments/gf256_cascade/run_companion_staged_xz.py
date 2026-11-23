#!/usr/bin/env python3
"""Full Pauli simulation: staged X recovery and direct G2-to-Delta Z BP."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy import sparse
from construct_companion_cascade import ROOT, read_outer, mul
from run_companion_z_bp import OuterBP, row_basis, member
from run_identity_outer_comparison import trace


def parity_conditioned_head_prior(head_prior, tail_flip_probability):
    """Return P(head pattern | measured nine-bit parity) for both parities."""
    head_parity = np.array([pattern.bit_count() & 1 for pattern in range(256)])
    tail_bits = np.arange(2)[:, None] ^ head_parity
    tail_likelihood = np.where(tail_bits, tail_flip_probability, 1 - tail_flip_probability)
    table = head_prior[None, :] * tail_likelihood
    return table / table.sum(axis=1, keepdims=True)


def decode_x(sx, bp, prior, inverse_dual, dual, n, iterations, parity_prior=None):
    # Only measured syndrome and channel probabilities enter this function.
    upper_size = 8 * len(bp.cols)
    upper = sx[:upper_size].reshape(-1, 8) @ (1 << np.arange(8))
    lower = sx[upper_size:]
    symbol_prior = prior[dual] if parity_prior is None else parity_prior[lower][:, dual]
    symbols, used, valid = bp.decode(inverse_dual[upper], symbol_prior, iterations)
    estimate = np.zeros((n, 9), dtype=np.int64)
    estimate[:, :8] = (dual[symbols, None] >> np.arange(8)) & 1
    estimate[:, 8] = lower ^ (estimate[:, :8].sum(axis=1) % 2)
    return estimate.ravel(), used, valid



def decode_binary_bp_fallback(h, syndrome, p, max_iterations):
    """Full binary Tanner-graph sum-product BP; no truth input."""
    h = h.tocsr()
    edge_rows = np.repeat(np.arange(h.shape[0]), np.diff(h.indptr))
    edge_vars = h.indices.copy()
    by_row = [np.arange(h.indptr[r], h.indptr[r+1]) for r in range(h.shape[0])]
    llr0 = np.log((1-p)/p)
    v2c = np.full(len(edge_vars), llr0, dtype=float)
    c2v = np.zeros_like(v2c)
    syndrome = np.asarray(syndrome, dtype=np.int64)
    for iteration in range(1, max_iterations+1):
        for r, edges in enumerate(by_row):
            vals = np.tanh(v2c[edges]/2)
            prefix = np.concatenate(([1.0], np.cumprod(vals[:-1])))
            suffix = np.concatenate((np.cumprod(vals[:0:-1])[::-1], [1.0]))
            prod = prefix*suffix*(-1.0 if syndrome[r] else 1.0)
            c2v[edges] = 2*np.arctanh(np.clip(prod, -1+1e-12, 1-1e-12))
        posterior = llr0 + np.bincount(edge_vars, weights=c2v, minlength=h.shape[1])
        estimate = (posterior < 0).astype(np.int64)
        parity = np.bincount(edge_rows, weights=estimate[edge_vars], minlength=h.shape[0]).astype(np.int64)&1
        if np.array_equal(parity, syndrome):
            return estimate, iteration, True
        v2c = posterior[edge_vars]-c2v
    return estimate, max_iterations, False


def decode_z(sz, bp, prior, n, iterations):
    symbols, used, valid = bp.decode(sz.reshape(-1, 8) @ (1 << np.arange(8)), prior, iterations)
    estimate = np.zeros((n, 9), dtype=np.int64)
    estimate[:, :8] = (symbols[:, None] >> np.arange(8)) & 1
    return estimate.ravel(), used, valid



def z_bit_prior(x_estimate, p_total, alpha):
    """P(Z component=1 | estimated X component) for pX=pY and pZ=alpha*pX."""
    x_estimate = np.asarray(x_estimate).reshape(-1, 9)
    if not np.all((x_estimate == 0) | (x_estimate == 1)):
        raise ValueError('X estimate must be binary')
    p_x = p_total / (2 + alpha)
    p_z = alpha * p_x
    return np.where(x_estimate == 1, .5, p_z / (1 - 2*p_x))


def conditioned_z_prior(x_estimate, p_depolar, alpha=1.0):
    """G2=[I8|1] pushforward conditioned on an estimated X block, not truth.

    This is a hard plug-in channel model, not the exact posterior given X decoding.
    """
    p = z_bit_prior(x_estimate, p_depolar, alpha)
    bits = (np.arange(256)[:, None] >> np.arange(8)) & 1
    head = np.prod(np.where(bits[None, :, :], p[:, None, :8], 1-p[:, None, :8]), axis=-1)
    prior = (1-p[:, 8, None])*head + p[:, 8, None]*head[:, ::-1]
    if not np.allclose(prior.sum(axis=1), 1):
        raise ValueError('Conditional prior normalization')
    return prior


def likelihood_lift(zero_tail_estimate, probabilities):
    """Choose between the two binary words with the decoded projection."""
    blocks = zero_tail_estimate.reshape(-1, 9)
    llr = np.log1p(-probabilities) - np.log(probabilities)
    prefer_complement = np.sum((1 - 2*blocks) * llr, axis=1) < 0
    return (blocks ^ prefer_complement[:, None]).ravel()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--x-fallback', choices=['none', 'binary_bp', 'parity_gf256'], default='none',
                        help='Retry X with full binary BP or GF(256) BP conditioned on measured block parity')
    parser.add_argument('--z-prior', choices=['marginal', 'conditioned'], default='marginal')
    parser.add_argument('--frames', type=int, default=20)
    parser.add_argument('--p-depolar', '--p-total', type=float, default=.07125)
    parser.add_argument('--alpha', type=float, default=1.0, help='pZ/pX with pY=pX')
    parser.add_argument('--compare-lifts', action='store_true',
                        help='Evaluate the likelihood representative alongside the historical zero-tail output')
    parser.add_argument('--seed', type=int, default=20261120)
    parser.add_argument('--max-iters', type=int, default=100)
    parser.add_argument('--matrix-dir', type=Path, default=ROOT/'results/companion_binary_n9_k8_matrix20261019')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.frames < 1 or args.max_iters < 1 or not 0 < args.p_depolar < 1
            or not np.isfinite(args.alpha) or args.alpha <= 0):
        raise ValueError('Invalid simulation parameters')
    n, delta = read_outer(args.matrix_dir/'HZ_TIT_Outer_For_HX_GF256.txt')
    nb, gamma = read_outer(args.matrix_dir/'HX_TIT_Outer_For_HZ_GF256.txt')
    hx = sparse.load_npz(args.matrix_dir/'HX_binary.npz').astype(np.int64)
    hz = sparse.load_npz(args.matrix_dir/'HZ_binary.npz').astype(np.int64)
    g = np.loadtxt(args.matrix_dir/'G2_binary.txt', dtype=int)
    h = np.loadtxt(args.matrix_dir/'H2_binary.txt', dtype=int, ndmin=2)
    if nb != n or not np.array_equal(g, np.hstack((np.eye(8,dtype=int),np.ones((8,1),dtype=int)))) or not np.array_equal(h,np.ones((1,9),dtype=int)):
        raise ValueError('Requires G2=[I8|1], H2=1x9 ones and matched outer lengths')
    if hx.shape != (8*len(delta),9*n) or hz.shape != (8*len(gamma)+n,9*n):
        raise ValueError('Unexpected physical check dimensions')
    product = hx@hz.T
    if np.any(product.data % 2): raise ValueError('CSS orthogonality failed')
    # These bases are used only by the evaluator, never by the decoders.
    x_stabilizers, z_stabilizers = row_basis(hx), row_basis(hz)
    dual = np.array([sum(trace(mul(v,1<<i))<<i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)
    bx, bz = OuterBP(gamma,n), OuterBP(delta,n)
    pd = args.p_depolar
    pure_x = pd / (2 + args.alpha)
    pure_z = args.alpha * pure_x
    fm = 2*pure_x
    z_flip = pure_x + pure_z
    weight = np.array([v.bit_count() for v in range(256)])
    px = fm**weight * (1-fm)**(8-weight)
    parity_table = parity_conditioned_head_prior(px, fm) if args.x_fallback == 'parity_gf256' else None
    head_z = z_flip**weight * (1-z_flip)**(8-weight)
    pz = (1-z_flip)*head_z + z_flip*head_z[::-1]
    rng = np.random.default_rng(args.seed)
    records = []
    start = time.perf_counter()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for frame in range(args.frames):
        draw = rng.random((n,9))
        # X, Z, Y, I intervals. Keep the original symmetric-frame replay order.
        x = ((draw < pure_x) | ((draw >= z_flip) & (draw < pd))).astype(np.int64)
        z = ((draw >= pure_x) & (draw < pd)).astype(np.int64)
        sx,sz = (hz@x.ravel())%2,(hx@z.ravel())%2
        # Audit physical/field mappings independently of decoder selection.
        ux = inverse_dual[x[:,:8]@(1<<np.arange(8))]
        uz = ((z@g.T)%2)@(1<<np.arange(8))
        if not np.array_equal(bx.syndrome(ux),inverse_dual[sx[:8*len(gamma)].reshape(-1,8)@(1<<np.arange(8))]):
            raise ValueError('X upper mapping mismatch')
        if not np.array_equal(bz.syndrome(uz),sz.reshape(-1,8)@(1<<np.arange(8))):
            raise ValueError('Z G2 mapping mismatch')
        t=time.perf_counter()
        xhat,xi,xvalid=decode_x(sx,bx,px,inverse_dual,dual,n,args.max_iters)
        x_fallback_attempted = args.x_fallback != 'none' and not xvalid
        x_fallback_succeeded = False
        x_fallback_iterations = 0
        if x_fallback_attempted:
            if args.x_fallback == 'binary_bp':
                candidate, x_fallback_iterations, valid = decode_binary_bp_fallback(hz, sx, fm, args.max_iters)
            else:
                candidate, x_fallback_iterations, valid = decode_x(
                    sx, bx, px, inverse_dual, dual, n, args.max_iters, parity_prior=parity_table)
            x_fallback_succeeded = valid and np.array_equal((hz @ candidate) % 2, sx)
            if x_fallback_succeeded:
                xhat, xi, xvalid = candidate, x_fallback_iterations, True
        condition_used = args.z_prior == 'conditioned' and xvalid
        zprior = conditioned_z_prior(xhat, pd, args.alpha) if condition_used else pz
        zhat,zi,zvalid=decode_z(sz,bz,zprior,n,args.max_iters)
        if xvalid != np.array_equal((hz@xhat)%2,sx) or zvalid != np.array_equal((hx@zhat)%2,sz):
            raise ValueError('Decoded physical syndrome mismatch')
        xok=xvalid and member(x.ravel()^xhat,x_stabilizers)
        zok=zvalid and member(z.ravel()^zhat,z_stabilizers)
        x_residual = x.ravel() ^ xhat
        z_residual = z.ravel() ^ zhat
        pauli_residual = x_residual | z_residual
        r=dict(frame=frame,x_fallback_attempted=x_fallback_attempted,
               x_fallback_succeeded=x_fallback_succeeded,x_fallback_iterations=x_fallback_iterations,
               conditioned_prior_used=condition_used,x_syndrome_ok=bool(xvalid),z_syndrome_ok=bool(zvalid),
               x_logical_ok=bool(xok),z_logical_ok=bool(zok),logical_ok=bool(xok and zok),
               exact_x=bool(np.array_equal(x.ravel(),xhat)),exact_z=bool(np.array_equal(z.ravel(),zhat)),
               exact_z_projection=bool(np.array_equal((zhat.reshape(n,9)@g.T)%2,(z@g.T)%2)),
               x_residual_weight=int(x_residual.sum()),z_residual_weight=int(z_residual.sum()),
               pauli_residual_weight=int(pauli_residual.sum()),
               x_physical_block_error=bool(np.any(x_residual)),z_physical_block_error=bool(np.any(z_residual)),
               physical_block_error=bool(np.any(pauli_residual)),
               x_iterations=xi,z_iterations=zi,seconds=time.perf_counter()-t)
        if args.compare_lifts:
            probabilities = z_bit_prior(xhat, pd, args.alpha) if condition_used else z_flip
            lifted = likelihood_lift(zhat, probabilities)
            delta = lifted ^ zhat
            if np.any((hx @ delta) % 2) or not member(delta, z_stabilizers):
                raise ValueError('Lift changed syndrome or stabilizer class')
            lift_z_residual = lifted ^ z.ravel()
            lift_pauli = x_residual | lift_z_residual
            r.update(lift_z_weight=int(lift_z_residual.sum()), lift_pauli_weight=int(lift_pauli.sum()),
                     lift_complemented_blocks=int(delta.reshape(n, 9)[:, 8].sum()))
        r.update(physical_x_only=int(np.count_nonzero(x & (1-z))),
                 physical_y=int(np.count_nonzero(x & z)), physical_z_only=int(np.count_nonzero((1-x) & z)))
        records.append(r)
        print(f'{frame+1}/{args.frames}: X syndrome={xvalid}, Z syndrome={zvalid}, logical={r["logical_ok"]}',flush=True)
    summary=dict(frames=args.frames,p_D=pd if args.alpha == 1 else None,p_total=pd,
                 f_m=fm if args.alpha == 1 else None, alpha=args.alpha,
                 f_x=fm,f_z=z_flip,seed=args.seed,max_iters=args.max_iters,
                 mode=f'staged_X_plus_G2_outer_Z; x_fallback={args.x_fallback}',
                 z_prior=args.z_prior,
                 x_fallback=args.x_fallback,
                 conditioned_frames=sum(r['conditioned_prior_used'] for r in records),
                 x_syndrome_failures=sum(not r['x_syndrome_ok'] for r in records),
                 x_fallback_attempts=sum(r['x_fallback_attempted'] for r in records),
                 x_fallback_successes=sum(r['x_fallback_succeeded'] for r in records),
                 z_syndrome_failures=sum(not r['z_syndrome_ok'] for r in records),
                 x_logical_failures=sum(not r['x_logical_ok'] for r in records),
                 z_logical_failures=sum(not r['z_logical_ok'] for r in records),
                 logical_failures=sum(not r['logical_ok'] for r in records),
                 undetected_logical_failures=sum(r['x_syndrome_ok'] and r['z_syndrome_ok'] and not r['logical_ok'] for r in records),
                 exact_x_frames=sum(r['exact_x'] for r in records),
                 qber_x=sum(r['x_residual_weight'] for r in records)/(args.frames*9*n),
                 qber_z=sum(r['z_residual_weight'] for r in records)/(args.frames*9*n),
                 qber_pauli=sum(r['pauli_residual_weight'] for r in records)/(args.frames*9*n),
                 wer_x=sum(r['x_physical_block_error'] for r in records)/args.frames,
                 wer_z=sum(r['z_physical_block_error'] for r in records)/args.frames,
                 wer_physical=sum(r['physical_block_error'] for r in records)/args.frames,
                 elapsed_seconds=time.perf_counter()-start,
                 source_sha256={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'run_companion_z_bp.py',ROOT/'run_identity_outer_comparison.py']},
                 matrix_manifest=json.loads((args.matrix_dir/'manifest.json').read_text()),records=records)
    summary['logical_fer']=summary['logical_failures']/args.frames
    if args.compare_lifts:
        summary['likelihood_lift'] = dict(
            qber_z=sum(r['lift_z_weight'] for r in records)/(args.frames*9*n),
            qber_pauli=sum(r['lift_pauli_weight'] for r in records)/(args.frames*9*n),
            wer_z=sum(r['lift_z_weight'] > 0 for r in records)/args.frames,
            wer_physical=sum(r['lift_pauli_weight'] > 0 for r in records)/args.frames)
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['records','matrix_manifest','source_sha256']},indent=2))

if __name__=='__main__':main()
