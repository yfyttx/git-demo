#!/usr/bin/env python3
"""Compare the [9,8] inner code with independent tail qubits on paired Pauli frames."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy import sparse
from scipy.stats import binomtest

from construct_companion_cascade import ROOT, expand, mul, read_outer
from run_companion_local_hadamard import simulate as simulate_cascade
from run_companion_swapped import BIT_WEIGHTS, HEAD_BITS, trace
from run_companion_z_bp import OuterBP, member, row_basis
from run_pra_batches import MATRICES
from run_pra_fig12_wer5000 import upper95
from run_pra_local_h_batches import MASK


def head_prior(probabilities):
    """Distribution over the 256 binary patterns on the eight data qubits."""
    probabilities = np.asarray(probabilities)[..., None, :8]
    return np.prod(np.where(HEAD_BITS, probabilities, 1 - probabilities), axis=-1)


def baseline_checks(delta_rows, gamma_rows, n):
    selector = np.column_stack((np.eye(8, dtype=np.int64), np.zeros(8, dtype=np.int64)))
    tail = np.eye(9, dtype=np.int64)[-1:]
    hx = expand(delta_rows, n, selector)
    hz = sparse.vstack(
        (expand(gamma_rows, n, selector, transpose=True),
         sparse.kron(sparse.eye(n, dtype=np.int64), tail)),
        format='csr',
    )
    if np.any((hx @ hz.T).data % 2):
        raise ValueError('Baseline CSS checks do not commute')
    return hx, hz


def simulate_baseline(args):
    n, delta_rows = read_outer(args.matrix_dir/'HZ_TIT_Outer_For_HX_GF256.txt')
    n_gamma, gamma_rows = read_outer(args.matrix_dir/'HX_TIT_Outer_For_HZ_GF256.txt')
    if n != n_gamma:
        raise ValueError('Mismatched outer lengths')
    hx, hz = baseline_checks(delta_rows, gamma_rows, n)
    x_basis, z_basis = row_basis(hx), row_basis(hz)
    source_matrix = json.loads((args.matrix_dir/'manifest.json').read_text())
    if (hx.shape[1], len(x_basis), len(z_basis)) != (
        source_matrix['physical_qubits'], source_matrix['hx_rank'], source_matrix['hz_rank']
    ):
        raise ValueError('Baseline does not match the cascade length and check ranks')
    delta, gamma = OuterBP(delta_rows, n), OuterBP(gamma_rows, n)
    dual = np.array([sum(trace(mul(v, 1 << i)) << i for i in range(8)) for v in range(256)])
    inverse_dual = np.argsort(dual)

    mask = np.isin(np.arange(1, 10), args.hadamards)
    px = py = args.p_total/(args.alpha + 2)
    pz = args.alpha*px
    x_marginal = np.where(mask, pz + py, px + py)
    z_marginal = np.where(mask, px + py, pz + py)
    z_only = np.where(mask, px, pz)
    x_prior = head_prior(x_marginal)[dual]
    z_prior = head_prior(z_marginal)
    conditional = np.stack((z_only/(1 - x_marginal), py/x_marginal))
    patterns = np.column_stack((HEAD_BITS, np.zeros(len(HEAD_BITS), dtype=int)))
    conditional_prior = head_prior(conditional[patterns, np.arange(9)])

    rng = np.random.default_rng(args.seed)
    noise_hash = hashlib.sha256()
    records = []
    for frame in range(args.frames):
        # Match the physical Pauli draws used by the existing cascade simulator.
        draw = rng.random((n, 9))
        physical_x = (draw < px) | ((draw >= pz + py) & (draw < args.p_total))
        physical_z = (draw >= px) & (draw < args.p_total)
        noise_hash.update((physical_x + 2*physical_z).astype(np.uint8).tobytes())
        x = np.where(mask, physical_z, physical_x).ravel().astype(np.int64)
        z = np.where(mask, physical_x, physical_z).ravel().astype(np.int64)
        sx, sz = (hz @ x) % 2, (hx @ z) % 2

        start = perf_counter()
        sx_outer = inverse_dual[sx[:-n].reshape(-1, 8) @ BIT_WEIGHTS]
        sz_outer = sz.reshape(-1, 8) @ BIT_WEIGHTS
        x_symbols, xi, xv = gamma.decode(sx_outer, x_prior, args.max_iters)
        x_head = dual[x_symbols]
        # The last syndrome rows measure each independent auxiliary qubit.
        x_tail = sx[-n:]
        xhat = np.column_stack((HEAD_BITS[x_head], x_tail)).ravel()
        z_symbols, zi, zv = delta.decode(
            sz_outer, conditional_prior[x_head] if xv else z_prior, args.max_iters
        )
        zhat = np.column_stack((HEAD_BITS[z_symbols], conditional[x_tail, 8] > 0.5)).ravel()
        seconds = perf_counter() - start

        if xv != np.array_equal((hz @ xhat) % 2, sx) or zv != np.array_equal((hx @ zhat) % 2, sz):
            raise ValueError('Baseline field decoder disagrees with the binary checks')
        xr, zr = x ^ xhat, z ^ zhat
        logical_ok = bool(xv and zv and member(xr, x_basis) and member(zr, z_basis))
        records.append(dict(frame=frame, staged_syndrome_ok=bool(xv), projected_syndrome_ok=bool(zv),
                            logical_ok=logical_ok, undetected=bool(xv and zv and not logical_ok),
                            physical_residual_weight=int((xr | zr).sum()), seconds=seconds,
                            staged_iterations=xi, projected_iterations=zi))

    return dict(noise_sha256=noise_hash.hexdigest(), records=records)


def write_json(path, data):
    tmp = path.with_suffix('.partial')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def run_pair(matrix, seed, settings, out, source_sha256):
    args = argparse.Namespace(**settings, seed=seed, matrix_dir=MATRICES[matrix - 1])
    settings_for_pair = dict(settings, matrix=matrix, seed=seed)
    output = out/f'm{matrix}_s{seed}.json'
    metadata = dict(settings=settings_for_pair, source_sha256=source_sha256,
                    matrix_manifest=json.loads((args.matrix_dir/'manifest.json').read_text()))
    if output.is_file():
        data = json.loads(output.read_text())
        if any(data[name] != value for name, value in metadata.items()):
            raise ValueError(f'Existing batch has different settings, matrix or source: {output}')
        return data

    runs = {}
    # Alternate order to reduce systematic timing differences from startup.
    methods = [('baseline', simulate_baseline), ('cascade', simulate_cascade)]
    for name, simulate in methods[::1 if seed % 2 == 0 else -1]:
        start = process_time()
        run = simulate(args)
        runs[name] = dict(cpu_seconds=process_time() - start,
                          records=run['records'], noise_sha256=run['noise_sha256'])
    if runs['baseline']['noise_sha256'] != runs['cascade']['noise_sha256']:
        raise ValueError(f'Physical noise differs between paired batches: {output}')
    data = dict(metadata, **runs)
    write_json(output, data)
    return data


def summarize(batches):
    arms = ('baseline', 'cascade')
    records = {arm: [row for batch in batches for row in batch[arm]['records']] for arm in arms}
    summary = {}
    frames = sum(batch['settings']['frames'] for batch in batches)
    for arm in arms:
        rows = records[arm]
        failures = sum(not row['logical_ok'] for row in rows)
        summary[arm] = dict(frames=frames, logical_failures=failures, logical_fer=failures/frames,
                            logical_fer_upper95=upper95(failures, frames),
                            x_syndrome_failures=sum(not row['staged_syndrome_ok'] for row in rows),
                            z_syndrome_failures=sum(not row['projected_syndrome_ok'] for row in rows),
                            undetected_logical_failures=sum(row['undetected'] for row in rows),
                            physical_error_frames=sum(row['physical_residual_weight'] > 0 for row in rows),
                            batch_cpu_seconds_per_frame=sum(b[arm]['cpu_seconds'] for b in batches)/frames,
                            mean_x_bp_iterations=sum(row['staged_iterations'] for row in rows)/frames,
                            mean_z_bp_iterations=sum(row['projected_iterations'] for row in rows)/frames)
    pairs = list(zip(records['baseline'], records['cascade'], strict=True))
    rescued = sum(not a['logical_ok'] and b['logical_ok'] for a, b in pairs)
    deteriorated = sum(a['logical_ok'] and not b['logical_ok'] for a, b in pairs)
    summary.update(cascade_only_success=rescued, baseline_only_success=deteriorated,
                   both_failed=sum(not a['logical_ok'] and not b['logical_ok'] for a, b in pairs),
                   mcnemar_exact_two_sided=binomtest(rescued, rescued + deteriorated).pvalue
                   if rescued + deteriorated else 1.0)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p-total', type=float, default=0.120)
    parser.add_argument('--alpha', type=float, default=100)
    parser.add_argument('--frames-per-matrix', type=int, default=300)
    parser.add_argument('--batch-size', type=int, default=100)
    parser.add_argument('--seed-start', type=int, default=2026112200)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if (min(args.batch_size, args.frames_per_matrix, args.workers) < 1
            or args.frames_per_matrix % args.batch_size or args.seed_start < 0
            or not 0 < args.p_total < 1 or not np.isfinite(args.alpha) or args.alpha <= 0):
        parser.error('Invalid channel, seed or batch sizes')
    settings = dict(frames=args.batch_size, alpha=args.alpha, p_total=args.p_total,
                    hadamards=MASK, max_iters=100)
    seeds = range(args.seed_start, args.seed_start + args.frames_per_matrix//args.batch_size)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = ('run_inner_code_comparison.py', 'run_companion_local_hadamard.py', 'run_companion_swapped.py',
             'run_companion_z_bp.py', 'run_companion_staged_xz.py', 'construct_companion_cascade.py',
             'run_identity_outer_comparison.py')
    source_sha256 = {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(run_pair, m, s, settings, args.output_dir, source_sha256): (m, s)
                for s in seeds for m in (1, 2)}
        for done, job in enumerate(as_completed(jobs), 1):
            matrix, seed = jobs[job]
            result = results[matrix, seed] = job.result()
            counts = [sum(not row['logical_ok'] for row in result[key]['records'])
                      for key in ('baseline', 'cascade')]
            print(f'[{done}/{len(jobs)}] matrix={matrix} seed={seed} '
                  f'baseline={counts[0]}/{args.batch_size} cascade={counts[1]}/{args.batch_size}', flush=True)
    for seed in seeds:
        if results[1, seed]['baseline']['noise_sha256'] != results[2, seed]['baseline']['noise_sha256']:
            raise ValueError(f'Physical noise differs across matrices for seed {seed}')
    summary = dict(settings=settings, seeds=list(seeds), matrices={})
    for matrix in (1, 2):
        batches = [results[matrix, seed] for seed in seeds]
        manifest = batches[0]['matrix_manifest']
        summary['matrices'][str(matrix)] = dict(
            physical_qubits=manifest['physical_qubits'], logical_qubits=manifest['logical_qubits'],
            quantum_rate=manifest['quantum_rate'], **summarize(batches))
    write_json(args.output_dir/'summary.json', summary)


if __name__ == '__main__':
    main()
