#!/usr/bin/env python3
"""Validate one Fig. 12 working point at alpha=100 with 5000 frames per matrix."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
from functools import partial
import json
from pathlib import Path

from scipy.stats import beta

from run_pra_local_h_batches import MASK, run


OUT = Path(__file__).resolve().parent/'results/pra_fig12_wer5000'
DEFAULT_P_TOTAL = 0.08
ALPHA = 100
DEFAULT_SEED_START = 20261900
BATCHES_PER_MATRIX = 50
TARGET_WER = 1e-3


def upper95(events, frames):
    """Exact one-sided Clopper–Pearson upper bound."""
    return 1.0 if events == frames else float(beta.ppf(0.95, events+1, frames-events))


def summarize(out_dir, p_total, seeds):
    results = []
    paired_noise = {}
    source_hashes = None
    for matrix in (1, 2):
        records = []
        batches = []
        for seed in seeds:
            batch_path = out_dir/f'fig12_m{matrix}_a{ALPHA}_p{p_total}_s{seed}.json'
            batch = json.loads(batch_path.read_text())
            batch_records = batch['records']
            if len(batch_records) != 100 or any(r['frame'] != i for i, r in enumerate(batch_records)):
                raise ValueError(f'Invalid frame records in {batch_path}')
            noise_hash = batch['noise_sha256']
            if paired_noise.setdefault(seed, noise_hash) != noise_hash:
                raise ValueError(f'Noise differs across matrices for seed {seed}')
            if source_hashes is None:
                source_hashes = batch['source_sha256']
            elif source_hashes != batch['source_sha256']:
                raise ValueError('Simulator source differs across batches. Use one frozen simulator version.')
            records.extend(batch_records)
            batches.append(batch_path.name)
        frames = len(records)
        n = batch['matrix_manifest']['physical_qubits']
        wer_events = sum(r['physical_residual_weight'] > 0 for r in records)
        logical_events = sum(not r['logical_ok'] for r in records)
        wer_upper = upper95(wer_events, frames)
        results.append(dict(matrix=matrix, physical_qubits=n,
                            quantum_rate=batch['matrix_manifest']['quantum_rate'],
                            frames=frames, physical_error_frames=wer_events,
                            physical_wer=wer_events/frames, physical_wer_upper_95=wer_upper,
                            supports_wer_below_target=wer_upper < TARGET_WER,
                            physical_qber=sum(r['physical_residual_weight'] for r in records)/(n*frames),
                            logical_failures=logical_events, logical_fer=logical_events/frames,
                            logical_fer_upper_95=upper95(logical_events,frames),
                            staged_syndrome_failures=sum(not r['staged_syndrome_ok'] for r in records),
                            projected_syndrome_failures=sum(not r['projected_syndrome_ok'] for r in records),
                            undetected_logical_failures=sum(r['undetected'] for r in records),
                            batches=batches))
    return dict(p_total=p_total, alpha=ALPHA, hadamards=MASK, max_iters=100,
                target_physical_wer=TARGET_WER, noise_seeds=list(seeds),
                source_sha256=source_hashes, matrices=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p-total', type=float, default=DEFAULT_P_TOTAL,
                        help='total physical Pauli error probability (default: %(default)s)')
    parser.add_argument('--seed-start', type=int, default=DEFAULT_SEED_START,
                        help='first of 50 consecutive noise seeds (default: %(default)s)')
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if not 0 < args.p_total < 1:
        parser.error('--p-total must be between 0 and 1')
    if args.seed_start < 0:
        parser.error('--seed-start must be nonnegative')
    seeds = range(args.seed_start, args.seed_start+BATCHES_PER_MATRIX)
    out_dir = OUT
    if (args.p_total, args.seed_start) != (DEFAULT_P_TOTAL, DEFAULT_SEED_START):
        out_dir = OUT.with_name(f'{OUT.name}_p{args.p_total}_s{args.seed_start}')
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(ALPHA, args.p_total, seed, 'fig12', matrix) for seed in seeds for matrix in (1, 2)]
    with (out_dir/'.run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error(f'A run already holds the lock in {out_dir}')
        pool = ThreadPoolExecutor(max_workers=args.workers)
        try:
            submit = partial(run, out_dir=out_dir)
            pending = [pool.submit(submit, job) for job in jobs]
            for done, future in enumerate(as_completed(pending), 1):
                future.result()
                print(f'[{done}/{len(jobs)}] validated batches', flush=True)
        except KeyboardInterrupt:
            print('Interrupted. Run the same command again to resume completed batches.', flush=True)
            raise SystemExit(130)
        finally:
            pool.shutdown(wait=True,cancel_futures=True)
        summary = summarize(out_dir, args.p_total, seeds)
        temporary = out_dir/'summary.json.partial'
        temporary.write_text(json.dumps(summary,indent=2)+'\n')
        temporary.replace(out_dir/'summary.json')
        for result in summary['matrices']:
            print(f"Matrix {result['matrix']}: WER {result['physical_error_frames']}/{result['frames']}, "
                  f"one-sided 95% upper {result['physical_wer_upper_95']:.6%}, "
                  f"target met: {result['supports_wer_below_target']}")
        print(f"Summary: {out_dir/'summary.json'}")


if __name__ == '__main__':
    main()
