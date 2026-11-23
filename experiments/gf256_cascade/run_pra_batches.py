"""Run fixed-size independent Pauli-noise batches for the Fig. 10 comparison."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results/pra_alpha100_20261118'
MATRICES = [ROOT/'results'/name for name in
            ('companion_binary_n9_k8_matrix20261019', 'companion_binary_n9_k8_outer_path2')]


def run(job):
    path, command, count = job
    if path.is_file():
        old = json.loads(path.read_text())
        if old['frames'] != count:
            raise ValueError(f'Existing file has incompatible count: {path}')
        return
    tmp = path.with_suffix('.partial')
    with path.with_suffix('.log').open('w') as out:
        subprocess.run(command + ['--output', str(tmp)], check=True, stdout=out,
                       stderr=subprocess.STDOUT, env={**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'})
    tmp.replace(path)
    print(f'Finished {path.name}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['original','swapped','swapped_plain'])
    parser.add_argument('--p-total', type=float, required=True)
    parser.add_argument('--frames', type=int, required=True, help='Per matrix')
    parser.add_argument('--batch', type=int, default=100)
    parser.add_argument('--seed-start', type=int, required=True)
    parser.add_argument('--workers', type=int, default=4)
    a = parser.parse_args()
    if min(a.frames,a.batch,a.workers)<1 or a.frames % a.batch:
        parser.error('Require positive sizes and an integer number of complete batches')
    directory = OUT/a.mode
    directory.mkdir(parents=True, exist_ok=True)
    jobs = []
    for b in range(a.frames//a.batch):
        seed = a.seed_start + b
        for m, matrix in enumerate(MATRICES,1):
            path = directory/f'm{m}_p{a.p_total}_s{seed}.json'
            if a.mode == 'original':
                command = [sys.executable,str(ROOT/'run_companion_staged_xz.py'),
                           '--compare-lifts','--z-prior','conditioned']
            else:
                command = [sys.executable,str(ROOT/'run_companion_swapped.py')]
                if a.mode=='swapped':
                    command.append('--parity-prior')
            command += ['--frames',str(a.batch),'--p-total',str(a.p_total),'--alpha','100',
                        '--seed',str(seed),'--max-iters','100','--matrix-dir',str(matrix)]
            jobs.append((path,command,a.batch))
    with ThreadPoolExecutor(a.workers) as pool:
        list(pool.map(run,jobs))


if __name__=='__main__':
    main()
