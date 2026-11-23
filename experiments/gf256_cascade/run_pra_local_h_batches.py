"""Run the frozen Hadamard-mask comparison in resumable 100-frame batches."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

from run_pra_batches import MATRICES

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results/pra_alpha100_local_h_20261119'
PROTOCOL = (
    (100, .07, 1000, 20261500, 'paired'),
    (100, .07, 1000, 20261700, 'independent'),
    (100, .08, 1000, 20261710, 'curve'),
    (100, .06, 400, 20261720, 'curve'),
    (100, .075, 400, 20261724, 'curve'),
    (100, .09, 400, 20261728, 'curve'),
    (100, .10, 400, 20261732, 'curve'),
    (10, .07, 400, 20261736, 'bias'),
    (1, .07, 400, 20261740, 'bias'),
)
MASK = [1,3,5,7,9]


def run(job, *, out_dir=OUT):
    alpha, p, seed, label, m = job
    output = out_dir/f'{label}_m{m}_a{alpha}_p{p}_s{seed}.json'
    if output.is_file():
        saved = json.loads(output.read_text())
        actual = [saved[k] for k in ('frames','alpha','p_total','seed','hadamards','max_iters')]
        if actual != [100,alpha,p,seed,MASK,100]:
            raise ValueError(f'Existing batch has incompatible settings: {output}')
        expected = json.loads((MATRICES[m-1]/'manifest.json').read_text())
        if saved['matrix_manifest'] != expected:
            raise ValueError(f'Existing batch has a different matrix: {output}')
        return
    cmd = [sys.executable, str(ROOT/'run_companion_local_hadamard.py'), '--alpha',str(alpha),
           '--p-total',str(p),'--seed',str(seed),'--frames','100','--max-iters','100',
           '--matrix-dir',str(MATRICES[m-1]),'--output',str(output),
           '--hadamards',*map(str,MASK)]
    with output.with_suffix('.log').open('w') as log:
        subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT,
                       env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    print(f'Finished {output.name}',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('Require a positive worker count')
    OUT.mkdir(parents=True,exist_ok=True)
    jobs = [(alpha,p,seed,label,m) for alpha,p,frames,start,label in PROTOCOL
            for seed in range(start,start+frames//100) for m in (1,2)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run,jobs))


if __name__ == '__main__':
    main()
