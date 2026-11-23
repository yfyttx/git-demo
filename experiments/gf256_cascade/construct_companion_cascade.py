#!/usr/bin/env python3
"""Binary single-parity-check inner [9,8,2] + GF(256) outer companion-block construction.
Separate format/experiment: do not pass these binary files to the old GF256 cascade runner.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parent

def mul(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        if a & 256:
            a ^= 0x11d
        b >>= 1
    return value

EXP = [1]
for _ in range(254):
    EXP.append(mul(EXP[-1], 2))

def read_outer(path):
    lines = path.read_text().splitlines()
    n, m, q = map(int, lines[0].split())
    if q != 256:
        raise ValueError('Expected GF(256) ESSAI file')
    rows = []
    for line in lines[3:]:
        tokens = list(map(int, line.split()))
        if len(tokens) % 2:
            raise ValueError('Odd ESSAI row length')
        rows.append([(tokens[i]-1, EXP[tokens[i+1] % 255]) for i in range(0,len(tokens),2)])
    if len(rows) != m:
        raise ValueError('Wrong row count')
    return n, rows

def multiplication_matrix(a):
    # Column coordinates in basis (1, alpha, ..., alpha^7).
    return np.array([[(mul(a, 1 << j) >> i) & 1 for j in range(8)]
                     for i in range(8)], dtype=np.int64)

def expand(rows, n, inner, transpose=False):
    rr, cc = [], []
    for i, row in enumerate(rows):
        for j, a in row:
            block = multiplication_matrix(a)
            if transpose:
                block = block.T
            block = (block @ inner) % 2
            r, c = np.nonzero(block)
            rr.extend((8*i+r).tolist())
            cc.extend((inner.shape[1]*j+c).tolist())
    return sparse.csr_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc)),
                             shape=(8*len(rows),inner.shape[1]*n))

def rank_binary(matrix):
    pivots = {}
    for i in range(matrix.shape[0]):
        value = 0
        for col in matrix.indices[matrix.indptr[i]:matrix.indptr[i+1]]:
            value ^= 1 << int(col)
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer-dir', type=Path, default=ROOT/'results/matrix_seed_20261019_path_n12/path/build')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'results/companion_binary_n9_k8_matrix20261019')
    parser.add_argument("--identity-inner", action="store_true", help="Use G2=I8 and no inner checks")
    args = parser.parse_args()
    if args.identity_inner and args.output_dir == ROOT/"results/companion_binary_n9_k8_matrix20261019":
        args.output_dir = ROOT/"results/companion_identity_n8_k8_matrix20261019"
    a_path = args.outer_dir/'HZ_TIT_Outer_For_HX_GF256.txt'
    b_path = args.outer_dir/'HX_TIT_Outer_For_HZ_GF256.txt'
    n, a = read_outer(a_path)
    nb, b = read_outer(b_path)
    if n != nb:
        raise ValueError('Outer lengths differ')
    # H2 is one all-one row. Nullity = 9 - 1 = 8.
    # Systematic nullspace basis: each row has an information bit and tail bit.
    g = np.hstack((np.eye(8,dtype=np.int64),np.ones((8,1),dtype=np.int64)))
    h = np.ones((1,9),dtype=np.int64)
    selector = np.hstack((np.eye(8,dtype=np.int64),np.zeros((8,1),dtype=np.int64)))
    if args.identity_inner:
        g = np.eye(8,dtype=np.int64)
        h = np.zeros((0,8),dtype=np.int64)
        selector = g.copy()
    hx = expand(a,n,g)
    hz_upper = expand(b,n,selector,transpose=True)
    hz_lower = sparse.kron(sparse.eye(n,format='csr',dtype=np.int64),sparse.csr_matrix(h),format='csr')
    hz = sparse.vstack((hz_upper,hz_lower),format='csr').astype(np.int64)
    css = hx @ hz.T
    css_bad = int(np.count_nonzero(css.data % 2))
    if np.any((g @ h.T) % 2) or css_bad:
        raise ValueError('Orthogonality failed; refusing export')
    d = min(int(np.count_nonzero((np.array([(v >> i)&1 for i in range(8)]) @ g)%2)) for v in range(1,256))
    rx, rz = rank_binary(hx), rank_binary(hz)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for name, matrix in [('G2_binary',g),('H2_binary',h),('selector_binary',selector)]:
        np.savetxt(args.output_dir/(name+'.txt'),matrix,fmt='%d')
    for name, matrix in [('HX_binary',hx),('HZ_binary',hz),('HZ_upper_binary',hz_upper),('HZ_lower_binary',hz_lower)]:
        sparse.save_npz(args.output_dir/(name+'.npz'),matrix)
    inputs = {}
    for path in (a_path,b_path):
        content = path.read_bytes()
        (args.output_dir/path.name).write_bytes(content)
        inputs[path.name] = hashlib.sha256(content).hexdigest()
    report = dict(construction='binary companion-block multiplication, not GF256 scalar Kronecker',
                  polynomial='0x11d', primitive_element=2,
                  companion_convention='M(a) is column multiplication; paper [a]=M(a).T',
                  hx_block='M(a) @ G2 mod 2',hz_upper_block='M(b).T @ selector mod 2',
                  inner_n=g.shape[1],inner_k=g.shape[0],inner_distance=d,outer_columns=n,
                  hx_shape=list(hx.shape),hz_shape=list(hz.shape),hx_rank=rx,hz_rank=rz,
                  physical_qubits=hx.shape[1],logical_qubits=hx.shape[1]-rx-rz,
                  quantum_rate=(hx.shape[1]-rx-rz)/hx.shape[1],css_nonzero_entries=css_bad,
                  source_sha256=inputs,fer=None)
    (args.output_dir/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    print('Saved:',args.output_dir)

if __name__ == '__main__':
    main()
