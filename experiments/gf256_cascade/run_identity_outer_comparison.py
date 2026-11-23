#!/usr/bin/env python3
"""Same-frame identity-inner vs direct outer syndrome comparison (both components)."""
import argparse,json,time
import numpy as np
from scipy import sparse
from construct_companion_cascade import ROOT,read_outer,mul,expand
from run_companion_z_bp import OuterBP

def trace(v):
    t=0
    for _ in range(8):t^=v;v=mul(v,v)
    return t

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frames',type=int,default=100)
    p.add_argument('--fm',type=float,default=.0475)
    p.add_argument('--seed',type=int,default=20261118)
    p.add_argument('--max-iters',type=int,default=50)
    a=p.parse_args()
    if a.frames<1 or not 0<a.fm<2/3 or a.max_iters<1:raise ValueError('Invalid parameters')
    d=ROOT/'results/companion_identity_n8_k8_matrix20261019'
    n,delta=read_outer(d/'HZ_TIT_Outer_For_HX_GF256.txt')
    _,gamma=read_outer(d/'HX_TIT_Outer_For_HZ_GF256.txt')
    hx=sparse.load_npz(d/'HX_binary.npz');hz=sparse.load_npz(d/'HZ_binary.npz')
    # Independently expand outer blocks without inner coding.
    ident=np.eye(8,dtype=np.int64)
    if (hx!=expand(delta,n,ident)).nnz or (hz!=expand(gamma,n,ident,True)).nnz:raise ValueError('Matrix mismatch')
    dual=np.array([sum(trace(mul(v,1<<i))<<i for i in range(8)) for v in range(256)])
    inv=np.argsort(dual)
    weight=np.array([v.bit_count() for v in range(256)])
    prior=a.fm**weight*(1-a.fm)**(8-weight)
    bpz=OuterBP(delta,n);bpx=OuterBP(gamma,n)
    rng=np.random.default_rng(a.seed); records=[];start=time.perf_counter()
    for frame in range(a.frames):
        r=rng.random((n,8));pd=1.5*a.fm
        x=((r<pd/3)|((r>=2*pd/3)&(r<pd))).astype(int)
        z=((r>=pd/3)&(r<pd)).astype(int)
        zs=z@(1<<np.arange(8));xs=inv[x@(1<<np.arange(8))]
        sz=((hx@z.ravel())%2).reshape(-1,8)@(1<<np.arange(8))
        sx=inv[((hz@x.ravel())%2).reshape(-1,8)@(1<<np.arange(8))]
        directz=bpz.syndrome(zs);directx=bpx.syndrome(xs)
        if not np.array_equal(sz,directz) or not np.array_equal(sx,directx):raise ValueError('Syndrome mismatch')
        # Explicitly execute both deterministic paths on the same physical frame.
        zh,zi,zok=bpz.decode(sz,prior,a.max_iters)
        xh,xi,xok=bpx.decode(sx,prior[dual],a.max_iters)
        zd,_,_=bpz.decode(directz,prior,a.max_iters)
        xd,_,_=bpx.decode(directx,prior[dual],a.max_iters)
        if not np.array_equal(zh,zd) or not np.array_equal(xh,xd):raise ValueError('Decoder mismatch')
        records.append(dict(frame=frame,delta_block_error=bool(np.any(zh!=zs)),gamma_block_error=bool(np.any(xh!=xs)),delta_syndrome_failure=not zok,gamma_syndrome_failure=not xok,z_iters=zi,x_iters=xi))
        if (frame+1)%10==0:print(f'{frame+1}/{a.frames}',flush=True)
    result=dict(frames=a.frames,f_m=a.fm,p_D=1.5*a.fm,seed=a.seed,max_iters=a.max_iters,
                convention='current cascade: Delta polynomial/Z, Gamma transposed/dual-X',
                matrix_mismatches=0,syndrome_mismatches=0,decoder_mismatches=0,
                elapsed_seconds=time.perf_counter()-start,records=records)
    for key in ('delta_block_error','gamma_block_error','delta_syndrome_failure','gamma_syndrome_failure'):
        result[key+'s']=sum(r[key] for r in records)
    out=ROOT/f'results/identity_outer_comparison_fm{a.fm}_seed{a.seed}'
    out.mkdir(parents=True,exist_ok=True);(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()
