#!/usr/bin/env python3
"""Z-only syndrome BP for binary companion cascade; no true errors enter BP."""
import argparse, json, time, hashlib
from pathlib import Path
import numpy as np
from scipy import sparse
from construct_companion_cascade import ROOT, read_outer, mul

def fwht(x):
    y=x.copy()
    width=1
    while width<256:
        v=y.reshape(*y.shape[:-1],-1,2*width)
        a=v[...,:width].copy(); b=v[...,width:].copy()
        v[...,:width]=a+b; v[...,width:]=a-b
        width*=2
    return y

def row_basis(h):
    piv={}
    for row in range(h.shape[0]):
        v=sum(1<<int(c) for c in h.indices[h.indptr[row]:h.indptr[row+1]])
        while v:
            k=v.bit_length()-1
            if k not in piv:
                piv[k]=v; break
            v^=piv[k]
    return piv

def member(bits,piv):
    v=sum(1<<int(c) for c in np.flatnonzero(bits))
    while v:
        k=v.bit_length()-1
        if k not in piv:return False
        v^=piv[k]
    return True

class OuterBP:
    def __init__(self, rows,n):
        if len(set(map(len,rows)))!=1: raise ValueError('Expected regular check degree')
        self.cols=np.array([[c for c,a in row] for row in rows])
        self.h=np.array([[a for c,a in row] for row in rows])
        self.n=n
        self.products=np.array([[[mul(int(a),x) for x in range(256)] for a in row] for row in self.h])
        self.invperm=np.argsort(self.products,axis=-1)
    def syndrome(self,u):
        vals=np.take_along_axis(self.products,u[self.cols][...,None],axis=-1)[...,0]
        return np.bitwise_xor.reduce(vals,axis=1)
    def decode(self,syn,prior,iters,*,return_log_posterior=False):
        prior = np.asarray(prior, dtype=float)
        if prior.shape not in ((256,), (self.n, 256)):
            raise ValueError('Prior must have shape (256,) or (n,256)')
        if not np.all(np.isfinite(prior)) or np.any(prior < 0) or not np.allclose(prior.sum(axis=-1), 1):
            raise ValueError('Prior must be finite, nonnegative and normalized')
        logprior=np.log(np.maximum(prior,1e-300))
        msg=(np.broadcast_to(prior,(*self.cols.shape,256)).copy()
             if prior.ndim == 1 else prior[self.cols].copy())
        for it in range(1,iters+1):
            weighted=np.take_along_axis(msg,self.invperm,axis=-1)
            ft=fwht(weighted)
            pre=np.concatenate((np.ones_like(ft[:,:1]),np.cumprod(ft[:,:-1],axis=1)),axis=1)
            suf=np.concatenate((np.cumprod(ft[:, :0:-1],axis=1)[:,::-1],np.ones_like(ft[:,:1])),axis=1)
            conv=fwht(pre*suf)/256
            lookup=self.products^syn[:,None,None]
            replies=np.maximum(np.take_along_axis(conv,lookup,axis=-1),1e-300)
            logs=np.log(replies)
            post=np.broadcast_to(logprior,(self.n,256)).copy()
            np.add.at(post,self.cols.ravel(),logs.reshape(-1,256))
            est=np.argmax(post,axis=1)
            if np.array_equal(self.syndrome(est),syn):
                return (est,it,True,post) if return_log_posterior else (est,it,True)
            outgoing=post[self.cols]-logs
            outgoing-=outgoing.max(axis=-1,keepdims=True)
            msg=np.exp(outgoing);msg/=msg.sum(axis=-1,keepdims=True)
        return (est,it,False,post) if return_log_posterior else (est,it,False)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix-dir',type=Path,default=ROOT/'results/companion_binary_n9_k8_matrix20261019')
    p.add_argument('--frames',type=int,default=20)
    p.add_argument('--p-depolar',type=float,default=.07125)
    p.add_argument('--seed',type=int,default=20261117)
    p.add_argument('--max-iters',type=int,default=100)
    p.add_argument('--output',type=Path,default=ROOT/'results/companion_z_bp_pilot20/summary.json')
    args=p.parse_args()
    if args.frames<1 or args.max_iters<1 or not 0<args.p_depolar<1:raise ValueError('Invalid parameters')
    n,rows=read_outer(args.matrix_dir/'HZ_TIT_Outer_For_HX_GF256.txt')
    hx=sparse.load_npz(args.matrix_dir/'HX_binary.npz');hz=sparse.load_npz(args.matrix_dir/'HZ_binary.npz')
    g=np.loadtxt(args.matrix_dir/'G2_binary.txt',dtype=int)
    if not np.array_equal(g,np.hstack((np.eye(8,dtype=int),np.ones((8,1),dtype=int)))):raise ValueError('Requires systematic [9,8] G2')
    basis=row_basis(hz); bp=OuterBP(rows,n)
    f=2*args.p_depolar/3
    w=np.array([x.bit_count() for x in range(256)])
    prior=(1-f)*f**w*(1-f)**(8-w)+f*f**(8-w)*(1-f)**w
    if not np.isclose(prior.sum(),1):raise ValueError('Prior normalization')
    rng=np.random.default_rng(args.seed); records=[];start=time.perf_counter()
    for frame in range(args.frames):
        # Correct marginal of depolarizing noise. No X decoding/conditioning in this experiment.
        true_z=(rng.random((n,9))<f).astype(np.int64)
        synbits=np.asarray(hx@true_z.ravel()).ravel()%2
        syn=synbits.reshape(-1,8)@ (1<<np.arange(8))
        # Simulation-only identity audit, not decoder input.
        true_u=((true_z@g.T)%2)@(1<<np.arange(8))
        if not np.array_equal(bp.syndrome(true_u),syn):raise ValueError('Block mapping mismatch')
        t=time.perf_counter();est,it,valid=bp.decode(syn,prior,args.max_iters)
        recovered=np.zeros((n,9),dtype=np.int64)
        recovered[:,:8]=(est[:,None]>>np.arange(8))&1
        physical_valid=np.array_equal((hx@recovered.ravel())%2,synbits)
        if physical_valid!=valid:raise ValueError('Recovery syndrome mismatch')
        logical_ok=physical_valid and member(true_z.ravel()^recovered.ravel(),basis)
        record=dict(frame=frame,iterations=it,syndrome_ok=bool(valid),z_logical_ok=bool(logical_ok),
                    exact_outer_symbols=bool(np.array_equal(est,true_u)),seconds=time.perf_counter()-t)
        records.append(record)
        print(f'{frame+1}/{args.frames}: Z syndrome={valid}, Z logical={logical_ok}, iterations={it}',flush=True)
    summary=dict(scope='Z component only; no X-side decoding or conditioning; not full quantum FER',
                 frames=args.frames,p_D=args.p_depolar,f_m=f,seed=args.seed,max_iters=args.max_iters,
                 physical_qubits=hx.shape[1],prior='exact G2 pushforward of independent Bernoulli(f_m) bits',
                 z_logical_failures=sum(not r['z_logical_ok'] for r in records),
                 z_syndrome_failures=sum(not r['syndrome_ok'] for r in records),
                 z_undetected_logical_failures=sum(r['syndrome_ok'] and not r['z_logical_ok'] for r in records),
                 exact_outer_frames=sum(r['exact_outer_symbols'] for r in records),
                 elapsed_seconds=time.perf_counter()-start,
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 matrix_manifest=json.loads((args.matrix_dir/'manifest.json').read_text()),records=records)
    summary['z_logical_fer']=summary['z_logical_failures']/args.frames
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('records','matrix_manifest')},indent=2))
if __name__=='__main__':main()
