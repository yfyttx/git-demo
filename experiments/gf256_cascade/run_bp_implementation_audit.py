#!/usr/bin/env python3
"""Controlled C++/Python BP audit and fixed-kernel G2 channel comparison."""
import argparse,json,subprocess,time,hashlib
from pathlib import Path
import numpy as np
from construct_companion_cascade import ROOT, read_outer,mul
from run_companion_z_bp import OuterBP
from run_identity_outer_comparison import trace

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bridge',type=Path,default=Path('/tmp/git-demo-kasai-direct-build/audit_bp_bridge'))
p.add_argument('--frames',type=int,default=100)
p.add_argument('--cross-frames',type=int,default=20)
p.add_argument('--max-iters',type=int,default=100)
p.add_argument('--fm',type=float,default=.0475)
p.add_argument('--seed',type=int,default=20261121)
p.add_argument('--output',type=Path,default=ROOT/'results/bp_implementation_audit_20261117')
a=p.parse_args()
if not 0<a.fm<.5 or not 0<a.cross_frames<=a.frames:raise ValueError('invalid parameters')
a.output.mkdir(parents=True,exist_ok=True)
d=ROOT/'results/companion_binary_n9_k8_matrix20261019'
f=a.fm
weights=np.array([v.bit_count() for v in range(256)])
raw=f**weights*(1-f)**(8-weights)
projected=(1-f)*raw+f*raw[::-1]
dual=np.array([sum(trace(mul(v,1<<i))<<i for i in range(8)) for v in range(256)])
inv=np.argsort(dual)
start=time.perf_counter();allresults={}
for name,filename,isdual,nframes in [
 ('Delta_polynomial','HZ_TIT_Outer_For_HX_GF256.txt',False,a.frames),
 ('Gamma_polynomial','HX_TIT_Outer_For_HZ_GF256.txt',False,a.cross_frames),
 ('Delta_dual','HZ_TIT_Outer_For_HX_GF256.txt',True,a.cross_frames)]:
 n,rows=read_outer(d/filename);bp=OuterBP(rows,n)
 rng=np.random.default_rng(a.seed)
 cases=[];truth=[];priors=[];labels=[]
 for frame in range(nframes):
  noise=(rng.random((n,9))<f).astype(int)
  direct=noise[:,:8]@(1<<np.arange(8))
  if isdual:direct=inv[direct]
  variants=[('original',direct,raw[dual] if isdual else raw,1)]
  if name=='Delta_polynomial':
   u=(noise[:,:8]^noise[:,8,None])@(1<<np.arange(8))
   variants.extend([('G2_exact_prior',u,projected,0),('G2_raw_prior_control',u,raw,0)])
  for mode,t,prior,legacy in variants:
   cases.append((bp.syndrome(t),legacy));truth.append(t);priors.append(prior);labels.append((frame,mode))
 inp=a.output/(name+'.input.txt');out=a.output/(name+'.cpp.txt')
 with inp.open('w') as stream:
  stream.write(f'{len(rows)} {n} {sum(map(len,rows))} {len(cases)} {a.max_iters}\n')
  for r,row in enumerate(rows):
   for c,v in row:stream.write(f'{r} {c} {v}\n')
  for (syn,legacy),prior in zip(cases,priors):
   stream.write(f'{f:.17g} {int(isdual)} {legacy}\n')
   stream.write(' '.join(format(v,'.17g') for v in prior)+'\n')
   stream.write(' '.join(map(str,syn))+'\n')
 np.savez_compressed(a.output/(name+'.truth.npz'),truth=np.array(truth),frame_mode=np.array(labels,dtype=str))
 subprocess.run([str(a.bridge),str(inp),str(out)],check=True)
 cpp=np.loadtxt(out,dtype=int,ndmin=2);records=[]
 for k,((syn,_),t,prior,(frame,mode)) in enumerate(zip(cases,truth,priors,labels)):
  ci,cv=cpp[k,:2];ce=cpp[k,2:]
  rec=dict(frame=frame,mode=mode,cpp_iterations=int(ci),cpp_syndrome_ok=bool(cv),cpp_block_error=bool(np.any(ce!=t)))
  if frame<a.cross_frames:
   pe,pi,pv=bp.decode(syn,prior,a.max_iters)
   rec.update(python_iterations=pi,python_syndrome_ok=bool(pv),python_block_error=bool(np.any(pe!=t)),
              estimate_equal=bool(np.array_equal(ce,pe)),iterations_equal=bool(ci==pi),syndrome_status_equal=bool(cv==pv),different_symbols=int(np.count_nonzero(ce!=pe)))
  records.append(rec)
  if (k+1)%20==0:print(name,k+1,'/',len(cases),flush=True)
 modes={}
 for mode in sorted(set(r['mode'] for r in records)):
  rs=[r for r in records if r['mode']==mode];cross=[r for r in rs if 'estimate_equal' in r]
  modes[mode]=dict(frames=len(rs),cpp_block_errors=sum(r['cpp_block_error'] for r in rs),cpp_syndrome_failures=sum(not r['cpp_syndrome_ok'] for r in rs),
                  cross_frames=len(cross),estimate_mismatches=sum(not r['estimate_equal'] for r in cross),iteration_mismatches=sum(not r['iterations_equal'] for r in cross),status_mismatches=sum(not r['syndrome_status_equal'] for r in cross),python_block_errors=sum(r['python_block_error'] for r in cross))
 allresults[name]=dict(modes=modes,records=records)
 (a.output/(name+'.json')).write_text(json.dumps(allresults[name],indent=2)+'\n')
 print(name,json.dumps(modes),flush=True)
summary=dict(f_m=f,p_D=1.5*f,seed=a.seed,max_iters=a.max_iters,elapsed_seconds=time.perf_counter()-start,
             results={name:r['modes'] for name,r in allresults.items()},
             note='G2_raw_prior_control deliberately mismatches the projected channel; not proposed decoder')
summary['source_sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),ROOT/'run_companion_z_bp.py',ROOT/'../../src/decoder/NB_BP_Decoder.cpp',ROOT/'audit_bp_bridge.cpp']}
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
