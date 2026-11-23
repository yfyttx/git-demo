"""Selected-failure mechanism replay only; never a FER estimator."""
import concurrent.futures,csv,hashlib,json,os,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
AUDIT=ROOT/'experiments/gf256_cascade/results/z_mechanism_audit_20261116'
STUDY=ROOT/'experiments/gf256_cascade/results/matrix20261019_tuned_vs_syndrome_flip_holdout1000_pD0p07125_seed20353000'
cfg=json.loads((STUDY/'config.json').read_text());exe=Path('/tmp/ldpc-h2-build/gf256_cascade_simulate');matrix=ROOT/cfg['matrix_dir']/'build'
for name,digest in cfg['matrix_sha256'].items():assert hashlib.sha256((matrix/name).read_bytes()).hexdigest()==digest
selected=json.loads((AUDIT/'mechanism_summary.json').read_text())['selected'];selected=[x for x in selected if x['group']=='persistent_z']
out=AUDIT/'budget_diagnostic';out.mkdir(exist_ok=True)
settings=[('bits16_weight1',16,1),('bits8_weight2',8,2)]
(out/'config.json').write_text(json.dumps({'selected':selected,'settings':settings,'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'matrix_sha256':cfg['matrix_sha256']},indent=2)+'\n')
def run(job):
 item,(label,bits,weight)=job;stem=out/f"{label}_seed{item['seed']}_trial{item['trial']}"
 env=os.environ.copy();env.update({'GF256_CASCADE_SKIP_TARGET_DIAG':'1','GF256_CASCADE_SF_BITS':str(bits),'GF256_CASCADE_SF_MAX_WEIGHT':str(weight),'GF256_FRAME_RESULTS':'1'})
 cmd=[str(exe),'5',str(item['seed']),'.07125','conditioned_robust_tuned_syndrome_flip',str(item['trial']),'8',str(matrix),'12']
 completed=all(path.exists() and path.stat().st_size>0 for path in
               (stem.with_suffix('.tsv'),stem.with_suffix('.log')))
 if not completed:
  with stem.with_suffix('.tsv').open('w') as stdout,stem.with_suffix('.log').open('w') as stderr:subprocess.run(cmd,env=env,cwd=ROOT,stdout=stdout,stderr=stderr,check=True)
 log=stem.with_suffix('.log').read_text();s=re.search(r'\[syndrome-flip\] trial=\d+ attempts=(\d+).*?valid_candidates=(\d+)',log);r=re.search(r'\[frame-result\] trial=\d+ logical_failure=(\d+) x_syndrome_ok=(\d+) z_syndrome_ok=(\d+) exact_x=(\d+)',log)
 if r: failure,x_syn,z_syn,exact_x=map(int,r.groups())
 else:
  f=re.search(r'\[failure\] trial=\d+ Xstage=\d+ Xexact=(\d+) Xsyn=(\d+) Zsyn=(\d+)',log)
  failure,x_syn,z_syn,exact_x=(1,int(f[2]),int(f[3]),int(f[1])) if f else (0,1,1,1)
 return {'seed':item['seed'],'trial':item['trial'],'setting':label,'attempts':int(s[1]),'valid_candidates':int(s[2]),'logical_failure':failure,'x_syndrome_ok':x_syn,'z_syndrome_ok':z_syn,'exact_x':exact_x}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,[(x,z) for x in selected for z in settings]))
with (out/'selected_replays.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
