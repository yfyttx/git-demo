"""Audit completed paired logs; select deterministic strata for offline replay.
Does not send truth to the decoder, and does not estimate population FER from
selected replay frames.
"""
import argparse,csv,importlib.util,json,re,hashlib,os,subprocess,concurrent.futures
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('study',Path(__file__).with_name('run_paired_decoder_study.py'))
study=importlib.util.module_from_spec(spec);spec.loader.exec_module(study)
p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--replay',action='store_true');args=p.parse_args()
args.out.mkdir(parents=True,exist_ok=True)
cfg=json.loads((args.study/'config.json').read_text());modes=cfg['modes']; records=[]
for seed in range(cfg['seed_start'],cfg['seed_start']+200):
 data=[];logs=[]
 for mode in modes:
  tsv,log,_=study.task_paths(args.study,.07125,seed,mode)
  data.append(study.read_result(tsv,log));logs.append(log.read_text())
 sf={int(m[0]):{'attempts':int(m[1]),'valid_candidates':int(m[2])} for m in re.findall(r'\[syndrome-flip\] trial=(\d+) attempts=(\d+).*?valid_candidates=(\d+)',logs[1])}
 for t in range(5):
  a=data[0]['failure_types'].get(t,'success');b=data[1]['failure_types'].get(t,'success')
  group=('persistent_z' if a==b=='z_syndrome' else 'rescued' if a=='z_syndrome' and b=='success' else 'undetected' if b=='undetected_logical' else 'success' if a==b=='success' else 'other')
  records.append({'seed':seed,'trial':t,'group':group,'baseline':a,'syndrome_flip':b,**sf.get(t,{'attempts':0,'valid_candidates':0})})
with (args.out/'all_1000_frame_mechanisms.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
summary={g:{'frames':sum(r['group']==g for r in records),'with_valid_sf_candidate':sum(r['group']==g and r['valid_candidates']>0 for r in records),'total_attempts':sum(r['attempts'] for r in records if r['group']==g)} for g in sorted({r['group'] for r in records})}
selected=[]
for g,n in [('persistent_z',3),('rescued',1),('undetected',1),('success',1)]:selected += [r for r in records if r['group']==g][:n]
(args.out/'mechanism_summary.json').write_text(json.dumps({'groups':summary,'selection_rule':'First three persistent Z failures and first one from each other category, sorted seed/trial; diagnostic only','selected':selected},indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
if args.replay:
 exe=Path('/tmp/ldpc-h2-build/gf256_cascade_simulate'); matrices=ROOT/cfg['matrix_dir']/'build'
 for name,digest in cfg['matrix_sha256'].items():assert study.sha256(matrices/name)==digest
 manifest={'executable_sha256':study.sha256(exe),'matrix_sha256':cfg['matrix_sha256'],'selected':selected,'mode':modes[1],'p_D':.07125}
 (args.out/'replay_config.json').write_text(json.dumps(manifest,indent=2)+'\n')
 def run(r):
  stem=args.out/f"{r['group']}_seed{r['seed']}_trial{r['trial']}"
  env=os.environ.copy();env['GF256_CASCADE_SKIP_TARGET_DIAG']='1';env['GF256_CASCADE_PROJECTION_RISK_AUDIT']='1';env['GF256_FRAME_EXPORT_PREFIX']=str(stem.resolve());env['GF256_FRAME_RESULTS']='1'
  cmd=[str(exe),'5',str(r['seed']),'.07125',modes[1],str(r['trial']),'8',str(matrices),'12']
  with stem.with_suffix('.tsv').open('w') as out,stem.with_suffix('.log').open('w') as log:subprocess.run(cmd,env=env,cwd=ROOT,stdout=out,stderr=log,check=True)
  print('replayed',stem.name,flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run,selected))
