#!/usr/bin/env python3
"""Summarize paired logical FER and X-side fallback outcomes."""
import argparse,json
from pathlib import Path
from scipy.stats import binomtest

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('no_fallback',type=Path);p.add_argument('fallback',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 x=json.loads(a.no_fallback.read_text());y=json.loads(a.fallback.read_text())
 for k in ['frames','p_D','seed','max_iters','z_prior','matrix_manifest']:
  if x[k]!=y[k]:raise ValueError(f'unmatched {k}')
 rescue=regress=both=0
 for u,v in zip(x['records'],y['records']):
  if u['frame']!=v['frame']:raise ValueError('unmatched frame order')
  rescue+=not u['logical_ok'] and v['logical_ok'];regress+=u['logical_ok'] and not v['logical_ok'];both+=not u['logical_ok'] and not v['logical_ok']
 summary={'frames':x['frames'],'p_D':x['p_D'],'f_m':x['f_m'],'seed':x['seed'],'max_iters':x['max_iters'],
  'no_fallback':{k:x[k] for k in ['logical_failures','logical_fer','x_syndrome_failures','z_syndrome_failures','x_logical_failures','z_logical_failures','undetected_logical_failures','exact_x_frames','conditioned_frames','elapsed_seconds']},
  'binary_HZ_BP_fallback':{k:y[k] for k in ['logical_failures','logical_fer','x_syndrome_failures','z_syndrome_failures','x_logical_failures','z_logical_failures','undetected_logical_failures','exact_x_frames','conditioned_frames','x_fallback_attempts','x_fallback_successes','elapsed_seconds']},
  'rescued':rescue,'regressed':regress,'both_failed':both,
  'mcnemar_exact_two_sided_p':float(binomtest(rescue,rescue+regress).pvalue) if rescue+regress else 1.}
 a.output.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
