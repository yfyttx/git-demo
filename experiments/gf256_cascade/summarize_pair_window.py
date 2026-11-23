"""Summarize one offline joint-physical-window candidate probe."""
import argparse
import csv
import json
import re
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    folder=args.directory
    cfg=json.loads((folder/'config.json').read_text())
    outcomes=list(csv.DictReader((folder/'replay_outcomes.csv').open()))
    rows=list(csv.DictReader((folder/'projection_coverage.csv').open()))
    assert len(outcomes)==cfg['limit'] and all(x['replay_matches']=='True' for x in outcomes)
    report=[]
    for frame in outcomes:
        seed,trial=frame['seed'],frame['trial']
        rr=[r for r in rows if r['seed']==seed and r['trial']==trial]
        pair=[r for r in rr if r['stage'].endswith('_pair')]
        assert len(pair)==11,(seed,trial,len(pair))
        base=[r for r in rr if r['stage']=='initial']
        assert len(base)==11
        log=(folder/f'seed{seed}_trial{trial}.log').read_text()
        meta=re.findall(r'\[pair-candidate-audit\] trial=\d+ stage=(\w+) left=(-?\d+) physical_candidates=(\d+)',log)
        assert len(meta)==1,meta
        _,left,candidates=meta[0];left=int(left)
        changes=[]
        for a in (range(11) if left == -2 else (left,left+1)):
            before=next(r for r in base if int(r['coordinate'])==a)
            after=next(r for r in pair if int(r['coordinate'])==a)
            assert int(after['coset_covered'])>=int(before['coset_covered'])
            changes.append(dict(coordinate=a,base_coset=int(before['coset_covered']),pair_coset=int(after['coset_covered']),
                                base_candidates=int(before['candidates']),pair_candidates=int(after['candidates'])))
        report.append(dict(seed=int(seed),trial=int(trial),left=left,
                           physical_candidates=int(candidates),changes=changes))
    output=dict(frames=report,limitation='Selected failure frames; diagnostic candidate lists only; no online FER estimate')
    (folder/'pair_window_summary.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__': main()
