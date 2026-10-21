#!/usr/bin/env python3
"""Bounded floating pair-screen; exact selected pair is bound separately."""
import argparse,gzip,json
from hashlib import sha256
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_pair_20261015_pairs'
def run():
    sources=[f'scripts/{TAG}.py','results/fresh_h9_joint_pair_20261015_twocaps.json'];rows=[]
    for suffix in ('escape_right','escape'):
        path=f'results/fresh_h9_joint_pair_20261015_scout_{suffix}.json';d=json.loads((ROOT/path).read_text());bank=d['candidate_bank']['path'];blob=(ROOT/bank).read_bytes()
        assert sha256(blob).hexdigest()==d['candidate_bank']['sha256'];b=json.loads(gzip.decompress(blob));assert len(b)==d['candidate_count'];sources.extend([path,bank])
        N=np.array([r['physical_normal']for r in b]);metric=np.diag([-1.]+[1.]*9);values=-2*np.cos(np.pi/np.arange(2,257));pairs=[];dotted=0
        for i in range(len(b)):
            vals=(N[i]@metric)@N[i+1:].T;dotted+=int(sum(vals<-2-1e-7))
            for k in np.flatnonzero((vals>-2+1e-6)&(vals<1e-7)):
                m=int(np.argmin(abs(vals[k]-values)))+2;error=abs(vals[k]-values[m-2]);j=i+1+int(k)
                if error>=1e-7:continue
                labels=[m]+[x[1]for x in b[i]['anchors']+b[j]['anchors']+b[i]['support_contacts']+b[j]['support_contacts']]
                pairs.append(dict(ids=[i,j],label=m,error=float(error),pairing=float(vals[k]),max_label=max(labels)))
        pairs.sort(key=lambda r:(r['max_label'],r['error']))
        rows.append(dict(scout=suffix,bank_count=len(b),tested_pairs=len(b)*(len(b)-1)//2,dotted_pair_events=dotted,fitted_finite_pairs=pairs,
                         low_label_recipe_records=[dict(pair=r,recipes=[b[i]for i in r['ids']])for r in pairs if r['max_label']<=5]))
    return dict(status='BOUNDED_FLOATING_MUTUAL_PAIR_SCREEN_WITH_SEPARATE_EXACT_PAIR',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},banks=rows,
        scope='Numericalscreen only. Rounded proposal bank may contain geometric duplicates. Finite anglefits and nearparallel labels can be falsepositives. Only the separately bound two-cap certificate establishes an exact pi3 pair. No exhaustive joint-completion classification.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],[(r['bank_count'],len(r['fitted_finite_pairs']))for r in d['banks']])
