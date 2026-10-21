#!/usr/bin/env python3
"""Bounded nonincident angle events for the14exact affine local profiles."""
import argparse,json,math
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_h9_h8_disjoint_walk_20261014_scout as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_walk_20261014_affine_scout'
def run():
 d,f,profiles,Q,N,T,eta,h,W,G,E=S.data()
 apath='results/fresh_h9_h8_disjoint_walk_20261014_affine.json';a=json.loads((ROOT/apath).read_text())
 labels=[2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64,0]
 pool=-2*np.cos(np.pi/np.arange(2,257));out=[];trials=0
 for row in a['supports']:
  eid=row['germ'];queue=f['queue'][eid];block=row['support'];outside=[i for i in range(15)if i not in block]
  P=S.vec(queue['start_point'])@T.T;D=S.vec(queue['direction'])@T.T
  U=P/math.sqrt(-P@eta@P);V=D/math.sqrt(D@eta@D);Z=U+V
  assert abs(Z@eta@Z)<1e-8 and max(Z@eta@W[outside].T)<-1e-8
  J=np.linalg.inv(G[np.ix_(block,block)])
  for k,prof in enumerate(row['affine_profiles']):
   c=np.zeros(8)
   for i,m in prof['contacts']:c[block.index(i)]=-2*math.cos(math.pi/m)
   p=c@J@W[block];assert abs(p@eta@p-2)<1e-7
   for j in outside:
    for m in labels:
     target=-2*math.cos(math.pi/m)if m else -2
     t=(target-p@eta@W[j])/(Z@eta@W[j]);trials+=1
     if t<=0:continue
     r=p+t*Z;vals=r@eta@W.T
     err=np.min(abs(vals[:,None]-pool[None,:]),axis=1)
     if not np.all((vals<=-2+1e-8)|(err<1e-8)):continue
     assert abs(r@eta@r-2)<1e-6 and r@eta@P<0
     out.append(dict(germ=eid,profile=k,contacts=prof['contacts'],anchor=j,anchor_label=m,null_coefficient=float(t),
        pairings=vals.tolist(),escape_pairing=float(r@eta@E),cuts_escape=bool(r@eta@E>1e-8)))
 files=S.FILES+[apath,'scripts/fresh_h9_h8_disjoint_walk_20261014_scout.py',f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_AFFINE_ENDPOINT_EVENTS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
   anchor_labels=labels,trial_count=trials,accepted_event_count=len(out),escape_cutters=sum(r['cuts_escape']for r in out),events=out,
   scope='Bounded numerical finite/parallel-anchor events for14complete local affine profiles. Label0 denotes parallel; global finite-angle matching uses labels2..256,tolerance1e-8. Continuous all-dotted parameter tails and unsampled anchor labels are not excluded. No exact newcap or finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['accepted_event_count'],d['escape_cutters'],flush=True)
