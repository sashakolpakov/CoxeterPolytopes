#!/usr/bin/env python3
"""Global halfspace geometry after conditioning on a new branch."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_global_h5_20261013_scout as S
def compatible(pair,tol=2e-4):
 mask=pair<-2+tol;ii=np.flatnonzero((pair>=-2+tol)&(pair<=tol));v=pair[ii];theta=np.arccos(np.clip(-v/2,-1,1));m=np.maximum(2,np.rint(np.pi/np.maximum(theta,1e-10)));mask[ii]=abs(v+2*np.cos(np.pi/m))<=tol;return mask
def run():
 R,G,c,names,d,p,C,NR,refs,sources,counts=S.load();bank=json.loads((S.P/'data/enumerations'/f'{S.TAG}_opposite_scout.json').read_text());rows=bank['candidates'];new=np.array([r['normal']for r in rows]);W=np.vstack([R,NR,new]);cov=W@G;A=cov@C;b=cov@p;nonzero=np.linalg.norm(A,axis=1)>1e-7;H=np.c_[A/np.maximum(-b,1e-100)[:,None],-np.ones(len(b))];hard=np.array([i for i,n in enumerate(names)if n.startswith('seed')or n in('compact1','compact2','compact3')]);checks=[];prior=json.loads((S.P/'data/enumerations'/f'{S.TAG}_backtrack.json').read_text())
 for i,r in enumerate(new):
  ps=cov@r;mask=compatible(ps);mask[len(R)+len(NR)+i]=True;assert all(mask[hard]);idx=np.flatnonzero(mask&nonzero)
  oldcase=next(x for x in prior['coherent_cases']if x['index']==i);directions=[C.T@G@np.array(oldcase['leading_point'])]+list(np.eye(5))+list(-np.eye(5));found=None;trials=[]
  for direction in directions:
   ans=linprog(-direction,A_ub=H[idx,:5],b_ub=np.ones(len(idx)),bounds=[(-10,10)]*5,method='highs');assert ans.success,(i,ans.status);z=ans.x;q=float(z@z-1);trials.append(q)
   if q>1e-4:
    z=.999*z;residual=float(max(H[idx,:5]@z-1));assert residual<-1e-6;found={'norm':float(z@z-1),'slice_coordinates':z.tolist(),'point':(p+C@z).tolist(),'maximum_normalized_residual':residual};break
  checks.append({'candidate':i,'allowed_halfspaces':len(idx),'strict_spacelike_witness':found,'lp_trials':trials})
  print(i,len(idx),'escape',None if found is None else found['norm'],flush=True)
 checks.sort(key=lambda r:r['candidate']);out={'status':'NUMERICAL_OPPOSITE_BRANCH_CONDITIONED_GLOBAL_RELAXATIONS','bank_size':len(W),'new_bank_candidates':len(new),'checks':checks,'scope':'Each numerical superset contains all walls judged Coxeter-compatible with its required opposite cap; all remaining mutual conflicts are ignored. All integer labels are allowed numerically through rounded inverse cosine, and parallel pairs remain. No exact graph or obstruction claim.'}
 (S.P/'data/enumerations'/f'{S.TAG}_conditioned.json').write_text(json.dumps(out,indent=2)+'\n');print('escaped',sum(r['strict_spacelike_witness']is not None for r in checks),'of',len(checks),flush=True)
if __name__=='__main__':run()
