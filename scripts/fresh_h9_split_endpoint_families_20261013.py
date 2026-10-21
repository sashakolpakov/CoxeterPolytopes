#!/usr/bin/env python3
"""Numerical discovery in continuous finite-endpoint cap families.

The seed pairings at its eight orthogonal walls have sum of squares <4.
Only the explicitly listed finite labels are sampled; no global exclusion.
"""
import argparse,itertools,json,os
from pathlib import Path
import numpy as np
from scripts.fresh_h9_split_leaf_caps_20261013 import seed,PHI
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_split_endpoint_families_20261013_scout.json'

def profiles(limit):
 values=[2*np.cos(np.pi/m) for m in range(3,limit+1)]
 yield np.zeros(8),[]
 for i in range(8):
  for m,x in zip(range(3,limit+1),values):
   p=np.zeros(8);p[i]=x;yield p,[(i,m)]
 for i,j in itertools.combinations(range(8),2):
  for mi,mj in [(3,3),(3,4),(4,3),(3,5),(5,3)]:
   p=np.zeros(8);p[i]=2*np.cos(np.pi/mi);p[j]=2*np.cos(np.pi/mj);yield p,[(i,mi),(j,mj)]
 for inds in itertools.combinations(range(8),3):
  p=np.zeros(8);p[list(inds)]=1;yield p,[(i,3)for i in inds]

def run(limit):
 Path('/tmp/fresh_h9_split_endpoint_families_20261013.pid').write_text(str(os.getpid()))
 a=PHI+.5;b=-a/(12*a+1);G=seed(a,b);eps=np.r_[np.ones(4),-np.ones(4)];S=np.ones(16)
 H=128*a*(1+2*a);L=-128*b*(1+2*b);w=np.array([8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]);gw=G@w
 vals=np.array([0]+[2*np.cos(np.pi/m)for m in range(3,limit+1)])
 rows=[];seen=set();maximum=(-1e100,None);counts={'profiles_per_edge':0,'tested_parameters':0,'accepted_discrete':0,'accepted_continuous_boundaries':0}
 def keep(side,sg,x,labels,t,kind,q):
  if t<=0:return
  y=(np.sqrt(H)*t-sg*eps*np.sqrt(L*(q+t*t)))/8-x-(a*sum(x)+b*np.dot(x,eps)*eps)/2
  if np.min(y)<-1e-7:return
  if kind=='finite':
   err=np.min(abs(y[:,None]-vals[None,:]),axis=1)
   if any((y<2-1e-7)&(err>1e-7)):return
  P=S.copy();P[side*8:side*8+8]+=4*a
  R=np.tile(eps,2);R[side*8:side*8+8]+=4*b*eps
  r=t/np.sqrt(H)*P+sg*np.sqrt((q+t*t)/L)*R;r[side*8:side*8+8]-=x/2
  cov=G@r;assert abs(r@cov-2)<1e-6 and np.max(cov)<1e-6 and sum(cov)<0
  key=tuple(np.round(cov,7))
  if key in seen:return
  seen.add(key);pair=float(r@gw)
  row={'seed_side':side,'endpoint_sign':sg,'finite_profile':labels,'parameter_sqrt_c':float(t),'projection_norm':float(q),'kind':kind,'normal_coefficients':r.tolist(),'seed_covector':cov.tolist(),'old_universal_escape_pairing':pair}
  rows.append(row);counts['accepted_discrete' if kind=='finite' else 'accepted_continuous_boundaries']+=1
  nonlocal maximum
  if pair>maximum[0]:maximum=(pair,len(rows)-1)
 for x,labels in profiles(limit):
  counts['profiles_per_edge']+=1;q=2-np.dot(x,x)/2;assert q>0
  for side in range(2):
   for sg in [-1,1]:
    targets=[]
    for j in range(8):
     B=sg*eps[j]*np.sqrt(L);K=H-L
     for v in vals:
      D=8*(v+x[j]+(a*sum(x)+b*np.dot(x,eps)*eps[j])/2)
      t=(np.sqrt(H)*D+B*np.sqrt(D*D+q*K))/K
      targets.append(float(t))
    for t in sorted(set(round(t,12)for t in targets)):
     counts['tested_parameters']+=1;keep(side,sg,x,labels,t,'finite',q)
    boundaries=[]
    for j in range(8):
     B=sg*eps[j]*np.sqrt(L);K=H-L;D=8*(2+x[j]+(a*sum(x)+b*np.dot(x,eps)*eps[j])/2)
     boundaries.append((np.sqrt(H)*D+B*np.sqrt(D*D+q*K))/K)
    keep(side,sg,x,labels,max(boundaries),'all_dotted_boundary',q)
 print({'counts':counts,'unique_caps':len(rows),'positive_old_escape_cutters':sum(r['old_universal_escape_pairing']>1e-7 for r in rows),'max_pairing':maximum[0]},flush=True)
 out={'status':'NUMERICAL_ENDPOINT_FAMILY_SCOUT','scope':'Finite endpoint profiles use labels3..'+str(limit)+'; discrete seed-compatible heights and one all-dotted threshold per continuous family are tested. New cap pairs/volume and exactness require separate checks. Thresholds may include parallel seed pairs.','counts':counts,'max_old_escape_pairing':maximum[0],'best_row':maximum[1],'caps':rows}
 OUT.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-label',type=int,default=30);a=p.parse_args();run(a.max_label)
