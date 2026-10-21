#!/usr/bin/env python3
"""Replace the last cap using a five-direction coupled profile at its old target."""
import json,itertools,numpy as np
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013_orbit_exact import OUT as OLD
from fresh_h20_coupled_completion_20261013 import read,num
Q,R,h,ids,Qe,he,hc,roots,f,g=setup();j=read(OLD);p=np.array([num(f,x)for x in map(f.parse,j['final_hyperideal_point']['vector'])]);p=-p/(p@Q@R[18]);R=R[:175];ns=np.einsum('ij,jk,ik->i',R,Q,R);leaves=np.array([R[18],R[30]/3,R[46]/3,R[1],R[10]]);values=[0.,1.,np.sqrt(2),np.sqrt(3),2.];cosines=np.array([0]+[np.cos(np.pi/m)for m in range(3,121)]);out=[];best=[];pq=p@Q@p;seen=set()
for tags in itertools.product(range(5),repeat=5):
 a=np.array([values[i]for i in tags]);u=a@leaves;bb=2*p@Q@u;cc=u@Q@u-2;D=bb*bb-4*pq*cc
 if D<-1e-8:continue
 for branch,tt in enumerate([(-bb+sg*np.sqrt(max(0,D)))/(2*pq)for sg in[-1,1]]):
  r=u+tt*p
  if r@Q@h>=0 or r@Q@p<=1e-7:continue
  z=-R@Q@r/np.sqrt(2*ns);good=(z>=1-1e-7)|(np.min(abs(z[:,None]-cosines),axis=1)<1e-7);bad=np.flatnonzero(~good).tolist();row={'tags':tags,'coefficients':a.tolist(),'polar_coefficient':tt,'branch':branch,'normal':r.tolist(),'bad':bad,'bad_ratios':[float(z[i])for i in bad],'finite_neighbors':np.flatnonzero(z<1-1e-7).tolist(),'parallel_neighbors':np.flatnonzero(abs(z-1)<1e-7).tolist()}
  if len(bad)<=2:best.append(row)
  if not bad:
   key=tuple(np.round(r,8))
   if key not in seen:seen.add(key);out.append(row)
out.sort(key=lambda x:(-len(x['finite_neighbors']),len(x['parallel_neighbors'])));print('accepted',len(out),'top',[(x['tags'],len(x['finite_neighbors']),x['parallel_neighbors'])for x in out[:20]],flush=True)
(ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_replacement_scout.json').write_text(json.dumps({'status':'NUMERICAL_FIVE_DIRECTION_REPLACEMENT_SCOUT','source':OLD,'tested_profiles':3125,'caps':out,'near_misses':best},indent=2)+'\n')
