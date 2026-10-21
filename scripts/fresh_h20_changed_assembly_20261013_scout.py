#!/usr/bin/env python3
"""Bounded numerical proposals for new polar-reflected H20 intersections.

This changes the reflecting hyperplane, rather than reusing the earlier bank.
Any accepted exact bank would be coherent internally because a single Lorentz
reflection transports it. Only the cross-pair gates need to be tested anew.
"""
import argparse,json
import numpy as np
from scipy.optimize import linprog
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013_second_exact import OUT as SOURCE
from fresh_h20_coupled_completion_20261013 import num,read
from fresh_dimension_bridge_two_angle_20261013_independent import B

OUT='data/enumerations/fresh_h20_changed_assembly_20261013_scout.json'

def run(limit):
 Q,R,h,*_=setup();d=read(SOURCE)
 q=np.einsum('ij,jk,ik->i',R,Q,R);sqrtq=np.sqrt(q);G=R@Q@R.T
 cs=np.array([0.]+[np.cos(np.pi/m)for m in range(3,121)])
 hn=h/np.sqrt(-h@Q@h);cov=R@Q/sqrtq[:,None];rng=np.random.default_rng(2026101337)
 polars=[{'kind':kind,'vector':[num(B,B.parse(x))for x in d[kind]['vector']]}for kind in('final_hyperideal_point','strict_escape')]
 seen=set();trials=[];target=np.array(polars[1]['vector'])
 for i in range(limit):
  ans=linprog(rng.normal(size=21),A_ub=cov,b_ub=np.zeros(len(R)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if ans.success and ans.x@Q@ans.x>1e-6:
   active=tuple(np.flatnonzero(abs(cov@ans.x)<1e-7).tolist())
   if active not in seen:seen.add(active);polars.append({'kind':'LP_gap','trial':i,'active_indices':list(active),'vector':ans.x.tolist()})
 for i in range(limit):
  ans=linprog(rng.normal(size=21),A_ub=cov,b_ub=np.zeros(len(R)),A_eq=np.vstack([hn@Q,cov[175]]),b_eq=[-1,0],bounds=[(None,None)]*21,method='highs')
  if ans.success and ans.x@Q@ans.x>1e-6:
   active=tuple(np.flatnonzero(abs(cov@ans.x)<1e-7).tolist())
   if active not in seen:seen.add(active);polars.append({'kind':'LP_last_facet_gap','trial':i,'active_indices':list(active),'vector':ans.x.tolist()})
 for t,polar in enumerate(polars):
  p=np.array(polar['vector']);pq=float(p@Q@p);pr=R@Q@p
  assert pq>0 and max(pr)<1e-6 and p@Q@h<0
  X=R-2*pr[:,None]*p/pq
  cross=G-2*pr[:,None]*pr[None,:]/pq
  z=-cross/np.outer(sqrtq,sqrtq)
  ok=(z>=1-1e-7)|(np.min(abs(z[...,None]-cs),axis=2)<1e-7)
  unchanged=np.linalg.norm(X-R,axis=1)<1e-6
  accepted=np.flatnonzero(np.all(ok,axis=0)&~unchanged).tolist()
  cp=h-(h@Q@p)/pq*p
  assert max(R@Q@cp)<0 and max(X@Q@cp)<0 and cp@Q@cp<0
  U=np.vstack([R,X[accepted]])
  uq=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(uq[:,None]);ch=cp/np.sqrt(-cp@Q@cp)
  escape=None;volume_trials=[]
  for k in range(40):
   res=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(ch@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
   if res.success:
    qq=float(res.x@Q@res.x);volume_trials.append(qq)
    if qq>1e-6:
     escape={'vector':res.x.tolist(),'norm':qq,'active_indices':np.flatnonzero(abs(A@res.x)<1e-7).tolist()};break
  row=dict(polar,norm=pq,accepted_reflected_originals=accepted,unchanged_originals=np.flatnonzero(unchanged).tolist(),
           rejected_originals=np.flatnonzero(~np.all(ok,axis=0)).tolist(),
           rejection_counts=np.sum(~ok,axis=0).tolist(),strict_center=cp.tolist(),
           reflected_target_pairings=(X@Q@target).tolist(),
           accepted_meeting_last_facet=[i for i in accepted if z[175,i]<=1+1e-7],
           old_target_cut_by_accepted=bool(accepted and max(X[accepted]@Q@target)>1e-7),
           volume_trials=volume_trials,escape=escape)
  trials.append(row)
  print(t,polar['kind'],'q',pq,'accepted',len(accepted),'cuts',row['old_target_cut_by_accepted'],'escape',None if escape is None else escape['norm'],flush=True)
 out={'status':'NUMERICAL_NEW_POLAR_REFLECTED_ASSEMBLIES_ONLY','source':SOURCE,'LP_objectives':limit,
      'polars':trials,'scope':'All pair gates and volume searches in this scout are floating proposals; no exact construction or finite-volume conclusion.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=32);a=p.parse_args();run(a.limit)
