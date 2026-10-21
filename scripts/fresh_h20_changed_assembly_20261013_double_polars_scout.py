#!/usr/bin/env python3
"""Bounded rational polar proposals after the doubled-source core replacement."""
import json,numpy as np
from fractions import Fraction as F
from functools import reduce
from math import gcd,lcm
from scipy.optimize import linprog
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013 import read
OUT='data/enumerations/fresh_h20_changed_assembly_20261013_double_polars_scout.json'
Q,old,h,*_=setup();data=read('data/enumerations/fresh_h20_changed_assembly_20261013_doubled_lift_scout.json');U=np.array([[float(F(x))for x in row]for row in data['lift_normals']]);initial=read('data/enumerations/fresh_h20_changed_assembly_20261013_double_volume_scout.json');ids=initial['retained_old'];R=np.vstack([U,old[ids]]);h=np.array(initial['center']);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(2026101341);allowed=np.array([0.]+[np.cos(np.pi/m)for m in range(3,121)]);added=[];rejected=[];stages=[];seen=set()
first=np.array([52,-22,-14,7,8,13,29,36,39,-13,-20,-11,-7,-8,-13,18,22,26,-3,-6,-1])
firstq=float(first@Q@first);firstns=np.einsum('ij,jk,ik->i',R,Q,R);firstz=-R@Q@first/np.sqrt(firstq*firstns)
assert firstq==1 and np.all((firstz>=1-1e-7)|(np.min(abs(firstz[:,None]-allowed),axis=1)<1e-7))
added.append({'stage':-1,'trial':None,'integer_normal':first.tolist(),'integer_norm':firstq,'active_indices':np.flatnonzero(abs(firstz)<1e-7).tolist(),'finite_neighbors':np.flatnonzero(firstz<1-1e-7).tolist(),'parallel_neighbors':np.flatnonzero(abs(firstz-1)<1e-7).tolist(),'bad':[]})
R=np.vstack([R,first])
for stage in range(9):
 ns=np.einsum('ij,jk,ik->i',R,Q,R);cov=R@Q/np.sqrt(ns[:,None]);found=None;escape=None;trials=[]
 for trial in range(100):
  ans=linprog(rng.normal(size=21),A_ub=cov,b_ub=np.zeros(len(R)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if not ans.success:continue
  q=float(ans.x@Q@ans.x);trials.append(q)
  if q<1e-6:continue
  active=tuple(np.flatnonzero(abs(cov@ans.x)<1e-7).tolist())
  row={'stage':stage,'trial':trial,'vector':ans.x.tolist(),'norm':q,'active_indices':list(active)}
  if escape is None:escape=row
  if active in seen:continue
  seen.add(active);vals=ans.x/np.min(np.abs(ans.x[np.abs(ans.x)>1e-6]));f=[F(float(x)).limit_denominator(1000)for x in vals];den=lcm(*(x.denominator for x in f));v=[int(x*den)for x in f];gg=reduce(gcd,v);v=[x//gg for x in v]
  vv=np.array(v,dtype=float);scale=(vv@ans.x)/(ans.x@ans.x)
  if max(abs(vv-scale*ans.x))>1e-6 or max(abs(vv))>10**7:
   rejected.append(dict(row,rejection='nonrational approximation'));continue
  qq=float(vv@Q@vv);z=-R@Q@vv/np.sqrt(ns*qq);bad=np.flatnonzero(~((z>=1-1e-7)|(np.min(abs(z[:,None]-allowed),axis=1)<1e-7))).tolist()
  row.update(integer_normal=v,integer_norm=qq,bad=bad,bad_cosines=[float(z[i])for i in bad],finite_neighbors=np.flatnonzero(z<1-1e-7).tolist(),parallel_neighbors=np.flatnonzero(abs(z-1)<1e-7).tolist())
  if bad:rejected.append(row);continue
  if len(row['finite_neighbors'])<20:
   rejected.append(dict(row,rejection='fewer than20 finite neighbors'));continue
  found=row;break
 stages.append({'stage':stage,'walls':len(R),'trials':len(trials),'largest_trial_square':max(trials,default=None),'escape':escape})
 print('stage',stage,'walls',len(R),'found',None if found is None else(found['integer_norm'],len(found['finite_neighbors']),found['integer_normal']),'escape',None if escape is None else escape['norm'],flush=True)
 if found is None or stage==8:break
 added.append(found);R=np.vstack([R,np.array(found['integer_normal'])])
(ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_CHANGED_CORE_BOUNDED_POLAR_EXTENSIONS','retained_old':ids,'added_polars':added,'rejected_polar_rays':rejected,'stages':stages,'final_wall_count':len(R)},indent=2)+'\n')
