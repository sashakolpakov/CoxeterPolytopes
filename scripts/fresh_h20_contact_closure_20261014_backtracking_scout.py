#!/usr/bin/env python3
"""Bounded changed-profile proposals at the new pi/3 finite-edge endpoints.
Old60 retained, sixteen bridge caps optional; numerical proposals only.
"""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.linalg import null_space
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014'
def run():
 D=json.loads((ROOT/f'results/{TAG}.json').read_text());val=lambda x:float(F(x['a']))+np.sqrt(2)*float(F(x['b']));Q=np.array([[val(x)for x in r]for r in D['form']]);N=np.array([[val(x)for x in r]for r in D['normals']]);H=np.array([val(x)for x in D['strict_center']]);E=np.array([val(x)for x in D['strict_escape']['point']]);G=N@Q@N.T;norms=np.diag(G);caps=[];records=[];trials=0;supports=set();profilecount=0
 for vertex in D['finite_vertices']:
  for released in vertex['incident']:
   S=tuple(sorted(set(vertex['incident'])-{released}))
   if S in supports:continue
   supports.add(S);B=G[np.ix_(S,S)];I=np.linalg.inv(B);X=null_space(N[list(S)]@Q);Q2=X.T@Q@X;others=[i for i in range(76)if i not in S];anchors=[i for i in others if norms[i]<3];goodids=[j for j,i in enumerate(S)if abs(norms[i]-2)<1e-6];profiles=[[]]
   for n in[1,2,3]:
    for J in combinations(goodids,n):
     c=np.zeros(19);c[list(J)]=-1;rho=2-c@I@c
     if rho>1e-7:profiles.append(list(J))
   profilecount+=len(profiles)
   for J in profiles:
    c=np.zeros(19);c[J]=-1;P=N[list(S)].T@I@c;rho=2-c@I@c
    for a in anchors:
     l=X.T@Q@N[a];ll=l@l
     if ll<1e-15:continue
     tangent=np.array([l[1],-l[0]]);qa=tangent@Q2@tangent
     for av in[0.,-1.]:
      z0=l*(av-P@Q@N[a])/ll;qb=2*tangent@Q2@z0;qc=z0@Q2@z0-rho
      if abs(qa)<1e-12:roots=[]if abs(qb)<1e-12 else[-qc/qb]
      else:
       dis=qb*qb-4*qa*qc
       if dis<0:continue
       roots=[(-qb+np.sqrt(dis))/(2*qa),(-qb-np.sqrt(dis))/(2*qa)]
      for t in roots:
       r=P+X@(z0+t*tangent);trials+=1
       if r@Q@H>=-1e-8 or r@Q@E<=1e-7:continue
       pp=r@Q@N.T/np.sqrt(2*norms);bad=(pp>1e-6)|((pp>-1+1e-6)&(abs(pp)>1e-6)&(abs(pp+.5)>1e-6)&(abs(pp+2**-.5)>1e-6)&(abs(pp+np.sqrt(3)/2)>1e-6))
       if np.any(bad[:60]):continue
       if any(np.max(abs(r-x))<1e-5 for x in caps):continue
       caps.append(r);records.append(dict(support=S,profile=[S[j]for j in J],anchor=a,anchor_pairing=av,normal=r.tolist(),cut=float(r@Q@E),max_norm_error=float(abs(r@Q@r-2)),conflicting_optional_walls=np.flatnonzero(bad[60:]).__add__(60).tolist()))
 print('supports',len(supports),'profiles',profilecount,'trials',trials,'gapcutters',len(caps),flush=True)
 return dict(status='BOUNDED_NUMERICAL_CHANGED_PROFILE_SCOUT_ONLY',retained_old_walls=60,optional_bridge_walls=16,supports=len(supports),profiles=profilecount,quadratic_branches=trials,candidates=records,scope='Only outgoing supports of16 specified vertices; profiles right except at most3 pi3 contacts to norm2 walls, positive Schur residual; one extra anchor with pairing0 or−1 and oldnorm1/2. Complete all-label/all-parameter search not claimed. Every survivor still needs independent exact reconstruction and global volume testing.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}_backtracking_scout.json').write_text(json.dumps(d,indent=2)+'\n')
