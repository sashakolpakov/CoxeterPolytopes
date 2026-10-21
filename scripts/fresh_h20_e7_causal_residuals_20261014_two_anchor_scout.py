#!/usr/bin/env python3
"""Bounded mixed-label E7 contacts with two off-face finite anchors."""
import json,math,time
from itertools import combinations,product
from collections import Counter
import numpy as np
import fresh_h20_e7_causal_residuals_20261014_scout as B
def run():
 start=time.time();d,f,S,J,N,V,Q,VG,VP,VE,SP,SE,edges=B.data();A=f['face_facet_walls'][:13];I=np.linalg.inv(N[A]@Q@N[A].T);W=I@N[A]
 val=lambda v:sum(float(B.F(a))*math.sqrt(b)for a,b in zip(v,[1,2,79,158,11,22,869,1738]))
 T=np.array([val(v)for v in f['face_basis'][-1]])/math.sqrt(2);E=np.array([val(v)for v in d['strict_escape']['point']]);frame=np.array([-W[0],-W[12],T]);C=frame@Q@frame.T;C[np.abs(C)<1e-9]=0
 AP=frame@Q@N[[61,75]].T;U,sv,Vh=np.linalg.svd(AP.T);direction=Vh[-1];P=np.linalg.pinv(AP.T)
 labels=[2,3,4,6];cosines=np.cos(np.pi/np.arange(2,65));norms=np.array(d['norms'],float);counts=Counter();candidates=[];near=[]
 profiles=[]
 for size in [1,2,3]:
  for support in combinations(range(7),size):
   choices=product([3,4,6],repeat=size)if size<3 else[(3,3,3)]
   for ll in choices:
    c=np.zeros(7);c[list(support)]=[2*math.cos(math.pi/l)for l in ll];profiles.append((support,ll,c))
 for support,ll,c in profiles:
  wc=c@J;w=wc@N[S];wq=wc@c
  if wq<2+1e-8:continue
  for mm in product(labels,repeat=2):
   rhs=np.array([w@Q@N[i]-math.sqrt(2*norms[i])*math.cos(math.pi/m)for i,m in zip([61,75],mm)])
   v=P@rhs;aa=direction@C@direction;bb=2*v@C@direction;cc=v@C@v+wq-2
   roots=np.roots([aa,bb,cc])if abs(aa)>1e-12 else[-cc/bb]if abs(bb)>1e-12 else[]
   for z in roots:
    if abs(np.imag(z))>1e-7:continue
    coeff=v+float(np.real(z))*direction;counts['real_branches']+=1
    if min(coeff)<-1e-8:continue
    r=-w+coeff@frame;cut=r@Q@E
    if cut<=1e-7:continue
    counts['escape_cutters']+=1;p=r@Q@N.T
    if np.max(p)>1e-7:continue
    counts['acute_cutters']+=1;cos=-p/np.sqrt(2*norms);finite=(cos>1e-7)&(cos<1-1e-7)
    err=np.min(abs(cos[finite,None]-cosines[None,:]),axis=1);bad=np.flatnonzero(finite)[err>1e-7]
    row=dict(E7_support=list(support),E7_labels=list(ll),anchor_walls=[61,75],anchor_labels=list(mm),source_coefficients=coeff.tolist(),normal=r.tolist(),escape_pairing=float(cut))
    if len(bad):
     if len(near)<30 or len(bad)<near[-1]['bad_count']:
      near.append(dict(row,bad_count=len(bad),bad_walls=bad.tolist(),bad_cosines=cos[bad].tolist()));near.sort(key=lambda a:a['bad_count']);near=near[:30]
    else:assert abs(r@Q@r-2)<1e-5;candidates.append(dict(row,pairings=p.tolist()))
 paths=[B.BASE,B.FACE,f'scripts/{B.TAG}_scout.py',f'scripts/{B.TAG}_two_anchor_scout.py']
 out=dict(status='BOUNDED_NUMERICAL_TWO_OFF_FACE_ANCHOR_SCOUT',sources={p:B.sha(p)for p in paths},profile_count=len(profiles),E7_support_sizes=[1,2,3],E7_labels='all3/4/6 for one/two support; all3 for triples',source_spatial_support=[0,12],anchor_walls=[61,75],anchor_labels=labels,counts=dict(counts),candidates=candidates,near_misses=near,elapsed_seconds=time.time()-start,scope='Two exact-definition linear anchor equations and norm quadratic, numerically evaluated. No continuous/all-label exclusion and no exact acceptance.')
 (B.ROOT/f'results/{B.TAG}_two_anchor_scout.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 d=run();print(d['status'],d['counts'],len(d['candidates']))
