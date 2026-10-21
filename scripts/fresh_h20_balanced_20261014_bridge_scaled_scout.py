#!/usr/bin/env python3
"""Bounded scaled-product-face central bridges; numerical proposals only."""
import json,math,time
from itertools import product
import numpy as np
from fresh_h20_balanced_20261014 import ROOT,data,sha
from fresh_h20_balanced_20261014_bridge_scout import reps,allowed
TAG='fresh_h20_balanced_20261014_bridge_scaled_scout'
def run():
 start=time.time();D=data();V=np.array([v['point']for v in D['vertices']]);fv=D['factor_vertices'];factors=[reps(h)for h in[7,7,8]]
 cs=np.array([math.cos(math.pi/m)for m in range(2,61)]+[1.]);labels=list(range(2,61))+[0]
 profiles=[]
 for n,m in product([3,4,5,6],[2,3,4,5,6,0]):
  a=math.cos(math.pi/n);b=1 if m==0 else math.cos(math.pi/m);N=196
  aa=N-1/a**2;bb=2/a**2+4*b/a;cc=-4-1/a**2-4*b/a
  g=(-bb+math.sqrt(bb*bb-4*aa*cc))/(2*aa);r=(1-g)/a
  profiles.append((n,m,g,r))
 radialC=np.concatenate([g*V for n,m,g,r in profiles]);radialR=np.concatenate([np.full(392,r)for n,m,g,r in profiles])
 accepted=[];contacts=[];count=0
 for ia,ib,ic in product(*[range(len(f))for f in factors]):
  if ia>ib:continue
  ss=[factors[k][j][0]for k,j in enumerate([ia,ib,ic])]
  center0=sum(np.mean(fv[k][list(s)],axis=0)for k,s in enumerate(ss));N=float(center0@center0)
  if N<=4+1e-9:continue
  count+=1
  for n,m,_,_ in profiles:
   a=math.cos(math.pi/n);b=1 if m==0 else math.cos(math.pi/m)
   aa=N-1/a**2;bb=2/a**2+4*b/a;cc=-4-1/a**2-4*b/a
   if aa<=0:continue
   g=(-bb+math.sqrt(bb*bb-4*aa*cc))/(2*aa);r=(1-g)/a
   if not 0<g<1 or r<1e-6:continue
   center=g*center0;gaps=D['offsets']-D['vectors']@center
   if not allowed(gaps/r,cs).all():continue
   co=(np.sum((radialC-center)**2,axis=1)-radialR**2-r*r)/(2*radialR*r)
   pos=np.searchsorted(cs[:-1],co);j=np.clip(pos,0,len(cs)-2);jm=np.clip(pos-1,0,len(cs)-2)
   finite=np.flatnonzero(np.minimum(abs(co-cs[j]),abs(co-cs[jm]))<1e-8)
   row=dict(subsets=[list(s)for s in ss],unscaled_norm=N,seed_label=n,central_label=m,gamma=g,radius=r,radial_finite_contacts=finite.tolist())
   accepted.append(row)
   if len(finite):contacts.append(row)
 out=dict(status='BOUNDED_NUMERICAL_SCALED_FACE_BRIDGES',sources={p:sha(p)for p in['scripts/fresh_h20_balanced_20261014.py','scripts/fresh_h20_balanced_20261014_bridge_scout.py',f'scripts/{TAG}.py']},product_representatives_with_norm_over4=count,labels_seed=[3,4,5,6],labels_central=[2,3,4,5,6,'parallel'],finite_pair_labels_tested=list(range(2,61)),accepted_representatives=len(accepted),with_radial_finite_contacts=len(contacts),recipes=accepted,elapsed_seconds=time.time()-start,scope='Numerical necessary seed/central gates; finite contacts to the24 radial families tested. No exhaustive labels, pair classification proof, clique, or volume certificate.')
 (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');print('accepted',len(accepted),'with radial finite contacts',len(contacts));print(contacts[:5])
if __name__=='__main__':run()
