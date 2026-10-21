#!/usr/bin/env python3
"""Bounded direct integral Gale proposals; every retained matrix is checked exactly.

Numerical SVD only constructs starting proposals. No numerical determinant or
incidence test is used. This is not an exhaustive Gale or oriented-matroid census.
"""
from pathlib import Path
import random,itertools,sys,collections,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
N=14;B=list(itertools.combinations(range(N),4));BI={b:i for i,b in enumerate(B)};C=list(itertools.combinations(range(N),5));IX=np.array(B);CI=np.array([[BI[c[:i]+c[i+1:]]for i in range(5)]for c in C]);CS=np.array([1,-1,1,-1,1]);MASK=np.array([sum(1<<i for i in c)for c in C]);r=random.Random(20261014);counts=collections.Counter();ex=[];t=time.time()
def calc(V):
 assert V.shape==(14,4) and np.max(np.abs(V))<=8192
 assert 24*8192**4 < 2**63
 a=V[IX];out=np.zeros(len(B),dtype=np.int64)
 for perm in itertools.permutations(range(4)):
  sign=(-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4));out+=sign*a[:,0,perm[0]]*a[:,1,perm[1]]*a[:,2,perm[2]]*a[:,3,perm[3]]
 return out
import json, importlib, scipy.linalg as la
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
seen=set();rows=[];seeds=[]
def accept(V):
 det=calc(V)
 if np.any(det==0):return None
 co=np.sign(det)[CI]*CS;pc=MASK[np.all(co==co[:,[0]],axis=1)]
 if not len(pc)or np.bitwise_or.reduce(pc)!=16383:return None
 vf=[16383^int(m)for m in pc]
 if any(not any(v>>i&1 for v in vf)for i in range(14)):return None
 mm=clutter(vf,14)
 if any(m.bit_count()>5 or m.bit_count()==3 for m in mm)or sum(m.bit_count()==2 for m in mm)<2:return None
 return vf,mm
for k in [0,1,3,4,5,7,9,31]:
 mod=importlib.import_module('coxeter_polytopes.combinatorics.h9_14_record'+str(k));coords=mod.BASE_COORDINATES;labels=sorted(coords)
 A=np.array([[1.]+list(map(float,coords[i]))for i in labels]).T;g=la.null_space(A)
 if hasattr(mod,'BASE_BLOCKS'):blocks=mod.BASE_BLOCKS
 else:blocks=getattr(mod,'CLONE_BLOCKS',None)or mod.CONSTRUCTION_CLONE_BLOCKS
 target=getattr(mod,'TARGET_RELABEL',None)or getattr(mod,'CONSTRUCTION_TO_TARGET',None)
 V=np.zeros((14,4))
 for row,label in enumerate(labels):
  for j in blocks.get(label,(label,)):V[target[j]if target else j]=g[row]
 V=V/np.linalg.norm(V,axis=1)[:,None]
 base=np.rint(2048*V).astype(np.int64)
 for z in range(1000):
  trial=base+np.array([[r.randint(-4,4)for j in range(4)]for i in range(14)])
  ac=accept(trial)
  if ac is not None:seeds.append((k,trial.tolist()));break
 else:print('BAD SEED',k);continue
 V=trial.copy();vf,mm=ac;cc,pp=canonical(mm,14);print('seed',k,len(vf),collections.Counter(m.bit_count()for m in mm),flush=True)
 for step in range(10000):
  tr=V.copy();i=r.randrange(14);j=r.randrange(4);tr[i,j]+=r.choice([-1,1])*r.choice([2,8,32,128,512]);
  if np.max(np.abs(tr))>8192:counts['coordinate_bound']+=1;continue
  ac=accept(tr)
  if ac is None:counts['rejected']+=1;continue
  V=tr;vf,mm=ac;counts['accepted']+=1
  can,perm=canonical(mm,14)
  if can in seen:continue
  seen.add(can);rows.append(dict(seed=k,step=step,V=V.tolist(),vertices=vf,missing=mm,canonical=can,permutation=perm))
  
  if len(rows)%25==0:print('NEW',len(rows),k,step,len(vf),flush=True)
json.dump(dict(seeds=seeds,counts=dict(counts),rows=rows),open(Path(__file__).resolve().parents[1]/'results/fresh_h9_gale_frontier_20261014_walk.json','w'))
print(counts,'distinct',len(rows),'elapsed',time.time()-t,flush=True)
