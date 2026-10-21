#!/usr/bin/env python3
"""Second bounded Gale walk, biased toward few missing four-faces.

The exponential acceptance rule only guides proposals; saved integral matrices
receive exact incidence and Coxeter checks separately.
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
import json, math
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
D=json.load(open(Path(__file__).resolve().parents[1]/'results/fresh_h9_gale_frontier_20261014_walk.json'));seen={tuple(row['canonical'])for row in D['rows']};rows=[];seeds=[];r=random.Random(2026101402)
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
for seed in [0,1,3,4,5,7,9,31]:
 source=min((a for a in D['rows']if a['seed']==seed),key=lambda a:(sum(m.bit_count()==4 for m in a['missing']),sum(m.bit_count()==2 for m in a['missing']),-len(a['vertices'])))
 V=np.array(source['V'],dtype=np.int64);vf,mm=accept(V);score=sum(m.bit_count()==4 for m in mm);seeds.append(dict(seed=seed,V=V.tolist(),score=score))
 for step in range(12500):
  tr=V.copy();i=r.randrange(14);j=r.randrange(4);tr[i,j]+=r.choice([-1,1])*r.choice([4,16,64,256,1024]);
  if np.max(np.abs(tr))>8192:counts['coordinate_bound']+=1;continue
  ac=accept(tr)
  if ac is None:counts['incidence_rejected']+=1;continue
  vf,mm=ac;sc=sum(m.bit_count()==4 for m in mm)
  if sc>score and r.random()>math.exp((score-sc)/.6):counts['objective_rejected']+=1;continue
  V=tr;score=sc;counts['accepted']+=1
  can,perm=canonical(mm,14)
  if can in seen:continue
  seen.add(can);rows.append(dict(seed=seed,step=step,V=V.tolist(),vertices=vf,missing=mm,canonical=can,permutation=perm));counts['new_score_'+str(score)]+=1
  if len(rows)%20==0 or score==0:print('NEW',len(rows),seed,step,len(vf),collections.Counter(m.bit_count()for m in mm),flush=True)
json.dump(dict(seeds=seeds,counts=dict(counts),rows=rows),open(Path(__file__).resolve().parents[1]/'results/fresh_h9_gale_frontier_20261014_low4.json','w'));print(counts,'distinct',len(rows),'elapsed',time.time()-t,flush=True)
