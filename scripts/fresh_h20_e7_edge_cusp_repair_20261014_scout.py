#!/usr/bin/env python3
"""Bounded nonincident contact events for five exact quotient cusp roots."""
import json,math
from pathlib import Path
from hashlib import sha256
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_edge_cusp_repair_20261014_scout';mp.mp.dps=70
BASE='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json';FRONT='results/fresh_h20_e7_finite_edgewalk_20261014_frontier.json'
def dec(x):
 a=str(x).split('/');return mp.mpf(a[0])/(mp.mpf(a[1])if len(a)>1 else 1)
def run():
 d=json.loads((ROOT/BASE).read_text());f=json.loads((ROOT/FRONT).read_text());rad=d['field']['base_square_root_basis']
 def base(x):return sum(dec(c)*mp.sqrt(q)for c,q in zip(x,rad))
 disc=base(d['field']['discriminant'])
 def ev(v):return np.array([float(base(a)+base(b)*mp.sqrt(disc))for a,b in v])
 def kv(v):return np.array([float(base(a))for a in v])
 N=np.array(list(map(ev,d['normals'])));H=ev(d['strict_center']);E=ev(d['strict_escape']['point']);Q=np.array([[float(dec(q))for q in row]for row in d['form']]);norms=np.array(list(map(float,d['norms'])))
 delta=kv(f['new_parallel_cusp']['point']);line=3*math.sqrt(3)/7*kv(f['new_parallel_cusp']['lineality_vector']);inc=set(f['new_parallel_cusp']['incident']);slopes=N@Q@delta
 bases=[-N[[12,13,18,31,65]].sum(axis=0),-N[15],-math.sqrt(2)*(N[17]+N[79])-N[83],line,-line]
 names=['A5','A1','C3_short','line_plus','line_minus'];vertices=np.array([kv(v['point'])for v in d['preserved_old_vertices']+d['new_cap_finite_vertices']])
 results=[]
 for name,a in zip(names,bases):
  assert abs(a@Q@a-2)<1e-8 and abs(a@Q@delta)<1e-8
  offset=N@Q@a;candidates=[]
  for j in range(90):
   if j in inc:continue
   assert slopes[j]<-1e-8
   for label in [2,3,4,6,'parallel']:
    co=1 if label=='parallel'else math.cos(math.pi/label)
    k=(-math.sqrt(2*norms[j])*co-offset[j])/slopes[j]
    if k<=0:continue
    r=a+k*delta;g=offset+k*slopes;c=-g/np.sqrt(2*norms)
    if r@Q@H>=-1e-8 or min(c)<-1e-7:continue
    if any(min(abs(z),abs(z-.5),abs(z-math.sqrt(.5)),abs(z-math.sqrt(.75)))>2e-7 and z<1-2e-7 for z in c):continue
    vp=vertices@Q@r
    if max(vp)>=-1e-8:continue
    finite=[i for i,z in enumerate(c)if z<1-1e-7]
    candidates.append(dict(anchor=j,label=label,k=float(k),nonincident_finite=[i for i in finite if i not in inc],escape_pair=float(r@Q@E),max_vertex_pair=float(max(vp))))
  candidates.sort(key=lambda x:(x['k'],-len(x['nonincident_finite']),x['anchor']))
  k=1
  while True:
   r=a+k*delta;c=-(offset+k*slopes)/np.sqrt(2*norms)
   if r@Q@H<0 and max(vertices@Q@r)<0 and all(c[j]>1+1e-6 for j in range(90)if j not in inc):break
   k*=2
  results.append(dict(name=name,events=candidates,safe_integer_shift=k,safe_escape_pair=float(r@Q@E)));print(name,len(candidates),candidates[:3],'fallback',k,flush=True)
 paths=[BASE,FRONT,f'scripts/{TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_EVENT_PROPOSALS_NOT_EXACT_ACCEPTANCE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},families=results,scope='Five fixed quotient roots; nonincident anchor labels2,3,4,6 orparallel, retaining84 specified vertices. Numerical event screen only, with safe largeinteger proposals. No all-angle/all-parameter exclusion.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n')
