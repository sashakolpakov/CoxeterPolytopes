#!/usr/bin/env python3
"""Exact incidence/gate helpers for direct fifteen-point rank-five Gale data."""
from collections import Counter
from itertools import combinations,permutations
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_rank5_gale_20261014'
N=15;R=5;FULL=(1<<N)-1
BASES=list(combinations(range(N),R));BI={b:i for i,b in enumerate(BASES)};IX=np.array(BASES)
CIRCUITS=list(combinations(range(N),R+1));CM=np.array([sum(1<<i for i in c)for c in CIRCUITS],dtype=np.int64)
CI=np.array([[BI[c[:i]+c[i+1:]]for i in range(R+1)]for c in CIRCUITS])
ALT=np.array([1,-1,1,-1,1,-1]);PERMS=[(p,(-1)**sum(p[i]>p[j]for i in range(R)for j in range(i+1,R)))for p in permutations(range(R))]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def checked(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def exact_minors(V):
 V=np.asarray(V,dtype=np.int64);assert V.shape==(N,R)and np.max(np.abs(V))<=2048
 assert 120*2048**5<2**63
 A=V[IX];out=np.zeros(len(BASES),dtype=np.int64)
 for p,s in PERMS:
  term=np.ones(len(BASES),dtype=np.int64)
  for i in range(R):term*=A[:,i,p[i]]
  out+=s*term
 return out
def incidence(V,require_small=True):
 det=exact_minors(V)
 if np.any(det==0):return None,'nonuniform'
 signs=np.sign(det)[CI]*ALT
 pc=CM[np.all(signs==signs[:,[0]],axis=1)]
 if not len(pc)or int(np.bitwise_or.reduce(pc))!=FULL:return None,'no_positive_kernel_cover'
 facets=sorted(FULL^int(c)for c in pc)
 if not all(any(f>>i&1 for f in facets)for i in range(N)):return None,'unused_vertex'
 mm=clutter(facets,N)
 if require_small and max(map(int.bit_count,mm))>5:return None,'large_minimal_nonface'
 return dict(facets=facets,minimal_nonfaces=list(mm),det=det),None
def simplex_link_factors(facets):
 """Literal sufficient test: an isolated one-edge component in link clutter."""
 found=[]
 for v in range(N):
  ff=[f^(1<<v)for f in facets if f>>v&1];mask=0
  for f in ff:mask|=f
  vertices=[i for i in range(N)if mask>>i&1]
  reduced=[sum(1<<k for k,i in enumerate(vertices)if f>>i&1)for f in ff]
  mm=clutter(reduced,len(vertices))
  for m in mm:
   if m.bit_count()>=2 and all(m==k or not(m&k)for k in mm):
    found.append(dict(vertex=v,factor=sum(1<<vertices[i]for i in range(len(vertices))if m>>i&1)))
 return found
def realize(V,data):
 V=np.asarray(V,dtype=np.int64);det=data['det'];weights=[0]*N
 for cc,indices in zip(CIRCUITS,CI):
  cof=[int(det[j])*int(s)for j,s in zip(indices,ALT)]
  if all(x>0 for x in cof)or all(x<0 for x in cof):
   for i,x in zip(cc,cof):weights[i]+=abs(x)
 assert all(w>0 for w in weights)
 assert all(sum(weights[i]*int(V[i,j])for i in range(N))==0 for j in range(R))
 facets=data['facets'];mm=data['minimal_nonfaces'];can,order=canonical(mm,N);inverse={v:i for i,v in enumerate(order)}
 cv=sorted(sum(1<<inverse[i]for i in range(N)if f>>i&1)for f in facets)
 ridges=Counter(f^(1<<i)for f in facets for i in range(N)if f>>i&1)
 assert set(ridges.values())=={2}
 factors=simplex_link_factors(facets);assert not factors
 return dict(gale=V.tolist(),positive_kernel=weights,nonzero_five_minors=len(det),minimum_absolute_minor=int(np.min(np.abs(det))),
  facets=facets,minimal_nonfaces=mm,canonical_minimal_nonfaces=list(can),canonical_permutation=list(order),canonical_vertices=cv,
  missing_size_counts=dict(Counter(m.bit_count()for m in mm)),ridge_count=len(ridges),simplex_boundary_link_factors=factors)
class LocalGates:
 def __init__(self):
  self.cache={};self.bad={}
  self.paths=['results/fresh_h9_newfacet_frontier_20261014_patterns.json','results/fresh_h9_newfacet_frontier_20261014_local_independent.json',
   'results/fresh_h9_newfacet_frontier_20261014_relations.json','results/fresh_h9_bugaenko_remaining_20261014_block_gate.json',
   'results/fresh_h9_bugaenko_remaining_20261014_block_gate_review.json','results/fresh_h9_fifteen_facet_20261014_patterns.json']
  for path in(self.paths[0],self.paths[-1]):
   for row in checked(path)['patterns']:
    if row['verdict']=='EXACT_LOCAL_OBSTRUCTION':self.bad[row['order'],tuple(row['minimal_nonfaces'])]=(path,row['pattern'])
  self.bad[7,(31,99)]=(self.paths[2],'seven_node_relation')
  for path in self.paths[1:-1]:checked(path)
 def check(self,mm):
  for a,b in combinations(mm,2):
   sa,sb=a.bit_count(),b.bit_count()
   if min(sa,sb)<3:continue
   u=a|b;inside=[c for c in mm if c&u==c]
   if not(a&b)and set(inside)=={a,b}and sorted((sa,sb))in([3,4],[3,5],[4,4],[4,5],[5,5]):
    return dict(gate='two_lanner_blocks',blocks=[a,b],support=u)
   if u.bit_count()>9 or any(m.bit_count()==2 for m in inside):continue
   vv=[i for i in range(N)if u>>i&1]
   red=tuple(sorted(sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in inside));key=len(vv),red
   if key not in self.cache:self.cache[key]=canonical(red,len(vv))
   can,order=self.cache[key]
   if(len(vv),can)in self.bad:
    path,ident=self.bad[len(vv),can]
    return dict(gate='literal_induced_pattern',source=path,pattern=ident,support=u,permutation=[vv[j]for j in order])
  return None
