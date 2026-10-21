#!/usr/bin/env python3
"""Exact bounded hyperplane slices of products of simplices.

Proposal families have explicit rational convex realizations. Literal local
Coxeter obstructions are necessary tests, with missing pairs left dotted and
unassigned. This is not a census of all product sections or all H9 nerves.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
from hashlib import sha256
import numpy as np
from fresh_h9_rank5_gale_20261014_core import LocalGates,canonical
from fresh_h9_nonproduct_incidence_20261014_multiwedge import join
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_structured_incidence_20261014'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return [i for i in range(m.bit_length())if m>>i&1]

class Gates(LocalGates):
 def check(self,mm):
  for a,b in combinations(mm,2):
   sa,sb=a.bit_count(),b.bit_count()
   if min(sa,sb)<3:continue
   u=a|b;inside=[c for c in mm if c&u==c]
   if not(a&b)and set(inside)=={a,b}and sorted((sa,sb))in([3,4],[3,5],[4,4],[4,5],[5,5]):
    return dict(gate='two_lanner_blocks',blocks=[a,b],support=u)
   if u.bit_count()>9 or any(m.bit_count()==2 for m in inside):continue
   vv=bits(u);red=tuple(sorted(sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in inside));key=len(vv),red
   if key not in self.cache:self.cache[key]=canonical(red,len(vv))
   can,order=self.cache[key]
   if(len(vv),can)in self.bad:
    path,ident=self.bad[len(vv),can]
    return dict(gate='literal_induced_pattern',source=path,pattern=ident,support=u,permutation=[vv[j]for j in order])
  return None


def recipes():
 # Permuting equal-dimensional factors and reversing the entire affine
 # functional are safe symmetries. No claim of complete type deduplication.
 for w in combinations_with_replacement(range(1,5),10):
  yield 'ten_cube',[[0,x]for x in w]
 shapes=[(0,a,b)for b in range(2,5)for a in range(1,b)]
 for c in combinations_with_replacement(shapes,5):yield 'five_triangles',list(map(list,c))
 for dims in [(1,1,1,3,4),(1,1,2,2,4),(1,1,2,3,3),(1,2,2,2,3)]:
  groups=[(d,dims.count(d))for d in sorted(set(dims))]
  choices=[list(combinations_with_replacement(range(1,4),n))for d,n in groups]
  for ws in product(*choices):
   co=[list(range(0,(d+1)*w,w))for(d,n),v in zip(groups,ws)for w in v]
   yield 'mixed_simplex_products',co


def interval(co,mask):
 lo=hi=0;p=0
 for c in co:
  v=[x for j,x in enumerate(c)if not(mask>>(p+j)&1)]
  if not v:return None
  lo+=min(v);hi+=max(v);p+=len(c)
 return lo,hi


def triangle_data(co,h):
 n=sum(map(len,co));X=np.arange(1<<n,dtype=np.int64)
 low=np.zeros(len(X),dtype=np.int64);high=low.copy();valid=np.ones(len(X),dtype=bool);p=0
 for c in co:
  k=len(c);ii=(X>>p)&((1<<k)-1);ls=[];hs=[]
  for m in range(1<<k):
   v=[x for j,x in enumerate(c)if not(m>>j&1)]
   ls.append(min(v)if v else 0);hs.append(max(v)if v else 0)
  low+=np.array(ls)[ii];high+=np.array(hs)[ii];valid&=ii!=(1<<k)-1;p+=k
 feasible=valid&(2*low<h)&(h<2*high)
 minimal=~feasible
 for i in range(n):minimal&=((X&(1<<i))==0)|feasible[X^(1<<i)]
 active=[i for i in range(n)if feasible[1<<i]]
 return active,[int(m)for m in X[minimal]if m.bit_count()>1]


def cube_data(co,h):
 w=[c[1]for c in co];W=sum(w);n=len(w);sums=[0]*(1<<n)
 for m in range(1,1<<n):b=m&-m;sums[m]=sums[m^b]+w[b.bit_length()-1]
 active=[2*i+j for i in range(n)for j in range(2)if (2*w[i]<h if j==0 else h<2*(W-w[i]))]
 aset=set(active);mm=[]
 for i in range(n):
  if 2*i in aset and 2*i+1 in aset:mm.append(3<<(2*i))
 # Removing barycentric vertex 0 forces x=1; removing vertex 1 forces x=0.
 for side,t in [(0,h),(1,2*W-h)]:
  for m in range(1,1<<n):
   if 2*sums[m]>t and all(2*sums[m^(1<<i)]<t for i in bits(m)):
    a=sum(1<<(2*i+side)for i in bits(m))
    if a.bit_count()>1 and all(i in aset for i in bits(a)):mm.append(a)
 return active,sorted(mm)


def vertices(co,h,active):
 offsets=[];n=0
 for c in co:offsets.append(n);n+=len(c)
 amap={a:i for i,a in enumerate(active)};out={};target=F(h,2)
 for k,c in enumerate(co):
  others=[j for j in range(len(co))if j!=k]
  for picks in product(*[range(len(co[j]))for j in others]):
   base=sum(co[j][v]for j,v in zip(others,picks))
   for a,b in combinations(range(len(c)),2):
    if not base+c[a]<target<base+c[b]:continue
    s=(target-base-c[a])/F(c[b]-c[a]);point=[F(0)]*n
    for j,v in zip(others,picks):point[offsets[j]+v]=1
    point[offsets[k]+a]=1-s;point[offsets[k]+b]=s
    zeros=[i for i,x in enumerate(point)if x==0]
    assert all(i in amap for i in zeros)
    v=sum(1<<amap[i]for i in zeros);assert v.bit_count()==9 and v not in out
    assert all(x>=0 for x in point)
    assert sum(F(a)*x for c0,p0 in zip(co,offsets)for a,x in zip(c0,point[p0:p0+len(c0)]))==target
    out[v]=dict(free_factor=k,free_pair=[a,b],other_vertices=list(picks),fraction=str(s))
 assert out
 assert set(Counter(v^(1<<i)for v in out for i in bits(v)).values())=={2}
 assert all(any(v>>i&1 for v in out)for i in range(len(active)))
 return [dict(incident=v,**out[v])for v in sorted(out)]


def components(mm,n):
 left=set(range(n));out=[]
 while left:
  hit={min(left)}
  while True:
   nh=hit|{i for m in mm if any(m>>j&1 for j in hit)for i in bits(m)}
   if nh==hit:break
   hit=nh
  out.append(sorted(hit));left-=hit
 return out


def run():
 gates=Gates();stats=Counter();rows=[];cases=[];seen={};byfamily={};newkey=0
 for rid,(family,co)in enumerate(recipes()):
  assert sum(len(c)-1 for c in co)==10 and all(c==sorted(set(c))for c in co)
  counts=byfamily.setdefault(family,Counter());W=sum(c[-1]for c in co)
  for h in range(1,W+1,2):
   counts['parameter_recipes']+=1;stats['parameter_recipes']+=1
   active,oldmm=(cube_data if family=='ten_cube'else triangle_data)(co,h)
   base=dict(recipe=rid,family=family,coefficients=co,twice_level=h,facet_count=len(active))
   if len(active)<14:
    row=base|dict(verdict='OUTSIDE_AT_LEAST_14_FACET_SCOPE',active_original_facets=active)
   else:
    big=next((m for m in oldmm if m.bit_count()>5),None)
    if big is not None:row=base|dict(verdict='EXACT_LARGE_MINIMAL_NONFACE',witness_original=big)
    else:
     amap={a:i for i,a in enumerate(active)};mm=sorted(sum(1<<amap[i]for i in bits(m))for m in oldmm)
     can,order=canonical(mm,len(active));key=(len(active),tuple(can))
     if key in seen:
      cid=seen[key];target=cases[cid]['canonical_order'];mp={order[k]:target[k]for k in range(len(order))}
      row=base|dict(verdict='DUPLICATE_NERVE',case=cid,isomorphism_to_case=[mp[i]for i in range(len(order))])
     else:
      cid=len(cases);seen[key]=cid;v=vertices(co,h,active);cc=components(mm,len(active));gate=gates.check(mm)
      row=base|dict(case=cid,active_original_facets=active,minimal_nonfaces=mm,canonical_minimal_nonfaces=list(can),canonical_order=list(order),
                   vertex_records=v,product_factors=cc,missing_size_counts=dict(Counter(map(int.bit_count,mm))))
      if gate:row.update(verdict='EXACT_INHERITED_OBSTRUCTION',gate=gate)
      else:
       trial=join([tuple(bits(m))for m in mm],[tuple(bits(z['incident']))for z in v],5000)
       row.update(verdict=trial['verdict'],label_join=trial)
      cases.append(row);print('case',cid,family,'facets',len(active),'vertices',len(v),'factors',list(map(len,cc)),row['verdict'],flush=True)
   counts[row['verdict']]+=1;stats[row['verdict']]+=1
   rows.append({k:v for k,v in row.items()if k not in('vertex_records','minimal_nonfaces','canonical_minimal_nonfaces','canonical_order','product_factors','missing_size_counts','gate','label_join')})
  if rid%100==0:print('recipe',rid,'counts',dict(stats),flush=True)
 files=[f'scripts/{TAG}.py','scripts/fresh_h9_rank5_gale_20261014_core.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py','scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py',
  'scripts/independent_h9_14_record30_a1_labels.py']+gates.paths
 return dict(status='TERMINAL_BOUNDED_EXACT_PRODUCT_SECTION_SEARCH',sources={p:sha(p)for p in files},
  settings=dict(cube_weights=[1,2,3,4],triangle_shapes=[[0,a,b]for b in range(2,5)for a in range(1,b)],mixed_factor_weights=[1,2,3],
   half_integer_levels='0 < twice_level <= sum(factor_maxima), twice_level odd',minimum_facets=14,max_join_states=5000),
  counts=dict(stats),families={k:dict(v)for k,v in byfamily.items()},records=rows,cases=cases,
  scope='Finite explicit rational slices of ten-dimensional products of simplices. Generic half-integer levels avoid every source vertex and give simple nine-polytopes. Cases below14facets are outside this search scope, not excluded. Missing pairs are never assigned finite labels by the inherited dictionary; they retain unspecified dotted distances. A label-join survivor/limit remains unresolved. This is not all product sections or all compact H9 incidences.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 d=json.loads(json.dumps(d))
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],flush=True)
