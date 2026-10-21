#!/usr/bin/env python3
"""Separate stdlib convex-section/nerve audit; no producer or label imports."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations,combinations_with_replacement
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_structured_incidence_20261014';SOURCE=f'results/{TAG}.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):
 while m:
  b=m&-m;yield b.bit_length()-1;m-=b

def interval(co,m):
 p=0;lo=hi=0
 for c in co:
  allowed=[a for k,a in enumerate(c)if not(m>>(p+k)&1)]
  if not allowed:return None
  lo+=min(allowed);hi+=max(allowed);p+=len(c)
 return lo,hi

def feasible(co,h,m):
 z=interval(co,m)
 return z is not None and 2*z[0]<h<2*z[1]

def recipe_list():
 for w in combinations_with_replacement(range(1,5),10):yield 'ten_cube',[[0,a]for a in w]
 ss=sorted((0,a,b)for b in range(2,5)for a in range(1,b))
 # Producer uses b-then-a ordering; equality below compares sets of tuples.
 for v in combinations_with_replacement(ss,5):yield 'five_triangles',[list(c)for c in v]
 for dims in [(1,1,1,3,4),(1,1,2,2,4),(1,1,2,3,3),(1,2,2,2,3)]:
  groups=[(i,dims.count(i))for i in sorted(set(dims))]
  for a in product(*[list(combinations_with_replacement(range(1,4),k))for i,k in groups]):
   yield 'mixed_simplex_products',[[j*w for j in range(i+1)]for(i,k),ws in zip(groups,a)for w in ws]

def recipe_key(f,c,h):return f,tuple(sorted(map(tuple,c))),h

def edge_vertices(co,h,active):
 off=[];n=0
 for c in co:off.append(n);n+=len(c)
 index={v:i for i,v in enumerate(active)};source=list(product(*[range(len(c))for c in co]));target=F(h,2)
 oldheight={v:sum(co[k][j]for k,j in enumerate(v))for v in source};out={}
 assert all(F(v)!=target for v in oldheight.values())
 for v in source:
  for k,c in enumerate(co):
   for j in range(v[k]+1,len(c)):
    w=v[:k]+(j,)+v[k+1:];a,b=oldheight[v],oldheight[w]
    if not a<target<b:continue
    lam=(target-a)/(b-a);point=[F(0)]*n
    for factor,(vi,wi)in enumerate(zip(v,w)):
     point[off[factor]+vi]+=1-lam;point[off[factor]+wi]+=lam
    assert all(x>=0 for x in point)
    assert all(sum(point[o:o+len(c)])==1 for o,c in zip(off,co))
    assert sum(a*x for c,o in zip(co,off)for a,x in zip(c,point[o:o+len(c)]))==target
    zero=[j for j,x in enumerate(point)if not x];assert all(j in index for j in zero)
    mask=sum(1<<index[j]for j in zero);assert mask.bit_count()==9 and mask not in out
    out[mask]=point
 assert all(any(v>>i&1 for v in out)for i in range(len(active)))
 return out

def missing_from_facets(vertices,n):
 # Full downward face complex, then check every possible minimal nonface.
 # No threshold-missing-face formula is used here.
 faces=set()
 for v in vertices:
  z=v
  while z:faces.add(z);z=(z-1)&v
 faces.add(0)
 boundary=set()
 for f in faces:
  for i in range(n):
   if not(f>>i&1):
    m=f|(1<<i)
    if m not in faces:boundary.add(m)
 return sorted(m for m in boundary if all(m^(1<<i)in faces for i in bits(m)))

def gate(r,d):
 g=r['gate'];mm=r['minimal_nonfaces'];u=g['support'];inside=[m for m in mm if m&u==m]
 if g['gate']=='two_lanner_blocks':
  a,b=g['blocks'];assert not(a&b)and a|b==u and set(inside)=={a,b}
  assert sorted([a.bit_count(),b.bit_count()])in [[3,4],[3,5],[4,4],[4,5],[5,5]]
 else:
  p=g['permutation'];assert set(p)==set(bits(u));rel=sorted(sum(1<<j for j,i in enumerate(p)if m>>i&1)for m in inside)
  src=json.loads((ROOT/g['source']).read_text())
  if g['pattern']=='seven_node_relation':assert rel==[31,99]
  else:
   rr=next(x for x in src['patterns']if x['pattern']==g['pattern']);assert rr['verdict']=='EXACT_LOCAL_OBSTRUCTION'
   assert rel==rr['minimal_nonfaces']

def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 expected={recipe_key(f,c,h)for f,c in recipe_list()for h in range(1,sum(x[-1]for x in c)+1,2)}
 actual=[recipe_key(r['family'],r['coefficients'],r['twice_level'])for r in d['records']]
 assert len(actual)==len(set(actual)) and set(actual)==expected
 counts=Counter();geom=[];checks=0;small={}
 for r in d['records']:
  c=r['coefficients'];h=r['twice_level'];n=sum(map(len,c));assert sum(len(x)-1 for x in c)==10
  act=[j for j in range(n)if feasible(c,h,1<<j)];assert len(act)==r['facet_count']
  v=r['verdict'];counts[v]+=1
  if v=='OUTSIDE_AT_LEAST_14_FACET_SCOPE':assert len(act)<14 and act==r['active_original_facets'];continue
  if v=='EXACT_LARGE_MINIMAL_NONFACE':
   m=r['witness_original'];assert m.bit_count()>5 and all(i in act for i in bits(m))
   assert not feasible(c,h,m)and all(feasible(c,h,m^(1<<i))for i in bits(m));checks+=1+m.bit_count();continue
  # Eligible and duplicate recipes receive a fresh full convex vertex/nerve check.
  verts=edge_vertices(c,h,act);mm=missing_from_facets(verts,len(act));assert max(map(int.bit_count,mm))<=5
  target=d['cases'][r['case']]
  if v=='DUPLICATE_NERVE':
   p=r['isomorphism_to_case'];assert sorted(p)==list(range(len(act)))
   assert sorted(sum(1<<p[i]for i in bits(m))for m in mm)==target['minimal_nonfaces']
  else:
   assert act==target['active_original_facets']and mm==target['minimal_nonfaces']
   assert sorted(verts)==[z['incident']for z in target['vertex_records']]
   assert set(Counter(v^(1<<i)for v in verts for i in bits(v)).values())=={2}
   small[r['case']]=verts
   if v=='EXACT_INHERITED_OBSTRUCTION':gate(target,d)
   else:
    trial=target['label_join'];assert v=='EXACT_MANDATORY_LANNER_VERTEX_OBSTRUCTION'
    assert trial['stages'][-1]['passing_states']==0 and not trial['state_limit_reached']
    for stage in trial['stages']:
     assert not stage['truncated'] and sum(1<<i for i in stage['missing_face'])in mm
  geom.append(dict(recipe=r['recipe'],twice_level=h,case=r['case'],facets=len(act),vertices=len(verts)))
 assert dict(counts)=={k:v for k,v in d['counts'].items()if k!='parameter_recipes'}
 paths=[SOURCE,f'scripts/{TAG}_independent.py']
 return dict(status='PASS_SEPARATE_EXACT_STRUCTURED_SECTION_GEOMETRY',sources={p:sha(p)for p in paths},
  parameter_recipes=len(actual),large_missing_face_feasibility_checks=checks,eligible_recipes=len(geom),
  exact_edge_intersection_vertices=sum(r['vertices']for r in geom),distinct_saved_cases=len(small),records=geom,
  scope='No producer imports. Independent rational source-edge intersections, complete downward face complexes, literal duplicate isomorphisms, large-missing-face minimality and local obstruction applicability. The inherited obstruction arithmetic and complete Lanner joins are not reimplemented here; primary saved replay covers the latter.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();p=ROOT/f'results/{TAG}_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['parameter_recipes'],d['eligible_recipes'],flush=True)
