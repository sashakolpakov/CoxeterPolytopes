#!/usr/bin/env python3
"""Exact rational finite-support census on the two new polar walls.

Actuality uses the all-wall acute/common-center certificate, once bound.
This script imports no floating point geometry or producer arithmetic.
"""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_polar_pair_20261014_census'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PRIMARY='results/fresh_h20_polar_pair_20261014.json'
IDS=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,24,25,26,27,28,31,32,37,38,47,48,49,50,55,56]
POLARS=[[7,14,21,18,15,12,9,6,3,1,-1,-3,-5,-6,-7,-8,-9,-10,10,0,0],
 [14,27,40,35,30,25,20,15,10,5,0,-5,-10,-12,-14,-16,-18,-20,20,0,0]]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(y)for y in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 old=json.loads((ROOT/BASE).read_text());primary=json.loads((ROOT/PRIMARY).read_text())
 for doc in [old,primary]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert primary['status']=='PASS_100_WALL_RATIONAL_POLAR_PAIR_STILL_INFINITE'
 assert primary['recipes']==POLARS
 Q=[list(map(F,row))for row in old['form']];N={}
 for i in IDS:
  v=[]
  for co in old['normals'][i]:
   assert not any(map(F,co[1])) and not any(map(F,co[0][1:]))
   v.append(F(co[0][0]))
  N[i]=v
 N.update({98:list(map(F,POLARS[0])),99:list(map(F,POLARS[1]))})
 sparse=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def dual(v):return [sum(q*v[j]for j,q in row)for row in sparse]
 def dot(a,b):return sum(x*y for x,y in zip(a,dual(b)))
 G={(i,j):dot(v,w)for i,v in N.items()for j,w in N.items()}
 assert G[98,98]==2 and G[99,99]==4 and G[98,99]==0
 assert all(G[i,j]<=0 for i,j in combinations(N,2))
 finite={i:{j for j in N if j!=i and G[i,j]**2<G[i,i]*G[j,j]}for i in N}
 neighbor={r:sorted(finite[r])for r in[98,99]}
 assert len(neighbor[98])==36 and len(neighbor[99])==25
 assert [neighbor[98],neighbor[99]]==primary['complete_finite_neighbors']
 # The complete full-wall pair certificate closes the rational-pool premise.
 polar_zero_ranks=[]
 for r,ids in zip([98,99],primary['old98_zero_sets']):
  assert set(ids)<=set(IDS) and all(G[r,i]==0 for i in ids)
  M=[N[i][:]for i in ids];k=0
  for j in range(21):
   p=next((i for i in range(k,len(M))if M[i][j]),None)
   if p is None:continue
   M[k],M[p]=M[p],M[k];v=M[k][j];M[k]=[x/v for x in M[k]]
   for i in range(k+1,len(M)):
    if M[i][j]:
     v=M[i][j];M[i]=[a-v*b for a,b in zip(M[i],M[k])]
   k+=1
  assert k==20;polar_zero_ranks.append(k)
 def maxcliques(vertices):
  result=[]
  def visit(R,P,X):
   if not P and not X:result.append(tuple(sorted(R)));return
   u=max(P|X,key=lambda v:(len(P&finite[v]),-v))
   for v in sorted(P-finite[u]):
    visit(R|{v},P&finite[v],X&finite[v]);P.remove(v);X.add(v)
  visit(set(),set(vertices),set());return sorted(result)
 cache={};counts=Counter()
 def pd_component(ids):
  ids=tuple(sorted(ids))
  if ids in cache:return cache[ids]['positive']
  M=[[G[i,j]for j in ids]for i in ids];piv=[];ok=True
  for k in range(len(ids)):
   p=M[k][k];piv.append(p)
   if p<=0:ok=False;break
   for i in range(k+1,len(ids)):
    if not M[i][k]:continue
    for j in range(k+1,len(ids)):
     if M[k][j]:M[i][j]-=M[i][k]*M[k][j]/p
  cache[ids]=dict(roots=ids,positive=ok,leading_pivots=piv)
  return ok
 def pd(ids):
  unseen=set(ids)
  while unseen:
   first=min(unseen);todo=[first];comp={first};unseen.remove(first)
   while todo:
    v=todo.pop();extra={j for j in unseen if G[v,j]!=0};comp|=extra;unseen-=extra;todo.extend(extra)
   if not pd_component(comp):return False
  return True
 cliques={r:maxcliques(neighbor[r])for r in[98,99]}
 supports=set();tested=set()
 for r in[98,99]:
  for cc in cliques[r]:
   for chosen in combinations(cc,19):
    ids=tuple(sorted(chosen+(r,)))
    if ids in tested:continue
    tested.add(ids)
    if pd(ids):supports.add(ids)
 print('PD supports',len(supports),'tested',len(tested),'components',len(cache),flush=True)
 H=[sum(v[k]for v in N.values())for k in range(21)];assert dot(H,H)<0
 def nullpoint(ids):
  M=[dual(N[i])for i in ids];piv=[];k=0
  for j in range(21):
   p=next((i for i in range(k,20)if M[i][j]),None)
   if p is None:continue
   M[k],M[p]=M[p],M[k];v=M[k][j];M[k]=[x/v for x in M[k]]
   for i in range(20):
    if i!=k and M[i][j]:
     v=M[i][j];M[i]=[a-v*b for a,b in zip(M[i],M[k])]
   piv.append(j);k+=1
   if k==20:break
  assert len(piv)==20
  free=next(j for j in range(21)if j not in piv);x=[F(0)]*21;x[free]=1
  for i,j in enumerate(piv):x[j]=-M[i][free]
  den=math.lcm(*(v.denominator for v in x));x=[int(v*den)for v in x];g=reduce(math.gcd,x);x=[F(v//g)for v in x]
  if dot(x,H)>0:x=[-v for v in x]
  assert dot(x,x)<0 and dot(x,H)<0 and all(dot(x,N[i])==0 for i in ids)
  assert all(dot(x,N[i])<0 for i in N if i not in ids)
  return x
 vertices=[]
 for ids in sorted(supports):
  x=nullpoint(ids);vertices.append(dict(incident=ids,primitive_rational_point=x,norm=dot(x,x),future_pairing=dot(x,H)))
 assert len({tuple(v['primitive_rational_point'])for v in vertices})==len(vertices)
 for v in vertices:
  for r in[98,99]:
   if r in v['incident']:counts['on_'+str(r)]+=1
  if 98 in v['incident'] and 99 in v['incident']:counts['on_both']+=1
 paths=[BASE,PRIMARY,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_RATIONAL_POLAR_FINITE_SUPPORT_CENSUS',sources={p:sha(p)for p in paths},
  new_normals=POLARS,new_norms=[2,4],old_zero_set_ranks=polar_zero_ranks,rational_old_wall_pool=IDS,complete_finite_neighbors=neighbor,
  maximal_finite_neighbor_cliques=cliques,distinct_twenty_subsets=len(tested),exact_component_checks=list(cache.values()),
  positive_support_count=len(vertices),vertex_counts=dict(counts),known_global_finite_vertices_at_least=primary['known_finite_vertex_count_at_least']+len(vertices),future_reference_point=H,future_reference_norm=dot(H,H),vertices=vertices,
  scope='Exact complete finite-vertex census on new walls98/99. The bound100-wall certificate supplies complete finite-angle neighbor lists, acuteness and a common strict center; projection to each positive20-support proves actuality and strictness at every other wall. The16 distinct primitive rational points are checked directly against every rational-pool wall. All old finite vertices are strictly preserved. No global old-vertex census or finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['positive_support_count'],d['vertex_counts'],flush=True)
