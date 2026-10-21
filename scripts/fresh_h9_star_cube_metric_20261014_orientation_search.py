#!/usr/bin/env python3
"""Continue the exact bridge search after the certified sixteen-orientation cut-family failure.

Two disjoint induced missing faces both have negative Gram directions.
If all their cross contacts are forced right by spherical-vertex tests,
their union already has two negative directions. Missing pairs are treated
as arbitrary strict dotted pairs; no distance is substituted.
"""
import argparse,json
from collections import defaultdict,Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from hashlib import sha256
from fresh_h9_highbase_construction_20261014_cell import codes
from independent_h9_14_record30_a1_labels import finite_component
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014_orientation_search';SOURCE='results/fresh_h9_tree_cube_wedges_20261014.json'
METRIC='results/fresh_h9_star_cube_metric_20261014_orientations.json'
FIRST='results/fresh_h9_tree_cube_wedges_20261014_case_00.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run_case(g,limit):
 metric=json.loads((ROOT/METRIC).read_text())
 for p,h in metric['sources'].items():assert sha(p)==h,p
 assert metric['status']=='PASS_ALL_SIXTEEN_ORIENTATIONS_EXCLUDED'
 ms=[tuple(m)for m in g['minimal_nonfaces']];ts=[m for m in ms if len(m)==3];vs=[tuple(z['incident'])for z in g['vertices']]
 pairs=list(combinations(range(18),2));index={p:i for i,p in enumerate(pairs)}
 T=[tuple(index[p]for p in combinations(t,2))for t in ts]
 V=[tuple(index[p]for p in combinations(v,2))for v in vs]
 pv=[tuple(k for k,v in enumerate(V)if i in v)for i in range(len(pairs))]
 clauses=[]
 for aa,bb in combinations(ms,2):
  A,B=set(aa),set(bb)
  if A&B or sum(set(c)<=A|B for c in ms)!=2:continue
  cross=tuple(sorted(index[tuple(sorted((i,j)))]for i in aa for j in bb))
  assert all(pv[i]for i in cross)
  clauses.append(dict(blocks=[aa,bb],cross_pairs=cross))
 assert clauses
 opts=sorted(codes(3),key=lambda c:(max(c),sum(c),c.count(2),c));domains=[list(opts)for t in ts]
 for a,b in [b for b in g['blocks']if len(b)==2]:
  pm=list(range(18));pm[a],pm[b]=b,a
  assert{tuple(sorted(pm[i]for i in v))for v in vs}==set(vs)
  k=next(k for k,t in enumerate(ts)if a in t and b in t);leaf=next(i for i in ts[k]if i not in(a,b))
  cp=list(combinations(ts[k],2));u=cp.index(tuple(sorted((leaf,a))));v=cp.index(tuple(sorted((leaf,b))))
  domains[k]=[c for c in domains[k]if c[u]<=c[v]]
 @lru_cache(2000000)
 def finite(edges):
  labels=defaultdict(lambda:2);adj=defaultdict(set)
  for k,m in edges:
   a,b=pairs[k];labels[a,b]=m;adj[a].add(b);adj[b].add(a)
  left=set(adj)
  while left:
   start=min(left);co={start};q=[start];left.remove(start)
   for a in q:
    for b in adj[a]&left:left.remove(b);co.add(b);q.append(b)
   if len(co)>2 and not finite_component(co,labels):return False
  return True
 nodes=0;prune=Counter();solution=None;stopped=False;arr=[0]*153;trace=sha256();samples=[];frontier=None
 def bridge_gate():
  tested={};vk=[tuple((i,arr[i])for i in v if arr[i]>2)for v in V]
  def can(k):
   if arr[k]:return(arr[k]>2,-1)
   if k not in tested:
    tested[k]=(True,-1)
    for j in pv[k]:
     ed=tuple(sorted((*vk[j],(k,3))))
     if not finite(ed):tested[k]=(False,j);break
   return tested[k]
  for j,c in enumerate(clauses):
   fail=[]
   for k in c['cross_pairs']:
    ok,vi=can(k)
    if ok:break
    fail.append((k,vi))
   else:return j,fail
  return None
 def family_gate(dom,remaining):
  S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
  if any(arr[index[S[2*i],S[2*i+1]]]!=4 for i in range(4)):return False
  for i in range(4):
   if not any(all(arr[index[tuple(sorted((leaf,S[2*i+j])))]]==3 for j in range(2))for leaf in[P[i],C[i]]):return False
  return True
 def dfs(dom,remaining,path):
  nonlocal nodes,solution,stopped,frontier
  if nodes>=limit:stopped=True;frontier=path;return False
  nodes+=1
  if nodes%10000==0:print('nodes',nodes,'depth',len(path),'prunes',dict(prune),'cache',finite.cache_info().currsize,flush=True)
  if family_gate(dom,remaining):
   prune['certified_sixteen_orientations']+=1;return False
  bad=bridge_gate()
  if bad is not None:
   prune['orthogonal_negative_blocks']+=1
   record=dict(path=path,clause=bad[0],forced_cross_witnesses=bad[1]);trace.update(json.dumps(record,sort_keys=True,separators=(',',':')).encode())
   if len(samples)<20:samples.append(record)
   return False
  if not remaining:
   candidate=[list(p)+[m]for p,m in zip(pairs,arr)if m]
   solution=dict(pair_labels=candidate,triangle_path=path);return True
  k=min(remaining,key=lambda k:(len(dom[k]),-sum(arr[i]!=0 for i in T[k]),-sum(len(set(T[k])&set(T[j]))for j in remaining if j!=k),k));nxt=remaining-{k}
  for c in dom[k]:
   assert all(not arr[i]or arr[i]==m for i,m in zip(T[k],c))
   changed=[i for i in T[k]if not arr[i]]
   for i,m in zip(T[k],c):arr[i]=m
   affected=set(j for i in changed if arr[i]>2 for j in pv[i])
   good=all(finite(tuple((i,arr[i])for i in V[j]if arr[i]>2))for j in affected)
   if not good:prune['spherical_vertex']+=1
   else:
    nd=dom[:]
    for j in nxt:
     nd[j]=[o for o in dom[j]if all(not arr[i]or arr[i]==m for i,m in zip(T[j],o))]
     if not nd[j]:good=False;prune['empty_triangle_domain']+=1;break
    if good and dfs(nd,nxt,path+[dict(triangle=k,labels=c)]):return True
   for i in changed:arr[i]=0
   if stopped:return False
  return False
 dfs(domains,set(range(len(ts))),[])
 paths=[SOURCE,f'scripts/{TAG}.py',METRIC,FIRST,'scripts/fresh_h9_star_cube_metric_20261014_orientations.py','scripts/fresh_h9_tree_cube_wedges_20261014_bridge_csp.py','scripts/fresh_h9_coupled_cube_wedges_20261014_bridge_csp.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='FULL_MANDATORY_BRIDGE_SURVIVOR'if solution else'NODE_LIMIT_UNRESOLVED'if stopped else'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY',
  sources={p:sha(p)for p in paths},settings=dict(node_limit=limit,stop_at_first_survivor=True),nodes=nodes,prunes=dict(prune),
  triangle_supports=ts,pairs=pairs,bridge_clauses=clauses,bridge_trace_sha256=trace.hexdigest(),prune_samples=samples,
  solution=solution,unexplored_frontier=frontier,
  scope='Continuation of star case0 after excluding the source-bound sixteen orientation families by their exact rank and actual-vertex certificates. A partial state is excluded only when every B2 internal label is4 and each block already has one leaf with both own contacts3. No condition is imposed on the opposite leaves. All mandatory Lanner triangles and every actual spherical vertex are necessary. Category7 means all finite integer labels>=7, only in these classifications. If adding the smallest nonright contact3 to a relaxed vertex fails, monotonicity forces that cross contact right in every completion. Two disjoint missing-face blocks with all cross contacts forced right require two negative directions, forbidden in Lorentz signature. Missing pairs are never assigned a distance. A limit remains unresolved; a surviving mandatory assignment still requires all omitted meeting contacts and a global metric.')
def run(case,limit):
 data=json.loads((ROOT/SOURCE).read_text())
 for p,h in data['sources'].items():assert sha(p)==h,p
 assert case==0
 row=data['representatives'][case];assert row['id']==case
 result=run_case(row['geometry'],limit)
 result['case']=case;result['tree']=row['tree'];result['recipe_multiplicity']=row['recipe_multiplicity']
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',type=int,required=True);p.add_argument('--limit',type=int,default=100000);p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run(a.case,a.limit)));out=ROOT/f'results/{TAG}_case_{a.case:02d}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print('CASE',a.case,d['status'],d['nodes'],d['prunes'],flush=True)
