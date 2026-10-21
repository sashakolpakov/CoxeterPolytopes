#!/usr/bin/env python3
"""Continue the exact bridge search after the T-case8 proved continuous metric families.

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
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_other_trees_20261014_limited_search';SOURCE='results/fresh_h9_tree_cube_wedges_20261014.json'
TRANSFERS='results/fresh_h9_other_trees_20261014_limited_transfers.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check_manifest(path,seen=None):
 if seen is None:seen=set()
 if path in seen:return
 seen.add(path);obj=json.loads((ROOT/path).read_text())
 for p,h in obj.get('sources',{}).items():
  assert sha(p)==h,p
  if p.endswith('.json'):check_manifest(p,seen)
 return obj
def run_case(g,limit,case):
 transfers=check_manifest(TRANSFERS)
 assert transfers['status']=='PASS_LITERAL_REVIEWED_THEOREM_TRANSFERS_TO_LIMITED_CASE5'
 transfer=next(r for r in transfers['cases']if r['case']==case)
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
 nodes=0;prune=Counter();solution=None;stopped=False;arr=[0]*153;trace=sha256();samples=[];frontier=None; sparse_samples=[];sparse_trace=sha256()
 families=[dict(id=f"{row['kind']}_{i}",fixed=row['fixed'],bounds=row['bounds'],witness_map=row['witness_map'])for i,row in enumerate(transfer['families'])]
 for family in families:
  family['conditions']=[(index[a,b],'eq',m)for a,b,m in family['fixed']]+[(index[a,b],'ge',m)for a,b,m in family['bounds']]
 keys={k for f in families for k,op,m in f['conditions']}
 incidence={k:[(j,T[j].index(k))for j in range(len(T))if k in T[j]]for k in keys}
 assert all(incidence.values())
 def sparse_gate(dom):
  for family in families:
   witnesses=[]
   for k,op,m in family['conditions']:
    if arr[k]:possible={arr[k]}
    else:
     possible=None
     for j,pos in incidence[k]:
      values={row[pos]for row in dom[j]}
      possible=values if possible is None else possible&values
     if not possible:break
    if not all(x==m if op=='eq'else x>=m for x in possible):break
    witnesses.append(dict(pair=pairs[k],relation=op,bound=m,possible_categories=sorted(possible),already_assigned=bool(arr[k])))
   else:return family['id'],witnesses
  return None
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
 def dfs(dom,remaining,path):
  nonlocal nodes,solution,stopped,frontier
  if nodes>=limit:stopped=True;frontier=path;return False
  nodes+=1
  if nodes%10000==0:print('nodes',nodes,'depth',len(path),'prunes',dict(prune),'cache',finite.cache_info().currsize,flush=True)
  sparse_hit=sparse_gate(dom)
  if sparse_hit is not None:
   family,witnesses=sparse_hit
   prune['certified_'+family+'_partial']+=1
   record=dict(path=path,family=family,conditions=witnesses,remaining_triangle_count=len(remaining))
   sparse_trace.update(json.dumps(record,sort_keys=True,separators=(',',':')).encode())
   if len(sparse_samples)<30:sparse_samples.append(record)
   return False
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
 paths=[TRANSFERS,SOURCE,f'scripts/{TAG}.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='FULL_MANDATORY_BRIDGE_SURVIVOR'if solution else'NODE_LIMIT_UNRESOLVED'if stopped else'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY',
  sources={p:sha(p)for p in paths},proved_family_conditions=families,settings=dict(node_limit=limit,stop_at_first_survivor=True),nodes=nodes,prunes=dict(prune),
  triangle_supports=ts,pairs=pairs,bridge_clauses=clauses,bridge_trace_sha256=trace.hexdigest(),prune_samples=samples,
  solution=solution,unexplored_frontier=frontier,sparse_partial_samples=sparse_samples,sparse_trace_sha256=sparse_trace.hexdigest(),
  scope='Bounded necessary-label search on one other literal tree incidence, using only locally verified transferred analytic theorems. All finite labels>=7 remain one exact category; dotted distances remain arbitrary. Full mandatory survivors still require omitted meeting contacts and a Lorentz metric. Empty completion excludes only this incidence, while a node limit remains unresolved.')

def run(case,limit):
 data=json.loads((ROOT/SOURCE).read_text())
 for p,h in data['sources'].items():assert sha(p)==h,p
 assert case==5
 row=data['representatives'][case];assert row['id']==case
 result=run_case(row['geometry'],limit,case)
 result['case']=case;result['tree']=row['tree'];result['recipe_multiplicity']=row['recipe_multiplicity']
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',type=int,required=True);p.add_argument('--limit',type=int,default=100000);p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run(a.case,a.limit)));out=ROOT/f'results/{TAG}_case_{a.case:02d}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print('CASE',a.case,d['status'],d['nodes'],d['prunes'],flush=True)
