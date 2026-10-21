#!/usr/bin/env python3
"""Depth-first exact necessary triangle/vertex CSP for the new cube wedges.

Only mandatory-triangle contacts are assigned. A terminal state still needs
all omitted meeting contacts and the nine continuous dotted distances.
"""
import argparse,json,time
from collections import defaultdict,Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from hashlib import sha256
from fresh_h9_highbase_construction_20261014_cell import codes
from independent_h9_14_record30_a1_labels import finite_component
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_coupled_cube_wedges_20261014';SOURCE=f'results/{TAG}.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run(limit):
 data=json.loads((ROOT/SOURCE).read_text());g=data['geometry']
 for p,h in data['sources'].items():assert sha(p)==h
 ts=[tuple(m)for m in g['minimal_nonfaces']if len(m)==3];vs=[tuple(z['incident'])for z in g['vertices']]
 pairs=sorted({p for t in ts for p in combinations(t,2)});index={p:i for i,p in enumerate(pairs)}
 T=[tuple(index[p]for p in combinations(t,2))for t in ts]
 masks={sum(1<<k for k,(i,j)in enumerate(pairs)if i in v and j in v)for v in vs}
 masks=sorted(m for m in masks if not any(m!=n and m&n==m for n in masks))
 V=[tuple(k for k in range(len(pairs))if m>>k&1)for m in masks]
 opts=list(codes(3));opts.sort(key=lambda c:(max(c),sum(c),c.count(2),c))
 domains=[opts[:]for t in ts];swaps=[]
 for a,b in [(1,2),(4,5),(7,8),(10,11)]:
  perm=list(range(18));perm[a],perm[b]=b,a
  assert {tuple(sorted(perm[i]for i in v))for v in vs}==set(vs)
  k=next(k for k,t in enumerate(ts)if a in t and b in t)
  leaf=next(i for i in ts[k]if i not in(a,b));pa=tuple(sorted((leaf,a)));pb=tuple(sorted((leaf,b)))
  ia=list(combinations(ts[k],2)).index(pa);ib=list(combinations(ts[k],2)).index(pb)
  domains[k]=[c for c in domains[k]if c[ia]<=c[ib]];swaps.append(dict(twin_pair=[a,b],triangle=k,ordered_contact_positions=[ia,ib]))
 @lru_cache(1000000)
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
 nodes=0;leaves=0;prune=Counter();solution=None;stopped=False;stats=[];a=[0]*len(pairs)
 def dfs(dom,remaining,path):
  nonlocal nodes,leaves,solution,stopped
  if nodes>=limit:stopped=True;return False
  nodes+=1
  if nodes%10000==0:print('nodes',nodes,'depth',len(path),'finitecache',finite.cache_info().currsize,flush=True)
  if not remaining:
   leaves+=1;solution=dict(pair_labels=[list(p)+[m]for p,m in zip(pairs,a)],triangle_path=path);return True
  k=min(remaining,key=lambda k:(len(dom[k]),-sum(a[i]!=0 for i in T[k]),-sum(len(set(T[k])&set(T[j]))for j in remaining if j!=k),k))
  nxt=remaining-{k}
  for c in dom[k]:
   assert all(not a[i]or a[i]==m for i,m in zip(T[k],c))
   changed=[i for i in T[k]if not a[i]]
   for i,m in zip(T[k],c):a[i]=m
   affected=[v for v in V if any(i in v and a[i]>2 for i in changed)]
   good=all(finite(tuple((i,a[i])for i in v if a[i]>2))for v in affected)
   if not good:prune['spherical_vertex']+=1
   else:
    nd=dom[:]
    for j in nxt:
     nd[j]=[o for o in dom[j]if all(not a[i]or a[i]==m for i,m in zip(T[j],o))]
     if not nd[j]:good=False;prune['empty_triangle_domain']+=1;break
    if good and dfs(nd,nxt,path+[dict(triangle=k,labels=c)]):return True
   for i in changed:a[i]=0
   if stopped:return False
  return False
 dfs(domains,set(range(len(ts))),[])
 paths=[SOURCE,f'scripts/{TAG}_csp.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='FULL_MANDATORY_LABEL_SURVIVOR'if solution else'NODE_LIMIT_UNRESOLVED'if stopped else'EXACT_ALL_MANDATORY_LABELS_EMPTY',
  sources={p:sha(p)for p in paths},settings=dict(node_limit=limit,stop_at_first_survivor=True),nodes=nodes,terminal_assignments=leaves,prunes=dict(prune),
  triangles=ts,mandatory_pairs=pairs,actual_vertices=len(vs),maximal_pair_restrictions=len(V),symmetry_orders=swaps,solution=solution,
  scope='Exact DFS over every shared Lanner3 category, with category7 meaning every integer label>=7. Unassigned meeting pairs are relaxed to2 only in necessary spherical tests. Nine missing pairs retain independent strictly dotted distances. The four twin swaps are explicitly verified incidence automorphisms. Only a complete exhausted empty search is an exclusion; a full mandatory survivor requires additional finite contacts and a Lorentzian metric.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=200000);p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run(a.limit)));out=ROOT/f'results/{TAG}_csp.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['nodes'],d['prunes'],flush=True)
