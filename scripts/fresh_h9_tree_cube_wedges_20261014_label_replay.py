#!/usr/bin/env python3
"""Saved completed searches with a separate exact Gram-positivity predicate.

DFS traversal is reused; the finite-Dynkin classifier is replaced by cycle
and large-label witnesses, or Fraction Schur elimination on weighted trees.
Node-limited searches are deliberately not promoted to exclusions.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_h9_tree_cube_wedges_20261014_bridge_csp as traversal
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tree_cube_wedges_20261014'
Z=(F(0),F(0));O=(F(1),F(0));COS2={3:(F(1,4),F(0)),4:(F(1,2),F(0)),5:(F(3,8),F(1,8)),6:(F(3,4),F(0))}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return-a[0],-a[1]
def mul(a,b):return a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def inv(a):
 d=a[0]*a[0]-5*a[1]*a[1];assert d;return a[0]/d,-a[1]/d
def sign(a):
 x,y=a;s=lambda x:int(x>0)-int(x<0)
 if not y:return s(x)
 if not x:return s(y)
 if s(x)==s(y):return s(x)
 return s(x)*s(x*x-5*y*y)
def run():
 data=json.loads((ROOT/f'results/{TAG}.json').read_text());paths=[f'scripts/{TAG}_label_replay.py',f'scripts/{TAG}_bridge_csp.py',f'results/{TAG}.json']
 expected=[]
 for r in data['representatives']:
  p=f'results/{TAG}_case_{r["id"]:02d}.json';e=json.loads((ROOT/p).read_text());paths.append(p)
  for q,h in e['sources'].items():assert sha(q)==h,q
  if e['status']!='NODE_LIMIT_UNRESOLVED':expected.append(e)
 assert set(traversal.codes(3))=={c for c in product(range(2,8),repeat=3)if sum(F(1,m)for m in c)<1}
 calls=Counter();table={}
 @lru_cache(None)
 def pivot(tree):
  q=O
  for m,child in tree:
   p=pivot(child)
   if p is None:return None
   q=add(q,neg(mul(COS2[m],inv(p))))
  table[tree]=q
  return q if sign(q)>0 else None
 def independent(component,labels):
  calls['component_calls']+=1
  edges=[(i,j,labels[i,j])for i in sorted(component)for j in sorted(component)if i<j and labels[i,j]>2]
  n=len(component);assert n>=3
  if len(edges)>=n:calls['cycle_nonpositive']+=1;return False
  assert len(edges)==n-1
  if any(m>=7 for i,j,m in edges):
   # An incident three-node path: cos(pi/m)^2>3/4, its other
   # bond has square>=1/4, so the principal determinant is negative.
   calls['large_label_path_negative']+=1;return False
  adj={i:[]for i in component}
  for i,j,m in edges:adj[i].append((j,m));adj[j].append((i,m))
  deg={i:len(adj[i])for i in component};live=set(component)
  while len(live)>2:
   leaves=[i for i in live if deg[i]==1];assert leaves
   for i in leaves:
    live.remove(i)
    for j,m in adj[i]:deg[j]-=1
  def rooted(i,parent):return tuple(sorted((m,rooted(j,i))for j,m in adj[i]if j!=parent))
  tree=min(rooted(i,None)for i in live);ok=pivot(tree)is not None
  calls['tree_positive'if ok else'tree_nonpositive']+=1;return ok
 traversal.finite_component=independent;records=[]
 for e in expected:
  got=json.loads(json.dumps(traversal.run(e['case'],e['settings']['node_limit'])))
  assert got==e
  records.append(dict(case=e['case'],status=e['status'],nodes=e['nodes'],prunes=e['prunes'],bridge_trace_sha256=e['bridge_trace_sha256']))
  print('independent labels',e['case'],e['status'],e['nodes'],flush=True)
 return dict(status='PASS_COMPLETED_DFS_WITH_SEPARATE_TREE_SCHUR_ARITHMETIC',sources={p:sha(p)for p in paths},
  completed_cases=records,component_counts=dict(calls),distinct_schur_trees=len(table),
  schur_table=[dict(rooted_tree=k,pivot=list(map(str,v)),positive=sign(v)>0)for k,v in sorted(table.items())],
  scope='Reuses the producer DFS and exact geometry, but replaces every Dynkin finite-component predicate by independent Fraction arithmetic in Q(sqrt5), cycle witnesses, or an all-m>=7 negative 3x3 path minor. Reproduces full branch hashes/counts for every saved exhausted case and every first partial-label survivor; node-limited cases remain unresolved and are not replayed. Surviving mandatory labels are not metric realizations.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_label_replay.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['completed_cases']),flush=True)
