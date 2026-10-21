#!/usr/bin/env python3
"""Replay the frozen DFS with an independent exact positivity predicate.

Traversal is reused and stated explicitly. The inherited finite-diagram
classification is replaced by cycle witnesses and Schur elimination of
weighted trees over Q(sqrt5), cached only by exact rooted-tree shape.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import fresh_h9_coupled_cube_wedges_20261014_bridge_csp as traversal
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_coupled_cube_wedges_20261014'
Z=(F(0),F(0));O=(F(1),F(0));COS2={3:(F(1,4),F(0)),4:(F(1,2),F(0)),5:(F(3,8),F(1,8)),6:(F(3,4),F(0))}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return-x[0],-x[1]
def mul(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def inv(x):
 q=x[0]*x[0]-5*x[1]*x[1];assert q;return x[0]/q,-x[1]/q
def sign(x):
 a,b=x;s=lambda a:int(a>0)-int(a<0)
 if not b:return s(a)
 if not a:return s(b)
 if s(a)==s(b):return s(a)
 return s(a)*s(a*a-5*b*b)
def run():
 expected=json.loads((ROOT/f'results/{TAG}_bridge_csp.json').read_text())
 for p,h in expected['sources'].items():assert sha(p)==h,p
 assert expected['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
 counts=Counter();table={}
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
  counts['component_calls']+=1
  edges=[(a,b,labels[a,b])for a in sorted(component)for b in sorted(component)if a<b and labels[a,b]>2]
  n=len(component);assert n>=3
  if len(edges)>=n:
   counts['cycle_nonpositive_witness']+=1;return False
  assert len(edges)==n-1
  if any(m>=7 for a,b,m in edges):
   # A three-node path contains that edge and another edge>=3.
   # cos(pi/m)>cos(pi/6)=sqrt3/2 for every m>=7.
   # The adjacent edge has cosine>=1/2, so the path determinant
   # is strictly less than 1-3/4-1/4=0.
   assert F(1)-F(3,4)-F(1,4)==0
   counts['large_label_path_witness']+=1;return False
  adj={i:[]for i in component}
  for a,b,m in edges:adj[a].append((b,m));adj[b].append((a,m))
  deg={i:len(adj[i])for i in component};live=set(component)
  while len(live)>2:
   leaves=[i for i in live if deg[i]==1]
   assert leaves
   for i in leaves:
    live.remove(i)
    for j,m in adj[i]:deg[j]-=1
  def rooted(i,parent):return tuple(sorted((m,rooted(j,i))for j,m in adj[i]if j!=parent))
  tree=min(rooted(i,None)for i in live);out=pivot(tree)is not None
  counts['tree_positive'if out else'tree_nonpositive']+=1
  return out
 traversal.finite_component=independent
 got=traversal.run(expected['settings']['node_limit']);got=json.loads(json.dumps(got))
 assert got==expected
 entries=[dict(rooted_tree=k,schur_pivot=[str(a),str(b)],positive=sign((a,b))>0)for k,(a,b)in sorted(table.items())]
 paths=[f'scripts/{TAG}_bridge_csp.py',f'results/{TAG}_bridge_csp.json',f'results/{TAG}_independent.json',f'scripts/{TAG}_bridge_independent.py']
 return dict(status='PASS_FULL_DFS_WITH_SEPARATE_TREE_SCHUR_ARITHMETIC',sources={p:sha(p)for p in paths},
  reproduced_nodes=got['nodes'],reproduced_prunes=got['prunes'],counts=dict(counts),distinct_schur_subtrees=len(entries),schur_table=entries,
  exact_large_label_bound='For m>=7, cos(pi/m)>cos(pi/6)=sqrt3/2; an incident label>=3 path therefore has determinant<1-3/4-1/4=0',
  scope='The entire saved DFS traversal and branch digest are reproduced, replacing every inherited finite-component classification with independent Fraction tree Schur elimination over Q(sqrt5), or explicit cycle/large-label path witnesses. The traversal is reused, not independently reimplemented. Geometry has a separate exhaustive primal replay. Labels>=7 are quantified categories and dotted distances remain continuous.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_bridge_independent.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['reproduced_nodes'],d['distinct_schur_subtrees'],flush=True)
