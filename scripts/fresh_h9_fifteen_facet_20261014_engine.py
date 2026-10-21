#!/usr/bin/env python3
"""Necessary fifteen-facet Coxeter label solver, including missing triangles.

Port of the separately reviewed14-facet propagation with the ambient facet
count changed to15. Finite vertex rank remains9 and allowed Gram rank10.
No fourteen-facet-specific restriction is used.
"""
import json
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import pynauty
from fresh_h9_highbase_construction_20261014_cell import finite,codes
from fresh_h9_bugaenko_closure_20261014 import residue
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical as nerve_canonical
ROOT=Path(__file__).resolve().parents[1]
N=15;NP=N*(N-1)//2;PAIRS=list(combinations(range(N),2));PI={p:i for i,p in enumerate(PAIRS)}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def tup(m):return tuple(i for i in range(N)if m>>i&1)
class Gates:
 def __init__(self,row):
  self.M=list(map(tup,row['canonical_minimal_nonfaces']));self.V=list(map(tup,row['canonical_vertices']))
  assert all(len(v)==9 for v in self.V)and all(2<=len(m)<=5 for m in self.M)
  self.vp=[tuple(PI[p]for p in combinations(v,2))for v in self.V]
  self.dotted={PI[m]for m in self.M if len(m)==2}
  self.cliques=[(a,tuple(PI[p]for p in combinations(a,2)))for a in combinations(range(N),11)if all(PI[p]not in self.dotted for p in combinations(a,2))]
  self.orthogonal=[(i,j,tuple(PI[tuple(sorted((a,b)))]for a in A for b in B))for i,A in enumerate(self.M)for j,B in enumerate(self.M)if i<j and not set(A)&set(B)]
  self.det_cache={}
 def bad_vertex(self,labels):
  return next((i for i,pp in enumerate(self.vp)if not finite(9,tuple(labels[p]for p in pp))),None)
 def proof(self,labels,known):
  bad=self.bad_vertex(labels)
  if bad is not None:return dict(kind='NONSPHERICAL_VERTEX',vertex=bad)
  for i,j,pp in self.orthogonal:
   if all(p in known and labels[p]==2 for p in pp):return dict(kind='ORTHOGONAL_LANNER',missing_faces=[i,j])
  for a,pp in self.cliques:
   if not all(p in known and labels[p]<=6 for p in pp):continue
   c=tuple(labels[p]for p in pp)
   if c not in self.det_cache:self.det_cache[c]=residue(c)
   if self.det_cache[c]:return dict(kind='RANK_AT_LEAST_ELEVEN',support=a,residue=self.det_cache[c])
  return None
from functools import lru_cache
@lru_cache(None)
def six_relation():
 p='results/fresh_h9_newfacet_frontier_20261014_relations.json'
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return tuple(tuple(r)for r in d['six_node_relation']['complete_labels'])
def search(row,limit):
 g=Gates(row);constraints=[(tuple(PI[p]for p in combinations(a,2)),codes(len(a)))for a in g.M if len(a)>=3]
 # Add every occurrence of the complete96-row six-node relation.
 mm=row['canonical_minimal_nonfaces']; seen=set(); extra=0
 for a,b in combinations(mm,2):
  if a.bit_count()!=5 or b.bit_count()!=5 or(a&b).bit_count()!=4:continue
  u=a|b
  if any(m not in(a,b)and m&u==m for m in mm):continue
  if u in seen:continue
  seen.add(u);vv=[i for i in range(N)if u>>i&1]
  red=[sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in(a,b)]
  canon,order=nerve_canonical(red,6);assert canon==(61,62)
  perm=[vv[j]for j in order]
  pp=tuple(PI[tuple(sorted((perm[i],perm[j])))]for i,j in combinations(range(6),2))
  constraints.append((pp,six_relation()));extra+=1
 affected=[tuple(i for i,v in enumerate(g.vp)if p in v)for p in range(NP)]
 nodes=visits=0;stats=Counter();proofs=[];candidate=None;dead=set()
 twins=defaultdict(set)
 for i in range(N):twins[tuple(j for j,m in enumerate(g.M)if i in m)].add(i)
 twin_blocks=list(twins.values());adj={i:[]for i in range(N+NP)}
 for k,(i,j)in enumerate(PAIRS,N):adj[i].append(k);adj[j].append(k);adj[k]=[i,j]
 def canonical(left,labels,known):
  values=[-1 if p in g.dotted else(labels[p]if p in known else 0)for p in range(NP)];ec=defaultdict(set)
  for p,c in enumerate(values,N):ec[c].add(p)
  graph=pynauty.Graph(N+NP,adjacency_dict=adj,vertex_coloring=twin_blocks+[ec[c]for c in sorted(ec)])
  order=pynauty.canon_label(graph)[:N];assert set(order)==set(range(N))
  pos=0
  for b in twin_blocks:assert set(order[pos:pos+len(b)])==b;pos+=len(b)
  # The key is the explicitly relabelled state, not only a graph digest.
  return left,tuple(values[PI[tuple(sorted((order[i],order[j])))]]for i,j in PAIRS)
 def force_right(labels,known):
  out=set(known)
  for p in range(NP):
   if p in known or p in g.dotted:continue
   labels[p]=3
   if any(not finite(9,tuple(labels[k]for k in g.vp[i]))for i in affected[p]):out.add(p)
   labels[p]=2
  stats['forced_right_occurrences']+=len(out)-len(known)
  return out
 def walk(left,labels,known):
  nonlocal nodes,visits,candidate
  visits+=1;key=canonical(left,labels,known)
  if key in dead:stats['cache_hits']+=1;return False
  nodes+=1
  if nodes>limit:return None
  bad=g.bad_vertex(labels)
  if bad is not None:stats['NONSPHERICAL_VERTEX']+=1;dead.add(key);return False
  known=force_right(labels,known)
  pr=g.proof(labels,known)
  if pr:
   stats[pr['kind']]+=1
   if len(proofs)<100:proofs.append(dict(labels=labels[:],known=sorted(known),proof=pr))
   dead.add(key);return False
  pending=[i for i in range(len(constraints))if left>>i&1]
  if pending:
   # Every table row must extend the current spherical relaxation.
   # Rejection is monotone: all as-yet-unassigned cosines are zero.
   options={}
   for i in pending:
    pp,pool=constraints[i];keep=[]
    for code in pool:
     if not all(p not in known or code[j]==labels[p]for j,p in enumerate(pp)):continue
     trial=labels[:]
     for p,m in zip(pp,code):trial[p]=m
     stats['table_rows_tested']+=1
     if g.bad_vertex(trial)is not None:stats['table_rows_nonspherical']+=1;continue
     keep.append(code)
    if not keep:stats['NO_LANNER_EXTENSION']+=1;dead.add(key);return False
    options[i]=keep
   # Pairwise relation consistency on overlapping missing faces. These are
   # necessary projections of the full join, never guesses of a shared label.
   changed=True
   while changed:
    changed=False
    for i in pending:
     pp=constraints[i][0]
     for j in pending:
      if i==j:continue
      qq=constraints[j][0];common=sorted(set(pp)&set(qq))
      if not common:continue
      ii=[pp.index(p)for p in common];jj=[qq.index(p)for p in common]
      allowed={tuple(c[k]for k in jj)for c in options[j]}
      keep=[c for c in options[i]if tuple(c[k]for k in ii)in allowed]
      if len(keep)!=len(options[i]):
       stats['arc_rows_removed']+=len(options[i])-len(keep);changed=True;options[i]=keep
      if not keep:stats['EMPTY_RELATION_JOIN']+=1;dead.add(key);return False
   # Values common to every row are forced in every possible completion.
   forced={}
   for i in pending:
    pp=constraints[i][0]
    for j,p in enumerate(pp):
     values={c[j]for c in options[i]}
     if p not in known and len(values)==1:forced[p]=next(iter(values))
   if forced:
    child=labels[:]
    for p,m in forced.items():child[p]=m
    stats['table_forced_pairs']+=len(forced)
    return walk(left,child,known|set(forced))
   k=min(pending,key=lambda i:(len(options[i]),-len(constraints[i][0]),-sum(p in known for p in constraints[i][0]),i))
   pp=constraints[k][0]
   for code in options[k]:
    child=labels[:]
    for p,m in zip(pp,code):child[p]=m
    r=walk(left&~(1<<k),child,known|set(pp))
    if r is not False:return r
  else:
   unknown=set(range(NP))-known-g.dotted
   if not unknown:candidate=labels[:];return True
   # Complete the nearest fixed minor rather than deferring rank checks.
   targets=[(sum(p not in known for p in pp),pp)for a,pp in g.cliques if all(labels[p]<=6 for p in pp if p in known)]
   target=min(targets,key=lambda t:t[0])[1]if targets else []
   p=min(unknown,key=lambda p:(p not in target,-len(affected[p]),p))
   for m in range(2,8):
    child=labels[:];child[p]=m;r=walk(left,child,known|{p})
    if r is not False:return r
  dead.add(key);return False
 result=walk((1<<len(constraints))-1,[2]*NP,set())
 out=dict(case=row['case'],six_node_relations=extra,verdict='FULL_CATEGORY_METRIC_REQUIRED'if result is True else('NODE_LIMIT_UNRESOLVED'if result is None else'EXACT_COMPLETE_OBSTRUCTION'),nodes=nodes,visits=visits,twin_blocks=[sorted(b)for b in twin_blocks],counters=dict(stats),sample_proofs=proofs,distinct_modular_minor_codes=len(g.det_cache))
 if candidate is not None:out.update(label_vector=candidate,dotted_pairs=[PAIRS[p]for p in sorted(g.dotted)])
 return out
