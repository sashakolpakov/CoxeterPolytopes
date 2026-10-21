#!/usr/bin/env python3
"""Separate local CSP traversal and finite-tree classification replay.

No producer local search, maximal-face, canonicalization or propagation code
is imported. The classical order4/5 Lanner catalog is a shared dependency.
"""
import argparse,json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from fresh_h9_bugaenko_surgery_20261014_rank_independent import spherical
from independent_h9_14_record30_a1_labels import lanner_codes
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_local_independent'
PAT='results/fresh_h9_newfacet_frontier_20261014_patterns.json';REL='results/fresh_h9_newfacet_frontier_20261014_relations.json';INC='results/fresh_h9_bugaenko_remaining_20261014_newfacet.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return tuple(i for i in range(m.bit_length())if m>>i&1)
@lru_cache(None)
def table(n):
 if n==3:return tuple(c for c in product(range(2,8),repeat=3)if sum(Fraction(1,m)for m in c)<1)
 return tuple(tuple(c)for c in lanner_codes(n))
@lru_cache(None)
def finite(k,code):
 pp=dict(zip(combinations(range(k),2),code))
 return spherical(tuple(pp.get((i,j),2)for i,j in combinations(range(9),2)))
def faces(n,missing):
 out=[]
 for mask in range(1<<n):
  if any(mask&m==m for m in missing):continue
  if any(all((mask|(1<<i))&m!=m for m in missing)for i in range(n)if not mask>>i&1):continue
  out.append(bits(mask))
 return out

def solve(n,missing,all_solutions=False):
 pairs=list(combinations(range(n),2));index={p:i for i,p in enumerate(pairs)}
 vv=[(len(f),tuple(index[p]for p in combinations(f,2)))for f in faces(n,missing)]
 tt=[(tuple(index[p]for p in combinations(bits(m),2)),table(m.bit_count()))for m in missing]
 solutions=[];nodes=0
 def walk(todo,assignment):
  nonlocal nodes
  nodes+=1
  labels=[assignment.get(i,2)for i in range(len(pairs))]
  if any(not finite(nf,tuple(labels[p]for p in pp))for nf,pp in vv):return
  if todo:
   domains=[]
   for j in todo:
    pp,pool=tt[j];available=[]
    for row in pool:
     if all(p not in assignment or assignment[p]==m for p,m in zip(pp,row)):available.append(row)
    if not available:return
    domains.append((len(available),j,available))
   _,j,available=min(domains,key=lambda x:(x[0],x[1]));pp=tt[j][0]
   for row in available:
    child=assignment.copy();child.update(zip(pp,row));walk(tuple(i for i in todo if i!=j),child)
    if solutions and not all_solutions:return
  else:
   free=[p for p in range(len(pairs))if p not in assignment]
   if not free:solutions.append(labels);return
   p=free[0]
   for m in range(2,8):
    child=assignment.copy();child[p]=m;walk((),child)
    if solutions and not all_solutions:return
 walk(tuple(range(len(tt))),{});return solutions,nodes

def run():
 data={p:json.loads((ROOT/p).read_text())for p in(PAT,REL,INC)}
 for d in data.values():
  for p,h in d['sources'].items():assert sha(p)==h,p
 primary,relation,inc=(data[p]for p in(PAT,REL,INC));checked=[];excluded=set();literal_faces=0
 todo=[r for r in primary['patterns']if r['verdict']=='EXACT_LOCAL_OBSTRUCTION']+[relation['seven_node_obstruction']]
 for r in todo:
  n=r['order'];missing=r['minimal_nonfaces'];answer,nodes=solve(n,missing)
  assert not answer
  ff=faces(n,missing);assert sorted(sum(1<<j for j in f)for f in ff)==r['maximal_faces']
  for case,perm in r['applications']:
   assert len(perm)==len(set(perm))==n
   row=inc['cases'][case];union=sum(1<<i for i in perm)
   induced={m for m in row['canonical_minimal_nonfaces']if m&union==m}
   assert induced=={sum(1<<perm[i]for i in bits(m))for m in missing}
   for f in ff:
    lifted=sum(1<<perm[i]for i in f);assert any(v&lifted==lifted for v in row['canonical_vertices']);literal_faces+=1
   excluded.add(case)
  checked.append(dict(order=n,minimal_nonfaces=missing,separate_states=nodes,applications=len(r['applications'])))
 six=relation['six_node_relation'];solutions,nodes=solve(6,six['minimal_nonfaces'],True)
 assert sorted(solutions)==six['complete_labels']and len(solutions)==96
 assert excluded==set(relation['combined_excluded_cases'])and len(excluded)==3112
 paths=[PAT,REL,INC,f'scripts/{TAG}.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='PASS_SEPARATE_44_LOCAL_OBSTRUCTIONS_AND96_ROW_RELATION',sources={p:sha(p)for p in paths},obstruction_count=len(checked),records=checked,literal_actual_face_checks=literal_faces,excluded_frontier_types=len(excluded),six_node_relation_rows=96,six_node_relation_separate_states=nodes,scope='Separate table-join traversal, independent finite-tree classifier with isolated-A1 padding, independent maximal-face reconstruction and literal face containment checks. Replays all44 local obstructions and the complete96-row relation without producer search/canonicalization imports. The classical order4/5 Lanner table and source incidence are inherited dependencies. Bounded failed attempts and unconstrained metric realization are outside this negative claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['excluded_frontier_types'],flush=True)
