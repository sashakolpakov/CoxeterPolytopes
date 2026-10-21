#!/usr/bin/env python3
"""Full-face replay of the fifteen-facet one-truncation family.

No new incidence producer, clutter formula or canonicalization is imported.
Literal operation coverage is checked using verified permutations. Only the
209 cross-recipe isomorphisms use the same external nauty library.
"""
import argparse,json,time
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import pynauty
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import downward,clutter
from fresh_h9_newfacet_frontier_20261014_local_independent import faces,solve
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_facet_20261014'
SRC='results/fresh_h9_gale_frontier_20261014_geometry.json'
INC=f'results/{TAG}_incidence.json';PAT=f'results/{TAG}_patterns.json'
OLD='results/fresh_h9_newfacet_frontier_20261014_patterns.json'
REL='results/fresh_h9_newfacet_frontier_20261014_relations.json'
OLD_CHECK='results/fresh_h9_newfacet_frontier_20261014_local_independent.json'
BLOCK='results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def bits(m):return tuple(i for i in range(m.bit_length())if m>>i&1)
def mask(a):return sum(1<<i for i in a)
def move(m,p):return mask(p[i]for i in bits(m))
def full_faces(vv):
 out=0
 for v in vv:out|=downward(v)
 return out
def reconstruct(old,F):
 f=set(bits(F));vv=set()
 for raw in old['facets']:
  v=set(bits(raw))
  if f<=v:
   for i in f:vv.add(mask((v-{i})|{14}))
  else:vv.add(raw)
 vv=sorted(vv);assert all(v.bit_count()==9 for v in vv)
 mm=clutter(vv,15);assert max(map(int.bit_count,mm))<=5
 assert mask(i for i in range(15)if any(v>>i&1 for v in vv))==32767
 return mm,vv
def graph(mm):
 adj={i:[]for i in range(15+len(mm))};colors=[set(range(15))]
 for k in range(2,6):
  c=set()
  for j,m in enumerate(mm,15):
   if m.bit_count()==k:c.add(j)
  if c:colors.append(c)
 for j,m in enumerate(mm,15):
  for i in bits(m):adj[i].append(j);adj[j].append(i)
 return pynauty.Graph(len(adj),adjacency_dict=adj,vertex_coloring=colors)
@lru_cache(None)
def maximal(n,mm):return tuple(mask(v)for v in faces(n,mm))
def witness(w,mm,ff,patterns):
 u=w['support'];inside={m for m in mm if m&u==m}
 if w['kind']=='two_block':
  a,b=w['blocks'];assert not a&b and a|b==u
  assert sorted((a.bit_count(),b.bit_count()))in([3,4],[3,5],[4,4],[4,5],[5,5])
  assert inside=={a,b}
  tests=[(a^(1<<i))|(b^(1<<j))for i in bits(a)for j in bits(b)]
 else:
  n,small=patterns[w['pattern']];p=w['permutation']
  assert len(p)==len(set(p))==n and mask(p)==u
  assert inside=={move(m,p)for m in small}
  tests=[move(f,p)for f in maximal(n,small)]
 assert all(ff>>f&1 for f in tests)
 return len(tests)
def run():
 start=time.monotonic();src=load(SRC);d=load(INC);pat=load(PAT)
 old=load(OLD);rel=load(REL);load(OLD_CHECK);load(BLOCK)
 patterns={r['pattern']:(r['order'],tuple(r['minimal_nonfaces']))for r in old['patterns']if r['verdict']=='EXACT_LOCAL_OBSTRUCTION'}
 seven=rel['seven_node_obstruction'];patterns['seven']=(seven['order'],tuple(seven['minimal_nonfaces']))
 assert len(src['cases'])==len(d['parents'])==770
 counts=Counter();recipe_keys=set()
 for parent in d['parents']:
  raw=src['cases'][parent['case']];mm=raw['minimal_nonfaces'];vv=raw['facets'];ff=full_faces(vv)
  assert list(clutter(vv,14))==mm
  common=16383
  for w in parent['reducing_obstructions']:
   counts['inherited_actual_face_checks']+=witness(w,mm,ff,patterns);common&=w['support']
  assert common==parent['common_support']
  allfaces={mask(c)for k in range(2,6)for c in combinations(range(14),k)if ff>>mask(c)&1}
  eligible={f for f in allfaces if f&common==f}
  assert len(allfaces)==parent['source_face_count']and sorted(eligible)==parent['eligible_cuts']
  counts['literal_source_faces']+=len(allfaces);counts['eligible_literal_cuts']+=len(eligible)
  counts['old_obstruction_excluded_cuts']+=len(allfaces)-len(eligible)
  gen=parent['verified_automorphism_generators']
  for g in gen:
   assert sorted(g)==list(range(14))
   assert {move(m,g)for m in mm}==set(mm)and {move(v,g)for v in vv}==set(vv)
  reached=set()
  for orb in parent['cut_orbits']:
   expected=set(orb);actual={orb[0]};todo=[orb[0]]
   for f in todo:
    for g in gen:
     h=move(f,g)
     if h not in actual:actual.add(h);todo.append(h)
   assert actual==expected and not actual&reached;reached|=actual
   recipe_keys.add((parent['case'],orb[0],len(orb)))
  assert reached==eligible
  counts['parents_with_eligible_cuts'if eligible else'parents_all_cuts_closed']+=1
 assert recipe_keys=={(r['parent'],r['cut_face'],r['orbit_size'])for r in d['recipes']}
 assert len(recipe_keys)==len(d['recipes'])
 expanded={};canonical_keys=set();new_closed=set()
 for row in d['cases']:
  i=row['case'];assert i==len(expanded)
  mm,vv=reconstruct(src['cases'][row['parent']],row['cut_face'])
  order=row['canonical_permutation'];assert sorted(order)==list(range(15));inv=[order.index(i)for i in range(15)]
  cc=tuple(sorted(move(m,inv)for m in mm));assert list(cc)==row['canonical_minimal_nonfaces']and cc not in canonical_keys
  canonical_keys.add(cc);expanded[i]=(mm,vv);counts['actual_vertices']+=len(vv)
  if row['new_obstruction']:
   counts['new_obstruction_actual_face_checks']+=witness(row['new_obstruction'],mm,full_faces(vv),patterns);new_closed.add(i)
  if i%5000==0:print('incidence',i,'seconds',round(time.monotonic()-start,1),flush=True)
 for r in d['recipes']:
  target=d['cases'][r['case']]
  if(r['parent'],r['cut_face'])==(target['parent'],target['cut_face']):continue
  mm,vv=reconstruct(src['cases'][r['parent']],r['cut_face'])
  assert pynauty.isomorphic(graph(mm),graph(expanded[r['case']][0]))
  counts['duplicate_recipe_isomorphisms']+=1
 counts['operation_orbits']=len(d['recipes']);counts['distinct_fifteen_facet_nerves']=len(expanded)
 counts['new_local_obstruction']=len(new_closed);counts['needs_label_search']=len(expanded)-len(new_closed)
 for k,v in d['counts'].items():assert counts[k]==v,(k,counts[k],v)
 exact=[];closed=set()
 for r in pat['patterns']:
  if r['verdict']!='EXACT_LOCAL_OBSTRUCTION':continue
  n=r['order'];missing=tuple(r['minimal_nonfaces']);solutions,nodes=solve(n,missing);assert not solutions
  assert sorted(maximal(n,missing))==r['maximal_faces']
  for i,p in r['applications']:
   assert len(p)==len(set(p))==n
   row=d['cases'][i];u=mask(p);cc=row['canonical_minimal_nonfaces']
   assert {m for m in cc if m&u==m}=={move(m,p)for m in missing}
   old_order=row['canonical_permutation'];vv=[move(v,[old_order.index(j)for j in range(15)])for v in expanded[i][1]]
   ff=full_faces(vv)
   for f in maximal(n,missing):assert ff>>move(f,p)&1;counts['new_pattern_actual_face_checks']+=1
   assert i not in new_closed;closed.add(i)
  exact.append(dict(pattern=r['pattern'],order=n,minimal_nonfaces=missing,independent_states=nodes,applications=len(r['applications'])))
 assert closed=={i for i,p in pat['excluded_cases']}
 assert set(range(len(expanded)))-new_closed-closed==set(pat['remaining_cases'])
 counts['additional_pattern_excluded']=len(closed);counts['frontier_before_global_search']=len(pat['remaining_cases'])
 paths=[SRC,INC,PAT,OLD,REL,OLD_CHECK,BLOCK,f'scripts/{TAG}_incidence_independent.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py',
  'scripts/fresh_h9_newfacet_frontier_20261014_local_independent.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='PASS_SEPARATE_FULL_FACE_OPERATION_COVERAGE_AND_LOCAL_GATES',sources={p:sha(p)for p in paths},counts=dict(counts),new_exact_patterns=exact,
  scope='Independent set-face truncation and full downward-complex minimal-nonface recovery for every saved nerve; all2,124,486 source faces and every explicit symmetry orbit are checked. Canonical maps cannot merge unrelated cases:209 duplicate recipes receive explicit nauty isomorphism tests. Every old/new gate is checked on the literal induced complex and actual faces. Two new local failures use a separate join traversal and finite classifier. Prior770-source rational realizations,44 local theorems, two-block theorem and classical Lanner tables are source-bound dependencies. The output is a complete incidence audit of this one-cut family, not all fifteen-facet polytopes.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}_incidence_independent.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],flush=True)
