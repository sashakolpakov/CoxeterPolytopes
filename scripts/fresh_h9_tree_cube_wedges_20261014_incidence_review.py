#!/usr/bin/env python3
"""Independent literal stellar/wedge incidence and every saved identification.

No producer imports and no graph canonicalization. This is not a second
primal-coordinate replay or an independent nonisomorphism classification.
"""
import argparse,json
from itertools import combinations,permutations,product
from collections import Counter
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tree_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 path=f'results/{TAG}.json';d=json.loads((ROOT/path).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 assert d['recipe_count']==240 and d['representative_count']==33
 expected={(name,order,u)for name,edges in d['trees'].items()for order in permutations(map(tuple,edges))for u in range(5)}
 literal={(r['tree'],tuple(map(tuple,r['cut_order'])),r['unwedged_coordinate'])for r in d['recipes']}
 assert literal==expected and len(d['recipes'])==len(literal)
 reps={r['id']:r for r in d['representatives']};counts=Counter();checked_vertices=0
 for rec in d['recipes']:
  base={frozenset(2*i+s for i,s in enumerate(ss))for ss in product([0,1],repeat=5)}
  for k,(i,j)in enumerate(rec['cut_order']):
   e={2*i,2*j};new=set()
   for v in base:
    if e<=v:
     new.update((v-{z})|{10+k}for z in e)
    else:new.add(v)
   base=new
   assert set(map(len,base))=={5}
   assert set(Counter(tuple(sorted(a))for v in base for a in combinations(v,4)).values())=={2}
  faces={0}
  for v in base:
   mask=sum(1<<i for i in v);s=mask
   while s:faces.add(s);s=(s-1)&mask
  missing=[]
  for m in range(1,1<<14):
   if m in faces:continue
   if all(m^(1<<i)in faces for i in range(14)if m>>i&1):missing.append(tuple(i for i in range(14)if m>>i&1))
  assert set(map(len,missing))=={2}
  chosen={2*i+1 for i in range(5)if i!=rec['unwedged_coordinate']};blocks=[];n=0
  for i in range(14):
   w=2 if i in chosen else 1;blocks.append(tuple(range(n,n+w)));n+=w
  V=set()
  for v in base:
   fixed=set(j for i in v for j in blocks[i]);free=sorted(chosen-v)
   for js in product(*(blocks[i]for i in free)):V.add(tuple(sorted(fixed|set(js))))
  M={tuple(sorted(j for i in m for j in blocks[i]))for m in missing}
  assert n==18 and set(map(len,V))=={9}and all(len(m)in(2,3)for m in M)
  g=reps[rec['representative']]['geometry'];pm=rec['facet_isomorphism'];assert sorted(pm)==list(range(18))
  assert{tuple(sorted(pm[i]for i in v))for v in V}=={tuple(z['incident'])for z in g['vertices']}
  assert{tuple(sorted(pm[i]for i in m))for m in M}==set(map(tuple,g['minimal_nonfaces']))
  assert len(V)==rec['vertex_count'];counts[rec['representative']]+=1;checked_vertices+=len(V)
 for k,r in reps.items():
  assert counts[k]==r['recipe_multiplicity']
  order=r['canonical_facet_order'];assert sorted(order)==list(range(18));pos={x:i for i,x in enumerate(order)}
  assert sorted(tuple(sorted(pos[i]for i in m))for m in r['geometry']['minimal_nonfaces'])==list(map(tuple,r['canonical_clutter']))
 paths=[f'scripts/{TAG}.py',path,f'scripts/{TAG}_incidence_review.py']
 return dict(status='PASS_240_LITERAL_RECIPES_AND_EXPLICIT_ISOMORPHISMS',sources={p:sha(p)for p in paths},
  ordered_recipes=240,saved_representatives=33,mapped_complete_vertices=checked_vertices,
  recipe_multiplicities=dict(sorted(counts.items())),
  scope='Independent combinatorial reconstruction of all four stellar cuts and four wedges for each literal recipe; full vertex sets and minimal-nonface clutters match every stored facet bijection. The rational small-cut/ordinary-wedge realization is source/proof reviewed, not arithmetically rerun here. All identifications are checked without pynauty, but the absence of additional isomorphisms among the33 representatives remains the primary canonicalization dependency. No Coxeter label or metric search is duplicated.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_incidence_review.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['ordered_recipes'],d['mapped_complete_vertices'])
