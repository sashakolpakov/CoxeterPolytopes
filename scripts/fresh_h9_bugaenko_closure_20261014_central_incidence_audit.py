#!/usr/bin/env python3
"""Independent vertex surgeries and explicit isomorphisms for208 recipes."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import networkx as nx
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_closure_20261014_central_incidence_audit'
def support(m):return frozenset(i for i in range(14)if m>>i&1)
def graph(M):
 g=nx.Graph();g.add_nodes_from((i,{'kind':0})for i in range(14))
 for j,m in enumerate(M,14):
  g.add_node(j,kind=1)
  for i in m:g.add_edge(i,j)
 return g
def cut(V,F,new):
 F=set(F);assert 2<=len(F)<=5 and any(F<=v for v in V)
 return {v for v in V if not F<=v}|{frozenset((v-{i})|{new})for v in V if F<=v for i in F}
def run():
 files=['results/fresh_h9_bugaenko_closure_20261014_central.json',
        'results/fresh_h9_bugaenko_closure_20261014_central_closure.json',
        'results/fresh_h9_bugaenko_surgery_20261014_incidence.json',
        'data/known_examples/bugaenko_compact_h8.json']
 d,final,inc,base=[json.loads((ROOT/p).read_text())for p in files]
 for data in[d,final,inc]:
  for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 assert final['status']=='PASS_COMPLETE_CENTRAL_SINGLETON_CARRIER_EXCLUSION'
 old={frozenset(v)for v in base['geometry']['vertices']};wedged=set()
 for v in old:
  if 4 in v:wedged.add(v|{11})
  else:wedged|={v|{4},v|{11}}
 records=[]
 for recipe in d['eligible_recipes']:
  V=cut(cut(wedged,recipe['first_face'],12),recipe['second_face'],13)
  M=[support(m)for m in clutter([sum(1<<i for i in v)for v in V],14)]
  row=inc['cases'][recipe['case']];target=[support(m)for m in row['canonical_minimal_nonfaces']]
  matcher=nx.algorithms.isomorphism.GraphMatcher(graph(M),graph(target),node_match=lambda a,b:a['kind']==b['kind'])
  assert matcher.is_isomorphic();p=[matcher.mapping[i]for i in range(14)]
  assert sorted(p)==list(range(14))
  assert {frozenset(p[i]for i in m)for m in M}==set(target)
  assert {sum(1<<p[i]for i in v)for v in V}==set(row['canonical_vertices'])
  records.append(dict(case=recipe['case'],first_face=recipe['first_face'],second_face=recipe['second_face'],explicit_facet_isomorphism=p,vertices=len(V)))
 assert len(records)==208 and len({r['case']for r in records})==128
 assert {r['case']for r in records}=={r['case']for r in final['exclusion_records']}
 files +=[str(Path(__file__).resolve().relative_to(ROOT)), 'scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py']
 return dict(status='PASS_SEPARATE_CENTRAL_RECIPE_ISOMORPHISM_AUDIT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},recipe_orbits=208,distinct_excluded_types=128,records=records,scope='All208 retained recipes are rebuilt from source vertices; minimal nonfaces come from complete downward complexes. NetworkX supplies an isomorphism, whose literal facet map is then directly checked on both nonfaces and all vertices. No nauty canonicalization is reused. Upstream completeness of the208 safe recipe orbits and each exact label exclusion remain explicitly bound dependencies.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
