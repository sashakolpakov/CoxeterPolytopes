#!/usr/bin/env python3
"""Bound the exact incidence family and separate all search outcomes."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tree_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for q,h in d.get('sources',{}).items():assert sha(q)==h,q
 return d
def run():
 paths=[f'results/{TAG}.json',f'results/{TAG}_primal.json',f'results/{TAG}_label_replay.json',f'scripts/{TAG}_summary.py']
 g=read(paths[0]);primal=read(paths[1]);replay=read(paths[2]);rows=[]
 assert g['recipe_count']==240 and g['representative_count']==33
 assert len(primal['cases'])==33 and primal['total_vertices']==14988
 for r in g['representatives']:
  p=f'results/{TAG}_case_{r["id"]:02d}.json';d=read(p);paths.append(p)
  assert d['case']==r['id'] and d['recipe_multiplicity']==r['recipe_multiplicity']
  assert (d['solution']is not None)==(d['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR')
  assert (d['unexplored_frontier']is not None)==(d['status']=='NODE_LIMIT_UNRESOLVED')
  rows.append(dict(case=r['id'],tree=r['tree'],cut_order=r['cut_order'],unwedged_coordinate=r['unwedged_coordinate'],
   recipe_multiplicity=r['recipe_multiplicity'],vertices=len(r['geometry']['vertices']),
   missing_size_counts=r['geometry']['missing_size_counts'],status=d['status'],nodes=d['nodes'],
   node_limit=d['settings']['node_limit'],partial_survivor=d['solution']))
 closed=[r['case']for r in rows if r['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY']
 survived=[r['case']for r in rows if r['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR']
 limited=[r['case']for r in rows if r['status']=='NODE_LIMIT_UNRESOLVED']
 assert closed==[1,2,3,4] and survived==[0,8,17,22,25,26] and len(limited)==23
 assert {r['case']for r in replay['completed_cases']}==set(closed+survived)
 count=Counter(r['status']for r in rows)
 weighted={s:sum(r['recipe_multiplicity']for r in rows if r['status']==s)for s in sorted(count)}
 assert sum(weighted.values())==240
 return dict(status='EXACT_NONPATH_TREE_INCIDENCES_WITH_BOUNDED_LABEL_FRONTIER',sources={p:sha(p)for p in paths},
  representative_counts=dict(count),recipe_counts=weighted,total_search_nodes=sum(r['nodes']for r in rows),
  excluded_cases=closed,necessary_partial_survivors=survived,unresolved_limited_cases=limited,cases=rows,
  metric_status='No hyperbolic metric or compact Coxeter polytope has been constructed. The six survivors assign only mandatory triangle labels; omitted meeting labels and every missing-pair distance remain to be solved, together with rank10/signature(9,1).',
  scope='Exactly the 240 ordered star/T recipes in the geometry artifact. Four entire incidence classes (96 recipes) are excluded for every finite label and every strict dotted distance. Six classes (52 recipes) have first necessary partial-label survivors. Twenty-three classes (92 recipes) hit explicit node limits and remain unresolved. No limit is a negative result. Separate exact primal and finite-positivity replays are bound.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}_summary.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['representative_counts'],flush=True)
