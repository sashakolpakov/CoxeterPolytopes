#!/usr/bin/env python3
"""Bind literal star0 theorem conditions; no new analytic obstruction claimed."""
import argparse,json
from hashlib import sha256
from pathlib import Path
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_gates'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
TARGETS='results/fresh_h9_dotted_central_20261014_next_targets.json'
ORI='results/fresh_h9_star_cube_metric_20261014_orientations.json'
ORIR='results/fresh_h9_star_cube_metric_20261014_orientations_review.json'
FIRST='results/fresh_h9_star_cube_metric_20261014.json'
NEXT='results/fresh_h9_star_cube_metric_20261014_next_metric.json'
STATE='results/fresh_h9_star_cube_metric_20261014_orientation_search_case_00.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 paths=[GEOM,TARGETS,ORI,ORIR,FIRST,NEXT,STATE]
 seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,h in d.get('sources',{}).items():
    assert sha(q)==h,q
    if q.endswith('.json'):read(q)
  return d
 geo,targets,ori,orir,first,nxt,state=[read(p)for p in paths]
 g=geo['representatives'][0]['geometry'];V=[set(v['incident'])for v in g['vertices']]
 ms={frozenset(m)for m in g['minimal_nonfaces']}
 assert len(V)==378
 assert {tuple(sorted(v))for v in V}=={v for v in combinations(range(18),9)if not any(m<=set(v)for m in ms)}
 target=targets['cases'][0];assert target['case']==0 and target['total_unique_conditions']==132 and not target['closed']
 families=[];counts={}
 for kind,record in target['gates'].items():
  for row in record['maps']:
   assert all(frozenset(t)in ms for t in row['missing_own_triples'])
   for s in row['supports']+row['threshold_pair_witnesses']:assert set(s['support'])<=V[s['actual_vertex']]
  for i,row in enumerate(record['conditions']):
   assert row['bounds']==record['maps'][row['witness_map']]['bounds']
   families.append(dict(id=kind+'_'+str(i),origin=kind,fixed=[],bounds=row['bounds']))
  counts[kind]=len(record['conditions'])
 assert counts==dict(dense_high=24,strong_contact=72,two_dense=36)
 assert ori['status']=='PASS_ALL_SIXTEEN_ORIENTATIONS_EXCLUDED'
 assert len(ori['orientations'])==16
 S=first['fixed_positive_support'];P=first['positive_leaves'];C=first['cut_leaves']
 assert S==[3,4,6,7,9,10,12,13] and P==[2,5,8,11] and C==[14,15,16,17]
 for r in ori['all_leaf_both_orientation_budget_guards']:assert set(r['support'])<=V[r['actual_vertex']]
 for o in ori['orientations']:
  assert o['status']=='EXACT_ORIENTATION_EXCLUDED'
  fixed=[S[2*i:2*i+2]+[4]for i in range(4)]
  for i,leaf in enumerate(o['fixed_33_leaves']):
   assert leaf==(C[i]if o['orientation']&(1<<i)else P[i])
   fixed += [sorted((leaf,S[2*i+j]))+[3]for j in range(2)]
  families.append(dict(id='orientation_'+str(o['orientation']),origin='reviewed_star_orientations',fixed=sorted(fixed),bounds=[]))
 assert nxt['status']=='PASS_NEXT_ORIENTATION_ASSIGNMENT_RANK_ELEVEN' and not nxt['surviving_rank_cases']
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 families.append(dict(id='previous_next_assignment',origin='frozen_exact_rank_gate',fixed=state['solution']['pair_labels'],bounds=[]))
 paths += [f'scripts/{TAG}.py']
 return dict(status='PASS_LITERAL_STAR_ZERO_CONDITIONS_BOUND',sources={p:sha(p)for p in paths},
  recursive_manifests=len(seen),case=0,vertices=378,block_counts=counts,orientation_count=16,
  exact_assignment_count=1,families=families,
  scope='132 literal block predicates,16 source-bound orientation families and one previously exact-excluded assignment. Analytic proofs are inherited, not newly proved; omitted contacts and dotted distances remain free.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
