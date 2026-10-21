#!/usr/bin/env python3
"""Disjoint proof ledger for the fully closed5536-type changed-order family."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_checkpoint'
FILES={
'incidence':'results/fresh_h9_bugaenko_remaining_20261014_newfacet.json',
'old_closure':'results/fresh_h9_bugaenko_remaining_20261014_closure.json',
'old_block_audit':'results/fresh_h9_bugaenko_remaining_20261014_block_independent.json',
'pilot':'results/fresh_h9_bugaenko_remaining_20261014_newfacet_search.json',
'pilot_audit':'results/fresh_h9_bugaenko_remaining_20261014_newfacet_audit.json',
'patterns':'results/fresh_h9_newfacet_frontier_20261014_patterns.json',
'relations':'results/fresh_h9_newfacet_frontier_20261014_relations.json',
'local_audit':'results/fresh_h9_newfacet_frontier_20261014_local_independent.json',
'incidence_audit':'results/fresh_h9_newfacet_frontier_20261014_incidence_independent.json',
'search':'results/fresh_h9_newfacet_frontier_20261014_search.json',
'search_audit':'results/fresh_h9_newfacet_frontier_20261014_search_independent.json'}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 data={k:json.loads((ROOT/p).read_text())for k,p in FILES.items()}
 for d in data.values():
  for p,h in d['sources'].items():assert sha(p)==h,p
 rows=data['incidence']['cases'];assert len(rows)==5536
 classes={
  'already_excluded_old_catalogue':{r['case']for r in rows if r['old_catalogue_case']is not None},
  'additional_two_block_gate':{r['case']for r in rows if r['old_catalogue_case']is None and r['block_obstruction']is not None},
  'previous_complete_pilot':{r['case']for r in data['pilot']['cases']},
  'new_induced_nerve_obstructions':set(data['relations']['combined_excluded_cases']),
  'new_complete_global_searches':{r['case']for r in data['search']['cases']}}
 assert [len(v)for v in classes.values()]==[1040,1220,32,3112,132]
 used=set()
 for ids in classes.values():assert not used&ids;used|=ids
 assert used==set(range(5536))
 assert all(r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'for r in data['pilot']['cases']+data['search']['cases'])
 assert not data['search']['metric_scouts']
 assert set(data['relations']['remaining_case_ids'])==classes['new_complete_global_searches']
 assert data['local_audit']['excluded_frontier_types']==3112
 assert data['incidence_audit']['cases']==5536
 assert {r['case']for r in data['search_audit']['records']}==classes['new_complete_global_searches']
 assert data['search_audit']['total_states']==1465 and data['search_audit']['maximum_states']==67
 files=list(FILES.values())+[f'scripts/{TAG}.py']
 return dict(status='PASS_ALL5536_CHANGED_ORDER_COMPACT_H9_TYPES_EXCLUDED',sources={p:sha(p)for p in files},
  family='Truncate a codimension2..5 face of the fixed11-facet Bugaenko H8; wedge the new truncation facet; truncate a codimension2..5 face of that H9 polytope.',
  operation_representatives=11729,large_missing_face_failures=5503,admissible_recipes=6226,admissible_incidence_types=5536,
  proof_counts={k:len(v)for k,v in classes.items()},proof_classes={k:sorted(v)for k,v in classes.items()},
  new_local_obstructions=44,new_six_node_relation_rows=96,new_complete_search_states=1465,new_maximum_search_states=67,
  unresolved_types=0,metric_survivors=0,finite_volume_examples=0,
  scope='Complete incidence exclusion for this specified changed-order family. No label or dotted-distance bound enters the final exclusions. The160-trial local discovery included80 limited probes, but every final catalogue type has a separate complete proof path. No claim covers arbitrary15+-facet surgeries, a different source polytope, or all compactH9 polytopes. A complete label survivor would require metric solving; none survives.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True)
