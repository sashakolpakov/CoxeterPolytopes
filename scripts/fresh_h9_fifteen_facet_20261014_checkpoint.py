#!/usr/bin/env python3
"""Disjoint final ledger of the frozen fifteen-facet construction attempt."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_facet_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def run():
 suffixes=['incidence','patterns','search','metric_priority','incidence_independent','search_independent','metric_priority_independent','review']
 paths=[f'results/{TAG}_{s}.json'for s in suffixes];data={s:load(p)for s,p in zip(suffixes,paths)}
 inc=data['incidence'];pat=data['patterns'];searches=[data['search'],data['metric_priority']]
 old={r['case']for r in inc['cases']if r['new_obstruction']is not None}
 local={i for i,p in pat['excluded_cases']};assert not old&local
 done=set();states=0
 for s,a in zip(searches,[data['search_independent'],data['metric_priority_independent']]):
  good={r['case']for r in s['cases']if r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'}
  assert good=={r['case']for r in a['records']}and len(good)==len(s['cases'])
  assert not good&(old|local|done);done|=good;states+=a['states']
 allids=set(range(len(inc['cases'])));left=allids-old-local-done
 assert left==set(pat['remaining_cases'])-done
 assert(len(old),len(local),len(done),len(left),states)==(39931,4240,174,1334,2847)
 assert data['incidence_independent']['counts']['actual_vertices']==9817894
 paths.append(f'scripts/{TAG}_checkpoint.py')
 return dict(status='FROZEN_VERIFIED_PROGRESS_WITH_1334_OPEN_FIFTEEN_FACET_NERVES',sources={p:sha(p)for p in paths},
  incidence_counts=inc['counts'],distinct_nerves=45679,inherited_local_or_block_exclusions=len(old),new_local_pattern_exclusions=len(local),
  complete_global_label_exclusions=len(done),global_states=states,remaining_without_complete_global_search=len(left),
  excluded_by_global_search=sorted(done),open_case_ids=sorted(left),full_category_survivors=0,metric_realizations=0,
  limits=dict(local_patterns_tried=30,local_state_limit=10000,local_label_survivors=19,local_node_limits=9,
   initial_global_cases=50,initial_state_limit=3000,diverse_parent_cases=100,diverse_parent_state_limit=2000,
   metric_priority_cases=24,metric_priority_state_limit=5000,global_node_limited_cases=0),
  scope='Complete codimension2..5 one-cut incidence coverage of the770 saved realizable14-facet parents. Exact exclusions cover44,345 distinct resulting nerves. The1,334 listed cases have no completed global label search and remain open. The30 local trials and174 global selections are bounded, and no full labeling reached metric fitting. This is not all15-facet combinatorics and proves no general compact-H9 impossibility.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}_checkpoint.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['remaining_without_complete_global_search'],flush=True)
