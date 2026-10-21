#!/usr/bin/env python3
"""Freeze the bounded remaining-T search scope without claiming an empty incidence."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tremaining_20261015_ledger'
PREFIX='fresh_h9_tremaining_20261015'


def run():
    seen=set();references=0
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    def read(p):
        nonlocal references
        d=json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q,h in d.get('sources',{}).items():
                assert sha(q)==h,q;references+=1
                if q.endswith('.json'):read(q)
        return d
    priorpath='results/fresh_h9_tnext_20261015_ledger.json';prior=read(priorpath)
    paths=[priorpath,f'scripts/{TAG}.py']
    prepared=[]
    for case in [10,19,29]:
        path=f'results/{PREFIX}_gates_case_{case:02d}.json';paths.append(path);d=read(path)
        assert case in prior['remaining_representatives']and d['case']==case
        prepared.append(dict(case=case,vertices=d['vertices'],inherited_predicates=d['inherited_predicates'],
                             three_large_predicates=d['three_large_predicates'],two_high_predicates=d['two_high_predicates'],
                             initial_total_predicates=d['total_predicates'],searched=case==10))
    p1=f'results/{PREFIX}_search_case_10_batch_001.json'
    p2=f'results/{PREFIX}_case10_metric_search_case_10_batch_002.json'
    p3=f'results/{PREFIX}_case10_metric_gates.json'
    p4=f'results/{PREFIX}_case10_rank_gate.json'
    p5=f'results/{PREFIX}_case10_rank_gate_review.json'
    p6=f'results/{PREFIX}_case10_high_bound_review.json'
    p7=f'results/{PREFIX}_case10_batch002_timing.json'
    paths += [p1,p2,p3,p4,p5,p6,p7]
    first,last,binder,gate,review,highreview,timing=map(read,(p1,p2,p3,p4,p5,p6,p7))
    assert first['status']==last['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert first['new_nodes']==first['cumulative_additional_nodes']==9008
    assert last['new_nodes']==14 and last['cumulative_additional_nodes']==9022
    assert len(last['pending_frames'])==10 and last['settings']['resume_input']==p1
    assert last['settings']['node_budget']==100000
    assert binder['total_predicates']==153 and binder['new_metric_predicates']==2
    assert gate['status']=='PASS_T10_FIRST_MANDATORY_FAMILY_ALL_HIGH_LABELS_EXCLUDED'
    assert review['status']=='PASS_INDEPENDENT_T10_RANK_GATE_REVIEW'and highreview['status'].startswith('PASS')
    a={(x,y):m for x,y,m in first['solution']['pair_labels']}
    b={(x,y):m for x,y,m in last['solution']['pair_labels']}
    changed=[dict(pair=list(p),old=a.get(p),new=m)for p,m in b.items()if a.get(p)!=m]
    assert changed==[dict(pair=[12,13],old=4,new=5)]
    assert timing['exit_code']==0 and timing['output_sha256']==sha(p2)
    assert not any(all(b.get((x,y))==m for x,y,m in f['fixed'])and
                   all(b.get((x,y),0)>=m for x,y,m in f['bounds'])for f in binder['families'])
    return dict(status='PASS_BOUNDED_T10_METRIC_FAMILY_PROGRESS_WITH_OPEN_CURSOR',
                sources={p:sha(p)for p in paths},prepared_cases=prepared,
                new_closed_representatives=[],closed_representatives=prior['all_closed_representatives'],
                remaining_representatives=prior['remaining_representatives'],
                closed_recipes=prior['closed_recipes'],remaining_recipes=prior['remaining_recipes'],total_recipes=240,
                first_survivor_nodes=9008,first_family_exactly_excluded=True,
                independent_family_review=True,metric_predicates=153,
                resumed_nodes=14,cumulative_nodes=9022,pending_frames=10,
                configured_new_node_budget=100000,unspent_new_node_budget=99986,
                stopped_at_first_outside_survivor=True,elapsed_seconds=timing['elapsed_seconds'],
                next_mandatory_family=last['solution']['pair_labels'],changed_labels=changed,
                next_internal_labels=[4,4,7,5],symbolic_category7='Every finite integer m>=7',
                unsearched_prepared_cases=[19,29],checked_json_manifests=len(seen),checked_hash_references=references,
                scope='No new whole-incidence exclusion or metric realization is claimed. T10 remains among18 open representatives. The first saved mandatory family is independently excluded for all high labels and omitted contacts; the next necessary family differs by Dlabel4to5 and remains open. The resumable DFS has preserved all pending siblings. No empty-start exclusion traversal is claimed for this nonempty search, and cases19/29 were prepared but not searched. The prior15 representative/164recipe exclusions are inherited without modification.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['cumulative_nodes'],d['remaining_representatives'])
