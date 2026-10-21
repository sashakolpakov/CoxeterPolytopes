#!/usr/bin/env python3
"""Compose the independently reviewed T10 metric-family exclusion with the old gates."""
import argparse
import json
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tremaining_20261015_case10_metric_gates'
OLD='results/fresh_h9_tremaining_20261015_gates_case_10.json'
THEOREM='results/fresh_h9_tremaining_20261015_case10_rank_gate.json'
REVIEW='results/fresh_h9_tremaining_20261015_case10_rank_gate_review.json'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json'


def run():
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    seen=set()
    def read(p):
        d=json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q,h in d.get('sources',{}).items():
                assert sha(q)==h,q
                if q.endswith('.json'):read(q)
        return d
    old,theorem,review,geo,state=map(read,(OLD,THEOREM,REVIEW,GEOM,STATE))
    assert old['case']==10 and old['total_predicates']==151
    assert theorem['status']=='PASS_T10_FIRST_MANDATORY_FAMILY_ALL_HIGH_LABELS_EXCLUDED'
    assert review['status'].startswith('PASS')and review['sources'][THEOREM]==sha(THEOREM)
    assert theorem['fixed_mandatory_family']==state['solution']['pair_labels']
    g=geo['representatives'][10]['geometry'];vertices={tuple(v['incident'])for v in g['vertices']}
    missing={tuple(m)for m in g['minimal_nonfaces']};blocks=[b for b in g['blocks']if len(b)==2]
    assert blocks==[[3,4],[6,7],[9,10],[12,13]]
    unique={};maps=[]
    for bits in product([0,1],repeat=4):
        p=list(range(18))
        for bit,(a,b)in zip(bits,blocks):
            if bit:p[a],p[b]=b,a
        assert {tuple(sorted(p[i]for i in v))for v in vertices}==vertices
        assert {tuple(sorted(p[i]for i in m))for m in missing}==missing
        fixed=[];bounds=[]
        for a,b,m in theorem['fixed_mandatory_family']:
            x,y=sorted((p[a],p[b]));(bounds if m==7 else fixed).append([x,y,m])
        fixed.sort();bounds.sort();assert len(bounds)==1 and bounds[0][2]==7
        key=(tuple(map(tuple,fixed)),tuple(map(tuple,bounds)))
        maps.append(dict(swaps=list(bits),permutation=p,fixed=fixed,bounds=bounds))
        unique.setdefault(key,len(unique))
    families=list(old['families'])
    for (fixed,bounds),i in unique.items():
        families.append(dict(id=f'case10_exact_residual_family_{i}',
                             origin='independently_reviewed_T10_rank_and_inertia_gate',
                             fixed=list(map(list,fixed)),bounds=list(map(list,bounds))))
    assert len(unique)==2 and len(families)==153
    arr={(a,b):m for a,b,m in state['solution']['pair_labels']}
    hits=[f['id']for f in families[151:]if all(arr.get((a,b))==m for a,b,m in f['fixed'])
          and all(arr.get((a,b),0)>=m for a,b,m in f['bounds'])]
    assert len(hits)==1
    paths=[OLD,THEOREM,REVIEW,GEOM,STATE,f'scripts/{TAG}.py']
    return dict(status='PASS_TNEXT_LITERAL_REVIEWED_CONDITIONS',package=TAG,case=10,
                sources={p:sha(p)for p in paths},vertices=len(vertices),
                inherited_predicates=151,new_metric_predicates=2,total_predicates=153,
                families=families,literal_twin_transports=maps,saved_survivor_hits=hits,
                scope='Adds only the independently reviewed full saved-family metric exclusion and its literal root-swap image. All16 root-pair swaps preserve the complete444-vertex nerve and minimal-nonface clutter. Category7 still means every integer>=7; the composed theorem covers its whole tail. This binder does not itself exclude T10 or consume any pending DFS frame.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['total_predicates'],d['saved_survivor_hits'])
