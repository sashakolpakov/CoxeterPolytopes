#!/usr/bin/env python3
"""Exact omitted-contact domains at the first saved T10 necessary survivor.

Category7 remains every finite integer label>=7. This is a necessary-domain
screen, not a metric or whole-incidence exclusion. Dotted distances stay free.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_tremaining_20261015_search as runner

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tremaining_20261015_case10_domains'
STATE='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json'


def run():
    E=runner.engine
    env=E.setup(10)
    saved=E.read_checked(STATE)
    assert saved['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR' and saved['case']==10
    pairs=env['pairs'];ix={p:i for i,p in enumerate(pairs)}
    assignment=[0]*153
    for a,b,m in saved['solution']['pair_labels']:
        assignment[ix[a,b]]=m
    geometry=env['g'];vertices=[r['incident'] for r in geometry['vertices']]
    vertex_edges=[tuple(ix[p] for p in combinations(v,2))for v in vertices]
    masks=[tuple((i,assignment[i])for i in es if assignment[i]>2)for es in vertex_edges]
    assert all(env['finite'](m)for m in masks)
    touching=[tuple(j for j,es in enumerate(vertex_edges)if i in es)for i in range(153)]
    missing_pairs={tuple(m) for m in geometry['minimal_nonfaces']if len(m)==2}
    assert {p for p,t in zip(pairs,touching)if not t}==missing_pairs
    rows=[]
    for i,pair in enumerate(pairs):
        if assignment[i] or not touching[i]:continue
        allowed=[2];failures=[]
        for label in range(3,8):
            bad=next((j for j in touching[i]if not env['finite'](tuple(sorted(masks[j]+((i,label),))))),None)
            if bad is None:allowed.append(label)
            else:failures.append(dict(label=label,actual_vertex=bad,support=vertices[bad]))
        assert allowed==list(range(2,max(allowed)+1))
        rows.append(dict(pair=list(pair),domain=allowed,failures=failures))
    forced=[r['pair']for r in rows if r['domain']==[2]]
    optional=[r for r in rows if len(r['domain'])>1]
    blocks=[b for b in geometry['blocks']if len(b)==2]
    cross=[list(p)for p in combinations(sum(blocks,[]),2)if list(p)not in blocks]
    assert len(cross)==24
    cross_right=all(p in forced or assignment[ix[tuple(p)]]==2 for p in cross)
    witnesses=[dict(pair=p,witness=next(r['failures'][0]for r in rows if r['pair']==p))for p in cross if p in forced]
    paths=[STATE,runner.engine.GATE.format(case=10),
           'scripts/fresh_h9_tremaining_20261015_search.py',
           'scripts/fresh_h9_tnext_20261015_search.py',
           'scripts/independent_h9_14_record30_a1_labels.py',f'scripts/{TAG}.py']
    return dict(status='PASS_EXACT_T10_SAVED_FAMILY_CONTACT_DOMAINS',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                vertices=len(vertices),fixed_meeting_pairs=saved['solution']['pair_labels'],
                dotted_pairs=sorted(map(list,missing_pairs)),unknown_meeting_pairs=len(rows),
                forced_right_pairs=forced,optional_meeting_pairs=optional,
                all_domains=rows,block_roots=blocks,cross_block_right=cross_right,
                cross_block_guards=witnesses,
                scope='Exact necessary finite-vertex contact domains for this saved mandatory-label '
                      'family, including its unrestricted >=7 block. Unknown meeting contacts are '
                      'otherwise relaxed to right. This does not jointly realize the domains or solve '
                      'rank/signature; all dotted distances remain arbitrary. No continuation is run.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'forced',len(d['forced_right_pairs']),'optional',len(d['optional_meeting_pairs']),
          'cross_right',d['cross_block_right'])
