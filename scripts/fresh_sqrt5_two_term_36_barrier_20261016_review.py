#!/usr/bin/env python3
"""Independent literal-graph and two-inequality review of a payment-model barrier.

This checks one model obligation, not its realizability as a global compact
polytope and not the full marked-forest dynamic program.
"""
import argparse
import json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_two_term_36_barrier_20261016_review'
PRIMARY='results/fresh_sqrt5_marked_h4_remainder_free_bonus_20261016.json'


def graph(vertices,edges):
    out={v:set() for v in vertices}
    for a,b,m in edges:
        if a in out and b in out and m!=2:
            assert m==3
            out[a].add(b);out[b].add(a)
    return out


def components(g):
    unseen=set(g);out=[]
    while unseen:
        seed=min(unseen);todo=[seed];seen={seed};unseen.remove(seed)
        while todo:
            v=todo.pop()
            for w in g[v]&unseen:seen.add(w);todo.append(w)
            unseen-=g[v]
        out.append({v:g[v]&seen for v in sorted(seen)})
    return out


def kind(g):
    n=len(g);degrees=sorted(map(len,g.values()))
    assert sum(degrees)==2*(n-1)
    if n>=2 and degrees==[1,1]+[2]*(n-2):return f'A{n}'
    if n==4 and degrees==[1,1,1,3]:return 'D4'
    if n==5 and degrees==[1,1,1,2,3]:return 'D5'
    raise AssertionError((n,degrees))


def distance(g,a,b):
    todo=deque([(a,0)]);seen={a}
    while todo:
        v,d=todo.popleft()
        if v==b:return d
        for w in g[v]-seen:seen.add(w);todo.append((w,d+1))
    raise AssertionError('disconnected')


def pair_isomorphism(g,pair,h,target):
    ids=sorted(g);hs=sorted(h)
    for order in permutations(hs):
        mp=dict(zip(ids,order))
        if {mp[x] for x in pair}!=set(target):continue
        if all((b in g[a])==(mp[b] in h[mp[a]]) for a,b in combinations(ids,2)):
            return [[a,mp[a]] for a in ids]
    raise AssertionError('no marked-pair isomorphism')


def run():
    d=json.loads((ROOT/PRIMARY).read_text());barrier=d['decisive_two_term_barrier']
    cert=d['certificate'];assert F(cert['lower_bound'])==36 and len(cert['terms'])==2
    assert all(F(x)==0 for x in cert['discount_upper_multipliers'])
    quad=next(t for t in cert['terms'] if t['metadata']['kind']=='quad')
    vertex=next(t for t in cert['terms'] if t['metadata']['kind']=='vertex')
    w=barrier['variable_index'];budget=d['variables']
    assert quad['metadata']==dict(kind='quad',family='prior',index=64,external_counts=[1,0,0])
    assert dict((i,F(v)) for i,v in quad['row'])=={w:F(-2)} and F(quad['rhs'])==-1
    assert dict((i,F(v)) for i,v in vertex['row'])=={w:F(8),budget:F(-1)} and F(vertex['rhs'])==-32
    assert F(quad['multiplier'])==4 and F(vertex['multiplier'])==1
    total={}
    for term in (quad,vertex):
        for i,z in term['row']:total[i]=total.get(i,F(0))+F(term['multiplier'])*F(z)
    assert {i:z for i,z in total.items() if z}=={budget:F(-1)}
    assert sum(F(t['multiplier'])*F(t['rhs']) for t in (quad,vertex))==-36
    row=barrier['quadrilateral']['literal_row'];assert row['nodes']==10
    edges=[(a,b,m) for (a,b),m in zip(combinations(range(10),2),row['labels'])]
    assert len(edges)==45 and row['white']==list(range(4,10))
    D5=graph(range(5),[(0,1,3),(0,3,3),(0,4,3),(1,2,3)])
    pair_records=[]
    for a,b in combinations(range(5),2):
        dist=distance(D5,a,b)
        pair_records.append(dict(pair=[a,b],distance=dist,baseline=int(dist<=2),w_coefficient=int(dist==3)))
    assert sum(p['baseline'] for p in pair_records)==8
    unpaid=[p['pair'] for p in pair_records if p['w_coefficient']]
    assert unpaid==[[2,3],[2,4]]
    pair_isomorphism(D5,unpaid[0],D5,unpaid[1])
    for actual,raw in zip(pair_records,barrier['marked_vertex']['D5_pairs']):
        assert actual['pair']==raw['pair'] and F(raw['baseline'])==actual['baseline']
        assert raw['bonus']==('w' if actual['w_coefficient'] else '0')
    corners=[]
    for a,b in product((0,2),(1,3)):
        g=graph(row['white']+[a,b],edges);pieces=components(g)
        types=sorted(kind(p) for p in pieces)
        containing=next(p for p in pieces if a in p)
        if b in containing:
            assert kind(containing)=='D5' and distance(containing,a,b)==3
            iso=pair_isomorphism(containing,[a,b],D5,unpaid[0]);payment=1
        else:
            assert types in (['A4','A4'],['D4','D4'])
            iso=None;payment=0
        assert len(g)+1+12==21
        corners.append(dict(pair=[a,b],components=types,w_coefficient=payment,D5_pair_isomorphism=iso))
    assert sum(c['w_coefficient'] for c in corners)==2
    assert 1+4*len(D5)==21
    baseline=4*sum(p['baseline'] for p in pair_records)
    coefficient=4*sum(p['w_coefficient'] for p in pair_records)
    assert (baseline,coefficient)==(32,8)
    lower=baseline+coefficient*F(1,2);assert lower==36
    paths=[f'scripts/{TAG}.py',PRIMARY]
    return dict(status='PASS_INDEPENDENT_LITERAL_TWO_TERM_WEIGHT_MODEL_BARRIER_36',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                quad_corners=corners,common_spectators=['marked A1','A12'],D5_pairs=pair_records,
                vertex_forest=['marked A1']+['D5']*4,vertex_baseline=baseline,vertex_w_coefficient=coefficient,
                model_lower_bound=str(lower),target=35,
                proof='The literal prior64 quadrilateral has two connected D5 distance3 black pairs and two disconnected A4/A4 or D4/D4 black pairs. The current family pays only the former, both by the same w, forcing2w>=1. In markedA1+4D5, everyD5 has8 baseline pairs and2copies ofthat same marked-pair orbit, yielding32+8w>=36. The two exact dual rows cancelw andgive the identical bound.',
                scope='Independent review of this specific restricted payment-family obstruction and its literal graphs, assuming the stated baseline, absent A4/D4 cross families and rank>=6 discounts. No claim that the model quad is globally realizable, no full DP/facecoverage audit, and no geometric dimension bound. Adding cross payments or changing D5 baseline escapes this obstruction.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
