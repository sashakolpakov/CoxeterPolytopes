#!/usr/bin/env python3
"""Exact affine graph-cost lemma for distance8 plus feature-near distance9.

The all-rank proof uses explicit spine intervals, not a finite-rank cutoff.
This package does not prove triangle/quad payment or a dimension bound.
"""
import argparse
from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_private_affine_tree_exclusion_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_distance8_feature_payment_20261017'
RADIUS=14


def graph_data(m,edges):
    adj=[set()for _ in range(m)]
    for a,b,label in edges:
        assert a!=b and b not in adj[a]
        adj[a].add(b);adj[b].add(a)
    feature={i for i in range(m)if len(adj[i])>=3}
    for a,b,label in edges:
        if label=='parallel' or label>=4:feature.update((a,b))
    distance=[]
    for a in range(m):
        ds={a:0};queue=deque([a])
        while queue:
            b=queue.popleft()
            for c in sorted(adj[b]):
                if c not in ds:ds[c]=ds[b]+1;queue.append(c)
        assert len(ds)==m
        distance.append([ds[b]for b in range(m)])
    return adj,feature,distance


def payment(kind,m,edges,radius=RADIUS,whole_small=False):
    adj,feature,d=graph_data(m,edges)
    near={i for i in range(m)if feature and min(d[i][j]for j in feature)<=radius}
    histogram=Counter();short=[];extra=[]
    for a,b in combinations(range(m),2):
        distance=d[a][b];histogram[distance]+=1
        if distance<=8:short.append((a,b))
        elif distance==9 and (a in near or b in near):extra.append((a,b))
    assert not(set(short)&set(extra))
    count=len(short)+len(extra)
    paid=m*(m-1)//2 if whole_small and m<=17 else count
    spine=m if kind in('A','C')else m-2 if kind=='B'else m-4 if kind=='D'else None
    classes=Counter()
    if spine is not None:
        for a,b in extra:
            classes['spine_spine'if a<spine and b<spine else'leaf_leaf'if a>=spine and b>=spine else'leaf_spine']+=1
    return dict(kind=kind,order=m,feature_vertices=sorted(feature),near_feature_vertices=sorted(near),
        pair_distance_histogram=sorted(histogram.items()),distance_through8=len(short),
        feature_near_distance9=len(extra),extra_categories=dict(sorted(classes.items())),
        unaugmented_cost=count,whole_small_component_enabled=whole_small,paid_pairs=paid,
        eight_m=8*m,all_pairs=m*(m-1)//2,budget=min(8*m,m*(m-1)//2),radius=radius)


def classical_formula(kind,m):
    assert m>=18
    if kind=='A':return 8*m,0,dict()
    if kind=='C':return 8*m-36,min(m-9,32),dict(spine_spine=min(m-9,32))
    if kind=='B':return 8*m-35,2+min(m-11,31),dict(leaf_spine=2,spine_spine=min(m-11,31))
    assert kind=='D'
    return 8*m-34,4+min(m-13,30),dict(leaf_spine=4,spine_spine=min(m-13,30))


def small_catalogue():
    yield 'A1',2,[(0,1,'parallel')]
    for kind,lo in [('A',3),('B',4),('C',3),('D',5)]:
        for m in range(lo,18):yield kind,m,P.graph(kind,m)[0]
    for kind,arms in [('E6',(2,2,2)),('E7',(1,3,3)),('E8',(1,2,5))]:
        m,es=P.arms(arms);yield kind,m,es
    yield 'F4',5,[(0,1,3),(1,2,3),(2,3,4),(3,4,3)]
    yield 'G2',3,[(0,1,3),(1,2,6)]


def interval_union_count(length,left,right):
    """Two endpoint intervals covering left/right slots among length slots."""
    assert length>=0 and left>=0 and right>=0
    x=set(range(min(length,left)))
    y=set(range(max(0,length-right),length))
    assert len(x|y)==min(length,left+right)
    return len(x|y)


def run():
    # Exact finite checks corroborate the separately stated all-rank spine proof.
    rows=[]
    for m in range(18,129):
        for kind in 'ABCD':
            es,_=P.graph(kind,m)
            row=payment(kind,m,es)
            short,extra,categories=classical_formula(kind,m)
            assert row['distance_through8']==short
            assert row['feature_near_distance9']==extra
            assert row['extra_categories']==categories
            assert row['paid_pairs']==short+extra<=row['budget']
            if kind=='A':assert row['feature_vertices']==[] and row['paid_pairs']==8*m
            if kind=='C':assert row['feature_vertices']==[0,1,m-2,m-1]
            if kind=='B':assert row['feature_vertices']==[0,m-4,m-3]
            if kind=='D':assert row['feature_vertices']==[0,m-5]
            rows.append(row)
    small=[]
    for kind,m,es in small_catalogue():
        row=payment(kind,m,es,whole_small=True)
        assert row['unaugmented_cost']<=row['paid_pairs']==row['all_pairs']==row['budget']
        row['edges']=es;small.append(row)
    assert len(small)==63 and len(rows)==444
    # Independent interval enumeration checks exactly the counts appearing in
    # the proof. Arbitrary length follows by the explicit endpoint-union formula.
    intervals=[]
    for length in range(0,150):
        for kind,left,right in [('C',16,16),('B',15,16),('D',15,15)]:
            intervals.append([kind,length,interval_union_count(length,left,right)])
    # Algebraic branch certificates: deficit from8m is max(small branch,limit).
    deficits={kind:dict(first_deficit_constant=c,limiting_deficit=l,
        formula=f'max({c}-m,{l})',saturation_order=c-l)
        for kind,c,l in [('C',45,4),('B',44,2),('D',43,0)]}
    for row in rows:
        kind,m=row['kind'],row['order']
        if kind in deficits:
            z=deficits[kind];assert 8*m-row['paid_pairs']==max(z['first_deficit_constant']-m,z['limiting_deficit'])
    # Optional sharpness diagnostic: radius15 already exceeds8m in affineD.
    sharp=payment('D',100,P.graph('D',100)[0],radius=15)
    assert sharp['paid_pairs']==8*100+2
    paths=[Path(__file__),Path(P.__file__),
        ROOT/'scripts/fresh_simple_cusp_distance_bound_20261017.py',
        ROOT/'results/fresh_simple_cusp_distance_bound_20261017.json',
        ROOT/'results/fresh_simple_cusp_distance_bound_20261017.md',
        ROOT/'scripts/fresh_h20_parabolic_endpoint_profiles_20261017.py']
    out=dict(status='PASS_ALL_RANK_AFFINE_DISTANCE8_FEATURE_PAIR_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        rule=dict(short_distance=8,extra_distance=9,feature_neighborhood_radius=14,
            feature_vertices='degree>=3 or endpoint of a Coxeter bond with label>=4 (parallel bond also counts)',
            extra_condition='at least one endpoint is within14 graph steps of the feature set',
            optional_whole_component_order_at_most=17),
        conclusion='For every connected affine Coxeter diagram of order m, paid unordered pairs<=min(8m,m(m-1)/2).',
        classical_checks=rows,small_affine_catalogue=small,interval_union_checks=intervals,
        all_rank_formulas=dict(A=dict(short='8m',extra='0',total='8m'),
            C=dict(short='8m-36',extra='min(m-9,32)',total='min(9m-45,8m-4)'),
            B=dict(short='8m-35',extra='2+min(m-11,31)',total='min(9m-44,8m-2)'),
            D=dict(short='8m-34',extra='4+min(m-13,30)',total='min(9m-43,8m)')),
        symbolic_deficits=deficits,radius15_failure=sharp,
        counts=dict(classical_bfs_checks=len(rows),small_affine_checks=len(small),endpoint_interval_checks=len(intervals)),
        affine_classification_input='The standard connected affine Coxeter classification: cycles and classical B/C/D trees in unbounded order; remaining exceptions have order<=9. This classification is inherited, not newly re-enumerated.',
        graph_cost_lemma_complete=True,face_payment_proved_here=False,ambient_dimension_bound=None,
        scope='All-rank affine-factor graph cost only. Finite-vertex payment and every triangle/quadrilateral paying the proposed rule require separate arguments. No97 or other new dimension theorem is asserted.')
    return json.loads(json.dumps(out))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
