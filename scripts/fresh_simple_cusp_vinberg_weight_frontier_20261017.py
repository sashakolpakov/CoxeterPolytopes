#!/usr/bin/env python3
"""Exact affine 3n-cost interface; quadrilateral payment is explicitly open."""
import argparse
from collections import Counter, deque
import hashlib
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_vinberg_weight_frontier_20261017'
AFFINE='fresh_simple_cusp_affine_exclusion_20261017'
COMPACT='fresh_compact_dimension_bound_20261015'


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def distances(n,edges):
    g=[set() for _ in range(n)]
    for a,b in edges:g[a].add(b);g[b].add(a)
    out={}
    for a in range(n):
        ds={a:0};queue=deque([a])
        while queue:
            x=queue.popleft()
            for y in g[x]:
                if y not in ds:ds[y]=ds[x]+1;queue.append(y)
        assert len(ds)==n
        for b in range(a+1,n):out[a,b]=ds[b]
    return out


def model(kind,n):
    if kind=='cycle':
        return [(i,(i+1)%n) for i in range(n)],[]
    if kind=='path':
        return [(i,i+1) for i in range(n-1)],[list(range(6)),list(range(n-6,n))]
    if kind=='single_fork':
        m=n-2
        return ([(i,i+1) for i in range(m-1)]+[(0,m),(0,m+1)],
                [list(range(4))+[m,m+1],list(range(m-6,m))])
    assert kind=='double_fork'
    m=n-4
    return ([(i,i+1) for i in range(m-1)]+[(0,m),(0,m+1),(m-1,m+2),(m-1,m+3)],
            [list(range(4))+[m,m+1],list(range(m-4,m))+[m+2,m+3]])


def run():
    affine=json.loads((ROOT/f'results/{AFFINE}.json').read_text())
    assert affine['status']=='PASS_EXACT_ALL_SIMPLE_AFFINE_EXCLUSION_DIMENSION61'
    for p,h in affine['sources'].items():assert sha(p)==h,p
    rows=[]
    for n in range(14,101):
        for kind,deficit in [('cycle',0),('path',6),('single_fork',5),('double_fork',4)]:
            edges,ends=model(kind,n);d=distances(n,edges)
            short={p for p,k in d.items() if k<=3}
            endpairs={tuple(sorted(p)) for end in ends for p in combinations(end,2)}
            extra=endpairs-short;weighted=short|endpairs
            assert len(short)==3*n-deficit and len(extra)==deficit and len(weighted)==3*n
            rows.append(dict(type=kind,order=n,short_distance_cost=len(short),
                extra_end_pairs=sorted(extra),weight=len(weighted),budget=3*n))
    envelopes=[]
    for n in range(14,41):
        edges,_=model('cycle',n);d=distances(n,edges);union=set()
        for cut in range(n):
            path=[(cut+j)%n for j in range(1,n)]
            pairs={tuple(sorted((path[j],path[j+1]))) for j in range(n-2)}
            for end in (path[:6],path[-6:]):pairs.update(tuple(sorted(p)) for p in combinations(end,2))
            union.update(pairs)
        expected={p for p,k in d.items() if k<=5}
        assert union==expected and len(union)==5*n
        envelopes.append(dict(order=n,max_spherical_deletion_cost=len(union),three_n=3*n))
    paths=[f'scripts/{TAG}.py',f'scripts/{AFFINE}.py',f'results/{AFFINE}.json',
        f'results/{AFFINE}.md',f'results/{AFFINE}_review.md',
        f'scripts/{COMPACT}.py',f'results/{COMPACT}.json',f'results/{COMPACT}.md',
        f'results/{COMPACT}_review.md']
    return dict(status='PASS_AFFINE_3N_WEIGHT_INTERFACE_QUADRILATERAL_PAYMENT_OPEN',
        sources={p:sha(p) for p in paths},minimum_dimension=14,
        affine_graph_checks=rows,affine_cost_formula='3n in every classical affine type, n>=14',
        all_rank_short_costs=dict(cycle='3n',path='3n-6',single_fork='3n-5',double_fork='3n-4'),
        all_rank_extra_end_costs=dict(cycle=0,path=6,single_fork=5,double_fork=4),
        cycle_deletion_envelope_checks=envelopes,
        triangle_payment='At least2: inherited literal Lanner cases; newly affine corners use distance<=3.',
        quadrilateral_payment='OPEN: full affine corner cycles invalidate the old elliptic tree decomposition.',
        new_ambient_dimension_bound=None,
        scope='Primary-literature scope audit plus exact all-rank affine budget and triangle-transfer proof interface. Not a proof of all-simple29 or31. The maximum spherical-deletion envelope costs5n on affineA and does not establish face payment. Frozen61 remains valid.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['affine_graph_checks']),'affine checks;',
          len(result['cycle_deletion_envelope_checks']),'cycle envelopes; no new bound')
