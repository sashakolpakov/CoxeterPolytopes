#!/usr/bin/env python3
"""Small exact checks for global affine exclusion in the all-simple class."""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_exclusion_20261017'
PREV='fresh_simple_cusp_distance_bound_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def lower(n):
    return F(n*(n-6),8) if n%2==0 else F((n-1)*(n-5),8)


def distances(g,starts):
    d={i:0 for i in starts};queue=deque(starts)
    while queue:
        i=queue.popleft()
        for j in g[i]:
            if j not in d:d[j]=d[i]+1;queue.append(j)
    return d


def run():
    prev=json.loads((ROOT/f'results/{PREV}.json').read_text())
    assert prev['status']=='PASS_EXACT_DISTANCE15_SIMPLE_CUSP_AND_COMPACT_FACET_BOUNDS'
    for p,h in prev['sources'].items(): assert sha(p)==h,p
    census=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    assert census['status']=='PASS_EXACT_HIGH_RANK_QUASILANNER_ATTACHMENT_CENSUS'
    for p,h in census['sources'].items():assert sha(p)==h,p
    affine_rows=[]
    for row in census['diagrams']:
        assert row['ideal_deletions'] and 6<=row['order']<=10
        for i in row['ideal_deletions']:
            deletion=row['deletions'][i]
            assert len(deletion)==1 and deletion[0]['kind']=='affine'
        affine_rows.append(dict(diagram=row['index'],order=row['order'],
            affine_order=row['order']-1,affine_deletions=row['ideal_deletions']))
    assert len(affine_rows)==26 and max(x['affine_order'] for x in affine_rows)==9
    # This finite confirmation uses ALL connected simple graphs of orders
    # 2..5, not an assumed list of low-order Coxeter diagrams. The all-rank
    # argument for an N-node graph is the elementary radius <= N-2 proof.
    small=[]
    for n in range(2,6):
        edges=list(combinations(range(n),2));connected=0;pairs=0;hist=Counter()
        for mask in range(1<<len(edges)):
            g=[[] for _ in range(n)]
            for k,(a,b) in enumerate(edges):
                if mask>>k&1:g[a].append(b);g[b].append(a)
            if len(distances(g,[0]))!=n:continue
            connected+=1
            for black in combinations(range(n),2):
                radius=max(distances(g,black).values())
                assert radius<=n-2<=3
                hist[radius]+=1;pairs+=1
        small.append(dict(order=n,connected_graphs=connected,marked_pairs=pairs,
                          radius_histogram=sorted(hist.items())))
    graph_rows=[]
    for row in prev['graph_formula_checks']:
        n=row['order'];counts=dict(row['counts'])
        assert all(a<=n for a in counts.values())
        cost=sum(a for k,a in counts.items() if k<=7)
        assert cost<=7*n
        graph_rows.append(dict(type=row['type'],order=n,cost=cost,budget=7*n))
    assert len(graph_rows)==245
    exceptions=[]
    for row in prev['exceptional_graph_checks']:
        n=row['order'];total=sum(a for _,a in row['counts'])
        assert n<=9 and total==n*(n-1)//2<=4*n<=7*n
        exceptions.append(dict(type=row['type'],order=n,all_pairs=total,budget=7*n))
    threshold=[]
    for t in range(100):
        for n,poly in [(62+2*t,4*t*t+124*t),(63+2*t,4*t*t+128*t+68)]:
            diff=lower(n)-7*n
            assert 8*diff==poly and diff>=0
            threshold.append([n,str(diff)])
    assert lower(62)==7*62 and lower(61)<7*61
    paths=[f'scripts/{TAG}.py',f'scripts/{PREV}.py',f'results/{PREV}.json',
        f'results/{PREV}.md',f'results/{PREV}_review.md',f'scripts/{CENSUS}.py',
        f'results/{CENSUS}.json','results/fresh_multi_ideal_triangle_scope_20261017.md',
        'results/fresh_nikulin_narrow_20261015_compact_facet.md']
    return dict(status='PASS_EXACT_ALL_SIMPLE_AFFINE_EXCLUSION_DIMENSION61',
        sources={p:sha(p) for p in paths},affine_exclusion_hypothesis='Every vertex is simple; n>=10.',
        global_affine_order='Every connected affine subdiagram has exactly n nodes.',
        minimal_indefinite_core_order_bound=5,simplex_face_dimension_bound=4,
        high_rank_quasilanner_affine_checks=affine_rows,small_graph_radius_checks=small,
        inherited_classical_distance7_checks=graph_rows,exceptional_total_pair_checks=exceptions,
        triangle_path_bound=3,quadrilateral_path_bound=7,vertex_cost_per_dimension=7,
        threshold_polynomials_times8=dict(even=[4,124,0],odd=[4,128,68]),
        threshold_checks=threshold,simple_cusp_dimension_bound=61,
        scope='All finite-volume hyperbolic Coxeter polytopes whose ideal vertices are simple; finite vertices are automatically simple. No arithmetic hypothesis. No compact-facet127 change, unrestricted163 change, or novelty/current-best claim. Graph rank checks confirm inherited all-rank formulas, not a rank-cutoff proof.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],sum(r['marked_pairs'] for r in result['small_graph_radius_checks']),
          'small marked graphs;',len(result['inherited_classical_distance7_checks']),
          'classical graph costs; bound61')
