#!/usr/bin/env python3
"""Exact radius improvement and classical finite/affine graph counts."""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_distance_bound_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
FACET='fresh_nikulin_narrow_20261015_compact_facet'


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def graph(n,edges):
    g=[set() for _ in range(n)]
    for a,b in edges:g[a].add(b);g[b].add(a)
    return g


def distances(g,starts):
    d={v:0 for v in starts};q=deque(starts)
    while q:
        v=q.popleft()
        for w in sorted(g[v]):
            if w not in d:d[w]=d[v]+1;q.append(w)
    assert len(d)==len(g)
    return d


def pair_counts(g):
    c=Counter()
    for a in range(len(g)):
        for b,d in distances(g,[a]).items():
            if a<b:c[d]+=1
    assert sum(c.values())==len(g)*(len(g)-1)//2
    return c


def path(m):return [(i,i+1) for i in range(m-1)]


def star_arms(lengths):
    edges=[];n=1
    for length in lengths:
        last=0
        for _ in range(length):edges.append((last,n));last=n;n+=1
    return graph(n,edges)


def lower(d):
    return F(d*(d-6),8) if d%2==0 else F((d-1)*(d-5),8)


def run():
    census=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    assert census['status']=='PASS_EXACT_HIGH_RANK_QUASILANNER_ATTACHMENT_CENSUS'
    for p,h in census['sources'].items():assert sha(p)==h,p
    assert census['orders']==[[6,12],[7,3],[8,4],[9,4],[10,3]]
    rows=[];hist=Counter();maximizers=[]
    for L in census['diagrams']:
        n=L['order'];g=graph(n,[(a,b) for a,b,_ in L['edges']])
        finite=[i for i,deletion in enumerate(L['deletions'])
                if all(c['kind']=='finite' for c in deletion)]
        assert finite==L['finite_deletions']
        for black in combinations(finite,2):
            d=distances(g,black);radius=max(d.values())
            assert radius<=7
            hist[radius]+=1
            row=dict(diagram=L['index'],order=n,black=list(black),
                     distances=[d[v] for v in range(n)],radius=radius,
                     farthest=[v for v in range(n) if d[v]==radius])
            rows.append(row)
            if radius==7:maximizers.append(row)
    assert len(rows)==384
    assert hist=={1:10,2:123,3:143,4:61,5:31,6:12,7:4}
    assert [(r['diagram'],r['black'],r['farthest']) for r in maximizers]==[
        (23,[0,1],[8]),(25,[1,2],[9]),(25,[1,3],[9]),(25,[2,3],[9])]
    checks=[]
    # Graph formula checks are independent breadth-first enumerations. The
    # formulas themselves, not this finite range, establish every rank.
    for n in range(2,65):
        g=graph(n,path(n));c=pair_counts(g)
        assert all(c[k]==max(n-k,0) for k in range(1,n+1))
        assert sum(v for k,v in c.items() if k<=15)<=15*n
        checks.append(dict(type='path',order=n,counts=sorted(c.items())))
    for n in range(3,65):
        c=pair_counts(graph(n,path(n)+[(n-1,0)]))
        assert all(c[k]==(n if 2*k<n else n//2 if 2*k==n else 0)
                   for k in range(1,n+1))
        assert all(v<=n for v in c.values())
        checks.append(dict(type='cycle',order=n,counts=sorted(c.items())))
    for m in range(2,63):
        n=m+2;c=pair_counts(graph(n,path(m)+[(0,m),(0,m+1)]))
        assert all(c[k]==max(m-k,0)+2*(k<=m)+(k==2) for k in range(1,n+1))
        assert all(v<=n for v in c.values())
        checks.append(dict(type='single_fork',order=n,counts=sorted(c.items())))
    for m in range(2,61):
        n=m+4;c=pair_counts(graph(n,path(m)+[(0,m),(0,m+1),(m-1,m+2),(m-1,m+3)]))
        assert all(c[k]==max(m-k,0)+4*(k<=m)+2*(k==2)+4*(k==m+1)
                   for k in range(1,n+1))
        assert all(v<=n for v in c.values())
        checks.append(dict(type='double_fork',order=n,counts=sorted(c.items())))
    exceptions=[]
    for name,lengths in [('D4_affine',[1,1,1,1]),('E6',[1,2,2]),('E7',[1,2,3]),
                         ('E8',[1,2,4]),('E6_affine',[2,2,2]),
                         ('E7_affine',[1,3,3]),('E8_affine',[1,2,5])]:
        g=star_arms(lengths);n=len(g);c=pair_counts(g)
        assert n<=9 and sum(c.values())<=4*n<=15*n
        exceptions.append(dict(type=name,order=n,counts=sorted(c.items())))
    # Strict h-vector lower bounds contradict an upper budget even when
    # their scalar difference is zero (the first even whole-polytope case).
    whole=[];facet=[]
    for t in range(100):
        for n,expected in [(126+2*t,F(4*t*t+252*t,8)),
                           (127+2*t,F(4*t*t+256*t+132,8))]:
            diff=lower(n)-15*n
            assert diff==expected and diff>=0
            whole.append([n,str(diff)])
        for d,expected in [(127+2*t,F(4*t*t+256*t+12,8)),
                           (128+2*t,F(4*t*t+260*t+136,8))]:
            diff=lower(d)-15*(d+1)
            assert diff==expected and diff>0
            facet.append([d,str(diff)])
    assert lower(142)==17*142 and lower(143)>17*143
    paths=[f'scripts/{TAG}.py',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           f'results/{CENSUS}.md',f'scripts/{FACET}.py',f'results/{FACET}.json',
           f'results/{FACET}.md','results/fresh_multi_ideal_triangle_scope_20261017.md',
           'results/fresh_finite_volume_simplex_truncation_20261017.md']
    return dict(status='PASS_EXACT_DISTANCE15_SIMPLE_CUSP_AND_COMPACT_FACET_BOUNDS',
        sources={p:sha(p) for p in paths},black_pair_checks=rows,
        radius_histogram=sorted(hist.items()),radius_maximizers=maximizers,
        low_order_argument='Any connected diagram of order<=5 and two marked nodes has distance to their set<=3; a dotted two-node core has radius0.',
        triangle_path_bound=8,quadrilateral_path_bound=15,
        graph_formula_checks=checks,exceptional_graph_checks=exceptions,
        all_rank_graph_formulas=dict(path='max(N-k,0)',
            cycle='N if 2k<N; N/2 if 2k=N; 0 otherwise',
            single_fork='max(m-k,0)+2*[k<=m]+[k=2], N=m+2,m>=2',
            double_fork='max(m-k,0)+4*[k<=m]+2*[k=2]+4*[k=m+1], N=m+4,m>=2',
            exceptions='Every remaining finite/affine connected underlying graph has N<=9, so all pairs<=4N<=15N.'),
        whole_simple_threshold_checks=whole,compact_facet_threshold_checks=facet,
        threshold_polynomials=dict(whole_even=[4,252,0],whole_odd=[4,256,132],
            facet_odd=[4,256,12],facet_even=[4,260,136]),
        simple_cusp_dimension_bound=125,compact_facet_ambient_dimension_bound=127,
        coarse_distance17_simple_cusp_bound=141,
        scope='All-simple means every finite and ideal vertex is simple; finite-volume Coxeter polytopes automatically have simple finite vertices. The independent compact-facet conclusion permits arbitrary nonsimple ideal vertices elsewhere. No claim that every chamber lies in either class. The full unrestricted163 bound is unchanged. No novelty or current-best claim; no re-enumeration of the frozen26-type catalogue.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['black_pair_checks']),'black pairs;',
          len(result['graph_formula_checks']),'classical graph checks; bounds125/127')
