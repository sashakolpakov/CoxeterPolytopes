#!/usr/bin/env python3
"""Finite8 plus marked-leaf9 and broad half-coupon budget, independent of face coverage."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from collections import deque

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite8_leaf9_half_coupon_budget_20261017'

def graph(kind,n):
    g=[set()for _ in range(n)]
    edges=[(i,i+1)for i in range(n-1)]if kind=='path'else[(i,i+1)for i in range(n-2)]+[(n-3,n-1)]
    for i,j in edges:g[i].add(j);g[j].add(i)
    return g

def graph_counts(kind,n):
    g=graph(kind,n);short=leaf=0
    for i in range(n):
        d={i:0};q=deque([i])
        while q:
            k=q.popleft()
            for j in g[k]-d.keys():d[j]=d[k]+1;q.append(j)
        for j in range(i+1,n):
            short+=d[j]<=8
            leaf+=d[j]==9 and(len(g[i])==1 or len(g[j])==1)
    return short,leaf

def run():
    # These samples check the all-rank graph identities proved in the report;
    # they do not supply the finite-type classification or extrapolate it.
    rows=[]
    for n in range(18,81):
        for kind,constant,leaves in(('path',36,2),('D',35,3)):
            short,leaf=graph_counts(kind,n)
            assert short==8*n-constant and leaf==leaves
            assert short+leaf<=8*n
            rows.append(dict(type=kind,rank=n,short_pairs=short,leaf9_pairs=leaf,
                base_cost=short+leaf,slack=8*n-short-leaf))
    small=[]
    for n in range(1,18):
        base=F(n*(n-1),2);slack=8*n-base;assert slack>=0
        if n<=4:assert slack>=F(13*n,2)
        small.append(dict(rank=n,all_pair_cost=str(base),slack=str(slack)))
    large=[]
    for name,n,marks in[('D'+str(n),n,3)for n in range(12,17)]+[('A13',13,1),('A15',15,1)]:
        slack=F(8*n)-F(n*(n-1),2);assert slack>=8 and marks<=3
        large.append(dict(type=name,rank=n,eligible_marks=marks,slack=str(slack)))
    cases=[];endpoints=[]
    for t in range(7):
        c=F(3*t-13,2);a=80-12*t if c>0 else 0
        val=c*a-8*t;assert val<=0
        cases.append(dict(large_count=t,maximum_small_rank=80-12*t,
            maximizing_small_rank=a,maximum_excess=str(val)))
    for n in range(1,81):
        for t in range(n//12+1):
            for a in sorted({0,n-12*t}):
                val=F(a*(3*t-13),2)-8*t;assert val<=0
                endpoints.append(dict(dimension=n,small_rank=a,large_count=t,excess=str(val)))
    stems=['fresh_broad_finite_half_coupons_20261017','fresh_selective_distance9_dimension108_20261017']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE8_LEAF9_HALF_COUPON_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        dimension_maximum=80,finite_vertex_budget='8n',all_rank_base_formulas={
            'A_or_B_rank_at_least18':'8r-36+2=8r-34',
            'D_rank_at_least18':'8r-35+3=8r-32',
            'rank_at_most17':'binom(r,2)<=8r'},
        graph_identity_checks=rows,small_rank_slack=small,large_coupon_components=large,
        analytic_count_cases=cases,dimension_endpoints=endpoints,
        surplus_bound='a*(3t/2-13/2)-8t, with a+12t<=n<=80',
        rule='Unit on same-component distance<=8, or distance9 with at least one actual leaf mark, or whole component rank<=17. Add the frozen broad small<=4 x D12..D16/A13/A15 half-coupons.',
        scope='Universal finite-vertex counting theorem only. The new finite face-payment coverage and whole-leaf/whole-component geometric premises are separate obligations; no dimension improvement asserted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['dimension_endpoints']),'endpoint checks',flush=True)
