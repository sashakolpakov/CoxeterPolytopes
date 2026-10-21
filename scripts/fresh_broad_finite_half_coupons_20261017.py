#!/usr/bin/env python3
"""An exact finite-vertex budget for broad small-component half-coupons."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_broad_finite_half_coupons_20261017'


def run():
    small=[]
    for r in range(1,5):
        base=F(r*(r-1),2)
        slack=9*r-base
        assert slack>=F(15*r,2)
        small.append(dict(rank=r,base_cost=str(base),slack=str(slack),
                          slack_per_rank=str(slack/r)))
    large=[]
    for name,r,marks in [('D'+str(r),r,3)for r in range(12,17)]+[('A13',13,1),('A15',15,1)]:
        base=F(r*(r-1),2);slack=9*r-base
        assert r>=12 and slack>=24 and marks<=3
        large.append(dict(type=name,rank=r,paid_mark_count=marks,
                          base_cost=str(base),slack=str(slack)))
    cases=[]
    for t in range(10):
        coeff=F(3*t-15,2)
        a=108-12*t if coeff>0 else 0
        excess=coeff*a-24*t
        assert excess<=0
        cases.append(dict(large_count=t,maximum_small_rank=108-12*t,
                          coefficient=str(coeff),maximizing_small_rank=a,
                          maximum_excess=str(excess)))
    rows=[]
    for n in range(1,109):
        endpoints=[]
        for t in range(n//12+1):
            for a in sorted({0,n-12*t}):
                excess=F(3*a*t,2)-F(15*a,2)-24*t
                assert excess<=0
                endpoints.append([a,t,str(excess)])
        rows.append(dict(dimension=n,endpoint_checks=endpoints))
    stems=['fresh_selective_distance9_dimension108_20261017',
           'fresh_order16_finite_cross_payments_20261017']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
                             for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_BROAD_FINITE_HALF_COUPONS_COST_AT_MOST_9N',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
                dimension_maximum=108,small_components=small,large_components=large,
                all_rank_reduction='a=total rank of components of rank<=4; t=number of D12..D16,A13,A15 components; a+12t<=n. Cost<=9n+a*(3t/2-15/2)-24t.',
                cases_at_maximum_dimension=cases,rows=rows,
                coupon='One half for a marked pair in different whole finite components, one of rank<=4 and the other D12..D16 with its mark a fork center or short-arm leaf, or A13/A15 with its mark the unique middle.',
                base_rule='Unit weight for same-component distance<=9 or whole-component order<=17. Retain the inherited all-rank finite cost<=9rank.',
                scope='Universal finite-vertex counting bound only. Whether a given unpaid face admits one of these whole-component coupons requires a separate geometric proof; no dimension conclusion is asserted here.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],'108 dimensions,10 analytic count cases',flush=True)
