#!/usr/bin/env python3
"""Finite8 cost with an additional whole A4/D13-neighbor half coupon."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'
BASE='fresh_finite8_a17_d11_coupon_budget_20261017'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run():
    prior=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert prior['finite_vertex_budget']=='8n' and prior['dimension_maximum']==80
    for path,h in prior['sources'].items():
        assert digest(ROOT/path)==h,path
    old_endpoints=0
    for row in prior['dimension_endpoints']:
        n,a,t,u=(row[k]for k in('dimension','small_rank','old_large_count','A17_count'))
        assert a+11*t+17*u<=n<=80
        assert F(row['excess'])==a*(F(3*t,2)+u-6)-8*t<=0
        old_endpoints+=1
    assert old_endpoints==1487
    # The new D13 mark is the unique non-leaf neighbor of the fork center.
    # It differs from all three old D13 marks. Whole A4 has two endpoints;
    # two marked pairs of weight1/2 give one unit per component pair.
    small_rank=4
    new_cost_per_component_pair=F(2,2)
    assert new_cost_per_component_pair==1
    d13_rank=13;d13_unit_cost=F(d13_rank*(d13_rank-1),2)
    d13_slack=8*d13_rank-d13_unit_cost
    assert d13_unit_cost==78 and d13_slack==26
    extra_slack=d13_slack-8
    assert extra_slack==18 and d13_rank-11==2
    # Analytic dominance: b>0 implies t>=1 and a<=80-11-2=67.
    maximum_small_rank_with_D13=80-13
    per_D13_delta_upper=F(maximum_small_rank_with_D13,4)-extra_slack
    assert maximum_small_rank_with_D13==67 and per_D13_delta_upper==F(-5,4)<0
    endpoints=[];positive80=[];counts={};total=0;at80=0;raw80=0
    for n in range(1,81):
        cases=0
        for t in range(n//11+1):
            for b in range(t+1):
                for u in range((n-11*t-2*b)//17+1):
                    maximum_a=n-11*t-2*b-17*u
                    assert maximum_a>=0
                    coefficient=F(3*t,2)+u-6+F(b,4)
                    raw80+=2*(n==80)
                    for a in sorted({0,maximum_a}):
                        old=a*(F(3*t,2)+u-6)-8*t
                        delta=b*(F(a,4)-extra_slack)
                        excess=a*coefficient-8*t-extra_slack*b
                        assert old<=0 and excess==old+delta<=0
                        assert 0<=b<=t and a+11*t+2*b+17*u<=n<=80
                        if b:
                            assert a<=67 and delta<=F(-5*b,4)<0
                        else:
                            assert delta==0
                        endpoints.append(dict(dimension=n,small_rank=a,old_large_count=t,
                            D13_count=b,A17_count=u,old_excess=str(old),
                            additional_excess=str(delta),excess=str(excess)))
                        cases+=1;total+=1;at80+=(n==80)
                    if n==80 and coefficient>0:
                        positive80.append(dict(old_large_count=t,D13_count=b,A17_count=u,
                            maximum_small_rank=maximum_a,coefficient=str(coefficient),
                            maximum_excess=str(maximum_a*coefficient-8*t-extra_slack*b)))
        counts[str(n)]=cases
    assert total==len(endpoints)==3465 and at80==138 and raw80==140
    assert len(positive80)==38 and all(F(r['maximum_excess'])<0 for r in positive80)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{BASE}.{ext}'for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE8_A4_D13_NEIGHBOR_HALF_COUPON_BUDGET',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_maximum=80,finite_vertex_budget='8n',
        retained_rule=prior['rule'],
        additional_rule='Add one half to an unordered pair in distinct whole components when one is an endpoint of whole A4 and the other is the unique non-leaf neighbor of the fork center of whole D13. These D13 marks differ from the inherited fork-center/short-leaf marks; count each pair once.',
        component_count=dict(A4_endpoints=2,D13_new_marks=1,half_weight='1/2',
            cost_per_A4_D13_component_pair='1',A4_count_bound='a/4',
            D13_unit_cost='78',D13_base_slack='26',additional_D13_slack='18',additional_D13_rank=2),
        variables='a is the total rank of components of rank<=5; t counts D11..D16/A13/A15; b counts whole D13 and satisfies0<=b<=t; u counts whole A17.',
        rank_constraint='a+11t+2b+17u<=n<=80',
        surplus_bound='a*(3t/2+u-6+b/4)-8t-18b',
        analytic_dominance=dict(old_surplus='a*(3t/2+u-6)-8t',difference='b*(a/4-18)',
            maximum_small_rank_if_b_positive=67,per_D13_difference_upper='-5/4',
            proof='If b=0 use the frozen budget. If b>=1 then t>=1, a<=67, and the new surplus is at most the old surplus minus5b/4. The old rank constraint remains valid.'),
        inherited_endpoint_checks=old_endpoints,dimension_endpoint_counts=counts,
        dimension_endpoints=endpoints,positive_coefficient_cases_at_80=positive80,
        counts=dict(distinct_endpoints_all_dimensions=3465,distinct_endpoints_at_80=138,
            endpoints_at_80_with_duplicate_zero_endpoints=140,positive_coefficient_cases_at_80=38),
        scope='Universal finite-vertex counting theorem for n<=80. Retains the entire old finite8 rule and adds the stated whole A4/D13-neighbor half coupon. Application to an actual deficient face requires a separate local transfer. No ideal-cost change or new ambient dimension bound is asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
