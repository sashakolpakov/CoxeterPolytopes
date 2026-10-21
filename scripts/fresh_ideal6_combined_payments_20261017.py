#!/usr/bin/env python3
"""Bind the finite9/ideal6 payment theorem and its piecewise vertex budgets."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal6_combined_payments_20261017'
TRI='fresh_ideal6_triangle_payments_20261017'
MIXED='fresh_ideal_distance6_mixed_replacement_20261017'
D14='fresh_affine_d14_long_pair_transfer_20261017'
FINITE='fresh_broad_finite_half_coupons_20261017'
EMBED='fresh_ideal6_small_factor_frontier_20261017'
OLD='fresh_order16_finite_cross_payments_20261017'
ENVELOPE={'1_to_13':'r/2','14':'20/3','15':'111/16','16':'127/17','17_and_larger':'95/13'}


def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def run():
    stems=(TRI,MIXED,D14,FINITE,EMBED,OLD,'fresh_ideal7_finite_cross_payments_20261017',
           'fresh_cycle16_half_coupon_transfer_20261017','fresh_selective_distance9_dimension108_20261017',
           'fresh_distance9_low_core_20261017')
    inputs={stem:read(stem)for stem in stems}
    for stem,d in inputs.items():
        assert d['status'].startswith('PASS')
        for rel,h in d['sources'].items():assert sha256((ROOT/rel).read_bytes()).hexdigest()==h,(stem,rel)
    assert len(inputs[TRI]['rows'])==920
    assert len(inputs[EMBED]['patterns'])==437
    assert inputs[D14]['status']=='PASS_UNIVERSAL_AFFINE_D14_LONG_PAIR_HALF_COUPON_TRANSFER'
    assert inputs[D14]['counts']==dict(marked_pairs=9,end_formula_checks=96,
        separated_cross_checks=256,center_affine13_payments=2,scalar_cost_rows=14,
        strict_nonextreme_end_pairs=4,equality_cases=72,equality_finite_corners=216,
        two_half_coupon_cases=32)
    d14_pairs={tuple(sorted(r['marks']))for r in inputs[D14]['equality']['rows']}
    assert d14_pairs=={tuple(sorted((a,b)))for a in(0,11,12)for b in(10,13,14)}
    assert all(r['finite_units']+F(r['coupon_total'])>=1
               for r in inputs[D14]['equality']['rows'])
    assert inputs[FINITE]['dimension_maximum']==108
    assert all(F(p[2])<=0 for r in inputs[FINITE]['rows']for p in r['endpoint_checks'])
    rows=[];costs={}
    for row in inputs[EMBED]['records']:
        m,kind=row['order'],row['kind'];paid=set(map(tuple,row['paid_pairs']))
        removed=set(map(tuple,row['proposed_unproved_order15_drop']))if m==15 and kind=='D'else set()
        if removed:assert removed==d14_pairs and removed<=paid
        final=paid-removed
        rows.append(dict(order=m,kind=kind,paid_pairs=[list(p)for p in sorted(final)],
                         newly_removed_pairs=[list(p)for p in sorted(removed)],cost=len(final)))
        costs.setdefault(str(m),{})[kind]=len(final)
    assert costs=={'15':dict(A=90,B=100,C=95,D=96),
                  '16':dict(A=96,B=107,C=103,D=111),
                  '17':dict(A=102,B=111,C=108,D=127)}
    maxima={m:str(F(max(v.values()),int(m)))for m,v in costs.items()}
    assert maxima=={'15':'20/3','16':'111/16','17':'127/17'}
    large=[]
    for r in inputs[MIXED]['basic_affine_cost']['all_rank_branches']:
        slope=F(95,13)-r['slope'];at=r['upper_order']if slope<0 else r['lower_order']
        assert at is not None and( r['upper_order'] is not None or slope>=0)
        slack=slope*at-r['intercept'];assert slack>=0
        large.append(dict(**r,common_bound='95/13',common_minimum_at=at,common_minimum_slack=str(slack)))
    assert F(95,13)>6 and len(large)==12
    assert max(r['selected_corner']['component_order_bound']for r in inputs['fresh_distance9_low_core_20261017']['extensions']['cases'])<=12
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
                             for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_COMPLETE_FINITE9_IDEAL6_LOCAL_PAYMENT_AND_BUDGET',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
                small_order_pair_sets=rows,small_order_costs=costs,small_order_normalized_maxima=maxima,
                large_order_linear_branches=large,ideal_rank_envelope=ENVELOPE,
                finite_vertex_budget='9n',finite_budget_dimension_maximum=108,
                triangle_minimum='2',quadrilateral_minimum='1',
                scope='Complete composition of the stated triangle and quadrilateral transfers, full preferred-pattern classification, broad whole-component finite coupons and all-rank ideal budgets. Actual cusp multiplicities are inherited. Dimension comparison is separate; external/formal review is not asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['small_order_costs'],flush=True)
