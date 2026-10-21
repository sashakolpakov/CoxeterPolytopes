#!/usr/bin/env python3
"""Remove the retained D13 center half by upgrading its whole finite coupons.

This is an exact composition of the frozen universal half theorem, not a
new census or a finite-forest budget. The written proof treats every branch.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order13_d_center_pair_zero_transfer_20261018'
HALF='fresh_order13_d_center_pair_half_transfer_20261018'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
BUDGET='fresh_sparse_coupon_rank63_forest_budget_20261018'
STEMS=(HALF,CURRENT,BUDGET)


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,doc in data.items():
        assert doc['status'].startswith('PASS'),tag
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old=data[HALF];assert data[CURRENT]['dimension_bound']==62
    assert old['results']['new_affine_total']=='147/2' and old['new_weight']=='1/2'
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert weights[0,8]==F(1,2) and sum(weights.values(),F())==F(147,2)
    weights[0,8]=F();assert sum(weights.values(),F())==73
    coupons=[]
    for row,count,old_total,new_total in zip(old['new_coupons'],(2,1),(F(1),F(1,2)),(F(2),F(1))):
        assert row['weight']=='1/2' and F(row['component_pair_increment'])==old_total
        assert count*F(1,2)==old_total and count==new_total
        coupons.append(dict(small=row['small'],small_role=row['small_role'],large=row['large'],
                            large_role=row['large_role'],old_weight='1/2',new_weight='1',
                            eligible_marked_pairs=count,old_component_pair_total=str(old_total),
                            new_component_pair_total=str(new_total),increment=str(new_total-old_total)))
    assert [(r['small'],r['small_role']) for r in coupons]==[('whole_A2','either_endpoint'),('whole_A3','unique_middle')]
    branches=[]
    for row in old['immediate_affine_unit_corners']:
        assert row['retained_payment']=='1' and 13 in row['marks']
        branches.append(dict(branch='immediate_affine',side=row['side'],profile=row['profile'],
                             type=row['whole_factor'],canonical_pair=row['canonical_pair'],distinct_payment='1'))
    for row in old['finite_coupon_corners']:
        assert row['original_ideal_payment']==row['coupon_payment']=='1/2' and row['total']=='1'
        assert 13 in row['marks'] and set(row['marks'])!={0,8}
        expected='middle' if row['whole_small_type']=='A3' else 'endpoint'
        assert row['small_mark_role']==expected and row['whole_large_type']=='D10' and row['large_mark_role']=='fork_center'
        branches.append(dict(branch='whole_finite_coupon',side=row['side'],old_leaf=row['attached_old_leaf'],
                             small=row['whole_small_type'],small_role=expected,large='D10',large_role='fork_center',
                             original_ideal_payment='0',distinct_payment='1',new_total='1'))
    for row in old['mixed_energy_one_cases']:
        assert row['opposite_black_path']==4 and row['payment']=='1'
        branches.append(dict(branch='mixed_energy_one',type=row['type'],path_length=4,distinct_payment='1'))
    assert len(branches)==17 and all(r['distinct_payment']=='1' for r in branches)
    families={r['kind']:F(r['total_cost']) for r in data[CURRENT]['new_factor_interfaces'] if r['order']==13}
    assert families==dict(A=F(65),B=F(73),C=F(73),D=F(147,2))
    families['D']=F(73);assert max(families.values())==73
    assert old['formulas']['remaining_nonextreme_maximum']==3
    assert old['formulas']['one_extreme_weighted_energy_maximum']==1
    assert old['formulas']['both_extreme_unscaled_energy']==1
    paths=[Path(__file__)]+[ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_UNIVERSAL_ORDER13_D_CENTER_PAIR_ZERO_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='D12',order=13,changed_pairs=[[0,8]],old_weight='1/2',new_weight='0',
                old_affine_total='147/2',new_affine_total='73',new_normalized_envelope='73/13',
                all_family_order13_costs={k:str(v) for k,v in families.items()},
                weighted_pairs=[dict(pair=list(p),weight=str(w)) for p,w in sorted(weights.items())],
                upgraded_coupons=coupons,additional_finite_charge='a2*d10 + (a3*d10)/2',
                inherited_branches_with_independent_units=branches,
                triangle_minimum=2,quadrilateral_minimum=1,
                actual_corner_priority='Apply this independent unit-alternative theorem to any actual changed D13 center corner. If none exists, all ideal weights equal the completed62 rule on that face and finite coupons only increase.',
                finite_budget_recomputed=False,new_dimension_bound=None,
                scope='Complete local zero-transfer composition using every branch of the frozen arbitrary-spectator half theorem. Coupon branches now pay a distinct full unit without the original ideal half. Finite optimization and any ambient theorem remain separate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'order13 maximum',out['new_affine_total'],'coupon delta',out['additional_finite_charge'],flush=True)
