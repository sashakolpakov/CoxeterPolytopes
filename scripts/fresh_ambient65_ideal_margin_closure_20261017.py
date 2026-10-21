#!/usr/bin/env python3
"""Complete ambient65 cusp margins with proved order15/17 local rules."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ambient65_ideal_margin_closure_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
OLD='fresh_order13_b73_scalar_frontier_20261017'
D13='fresh_order13_affine_b_fork_projection_transfer_20261017'
B15='fresh_order15_affine_b_residual_terminal_transfer_20261017'
C15='fresh_order15_affine_c_terminal_projection_transfer_20261017'
D15='fresh_order15_affine_d_half_coupon_transfer_20261017'
D17='fresh_order17_affine_d_three_quarter_transfer_20261017'
FINITE='fresh_e6_d11_quarter_coupon_forest_budget_20261017'
SHAPE='fresh_cusp_link_shape_curvature_20261017'
FIBER='fresh_cusp_factor_payment_framework_20261017'
STEMS=(BASE,OLD,D13,B15,C15,D15,D17,FINITE,SHAPE,FIBER)


def profiles(left,minimum=1,prefix=(),low=False,high=False):
    if not left:yield prefix;return
    for r in range(minimum,left+1):
        if (r<10 and low) or (r>=18 and high):continue
        yield from profiles(left-r,r,prefix+(r,),low or r<10,high or r>=18)


def cost(r):
    if r<=10:return F(r,2)
    return {11:F(16,3),12:F(74,13),13:F(6),14:F(94,15),15:F(107,16),16:F(7)}[r] if r<=16 else F(163,23)


def scalar(rs):
    h=[1]
    for r in rs:
        new=[0]*(len(h)+r)
        for i,x in enumerate(h):
            for j in range(r+1):new[i+j]+=x
        h=new
    d=sum(rs);V=prod(r+1 for r in rs);G=sum(x*min(i,d-i) for i,x in enumerate(h))
    assert d==64 and sum(h)==V and h==h[::-1]
    charge=sum(map(cost,rs),F());curvature=F(5*G-64*V,2)
    margin=curvature/V-charge;assert margin>0
    return dict(ranks=rs,h_polynomial=h,link_vertices=V,G=G,
        normalized_cost=str(charge),curvature=str(curvature),normalized_margin=str(margin),
        absolute_margin=str(V*margin))


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert data[D13]['new_order13_costs']==dict(A=65,B=73,C=73,D=74)
    assert data[D15]['simultaneous_order15_costs']==dict(A=75,B=94,C=91,D=94)
    assert data[D17]['new_order17_costs']==dict(A=85,B=111,C=108,D=119)
    assert data[D17]['added_finite_coupon']['added_total_per_component_pair']=='1'
    assert data[FINITE]['finite_max65']=='487' and data[FINITE]['curvature65']=='480'
    rows=[scalar(rs) for rs in profiles(64)];assert len(rows)==1650
    prior={tuple(r['ranks']):r for r in data[OLD]['rows']}
    conditional={tuple(r['ranks']):r for r in data[OLD]['conditional_scenarios'][2]['rows']}
    assert set(prior)==set(conditional)=={tuple(r['ranks']) for r in rows}
    for r in rows:
        rs=r['ranks'];old=prior[rs]
        for key in ('h_polynomial','link_vertices','G'):assert r[key]==old[key]
        assert F(r['normalized_margin'])==F(old['margin'])+F(2*rs.count(14),15)+F(rs.count(16),17)
        assert r['normalized_margin']==conditional[rs]['margin']
    minimum=min(rows,key=lambda r:F(r['normalized_margin']))
    absolute=min(rows,key=lambda r:F(r['absolute_margin']))
    assert minimum['ranks']==(9,10,10,10,10,15) and minimum['normalized_margin']=='61/1936'
    assert absolute['ranks']==(64,) and absolute['absolute_margin']=='445/23'
    # Every legal concentration lowers product(r+1), as follows from
    # (a+1)(b+1)-(a-delta+1)(b+delta+1)=delta*(b-a+delta)>=0.
    extremum_gain=sum(2*j-64 for j in range(33,65))
    assert F(487)-480==7 and extremum_gain==1056
    assert 5*extremum_gain==5280==5*(65*65//4)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))]
    paths.append(ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md')
    return dict(status='PASS_COMPLETE_AMBIENT65_POSITIVE_IDEAL_MARGINS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        ambient_dimension=65,complete_terminal_profiles=1650,rows=rows,
        current_order15_costs=data[D15]['simultaneous_order15_costs'],
        current_order17_costs=data[D17]['new_order17_costs'],
        minimum_normalized_margin=minimum,minimum_absolute_margin=absolute,
        finite_maximum='487',finite_curvature='480',finite_deficit='7',
        inherited_global_constant='320',strengthened_global_additive_constant='5280',
        extremum_excess_lower_bound='1056',
        necessary_vertex_inequality='7*v_f >= 5280 + (445/23)*v_i',
        necessary_finite_vertex_count=755,necessary_finite_vertex_count_if_any_ideal=758,
        all_ideal_dimension_ceiling=64,unrestricted_dimension_ceiling=65,
        concentration='Arbitrary10..17 multiplicities; at most one1..9 and one>=18. Only linear-low and constant-high moves; each decreases curvature/V and V, preserving cost. Positive final margins imply absolute-margin monotonicity too.',
        scope='All cusp ranks covered by the inherited all-rank convex-order reduction and complete1650 terminal domain. Proved compatible order15/17 transfers replace the earlier hypothetical targets; new E6/D11 cost is counted exactly. All-ideal ceiling64 follows by excluding65 and using the completed unrestricted65 higher tail. General finite-volume ceiling remains65 because the finite deficit7 persists. No new example or arithmetic bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['complete_terminal_profiles'],out['minimum_normalized_margin']['normalized_margin'],out['minimum_absolute_margin']['absolute_margin'],flush=True)
