#!/usr/bin/env python3
"""Bind completed local, finite and ideal proofs into a dimension65 restriction.

This is an exact composition certificate, not another geometry or partition run.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_dimension65_large_finite_component_obstruction_20261017'
LOCAL='fresh_finite_d16_d17_long_pair_half_transfer_20261017'
BUDGET='fresh_d16_d17_d14_coupon_forest_budget_20261017'
IDEAL='fresh_ambient65_ideal_margin_closure_20261017'
REVIEW=IDEAL+'_review'
TAIL='fresh_ideal5_weighted_dimension65_20261017'
STEMS=(LOCAL,BUDGET,IDEAL,REVIEW,TAIL)
FROZEN_BUDGET={
    'scripts/'+BUDGET+'.py':'b69a37086bebb081fc007d5e994f91298cf13ded2c6ced98b704a9e599c40739',
    'results/'+BUDGET+'.json':'c3ac67cd54274fafcf62ccb2052dbff4eb0617c3f1b48a785f0e17e4c2435a69',
    'results/'+BUDGET+'.md':'2c109ca38bb38e69a9b28d971040e2029358a671bb1c961da35138ee321ad45d',
}


def digest(p):return sha256(p.read_bytes()).hexdigest()
def ceil(q):return -(-q.numerator//q.denominator)


def run():
    for p,h in FROZEN_BUDGET.items():assert digest(ROOT/p)==h,p
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    checked=0
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(ROOT/p)==h,(s,p);checked+=1
    local,budget,ideal,review,tail=(data[s] for s in STEMS)
    costs=local['new_finite_component_costs']
    assert costs==dict(A16=116,B16=104,D16=118,D17=125)
    assert all(budget['conditional_costs'][k]==v for k,v in costs.items())
    assert budget['conditional_costs']['A17']==114 and budget['conditional_costs']['B17']==105
    assert local['added_coupon']==dict(weight='1/2',small_components=['whole A1','whole A2'],
        small_mark='any node',large_component='whole D14',
        large_mark='the unique long-arm neighbor of its fork center',canonical_D14_node=10,
        added_A1_D14_cost='1/2',added_A2_D14_cost='1')
    assert budget['new_coupon']==('Weight1/2 between every mark of whole A1 or A2 and the unique '
        'long-arm neighbor of the fork in separate whole D14; added totals1/2 and1 per component pair respectively.')
    domains={r['small_components_only']:r for r in budget['results']}
    assert set(domains)=={False,True}
    large_row=domains[False]['dimension_maxima'][65]
    small_row=domains[True]['dimension_maxima'][65]
    assert large_row==dict(dimension=65,doubled_cost=974,forest=['D65'],curvature='480',slack='-7')
    assert small_row==dict(dimension=65,doubled_cost=958,forest=['D16','D16','D16','D17'],curvature='480',slack='1')
    beta=F((65-1)*(65-5),8);large=F(large_row['doubled_cost'],2);small=F(small_row['doubled_cost'],2)
    assert beta==480 and large==487 and small==479
    assert budget['max65_all']=='487' and budget['max65_components_at_most17']=='479'
    assert ideal['complete_terminal_profiles']==review['profiles']==1650
    assert ideal['minimum_normalized_margin']['normalized_margin']=='61/1936'
    assert ideal['minimum_absolute_margin']['absolute_margin']=='445/23'
    assert ideal['minimum_absolute_margin']['ranks']==[64]
    assert review['minimum_normalized']==ideal['minimum_normalized_margin']
    assert review['minimum_absolute']==ideal['minimum_absolute_margin']
    # Rebind the independently checked complete scalar rows, without
    # re-enumerating partitions or recalculating their curvature/payment data.
    canonical=json.dumps(sorted(ideal['rows'],key=lambda r:r['ranks']),sort_keys=True,separators=(',',':')).encode()
    assert sha256(canonical).hexdigest()==review['producer_scalar_rows_sha256']==review['independent_rows_sha256']
    surplus=F(ideal['minimum_absolute_margin']['absolute_margin'])
    constant=F(ideal['strengthened_global_additive_constant'])
    assert constant==5*(65*65//4)==5280 and surplus==F(445,23)
    excess=dict(v_L=beta-large,v_S=beta-small,v_i=surplus,constant=constant)
    assert excess==dict(v_L=F(-7),v_S=F(1),v_i=F(445,23),constant=F(5280))
    integer={k:int(23*v) for k,v in excess.items()}
    assert all(F(integer[k])==23*v for k,v in excess.items())
    assert integer==dict(v_L=-161,v_S=23,v_i=445,constant=121440)
    thresholds=[]
    for vi,want in ((0,755),(1,758)):
        lower=(constant+surplus*vi)/(large-beta);answer=ceil(lower)
        assert answer==want and F(answer-1)<lower<=answer
        thresholds.append(dict(ideal_vertex_lower_bound=vi,small_finite_vertex_lower_bound=0,
            rational_large_finite_lower_bound=str(lower),integer_large_finite_lower_bound=answer))
    assert tail['dimension_bound']==ideal['unrestricted_dimension_ceiling']==65
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_DIMENSION65_LARGE_FINITE_COMPONENT_OBSTRUCTION',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},checked_input_bindings=checked,
        ambient_dimension=65,large_finite_vertex_definition='A finite vertex with at least one irreducible finite Coxeter component of rank>=18.',
        small_finite_vertex_definition='A finite vertex all of whose irreducible finite Coxeter components have rank<=17.',
        local_component_costs=costs,new_coupon=local['added_coupon'],
        finite_curvature=str(beta),arbitrary_finite_budget=str(large),small_component_finite_budget=str(small),
        exact_finite_budget_rows=dict(arbitrary=large_row,small_components=small_row),
        inherited_complete_ideal_profiles=1650,inherited_normalized_ideal_margin='61/1936',
        inherited_absolute_ideal_margin=str(surplus),scalar_rows_sha256=sha256(canonical).hexdigest(),
        strengthened_global_constant=str(constant),lower_minus_upper_coefficients={k:str(v) for k,v in excess.items()},
        necessary_inequality='7*v_L >= v_S + 5280 + (445/23)*v_i',
        equivalent_integer_inequality='161*v_L >= 23*v_S + 121440 + 445*v_i',
        finite_count_thresholds=thresholds,minimum_large_finite_vertices=755,
        minimum_large_finite_vertices_if_any_ideal=758,
        dimension_ceiling_if_every_finite_vertex_has_components_at_most17=64,
        unrestricted_dimension_ceiling=65,
        scope='Exact composition only. The frozen local theorem supplies geometric coverage for the finite budget; current ideal payments are unchanged. Scalar rows are hash-bound to their independent reconstruction, not rerun. Complete geometric censuses, finite Bellman verification, all-rank concentration and the strengthened extremum argument are inherited proved inputs. No new unrestricted bound or chamber construction.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['necessary_inequality'],flush=True)
