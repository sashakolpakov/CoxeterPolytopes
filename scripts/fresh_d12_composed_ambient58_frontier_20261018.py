#!/usr/bin/env python3
"""Bind the complete D12 face theorem to its exact typed forest count.

This discharges the count's geometric hypothesis. It supplies a necessary
vertex inequality in dimension58, not a contradiction or dimension57 bound.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_d12_half_face_composition_20261018 as G
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_d12_composed_ambient58_frontier_20261018'
BASE=C.L.TAG
TAIL='fresh_a15_d15_dimension58_20261018'
FINAL_REPORTS={
 G.TAG:'9f697fed9259db6a96ca94376ea69691794f9376995bdedbeb968a7851ba0a47',
 C.TAG:'4e6bf11e259af23023e43a10f7525ec160d00b0e3266d55ccd35295b7c30e7f8',
}


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    tags=(G.TAG,C.TAG,BASE,TAIL)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in tags}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h,s
    geo,count,prior,tail=(docs[s]for s in tags)
    assert G.B.TAG==C.L.TAG==BASE
    for folder,ext in(('scripts','py'),('results','json'),('results','md')):
        p=f'{folder}/{BASE}.{ext}'
        assert geo['sources'][p]==count['sources'][p]==digest(p)
    assert count['conditional_on_complete_D12_whole_face_composition']
    assert geo['triangle_minimum']=='2'and geo['quadrilateral_minimum']=='1'
    assert geo['ideal_weights_changed']is False and geo['contextual_coupons_changed']is False
    assert geo['finite_pair_rule']==dict(type='whole D12',pairs=[[0,10],[0,11]],old='1',new='1/2',old_cost='66',new_cost='65')
    assert count['ordinary_coupon_increments']=={'whole_A1_A11_endpoints':'1','whole_A2_A11_endpoints':'2'}
    assert [r['pair_increment']for r in geo['added_coupons']]==['1','2']
    types=tuple(C.TYPES)
    assert set(types)=={t for t in G.B.C.C.C.P.Z.B.B.types()if t[1]<=58}
    for a in types:
        assert G.actual_base(a)==C.actual_base(a)
        for b in types:assert G.actual_cross(a,b)==C.actual_cross(a,b)
    # Both wrappers inherit this exact same contextual rule, with no delta.
    assert G.B is C.L
    ff=count['finite'];ref=count['actual_typed_refinement'];ideal=count['ideal']
    assert ff['exact_actual_maximum']==ff['certified_upper']==ref['exact_actual_maximum']=='1527/4'
    assert ff['actual_maximum_claimed']and ff['actual_attainment_of_global_upper_claimed']
    classes=ref['classes'];assert len(classes)==ref['complete_expanded_upper_classes']==26
    retained=[r for r in classes if F(r['actual_maximum'])>377]
    assert len(retained)==ref['classes_still_above']==23
    assert ref['all_other_classes_at_most_curvature']and ref['both_context_guards_exact']
    for row in classes:
        assert max(r['thirtyseconds_actual_maximum']for r in row['rows'])==32*F(row['actual_maximum'])
        assert row['actual_obstruction_retained']==(F(row['actual_maximum'])>377)
        for forest in row['maximizers']:
            fs=list(map(tuple,forest))
            assert sum(n for k,n in fs)==58 and C.value(fs)==32*F(row['actual_maximum'])
    witness=[('A',16)]+[('D',14)]*3
    assert ref['actual_attaining_forest']==[list(t)for t in witness]
    assert F(C.value(witness),32)==F(1527,4)
    assert ideal['complete_profiles']==2552 and ideal['all_positive']
    assert ideal['weights_and_envelopes_unchanged']and ideal['inherited_positive_low8_concentration_and_factor_fiber']
    assert ideal['absolute_surplus']==prior['ideal']['absolute_surplus']=='36'
    assert ideal['normalized_minimum']=='87/10240'
    d=57;g0=sum(min(j,d-j)for j in range(d+1))
    curvature=F(5*g0-d*(d+1),2)
    strict=5*sum(2*j-d for j in range(d+1)if 2*j>d)
    excess=F(ff['exact_actual_maximum'])-curvature
    assert (g0,curvature,strict,excess)==(812,F(377),4205,F(19,4))
    assert count['finite_curvature']==str(curvature)and count['strict_global_constant']==strict
    lower=F(strict,1)/excess;ceil_lower=-(-lower.numerator//lower.denominator)
    assert ceil_lower==886 and tail['dimension_bound']==58
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in tags
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return dict(status='PASS_COMPOSED_D12_AMBIENT58_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},
        hypotheses=tail['hypotheses'],ambient_dimension=58,completed_dimension_bound=58,new_dimension_bound=None,
        conditional_count_geometric_hypothesis_discharged=True,
        triangle_minimum=2,quadrilateral_minimum=1,
        finite_base_identities=len(types),finite_cross_identities=len(types)**2,
        contextual_rule_identical=True,ideal_rule_unchanged=True,
        finite_exact_maximum=str(F(1527,4)),finite_curvature=str(curvature),finite_excess=str(excess),
        finite_attaining_formal_forest=[list(t)for t in witness],
        above_curvature_upper_classes=26,above_curvature_actual_classes=23,
        every_obstructing_forest_enumerated=False,geometric_realization_asserted=False,
        positive_ideal_profiles=2552,actual_cusp_surplus='36',strict_global_constant=str(strict),
        necessary_ambient58_inequality='(19/4)*v_f >= 4205 + 36*v_i',
        integer_vertex_requirement='v_f >= ceil((16820 + 144*v_i)/19)',minimum_finite_vertices=ceil_lower,
        contradiction_obtained=False,DP_replayed=False,geometric_census_repeated=False,
        scope='Complete local face composition plus frozen exact typed forest count and unchanged positive-cusp theorem. Fresh identical-rule and saved witness checks discharge the conditional count. This is a necessary vertex inequality, not a new dimension bound or polytope.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_exact_maximum'],out['above_curvature_actual_classes'],out['minimum_finite_vertices'],flush=True)
