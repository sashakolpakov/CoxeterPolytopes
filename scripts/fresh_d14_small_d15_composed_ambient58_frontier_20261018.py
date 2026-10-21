#!/usr/bin/env python3
"""Bind the joint local rule to its complete exact ambient58 forest count.

This discharges the count's geometric hypothesis, yielding a necessary
vertex inequality. It does not exclude dimension58 or construct a polytope.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_d14_a16_d15_joint_face_composition_20261018 as G
import fresh_d14_small_d15_ambient58_frontier_20261018 as C

ROOT=G.ROOT
TAG='fresh_d14_small_d15_composed_ambient58_frontier_20261018'
BASE=C.C.TAG
IDEAL_BASE=C.C.L.TAG
TAIL='fresh_a15_d15_dimension58_20261018'
FINAL_REPORTS={
 G.TAG:'ee6c8bfd412f0502abf10f3f553673b8ec06b327de2cd0488e860dcd40c9363f',
 C.TAG:'36c5a69fc485405d4e794bd304cdd432c803464060d77ae57cee421c6f83a83e',
}


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    tags=(G.TAG,C.TAG,BASE,IDEAL_BASE,TAIL)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in tags}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h,s
    geo,count,prior,idealbase,tail=(docs[s]for s in tags)
    assert G.C is C.C and G.U is C.P and G.C.context is C.context
    for folder,ext in(('scripts','py'),('results','json'),('results','md')):
        p=f'{folder}/{BASE}.{ext}'
        assert geo['sources'][p]==count['sources'][p]==digest(p)
    assert count['conditional_on_joint_D14_A16_and_D15_center_whole_face_composition']
    assert geo['triangle_minimum']=='2'and geo['quadrilateral_minimum']=='1'
    assert geo['ideal_weights_changed'] is False and geo['contextual_rules_changed'] is False
    assert geo['D14_cost']=='179/2'and geo['A16_cost']=='443/4'and geo['small_D15_center_weight']=='3/8'
    assert count['actual_base_deltas']=={'D14':'-1','A16':'1/2'}
    for a in C.TYPES:
        assert G.actual_base(a)==C.actual_base(a)
        for b in C.TYPES:assert G.actual_cross(a,b)==C.actual_cross(a,b)
    ff=count['finite'];ref=count['actual_typed_refinement'];ideal=count['ideal']
    assert ff['exact_actual_maximum']==ff['certified_upper']==ref['exact_actual_maximum']=='763/2'
    assert ff['actual_maximum_claimed']and ff['actual_attainment_of_global_upper_claimed']
    classes=ref['classes'];assert len(classes)==ref['complete_expanded_upper_classes']==20
    assert len([r for r in classes if F(r['actual_maximum'])>377])==ref['classes_still_above']==16
    assert ref['all_other_classes_at_most_curvature']and ref['both_context_guards_exact']
    for row in classes:
        assert max(r['thirtyseconds_actual_maximum']for r in row['rows'])==32*F(row['actual_maximum'])
        assert row['actual_obstruction_retained']==(F(row['actual_maximum'])>377)
        for forest in row['maximizers']:
            fs=list(map(tuple,forest))
            assert sum(n for k,n in fs)==58 and C.value(fs)==32*F(row['actual_maximum'])
    witness=[('A',16)]*2+[('D',26)]
    assert ref['actual_attaining_forest']==[list(t)for t in witness]
    assert F(C.value(witness),32)==F(763,2)
    assert prior['finite']['exact_actual_maximum']=='1527/4'
    assert F(C.C.value(witness),32)==F(761,2)
    assert ideal['complete_profiles']==2552 and ideal['all_positive']
    assert ideal['weights_and_envelopes_unchanged']and ideal['inherited_positive_low8_concentration_and_factor_fiber']
    assert ideal['absolute_surplus']==idealbase['ideal']['absolute_surplus']=='36'
    assert ideal['normalized_minimum']=='87/10240'
    d=57;g0=sum(min(j,d-j)for j in range(d+1))
    curvature=F(5*g0-d*(d+1),2)
    strict=5*sum(2*j-d for j in range(d+1)if 2*j>d)
    excess=F(ff['exact_actual_maximum'])-curvature
    assert (g0,curvature,strict,excess)==(812,F(377),4205,F(9,2))
    assert count['finite_curvature']==str(curvature)and count['strict_global_constant']==strict
    lower=F(strict)/excess;ceil_lower=-(-lower.numerator//lower.denominator)
    assert ceil_lower==935 and tail['dimension_bound']==58
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in tags for folder,ext in
        [('scripts','py'),('results','json'),('results','md')]}
    return dict(status='PASS_COMPOSED_D14_SMALL_D15_AMBIENT58_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},hypotheses=tail['hypotheses'],
        ambient_dimension=58,completed_dimension_bound=58,new_dimension_bound=None,
        conditional_count_geometric_hypothesis_discharged=True,
        triangle_minimum=2,quadrilateral_minimum=1,finite_base_identities=len(C.TYPES),
        finite_cross_identities=len(C.TYPES)**2,contextual_rule_identical=True,ideal_rule_unchanged=True,
        finite_exact_maximum='763/2',previous_maximum='1527/4',finite_curvature=str(curvature),finite_excess=str(excess),
        finite_attaining_formal_forest=[list(t)for t in witness],
        above_curvature_upper_classes=20,above_curvature_actual_classes=16,
        every_obstructing_forest_enumerated=False,geometric_realization_asserted=False,
        positive_ideal_profiles=2552,actual_cusp_surplus='36',strict_global_constant=str(strict),
        necessary_ambient58_inequality='(9/2)*v_f >= 4205 + 36*v_i',
        integer_vertex_requirement='v_f >= ceil((8410 + 72*v_i)/9)',minimum_finite_vertices=ceil_lower,
        contradiction_obtained=False,DP_replayed=False,geometric_census_repeated=False,
        scope='Complete joint local face composition plus frozen exact typed forest count and unchanged positive-cusp theorem. Identical-rule and saved witness checks discharge the conditional count. A necessary vertex inequality, not a new dimension bound or polytope.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_exact_maximum'],out['above_curvature_actual_classes'],out['minimum_finite_vertices'],flush=True)
