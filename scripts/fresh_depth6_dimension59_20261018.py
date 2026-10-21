#!/usr/bin/env python3
"""Compose the completed depth-six local theorem and exact finite404 count.

Read-only theorem/certificate binding and exact saved-row checks. No geometric
census, convolution, Bellman optimization or new construction is performed.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_depth6_dimension59_20261018'
LOCAL='fresh_depth6_finite7_adaptive_transfer_20261018'
COUNT='fresh_depth6_sparse_ambient60_forest_frontier_20261018'
BASE='fresh_small_d16_coupon_and_ideal17_rebalance_20261018'
TAIL='fresh_order11_order16_dimension60_20261018'
CLARIFICATION='results/fresh_depth6_sparse_ambient60_forest_frontier_20261018_domain_clarification.md'
FINAL_REPORTS={
    LOCAL:'65e3c3dbae70ae1db09293710ff888e1ee8f18b4f30f468e14e7f0cafc4441bb',
    COUNT:'75f2d0a56f759815644f298f060697a8d3be6afe12f38cc3fe1dff359ecca399',
    BASE:'c0c10822e322c31689b21a53b2bc64b03522dd7077fbcac1758dfdf19ddacb28',
}


def digest(p):return sha256(p.read_bytes()).hexdigest()


def terminal_domain(total):
    def high(left,least=9):
        if left==0:
            yield ()
        for r in range(least,left+1):
            for p in high(left-r,r):yield (r,)+p
    return sorted((low,)+p if low else p for low in range(9) for p in high(total-low))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in (LOCAL,COUNT,BASE,TAIL)}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert digest(ROOT/p)==h,(stem,p)
    for stem,h in FINAL_REPORTS.items():assert digest(ROOT/f'results/{stem}.md')==h
    local,count,base,tail=(docs[s] for s in (LOCAL,COUNT,BASE,TAIL))
    base_source=f'results/{BASE}.json';base_report=f'results/{BASE}.md'
    assert local['sources'][base_source]==count['sources'][base_source]==digest(ROOT/base_source)
    assert local['sources'][base_report]==count['sources'][base_report]==FINAL_REPORTS[BASE]
    assert local['status']=='PASS_DEPTH6_FINITE7_ADAPTIVE_TRANSFER'
    assert count['status']=='PASS_CONDITIONAL_DEPTH6_SPARSE_AMBIENT60_FOREST_FRONTIER'
    for d in (local,base):assert (d['triangle_minimum'],d['quadrilateral_minimum'])==(2,1)
    assert not local['new_coupons'] and not local['new_forest_count']
    assert count['local_payment_geometry_proved_here'] is False
    assert count['geometric_hypothesis_and_full_payment_composition_required'] is True
    assert local['counts']['formulas']==count['conditional_all_rank_deleted_pair_formulas']
    assert local['counts']['formulas']==dict(A='max(0,n-19)',B='max(0,n-20)',D='max(0,n-21)')
    assert count['old_all_rank_deleted_pair_formulas']==dict(A='max(0,n-23)',B='max(0,n-24)',D='max(0,n-25)')
    assert local['profiles']['candidate_count']==129 and local['profiles']['accepted_count']==18
    assert local['minimal_core_reduction']['catalogue_diagrams']==26
    assert [(r['diagram'],r['a'],r['c']) for r in local['minimal_core_reduction']['survivors']]==[(25,8,1),(25,8,2)]
    assert (local['schur']['branch_count'],local['schur']['positive_branches'])==(766,765)
    assert len(local['schur']['central_exception'])==1
    assert local['endpoint_seal']['norm2_projection_lower']=='5/2'
    assert local['endpoint_seal']['retained_whole_alternatives']==[
        'finite_A8','finite_D9','affine_D9_order10','affine_E8_order9']

    finite=count['finite']
    assert finite['all_cross_component_costs_unchanged'] and finite['new_residual_tables']
    assert (finite['complete_fixed_vectors'],finite['interaction_keys'])==(1044,1042)
    assert len(finite['alternatives'])==1044 and len(finite['tables'])==1042
    assert len({tuple(r['counts']) for r in finite['alternatives']})==1044
    assert len({tuple(t['key']) for t in finite['tables']})==1042
    assert max(r['octupled_upper'] for r in finite['alternatives'])==3232
    assert finite['exact_maximum']=='404' and finite['curvature']=='405' and finite['finite_slack']=='1'
    assert finite['above_curvature']==0
    assert finite['actual_attainer']==['A15']*4 and finite['actual_attaining_cost']=='404'
    assert finite['actual_rank15_costs']==dict(A15='101',B15='97',D15='201/2')
    assert finite['rank15_representative_upper']=='101'
    assert 4*F(finite['actual_rank15_costs']['A15'])==F(finite['exact_maximum'])
    # The full recurrence certificate is proved/checked in COUNT; no DP is
    # rerun here. In particular residual D17 is not omitted: fixed17 is A17.
    assert all(len(t['values'])==len(t['last'])==61 for t in finite['tables'])

    ideal=count['ideal'];rows=ideal['rows']
    assert ideal['ideal_weights_unchanged'] and ideal['low9_concentration_only']
    assert ideal['inherited_order17_cost']=='229/2'
    assert rows==base['ideal']['rows'] and ideal['summary']==base['ideal']['summary']
    domain=terminal_domain(59)
    assert len(domain)==len(rows)==ideal['complete_profiles']==2113
    assert sorted(tuple(r['ranks']) for r in rows)==domain
    assert len({tuple(r['ranks']) for r in rows})==2113
    for row in rows:
        h=row['h_polynomial'];v=prod(r+1 for r in row['ranks'])
        assert len(h)==60 and h==h[::-1] and all(x>=1 for x in h)
        assert v==sum(h)==row['link_vertices']
        g=sum(min(j,59-j)*x for j,x in enumerate(h));assert g==row['G']
        kappa=F(5*g-59*v,2);assert kappa==F(row['curvature'])
        margin=kappa/v-F(row['normalized_cost']);assert margin==F(row['normalized_margin'])>0
        assert v*margin==F(row['absolute_margin'])>=54
    norm=min(F(r['normalized_margin']) for r in rows)
    absolute=min(F(r['absolute_margin']) for r in rows)
    assert norm==F(31,6000) and absolute==54
    assert ideal['positive_terminal_absolute_surplus_extends_to_all_actual_cusps']

    # Bind the inherited global shape-curvature constant, and independently
    # check the dimension/parity arithmetic. The extremum correction is once.
    inherited=base['necessary_ambient60_inequality']
    assert inherited==dict(finite_vertex_coefficient='5/2',constant='4500',ideal_vertex_coefficient='54')
    d=59;g0=sum(min(j,d-j) for j in range(d+1))
    beta=F(5*g0-d*(d+1),2)
    extremum_difference=sum(2*j-d for j in range(d+1) if 2*j>d)
    constant=5*extremum_difference
    assert (g0,beta,extremum_difference,constant)==(870,F(405),900,4500)
    assert constant==F(inherited['constant'])==F(count['strict_global_constant'])
    assert F(base['finite']['exact_maximum'])-beta==F(inherited['finite_vertex_coefficient'])
    assert absolute==F(inherited['ideal_vertex_coefficient'])
    assert beta-F(finite['exact_maximum'])==1
    assert tail['dimension_bound']==60 and tail['hypotheses']==[
        'finite-sided real hyperbolic Coxeter polytope','finite volume']
    paths=[Path(__file__),ROOT/CLARIFICATION]+[ROOT/f'{folder}/{s}.{ext}' for s in docs
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_COMPLETE_DEPTH6_DIMENSION59_COMPOSITION',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        hypotheses=tail['hypotheses'],excluded_dimension=60,inherited_dimension_bound=60,dimension_bound=59,
        arithmeticity_required=False,compactness_required=False,simple_cusps_required=False,
        facet_count_bound_required=False,triangle_minimum=2,quadrilateral_minimum=1,
        feature_depth=6,changed_distance=7,all_cross_component_coupons_unchanged=True,
        finite_maximum='404',finite_curvature='405',finite_slack='1',actual_attaining_forest=['A15']*4,
        complete_fixed_vectors=1044,interaction_keys=1042,complete_terminal_ideal_profiles=2113,
        normalized_ideal_minimum=str(norm),actual_ideal_absolute_surplus=str(absolute),
        normalized_minimizers=[r['ranks'] for r in rows if F(r['normalized_margin'])==norm],
        absolute_minimizers=[r['ranks'] for r in rows if F(r['absolute_margin'])==absolute],
        low9_concentration_unchanged=True,high_factor_concentration_used=False,
        strict_global_constant=str(constant),constant_counted_once=True,
        contradiction='0 >= 1*v_f + 54*v_i + 4500 > 0',
        local_composition='If no actual corner has a newly deleted finite distance7 unit, use the frozen sparse rule pointwise. Otherwise the depth6 theorem gives an independent short6 payment or a retained whole A8/D9/affineD9/E8 unit at a distinct corner. Triangles retain their protected payments; all typed coupons and ideal weights are unchanged.',
        conditional_count_hypothesis_discharged=True,
        conditional_report_erratum='Residual ranks are1 through8 and17 through60, including D17; fixed rank17 is A17. The conditional source and complete certificate already include residual17, so this is a report-only correction.',
        higher_dimension_tail='The completed bound60 excludes every dimension at least61; the present contradiction excludes60.',
        new_finite_volume_example=False,external_refereeing=False,formal_proof_assistant_verified=False,
        scope='Complete mathematical composition with frozen local proof and finite/ideal certificates, independent terminal-domain and saved-polynomial checks, and exact inherited curvature arithmetic. No geometric census, convolution, Bellman optimization or construction is rerun.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'dimension',out['dimension_bound'],'finite',out['finite_maximum'],
        'ideal absolute',out['actual_ideal_absolute_surplus'],'constant',out['strict_global_constant'],flush=True)
