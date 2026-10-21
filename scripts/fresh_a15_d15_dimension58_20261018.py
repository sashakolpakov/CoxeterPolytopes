#!/usr/bin/env python3
"""Compose frozen rank15 local transfers and the complete ambient59 count.

Independent saved-domain, polynomial and curvature checks; no new geometric
census, convolution, optimizer or Bellman recurrence is run in this binder.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_a15_d15_dimension58_20261018'
A15='fresh_finite_a15_three_coupon_zero_transfer_20261018'
D15='fresh_finite_d15_a4_a11_half_zero_transfer_20261018'
COUNT='fresh_a15_d15_ambient59_forest_closure_20261018'
BASE='fresh_a14_cd16_ambient59_frontier_20261018'
TAIL='fresh_depth6_dimension59_20261018'
FINAL_REPORTS={A15:'90b45e20ba0835656a32481f8d42eb47bf2e96d82c38da1c986e1f0fd30efaa2',
 D15:'5fc4e6622636fbce7dcc3f349028b657d2f868d2593ebfadaf51488c97adfe61',
 COUNT:'021b4e474b3a15b805788d8736ff732bd102c0f93d8f5921c148a45720a2607c',
 BASE:'257219a7c1e324f44d4c7ed1c39d846a0994da807d2328f07b2387bad9ddb6cf',
 TAIL:'9bc6ca854bc7bceb8e8f5a50915f391c0cec345f901be09e251d09917069cdd8'}


def digest(p):return sha256(p.read_bytes()).hexdigest()
def fixed_domain(left=59,r=9,prefix=()):
    if r==18:
        yield prefix;return
    for k in range(left//r+1):yield from fixed_domain(left-r*k,r+1,prefix+(k,))
def high_parts(left,least=9):
    if left==0:yield()
    for r in range(least,left+1):
        for p in high_parts(left-r,r):yield(r,)+p


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(A15,D15,COUNT,BASE,TAIL)}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(ROOT/p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(ROOT/f'results/{s}.md')==h,s
    aa,dd,count,base,tail=(docs[s]for s in(A15,D15,COUNT,BASE,TAIL))
    for ext in('json','md'):
        path=f'results/{BASE}.{ext}';h=digest(ROOT/path)
        assert aa['sources'][path]==dd['sources'][path]==count['sources'][path]==h
    for local in(A15,D15):
        for ext in('json','md'):
            path=f'results/{local}.{ext}';assert count['sources'][path]==digest(ROOT/path)
    assert aa['changed_pairs']==[[0,12],[2,14]]and dd['changed_pairs']==[[0,11]]
    assert aa['new_A15_cost']==dd['new_cost']=='100'
    assert aa['retained_A15_unit']==[1,13]and aa['retained_A16_halves']==[[1,13],[2,14]]
    assert dd['retained_D16_half']==[1,12]
    assert aa['complete_selected_support']['complete_patterns']==dd['complete_selected_support']['complete_patterns']==536
    assert aa['complete_selected_support']['complete_maps']==482 and dd['complete_selected_support']['complete_maps']==870
    assert aa['complete_selected_support']['affected_maps']==8 and dd['complete_selected_support']['affected_maps']==10
    assert aa['all_ideal_weights_unchanged']and dd['ideal_weights_unchanged']
    for d in(aa,dd,count,base):assert(d['triangle_minimum'],d['quadrilateral_minimum'])==(2,1)
    assert[cp['component_pair_increment']for cp in aa['new_coupons']]==['1/2','1/2','2']
    assert dd['new_coupon']['component_pair_increment']=='2'
    ll=count['combined_literal_interface']
    assert ll['configurations']==len(ll['rows'])==216
    assert(ll['changed_debit_corners'],ll['changed_credit_corners'])==(0,0)
    assert len({(r['source'],json.dumps(r['case'],sort_keys=True))for r in ll['rows']})==216
    for r in ll['rows']:
        assert sum(F(c['retained_payment_lower_bound'])for c in r['corners'])==F(r['total_lower_bound'])>=1

    ff=count['finite'];domain=sorted(fixed_domain())
    assert len(domain)==len(set(domain))==ff['fixed_large_vectors']==961
    assert sorted(tuple(r['counts'])for r in ff['alternatives'])==domain
    expanded={(cs,k5)for cs in domain for k5 in range((59-sum(r*k for r,k in zip(range(9,18),cs)))//5+1)}
    assert len(expanded)==ff['fixed_large_and_A5_vectors']==2574
    assert {(tuple(r['counts']),r['fixed_A5_count'])for r in ff['expanded_alternatives']}==expanded
    assert len(ff['expanded_alternatives'])==2574
    assert ff['rank15_actual_costs']==dict(A15='100',B15='97',D15='100')
    assert ff['rank15_representative_upper']=='100'
    assert ff['certified_upper']==ff['curvature']=='783/2'and ff['upper_slack']=='0'
    assert not ff['actual_maximum_claimed']and not ff['actual_attainment_of_global_upper_claimed']
    assert max(r['sixteenths_upper']for r in ff['expanded_alternatives'])==6264
    assert all(r['sixteenths_actual']<=r['sixteenths_upper']<=6264 for r in ff['expanded_alternatives'])
    assert ff['above_curvature']==[] and len(ff['tables'])==5148
    assert ff['interaction_keys']==2574 and ff['simultaneous_type_pairs']==32400
    assert ff['residual_item_checks']==128700 and ff['Bellman_upper_inequalities']==6723288
    assert ff['attaining_predecessors']==303732 and ff['context_options']==9191
    assert ff['residual_rank17_D17_fixed_rank17_A17']and ff['contextual_witness_multiplicity_not_counted']
    assert len({(tuple(t['key']),t['fixed_A5_count'],t['exclude2'])for t in ff['tables']})==5148
    assert all(len(t['values'])==len(t['last'])==60 and t['values'][0]==t['last'][0]==0 for t in ff['tables'])

    ideal=count['ideal'];assert ideal==base['ideal']
    rows=ideal['rows'];domain=sorted((low,)+p if low else p for low in range(9)for p in high_parts(58-low))
    assert len(domain)==len(set(domain))==len(rows)==ideal['complete_profiles']==1890
    assert sorted(tuple(r['ranks'])for r in rows)==domain
    env={int(k):F(v)for k,v in ideal['current_envelopes'].items()}
    for r in rows:
        h=r['h_polynomial'];rs=r['ranks'];v=prod(x+1 for x in rs)
        assert sum(rs)==58 and len(h)==59 and h==h[::-1]and all(isinstance(x,int)and x>0 for x in h)
        assert v==sum(h)==r['link_vertices']
        g=sum(min(j,58-j)*x for j,x in enumerate(h));assert g==r['G']
        assert 2*g==58*v-sum(abs(2*j-58)*x for j,x in enumerate(h))
        curv=F(5*g-58*v,2);cost=sum(env[n]for n in rs)
        assert curv==F(r['curvature'])and cost==F(r['normalized_cost'])
        margin=curv/v-cost;assert margin==F(r['normalized_margin'])>0
        assert v*margin==F(r['absolute_margin'])>=F(91,2)
    norm=min(F(r['normalized_margin'])for r in rows);absolute=min(F(r['absolute_margin'])for r in rows)
    assert norm==F(1,340)and absolute==F(91,2)
    assert ideal['positive_terminal_extension_to_actual_cusps']and ideal['actual_charge_bounded_by_envelope']

    d=58;g0=sum(min(j,d-j)for j in range(d+1));beta=F(5*g0-d*(d+1),2)
    correction=sum(2*j-d for j in range(d+1)if 2*j>d);constant=5*correction
    assert(g0,beta,correction,constant)==(841,F(783,2),870,4350)
    assert beta==F(count['finite_curvature'])==F(base['finite_curvature'])
    assert constant==int(count['strict_global_constant'])==int(base['strict_global_constant'])
    assert base['necessary_ambient59_inequality']=='2*v_f >= 4350 + (91/2)*v_i'
    assert F(base['finite']['exact_maximum'])-beta==2
    assert count['necessary_ambient59_inequality']=='0 >= 4350 + (91/2)*v_i'
    assert tail['dimension_bound']==59 and tail['hypotheses']==[
        'finite-sided real hyperbolic Coxeter polytope','finite volume']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in docs for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_COMPLETE_A15_D15_DIMENSION58_COMPOSITION',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},hypotheses=tail['hypotheses'],
        excluded_dimension=59,inherited_dimension_bound=59,dimension_bound=58,
        arithmeticity_required=False,compactness_required=False,simple_cusps_required=False,facet_count_bound_required=False,
        triangle_minimum=2,quadrilateral_minimum=1,finite_certified_upper='783/2',finite_curvature='783/2',finite_slack='0',
        actual_finite_maximum_claimed=False,complete_fixed_vectors=961,complete_expanded_vectors=2574,residual_tables=5148,
        complete_terminal_ideal_profiles=1890,normalized_ideal_minimum=str(norm),actual_ideal_absolute_surplus=str(absolute),
        strict_global_constant=str(constant),constant_counted_once=True,
        contradiction='0 >= 4350 + (91/2)*v_i > 0',
        local_composition='Independent prior finite6 and ideal/short-path actual-corner priorities first; remaining changed selected A15/D15 contexts use their complete arbitrary-support proofs and distinct retained A16/D16 halves plus explicitly counted whole-component coupons. All other selected/nonselected and triangle outputs retain their weights.',
        finite_representative_upper_not_actual_attainment=True,low9_concentration_unchanged=True,high_factor_concentration_used=False,
        higher_dimension_tail='The completed bound59 excludes every dimension at least60; the present contradiction excludes59.',
        new_finite_volume_example=False,external_refereeing=False,formal_proof_assistant_verified=False,
        scope='Complete mathematical composition with frozen local theorems and complete finite/ideal certificate. Independent saved-domain, polynomial and inherited curvature arithmetic checks; no geometric census, convolution, Bellman optimization or construction replay.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'dimension',out['dimension_bound'],'finite upper',out['finite_certified_upper'],
        'ideal surplus',out['actual_ideal_absolute_surplus'],'constant',out['strict_global_constant'],flush=True)
