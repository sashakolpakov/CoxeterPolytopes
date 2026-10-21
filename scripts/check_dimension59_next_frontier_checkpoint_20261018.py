#!/usr/bin/env python3
"""Integrity, links and saved-headline audit; no research computation replay."""
import argparse,ast,json,re
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension59_next_frontier_and_link_flip_20261018'
STEMS=(
 'fresh_depth6_dimension59_20261018',
 'fresh_ambient59_payment_frontier_20261018',
 'fresh_finite_rank14_long_pair_support_scout_20261018',
 'fresh_finite_rank14_existing_payment_transfer_20261018',
 'fresh_order14_c_basic_distance9_transfer_20261018',
 'fresh_raw76_complete_family_full_model_20261018',
 'fresh_raw76_complete_family_support_20261018',
 'fresh_joint44_65_compact_chain_20261018',
 'fresh_joint44_63_65_link_flip_20261018')


def run():
    hashes={}
    def digest(p):
        p=Path(p)
        if p not in hashes:
            h=sha256()
            with p.open('rb') as f:
                for block in iter(lambda:f.read(1<<20),b''):h.update(block)
            hashes[p]=h.hexdigest()
        return hashes[p]
    docs={}
    for s in STEMS:
        p=ROOT/f'results/{s}.json';docs[s]=json.loads(p.read_text());digest(p)
    bindings=0
    for stem,d in docs.items():
        assert d['status'].startswith(('PASS','BOUNDED')),stem
        maps=[d.get('sources',{}),d.get('prepared',{}).get('sources',{}),d.get('certification_sources',{})]
        assert any(maps)
        for mapping in maps:
            for p,h in mapping.items():assert digest(ROOT/p)==h,(stem,p);bindings+=1
    proof,front,scout,finite,ideal,model,support,construction,flip=[docs[s] for s in STEMS]
    assert proof['dimension_bound']==59 and proof['conditional_count_hypothesis_discharged']
    assert front['finite']['exact_maximum']=='795/2' and front['finite_curvature']=='783/2'
    assert front['finite']['complete_fixed_vectors']==961 and front['finite']['interaction_keys']==960
    assert front['ideal']['complete_profiles']==1890 and front['ideal']['summary']['negative_count']==10
    assert front['ideal']['summary']['normalized_minimum']=='-1537/14000'
    assert front['strict_global_constant']=='4350' and front['new_dimension_bound'] is None
    assert [r['complete_eligible_patterns'] for r in scout['records']]==[600]*3
    assert [r['complete_maps'] for r in scout['records']]==[552,260,1104]
    assert scout['new_weight_rule_proved'] is False
    assert finite['new_costs']==dict(A14='91',B14='87',D14='181/2')
    assert finite['selected_contexts']==13 and finite['saturated_extra_white_branches']==4
    assert finite['new_cross_component_coupons']==[] and finite['rank14_representative_upper_remains']=='91'
    assert finite['rejected_A14_shortcut']['retained_A14_cost']=='91'
    assert finite['combined_literal_interface']['changed_corners']=={'D':8}
    assert ideal['new_cost']=='80' and ideal['current_order14_family_costs']==dict(A=70,B=82,C=80,D=82)
    assert ideal['current_order14_maximum']=='82' and ideal['affected_selected_maps']==22 and ideal['affected_protected_maps']==4
    assert ideal['finite_weights_unchanged'] and ideal['new_cross_component_coupons']==[]
    for theorem in (finite,ideal):
        assert theorem['triangle_minimum']==2 and theorem['quadrilateral_minimum']==1
        assert theorem['new_dimension_bound'] is None
    a=finite['combined_literal_interface']['rows'];b=ideal['current_literal_interface']['rows']
    assert len(a)==len(b)==216
    for x,y in zip(a,b):
        assert (x['source'],x['case'])==(y['source'],y['case'])
        assert F(x['total_lower_bound'])>=1 and F(y['total_lower_bound'])>=1
        assert len(x['corners'])==len(y['corners'])
        for p,q in zip(x['corners'],y['corners']):
            assert p['marks']==q['marks'] and q['new_debit']=='0'
            assert F(p['previous_lower'])==F(q['retained_payment_lower_bound'])
    # The old exact upper remains valid after decreases. Its actual attainers
    # use no changed rank14 component, so the upper is still attained.
    witnesses=front['finite']['actual_obstructing_forests'][:3]
    assert len(witnesses)==3
    assert all(F(w['cost'])==F(795,2) and all(not t.endswith('14') for t in w['forest']) for w in witnesses)
    assert 3*F(101)+F(finite['new_costs']['D14'])==F(787,2)>F(783,2)
    assert model['solver_invocations']==1 and model['solver']['status']==1 and model['solver']['objective'] is None
    assert model['optimization_resolved'] is False and model['changed_face_masses'] is False
    lower=F(model['exact_interval']['lower']);upper=F(model['exact_interval']['upper'])
    assert lower==F(102699323852930280,8367588524424317)>12
    assert upper==F(6174801959947,500000000000)>lower
    assert F(model['exact_interval']['gap'])==upper-lower
    counts=model['prepared']['counts']
    assert counts['variables']==440098 and counts['appended_family_columns']==242
    assert counts['preserved_multi_owner_columns']==14 and counts['zero_budget_signed_variables']==96176
    assert support['counts']['positive']==4624 and support['counts']['unambiguous_missing_classes']==1288
    assert construction['summary']['changed_triples']==0 and construction['summary']['native_pairs']==345
    assert construction['summary']['unresolved']==457
    assert flip['summary']['selected_prefixes']==21 and flip['summary']['X_branches']==1584 and flip['summary']['X_native']==270
    assert flip['summary']['X_L64_finite']==flip['summary']['R_banks']==flip['summary']['completed_disks']==0
    assert flip['summary']['unresolved']==292
    assert construction['additional63_deletion_debit']==flip['exact_deletion_debit']
    assert flip['exact_deletion_debit']['lost_finite_count']==19
    assert flip['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    sources=[ROOT/f'scripts/{s}.py' for s in STEMS]
    sources.extend(ROOT/f'scripts/{s}.py' for s in (
        'certify_raw76_complete_family_fallback_20261018',
        'check_raw76_complete_family_model_saved_20261018',
        'check_joint44_65_and63_flip_saved_20261018'))
    sources.append(Path(__file__))
    for p in sources:ast.parse(p.read_text());digest(p)
    reports=[ROOT/f'results/{s}.md' for s in STEMS if s!=STEMS[6]]+[ROOT/f'results/{TAG}.md']
    links=0
    def check_links(p,body):
        nonlocal links
        for target in re.findall(r'\]\(([^)]+)\)',body):
            if target.startswith(('http:','https:','#')):continue
            dest=(p.parent/target.split('#')[0]).resolve()
            if dest==ROOT/f'results/{TAG}_audit.json':continue
            assert dest.exists(),(p,target);links+=1
    for p in reports:check_links(p,p.read_text());digest(p)
    for name,stop in (('README.md','The [preceding bound61 checkpoint]'),
                      ('RESEARCH_PLAN.md','The [preceding bound62 checkpoint]'),('GOAL.md','## Objective')):
        p=ROOT/name;body=p.read_text().split(stop,1)[0]
        assert TAG in body and '59' in body and '21' in body
        check_links(p,body);digest(p)
    return dict(status='PASS_DIMENSION59_NEXT_FRONTIER_CHECKPOINT_INTEGRITY_AND_HEADLINES',
        packets=len(STEMS),direct_bindings=bindings,files_hashed=len(hashes),source_AST_checks=len(sources),
        reports_and_openings=len(reports)+3,local_links=links,
        sources={str(p.relative_to(ROOT)):h for p,h in sorted(hashes.items())},
        saved_headlines=dict(dimension_bound=59,finite_rank59_maximum='795/2',curvature='783/2',
            negative_ideal_terminal_profiles=10,finite_rank14=finite['new_costs'],ideal_order14=ideal['current_order14_family_costs'],
            arithmetic_lower=str(lower),arithmetic_upper=str(upper),native_flip_caps=270,new_finite_volume_example=False),
        scope='Integrity, syntax, local links and saved headline arithmetic only. Does not replay geometric census, native searches, inertia/radical computation, finite Bellman certificates, ideal convolution, numerical cone reconstruction, full-matrix optimization/replay or formal proof verification.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}_audit.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],{k:v for k,v in out.items() if k not in ('sources','scope','saved_headlines')})
