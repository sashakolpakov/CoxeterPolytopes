#!/usr/bin/env python3
"""Bind proved local transfers to the exact rank64 finite maximum464.

The remaining all-family ideal deficits are explicit; no dimension63 claim.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ambient64_finite_closure_20261018'
FRONT='fresh_ambient64_payment_frontier_20261018'
ZERO='fresh_depth11_finite8_zero_transfer_20261018'
EXCHANGE='fresh_finite_d16_d17_two_unit_exchange_20261018'
CYCLE='fresh_cycle11_half_pair_transfer_20261018'
CURRENT='fresh_deep_finite8_dimension64_20261018'
STEMS=(FRONT,ZERO,EXCHANGE,CYCLE,CURRENT)


def run():
    docs={stem:json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem,data in docs.items():
        assert data['status'].startswith('PASS'),stem
        for name,h in data['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    zero,exchange,cycle,front=(docs[s] for s in (ZERO,EXCHANGE,CYCLE,FRONT))
    assert zero['changed_weight']=='0' and zero['corridor']['feature_depth']==11
    assert zero['triangle_minimum']==2 and zero['quadrilateral_minimum']==1
    assert zero['ideal_weights_unchanged']
    assert exchange['new_finite_component_costs']==dict(A16=116,B16=104,D16=116,D17=127)
    assert exchange['depth11_zero_compatible'] and exchange['all_coupons_unchanged']
    assert exchange['all_ideal_weights_unchanged']
    assert exchange['literal_interfaces']['minimum_retained_total']=='1'
    assert cycle['added_finite_cost']=='0' and cycle['new_total_cost']=='99/2'
    assert cycle['other_ideal_types_unchanged'] and cycle['order11_all_family_envelope']=='55'
    scenarios=[r for r in front['finite_conditional_scenarios'] if
               (r['feature_depth'],r['deep_pair_weight'],r['D16_drop'],r['D17_drop'])==(11,'0',2,-2)]
    assert len(scenarios)==1; budget=scenarios[0]
    assert budget['exact_maximum']=='464' and budget['finite_slack']=='0'
    assert budget['attaining_forest']==['D16']*4 and F(4)*116==464
    assert budget['large_vectors']==753 and budget['interaction_keys']==348
    assert budget['checks']==dict(bellman_upper=597168,attaining_predecessors=22272,
                                 literal_interactions=42921,global_alternatives=753)
    assert budget['type_domination_checks']==38025
    assert docs[CURRENT]['dimension_bound']==64
    ideal=front['ideal_current_summary']
    assert ideal['count']==1513 and ideal['negative']==24 and ideal['zero']==0
    assert ideal['minimum']['normalized_margin']=='-62665/419628'
    assert front['strengthened_global_constant']=='5120'
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
                            for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_AMBIENT64_FINITE_PAYMENT_CLOSURE',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=64,finite_maximum='464',finite_curvature='464',finite_slack='0',
                attaining_forest=budget['attaining_forest'],finite_certificate_checks=budget['checks'],
                finite_budget_hypotheses_discharged=True,
                simultaneous_rules=dict(deep_finite8_depth=11,deep_finite8_weight='0',
                    D16=116,D17=127,affine_cycle11='99/2',new_coupons=[]),
                face_minima=dict(triangle=2,quadrilateral=1),
                ideal_envelope_unchanged=True,ideal_profiles=1513,negative_ideal_profiles=24,
                worst_ideal_ranks=ideal['minimum']['ranks'],
                worst_ideal_margin=ideal['minimum']['normalized_margin'],
                global_additive_constant='5120',completed_dimension_ceiling=64,
                new_dimension_bound=None,
                scope='The written composition discharges the conditional finite budget hypotheses. Its exact all-forest rank64 maximum464 is now proved for this complete payment rule. Cycle11 halves add no finite cost but leave the all-family ideal envelope55. The ideal obstruction remains; no dimension63 theorem or finite-volume example.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite_maximum'],'ideal negative profiles',out['negative_ideal_profiles'],flush=True)
