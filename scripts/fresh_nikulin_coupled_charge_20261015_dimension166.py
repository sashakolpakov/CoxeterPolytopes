#!/usr/bin/env python3
"""Bind the reviewed common-ratio estimate and complete dimension166 split."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_coupled_charge_20261015_dimension166'
BASE='fresh_nikulin_unit_charge_20261015_dimension169'
PRIMARY='fresh_nikulin_coupled_charge_20261015'
REVIEW=PRIMARY+'_review'
FACET='fresh_nikulin_narrow_20261015_compact_facet'

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    tags=[BASE,PRIMARY,REVIEW,FACET]
    data={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in tags}
    assert data[BASE]['status']=='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_169_CANDIDATE'
    assert data[PRIMARY]['status']=='PASS_EXACT_BALANCED_SEVEN_SIXTHS_CHARGE_CANDIDATE166'
    assert data[REVIEW]['status']=='PASS_INDEPENDENT_GROUNDED_RATIO_DEPENDENT_CHARGE_AND_166_REVIEW'
    for p in (f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json'):
        assert data[REVIEW]['sources'][p]==digest(p)
    assert Q(data[PRIMARY]['charge'])==Q(data[REVIEW]['charge'])==Q(7,6)
    assert Q(data[PRIMARY]['common_height_ratio_limit'])==Q(data[REVIEW]['height_ratio_limit'])==Q(21,20)
    g=data[PRIMARY]['global_statement']
    assert g['total_endpoint_cost_budget']==4 and Q(g['avoid_penalty'])==Q(2,9)
    assert g['combined_candidate_at_most']==data[REVIEW]['combined_candidate_bound']==166
    assert data[FACET]['ambient_dimension_at_most']==143
    assert (Q(2*166,9)-Q(7,6)*4-Q(2,9))/2-2==14
    assert Q(21,20)*Q(166-43,9)==Q(287,20)>14
    assert data[REVIEW]['first_excluded_dimension']==167
    files=[f'scripts/{TAG}.py',f'results/{TAG}.md']
    for t in tags:files.extend([f'scripts/{t}.py',f'results/{t}.json',f'results/{t}.md'])
    seen,checked=set(),set()
    def recurse(path):
        if path in seen:return
        seen.add(path)
        for p,h in json.loads((ROOT/path).read_text()).get('sources',{}).items():
            assert digest(p)==h,(path,p)
            checked.add(p)
            if p.endswith('.json'):recurse(p)
    for t in tags:recurse(f'results/{t}.json')
    return dict(status='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_166_CANDIDATE',
                sources={p:digest(p) for p in files},recursive_manifests_checked=len(seen),
                recursive_source_files_checked=len(checked),
                candidate_statement='Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most166.',
                arithmetic_assumption=False,compactness_assumption=False,simple_cusp_assumption=False,
                chosen_narrow_facet_cases=[
                    dict(case='contains an ideal vertex',ambient_dimension_at_most=166,primary=PRIMARY,independent_review=REVIEW),
                    dict(case='only finite vertices, hence compact',ambient_dimension_at_most=143,primary=FACET,
                         independent_review='fresh_dimension225_hostile_20261015')],
                cost_accounting=dict(total_schur_cost_at_most=4,loss_per_cost='7/6 if M<=21/20; 4/3 otherwise',avoid_node_deficit='2/9'),
                common_height_ratio_split='21/20',first_excluded_dimension=167,
                near_balanced_pairing_lower='M*((n-1)-40)/9',other_pairing_lower='M*((n-1)-43)/9',
                combined_candidate_dimension_at_most=166,norm_two_narrow_threshold=14,frozen_baseline=BASE,
                primary_sources=data[BASE]['primary_sources'],
                scope='Frozen proof/review dependency and complete case binder. No producer arithmetic imported. Internal research proof candidate; external validation, novelty and current-best status unestablished. Published geometry and classification inputs are mathematical dependencies, not formal proof-assistant certification.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
