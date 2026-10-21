#!/usr/bin/env python3
"""Bind the reviewed four-thirds lemma and complete dimension169 case split."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_unit_charge_20261015_dimension169'
BASE='fresh_nikulin_dimension_181_candidate_20261015'
PRIMARY='fresh_nikulin_unit_charge_20261015'
REVIEW=PRIMARY+'_review'
FACET='fresh_nikulin_narrow_20261015_compact_facet'


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    tags=[BASE,PRIMARY,REVIEW,FACET]
    data={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in tags}
    assert data[BASE]['status']=='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_181_CANDIDATE'
    assert data[PRIMARY]['status']=='PASS_EXACT_UNIT_CHARGE_OBSTRUCTION_AND_FOUR_THIRDS_UNIFORM_LEMMA'
    assert data[REVIEW]['status']=='PASS_INDEPENDENT_GROUNDED_ENERGY_ALL_RANK_AND_FOUR_THIRDS_REVIEW'
    for p in (f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json'):
        assert data[REVIEW]['sources'][p]==digest(p)
    assert Q(data[PRIMARY]['charge_coefficient'])==Q(4,3)
    assert data[PRIMARY]['endpoint_cost_budget']==data[REVIEW]['total_schur_cost_bound']==4
    assert Q(data[PRIMARY]['avoid_penalty'])==Q(data[REVIEW]['avoid_node_penalty'])==Q(2,9)
    assert data[PRIMARY]['combined_candidate_dimension_at_most']==data[REVIEW]['combined_candidate_dimension_at_most']==169
    assert data[FACET]['ambient_dimension_at_most']==143
    assert (Q(2*169,9)-Q(4,3)*4-Q(2,9))/2-2==14
    assert data[REVIEW]['first_excluded_dimension']==170
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
    return dict(status='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_169_CANDIDATE',
                sources={p:digest(p) for p in files},recursive_manifests_checked=len(seen),
                recursive_source_files_checked=len(checked),
                candidate_statement='Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most169.',
                arithmetic_assumption=False,compactness_assumption=False,simple_cusp_assumption=False,
                chosen_narrow_facet_cases=[
                    dict(case='contains an ideal vertex',ambient_dimension_at_most=169,primary=PRIMARY,independent_review=REVIEW),
                    dict(case='only finite vertices, hence compact',ambient_dimension_at_most=143,primary=FACET,
                         independent_review='fresh_dimension225_hostile_20261015')],
                cost_accounting=dict(total_schur_cost_at_most=4,loss_per_cost='4/3',avoid_node_deficit='2/9'),
                negative_pairing_lower_bound='M*((n-1)-43)/9, M>=1',first_excluded_dimension=170,
                combined_candidate_dimension_at_most=169,norm_two_narrow_threshold=14,frozen_baseline=BASE,
                primary_sources=data[BASE]['primary_sources'],
                scope='Frozen proof/review dependency and complete case binder. No producer arithmetic imported. Internal research proof candidate; external validation, novelty and current-best status unestablished. Published geometry and classification inputs are mathematical dependencies, not formal proof-assistant certification.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
