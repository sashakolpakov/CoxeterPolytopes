#!/usr/bin/env python3
"""Final dependency/case binder for the reviewed 189 research proof candidate."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_dimension_189_candidate_20261015'
BASE='fresh_nikulin_dimension_225_candidate_20261015'
CHARGE='fresh_nikulin_energy_charge_20261015'
REVIEW='fresh_nikulin_energy_charge_20261015_review'
FACET='fresh_nikulin_narrow_20261015_compact_facet'


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    tags=[BASE,CHARGE,REVIEW,FACET]
    data={tag:json.loads((ROOT/'results'/f'{tag}.json').read_text()) for tag in tags}
    assert data[BASE]['status']=='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_225_CANDIDATE'
    assert data[CHARGE]['status']=='PASS_EXACT_CONTACT_COST_CHARGE_FOR_DIMENSION_189_CANDIDATE'
    assert data[REVIEW]['status']=='PASS_INDEPENDENT_FRACTION_ENERGY_CHARGE_AND_HOSTILE_SCOPE_REVIEW'
    assert data[REVIEW]['checked_candidate_dimension_at_most']==189
    for path in [f'scripts/{CHARGE}.py',f'results/{CHARGE}.json']:
        assert data[REVIEW]['sources'][path]==digest(path)
    assert data[CHARGE]['non_narrow_loss_per_unit_schur_cost']==2
    assert data[CHARGE]['total_two_endpoint_schur_cost_at_most']==4
    assert data[CHARGE]['separate_narrow_component_loss_at_most']==2
    assert data[CHARGE]['ideal_narrow_facet_dimension_at_most']==189
    assert data[CHARGE]['first_excluded_dimension']==190
    assert data[FACET]['ambient_dimension_at_most']==143
    assert Q(2*(190-1),9)-2*4-2==32
    assert Q(32,2)-2==14
    assert max(189,143)==189

    sources=[f'scripts/{TAG}.py',f'results/{TAG}.md']
    for tag in tags:
        sources.extend([f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md'])
    manifests=set();checked=set()
    def recurse(path):
        if path in manifests:return
        manifests.add(path)
        obj=json.loads((ROOT/path).read_text())
        for filename,expected in obj.get('sources',{}).items():
            assert digest(filename)==expected,(path,filename)
            checked.add(filename)
            if filename.endswith('.json'):recurse(filename)
    for tag in tags:recurse(f'results/{tag}.json')
    return {
        'status':'PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_189_CANDIDATE',
        'sources':{path:digest(path) for path in sources},
        'recursive_manifests_checked':len(manifests),
        'recursive_source_files_checked':len(checked),
        'candidate_statement':'Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most189.',
        'arithmetic_assumption':False,
        'review_status':data[REVIEW]['status'],
        'chosen_narrow_facet_cases':[
            {'case':'contains an ideal vertex','ambient_dimension_at_most':189,
             'primary':CHARGE,'independent_review':REVIEW},
            {'case':'has only finite vertices, hence is compact','ambient_dimension_at_most':143,
             'primary':FACET,'independent_review':'fresh_dimension225_hostile_20261015'},
        ],
        'energy_identity':'-x.y=(sum sigma_a*tau_a*E_a)/2-q-1/q',
        'cost_accounting':{'total_schur_cost_at_most':4,'loss_per_cost':2,'narrow_factor_extra_loss_at_most':2},
        'negative_pairing_lower_bound':'M*((n-1)/9-7), M>=1',
        'first_excluded_dimension':190,
        'norm_two_narrow_threshold':14,
        'frozen_baseline':BASE,
        'historical_pending_review_sentence':'The primary energy-charge report was frozen before its independent review; that pending-review statement describes the earlier checkpoint. The final review is bound here.',
        'literature_status':'Nikulin ICM1986 printed659 announces n<300. This is an explicit reviewed research proof candidate, without publication-priority or current-best assertion.',
        'primary_sources':data[BASE]['primary_sources'],
        'scope':'No producer arithmetic imported. This checker binds frozen proofs and independent reviews, validates recursive hashes, and checks the scalar/case combination. Published classifications and geometric arguments are mathematical dependencies, not a formal proof-assistant development.',
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();result=run();destination=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert json.loads(destination.read_text())==result
    else:destination.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
