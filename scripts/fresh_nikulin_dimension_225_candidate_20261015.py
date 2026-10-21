#!/usr/bin/env python3
"""Bind the two reviewed geometric branches; no producer arithmetic imports."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_dimension_225_candidate_20261015'
PACKAGES = [
    'fresh_nikulin_narrow_20261015',
    'fresh_nikulin_narrow_20261015_products',
    'fresh_nikulin_product_independent_20261015',
    'fresh_nikulin_narrow_20261015_compact_facet',
    'fresh_dimension225_hostile_20261015',
    'dimension_bound_prokhorov_20261015_proof_review',
]


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    data = {tag:json.loads((ROOT/'results'/f'{tag}.json').read_text()) for tag in PACKAGES}
    product = data['fresh_nikulin_narrow_20261015_products']
    facet = data['fresh_nikulin_narrow_20261015_compact_facet']
    assert product['status'] == 'PASS_CONDITIONAL_PRODUCT_CUSP_BOUND_225'
    assert product['first_excluded_dimension'] == 226
    assert product['bad_components_at_most'] == 9
    assert product['retained_nonright_per_endpoint'] == 4
    assert facet['status'] == 'PASS_EXACT_CHECKS_FOR_COMPACT_FACET_BOUND_143'
    assert facet['ambient_dimension_at_most'] == 143
    assert facet['compact_facet_dimension_at_most'] == 142
    assert facet['quadrilateral_path_bound'] == 17
    assert data['fresh_nikulin_product_independent_20261015']['status'] == 'PASS_INDEPENDENT_CONDITIONAL_PRODUCT_CUSP_ARITHMETIC'
    assert data['fresh_dimension225_hostile_20261015']['status'] == 'PASS_HOSTILE_AUDIT_WITH_COARSE_COMPACT_FACET_REPAIR'
    assert data['dimension_bound_prokhorov_20261015_proof_review']['status'] == 'PASS_PRIMARY_GEOMETRY_AND_SCOPE_AUDIT'
    assert Q(226-1,9)-11 == 14
    assert max(225, facet['ambient_dimension_at_most']) == 225

    bound = [f'scripts/{TAG}.py',f'results/{TAG}.md']
    for tag in PACKAGES:
        bound.extend([f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md'])
    visited = set(); checked_files = set()
    def verify_manifest(path):
        if path in visited:
            return
        visited.add(path)
        obj=json.loads((ROOT/path).read_text())
        for source, expected in obj.get('sources',{}).items():
            assert digest(source) == expected, (path,source)
            checked_files.add(source)
            if source.endswith('.json'):
                verify_manifest(source)
    for tag in PACKAGES:
        verify_manifest(f'results/{tag}.json')
    return {
        'status':'PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_225_CANDIDATE',
        'sources':{p:digest(p) for p in bound},
        'recursive_manifests_checked':len(visited),
        'recursive_source_files_checked':len(checked_files),
        'candidate_statement':'Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most225.',
        'arithmetic_assumption':False,
        'chosen_narrow_facet_cases':[
            {'case':'contains an ideal vertex','ambient_dimension_at_most':225,
             'proof':'fresh_nikulin_narrow_20261015_products',
             'independent_arithmetic':'fresh_nikulin_product_independent_20261015'},
            {'case':'has only finite vertices, hence is compact','ambient_dimension_at_most':143,
             'proof':'fresh_nikulin_narrow_20261015_compact_facet'},
        ],
        'geometric_review':'fresh_dimension225_hostile_20261015',
        'first_excluded_dimension':226,
        'narrow_pairing_threshold':'14 for outward normals of square2',
        'historical_status':'Nikulin ICM1986 printed659 announces n<300. This is an explicit research proof candidate, not a publication-priority/current-best assertion.',
        'primary_sources':[
            {'url':'https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf',
             'pages':'printed655 TheoremA; printed657-659 method and announcement'},
            {'url':'https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng',
             'pages':'printed407-408, section2 Propositions1-2; quasi-Lanner extraction/order10'},
            {'url':'https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt',
             'pages':'printed70-71 h-vector; printed77 Proposition10 spherical distance-pair bound'},
        ],
        'checker_scope':'Independent manifest closure and elementary case/scalar checks; imports no producer arithmetic. Geometric proofs and published classifications remain explicit mathematical dependencies, not a formal proof-assistant verification.',
    }


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true')
    args=parser.parse_args(); result=run(); destination=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(destination.read_text()) == result
    else:
        destination.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
