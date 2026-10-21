#!/usr/bin/env python3
"""Bind the reviewed aggregate estimate and complete dimension165 case split."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_height_charge_20261015_dimension165'
BASE = 'fresh_nikulin_coupled_charge_20261015_dimension166'
PRIMARY = 'fresh_nikulin_height_charge_20261015'
REVIEW = PRIMARY + '_review'
FACET = 'fresh_nikulin_narrow_20261015_compact_facet'


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def run():
    tags = [BASE, PRIMARY, REVIEW, FACET]
    data = {t: json.loads((ROOT / f'results/{t}.json').read_text()) for t in tags}
    assert data[BASE]['status'] == 'PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_166_CANDIDATE'
    assert data[PRIMARY]['status'] == 'PASS_EXACT_HEIGHT_CHARGE_AND_TWO_DEFECT_AGGREGATION_CANDIDATE165_PENDING_REVIEW'
    assert data[REVIEW]['status'] == 'PASS_INDEPENDENT_EXACT_HEIGHT_CHARGE_DEFECT_BUDGET_AND_165_REVIEW'
    for path in (f'scripts/{PRIMARY}.py', f'results/{PRIMARY}.json'):
        assert data[REVIEW]['sources'][path] == digest(path)
    assert Q(data[PRIMARY]['height_window']) == Q(data[REVIEW]['height_ratio_limit']) == Q(31, 30)
    p, r = data[PRIMARY]['aggregate'], data[REVIEW]['aggregate']
    assert p['total_cost_budget'] == r['total_cost'] == 4
    assert Q(r['minimum_defective_factor_cost']) == Q(3, 2)
    assert 3 * Q(3, 2) > 4
    assert p['maximum_defective_factors'] == r['maximum_defective_factors'] == 2
    assert Q(r['maximum_each_defect']) == Q(2, 9)
    assert Q(p['total_loss']) == Q(r['combined_cost_and_defect_loss']) == 4 + 2 * Q(2, 9)
    # N=n-1. The narrow-factor target loses 2/9 and the Lorentz formula 2M.
    near = (Q(2 * 165, 9) - Q(2, 9) - Q(40, 9)) / 2 - 2
    far = Q(31, 30) * Q(165 - 43, 9)
    assert near == Q(r['near_at_N165']) == 14
    assert far == Q(r['far_at_N165']) == 14 + Q(1, 135)
    assert Q(r['near_slope']) == Q(1, 9) > 0
    assert Q(r['far_slope']) == Q(31, 270) > 0
    assert p['combined_candidate'] == r['combined_candidate_bound'] == 165
    assert p['compact_facet_branch_at_most'] == r['compact_branch_inherited_bound'] == data[FACET]['ambient_dimension_at_most'] == 143
    assert (Q(2 * 164, 9) - Q(2, 9) - Q(40, 9)) / 2 - 2 < 14
    files = [f'scripts/{TAG}.py', f'results/{TAG}.md']
    for t in tags:
        files += [f'scripts/{t}.py', f'results/{t}.json', f'results/{t}.md']
    seen, checked = set(), set()

    def recurse(path):
        if path in seen:
            return
        seen.add(path)
        for source, expected in json.loads((ROOT / path).read_text()).get('sources', {}).items():
            assert digest(source) == expected, (path, source)
            checked.add(source)
            if source.endswith('.json'):
                recurse(source)

    for t in tags:
        recurse(f'results/{t}.json')
    return dict(
        status='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_165_CANDIDATE',
        sources={path: digest(path) for path in files},
        recursive_manifests_checked=len(seen), recursive_source_files_checked=len(checked),
        candidate_statement='Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most165.',
        arithmetic_assumption=False, compactness_assumption=False, simple_cusp_assumption=False,
        chosen_narrow_facet_cases=[
            dict(case='contains an ideal vertex', ambient_dimension_at_most=165,
                 primary=PRIMARY, independent_review=REVIEW),
            dict(case='only finite vertices, hence compact', ambient_dimension_at_most=143,
                 primary=FACET, independent_review='fresh_dimension225_hostile_20261015')],
        common_height_ratio_split='31/30', total_retained_cost_at_most=4,
        minimum_defective_factor_cost='3/2', maximum_defective_factors=2,
        maximum_each_defect='2/9', total_cost_and_defect_loss='40/9', avoid_node_deficit='2/9',
        near_pairing_lower='M*((n-1)-39)/9', far_pairing_lower='M*((n-1)-43)/9',
        first_excluded_dimension=166, combined_candidate_dimension_at_most=165,
        norm_two_narrow_threshold=14, frozen_baseline=BASE,
        primary_sources=data[BASE]['primary_sources'],
        scope='Reviewed dependency and complete case binder, not an independent proof of the inherited geometry. Internal research proof candidate; external validation, novelty and current-best status remain unestablished. All prior frozen artifacts are unchanged.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
