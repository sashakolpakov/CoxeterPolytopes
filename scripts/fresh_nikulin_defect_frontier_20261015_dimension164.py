#!/usr/bin/env python3
"""Bind the independently reviewed residual lemma and dimension164 split."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_defect_frontier_20261015_dimension164'
BASE = 'fresh_nikulin_height_charge_20261015_dimension165'
PRIMARY = 'fresh_nikulin_defect_frontier_20261015'
REVIEW = PRIMARY + '_review'
FACET = 'fresh_nikulin_narrow_20261015_compact_facet'


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def run():
    tags = [BASE, PRIMARY, REVIEW, FACET]
    data = {t: json.loads((ROOT / f'results/{t}.json').read_text()) for t in tags}
    assert data[BASE]['status'] == 'PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_165_CANDIDATE'
    assert data[PRIMARY]['status'] == 'PASS_EXACT_RESIDUAL_COST_ENVELOPE_CANDIDATE164_PENDING_INDEPENDENT_REVIEW'
    assert data[REVIEW]['status'] == 'PASS_INDEPENDENT_PHYSICAL_GRAM_KKT_SUBUNIT_ENVELOPE_AND_164_SCOPE_REVIEW'
    for path in (f'scripts/{PRIMARY}.py', f'results/{PRIMARY}.json'):
        assert data[REVIEW]['sources'][path] == digest(path)
    p = data[PRIMARY]
    assert Q(p['height_ratio_window']) == Q(21, 20)
    assert Q(p['total_loss_at_most']) == Q(38, 9)
    # Check the aggregate split independently of the symbolic producer.
    bounds = [4 + Q(2, 9), Q(28, 9) + 1,
              Q(30, 9) + Q(2, 3), Q(32, 9) + Q(2, 3)]
    assert bounds == [Q(c['loss_bound']) for c in p['aggregate_cases']]
    assert max(bounds) == Q(38, 9)
    assert 4 - 2 * Q(14, 9) < 1
    assert 4 - Q(14, 9) - Q(3, 2) < 1
    assert 3 * Q(3, 2) > 4
    # N=n-1; exactly one target loses 2/9 for the prescribed avoided node.
    near = (Q(2 * 164, 9) - Q(2, 9) - Q(38, 9)) / 2 - 2
    far = Q(21, 20) * Q(164 - 43, 9)
    t = p['threshold_checks']
    assert near == Q(t['near_lower']) == 14
    assert far == Q(t['far_lower']) == 14 + Q(7, 60)
    assert t['combined_candidate_bound'] == t['ideal_branch_bound'] == 164
    assert data[REVIEW]['aggregate']['combined_candidate'] == data[REVIEW]['aggregate']['ideal_branch_bound'] == 164
    assert data[REVIEW]['aggregate']['inherited_compact_branch'] == 143
    assert Q(data[REVIEW]['aggregate']['zero_or_one_defect_loss']) == Q(38, 9)
    assert t['compact_facet_branch_inherited_bound'] == data[FACET]['ambient_dimension_at_most'] == 143
    assert Q(1, 9) > 0 and Q(21, 20) / 9 == Q(7, 60) > 0
    assert (Q(2 * 163, 9) - Q(2, 9) - Q(38, 9)) / 2 - 2 < 14
    files = [f'scripts/{TAG}.py', f'results/{TAG}.md']
    for tag in tags:
        files += [f'scripts/{tag}.py', f'results/{tag}.json', f'results/{tag}.md']
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

    for tag in tags:
        recurse(f'results/{tag}.json')
    return dict(
        status='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_164_CANDIDATE',
        sources={path: digest(path) for path in files},
        recursive_manifests_checked=len(seen), recursive_source_files_checked=len(checked),
        candidate_statement='Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most164.',
        arithmetic_assumption=False, compactness_assumption=False, simple_cusp_assumption=False,
        chosen_narrow_facet_cases=[
            dict(case='contains an ideal vertex', ambient_dimension_at_most=164,
                 primary=PRIMARY, independent_review=REVIEW),
            dict(case='only finite vertices, hence compact', ambient_dimension_at_most=143,
                 primary=FACET, independent_review='fresh_dimension225_hostile_20261015')],
        common_height_ratio_split='21/20', total_retained_cost_at_most=4,
        strict_subunit_residual_loss_at_most='2/3', total_energy_loss_at_most='38/9',
        avoid_node_deficit='2/9', near_pairing_lower='M*((n-1)-38)/9',
        far_pairing_lower='M*((n-1)-43)/9', first_excluded_dimension=165,
        combined_candidate_dimension_at_most=164, norm_two_narrow_threshold=14,
        frozen_baseline=BASE, primary_sources=data[BASE]['primary_sources'],
        scope='Reviewed local-proof dependency and complete case binder. Inherited geometry remains a mathematical dependency. Internal research proof candidate; external validation, novelty and current-best status remain unestablished. No bound163 follows from this estimate.')


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
