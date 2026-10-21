#!/usr/bin/env python3
"""Compose the frozen dimension21 reduction with the reviewed F25 exclusion."""
import argparse
import json
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_single_field_20261016'
OLD = 'results/fresh_compact21_field_reduction_20261016.json'
PRIMARY = 'results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.json'
REVIEW = 'results/fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.json'


@lru_cache(None)
def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def run():
    old, primary, review = [json.loads((ROOT / p).read_text()) for p in (OLD, PRIMARY, REVIEW)]
    assert old['status'] == 'PASS_COMPACT_DIMENSION21_REDUCTION_TO_TWO_OPEN_FIELDS'
    assert old['remaining_fields'] == ['Q(sqrt5)', 'Q(sqrt2,sqrt5)']
    assert not old['dimension21_excluded']
    assert primary['status'] == 'PASS_COMPLETE_F25_LANNER_CATALOGUE_AND_COMPACT_BOUND20_CANDIDATE'
    assert review['status'] == 'PASS_INDEPENDENT_F25_COMPACT_BOUND20_REVIEW'
    assert primary['compact_dimension_upper_bound_candidate'] == review['compact_dimension_upper_bound_candidate'] == 20
    assert primary['triangles']['complete_triangles'] == review['triangles']['accepted'] == [[2, 5, 8], [4, 5, 5]]
    assert review['counting']['strict_lower'] == review['counting']['vertex_upper'] == 40
    branches = []
    for row in old['branches']:
        assert row['dimension'] == 21
        if row['field'] == 'Q(sqrt2,sqrt5)':
            assert row['status'] == 'OPEN'
            row = dict(field=row['field'], dimension=21,
                       status='EXCLUDED_RESEARCH_CANDIDATE', upper_bound=20,
                       evidence=[PRIMARY, REVIEW])
        branches.append(row)
    assert len(branches) == 13
    assert len({r['field'] for r in branches}) == 13
    remaining = [row['field'] for row in branches if row['status'] == 'OPEN']
    assert remaining == ['Q(sqrt5)']
    checked = set()
    references = 0

    def check(path):
        nonlocal references
        if not path.endswith('.json') or path in checked:
            return
        checked.add(path)
        doc = json.loads((ROOT / path).read_text())
        for key in ('sources', 'source_hashes'):
            for source, digest in doc.get(key, {}).items():
                assert isinstance(digest, str) and sha(source) == digest, (path, source)
                references += 1
                check(source)

    paths = [OLD, PRIMARY, REVIEW, f'scripts/{TAG}.py']
    for path in paths:
        check(path)
    return dict(status='PASS_COMPACT_DIMENSION21_REDUCTION_TO_SINGLE_OPEN_FIELD',
                sources={p: sha(p) for p in paths}, branches=branches,
                remaining_fields=remaining, dimension21_excluded=False,
                compact_arithmetic_upper_bound_unchanged=21,
                finite_volume_arithmetic_upper_bound_unchanged=21,
                new_field_bound20_candidate='Q(sqrt2,sqrt5)',
                recursive_manifest_count=len(checked), source_reference_count=references,
                scope='Composition of internally reviewed proof candidates. Published dimension>=14 field list and rational isotropy are inherited explicitly. No exclusion of compact dimension21 over Q(sqrt5), no compact bound20 for all fields, and no external validation or priority claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['remaining_fields'], result['recursive_manifest_count'])
