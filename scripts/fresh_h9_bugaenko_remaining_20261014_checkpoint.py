#!/usr/bin/env python3
"""Frozen theorem plus the explicitly unresolved changed-order frontier."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_checkpoint'
FILES = ['results/fresh_h9_bugaenko_remaining_20261014_' + name + '.json' for name in
         ['closure', 'block_gate_review', 'cut_normalization_review',
          'newfacet', 'newfacet_search', 'newfacet_audit']]


def run():
    data = {}
    for p in FILES:
        d = json.loads((ROOT / p).read_text())
        for q, h in d['sources'].items():
            assert sha256((ROOT / q).read_bytes()).hexdigest() == h, q
        data[p] = d
    old, review, cut_review, inc, search, audit = [data[p] for p in FILES]
    assert old['status'] == 'PASS_COMPLETE_BUGAENKO_WEDGE_TWO_TRUNCATION_EXCLUSION'
    assert old['distinct_eligible_types'] == 2953 and old['unresolved_types'] == 0
    eligible = set(inc['genuinely_new_unexcluded_types'])
    excluded = {r['case'] for r in search['cases'] if r['verdict'] == 'EXACT_COMPLETE_OBSTRUCTION'}
    assert len(eligible) == 3276 and len(excluded) == 32 and excluded <= eligible
    assert len(audit['records']) == 32 and {r['case'] for r in audit['records']} == excluded
    remaining = sorted(eligible - excluded)
    assert len(remaining) == 3244
    files = FILES + [str(Path(__file__).resolve().relative_to(ROOT))]
    return dict(status='PASS_CLOSED_OLD_FAMILY_CHANGED_ORDER_FRONTIER',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        old_family=dict(all_wedge_carriers=True, all_two_proper_face_truncations=True,
                        eligible_types=2953, exactly_excluded=2953, unresolved=0),
        changed_order=dict(operation='original H8 face cut; wedge its new facet; second face cut',
                           first_face_orbits=39, recipe_representatives=11729,
                           passing_recipes=6226, incidence_types=5536,
                           old_catalogue_duplicates=1040, additional_block_exclusions=1220,
                           pilot_exact_exclusions=32, metric_survivors=0,
                           unsearched_unexcluded=3244,
                           open_case_ids=remaining),
        scope='No new compact H9 polytope. The original operation family is completely '
              'excluded; the changed operation order remains open on the3244 explicitly '
              'listed incidence types. No numerical failure or node limit is promoted '
              'to an exact exclusion.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    result = run()
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['changed_order']['unsearched_unexcluded'])
