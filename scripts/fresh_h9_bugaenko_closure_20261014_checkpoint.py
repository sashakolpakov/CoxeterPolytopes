#!/usr/bin/env python3
"""Hash-bound scope ledger; it does not rerun or enlarge any search."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_closure_20261014'
BATCHES = ['search_0_25', 'search_25_250', 'search_500_650',
           'search_524_525', 'central', 'large_first_20_21', 'large_first_0_505']
CHECKS = ['', 'independent', 'central_closure', 'central_incidence_audit',
          'search_audit_0_25', 'search_audit_25_250', 'search_audit_500_650',
          'search_audit_524_525', 'search_audit_12_201',
          'large_first_audit_20_21', 'large_first_parallel_audit']


def path(suffix):
    return f'results/{TAG}{"_" + suffix if suffix else ""}.json'


def load_bound(p):
    data = json.loads((ROOT / p).read_text())
    for s, h in data['sources'].items():
        assert sha256((ROOT / s).read_bytes()).hexdigest() == h, s
    return data


def run():
    files = {path(s) for s in BATCHES + CHECKS}
    files.add('results/fresh_h9_bugaenko_surgery_20261014_incidence.json')
    data = {p: load_bound(p) for p in files}
    inc = data['results/fresh_h9_bugaenko_surgery_20261014_incidence.json']
    assert len(inc['cases']) == 2953
    tested = set()
    providers = {}
    batches = []
    for suffix in BATCHES:
        p = path(suffix)
        d = data[p]
        rows = d['cases']
        for r in rows:
            tested.add(r['case'])
            assert r['verdict'] in {'EXACT_COMPLETE_OBSTRUCTION', 'NODE_LIMIT_UNRESOLVED'}
            if r['verdict'] == 'EXACT_COMPLETE_OBSTRUCTION':
                assert r['nodes'] <= d['node_limit']
                providers[r['case']] = p
        batches.append(dict(source=p, cases=len(rows),
                            verdict_counts=dict(Counter(r['verdict'] for r in rows)),
                            node_limit=d['node_limit']))
    assert tested == set(range(650))
    unresolved = sorted(tested - providers.keys())
    assert unresolved == [22, 23, 318, 582, 629, 633, 636, 637]
    assert len(providers) == 642
    central = data[path('central_closure')]
    central_ids = {r['case'] for r in central['exclusion_records']}
    assert len(central_ids) == 128 and central_ids <= providers.keys()
    for suffix in CHECKS:
        assert data[path(suffix)]['status'].startswith('PASS'), suffix
    final_audit = data[path('large_first_parallel_audit')]
    assert final_audit['cases_replayed'] == 505 and final_audit['exact_exclusions'] == 502
    files.add(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(status='PASS_FROZEN_BOUNDED_CLOSURE_LEDGER',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in sorted(files)},
        catalogue_types=2953, tested_types=650, exact_exclusions=642,
        tested_but_unresolved=unresolved, untested_types=2303,
        untested_case_range=[650, 2953], total_unresolved=2311,
        first_carrier=dict(types=505, exact_exclusions=502, unresolved=[22, 23, 318]),
        central_carrier=dict(ordered_operation_orbits=15790,
            excluded_by_missing_order=15582, eligible_operations=208,
            distinct_eligible_types=128, exact_exclusions=128,
            newly_introduced_types=49, overlapping_first_carrier_types=79),
        batches=batches,
        exact_exclusion_sources={str(i): providers[i] for i in sorted(providers)},
        scope='Compact H9 incidence exclusions only. All finite integer labels and arbitrary '
              'dotted distances are allowed by the necessary gates. No node-limited case is '
              'excluded. No new Coxeter metric or finite-volume polytope was obtained. '
              'This ledger checks saved results and source hashes, not an independent search.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    target = ROOT / f'results/{TAG}_checkpoint.json'
    if args.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['exact_exclusions'], result['tested_but_unresolved'])
