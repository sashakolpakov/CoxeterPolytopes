#!/usr/bin/env python3
"""Bind every old surgery type to an exact, independently checked exclusion."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_closure'
INC = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'
OLD = 'results/fresh_h9_bugaenko_closure_20261014_checkpoint.json'
GATE = 'results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
BLOCK_AUDIT = 'results/fresh_h9_bugaenko_remaining_20261014_block_independent.json'
SEARCH_AUDIT = 'results/fresh_h9_bugaenko_remaining_20261014_search_audit.json'
NORMALIZATION = 'results/fresh_h9_bugaenko_remaining_20261014_cut_normalization.json'


def read_bound(p):
    d = json.loads((ROOT / p).read_text())
    for s, h in d['sources'].items():
        assert sha256((ROOT / s).read_bytes()).hexdigest() == h, s
    return d


def run():
    files = {INC, OLD, GATE, BLOCK_AUDIT, SEARCH_AUDIT, NORMALIZATION}
    inc, old, gate, ba, sa, normalization = map(read_bound,
        [INC, OLD, GATE, BLOCK_AUDIT, SEARCH_AUDIT, NORMALIZATION])
    assert len(inc['cases']) == 2953
    assert ba['status'] == 'PASS_SEPARATE_FRACTION_TWO_BLOCK_AUDIT'
    assert ba['applications_checked'] == 2365 and len(ba['schur_proofs']) == 42
    assert sa['status'] == 'PASS_SEPARATE_ARITHMETIC_REPLAY_476_NEEDED_EXCLUSIONS'
    assert normalization['status'] == 'PASS_NESTED_TRUNCATION_VERTEX_IDENTITY'
    providers = {r['case']: dict(kind='TWO_BLOCK_INERTIA', source=GATE)
                 for r in gate['applications']}
    inherited = 0
    for key, p in old['exact_exclusion_sources'].items():
        i = int(key)
        if i in providers:
            continue
        d = read_bound(p)
        r = next(r for r in d['cases'] if r['case'] == i)
        assert r['verdict'] == 'EXACT_COMPLETE_OBSTRUCTION' and r['nodes'] <= d['node_limit']
        providers[i] = dict(kind='PRIOR_COMPLETE_SEARCH', source=p)
        files.add(p)
        inherited += 1
    assert inherited == 112
    replayed = {r['case'] for r in sa['records']}
    assert len(replayed) == 476
    for p in sa['sources']:
        if p in [INC, OLD, GATE] or not p.endswith('.json'):
            continue
        d = read_bound(p)
        files.add(p)
        for r in d.get('cases', []):
            i = r['case']
            if i not in replayed or i in providers or r['verdict'] != 'EXACT_COMPLETE_OBSTRUCTION':
                continue
            assert r['nodes'] <= d['node_limit']
            providers[i] = dict(kind='NEW_COMPLETE_SEARCH', source=p)
    assert set(providers) == set(range(2953))
    assert inc['counts']['ordered_recipe_orbits'] == 73513
    assert inc['counts']['passes_missing_order'] == 5986
    files.add(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(status='PASS_COMPLETE_BUGAENKO_WEDGE_TWO_TRUNCATION_EXCLUSION',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in sorted(files)},
        dimension=9, facets=14, safe_recipe_orbits=73513,
        recipes_excluded_by_large_missing_faces=73513-5986,
        eligible_recipes=5986, distinct_eligible_types=2953,
        structural_exclusions=2365, inherited_search_exclusions=112,
        new_search_exclusions=476, unresolved_types=0,
        type_certificates={str(i): providers[i] for i in sorted(providers)},
        scope='No compact hyperbolic Coxeter polytope arises by wedging any facet of '
              'the specified11-facet Bugaenko H8 incidence and then truncating any two '
              'proper faces. The original exhaustive cuts of codimensions2..5 cover '
              'all admissible final nerves; the nested-cut identity extends to initial '
              'cuts of codimension6..9. All finite integer labels and arbitrary dotted '
              'distances are included. No arithmetic-field restriction or numerical '
              'metric-failure claim is used. Other operation orders and starting '
              'incidences are not excluded.')


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
    print(result['status'], result['distinct_eligible_types'])
