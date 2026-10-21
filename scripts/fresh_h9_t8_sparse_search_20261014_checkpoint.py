#!/usr/bin/env python3
"""Fast source-bound consistency check of the completed T8 search phase.

This checks terminal coverage/status and provenance.  It does not replace
the saved full DFS replay or the separately implemented theorem reviews.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_sparse_search_20261014_checkpoint'
PREFIX = 'results/fresh_h9_t8_sparse_search_20261014'
STAGES = [PREFIX + s for s in ('_case_08.json', '_next_case_08.json',
                               '_final_case_08.json', '_outside_case_08.json')]
PROOFS = [PREFIX + s for s in ('_all_blocks.json', '_two_high_blocks.json',
                               '_bridge_variants.json', '_twin.json')]
REVIEWS = ['results/fresh_h9_t8_all_blocks_20261014_review.json',
           PREFIX + '_two_high_blocks_review.json',
           'results/fresh_h9_t8_two_high_transfers_20261014_review.json']

def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

def run():
    seen = set()
    refs = {}
    def load(p):
        obj = json.loads((ROOT / p).read_text())
        if p in seen:
            return obj
        seen.add(p)
        for q, h in obj.get('sources', {}).items():
            assert sha(q) == h, q
            if q in refs:
                assert refs[q] == h
            refs[q] = h
            if q.endswith('.json'):
                load(q)
        return obj
    states = [load(p) for p in STAGES]
    proofs = [load(p) for p in PROOFS]
    reviews = [load(p) for p in REVIEWS]
    assert all('PASS' in r['status'] for r in proofs + reviews)
    assert [s['nodes'] for s in states] == [1942, 4232, 21331, 59423]
    assert all(s['case'] == 8 and s['settings']['node_limit'] == 100000 for s in states)
    assert all(s['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR' for s in states[:3])
    def satisfies(state, fixed, bounds):
        d = {(a, b): m for a, b, m in state['solution']['pair_labels']}
        return all(d[a, b] == m for a, b, m in fixed) and all(d[a, b] >= m for a, b, m in bounds)
    assert satisfies(states[0], proofs[0]['fixed_labels'], proofs[0]['finite_label_lower_bounds'])
    assert satisfies(states[1], [], proofs[1]['finite_label_lower_bounds'])
    assert satisfies(states[2], proofs[3]['fixed_labels'], proofs[3]['finite_label_lower_bounds'])
    end = states[-1]
    assert end['status'] == 'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
    assert end['solution'] is None and end['unexplored_frontier'] is None
    assert end['nodes'] < end['settings']['node_limit']
    assert len(end['triangle_supports']) == 10
    assert len(end['proved_family_conditions']) == 4
    assert len(end['bridge_clauses']) == 60
    assert end['prunes'] == {
        'spherical_vertex': 238748,
        'orthogonal_negative_blocks': 32665,
        'certified_all_blocks_sparse_partial': 3420,
        'certified_two_high_placement_1_partial': 3648,
        'certified_two_high_placement_0_partial': 7296,
        'certified_all_blocks_sparse_twin_partial': 3420,
    }
    paths = STAGES + PROOFS + REVIEWS + [f'scripts/{TAG}.py']
    return dict(status='PASS_COMPLETED_T8_NECESSARY_LABEL_EXCLUSION_CHECKPOINT',
        sources={p: sha(p) for p in paths},
        source_manifest_count=len(seen), checked_unique_source_count=len(refs),
        stages=[dict(source=p, status=s['status'], nodes=s['nodes'], prunes=s['prunes'])
                for p, s in zip(STAGES, states)],
        complete=True, node_limit_reached=False, surviving_assignments=0,
        scope='The literal 18-facet T-case8 incidence is excluded for compact H9 Coxeter '
              'realizations. This terminal checker validates saved search status, all prior '
              'survivor exclusions and recursive provenance, not an independent DFS. '
              'Full saved DFS replay and the separately implemented all-label proof reviews '
              'have distinct roles. No claim about other catalogue incidences or all compact H9.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['source_manifest_count'], result['checked_unique_source_count'])
