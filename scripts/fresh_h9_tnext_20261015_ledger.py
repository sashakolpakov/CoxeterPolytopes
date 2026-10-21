#!/usr/bin/env python3
"""Compose only independently completed per-case exclusions; no search here."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_tnext_20261015_ledger'
GEOMETRY = 'results/fresh_h9_tree_cube_wedges_20261014.json'
PREVIOUS_CHECKPOINT = 'results/narrow_facet_bound_and_star_closure_2026_10_15.md'
PREVIOUS_CLOSED = (0, 1, 2, 3, 4, 5, 6, 8, 17, 22, 25, 26)


def run(cases):
    sources = {}
    seen = set()
    references = 0

    def digest(p):
        return sha256((ROOT / p).read_bytes()).hexdigest()

    def read(p):
        nonlocal references
        d = json.loads((ROOT / p).read_text())
        if p not in seen:
            seen.add(p)
            for q, expected in d.get('sources', {}).items():
                references += 1
                assert digest(q) == expected, q
                if q.endswith('.json'):
                    read(q)
        return d

    geo = read(GEOMETRY)
    assert len(geo['representatives']) == 33
    assert len(cases) == len(set(cases)) and cases
    assert not set(cases) & set(PREVIOUS_CLOSED)
    rows = []
    paths = [GEOMETRY, PREVIOUS_CHECKPOINT, f'scripts/{TAG}.py']
    for case in cases:
        initial_gpath = f'results/fresh_h9_tnext_20261015_gates_case_{case:02d}.json'
        initial_gate = read(initial_gpath)
        gpath = initial_gpath
        ppath = f'results/fresh_h9_tnext_20261015_search_case_{case:02d}_batch_001.json'
        ipath = f'results/fresh_h9_tnext_20261015_independent_case_{case:02d}.json'
        if case == 31:
            gpath = 'results/fresh_h9_tnext_20261015_case31_gates.json'
            ppath = 'results/fresh_h9_tnext_20261015_resume_case_31_batch_002.json'
            ipath = 'results/fresh_h9_tnext_20261015_case31_independent.json'
        gate, primary, independent = map(read, (gpath, ppath, ipath))
        assert all(d['case'] == case for d in (gate, primary, independent))
        assert primary['status'] == 'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
        assert primary['solution'] is None and primary['pending_frames'] == []
        assert independent['status'] == 'PASS_INDEPENDENT_COMPLETE_TNEXT_EXCLUSION'
        assert independent['complete'] and independent['cutoff'] is None
        assert independent['survivor'] is None and independent['unexplored_frontier'] is None
        assert independent['sources'][ppath] == digest(ppath)
        assert independent['sources'][gpath] == digest(gpath)
        assert independent['analytic_conditions'] == gate['total_predicates']
        rep = geo['representatives'][case]
        assert rep['id'] == case and rep['tree'] == 'T'
        rows.append(dict(case=case, vertices=initial_gate['vertices'],
                         inherited_predicates=initial_gate['inherited_predicates'],
                         new_predicates=gate['total_predicates']-initial_gate['inherited_predicates'],
                         total_predicates=gate['total_predicates'],
                         three_large_placements=[p['chosen_blocks'] for p in initial_gate['three_large_transfers']],
                         uncovered_three_large_triples=initial_gate['uncovered_three_large_block_triples'],
                         primary_nodes=primary['cumulative_additional_nodes'],
                         final_batch_nodes=primary['new_nodes'], independent_nodes=independent['nodes'],
                         finite_masks=independent['finite_masks_checked'],
                         recipe_multiplicity=rep['recipe_multiplicity']))
        paths.extend((initial_gpath, gpath, ppath, ipath))
    remaining = sorted(set(range(33)) - set(PREVIOUS_CLOSED) - set(cases))
    assert all(geo['representatives'][i]['tree'] == 'T' for i in remaining)
    new_recipes = sum(r['recipe_multiplicity'] for r in rows)
    closed = sorted(set(PREVIOUS_CLOSED) | set(cases))
    return dict(status='PASS_COMPOSITION_OF_INDEPENDENT_T_CASE_EXCLUSIONS',
                sources={p: digest(p) for p in paths}, cases=rows,
                prior_closed_representatives=list(PREVIOUS_CLOSED), new_closed_representatives=cases,
                all_closed_representatives=closed, remaining_representatives=remaining,
                newly_excluded_recipes=new_recipes,
                closed_recipes=sum(geo['representatives'][i]['recipe_multiplicity'] for i in closed),
                remaining_recipes=sum(geo['representatives'][i]['recipe_multiplicity'] for i in remaining),
                total_recipes=geo['recipe_count'],
                checked_json_manifests=len(seen), checked_hash_references=references,
                scope='The prior closures are inherited from the named frozen checkpoint. Each new closure '
                      'requires an empty primary traversal and a separate no-cutoff empty-start exact traversal. '
                      'This excludes only the listed compact H9 incidence representatives and their already '
                      'verified recipe orbits. No other eighteen-facet or noncompact incidence is excluded.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--cases', type=int, nargs='+', required=True)
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    data = run(args.cases)
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['new_closed_representatives'], 'remaining', data['remaining_representatives'])
