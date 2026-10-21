#!/usr/bin/env python3
"""Literal metric-obstruction transfers to the saved T31 frontier; no symmetry assumption."""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_tnext_20261015_case31_transfers'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
OLD = 'results/fresh_h9_tnext_20261015_gates_case_31.json'
STATE = 'results/fresh_h9_tnext_20261015_search_case_31_batch_001.json'
THEOREM = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json'
REVIEW = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    seen = set()
    def read(p):
        d = json.loads((ROOT / p).read_text())
        if p not in seen:
            seen.add(p)
            for q, h in d.get('sources', {}).items():
                assert sha(q) == h, q
                if q.endswith('.json'):
                    read(q)
        return d
    geo, old, state, theorem, review = map(read, (GEOM, OLD, STATE, THEOREM, REVIEW))
    assert old['case'] == state['case'] == 31
    assert old['total_predicates'] == 153
    assert state['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert theorem['status'] == 'PASS_T8_TWO_HIGH_BLOCKS_ALL_OTHER_LABELS_FREE_EXCLUDED'
    assert review['status'].startswith('PASS_') and review['sources'][THEOREM] == sha(THEOREM)
    g = geo['representatives'][31]['geometry']
    vertices = [set(v['incident']) for v in g['vertices']]
    missing = {frozenset(m) for m in g['minimal_nonfaces']}
    blocks = [b for b in g['blocks'] if len(b) == 2]
    assert blocks == [[1, 2], [4, 5], [7, 8], [12, 13]]
    outside = sorted(set(range(18)) - set(sum(blocks, [])))
    supports = [('positive_support', theorem['positive_support'])]
    supports += [('meeting_pair', r['pair']) for r in theorem['all45_core_meeting_pairs']]
    for key in ('cross_right_guards', 'anchor_high_right_budgets',
                'third_other_blocks_right_budgets', 'third_anchor_upper_stars'):
        supports += [(key, r['support']) for r in theorem[key]]
    assert len(supports) == 80
    def actual(s):
        return next((i for i, v in enumerate(vertices) if set(s) <= v), None)
    maps = []
    for high in combinations(range(4), 2):
        low = sorted(set(range(4)) - set(high))
        for own in low:
            order = [*high, next(i for i in low if i != own), own]
            for bridge in outside:
                if not all(frozenset(blocks[i] + [bridge]) in missing for i in high):
                    continue
                for anchor in outside:
                    if anchor == bridge or frozenset(blocks[own] + [anchor]) not in missing:
                        continue
                    mp = dict(zip(theorem['positive_support'], sum((blocks[i] for i in order), [])))
                    mp.update({14: bridge, 11: anchor})
                    assert len(set(mp.values())) == len(mp) == 10
                    checks = []
                    for kind, source in supports:
                        target = [mp[x] for x in source]
                        vi = actual(target)
                        if vi is None:
                            break
                        checks.append(dict(kind=kind, source=source, target=target, actual_vertex=vi))
                    else:
                        triples = [sorted(mp[x] for x in s) for s in theorem['missing_triples']]
                        assert all(frozenset(s) in missing for s in triples)
                        bounds = sorted(sorted((mp[a], mp[b])) + [m]
                                        for a, b, m in theorem['finite_label_lower_bounds'])
                        assert not theorem['fixed_labels']
                        maps.append(dict(mapping=sorted(map(list, mp.items())), high_blocks=list(high),
                                         anchor_block=own, bridge=bridge, anchor=anchor,
                                         supports=checks, missing_triples=triples, bounds=bounds))
    conditions = []
    keys = set()
    old_keys = {(tuple(map(tuple, f['fixed'])), tuple(sorted(map(tuple, f['bounds'])))) for f in old['families']}
    for i, row in enumerate(maps):
        key = ((), tuple(map(tuple, row['bounds'])))
        if key in keys or key in old_keys:
            continue
        keys.add(key)
        conditions.append(dict(id=f'two_high_metric_transfer_{len(conditions)}',
                               origin='reviewed_literal_two_high_negative_plane',
                               fixed=[], bounds=row['bounds'], witness_map=i))
    labels = {(a, b): m for a, b, m in state['solution']['pair_labels']}
    hits = [i for i, c in enumerate(conditions) if all(labels[a, b] >= m for a, b, m in c['bounds'])]
    assert hits
    paths = [GEOM, OLD, STATE, THEOREM, REVIEW, f'scripts/{TAG}.py']
    return dict(status='PASS_LITERAL_T31_TWO_HIGH_METRIC_TRANSFERS', case=31,
                sources={p: sha(p) for p in paths}, inherited_predicates=153,
                maps=maps, new_conditions=conditions, saved_survivor_hits=hits,
                total_predicates=153 + len(conditions),
                literal_support_occurrences=sum(len(r['supports']) for r in maps),
                scope='Each map separately checks all 80 source-proof face premises and three missing triples. '
                      'The inherited negative-plane proof excludes the full stated high-block domain with all '
                      'other labels and dotted distances free, in every Lorentz ambient dimension. '
                      'Only these maps and the saved survivor are handled; no full-nerve isomorphism or '
                      'whole T31-incidence exclusion is claimed. No continuation is run.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    data = run()
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], len(data['maps']), len(data['new_conditions']), data['saved_survivor_hits'])
