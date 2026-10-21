#!/usr/bin/env python3
"""Rebuild the32 tested new nerves from vertices and re-arithmetize search."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
import fresh_h9_bugaenko_closure_20261014 as C
import fresh_h9_bugaenko_closure_20261014_large_first as P
import fresh_h9_bugaenko_closure_20261014_large_first_audit as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_newfacet_audit'
INC = 'results/fresh_h9_bugaenko_remaining_20261014_newfacet.json'
SEARCH = 'results/fresh_h9_bugaenko_remaining_20261014_newfacet_search.json'
RAW = 'data/known_examples/bugaenko_compact_h8.json'


def cut(vertices, face, new):
    face = frozenset(face)
    assert any(face <= v for v in vertices)
    return {v for v in vertices if not face <= v} | {
        frozenset((v - {i}) | {new}) for v in vertices if face <= v for i in face}


def run():
    inc = json.loads((ROOT / INC).read_text())
    search = json.loads((ROOT / SEARCH).read_text())
    for d in [inc, search]:
        for p, h in d['sources'].items():
            assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    original = {frozenset(v) for v in json.loads((ROOT / RAW).read_text())['geometry']['vertices']}
    C.finite = P.finite = A.finite
    C.residue = A.determinant
    records = []
    for saved in search['cases']:
        row = inc['cases'][saved['case']]
        first = cut(original, row['first_face'], 11)
        wedged = set()
        for v in first:
            if 11 in v:
                wedged.add(v | {12})
            else:
                wedged.update([v | {11}, v | {12}])
        vertices = cut(wedged, row['second_face'], 13)
        masks = [sum(1 << i for i in v) for v in vertices]
        M = clutter(masks, 14)
        assert max(map(int.bit_count, M)) <= 5
        order = row['canonical_order']
        assert sorted(order) == list(range(14))
        mapping = {old: new for new, old in enumerate(order)}
        move = lambda m: sum(1 << mapping[i] for i in range(14) if m >> i & 1)
        assert sorted(map(move, M)) == row['canonical_minimal_nonfaces']
        assert sorted(map(move, masks)) == row['canonical_vertices']
        assert all(len(v) == 9 for v in vertices)
        assert set(Counter(frozenset(v - {i}) for v in vertices for i in v).values()) == {2}
        assert row['old_catalogue_case'] is None and row['block_obstruction'] is None
        actual = json.loads(json.dumps(P.search(row, search['node_limit'])))
        assert actual == saved and actual['verdict'] == 'EXACT_COMPLETE_OBSTRUCTION'
        records.append(dict(case=row['case'], vertices=len(vertices), minimal_nonfaces=len(M),
                            evaluated_search_states=actual['nodes']))
    assert len(records) == 32
    files = [INC, SEARCH, RAW, str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py',
        'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
    return dict(status='PASS_SEPARATE32_CHANGED_NERVE_AND_ARITHMETIC_REPLAY',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        records=records,
        scope='Every tested changed recipe is independently rebuilt by literal vertex '
              'surgeries; its minimal nonfaces are recovered from the full downward '
              'complex, with no producer clutter formula or nauty reused. All32 exact '
              'label exclusions are replayed with separate finite/modular arithmetic, '
              'while reusing the DFS and canonicalizer. The other5504 incidence types '
              'and complete operation-orbit enumeration are not independently rerun here.')


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
    print(result['status'], len(result['records']))
