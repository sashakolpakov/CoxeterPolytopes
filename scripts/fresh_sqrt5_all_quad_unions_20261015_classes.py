#!/usr/bin/env python3
"""Exact marked-isomorphism quotient of the completed all-quad retained rows.

Uses all eight permutations of four black labels preserving the opposite pairs,
then every white permutation within canonical adjacency-signature groups.
Every member is accompanied by its literal complete-graph bijection.
No principal arithmetic or DFS is repeated here.
"""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_sqrt5_all_quad_unions_20261015'
TAG = BASE + '_classes'
OPPOSITES = {frozenset((0, 2)), frozenset((1, 3))}
BLACK_MAPS = tuple(p for p in permutations(range(4))
                  if {frozenset((p[0], p[2])), frozenset((p[1], p[3]))} == OPPOSITES)
assert len(BLACK_MAPS) == 8


def canonical(row):
    n = row['nodes']
    g = [[2] * n for _ in range(n)]
    for (i, j), m in zip(combinations(range(n), 2), row['labels']):
        g[i][j] = g[j][i] = m
    best = None
    best_order = None
    for p in BLACK_MAPS:
        groups = defaultdict(list)
        for w in range(4, n):
            signature = (tuple(g[w][b] for b in p),
                         tuple(sorted(g[w][v] for v in range(4, n) if v != w)))
            groups[signature].append(w)
        blocks = [groups[k] for k in sorted(groups)]
        for chosen in product(*(permutations(b) for b in blocks)):
            order = p + tuple(v for group in chosen for v in group)
            word = tuple(g[order[i]][order[j]] for i, j in combinations(range(n), 2))
            if best is None or word < best:
                best, best_order = word, order
    assert best is not None
    return (n, best), best_order


def check_map(row, target, mp):
    n = row['nodes']
    assert n == target['nodes'] and sorted(mp) == list(range(n))
    assert all((i < 4) == (mp[i] < 4) for i in range(n))
    assert {frozenset((mp[0], mp[2])), frozenset((mp[1], mp[3]))} == OPPOSITES
    tg = dict(zip(combinations(range(n), 2), target['labels']))
    for (i, j), m in zip(combinations(range(n), 2), row['labels']):
        assert tg[tuple(sorted((mp[i], mp[j])))] == m


def generate():
    data = json.loads((ROOT / 'results' / (BASE + '.json')).read_text())
    assert data['status'] == 'COMPLETE_ALL_QUAD_UNIONS'
    assert data['state']['complete'] and data['state']['left'] == 153
    rows = data['labelled_union_supports']
    classes, by_key, orders = [], {}, {}
    checked_pairs = 0
    for index, row in enumerate(rows):
        key, order = canonical(row)
        if key not in by_key:
            by_key[key] = len(classes)
            orders[key] = order
            classes.append(dict(representative=index, members=[],
                                canonical_labels=list(key[1])))
        rec = classes[by_key[key]]
        mp = [None] * row['nodes']
        for source, target in zip(order, orders[key]):
            mp[source] = target
        check_map(row, rows[rec['representative']], mp)
        checked_pairs += row['nodes'] * (row['nodes'] - 1) // 2
        rec['members'].append(dict(index=index, to_representative=mp))
    assert sum(len(x['members']) for x in classes) == len(rows)
    counts = Counter('dotted' if 0 in rows[x['representative']]['labels'] else 'finite'
                     for x in classes)
    node_counts = Counter(rows[x['representative']]['nodes'] for x in classes)
    paths = [f'scripts/{TAG}.py', f'scripts/{BASE}.py', f'results/{BASE}.json',
             f'results/{BASE}_review.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015_incidence_review.json']
    sources = {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    return dict(status='PASS_COMPLETE_RETAINED_MARKED_ISOMORPHISM_QUOTIENT',
                sources=sources, literal_rows=len(rows), marked_classes=len(classes),
                classes_by_parameter_type=dict(counts), classes_by_node_count=dict(node_counts),
                literal_pair_checks=checked_pairs, marked_isomorphism_classes=classes,
                scope='Exact quotient only of retained completed-union supports; preserves black colors and opposite pairs. No external-white completion, dotted-distance feasibility, or global dimension claim.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = generate()
    path = ROOT / 'results' / (TAG + '.json')
    # Round trip normalizes integer dictionary keys without changing arithmetic.
    out = json.loads(json.dumps(out))
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['literal_rows'], out['marked_classes'], out['classes_by_parameter_type'])


if __name__ == '__main__':
    main()
