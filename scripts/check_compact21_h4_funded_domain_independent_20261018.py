#!/usr/bin/env python3
"""Independent PSD-domain and actual-angle check of the H4 transfer.

Imports the old exact field helper, not the new H4 producer or its
combinatorial finite-diagram classifier. High-label completeness is prose.
"""
from collections import Counter
import fcntl
from hashlib import file_digest
from itertools import combinations, permutations, product
import json
from pathlib import Path
import signal
import fresh_compact29_small_faces_20261015 as P
from run_lean_batches_20261018 import resources

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_h4_a1_triangle_priority_20261018'
OLD = {(0, 1): 5, (1, 2): 3, (2, 3): 3}


def matrix(edges):
    result = [[P.K(2 if i == j else 0) for j in range(6)] for i in range(6)]
    for (i, j), label in edges.items():
        result[i][j] = result[j][i] = -P.COS[label]
    return result


def payment(nodes, edges, marks):
    adjacent = {i: set() for i in nodes}
    for (i, j), label in edges.items():
        if label > 2 and i in adjacent and j in adjacent:
            adjacent[i].add(j); adjacent[j].add(i)
    distance = {marks[0]: 0}
    queue = [marks[0]]
    for i in queue:
        for j in sorted(adjacent[i]):
            if j not in distance:
                distance[j] = distance[i]+1; queue.append(j)
    if marks[1] not in distance:
        return (0, 0, 0), queue
    d = distance[marks[1]]
    ordinary = int(d <= 2)
    h4 = len(queue) == 4 and any(label == 5 and i in distance and j in distance
                                for (i, j), label in edges.items())
    five_edge = edges.get(tuple(sorted(marks))) == 5
    extra = int(h4 and (d == 3 or five_edge))
    removed = int(h4 and five_edge)
    return (ordinary+extra, ordinary+extra-removed, removed), queue


def canonical_source(nodes, edges, marks):
    _, component = payment(nodes, edges, marks)
    size = len(component)
    assert size in (4, 5)
    heavy = any(label == 4 and i in component and j in component
                for (i, j), label in edges.items())
    kind = ('B' if heavy else 'A')+str(size)
    canonical = {(i, i+1): (4 if heavy and i == 0 else 3) for i in range(size-1)}
    images = set()
    for perm in permutations(range(size)):
        naming = dict(zip(sorted(component), perm))
        if all(edges.get(tuple(sorted((i, j))), 2) ==
               canonical.get(tuple(sorted((naming[i], naming[j]))), 2)
               for i, j in combinations(component, 2)):
            images.add(tuple(sorted(naming[i] for i in marks)))
    assert images
    return kind, images


def main():
    with (ROOT/f'results/{TAG}.json').open() as stream:
        saved = json.load(stream)
    for path, want in saved['sources'].items():
        with (ROOT/path).open('rb') as stream:
            assert file_digest(stream, 'sha256').hexdigest() == want, path
    inverse = P.inverse(P.sub(matrix(OLD), list(range(5))))
    expected_roots = {tuple(tuple(c) for c in r['contacts']): r
                      for r in saved['surviving_roots']}
    accepted = set(); dangerous = {}; counts = Counter(); pd_tests = 0
    for labels in product(range(2, 6), repeat=5):
        contacts = tuple((i, m) for i, m in enumerate(labels) if m > 2)
        edges = OLD | {(i, 5): label for i, label in contacts}
        g = matrix(edges)
        finite = True
        for omit in (0, 1):
            # Reverse order differs from the producer's graph catalogue.
            ids = [i for i in reversed(range(6)) if i != omit]
            pd_tests += 1
            if P.psd_status(P.sub(g, ids)) != 1:
                finite = False; break
        if not finite:
            continue
        h = sum((P.COS[mi]*inverse[i][j]*P.COS[mj]
                 for i, mi in contacts for j, mj in contacts), P.K(-2))
        if h.sign() <= 0:
            counts['nonpositive_height'] += 1; continue
        cycles = [P.COS[mi]*P.COS[mj]
                  for (i, mi), (j, mj) in combinations(contacts, 2) if i < 4 and j < 4]
        if any(v.v[1] or v.v[3] for v in cycles):
            counts['cycle_outside_Qsqrt5'] += 1; continue
        assert not h.v[1] and not h.v[3]
        assert P.K((h.v[0], 0, -h.v[2], 0)).sign() <= 0
        counts['necessary_arithmetic_triangle_root'] += 1
        assert contacts in expected_roots
        record = expected_roots[contacts]
        assert P.enc(h) == record['height_squared']
        accepted.add(contacts)
        totals = []
        for omit in (0, 1, 5):
            marks = [i for i in (0, 1, 5) if i != omit]
            values, _ = payment([i for i in range(6) if i != omit], edges, marks)
            saved_corner = next(c for c in record['corners'] if c['omitted'] == omit)
            assert values == (saved_corner['total'], saved_corner['after_deleting_five_edge_extra'],
                              saved_corner['five_edge_extra'])
            totals.append(values)
        if sum(v[1] for v in totals) < 2:
            assert sum(v[0] for v in totals) == 2
            assert sum(v[1] for v in totals) == sum(v[2] for v in totals) == 1
            assert totals[-1][2] == 1
            kind, images = canonical_source([1, 2, 3, 4, 5], edges, [1, 5])
            dangerous[contacts] = (kind, images)
    assert accepted == set(expected_roots) and len(accepted) == 34
    assert dict(counts) == saved['status_counts']
    expected = {
        ((3, 3),): ('A4', {(0, 3)}),
        ((3, 3), (4, 3)): ('A5', {(0, 3), (1, 4)}),
        ((3, 3), (4, 4)): ('B5', {(1, 4)}),
        ((3, 4),): ('B4', {(0, 3)}),
    }
    assert dangerous == expected
    assert {tuple(tuple(c) for c in saved['surviving_roots'][i]['contacts'])
            for i in saved['dangerous_triangle_indices']} == set(expected)
    print('PASS_INDEPENDENT_H4_PSD_AND_FUNDED_ANGLES', pd_tests, sum(counts.values()),
          len(accepted), len(dangerous), 'capacities=1,2,1,1', flush=True)


if __name__ == '__main__':
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        signal.alarm(59); resources(); main(); resources(start=False)
