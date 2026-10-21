#!/usr/bin/env python3
"""Independent standard-library replay; no producer or numerical imports.

Mathematical implications checked by the finite data:
1. A positive kernel vector in an affine principal Gram block gives a nonzero
   future null vector, because its product with a strict interior point is
   negative. Every acute-compatible wall pairs nonpositively with that vector.
   Thus even additional such walls cannot make the chamber compact.
2. Pairwise orthogonal blocks of one or two unit-length-squared-2 normals have
   full rank unless a two-wall block is parallel (weight 2). A compact chamber
   cannot have that parallel block by implication 1. More than ten leaf normals
   in those blocks therefore contradict rank at most ten in H9.
"""
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / 'data/enumerations/fresh_h9_cusp_filling_obstructions_20261013.json'
OUTPUT = ROOT / 'results/fresh_h9_cusp_filling_obstructions_verify_20261013.json'


def positive_definite(matrix):
    a = [row[:] for row in matrix]
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot > 0
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / pivot


def check(d):
    path = ROOT / d['source']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == d['source_sha256']
    src = json.loads(path.read_text())['certified_types']
    grams = []
    for name in ('T1', 'T2'):
        row = next(r for r in src if r['id'] == d[name]['source_id'])
        g = [[-Q(x) / 2 for x in row] for row in row['facet_gram']]
        assert len(g) == 18 and all(len(row) == 18 for row in g)
        assert all(g[i][i] == 2 for i in range(18))
        assert all(g[i][j] == g[j][i] and g[i][j] in (0, -1, -2)
                   for i, j in combinations(range(18), 2))
        grams.append(g)
    g = grams[0]
    circuits = d['T1']['circuits']
    assert len({tuple(sorted(c['nodes'])) for c in circuits}) == len(circuits)
    for c in circuits:
        nodes, weights = c['nodes'], c['kernel']
        assert len(nodes) == len(weights) and len(set(nodes)) == len(nodes)
        assert all(x > 0 for x in weights)
        a = [[g[i][j] for j in nodes] for i in nodes]
        assert all(a[i][j] != -2 for i, j in combinations(range(len(a)), 2))
        assert all(sum(a[i][j] * weights[j] for j in range(len(a))) == 0
                   for i in range(len(a)))
        positive_definite([row[:-1] for row in a[:-1]])
        degrees = sorted(sum(x == -1 for x in row) for row in a)
        assert (c['type'], degrees) in [('affine_A3', [2, 2, 2, 2]),
                                       ('affine_D4', [1, 1, 1, 1, 4])]
    rows = d['T1']['five_deletion_cover']
    assert {tuple(r['deleted']) for r in rows} == set(combinations(range(18), 5))
    assert len(rows) == 8568
    for row in rows:
        assert not set(row['deleted']).intersection(circuits[row['circuit']]['nodes'])
    six = d['T1']['six_deletion_circuit_transversal']
    assert len(set(six)) == 6 and set(six) <= set(range(18))
    assert all(set(six).intersection(c['nodes']) for c in circuits)
    h, t = grams[1], d['T2']
    centers, leaves = set(t['centers']), set(t['leaves'])
    assert len(centers) == 2 and len(leaves) == 16
    assert not centers.intersection(leaves) and centers | leaves == set(range(18))
    assert all(sum(x == -1 for x in h[i]) == 8 for i in centers)
    pairs = {tuple(sorted(pair)) for pair in t['leaf_matching']}
    assert len(pairs) == 8 and sorted(x for pair in pairs for x in pair) == sorted(leaves)
    assert all(h[i][j] == (-2 if (i, j) in pairs else 0)
               for i, j in combinations(sorted(leaves), 2))
    counts = Counter(len(leaves - set(s)) for s in combinations(range(18), 5))
    assert dict(counts) == {int(k): v for k, v in t['five_deletion_leaf_count_histogram'].items()}
    assert min(counts) == 11 > t['ambient_rank_bound'] == 10
    # The polynomial determinant of every possible changed pair is 4-w^2.
    # On w >= 0 its only rank-drop point is w=2, the forbidden affine pair.
    pair = [[Q(2), Q(-2)], [Q(-2), Q(2)]]
    assert all(sum(row) == 0 for row in pair)
    positive_definite([[Q(2)]])
    return {'T1_affine_circuits': len(circuits),
            'T1_five_deletion_sets': len(rows),
            'T1_circuit_transversal_minimum': 6,
            'T2_minimum_retained_leaves': min(counts),
            'T2_ambient_rank': 10}


def main():
    data = json.loads(ARCHIVE.read_text())
    checks = check(data)
    mutations = {}
    bad = copy.deepcopy(data)
    bad['T1']['circuits'][0]['kernel'][0] += 1
    cases = [('wrong_affine_kernel', bad)]
    bad = copy.deepcopy(data)
    bad['T1']['five_deletion_cover'].pop()
    cases.append(('missing_deletion_case', bad))
    bad = copy.deepcopy(data)
    bad['T2']['leaf_matching'][0][1] = bad['T2']['leaf_matching'][1][1]
    cases.append(('wrong_leaf_matching', bad))
    for label, bad in cases:
        try:
            check(bad)
        except (AssertionError, IndexError, KeyError, ValueError):
            mutations[label] = 'REJECTED'
        else:
            raise AssertionError('Mutation was accepted: ' + label)
    result = {'status': 'PASS', 'arithmetic': 'stdlib Fraction only',
              'checks': checks, 'mutation_checks': mutations,
              'scope': 'Compact fillings only, under the declared fixed-angle '
                       'or fixed-orthogonality assumptions; no finite-volume exclusion.'}
    OUTPUT.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
