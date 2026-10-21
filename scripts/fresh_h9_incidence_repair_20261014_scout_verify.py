#!/usr/bin/env python3
"""Replay saved incidence proposals and rank witnesses, without a SAT solver.

The exact spherical predicate is reused from the separate stdlib family
checker. Numerical proposal selection and its vertex rejections are NOT
promoted to an exact exhaustion theorem by this replay.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import argparse
import hashlib
import json
from fresh_h9_incidence_repair_20261014_verify import spherical

ROOT = Path(__file__).resolve().parents[1]
STEM = 'fresh_h9_incidence_repair_20261014'
SOURCE = ROOT / 'results/fresh_h9_incidence_surgery_20261013.json'
PATTERNS = ROOT / 'results/fresh_h9_incidence_surgery_20261013_dotted.json'
RELEASED = ROOT / f'results/{STEM}_released_sat.json'
METRIC = ROOT / f'results/{STEM}_metric.json'
OUT = ROOT / f'results/{STEM}_scout_verify.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    d = json.loads(SOURCE.read_text())
    row = next(r for r in d['cases'] if r['id'] == 'A16_14_1_C01')
    G = [[F(x) for x in r] for r in row['source_positive_gram']]
    old = json.loads(PATTERNS.read_text())
    assert old['source_sha256'] == sha(SOURCE)
    pats = next(r for r in old['cases'] if r['id'] == row['id'])['fixed_principal_witnesses']
    return row, G, pats


def labels(case, row, G, pats):
    pr = pats[case['pattern_index']]
    assert pr['support_index'] == case['support_index']
    assert pr['new_dotted_mask'] == case['new_dotted_mask']
    assert case['released_old_disjointness'] and not case['fixed_free_pairs']
    L = [[1 if i == j else {F(0): 2, F(1, 4): 3, F(1, 2): 4, F(1): 0}
          [G[i][j]**2 / (G[i][i]*G[j][j])] for j in range(15)] for i in range(15)]
    # Every selected weakening and every old parallel incidence is released.
    # Its subsequent value is right or one of the model's explicit dotted pairs.
    for i, j in row['parallel_pairs'] + row['finite_edge_cover']['minimum_supports'][pr['support_index']]:
        L[i][j] = L[j][i] = 2
    for i, j, m in case['added_finite_contacts']:
        assert L[i][j] == 2 and m in (3, 4)
        L[i][j] = L[j][i] = m
    assert all(L[i][j] in (2, 3, 4) for i, j in combinations(range(15), 2))
    return L


def complex_levels(L, D):
    levels = {0: [()]}
    for k in range(1, 11):
        levels[k] = []
        for B in levels[k-1]:
            for j in range(B[-1]+1 if B else 0, 15):
                if any((i, j) in D for i in B):
                    continue
                C = B+(j,)
                if spherical(L, C):
                    levels[k].append(C)
    return levels


def determinant_mod(L, B):
    p = 241
    roots = {d: next(x for x in range(p) if x*x % p == d) for d in (2, 3, 5)}
    assert all(x*x % p == d for d, x in roots.items())
    v = {1: 2, 2: 0, 3: -1, 4: -roots[2],
         5: -(1+roots[5])*pow(2, -1, p), 6: -roots[3]}
    A = [[v[L[i][j]] % p for j in B] for i in B]
    det = 1
    for k in range(len(B)):
        pivot = next((j for j in range(k, len(B)) if A[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            A[pivot], A[k] = A[k], A[pivot]
            det = -det
        q = A[k][k]
        det = det*q % p
        inv = pow(q, -1, p)
        for i in range(k+1, len(B)):
            c = A[i][k]*inv % p
            for j in range(k+1, len(B)):
                A[i][j] = (A[i][j]-c*A[k][j]) % p
    return det % p


def replay():
    row, G, pats = inputs()
    released = json.loads(RELEASED.read_text())
    metric = json.loads(METRIC.read_text())
    for path, expected in metric['sources'].items():
        assert sha(ROOT/path) == expected
    model_checks = []
    for ci, case in enumerate(released['cases']):
        L = labels(case, row, G, pats)
        for mi, model in enumerate(case['models']):
            D = set(map(tuple, model['extra_dotted_pairs']))
            assert len(D) == 3
            assert all(0 <= i < j < 15 and L[i][j] == 2 for i, j in D)
            levels = complex_levels(L, D)
            V = set(levels[9])
            assert len(V) == 216 and not levels[10]
            assert V == set(map(tuple, model['vertices']))
            faces = {k: set(B for W in V for B in combinations(W, k)) for k in range(1, 9)}
            for k in range(1, 9):
                assert set(levels[k]) == faces[k]
            ridge_counts = Counter(B for W in V for B in combinations(W, 8))
            assert set(ridge_counts.values()) == {2}
            model_checks.append({'case': ci, 'model': mi, 'spherical9': len(V),
                                 'spherical8': len(levels[8]), 'all_positive_faces_extend': True,
                                 'every_eight_face_has_two_vertices': True})
    assert len(model_checks) == 16
    rank_checks = []
    numerical = []
    for ci, case in enumerate(metric['cases']):
        assert case['base_case'] in released['cases']
        assert case['base_model'] in case['base_case']['models']
        L = labels(case['base_case'], row, G, pats)
        D = set(map(tuple, case['fixed_dotted_pairs']))
        assert D == set(map(tuple, case['base_model']['extra_dotted_pairs']))
        options = {(i, j): ms for i, j, ms in case['allowed_new_contacts']}
        for rejection in case['modular_rejections']:
            B = rejection['principal_nodes']
            assert len(B) == len(set(B)) == 11 and all(0 <= i < 15 for i in B)
            assert not any(i in B and j in B for i, j in D)
            Q = [r[:] for r in L]
            chosen = {(i, j): m for i, j, m in rejection['choices']}
            assert set(chosen) == {p for p in options if p[0] in B and p[1] in B}
            for (i, j), m in chosen.items():
                assert m in options[i, j] and m in range(2, 7)
                Q[i][j] = Q[j][i] = m
            det = determinant_mod(Q, B)
            assert det and det == rejection['determinant_mod241']
            rank_checks.append({'case': ci, 'round': rejection['round'],
                                'nodes': B, 'determinant_mod241': det})
        numerical.append(len(case['numerical_vertex_rejections']))
        assert not case['distance_trials']
    assert len(rank_checks) == 55 and sum(numerical) == 145
    paths = [SOURCE, PATTERNS, RELEASED, METRIC, Path(__file__).resolve(),
             ROOT/f'scripts/{STEM}_verify.py', ROOT/f'scripts/{STEM}.py',
             ROOT/f'scripts/{STEM}_metric.py']
    return {'status': 'PASS_SAVED_COMBINATORIAL_PROPOSALS_AND_EXACT_RANK_WITNESSES',
            'sources': {str(p.relative_to(ROOT)): sha(p) for p in paths},
            'models': model_checks, 'rank_witnesses': rank_checks,
            'numerical_vertex_rejections_not_promoted_to_exact': numerical,
            'distance_trials': 0,
            'scope': 'Sixteen saved incidence proposals satisfy the checked necessary compact face conditions; no metric realization or sphere classification follows. Fifty-five rank11 witnesses are exact and independent of dotted distances. Numerical contact-pool selection and 145 joint vertex failures remain numerical; solver UNSAT is not an exact unrestricted-family conclusion.'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    d = replay()
    if a.verify:
        assert d == json.loads(OUT.read_text())
        print('PASS saved proposal and modular-witness replay')
    else:
        OUT.write_text(json.dumps(d, indent=2)+'\n')
        print('WROTE', OUT)


if __name__ == '__main__':
    main()
