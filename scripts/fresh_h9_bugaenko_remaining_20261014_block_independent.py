#!/usr/bin/env python3
"""Separate Fraction/squarefree-radical replay of the structural block gate.

No producer or SymPy arithmetic is imported. The classical Lanner tables are
shared. Finite diagrams use the older separately written tree classifier.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path
import fresh_h22_multiquadratic_20261013_review as R
from fresh_h9_bugaenko_surgery_20261014_rank_independent import spherical
from independent_h9_14_record30_a1_labels import lanner_codes

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_block_independent'
INPUT = 'results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
ZERO = {}
ONE = R.C(1)


def neg(a):
    return R.scale(a, -1)


def sub(a, b):
    return R.add(a, neg(b))


def inverse(a):
    assert a
    numerator = ONE
    for mask in range(1, 8):
        changed = {d: c * (-1)**sum(bool(mask >> j & 1) and d % p == 0
                                     for j, p in enumerate([2, 3, 5])) for d, c in a.items()}
        numerator = R.mul(numerator, changed)
    denominator = R.mul(a, numerator)
    assert set(denominator) == {1} and denominator[1]
    return R.scale(numerator, 1 / denominator[1])


@lru_cache(None)
def canonical(n):
    pairs = list(combinations(range(n), 2))
    where = {p: i for i, p in enumerate(pairs)}
    maps = [[where[tuple(sorted((p[i], p[j])))] for i, j in pairs]
            for p in permutations(range(n))]
    rows = lanner_codes(n) if n != 3 else [c for c in product(range(2, 8), repeat=3)
                                          if sum(F(1, m) for m in c) < 1]
    return sorted({min(tuple(row[i] for i in move) for move in maps) for row in rows})


@lru_cache(None)
def finite(n, code):
    d = dict(zip(combinations(range(n), 2), code))
    return spherical(tuple(d.get(pair, 2) for pair in combinations(range(9), 2)))


def pinv(A):
    n = len(A)
    B = [row[:] + [ONE if i == j else ZERO for j in range(n)] for i, row in enumerate(A)]
    signs = []
    for k in range(n):
        pivot = B[k][k]
        signs.append(R.sign(pivot))
        B[k] = [R.mul(x, inverse(pivot)) for x in B[k]]
        for i in range(n):
            if i == k:
                continue
            a = B[i][k]
            B[i] = [sub(x, R.mul(a, y)) for x, y in zip(B[i], B[k])]
    assert signs == [1] * (n - 1) + [-1], signs
    return [row[n:] for row in B]


def trim(p):
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def padd(a, b):
    return trim([R.add(a[i] if i < len(a) else ZERO, b[i] if i < len(b) else ZERO)
                 for i in range(max(len(a), len(b)))])


def pmul(a, b):
    p = [ZERO for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i + j] = R.add(p[i + j], R.mul(x, y))
    return trim(p)


COS = {2: ZERO, 3: R.C(F(1, 2)), 4: {2: F(1, 2)},
       5: {1: F(1, 4), 5: F(1, 4)}, 6: {3: F(1, 2)}}


def matrix(n, d):
    G = [[[ONE] if i == j else [ZERO] for j in range(n)] for i in range(n)]
    for (i, j), m in d.items():
        G[i][j] = G[j][i] = [ZERO, neg(ONE)] if m == 7 else [neg(COS[m])]
    return G


def pdet(A):
    n = len(A)
    total = [ZERO]
    for p in permutations(range(n)):
        term = [R.C((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n)))]
        for i, j in enumerate(p):
            term = pmul(term, A[i][j])
        total = padd(total, term)
    return trim(total)


def schur_proof(G, n, AI):
    m = len(G) - n
    S = [[G[n + i][n + j][:] for j in range(m)] for i in range(m)]
    for i in range(m):
        for j in range(m):
            correction = [ZERO]
            for a in range(n):
                for b in range(n):
                    correction = padd(correction, pmul(pmul(G[n + i][a], [AI[a][b]]), G[b][n + j]))
            S[i][j] = padd(S[i][j], [neg(x) for x in correction])
    p = pdet(S)
    degree = len(p) - 1
    shifted = [ZERO for _ in p]
    for k, coefficient in enumerate(p):
        for j in range(k + 1):
            shifted[j] = R.add(shifted[j], R.scale(coefficient, F(comb(k, j)) * F(9, 10)**(k - j) / 10**j))
    bernstein = [R.psum(R.scale(shifted[j], F(comb(i, j), comb(degree, j)))
                         for j in range(i + 1)) for i in range(degree + 1)]
    assert all(R.sign(c) == -1 for c in bernstein), (p, bernstein)
    return dict(polynomial_coefficients=[R.enc(c) for c in p],
                negative_bernstein_coefficients=[R.enc(c) for c in bernstein])


def run():
    saved = json.loads((ROOT / INPUT).read_text())
    for p, h in saved['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    proofs = []
    type_counts = []
    for table in saved['tables']:
        n, m = table['orders']
        A_pairs = list(combinations(range(n), 2))
        B_pairs = list(combinations(range(n, n + m), 2))
        proper = [list(combinations([i for i in range(n + m) if i not in [a, b]], 2))
                  for a in range(n) for b in range(n, n + m)]
        def good(d):
            return all(finite(n + m - 2, tuple(d.get(p, 2) for p in pp)) for pp in proper)
        expected = {(tuple(r['first_code']), tuple(r['second_code'])): r
                    for r in table['internal_pairs_allowing_nonright']}
        observed = set()
        right = pairs = nonright = 0
        for ca in canonical(n):
            A = matrix(n, dict(zip(A_pairs, ca)))
            AI = pinv([[a[0] for a in row] for row in A])
            for cb in canonical(m):
                pairs += 1
                d = dict(zip(A_pairs, ca))
                d.update(zip(B_pairs, cb))
                active = []
                for pair in product(range(n), range(n, n + m)):
                    d[pair] = 3
                    if good(d):
                        active.append(pair)
                    del d[pair]
                if not active:
                    right += 1
                    continue
                observed.add((ca, cb))
                row = expected[ca, cb]
                assert active == list(map(tuple, row['potentially_nonright_cross_pairs']))
                values_seen = []
                for values in product(range(2, 8), repeat=len(active)):
                    if max(values) == 2:
                        continue
                    trial = dict(d)
                    trial.update(zip(active, values))
                    if not good(trial):
                        continue
                    assert max(values) <= 6 and sum(v == 7 for v in ca + cb) <= 1
                    values_seen.append(values)
                    proof = schur_proof(matrix(n + m, trial), n, AI)
                    proofs.append(dict(orders=[n, m], first_code=ca, second_code=cb,
                                       cross_values=values, **proof))
                    nonright += 1
                assert values_seen == [tuple(r['cross_values']) for r in row['completions']]
        assert observed == set(expected)
        assert pairs == table['internal_type_pairs'] and right == table['all_cross_pairs_forced_right']
        assert nonright == table['nonright_completions']
        type_counts.append(dict(orders=[n, m], internal_pairs=pairs, forced_right=right, nonright=nonright))
    assert len(proofs) == 42
    incpath = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'
    inc = json.loads((ROOT / incpath).read_text())
    for app in saved['applications']:
        row = inc['cases'][app['case']]
        a = sum(1 << i for i in app['first_missing_face'])
        b = sum(1 << i for i in app['second_missing_face'])
        assert not a & b and a in row['canonical_minimal_nonfaces'] and b in row['canonical_minimal_nonfaces']
        for i in app['first_missing_face']:
            for j in app['second_missing_face']:
                face = (a | b) ^ (1 << i) ^ (1 << j)
                assert any(face & v == face for v in row['canonical_vertices'])
    files = [INPUT, incpath, str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h22_multiquadratic_20261013_review.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py',
        'scripts/independent_h9_14_record30_a1_labels.py']
    return dict(status='PASS_SEPARATE_FRACTION_TWO_BLOCK_AUDIT',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        tables=type_counts, schur_proofs=proofs, applications_checked=len(saved['applications']),
        scope='Separate finite-tree classifier, Fraction squarefree radical arithmetic, '
              'Galois-conjugate inversion and permutation determinants. All cross-category '
              'completions and strict Bernstein signs are rebuilt; no producer or SymPy '
              'arithmetic is imported. The Lanner tables, catalogue and elementary '
              'cos(pi/7)>9/10 identity remain shared mathematical inputs.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    data = json.loads(json.dumps(run()))
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], len(data['schur_proofs']), data['applications_checked'])
