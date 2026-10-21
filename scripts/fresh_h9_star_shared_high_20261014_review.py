#!/usr/bin/env python3
"""Independent Fraction-basis affine Schur census for the shared high block."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product, combinations
from pathlib import Path
from fresh_h9_star_zero_20261014_own_three_review import (
    Z, O, COS, add, sub, neg, scale, mul, inverse, sign, total, dec_anp, determinant3)
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_shared_high_20261014'
S = [3, 4, 6, 7, 9, 10, 12, 13]; P = [2, 5, 8, 11]; C = [14, 15, 16, 17]
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    source = f'results/{TAG}.json'; data = json.loads((ROOT/source).read_text())
    pending, seen = [source], set()
    while pending:
        path = pending.pop()
        if path in seen: continue
        seen.add(path); doc = json.loads((ROOT/path).read_text())
        for p, h in doc.get('sources', {}).items():
            assert sha(p) == h, p
            if p.endswith('.json'): pending.append(p)
    geo = json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][0]['geometry']
    V = [v['incident'] for v in geo['vertices']]; sets = list(map(set, V))
    missing = {frozenset(m) for m in geo['minimal_nonfaces']}
    gates = json.loads((ROOT/'results/fresh_h9_star_zero_20261014_gates.json').read_text())
    own = {tuple(sorted((leaf, a))): 3 for i in range(3) for leaf in (P[i], C[i]) for a in S[2*i:2*i+2]}
    own.update({(11, 12): 2, (11, 13): 3, (12, 17): 2, (13, 17): 3})
    assert own == {(a, b): m for a, b, m in data['fixed_own_contacts']}
    for i in range(4):
        for leaf in (P[i], C[i]): assert frozenset(S[2*i:2*i+2]+[leaf]) in missing
    assert F(1, 2)+F(1, 3)+F(1, 6) == 1
    def face(row, nodes): assert set(row['support']) == set(nodes) <= sets[row['actual_vertex']]
    expected_cross = {tuple(sorted((a, b))) for i, j in combinations(range(4), 2)
                      for a in S[2*i:2*i+2] for b in S[2*j:2*j+2]}
    assert {tuple(r['pair']) for r in data['cross_right_guards']} == expected_cross
    for row in data['cross_right_guards']:
        a, b = row['pair']; i, j = S.index(a)//2, S.index(b)//2
        face(row, S[2*i:2*i+2]+S[2*j:2*j+2])
    assert 1-F(1, 2)-F(1, 4)-F(1, 2)+F(1, 4) == 0
    for row in data['first_three_high_reductions']:
        i = row['high_block']; labels = own | {tuple(S[2*j:2*j+2]): 6 if i == j else 7 if j == 3 else 4 for j in range(4)}
        gate = next(g for g in gates['families'] if g['id'] == row['reviewed_gate'])
        assert gate['origin'] == 'two_dense' and not gate['fixed']
        assert all(labels.get((a, b), 0) >= m for a, b, m in gate['bounds'])
    assert {r['high_block'] for r in data['first_three_high_reductions']} == {0, 1, 2}
    for row, (leaf, drop) in zip(data['projection_budgets'], ((11, 12), (17, 12), (14, 3))):
        assert row['leaf'] == leaf; face(row, (set(S)-{drop}) | {leaf})
    assert len(data['projection_budgets']) == 3
    assert F(1, 4)+2*F(1, 4)*2 > 1 and F(1, 4)+F(1, 2)*2 > 1
    face(data['high_block_third_right_guard'], (12, 13, 14))
    face(data['mutual_right_guard'], (11, 13, 17))
    assert 1-F(3, 4)-F(1, 4) == 0
    assert determinant3([[O, neg(COS[3]), neg(COS[3])],
                         [neg(COS[3]), O, neg(COS[3])], [neg(COS[3]), neg(COS[3]), O]]) == Z
    swap = data['literal_twin_swap']; assert swap == list(range(12))+[13, 12]+list(range(14, 18))
    assert {frozenset(swap[x] for x in v) for v in V} == set(map(frozenset, V))
    assert {frozenset(swap[x] for x in m) for m in missing} == missing
    base = own | {tuple(S[2*i:2*i+2]): 4 for i in range(4)}
    def vertex_failure(row, labels):
        nodes = V[row['actual_vertex']]
        M = [[O if a == b else neg(COS[labels.get(tuple(sorted((a, b))), 2)]) for b in nodes] for a in nodes]
        pivots = []
        for k in range(9):
            pivot = M[k][k]; pivots.append(pivot)
            if sign(pivot) <= 0: break
            for a in range(k+1, 9):
                for b in range(k+1, 9): M[a][b] = sub(M[a][b], mul(mul(M[a][k], M[k][b]), inverse(pivot)))
        assert sign(pivots[-1]) <= 0 and k+1 == row['prefix_order']
        assert pivots == [dec_anp(x) for x in row['pivots']]
    assert {tuple(r['pair']) for r in data['third_upper_guards']} == {(11, 14), (14, 17)}
    for row in data['third_upper_guards']:
        assert set(row['pair']) <= sets[row['actual_vertex']] and row['test_label'] == 4
        vertex_failure(row, base | {tuple(row['pair']): 4})
    observed, positive, rejected, minima = set(), 0, 0, []
    cache = {}
    for row in data['cases']:
        blocks = tuple(row['block_labels']); ia, ib, ic = row['foreign_contacts']; ma, mb = row['third_labels']
        key = (blocks, ia, ib, ic, ma, mb); assert key not in observed; observed.add(key)
        if blocks not in cache:
            J = [[Z]*6 for _ in range(6)]
            for i, m in enumerate(blocks):
                c = COS[m]; d = inverse(sub(O, mul(c, c))); assert sign(sub(d, scale(O, 2))) >= 0
                for a in range(2):
                    for b in range(2): J[2*i+a][2*i+b] = d if a == b else mul(d, c)
            cache[blocks] = J
        J = cache[blocks]
        def vector(extra, third=False):
            v = [Z]*6
            if third: v[0] = v[1] = neg(COS[3])
            if extra is not None: v[extra] = neg(COS[3])
            return v
        vectors = [vector(ia), vector(ib), vector(ic, True)]
        def inner(v, w): return total(mul(v[i], total(mul(J[i][j], w[j]) for j in range(6))) for i in range(6))
        raw = [[O, Z, neg(COS[ma])], [Z, O, neg(COS[mb])], [neg(COS[ma]), neg(COS[mb]), O]]
        M0 = [[sub(raw[i][j], inner(vectors[i], vectors[j])) for j in range(3)] for i in range(3)]
        def at(D): return [[sub(M0[i][j], scale(O, F(D, 4))) if i < 2 and j < 2 else M0[i][j] for j in range(3)] for i in range(3)]
        B = determinant3(M0); bound = determinant3(at(4)); slope = scale(sub(bound, B), F(1, 4))
        assert determinant3(at(5)) == add(B, scale(slope, 5))
        assert (slope, B, bound) == tuple(dec_anp(row[k]) for k in ('slope', 'intercept', 'boundary_D4'))
        assert sign(slope) > 0; minima.append(slope)
        if sign(bound) > 0:
            assert row['status'] == 'STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR'; positive += 1
        else:
            assert row['status'] == 'NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR'
            labels = base | {tuple(S[2*i:2*i+2]): blocks[i] for i in range(3)} | {(11, 14): ma, (14, 17): mb}
            for leaf, extra in zip((11, 17, 14), (ia, ib, ic)):
                if extra is not None: labels[tuple(sorted((leaf, S[extra])))] = 3
            vertex_failure(row['vertex_failure'], labels); rejected += 1
    expected = set(product(product((4, 5), repeat=3), [None]+list(range(6)), [None]+list(range(6)), [None, 2, 3, 4, 5], (2, 3), (2, 3)))
    assert observed == expected and positive == 7392 and rejected == 448
    minimum = minima[0]
    for a in minima[1:]:
        if sign(sub(a, minimum)) < 0: minimum = a
    assert minimum == dec_anp(data['minimum_slope'])
    paths = [source, f'scripts/{TAG}_review.py', 'scripts/fresh_h9_star_zero_20261014_own_three_review.py']
    return dict(status='PASS_INDEPENDENT_FRACTION_SHARED_HIGH_AFFINE_REVIEW',
                sources={p: sha(p) for p in paths}, complete_cases=7840,
                positive_affine_cases=positive, actual_vertex_failures=rejected,
                scope='No shared-high producer or symbolic algebra imports. Separate Fraction-basis Schur matrices and permutation determinants reconstruct all7840 affine cases atD=0,4,5, all448 vertex pivots and new domain guards. Rank-one dependence makes the determinant affine for every realD; D>4 covers all integer labels>=7. Prior two-dense gate proofs are hash-bound and their placements checked. Whole incidence remains open.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}_review.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['complete_cases'])
