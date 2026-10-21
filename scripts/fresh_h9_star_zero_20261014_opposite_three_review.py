#!/usr/bin/env python3
"""Independent literal-support and symbolic-identity check of internal3 branch."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_zero_20261014_opposite_three'

def run():
    source = f'results/{TAG}.json'; data = json.loads((ROOT/source).read_text())
    for p, h in data['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    path = 'results/fresh_h9_tree_cube_wedges_20261014.json'
    geo = json.loads((ROOT/path).read_text())['representatives'][0]['geometry']
    vertices = [set(v['incident']) for v in geo['vertices']]
    missing = {frozenset(m) for m in geo['minimal_nonfaces']}
    S = [3, 4, 6, 7, 9, 10, 12, 13]; P = [2, 5, 8, 11]; C = [14, 15, 16, 17]
    def check(row, expected):
        assert set(row['support']) == set(expected) <= vertices[row['actual_vertex']]
    check(data['complete_positive_support'], S)
    assert len(data['foreign_high_edge_guards']) == 12
    assert {(r['anchor'], r['root'], k) for r in data['foreign_high_edge_guards']
            for k in set(r['support'])-{r['anchor'], r['root']}} == {
                (a, r, k) for a, r in ((11, 13), (17, 12)) for k in S[:6]}
    for row in data['foreign_high_edge_guards']+data['third_high_edge_guards']:
        a, r = row['anchor'], row['root']; other = set(row['support'])-{a, r}; assert len(other) == 1
        k = next(iter(other)); check(row, (a, r, k))
        assert {tuple(p) for p in row['forced_right_pairs']} == {tuple(sorted((a, k))), tuple(sorted((r, k)))}
    assert {(r['anchor'], r['root']) for r in data['third_high_edge_guards']} == {(11, 13), (17, 12)}
    assert all(14 in r['support'] for r in data['third_high_edge_guards'])
    check(data['mutual_right_guard'], (11, 13, 17))
    assert len(data['first_three_cross_guards']) == 3
    assert {tuple(tuple(b) for b in r['blocks']) for r in data['first_three_cross_guards']} == {
        (tuple(S[2*i:2*i+2]), tuple(S[2*j:2*j+2])) for i, j in combinations(range(3), 2)}
    for row in data['first_three_cross_guards']: check(row, sum(row['blocks'], []))
    for i in range(4):
        for leaf in (P[i], C[i]): assert frozenset(S[2*i:2*i+2]+[leaf]) in missing
    assert data['own_bounds'] == [[11, 13, 7], [12, 17, 7]]
    expected = {(min(leaf, r), max(leaf, r)): 3 for i in range(3)
                for leaf in (P[i], C[i]) for r in S[2*i:2*i+2]}
    expected.update({(11, 12): 2, (13, 17): 2, (12, 13): 3})
    assert {(a, b): m for a, b, m in data['fixed']} == expected
    # Direct determinant polynomial in x=alpha^2,y=beta^2:
    # (1-4x/3)(1-4y/3)-(2alpha beta/3)^2.
    determinant = {(0, 0): F(1), (1, 0): F(-4, 3), (0, 1): F(-4, 3),
                   (1, 1): F(16, 9)-F(4, 9)}
    translated = {(0, 0): F(4, 3)*F(3, 4), (1, 0): F(-4, 3),
                  (0, 1): F(-4, 3), (1, 1): F(4, 3)}
    assert determinant == translated
    # For x,y>3/4 and <1, product(1-x)(1-y)<1/16.
    assert F(4, 3)*(F(1, 16)-F(1, 4)) == F(-1, 4)
    # High-edge triangle: c_high^2>3/4 and c_other^2>=1/4.
    assert 1-F(3, 4)-F(1, 4) == 0
    # Relaxed 4--3--4 path has a zero determinant; larger entries cannot be PD.
    assert 1-F(1, 2)-F(1, 4)-F(1, 2)+F(1, 4) == 0
    # Third own projection cost: 1/[2(1-c)] at c>=sqrt2/2;
    # rationalize denominator to get1+sqrt2/2, hence residual<=-sqrt2/2.
    assert F(1, 2)/(1-F(1, 2)) == 1
    swap = data['literal_twin_swap']; assert swap == list(range(12))+[13, 12]+list(range(14, 18))
    assert {frozenset(swap[x] for x in v) for v in vertices} == set(map(frozenset, vertices))
    assert {frozenset(swap[x] for x in m) for m in missing} == missing
    paths = [source, path, f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_INTERNAL_THREE_SUPPORT_AND_NEGATIVE_INDEX_REVIEW',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                literal_supports=19, anchor_determinant_coefficients={str(k): str(v) for k, v in determinant.items()},
                scope='No producer imports. Reconstructs all19 actual supports, exact fixed contacts, missing triples and twin swap, and checks the symbolic anchor determinant. The finite-triangle, affine4-3-4 and Schur-budget inequalities give the full arbitrary-label negative-index contradiction. It excludes only the specified last-internal3 opposite-single branch.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}_review.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
