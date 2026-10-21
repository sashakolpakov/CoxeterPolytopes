#!/usr/bin/env python3
"""Small separate Fraction check supporting the all-real absorption review."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import json
import fresh_h20_cusp_sections_20261013_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_cusp_sections_20261013_absorption_review'
SOURCE = ROOT/'results/fresh_h20_cusp_sections_20261013_absorption.json'
PRODUCER = ROOT/'scripts/fresh_h20_cusp_sections_20261013_absorption.py'
OUT = ROOT/'results'/f'{TAG}.json'

def inverse(A):
    n = len(A)
    rows = [[F(x) for x in row]+[F(i == j) for j in range(n)] for i, row in enumerate(A)]
    for j in range(n):
        k = next(k for k in range(j, n) if rows[k][j])
        rows[j], rows[k] = rows[k], rows[j]
        t = rows[j][j]; rows[j] = [x/t for x in rows[j]]
        for k in range(n):
            if k != j:
                t = rows[k][j]
                rows[k] = [x-t*y for x, y in zip(rows[k], rows[j])]
    return [row[n:] for row in rows]

def run():
    d = json.loads(SOURCE.read_text())
    for p, h in d['sources'].items():
        assert R.sha(ROOT/p) == h
    raw = json.loads(R.RAW.read_text())['coxeter_input']
    Q = [[-F(x) for x in row] for row in raw['basis_gram']]
    W = [list(map(F, row)) for row in raw['wall_coordinates_as_columns_transposed']]
    C = [R.mv(Q, w) for w in W]
    G = [[R.dot(w, c) for c in C] for w in W]
    pair = lambda x, y: R.dot(x, R.mv(Q, y))
    delta = list(map(F, d['cusp']))
    assert pair(delta, delta) == 0
    assert R.add(W[162], W[164]) == R.scale(2, delta)
    D = [4, 9, 14, 26, 29, 48]
    centers = {29, 48}
    marks = {i: 2 if i in centers else 1 for i in D}
    B = inverse([[G[i][j] for j in D[:-1]] for i in D[:-1]])
    forms = []
    # Every star centered at either trivalent node is independently rebuilt.
    # This includes both orientations of the center-center edge.
    for center in sorted(centers):
        neighbors = [i for i in D if G[i][center] < 0]
        assert len(neighbors) == 3
        for size in [1, 2, 3]:
            for leaves in combinations(neighbors, size):
                columns = []
                for leaf in leaves:
                    values = {i: F(0) for i in D}
                    values[leaf] = F(1)
                    values[center] = -F(marks[leaf], marks[center])
                    columns.append([values[i] for i in D[:-1]])
                M = [[R.dot(x, R.mv(B, y)) for y in columns] for x in columns]
                expected = [[(F(1) if i in centers else F(1, 2)) if i == j else F(0) for j in leaves] for i in leaves]
                assert M == expected
                forms.append({'center': center, 'leaves': list(leaves), 'norm_matrix': [[str(x) for x in row] for row in M]})
    assert len(forms) == 14
    # Exhaust signed supports independently; no sparsity is assumed.
    support_counts = []
    for row in d['A7_opposite_supports'] + [{'walls': D, 'rows': d['D5_star_supports']}]:
        ids = row['walls']; actual = set()
        for signs in product([-1, 0, 1], repeat=len(ids)):
            pos = tuple(i for i, s in zip(ids, signs) if s > 0)
            neg = tuple(i for i, s in zip(ids, signs) if s < 0)
            if pos and neg and all(G[i][j] < 0 for i in pos for j in neg):
                assert min(len(pos), len(neg)) == 1
                actual.add((pos, neg))
        saved = {(tuple(r['positive']), tuple(r['negative'])) for r in row['rows']}
        assert actual == saved
        support_counts.append(len(actual))
    assert support_counts == [32, 32, 26]
    for row in d['A7_unique_update_recipes']:
        ids = row['walls']; a, b = row['opposite_pair']
        coeff = list(map(F, row['coefficients']))
        e = [sum(t*W[i][k] for t, i in zip(coeff, ids)) for k in range(22)]
        assert e == list(map(F, row['vector']))
        assert pair(e, e) == F(7, 8)
        assert [pair(e, W[i]) for i in ids] == [F(i == a)-F(i == b) for i in ids]
        assert max(coeff)-min(coeff) == F(7, 8)
    tangent = R.add(W[162], R.scale(-1, delta))
    assert pair(tangent, tangent) == 4
    assert (F(7, 8)**2)/F(1, 8) == F(d['width_squared_ratio']) == F(49, 8)
    plane = d['pure_A1_plane']; e = list(map(F, plane['base_normal']))
    assert e == R.scale(F(1, 2), tangent)
    assert pair(e, e) == 1 and pair(e, delta) == 0
    assert R.add(e, R.scale(F(1, 2), delta)) == R.scale(F(1, 2), W[162])
    assert R.add(e, R.scale(F(-1, 2), delta)) == R.scale(F(-1, 2), W[164])
    z = []
    for key in ['left_face', 'right_face']:
        row = plane[key]; v = list(map(F, row['point']))
        assert pair(v, v) == F(row['norm']) < 0 and pair(v, delta) < 0
        vals = [R.dot(v, c) for c in C]
        assert max(vals) == 0 and [i for i, q in enumerate(vals) if q == 0] == row['walls']
        z.append(v)
    left, right = z
    A, BB = -pair(left, delta), -pair(right, delta)
    assert [A, BB] == [F(plane['A']), F(plane['B'])] and min(A, BB) > 0
    assert pair(left, right) < 0
    # Coefficientwise identity for x(t)=(1/2-t)B zL+(t+1/2)A zR.
    x0 = R.add(R.scale(BB/2, left), R.scale(A/2, right))
    x1 = R.add(R.scale(-BB, left), R.scale(A, right))
    assert pair(x0, e) == pair(x1, e)+pair(x0, delta) == pair(x1, delta) == 0
    assert G[45][45] == G[63][63] == 2 and G[45][63] == 0
    assert [pair(e, W[i]) for i in [45, 63]] == [F(-1, 2), F(1, 2)]
    assert [pair(delta, W[i]) for i in [45, 63]] == [-1, -1]
    assert list(map(F, plane['projected_pairing_coefficients'])) == [F(1, 4), 0, -1]
    assert list(map(F, plane['projected_Gram_determinant_coefficients'])) == [3, 0, -4]
    paths = [SOURCE, PRODUCER, R.RAW, Path(R.__file__), Path(__file__)]
    return {'status': 'PASS_SEPARATE_FRACTION_A1_ABSORPTION_REVIEW',
            'sources': {str(p.relative_to(ROOT)): R.sha(p) for p in paths},
            'arithmetic': 'Separate Fraction reconstruction; imports only our older independent rational helper, not the absorption producer.',
            'signed_support_counts': support_counts, 'D5_star_quadratic_forms': forms,
            'A7_width_updates': 16, 'width_squared_ratio': '49/8',
            'actual_ridge_all210_inequalities_rechecked': True,
            'scope': 'Source/theorem review plus separate finite rational identities, all sign-support patterns, D5 norm forms, A7 widths, and pure-A1 actual-ridge interpolation. The all-real product/facet and inequality arguments are reviewed in the accompanying note; this is not a general H20 obstruction.'}

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    result = run()
    if a.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
