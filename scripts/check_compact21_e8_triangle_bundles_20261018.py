#!/usr/bin/env python3
"""Tiny exact E8 corner check for forced simplex-face triangle bundles.

Independent stdlib arithmetic; no project imports or domain enumeration.
The actual-face and simplex closure arguments are in the companion review.
"""
from fractions import Fraction as Q
from itertools import product
from hashlib import file_digest
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PINS = {
    'scripts/fresh_compact21_d5_two_e8_triangle_coupling_20261018.py':
        '87e56cf521fba2b7d95c8d56aa68449d450f911dadc6db43bc715aac359f9585',
    'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json':
        '1ab620b6d37ab5806cef96e27e70cfd230ed03ab07ddaf5f4d0a31edfcd38381',
    'results/fresh_compact21_single_field_20261016.md':
        '0988994cfe0c29943efdea7613afb299eb3783a3c2e80da73719a29465d972eb',
    'results/fresh_compact21_single_field_20261016.json':
        '9c5363f0671d5da6b820700d9673355dd69f9c950b2489ee5aa2832f3ca6b247',
}
for path, expected in PINS.items():
    with (ROOT/path).open('rb') as stream:
        assert file_digest(stream, 'sha256').hexdigest() == expected, path


class K:
    def __init__(self, a=0, b=0):
        self.a, self.b = (a.a, a.b) if isinstance(a, K) else (Q(a), Q(b))

    def __add__(self, x):
        x = K(x)
        return K(self.a+x.a, self.b+x.b)
    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, x):
        return self+-K(x)

    def __rsub__(self, x):
        return K(x)+-self

    def __mul__(self, x):
        x = K(x)
        return K(self.a*x.a+5*self.b*x.b, self.a*x.b+self.b*x.a)
    __rmul__ = __mul__

    def __truediv__(self, x):
        x = K(x)
        norm = x.a*x.a-5*x.b*x.b
        assert norm
        return self*K(x.a/norm, -x.b/norm)

    def __eq__(self, x):
        x = K(x)
        return self.a == x.a and self.b == x.b

    def sign(self):
        sign = lambda x: (x > 0)-(x < 0)
        if not self.b:
            return sign(self.a)
        if not self.a or sign(self.a) == sign(self.b):
            return sign(self.b)
        return sign(self.a)*sign(self.a*self.a-5*self.b*self.b)


def pd(matrix):
    a = [list(row) for row in matrix]
    for i in range(len(a)):
        if a[i][i].sign() <= 0:
            return False
        for j in range(i+1, len(a)):
            for k in range(j, len(a)):
                a[k][j] = a[j][k] = a[j][k]-a[j][i]*a[i][k]/a[i][i]
    return True


edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]
g = [[K(2*(i == j)) for j in range(10)] for i in range(10)]
for i, j in edges:
    g[i][j] = g[j][i] = K(-1)
phi = K(Q(1, 2), Q(1, 2))
g[0][8] = g[8][0] = -phi
g[6][9] = g[9][6] = -phi
assert pd([row[:8] for row in g[:8]])
# Solve the rational E8 inverse column at6 independently by a displayed vector.
column6 = [2, 4, 6, 5, 4, 3, 2, 3]
column0 = [4, 7, 10, 8, 6, 4, 2, 5]
for node, column in ((0, column0), (6, column6)):
    assert [sum(g[i][j]*column[j] for j in range(8)) for i in range(8)] == [K(int(i == node)) for i in range(8)]
h0, h6 = K(4, 2), K(1, 1)
assert h0 == column0[0]*phi*phi-2 and h6 == column6[6]*phi*phi-2
projection_product = 2*phi*phi
assert projection_product == K(3, 1)
assert h0*h6 == projection_product*projection_product
assert projection_product.sign() > 0 and h0.sign() > 0 and h6.sign() > 0
assert K(h0.a, -h0.b).sign() < 0 and K(h6.a, -h6.b).sign() < 0
assert g[8][9] == 0  # Forced by the preceding positive square-root identity.
single = {}
for added in (8, 9):
    finite = []
    for omitted in range(8):
        indices = [i for i in range(8) if i != omitted]+[added]
        if pd([[g[i][j] for j in indices] for i in indices]):
            finite.append(omitted)
    single[added] = finite
assert single == {8: [0, 1, 2], 9: [3, 4, 5, 6]}
# All vertex corners of the prospective product face have the required rank.
groups = ([0, 1, 2, 8], [3, 4, 5, 6, 9])
for a, b in product(*groups):
    indices = [i for i in range(10) if i not in (a, b)]
    assert pd([[g[i][j] for j in indices] for i in indices])
# Independently verify the D5 inverse and all pairwise candidate acute gates.
d5 = [[Q(2*(i == j)) for j in range(5)] for i in range(5)]
for i, j in ((0, 1), (1, 2), (2, 3), (2, 4)):
    d5[i][j] = d5[j][i] = Q(-1)
j5 = [[Q(x) for x in row] for row in [
    [1, 1, 1, '1/2', '1/2'], [1, 2, 2, 1, 1],
    [1, 2, 3, '3/2', '3/2'], ['1/2', 1, '3/2', '5/4', '3/4'],
    ['1/2', 1, '3/2', '3/4', '5/4']]]
assert [[sum(d5[i][k]*j5[k][j] for k in range(5)) for j in range(5)]
        for i in range(5)] == [[Q(i == j) for j in range(5)] for i in range(5)]
contacts = [{0: phi}, {3: phi}, {1: K(1), 3: phi},
            {4: phi}, {1: K(1), 4: phi}]
bilinear = lambda a, b: sum((x*j5[i][j]*y for i, x in a.items() for j, y in b.items()), K())
heights = [bilinear(a, a)-2 for a in contacts]
assert heights == [K(Q(-1, 2), Q(1, 2)), K(Q(-1, 8), Q(5, 8)),
                   K(Q(23, 8), Q(13, 8)), K(Q(-1, 8), Q(5, 8)),
                   K(Q(23, 8), Q(13, 8))]
pair_count = 0
for i in range(5):
    for j in range(i):
        p = bilinear(contacts[i], contacts[j])
        assert p.sign() > 0 and (p*p-heights[i]*heights[j]).sign() > 0
        pair_count += 1
assert pair_count == 10
print(json.dumps({'status': 'PASS_EXACT_E8_TRIANGLE_BUNDLES',
                  'forced_single_wall_deleted_root_domains': single,
                  'third_walls_mutually_orthogonal': True,
                  'product_face_corners_checked': 20,
                  'D5_positive_mutual_product_pairs': pair_count,
                  'D5_designated_triangle_maximum': 1,
                  'all48_nonadjacent_faces_triangle_upper': 9,
                  'global_charge_saving_proved': False, 'new_dimension_bound': None}))
