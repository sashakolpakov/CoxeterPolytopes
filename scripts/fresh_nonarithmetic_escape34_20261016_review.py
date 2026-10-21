#!/usr/bin/env python3
"""Independent nested-radical review of the tetrahedral-face cap."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json

import fresh_nonarithmetic_wall3_cap_20261016_review as P

ROOT, R, K, Z, O, V = P.ROOT, P.R, P.K, P.Z, P.O, P.V
X = P.Radical
D = 24*O+72*R.BASIS[1]/5
TAG = 'fresh_nonarithmetic_escape34_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_escape34_20261016.json'
PREVIOUS = P.PRIMARY
PRIOR_REVIEW = 'results/fresh_nonarithmetic_wall3_cap_20261016_review.json'


class Y:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        if isinstance(a, Y):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = X(a), X(b)
    def __add__(self, value):
        value = Y(value)
        return Y(self.a+value.a, self.b+value.b)
    __radd__ = __add__
    def __neg__(self):
        return Y(-self.a, -self.b)
    def __sub__(self, value):
        return self+-Y(value)
    def __rsub__(self, value):
        return Y(value)+-self
    def __mul__(self, value):
        value = Y(value)
        return Y(self.a*value.a+X(D)*self.b*value.b,
                 self.a*value.b+self.b*value.a)
    __rmul__ = __mul__
    def __truediv__(self, value):
        value = Y(value)
        assert value.b == 0 and not value.a.b and value.a.a
        return Y(self.a/value.a, self.b/value.a)
    def __eq__(self, value):
        value = Y(value)
        return self.a == value.a and self.b == value.b
    def sign(self):
        left, right = self.a.sign(), self.b.sign()
        if not left:
            return right
        if not right or left == right:
            return left
        difference = (self.a*self.a-X(D)*self.b*self.b).sign()
        return 0 if not difference else left if difference > 0 else right


def dy(v):
    return Y(P.dx(v[0]), P.dx(v[1]))


def vy(v):
    return list(map(dy, v))


def base(v):
    v = Y(v)
    assert v.b == 0 and not v.a.b
    return v.a.a


def dot(a, b):
    return sum((Y(x)*Y(y) for x, y in zip(a, b)), Y())


def mv(a, x):
    return [dot(row, x) for row in a]


def run():
    saved, old, prior = map(V.read, (PRIMARY, PREVIOUS, PRIOR_REVIEW))
    assert prior['status'] == 'PASS_INDEPENDENT_WALL3_CAP_FINITE_CENSUS_AND_RENEWED_ESCAPE'
    assert prior['sources'][PREVIOUS] == V.digest(PREVIOUS)
    assert saved['status'] == 'PASS_EXACT_ESCAPE_CUTTER_AND_COMPACT_TETRAHEDRAL_FACE'
    Q = [list(map(V.dec, row)) for row in old['form']]
    H, oldE = list(map(V.dec, old['center'])), list(map(V.dec, old['strict_escape']['point']))
    oldN = [P.vx(row) for row in old['normals']]
    pair = lambda a, b: dot(a, mv(Q, b))
    oldcov = [mv(Q, n) for n in oldN]
    support = [0, 1, 2, 3, 5, 6, 7, 12]
    rows = [[base(x) for x in oldcov[i]] for i in support+[8]]
    p, k = list(map(V.dec, saved['cap']['particular'])), list(map(V.dec, saved['cap']['kernel']))
    assert V.rank(rows) == 9 and V.mv(rows, p) == [Z]*8+[-O] and V.mv(rows, k) == [Z]*9
    assert pair(p, p) == Y(14*O+6*R.BASIS[1])
    assert pair(k, k) == Y(-(5*O+R.BASIS[1])/16) and pair(p, k) == 0
    assert R.sg(D) > 0
    r = [Y(a, -b) for a, b in zip(p, k)]
    assert r == vy(saved['cap']['normal']) and pair(r, r) == 2
    normals = [[Y(x) for x in n] for n in oldN]+[r]
    cov = [mv(Q, n) for n in normals]
    gram = [[dot(n, c) for c in cov] for n in normals]
    hp = [dot(H, c) for c in cov]
    assert pair(H, H).sign() < 0 and all(x.sign() < 0 for x in hp)
    assert pair(r, oldE).sign() > 0
    neighbors, counts = [], Counter()
    for i in range(34):
        g = gram[34][i]
        assert g.sign() <= 0
        if g == 0:
            kind, label = 'right', 2
        elif i == 8:
            assert g == -1 and gram[i][i] == 2
            kind, label = 'finite', 3
        else:
            assert (4*g*g/(gram[i][i]*gram[34][34])-4).sign() > 0
            kind, label = 'dotted', None
        rr = saved['cap']['relations'][i]
        assert (rr['kind'], rr['label']) == (kind, label) and dy(rr['pairing']) == g
        counts[kind] += 1
        if kind != 'dotted':
            neighbors.append(i)
    assert counts == dict(right=8, finite=1, dotted=25) and neighbors == sorted(support+[8])
    assert all(gram[i][j].sign() <= 0 for i in range(35) for j in range(i))
    for i in range(35):
        norm = gram[i][i]
        assert R.sg(base(norm)) > 0
        point = [Y(h)-hp[i]/norm*n for h, n in zip(H, normals[i])]
        assert point == vy(saved['facet_witnesses'][i]['point'])
        assert (pair(H, H)-hp[i]*hp[i]/norm).sign() < 0
        values = [hp[j]-hp[i]*gram[i][j]/norm for j in range(35)]
        assert [x.sign() for x in values] == [0 if j == i else -1 for j in range(35)]
    print('independent35 facets,34 contacts,old escape cut PASS', flush=True)
    for row in old['finite_vertices']:
        assert dot([Y(x) for x in P.vx(row['point'])], cov[34]).sign() < 0
    census = saved['complete_new_cap_finite_census']
    candidates = {tuple(row['support']): row for row in census['candidate_supports']}
    new = {tuple(row['incident']): row for row in census['vertices']}
    assert len(candidates) == 9 and len(new) == 3
    new_incidence = []
    for sub in combinations(neighbors, 8):
        ids = list(sub)+[34]
        G = [[base(gram[i][j]) for j in ids] for i in ids]
        pd = P.positive(G)
        assert pd == candidates[tuple(ids)]['positive_definite']
        if not pd:
            continue
        coeff = mv(P.inverse(G), [hp[i] for i in ids])
        point = [Y(H[t])-sum((a*normals[i][t] for a, i in zip(coeff, ids)), Y()) for t in range(10)]
        assert point == vy(new[tuple(ids)]['point']) and pair(point, point).sign() < 0
        assert pair(point, H).sign() < 0
        assert [dot(point, c).sign() for c in cov] == [0 if i in ids else -1 for i in range(35)]
        new_incidence.append(ids)
    incidence = old['finite_vertex_incidence']+new_incidence
    assert incidence == saved['finite_vertex_incidence'] and len(incidence) == 98
    edges = []
    for i, j in combinations(range(98), 2):
        common = sorted(set(incidence[i]) & set(incidence[j]))
        if len(common) == 8:
            edges.append(dict(vertices=[i, j], support=common))
    assert edges == saved['finite_edges'] and len(edges) == 263
    assert 9*98-2*263 == saved['open_finite_germ_count'] == 356
    print('independent98 finite vertices,263 edges PASS', flush=True)
    face = [0, 1, 2, 3, 5, 12]
    boundary = [6, 7, 8, 34]
    FG = [[base(gram[i][j]) for j in face] for i in face]
    assert P.positive(FG)
    FI = P.inverse(FG)
    projected = []
    for i in boundary:
        coeff = V.mv(FI, [base(gram[i][j]) for j in face])
        projected.append([normals[i][t]-sum((Y(a)*normals[j][t] for a, j in zip(coeff, face)), Y()) for t in range(10)])
    BG = [[base(pair(a, b)) for b in projected] for a in projected]
    assert BG == [list(map(V.dec, row)) for row in saved['compact_three_face']['projected_boundary_gram']]
    assert all(P.positive([[BG[i][j] for j in sub] for i in sub]) for sub in combinations(range(4), 3))
    BI = P.inverse(BG)
    assert all(R.sg(x) < 0 for row in BI for x in row)
    duals = []
    for omitted in range(4):
        point = [-sum((Y(BI[i][omitted])*projected[i][t] for i in range(4)), Y()) for t in range(10)]
        active = sorted(face+[boundary[j] for j in range(4) if j != omitted])
        assert pair(point, point).sign() < 0 and pair(point, H).sign() < 0
        assert [dot(point, c).sign() for c in cov] == [0 if j in active else -1 for j in range(35)]
        assert [dot(point, cov[j]) for j in boundary] == [Y(-int(j == omitted)) for j in range(4)]
        duals.append(point)
    escape = saved['strict_escape']
    active = [0, 1, 2, 4, 5, 6, 7, 9, 17]
    rows = [[base(x) for x in cov[i]] for i in active]+[V.mv(Q, H)]
    ray = V.mv(P.inverse(rows), [Z]*9+[-O])
    assert ray == list(map(V.dec, escape['ray']))
    E = [x+K.convert(F(1, 8192))*h for x, h in zip(ray, H)]
    assert E == list(map(V.dec, escape['point']))
    assert pair(E, E).sign() > 0 and pair(E, H).sign() < 0
    assert all(dot(E, c).sign() < 0 for c in cov)
    coeff = V.mv(FI, [base(dot(E, cov[i])) for i in face])
    EP = [Y(E[t])-sum((Y(a)*normals[i][t] for a, i in zip(coeff, face)), Y()) for t in range(10)]
    assert EP == vy(escape['compact_face_projection']) and pair(EP, EP).sign() < 0
    bary = [-dot(EP, cov[i]) for i in boundary]
    assert all(a.sign() > 0 for a in bary)
    assert EP == [sum((a*v[t] for a, v in zip(bary, duals)), Y()) for t in range(10)]
    paths = [f'scripts/{TAG}.py', PRIMARY, PREVIOUS, PRIOR_REVIEW,
             'scripts/fresh_nonarithmetic_wall3_cap_20261016_review.py',
             'scripts/fresh_nonarithmetic_two_caps_20261016_review.py',
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_TETRAHEDRAL_CAP_AND_COMPACT_FACE_ESCAPE_GATE',
                sources={p: V.digest(p) for p in paths}, actual_facets=35,
                new_pair_counts=dict(counts), finite_vertices=98, finite_edges=263,
                open_finite_germs=356, compact_three_face_support=face,
                compact_three_face_vertices=[1, 95, 96, 97],
                old_escape_cut=True, renewed_strict_escape=True,
                all_real_perpendicular_escape_cutter_excluded=True,
                scope='Separate nested radical arithmetic and inherited independent octic engine. All new pairs/facets,95old vertices,complete9-support new-facet census,98vertices/263edges,complete compact tetrahedral3face and strict renewed escape checked. Positive dual decomposition supports the all-real fixed-six-wall-support cutter obstruction. No intrinsic Coxeter-face angle claim, complete exterior-ray census or finite-volume polytope.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
