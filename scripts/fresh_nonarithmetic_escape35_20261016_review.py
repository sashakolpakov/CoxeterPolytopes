#!/usr/bin/env python3
"""Independent exact cap, finite-corner and complete Lambert-face review."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json

import fresh_nonarithmetic_escape34_20261016_review as P

ROOT, R, K, O, Z, V = P.ROOT, P.R, P.K, P.O, P.Z, P.V
INNER = P.Y
D = 10*O+68*R.BASIS[1]/15
TAG = 'fresh_nonarithmetic_escape35_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_escape35_20261016.json'
PREVIOUS = P.PRIMARY
PRIOR_REVIEW = 'results/fresh_nonarithmetic_escape34_20261016_review.json'


class T:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        if isinstance(a, T):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = INNER(a), INNER(b)
    def __add__(self, value):
        value = T(value)
        return T(self.a+value.a, self.b+value.b)
    __radd__ = __add__
    def __neg__(self):
        return T(-self.a, -self.b)
    def __sub__(self, value):
        return self+-T(value)
    def __rsub__(self, value):
        return T(value)+-self
    def __mul__(self, value):
        value = T(value)
        return T(self.a*value.a+INNER(D)*self.b*value.b,
                 self.a*value.b+self.b*value.a)
    __rmul__ = __mul__
    def __truediv__(self, value):
        value = T(value)
        assert value.b == 0
        q = P.base(value.a)
        assert q
        return T(self.a/q, self.b/q)
    def __eq__(self, value):
        value = T(value)
        return self.a == value.a and self.b == value.b
    def sign(self):
        left, right = self.a.sign(), self.b.sign()
        if not left:
            return right
        if not right or left == right:
            return left
        difference = (self.a*self.a-INNER(D)*self.b*self.b).sign()
        return 0 if not difference else left if difference > 0 else right


def dt(v):
    return T(P.dy(v[0]), P.dy(v[1]))


def vt(v):
    return list(map(dt, v))


def base(v):
    v = T(v)
    assert v.b == 0
    return P.base(v.a)


def dot(a, b):
    return sum((T(x)*T(y) for x, y in zip(a, b)), T())


def mv(a, x):
    return [dot(row, x) for row in a]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def run():
    saved, old, prior = map(V.read, (PRIMARY, PREVIOUS, PRIOR_REVIEW))
    assert prior['status'] == 'PASS_INDEPENDENT_TETRAHEDRAL_CAP_AND_COMPACT_FACE_ESCAPE_GATE'
    assert prior['sources'][PREVIOUS] == V.digest(PREVIOUS)
    assert saved['status'] == 'PASS_EXACT_ESCAPE_CUTTER_AND_COMPACT_QUADRILATERAL_FACE'
    Q = [list(map(V.dec, row)) for row in old['form']]
    H, oldE = list(map(V.dec, old['center'])), list(map(V.dec, old['strict_escape']['point']))
    normals = [[T(x) for x in P.vy(row)] for row in old['normals']]
    pair = lambda a, b: dot(a, mv(Q, b))
    cov = [mv(Q, n) for n in normals]
    support = [0, 1, 2, 4, 5, 6, 7, 9]
    rows = [[base(x) for x in cov[i]] for i in support+[8]]
    assert V.rank(rows) == 9
    k = [Z]*9+[O]
    assert V.mv(rows, k) == [Z]*9
    p = V.mv(P.P.inverse(rows+[V.mv(Q, k)]), [Z]*8+[-O, Z])
    assert p == list(map(V.dec, saved['cap']['particular']))
    assert k == list(map(V.dec, saved['cap']['kernel']))
    assert pair(p, p) == T(-5*O-5*R.BASIS[1]/2)
    assert pair(k, k) == T(R.BASIS[1]/4) and pair(p, k) == 0 and R.sg(D) > 0
    r = [T(a, -b) for a, b in zip(p, k)]
    assert r == vt(saved['cap']['normal']) and pair(r, r) == T(K.convert(F(2, 3)))
    normals.append(r)
    cov.append(mv(Q, r))
    gram = [[dot(n, c) for c in cov] for n in normals]
    hp = [dot(H, c) for c in cov]
    assert pair(H, H).sign() < 0 and all(x.sign() < 0 for x in hp)
    assert pair(r, oldE).sign() > 0
    assert all(gram[i][j].sign() <= 0 for i in range(36) for j in range(i))
    neighbors, counts = [], Counter()
    for i in range(35):
        g = gram[35][i]
        cycle = 4*g*g/(gram[i][i]*gram[35][35])
        if g == 0:
            kind, label = 'right', 2
        elif i == 8:
            assert g == -1 and cycle == 3
            kind, label = 'finite', 6
        else:
            assert (cycle-4).sign() > 0
            kind, label = 'dotted', None
        rr = saved['cap']['relations'][i]
        assert (rr['kind'], rr['label']) == (kind, label) and dt(rr['pairing']) == g
        counts[kind] += 1
        if kind != 'dotted':
            neighbors.append(i)
    assert counts == dict(right=8, finite=1, dotted=26)
    assert neighbors == sorted(support+[8])
    for i in range(36):
        norm = gram[i][i]
        assert R.sg(base(norm)) > 0
        point = [T(h)-hp[i]/norm*n for h, n in zip(H, normals[i])]
        assert point == vt(saved['facet_witnesses'][i]['point'])
        assert (pair(H, H)-hp[i]*hp[i]/norm).sign() < 0
        assert [x.sign() for x in [hp[j]-hp[i]*gram[i][j]/norm for j in range(36)]] == [0 if j == i else -1 for j in range(36)]
    print('independent36 facets,35 new contacts and escape cut PASS', flush=True)
    for row in old['finite_vertices']:
        assert dot([T(x) for x in P.vy(row['point'])], cov[35]).sign() < 0
    census = saved['complete_new_cap_finite_census']
    candidates = {tuple(row['support']):row for row in census['candidate_supports']}
    new = {tuple(row['incident']):row for row in census['vertices']}
    assert len(candidates) == 9 and len(new) == 2
    added = []
    for sub in combinations(neighbors, 8):
        ids = list(sub)+[35]
        G = [[base(gram[i][j]) for j in ids] for i in ids]
        positive = P.P.positive(G)
        assert positive == candidates[tuple(ids)]['positive_definite']
        if not positive:
            continue
        coeff = mv(P.P.inverse(G), [hp[i] for i in ids])
        point = [T(H[t])-sum((a*normals[i][t] for a, i in zip(coeff, ids)), T()) for t in range(10)]
        assert point == vt(new[tuple(ids)]['point']) and pair(point, point).sign() < 0
        assert pair(point, H).sign() < 0
        assert [dot(point, c).sign() for c in cov] == [0 if i in ids else -1 for i in range(36)]
        added.append(ids)
    incidence = old['finite_vertex_incidence']+added
    assert incidence == saved['finite_vertex_incidence'] and len(incidence) == 100
    edges = []
    for i, j in combinations(range(100), 2):
        common = sorted(set(incidence[i]) & set(incidence[j]))
        if len(common) == 8:
            edges.append(dict(vertices=[i, j], support=common))
    assert edges == saved['finite_edges'] and len(edges) == 266
    assert 9*100-2*266 == saved['open_finite_germ_count'] == 368
    print('independent complete100 finite vertices,266 edges PASS', flush=True)
    face = [0, 1, 2, 4, 5, 6, 9]
    boundary = [7, 8, 11, 35]
    G = [[base(gram[i][j]) for j in face] for i in face]
    assert P.P.positive(G)
    J = P.P.inverse(G)
    projected = []
    for i in boundary:
        cs = V.mv(J, [base(gram[i][j]) for j in face])
        projected.append([normals[i][t]-sum((T(a)*normals[j][t] for a,j in zip(cs,face)), T()) for t in range(10)])
    BG = [[pair(a, b) for b in projected] for a in projected]
    assert BG == [[dt(x) for x in row] for row in saved['compact_two_face']['projected_boundary_gram']]
    frame = projected[:3]
    frameG = [[base(pair(a, b)) for b in frame] for a in frame]
    assert V.rank(frameG) == 3
    fcov = [[pair(n, b) for b in frame] for n in projected]
    feasible = []
    for i, j in combinations(range(4), 2):
        cs = cross(fcov[i], fcov[j])
        ray = [sum((a*v[t] for a,v in zip(cs,frame)), T()) for t in range(10)]
        assert any(a.sign() for a in ray)
        assert pair(ray, projected[i]) == pair(ray, projected[j]) == 0
        good = []
        for direction in (-1, 1):
            v = [direction*x for x in ray]
            if all(pair(v, n).sign() <= 0 for n in projected):
                good.append(v)
        if not good:
            assert pair(ray, ray).sign() > 0
            continue
        assert len(good) == 1
        ray = good[0]
        active = sorted(face+[boundary[i],boundary[j]])
        assert pair(ray, ray).sign() < 0 and pair(ray,H).sign() < 0
        assert [dot(ray,c).sign() for c in cov] == [0 if k in active else -1 for k in range(36)]
        cycle = 4*BG[i][j]*BG[i][j]/(BG[i][i]*BG[j][j])
        assert cycle == 0 or cycle == 3
        feasible.append(dict(active=active, boundary_pair=[boundary[i],boundary[j]], angle_label=2 if cycle == 0 else 6))
    assert len(feasible) == 4 and sorted(row['angle_label'] for row in feasible) == [2,2,2,6]
    assert {tuple(row['active']) for row in feasible} == {tuple(incidence[i]) for i in [6,8,98,99]}
    active = [0,1,2,4,5,6,7,10,27]
    rows = [[base(x) for x in cov[i]] for i in active]+[V.mv(Q,H)]
    ray = V.mv(P.P.inverse(rows), [Z]*9+[-O])
    E = [v+K.convert(F(1,8192))*h for v,h in zip(ray,H)]
    assert E == list(map(V.dec,saved['strict_escape']['point']))
    assert pair(E,E).sign() > 0 and pair(E,H).sign() < 0
    assert all(dot(E,c).sign() < 0 for c in cov)
    cs = V.mv(J,[base(dot(E,cov[i])) for i in face])
    EP = [T(E[t])-sum((T(a)*normals[i][t] for a,i in zip(cs,face)),T()) for t in range(10)]
    assert EP == vt(saved['strict_escape']['compact_face_projection'])
    assert pair(EP,EP).sign() < 0 and pair(EP,H).sign() < 0
    assert [dot(EP,c).sign() for c in cov] == [0 if i in face else -1 for i in range(36)]
    paths = [f'scripts/{TAG}.py', PRIMARY, PREVIOUS, PRIOR_REVIEW,
             'scripts/fresh_nonarithmetic_escape34_20261016_review.py',
             'scripts/fresh_nonarithmetic_wall3_cap_20261016_review.py',
             'scripts/fresh_nonarithmetic_two_caps_20261016_review.py',
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_LAMBERT_CAP_FINITE_CENSUS_AND_RENEWED_ESCAPE',
                sources={p:V.digest(p) for p in paths}, actual_facets=36,
                new_pair_counts=dict(counts), finite_vertices=100,finite_edges=266,
                open_finite_germs=368, complete_two_face_corners=feasible,
                renewed_strict_escape=True, compact_face_projection_inside=True,
                scope='Independent outer radical arithmetic over frozen earlier independent engines; all contacts/facets,98old vertices,complete9support new-facet census,complete6pair quadrilateral line census and strict escape. Compact Lambert face angles verified. Inherited complete finite census; no complete ambient exterior-ray census or finite-volume claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
