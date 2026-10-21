#!/usr/bin/env python3
"""Independent literal review of the new wall3 cap and finite-corner update.

No new-cap producer import. The frozen octic arithmetic is shared; formal
radical operations, linear algebra and all new geometric checks are here.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json

import fresh_nonarithmetic_two_caps_20261016_review as V

ROOT, R, K, Z, O = V.ROOT, V.R, V.K, V.Z, V.O
TAG = 'fresh_nonarithmetic_wall3_cap_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_wall3_cap_20261016.json'
OLD_REVIEW = 'results/fresh_h9_joint_corner_20261015_bridge_review.json'
D = O+R.BASIS[1]


class Radical:
    __slots__ = ('a', 'b')

    def __init__(self, a=Z, b=Z):
        if isinstance(a, Radical):
            self.a, self.b = a.a, a.b
        else:
            self.a = a if hasattr(a, 'to_list') else K.convert(a)
            self.b = b if hasattr(b, 'to_list') else K.convert(b)

    def __add__(self, other):
        other = Radical(other)
        return Radical(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Radical(-self.a, -self.b)

    def __sub__(self, other):
        return self+-Radical(other)

    def __rsub__(self, other):
        return Radical(other)+-self

    def __mul__(self, other):
        other = Radical(other)
        return Radical(self.a*other.a+D*self.b*other.b,
                       self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = Radical(scalar)
        assert not scalar.b and scalar.a
        return Radical(self.a/scalar.a, self.b/scalar.a)

    def __eq__(self, other):
        other = Radical(other)
        return self.a == other.a and self.b == other.b

    def sign(self):
        return V.sign_radical_sum(self.a, O, self.b, D)


def dx(x):
    return Radical(V.dec(x[0]), V.dec(x[1]))


def vx(x):
    return list(map(dx, x))


def dot(a, b):
    assert len(a) == len(b)
    return sum((Radical(x)*Radical(y) for x, y in zip(a, b)), Radical())


def mv(a, v):
    return [dot(row, v) for row in a]


def inverse(a):
    n = len(a)
    rows = [row[:]+[O if i == j else Z for j in range(n)] for i, row in enumerate(a)]
    for j in range(n-1, -1, -1):
        pivot = next(i for i in range(j+1) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x/divisor for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[j])]
    assert [row[:n] for row in rows] == [[O if i == j else Z for j in range(n)] for i in range(n)]
    return [row[n:] for row in rows]


def positive(a):
    a = [row[:] for row in a]
    while a:
        pivot = a[-1][-1]
        if R.sg(pivot) <= 0:
            return False
        a = [[a[i][j]-a[i][-1]*a[-1][j]/pivot for j in range(len(a)-1)] for i in range(len(a)-1)]
    return True


def run():
    saved, base, bridge, old = map(V.read, (PRIMARY, V.BASE, V.BRIDGE, OLD_REVIEW))
    assert saved['status'] == 'PASS_EXACT_ESCAPE_CUTTING_WALL3_CAP_AND_FINITE_CENSUS'
    assert old['status'] == 'PASS_INDEPENDENT_OCTIC_BRIDGE_AND_COMPLETE_FINITE_VERTEX_CENSUS'
    assert old['sources'][V.BRIDGE] == V.digest(V.BRIDGE)
    assert (old['finite_vertices'], old['finite_edges']) == (92, 253)
    Q, N, H = R.m4(base['form']), R.m4(base['normals']), R.v4(base['center'])
    N.append(R.v8(bridge['cap']['normal']))
    oldE = R.v4(base['strict_escape']['point'])
    pair = lambda a, b: dot(a, mv(Q, b))
    basepair = lambda a, b: V.dot(a, V.mv(Q, b))
    assert Q == [list(map(V.dec, row)) for row in saved['form']]
    assert N == [list(map(V.dec, row)) for row in saved['old_normals']]
    p = list(map(V.dec, saved['cap']['particular']))
    w = list(map(V.dec, saved['cap']['radical_vector']))
    ids = [0, 1, 2, 5, 6, 10, 30, 3, 7]
    cov_old = [V.mv(Q, n) for n in N]
    rows = [cov_old[i] for i in ids]
    assert V.rank(rows) == 9
    assert V.mv(rows, p) == [Z]*8+[-O] and V.mv(rows, w) == [Z]*9
    assert basepair(p, w) == Z
    assert basepair(p, p) == -8*O-4*R.BASIS[1]
    assert basepair(w, w) == (5*O+3*R.BASIS[1])/2
    assert R.sg(D) > 0
    r = [Radical(a, -b) for a, b in zip(p, w)]
    assert r == vx(saved['cap']['normal']) and pair(r, r) == 2
    assert pair(r, H).sign() < 0 and pair(r, oldE).sign() > 0
    normals = [[Radical(a) for a in n] for n in N]+[r]
    cov = [mv(Q, n) for n in normals]
    norms = [dot(n, c) for n, c in zip(normals, cov)]
    gram = [[dot(n, c) for c in cov] for n in normals]
    assert all(not q.b and R.sg(q.a) > 0 for q in norms)
    counts, neighbors = Counter(), []
    for i in range(33):
        g = gram[33][i]
        cycle = 4*g*g/(norms[i]*norms[33])
        assert g.sign() <= 0
        if g == 0:
            kind, label = 'right', 2
        elif i == 7:
            assert g == -1 and cycle == 1
            kind, label = 'finite', 3
        else:
            assert (cycle-4).sign() > 0
            kind, label = 'dotted', None
        row = saved['cap']['relations'][i]
        assert row['wall'] == i and (row['kind'], row['label']) == (kind, label)
        assert dx(row['pairing']) == g and dx(row['cartan_cycle']) == cycle
        counts[kind] += 1
        if kind != 'dotted':
            neighbors.append(i)
    assert counts == dict(right=8, finite=1, dotted=24) and neighbors == sorted(ids)
    hp = [dot(H, c) for c in cov]
    assert pair(H, H).sign() < 0 and all(x.sign() < 0 for x in hp)
    assert len(saved['facet_witnesses']) == 34
    for i, normal in enumerate(normals):
        point = [Radical(h)-hp[i]/norms[i]*a for h, a in zip(H, normal)]
        values = [dot(point, c) for c in cov]
        assert point == vx(saved['facet_witnesses'][i]['point'])
        assert pair(point, point).sign() < 0 and pair(point, H).sign() < 0
        assert [x.sign() for x in values] == [0 if j == i else -1 for j in range(34)]
    print('independent cap,33 contacts,34 actual facets PASS', flush=True)
    old_points = [R.v4(row['point']) for row in base['finite_vertices']]
    old_points += [R.v8(row['point']) for row in bridge['complete_new_cap_finite_census']['vertices']]
    assert len(old_points) == 92
    for point, raw in zip(old_points, saved['old_finite_vertex_preservation_pairings']):
        value = dot(point, cov[33])
        assert value.sign() < 0 and value == dx(raw)
    census = saved['complete_new_cap_finite_census']
    candidates = {tuple(row['support']): row for row in census['candidate_supports']}
    new_rows = {tuple(row['incident']): row for row in census['vertices']}
    assert len(candidates) == 9 and len(new_rows) == 3
    new_incidence = []
    for sub in combinations(neighbors, 8):
        support = list(sub)+[33]
        active = [[gram[i][j] for j in support] for i in support]
        assert all(not x.b for row in active for x in row)
        active = [[x.a for x in row] for row in active]
        pd = positive(active)
        assert pd == candidates[tuple(support)]['positive_definite']
        if not pd:
            assert tuple(support) not in new_rows
            continue
        coefficients = mv(inverse(active), [hp[i] for i in support])
        point = [Radical(H[k])-sum((a*normals[i][k] for a, i in zip(coefficients, support)), Radical()) for k in range(10)]
        values = [dot(point, c) for c in cov]
        assert pair(point, point).sign() < 0 and pair(point, H).sign() < 0
        assert [x.sign() for x in values] == [0 if i in support else -1 for i in range(34)]
        assert point == vx(new_rows[tuple(support)]['point'])
        new_incidence.append(support)
    incidence = bridge['finite_vertex_incidence']+new_incidence
    assert incidence == saved['finite_vertex_incidence'] and len(incidence) == 95
    edges = []
    for i, j in combinations(range(95), 2):
        common = sorted(set(incidence[i]) & set(incidence[j]))
        if len(common) == 8:
            edges.append(dict(vertices=[i, j], support=common))
    assert edges == saved['finite_edges'] and len(edges) == 257
    assert 9*len(incidence)-2*len(edges) == saved['open_finite_germ_count'] == 341
    print('independent92 preserved+3 new finite vertices,257 edges PASS', flush=True)
    escape = saved['strict_escape']
    support = [0, 1, 2, 3, 5, 6, 7, 12, 30]
    assert support == base['complete_census']['hyperideal_rays'][1] == escape['ray_active']
    rows = [cov_old[i] for i in support]+[V.mv(Q, H)]
    ray = V.mv(inverse(rows), [Z]*9+[-O])
    assert ray == list(map(V.dec, escape['inherited_hyperideal_ray']))
    assert R.sg(basepair(ray, ray)) > 0 and basepair(ray, H) == -O
    point = [v+K.convert(F(1, 8192))*h for v, h in zip(ray, H)]
    assert point == list(map(V.dec, escape['point']))
    assert R.sg(basepair(point, point)) > 0 and R.sg(basepair(point, H)) < 0
    assert all(dot(point, c).sign() < 0 for c in cov)
    paths = [f'scripts/{TAG}.py', PRIMARY, V.BASE, V.BRIDGE, OLD_REVIEW,
             'scripts/fresh_nonarithmetic_two_caps_20261016_review.py',
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_WALL3_CAP_FINITE_CENSUS_AND_RENEWED_ESCAPE',
                sources={p: V.digest(p) for p in paths}, actual_facets=34,
                new_pair_counts=dict(counts), old_finite_vertices_preserved=92,
                candidate_new_facet_supports=9, new_finite_vertices=3,
                finite_vertices=95, finite_edges=257, open_finite_germs=341,
                old_escape_cut=True, strict_spacelike_escape=True, finite_volume=False,
                scope='Separate formal radical arithmetic and linear algebra, with frozen octic field arithmetic/sign evaluator shared. All33 new contacts,34 facets,92 vertex preservation signs,complete9-support new-facet census,finite edges and strict renewed escape checked. Old33-wall complete finite census is inherited; no complete new exterior-ray census or finite-volume example.')


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
