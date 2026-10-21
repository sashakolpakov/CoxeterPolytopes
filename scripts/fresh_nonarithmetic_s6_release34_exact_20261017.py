#!/usr/bin/env python3
"""Exact three-cap patch after deleting34 and replacing52.

The free dotted product x=44.52 is an isolated algebraic root of the
polar-angle equation. Interval signs and a derivative bound certify that
root; interval overlap is never used to establish an angle equality.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_escape62_20261017_cap as C

P, U, I, D = C.P, C.U, C.I, C.D
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_s6_release34_exact_20261017'
OLD_CONE = 'fresh_nonarithmetic_escape62_20261017_cone'
S = [0, 1, 2, 3, 5, 6]
NEW = [52, 63, 64, 65]


class Jet:
    """First derivative with outward rational interval bounds."""
    def __init__(self, x=0, derivative=0):
        if isinstance(x, Jet):
            self.v, self.d = x.v, x.d
        else:
            self.v, self.d = I(x), I(derivative)

    def __add__(self, other):
        other = Jet(other)
        return Jet(U.rnd(self.v+other.v), U.rnd(self.d+other.d))
    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d)

    def __sub__(self, other):
        return self+-Jet(other)

    def __rsub__(self, other):
        return Jet(other)+-self

    def __mul__(self, other):
        other = Jet(other)
        return Jet(U.rnd(self.v*other.v), U.rnd(self.d*other.v+self.v*other.d))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Jet(other)
        return Jet(U.rnd(self.v/other.v),
                   U.rnd((self.d*other.v-self.v*other.d)/(other.v*other.v)))

    def __rtruediv__(self, other):
        return Jet(other)/self

    def sqrt(self):
        q = U.sqrt(self.v)
        return Jet(q, U.rnd(self.d/(2*q)))


def jdot(a, b):
    return sum((x*y for x, y in zip(a, b)), Jet(0))


def jmv(a, b):
    return [jdot(row, b) for row in a]


def jinverse(matrix):
    n = len(matrix)
    a = [[Jet(x) for x in row]+[Jet(int(i == j)) for j in range(n)]
         for i, row in enumerate(matrix)]
    for k in range(n):
        choices = [i for i in range(k, n) if a[i][k].v.l*a[i][k].v.h > 0]
        assert choices, ('no certified jet pivot', k)
        j = max(choices, key=lambda i:min(abs(a[i][k].v.l), abs(a[i][k].v.h)))
        a[k], a[j] = a[j], a[k]
        pivot = a[k][k]
        a[k] = [x/pivot for x in a[k]]
        for i in range(n):
            if i == k:
                continue
            factor = a[i][k]
            a[i] = [x-factor*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def polar_equation(x, b):
    x, b = Jet(x, 1), Jet(b)
    zero, two, t, a = Jet(0), Jet(2), Jet(D.theta(4)), Jet(D.theta(15))
    j = jinverse([[two, zero, zero], [zero, two, x], [zero, x, two]])
    contacts, old30 = [zero, -a, -t], [b, zero, -t]
    rho = 2-jdot(contacts, jmv(j, contacts))
    w2 = 2-jdot(old30, jmv(j, old30))
    z = jdot(contacts, jmv(j, old30))-(rho*w2).sqrt()
    g = [[two, zero, -t, z], [zero, two, x, -a],
         [-t, x, two, -t], [z, -a, -t, two]]
    jj = jinverse(g)
    value = 2*jj[1][2]/(jj[1][1]*jj[2][2]).sqrt()+Jet(D.theta(27))
    return value


def isolate(b):
    lo, hi = F('-4.098064174'), F('-4.098064173')
    left, right = polar_equation(I(lo), b), polar_equation(I(hi), b)
    whole = polar_equation(I(lo, hi), b)
    assert left.v.h < 0 < right.v.l
    assert whole.d.l > 0
    initial = dict(interval=I(lo, hi).enc(), left=left.v.enc(), right=right.v.enc(),
                   derivative=whole.d.enc())
    trace = []
    while hi-lo > F(1, 10**22):
        mid = (lo+hi)/2
        value = polar_equation(I(mid), b).v
        assert value.l > 0 or value.h < 0, ('root refinement unresolved', str(mid), value.enc())
        trace.append(dict(point=str(mid), value=value.enc()))
        if value.h < 0:
            lo = mid
        else:
            hi = mid
    x = I(lo, hi)
    q = polar_equation(x, b)
    assert q.v.l <= 0 <= q.v.h and q.d.l > 0
    return x, dict(initial=initial, bisection=trace, isolated_root=x.enc(),
                   final_left=polar_equation(I(lo), b).v.enc(),
                   final_right=polar_equation(I(hi), b).v.enc(), derivative=q.d.enc(),
                   equation='2*(G4(x)^-1)[1,2]/sqrt((G4(x)^-1)[1,1]*(G4(x)^-1)[2,2])+2*cos(pi/27)=0',
                   proof='All radicands and inverse pivots are nonzero throughout the initial interval. Opposite endpoint signs and a strictly positive derivative prove exactly one root there. Rational bisection encloses that same root. The constants and operations are algebraic, so this defines an exact algebraic normal parameter.')


def geometry():
    z = C.geometry()
    old, pair = z['normals'], z['pair']
    b = U.rnd(pair(old[8], old[30]))
    oldx = U.rnd(pair(old[44], old[52]))
    x, root = isolate(b)
    print('RELEASE34_ROOT_ISOLATED', x.enc(), flush=True)
    zero, two, t, a = I(0), I(2), D.theta(4), D.theta(15)
    gram = [[two, b, zero], [b, two, zero], [zero, zero, two]]
    inverse, pivots = P.C53.inverse(gram)
    frame = [old[i] for i in [8, 30, 44]]
    c0, cc = [zero, -t, oldx], [zero, -t, x]
    p0 = P.combine(U.imv(inverse, c0), frame)
    k = P.sub(old[52], p0)
    k2 = U.rnd(2-U.idot(c0, U.imv(inverse, c0)))
    rho = U.rnd(2-U.idot(cc, U.imv(inverse, cc)))
    assert k2.l > 0 and rho.l > 0
    p = P.combine(U.imv(inverse, cc), frame)
    replacement = [U.rnd(v+U.sqrt(rho/k2)*w) for v, w in zip(p, k)]

    gram3 = [[two, zero, zero], [zero, two, x], [zero, x, two]]
    j3, pivots3 = P.C53.inverse(gram3)
    ccontact, old30 = [zero, -a, -t], [b, zero, -t]
    basis3 = [old[8], old[44], replacement]
    w = P.sub(old[30], P.combine(U.imv(j3, old30), basis3))
    w2 = U.rnd(2-U.idot(old30, U.imv(j3, old30)))
    crho = U.rnd(2-U.idot(ccontact, U.imv(j3, ccontact)))
    assert w2.l > 0 and crho.l > 0
    cp = P.combine(U.imv(j3, ccontact), basis3)
    c = [U.rnd(v-U.sqrt(crho/w2)*ww) for v, ww in zip(cp, w)]
    c30 = U.rnd(U.idot(ccontact, U.imv(j3, old30))-U.sqrt(crho*w2))
    g4 = [[two, zero, -t, c30], [zero, two, x, -a],
          [-t, x, two, -t], [c30, -a, -t, two]]
    j4, pivots4 = P.C53.inverse(g4)
    assert j4[1][1].l > 0 and j4[2][2].l > 0
    basis4 = [old[30], old[44], replacement, c]
    aa = [U.rnd(-U.sqrt(2/j4[2][2])*v) for v in P.combine(j4[2], basis4)]
    bb = [U.rnd(-U.sqrt(2/j4[1][1])*v) for v in P.combine(j4[1], basis4)]
    normals = {i:n for i, n in enumerate(old) if i not in [34, 52]}
    normals.update({52:replacement, 63:aa, 64:bb, 65:c})
    z.update(normals=normals, root_certificate=root,
             parameter_data=dict(b=b.enc(), old44_52=oldx.enc(), new44_52=x.enc(),
                 replacement_gram=[P.enc(row) for row in gram], replacement_inverse_pivots=pivots,
                 replacement_kernel=P.enc(k), replacement_kernel_square=k2.enc(), replacement_residual=rho.enc(),
                 C_gram=[P.enc(row) for row in gram3], C_inverse_pivots=pivots3,
                 C_kernel_square=w2.enc(), C_residual=crho.enc(), C30=c30.enc(),
                 polar_gram=[P.enc(row) for row in g4], polar_inverse_pivots=pivots4,
                 polar_inverse=[P.enc(row) for row in j4]),
             common_center=C.C.decode(C.C.read(C.TAG)['common_center']))
    return z


def run():
    z = geometry()
    normals, pair, h = z['normals'], z['pair'], z['common_center']
    labels = sorted(normals)
    prescribed = {
        52:{**{i:2 for i in S+[8]}, 30:4},
        63:{**{i:2 for i in S+[30, 44, 65]}, 64:27},
        64:{**{i:2 for i in S+[30, 52, 65]}, 63:27},
        65:{**{i:2 for i in S+[8, 63, 64]}, 44:15, 52:4},
    }
    old_labels = D.E.labels()

    def literal(i, j):
        if i == j:
            return I(2)
        m = prescribed.get(i, {}).get(j, prescribed.get(j, {}).get(i))
        if m is not None:
            return -D.theta(m)
        if max(i, j) < 55 and old_labels[i][j]:
            return -D.theta(old_labels[i][j])
        return pair(normals[i], normals[j])

    gates, counts = [], Counter()
    for i, j in combinations(labels, 2):
        if i not in NEW and j not in NEW:
            continue
        value = U.rnd(pair(normals[i], normals[j]))
        m = prescribed.get(i, {}).get(j, prescribed.get(j, {}).get(i))
        if m is None:
            assert value.h < -2, ('new dotted gate', i, j, value.enc())
            label = 'dotted'
        else:
            residual = U.rnd(value+D.theta(m))
            assert residual.l <= 0 <= residual.h, ('exact-identity consistency', i, j)
            label = m
        counts[str(label)] += 1
        gates.append(dict(pair=[i, j], label=label, product=value.enc()))
    assert len(gates) == 250
    hh = U.rnd(pair(h, h))
    hp = {i:U.rnd(pair(h, normals[i])) for i in labels}
    assert hh.h < 0 and all(v.h < 0 for v in hp.values())
    norms = []
    for i in NEW:
        q = U.rnd(pair(normals[i], normals[i]))
        assert q.l <= 2 <= q.h
        norms.append(dict(wall=i, norm=q.enc()))
    facets = []
    for i in labels:
        point = [U.rnd(x-hp[i]*v/2) for x, v in zip(h, normals[i])]
        square = U.rnd(hh-hp[i]*hp[i]/2)
        assert square.h < 0
        # Exact acuteness and h.n_j<0 imply every nonincident product
        # h.n_j-(h.n_i)(n_i.n_j)/2<0, without interval cancellation.
        facets.append(dict(wall=i, point=P.enc(point), norm_and_future=square.enc()))

    corner_sides = [[8, 44, 65], [8, 52, 65], [30, 44, 63], [30, 52, 64],
                    [44, 63, 65], [52, 64, 65], [30, 63, 64], [63, 64, 65]]
    vertices = []
    for sides in corner_sides:
        support = S+sides
        gram = [[literal(i, j) for j in support] for i in support]
        pd = P.pd(gram)
        inverse, pivots = P.C53.inverse(gram)
        hs = [hp[i] for i in support]
        coeff = U.imv(inverse, hs)
        point = P.sub(h, P.combine(coeff, [normals[i] for i in support]))
        square = U.rnd(hh-U.idot(hs, coeff))
        assert square.h < 0
        products = []
        for j in labels:
            value = U.rnd(hp[j]-U.idot(coeff, [literal(i, j) for i in support]))
            if j in support:
                assert value.l <= 0 <= value.h
            else:
                assert value.h < 0, ('corner gate', sides, j, value.enc())
            products.append(dict(wall=j, product=value.enc()))
        vertices.append(dict(sides=sides, support=support, point=P.enc(point),
                             PD_pivots=pd, inverse_pivots=pivots, norm_and_future=square.enc(), products=products))

    faces = []
    for cap, cycle in [(63, [30, 44, 65, 64]), (64, [30, 52, 65, 63]), (65, [8, 44, 63, 64, 52])]:
        support = S+[cap]
        gram = [[literal(i, j) for j in support] for i in support]
        pd = P.pd(gram)
        inverse, pivots = P.C53.inverse(gram)
        projections = []
        for j in labels:
            if j in support:
                continue
            contacts = [literal(i, j) for i in support]
            q = U.rnd(2-U.idot(contacts, U.imv(inverse, contacts)))
            assert q.l > 0 if j in cycle else q.h < 0, ('face side classification', cap, j, q.enc())
            projections.append(dict(wall=j, projected_square=q.enc()))
        for a, b in zip(cycle, cycle[1:]+cycle[:1]):
            assert any(set(v['sides']) == {cap, a, b} for v in vertices)
        faces.append(dict(support=support, cyclic_sides=cycle, vertex_count=len(cycle),
                          PD_pivots=pd, inverse_pivots=pivots, all_other_wall_projection_squares=projections))
    oldcone = C.C.read(OLD_CONE)
    assert oldcone['computational_root_basis'] == z['basis_ids']
    retained_vertices = []
    deleted_vertices = []
    for row in oldcone['final_rays']:
        if row['kind'] != 'finite_vertices':
            continue
        if set(row['active']) & {34, 52}:
            deleted_vertices.append(row['active'])
            continue
        point = C.C.decode(row['normalized_computational_coordinate_bounds'])
        products = [U.rnd(pair(point, normals[i])) for i in NEW]
        assert all(q.h < 0 for q in products), ('old finite vertex', row['active'])
        retained_vertices.append(dict(support=row['active'], new_products=P.enc(products)))
    assert len(retained_vertices) == 188 and len(deleted_vertices) == 5
    corner_supports = {tuple(sorted(v['support'])) for v in vertices}
    found, census = set(), []
    for wall in NEW:
        neighbors = sorted(j for j in labels if j != wall and
                           prescribed.get(wall, {}).get(j, prescribed.get(j, {}).get(wall)) is not None)
        records = []
        for chosen in combinations(neighbors, 8):
            support = tuple(sorted((wall,)+chosen))
            blocker = None
            for i, j in combinations(support, 2):
                m = prescribed.get(i, {}).get(j, prescribed.get(j, {}).get(i))
                if m is None and max(i, j) < 55 and old_labels[i][j]:
                    m = old_labels[i][j]
                if m is None:
                    value = literal(i, j)
                    assert value.h < -2
                    blocker = dict(pair=[i, j], product=value.enc())
                    break
            if blocker:
                records.append(dict(support=list(support), status='DOTTED_PAIR', blocker=blocker))
            else:
                assert support in corner_supports, ('new unlisted all-finite candidate', wall, support)
                found.add(support)
                records.append(dict(support=list(support), status='CERTIFIED_FINITE_CORNER'))
        census.append(dict(wall=wall, complete_finite_neighbors=neighbors, support_tests=records))
    assert found == corner_supports
    escape = C.C.decode(oldcone['strict_escape']['computational_point'])
    ep = [dict(wall=i, product=U.rnd(pair(escape, normals[i])).enc()) for i in labels]
    assert pair(escape, escape).l > 0 and pair(escape, h).h < 0
    assert all(F(row['product'][1]) < 0 for row in ep), ('escape changed', ep)
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json',
             f'scripts/{D.TAG}.py', f'results/{OLD_CONE}.json']
    return dict(status='PASS_EXACT_RELEASE34_65_WALL_CHAMBER_WITH_CLOSED_THREE_CAP_PATCH',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                ambient_dimension=9, deleted_old_walls=[34, 52], retained_old_wall_count=61,
                total_wall_count=len(labels), wall_labels=labels, new_labels=NEW,
                computational_root_basis=z['basis_ids'], computational_gram=[P.enc(row) for row in z['exact_gram']],
                root_certificate=z['root_certificate'], parameter_data=z['parameter_data'],
                new_normals=[dict(index=i, computational=P.enc(normals[i]), physical=P.enc(P.combine(normals[i], z['basis']))) for i in NEW],
                norm_consistency=norms, new_contact_checks=gates, new_contact_counts=dict(counts),
                common_center=P.enc(h), common_center_norm=hh.enc(), common_center_products=[dict(wall=i, product=hp[i].enc()) for i in labels],
                facet_witnesses=facets, finite_corners=vertices, complete_compact_two_faces=faces,
                preserved_old_finite_vertices=retained_vertices, deleted_old_finite_supports=deleted_vertices,
                new_wall_finite_vertex_census=census, complete_finite_vertex_count=196,
                strict_escape=dict(computational_point=P.enc(escape), norm=pair(escape, escape).enc(),
                                   future=pair(escape, h).enc(), products=ep),
                exact_equalities='Replacement52 is the norm-two completion of its [8,30,44] contacts, in the same perpendicular direction as old52. C is the negative normalized completion of its [8,44,52] contacts. A and B are negative normalized duals to52 and44 in [30,44,52,C]. All roots are right to S6. These projection identities prove their norms and prescribed angles; the unique algebraic root x supplies A.B=-2cos(pi/27). Interval overlap checks only consistency.',
                face_proof='For each S6+new-cap support the Gram is positive definite. Exactly its listed cyclic side normals have positive projected square; every other retained normal projects to a timelike vector whose halfspace contains the future plane, by the projected strict center. The listed cyclic adjacent corners are finite and feasible with all other walls strict. Their convex polygon is the entire side-halfspace intersection, hence a whole compact quadrilateral or pentagon.',
                finite_census_proof='Every finite vertex incident to a new root needs eight finite neighbors; the exhaustive subsets either contain a certified dotted pair or equal one of the eight verified new corners. A finite vertex incident only to retained old roots already occurred in the old acute chamber: projecting its strict center off the positive-definite nine-root support makes every other old wall strict by the nonnegative inverse-Gram property. Thus only the188 old corners avoiding34 and52 can survive without a new root. All four new products are strictly negative at each. The complete finite-vertex count is188+8=196.',
                scope='Exact alternative after explicitly deleting34 and replacing52. Three entire two-faces close the8/30/44/replacement52 patch. The complete finite-vertex count is196, but no complete ambient exterior cone or persistence of faces involving deleted walls is asserted. The strict escape proves infinite volume; this is not a new finite-volume H9 example.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['new_contact_counts'], flush=True)
