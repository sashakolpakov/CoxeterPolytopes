#!/usr/bin/env python3
"""All seven actual S6-facet profiles miss the strict 55-wall escape.

Every matrix entry denotes an exact algebraic number. The frozen physical
normal recipes and their outward rational intervals are reused; no floating
point value is used for an acceptance or rejection in this certificate.
"""
import argparse
import itertools
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_nonright_s6_20261016_cone as CN
import fresh_nonarithmetic_escape52_20261016_cap as C53

U = C53.U
B, D, A, I = U.B, U.D, U.A, U.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_a2_escape_partition_20261017'
CAP52 = 'results/fresh_nonarithmetic_escape52_20261016_cap.json'
CAPS55 = 'results/fresh_nonarithmetic_nonright_s6_20261016_caps.json'
CONE55 = 'results/fresh_nonarithmetic_nonright_s6_20261016_cone.json'
S = [0, 1, 2, 3, 5, 6]
FRAME = [7, 10, 30, 47]


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def enc(v):
    return [x.enc() for x in v]


def sub(v, w):
    return [U.rnd(x - y) for x, y in zip(v, w)]


def combine(coefficients, vectors):
    return [U.idot(coefficients, [v[j] for v in vectors]) for j in range(10)]


def pd(matrix):
    work = [list(row) for row in matrix]
    pivots = []
    for k in range(len(work)):
        pivot = work[k][k]
        assert pivot.l > 0
        pivots.append(pivot.enc())
        for i in range(k + 1, len(work)):
            for j in range(i, len(work)):
                work[i][j] = work[j][i] = U.rnd(
                    work[i][j] - work[i][k] * work[j][k] / pivot)
    return pivots


def physical_data():
    old = B.read(U.CAP)
    cluster = B.read(C53.CL)
    cap52 = B.read(CAP52)
    caps55 = B.read(CAPS55)
    cone55 = B.read(CONE55)
    q = [[A.B.db(x) for x in row] for row in old['form']]
    qi = [[U.iv(D.lift(x)) for x in row] for row in q]
    raw = [tuple(U.dec(x) for x in row) for row in old['normals']]
    normals = []
    normal_scales = []
    for n in raw:
        scale = U.sqrt(U.iv(D.pair(n, n, q)) / 2)
        normal_scales.append(scale.enc())
        normals.append([U.rnd(U.iv(x) / scale) for x in n])
    normals += [[I(*x) for x in row['normal_interval']]
                for row in cluster['caps']]
    normals += [[I(*x) for x in cap52['normal_interval']]]
    normals += [[I(*x) for x in row['normal_interval']]
                for row in caps55['caps']]
    assert len(normals) == 55
    h = [U.iv(U.dec(x)) for x in old['center']]
    e = [CN.ii(CN.decode(x)) for x in cone55['strict_escape']['point']]
    e_scale = CN.power_scale(e)
    e = [U.rnd(e_scale * x) for x in e]
    return q, qi, normals, h, e, e_scale, normal_scales, cone55


def run():
    q, qi, normals, h, e, e_scale, normal_scales, cone = physical_data()
    pair = lambda v, w: U.idot(v, U.imv(qi, w))
    assert pair(h, h).h < 0 < pair(e, e).l
    assert pair(h, e).h < 0
    assert all(pair(e, n).h < 0 for n in normals)
    sg = [[pair(normals[i], normals[j]) for j in S] for i in S]
    support_pivots = pd(sg)
    sj, support_inverse_pivots = C53.inverse(sg)

    def project(v):
        return sub(v, combine(U.imv(sj, [pair(normals[i], v) for i in S]),
                              [normals[i] for i in S]))

    frame = [project(normals[i]) for i in FRAME]
    gram = [[pair(v, w) for w in frame] for v in frame]
    inverse, inverse_pivots = C53.inverse(gram)

    # Coordinates are xi_i=<z,proj(n_FRAME[i])>, so q(z)=xi^T G^-1 xi.
    def row(v):
        return U.imv(inverse, [pair(w, v) for w in frame])

    escape_row = row(e)
    r53 = row(normals[53])
    # The saved exact Coxeter inventory, together with positive rescaling,
    # identifies the literal S6 Gram as H4 + A2.
    inventory = B.read('results/fresh_h9_persistent_cusp_screen_20261016.json')
    assert [(i, j, inventory['pair_matrix'][i][j])
            for i, j in itertools.combinations(S, 2)
            if inventory['pair_matrix'][i][j] != 2] == [
                (0, 1, 5), (1, 2, 3), (2, 3, 3), (5, 6, 3)]
    h4i, h4_pivots = C53.inverse([r[:4] for r in sg[:4]])
    assert all(x.l > 0 for r in h4i for x in r)
    assert all(h4i[i][i].l > 2 for i in range(4))
    phi = U.theta(5)
    threshold = U.rnd(2 / (phi * phi) - I(F(2, 3)))
    assert threshold.l > F(9, 100)

    # Zero profile: E_perp is a positive combination of three old projected
    # normals and one actual finite vertex. The latter persists under every
    # acute completion with a common strict future interior, not just H.
    vertex_support = S + [8, 10, 44]
    assert vertex_support in cone['complete_census']['finite_vertices']
    vg = [[pair(normals[i], normals[j]) for j in vertex_support]
          for i in vertex_support]
    vertex_pd = pd(vg)
    vj, vertex_inverse_pivots = C53.inverse(vg)
    v = sub(h, combine(U.imv(vj, [pair(h, normals[i]) for i in vertex_support]),
                       [normals[i] for i in vertex_support]))
    vnorm, vfuture = pair(v, v), pair(v, h)
    assert vnorm.h < 0 and vfuture.h < 0
    vertex_signs = []
    for i, n in enumerate(normals):
        if i not in vertex_support:
            product = pair(v, n)
            assert product.h < 0
            vertex_signs.append(dict(wall=i, product=product.enc()))
    zero_rows = [row(normals[i]) for i in [7, 47, 51]] + [row(v)]
    zi, zero_inverse_pivots = C53.inverse([list(r) for r in zip(*zero_rows)])
    zw = U.imv(zi, escape_row)
    assert all(x.l > 0 for x in zw)
    zero = dict(old_walls=[7, 47, 51], vertex_support=vertex_support,
                vertex=enc(v), vertex_norm=vnorm.enc(), vertex_future=vfuture.enc(),
                vertex_pd_pivots=vertex_pd, vertex_inverse_pivots=vertex_inverse_pivots,
                vertex_nonincident_products=vertex_signs,
                positive_coefficients=enc(zw), inverse_pivots=zero_inverse_pivots)

    # Node6: r=c(p6+z), c>0. Equality E_perp=sum lambda_i n_i,perp
    # and acuteness give <p6+z,E> <= <p6,E>-sum lambda_i<p6,n_i><0.
    p6 = combine([I(F(-1, 3)), I(F(-2, 3))], [normals[5], normals[6]])
    node6_walls = [7, 30, 47, 52]
    n6rows = [row(normals[i]) for i in node6_walls]
    n6i, n6pivots = C53.inverse([list(r) for r in zip(*n6rows)])
    n6w = U.imv(n6i, escape_row)
    assert all(x.l > 0 for x in n6w)
    n6bound = U.rnd(pair(p6, e) - U.idot(n6w, [pair(p6, normals[i])
                                                            for i in node6_walls]))
    assert n6bound.h < 0
    node6 = dict(old_walls=node6_walls, particular=enc(p6),
                 positive_coefficients=enc(n6w), inverse_pivots=n6pivots,
                 escape_upper_bound_per_contact=n6bound.enc())

    # Node5: a six-halfspace relaxation suffices. The first four coefficient
    # rows are exact unit rows by definition of xi, not interval zeros.
    p5 = combine([I(F(-2, 3)), I(F(-1, 3))], [normals[5], normals[6]])
    rows = [[I(int(i == j)) for j in range(4)] for i in range(4)]
    rows += [r53, [-x for x in escape_row]]
    rhs = [-pair(p5, normals[i]) for i in FRAME + [53]] + [pair(p5, e)]
    # Strictly positive dependence of all six rows proves the recession cone
    # is {0}, because the first four rows span the whole dual space.
    t = F(1, 32768)
    recession_weights = [U.rnd(a - t * b) for a, b in zip(escape_row, r53)]
    assert all(x.l > 0 for x in recession_weights)
    recession_weights += [I(t), I(1)]
    bases, vertices = [], []
    for active in itertools.combinations(range(6), 4):
        mi, pivots = C53.inverse([rows[i] for i in active])
        point = U.imv(mi, [rhs[i] for i in active])
        slacks = [U.rnd(U.idot(r, point) - b) for r, b in zip(rows, rhs)]
        outside = [i for i in range(6) if i not in active and slacks[i].l > 0]
        record = dict(active=list(active), inverse_pivots=pivots,
                      point=enc(point), inactive_slacks=[dict(row=i, slack=slacks[i].enc())
                                                      for i in range(6) if i not in active])
        if outside:
            record.update(feasible=False, positive_obstruction_rows=outside)
        else:
            assert all(slacks[i].h < 0 for i in range(6) if i not in active)
            norm = U.idot(point, U.imv(inverse, point))
            assert norm.h < F(1, 50)
            record.update(feasible=True, vertex_index=len(vertices), norm=norm.enc())
            vertices.append((active, point, norm))
        bases.append(record)
    assert len(bases) == 15 and len(vertices) == 9
    barycenter = [U.rnd(sum((v[j] for _, v, _ in vertices), I()) / 9) for j in range(4)]
    interior_slacks = [U.rnd(U.idot(r, barycenter) - b) for r, b in zip(rows, rhs)]
    assert all(x.h < 0 for x in interior_slacks)
    edges = []
    for ai, bi in itertools.combinations(range(9), 2):
        aa, va, qa = vertices[ai]
        ab, vb, qb = vertices[bi]
        if len(set(aa) & set(ab)) != 3:
            continue
        direction = sub(vb, va)
        qd = U.idot(direction, U.imv(inverse, direction))
        derivative_half = U.idot(va, U.imv(inverse, direction))
        edge = dict(vertices=[ai, bi], common_active=sorted(set(aa) & set(ab)),
                    direction_norm=qd.enc(), initial_half_derivative=derivative_half.enc())
        if qd.l > 0:
            edge['maximum_location'] = 'endpoint_by_strict_convexity'
        else:
            assert qd.h < 0
            critical = U.rnd(-derivative_half / qd)
            edge['critical_parameter'] = critical.enc()
            if critical.h < 0 or critical.l > 1:
                edge['maximum_location'] = 'endpoint_by_monotonicity'
            else:
                maximum = U.rnd(qa - derivative_half * derivative_half / qd)
                assert maximum.h < F(1, 50)
                edge.update(maximum_location='interior_or_endpoint',
                            stationary_norm_upper_bound=maximum.enc())
        edges.append(edge)
    assert len(edges) == 18
    assert all(sum(i in edge['vertices'] for edge in edges) == 4 for i in range(9))
    assert sum(edge['maximum_location'] == 'endpoint_by_strict_convexity'
               for edge in edges) == 17
    node5 = dict(particular=enc(p5), inequality_walls=FRAME + [53, 'negative_escape'],
                 coefficient_rows=[enc(r) for r in rows], rhs=enc(rhs),
                 positive_recession_dependence=enc(recession_weights),
                 recession_parameter=str(t), all_four_row_bases=bases,
                 strict_interior_point=enc(barycenter), strict_interior_slacks=enc(interior_slacks),
                 vertex_count=9, edges=edges, edge_count=18,
                 certified_norm_upper_bound='1/50', required_norm_lower_bound='9/100',
                 exact_smallest_required_norm='2/phi^2-2/3',
                 smallest_required_norm_interval=threshold.enc())
    paths = [f'scripts/{TAG}.py', f'scripts/{CN.TAG}.py', f'scripts/{C53.TAG}.py',
             f'scripts/{U.TAG}.py', U.CAP, C53.CL, CAP52, CAPS55, CONE55,
             'results/fresh_h9_persistent_cusp_screen_20261016.json',
             'results/fresh_nonarithmetic_nonright_s6_20261016_caps_review.json',
             'results/fresh_nonarithmetic_nonright_s6_20261016_cone_review.json',
             'results/fresh_nonarithmetic_h3_a2_profiles_20261016.json']
    for path in paths:
        if path.endswith('.json'):
            B.read(path)
    return dict(status='PASS_EXACT_ALL_SEVEN_ACTUAL_S6_FACET_PROFILES_MISS_E55',
                sources={path: sha(path) for path in paths}, retained_walls=55, support=S,
                escape_source=CONE55, escape_positive_rational_scale=str(e_scale),
                scaled_escape=enc(e), native_to_norm_two_divisors=normal_scales,
                support_gram=[enc(r) for r in sg], support_pd_pivots=support_pivots,
                support_inverse_pivots=support_inverse_pivots,
                H4_inverse=[enc(r) for r in h4i], H4_inverse_pivots=h4_pivots,
                profile_count=7, nonzero_profiles=[[node, m] for node in [5, 6] for m in [3, 4, 5]],
                residual_frame=FRAME, residual_gram=[enc(r) for r in gram],
                residual_inverse=[enc(r) for r in inverse], residual_inverse_pivots=inverse_pivots,
                escape_covector=enc(escape_row), zero_profile=zero, node6_profiles=node6,
                node5_profiles=node5,
                arithmetic='Exact literal algebraic recipes, with outward 256-bit rational intervals. Shared frozen physical/sign arithmetic; this is a primary certificate, not an independent review.',
                conclusion='Every additional Coxeter wall inducing an actual facet of the S6 face in an acute completion with a common strict future interior has strictly negative pairing with the retained55 strict spacelike escape. All finite labels and all real coefficient fields are covered through the complete seven-profile projection lemma.',
                scope='This is a fixed-face escape obstruction. It does not exclude a finite-volume completion using further ambient walls which do not induce facets of S6; it does not assert that closing S6 alone is impossible. No search-sampling conclusion is used.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
        print('PASS saved all-seven-profile escape partition replay')
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
