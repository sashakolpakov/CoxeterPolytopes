#!/usr/bin/env python3
"""Three coherent added walls and an actual compact H3 tetrahedral face.

The final two walls realize literal pi5 template9 at the first m=8 endpoint.
All signs are certified by outward rational interval arithmetic; prescribed
equalities follow from the displayed projection and quadratic recipes.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_h3_edge730_domain_20261017 as D

P, U, I = D.P, D.U, D.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_tetrahedron_completion_20261017'
CAT = 'fresh_h9_pi5_tetrahedron_templates_20261017'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def geometry():
    endpoint = D.E.B.read(f'results/{D.TAG}.json')
    catalogue = D.E.B.read(f'results/{CAT}.json')
    assert catalogue['literal_templates'] == 18
    _, qi, normals, h, _, _, _, _ = P.physical_data()
    physical_pair = lambda v, w: U.idot(v, U.imv(qi, w))
    a = [I(*x) for x in endpoint['A']]
    b = [I(*x) for x in endpoint['B']]
    t8 = I(*endpoint['t8'])
    cap = [U.rnd(t8*x+y/t8) for x, y in zip(a, b)]
    normals.append(cap)
    physical_vertex = [I(*x) for x in endpoint['common_interior_vertex']]
    white, fixed = [0, 1, 2, 3, 51, 55], [0, 1, 2, 3, 51, 55, 7, 30]
    # A literal-root basis avoids the severe cancellation of a globally
    # boosted physical kernel. Its complete Gram has only finite contacts.
    basis_ids = fixed+[5, 6]
    basis = [normals[i] for i in basis_ids]
    old_labels = D.E.labels()
    def basis_label(i, j):
        if 55 in (i, j):
            return 8 if 51 in (i, j) else 2
        return old_labels[i][j]
    exact_gram = [[I(2) if i == j else -D.theta(basis_label(i, j)) for j in basis_ids] for i in basis_ids]
    full_inverse, full_inverse_pivots = P.C53.inverse(exact_gram)
    def coordinates(x):
        return U.imv(full_inverse, [physical_pair(b, x) for b in basis])
    old_physical = normals
    normals = [([I(int(j == basis_ids.index(i))) for j in range(10)] if i in basis_ids else coordinates(n))
               for i, n in enumerate(old_physical)]
    h, vertex = coordinates(h), coordinates(physical_vertex)
    pair = lambda v, w: U.idot(v, U.imv(exact_gram, w))
    frame = [normals[i] for i in fixed]
    gram = [[pair(x, y) for y in frame] for x in frame]
    pd_pivots = P.pd(gram)
    inverse, pivots = P.C53.inverse(gram)

    def project(v):
        return P.sub(v, P.combine(U.imv(inverse, [pair(v, n) for n in frame]), frame))

    timelike_choices = []
    for coordinate0 in (8, 9):
        seed = [I(int(i == coordinate0)) for i in range(10)]
        trial = project(seed)
        trial_square = pair(trial, trial)
        if trial_square.h < 0:
            timelike_choices.append((coordinate0, trial, trial_square))
    if timelike_choices:
        coordinate0, k0, q0 = timelike_choices[0]
        if pair(k0, h).l > 0:
            k0 = [-x for x in k0]
        assert pair(k0, h).h < 0
    else:
        coordinate0, k0 = 'old_center', project(h)
        q0 = pair(k0, k0)
        assert q0.h < 0
    scale = U.sqrt(-q0)
    k0 = [U.rnd(x/scale) for x in k0]
    for coordinate in range(10):
        e = [I(int(i == coordinate)) for i in range(10)]
        v = project(e)
        c = pair(v, k0)
        k1 = [U.rnd(x+c*y) for x, y in zip(v, k0)]
        q1 = pair(k1, k1)
        if q1.l > 0:
            break
    else:
        raise AssertionError('no positive kernel coordinate')
    scale = U.sqrt(q1)
    k1 = [U.rnd(x/scale) for x in k1]
    wedge = [normals[7], normals[30]]
    wedge_inverse, _ = P.C53.inverse([[pair(x, y) for y in wedge] for x in wedge])
    sq2 = U.sqrt(I(2))
    return locals()


def run():
    z = geometry()
    pair, normals = z['pair'], z['normals']
    wedge, wi = z['wedge'], z['wedge_inverse']
    k0, k1, sq2 = z['k0'], z['k1'], z['sq2']
    template = z['catalogue']['records'][9]
    labels = template['labels']
    assert labels == [5, 2, 3, 3, 2, 4]
    pu = P.combine(U.imv(wi, [I(0), I(-1)]), wedge)
    pv = P.combine(U.imv(wi, [I(-1), I(0)]), wedge)
    ru, rv = U.rnd(2-pair(pu, pu)), U.rnd(2-pair(pv, pv))
    d = U.rnd(-D.theta(4)-pair(pu, pv))
    gamma2 = U.rnd(d*d/ru-rv)
    assert ru.l > 0 and rv.l > 0 and d.h < 0 and gamma2.l > 0
    gamma, uscale = U.sqrt(gamma2), U.sqrt(ru/2)
    dscale = U.rnd(d/U.sqrt(2*ru))
    aa = [U.rnd((x+y)/sq2) for x, y in zip(k0, k1)]
    bb = [U.rnd((-x+y)/sq2) for x, y in zip(k0, k1)]
    ua, ub = [U.rnd(uscale*x) for x in aa], [U.rnd(uscale*x) for x in bb]
    va = [U.rnd((dscale+gamma/sq2)*x) for x in aa]
    vb = [U.rnd((dscale-gamma/sq2)*x) for x in bb]
    qa, qb, qc = [pair(v, normals[6]) for v in (ua, pu, ub)]
    assert qa.h < 0 < qc.l
    parameter = U.rnd((-qb-U.sqrt(qb*qb-4*qa*qc))/(2*qa))
    assert parameter.l > 0
    u = [U.rnd(p+parameter*a+b/parameter) for p, a, b in zip(pu, ua, ub)]
    v = [U.rnd(p+parameter*a+b/parameter) for p, a, b in zip(pv, va, vb)]
    normals += [u, v]
    assert len(normals) == 58
    white = z['white']
    prescribed = {
        55: {**{i: 2 for i in [0, 1, 2, 3, 5, 6, 7, 30]}, 51: 8},
        56: {**{i: 2 for i in white+[7, 6]}, 30: 3},
        57: {**{i: 2 for i in white+[30]}, 7: 3, 56: 4},
    }
    checks = []
    for i in (55, 56, 57):
        norm = pair(normals[i], normals[i])
        assert norm.l <= 2 <= norm.h
        contacts = []
        for j in range(i):
            value = U.rnd(pair(normals[i], normals[j]))
            if j in prescribed[i]:
                label = prescribed[i][j]
                residual = value+D.theta(label)
                assert residual.l <= 0 <= residual.h
                contacts.append(dict(wall=j, label=label, product=value.enc(),
                                     exact_zero_residual='projection/quadratic identity'))
            else:
                assert value.h < -2, (i, j, value.enc())
                contacts.append(dict(wall=j, label='dotted', product=value.enc()))
        checks.append(dict(wall=i, norm=norm.enc(), contacts=contacts))

    # The same persistent old finite vertex is strict for all three additions.
    old_vertex, old_center = z['vertex'], z['h']
    vertex_products = [U.rnd(pair(old_vertex, normals[i])) for i in (55, 56, 57)]
    assert all(x.h < 0 for x in vertex_products)
    for power in range(64):
        epsilon = I(1)/I(2**power)
        center = [U.rnd(a+epsilon*b) for a, b in zip(old_vertex, old_center)]
        center_products = [U.rnd(pair(center, n)) for n in normals]
        if all(x.h < 0 for x in center_products):
            break
    else:
        raise AssertionError('common center perturbation')
    center_square = pair(center, center)
    center_future = pair(center, old_center)
    assert center_square.h < 0 and center_future.h < 0
    facets = []
    for i, n in enumerate(normals):
        point = [U.rnd(x-center_products[i]*y/2) for x, y in zip(center, n)]
        square, future = pair(point, point), pair(point, center)
        assert square.h < 0 and future.h < 0
        products = [U.rnd(pair(point, other)) for other in normals]
        assert products[i].l <= 0 <= products[i].h
        assert all(x.h < 0 for j, x in enumerate(products) if j != i)
        facets.append(dict(wall=i, point=P.enc(point), norm=square.enc(),
                           future=future.enc(), products=P.enc(products)))

    # All four perpendicular sides lie in the Lorentz4 complement of white.
    # Its dual simplex vertices prove that no other old inequality truncates
    # the tetrahedron. This establishes the entire actual face, not just four
    # local corners or a formal Gram diagram.
    sides = [7, 30, 56, 57]
    simplex = [[I(2), -D.theta(5), I(0), I(-1)],
               [-D.theta(5), I(2), I(-1), I(0)],
               [I(0), I(-1), I(2), -D.theta(4)],
               [I(-1), I(0), -D.theta(4), I(2)]]
    si, sp = P.C53.inverse(simplex)
    vertices = []
    for missing in range(4):
        point = P.combine([-si[j][missing] for j in range(4)], [normals[i] for i in sides])
        incident = white+[sides[j] for j in range(4) if j != missing]
        exact_corner = [[(z['exact_gram'][z['basis_ids'].index(a)][z['basis_ids'].index(b)]
                          if a in white and b in white else
                          I(0) if (a in white) != (b in white) else
                          simplex[sides.index(a)][sides.index(b)])
                         for b in incident] for a in incident]
        pivots = P.pd(exact_corner)
        square, future = pair(point, point), pair(point, center)
        assert square.h < 0 and future.h < 0
        products = [U.rnd(pair(point, n)) for n in normals]
        for j, value in enumerate(products):
            if j in incident:
                assert value.l <= 0 <= value.h
            else:
                assert value.h < 0, (missing, j, value.enc())
        vertices.append(dict(omitted_side=sides[missing], incident=incident,
                             point=P.enc(point), norm=square.enc(), future=future.enc(),
                             corner_PD_pivots=pivots, products=P.enc(products)))
    # Convert the reproducible computational-root coordinates to the same
    # physical covector representation used by the frozen 55-wall chamber.
    physical = lambda vector: P.combine(vector, z['basis'])
    paths = [f'scripts/{TAG}.py', f'scripts/{D.TAG}.py', f'results/{D.TAG}.json',
             f'scripts/{CAT}.py', f'results/{CAT}.json', f'scripts/{P.TAG}.py',
             f'results/{P.TAG}.json', P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55, P.CONE55]
    return dict(status='PASS_EXACT_58_WALL_CHAMBER_WITH_COMPACT_PI8_TETRAHEDRAL_FACE',
                sources={p: sha(p) for p in paths}, old_walls=55, total_walls=58,
                computational_root_basis=z['basis_ids'],
                computational_gram=[P.enc(row) for row in z['exact_gram']],
                kernel_future_unit=P.enc(k0), kernel_positive_unit=P.enc(k1),
                endpoint_parameter=z['t8'].enc(), endpoint_formula='r8=t8*A+B/t8 from frozen endpoint locus',
                tetrahedron_template_index=9, tetrahedron_pair_order=z['catalogue']['pair_order'],
                tetrahedron_labels=labels, new_face_support=white, tetrahedron_sides=sides,
                parameter=parameter.enc(), parameter_equation=P.enc([qa, qb, qc]),
                projection_vectors=dict(u=P.enc(pu), v=P.enc(pv)),
                hyperbola_coefficients=dict(u=[P.enc(ua), P.enc(ub)], v=[P.enc(va), P.enc(vb)]),
                residual_norms=P.enc([ru, rv]), residual_product=d.enc(), gamma_square=gamma2.enc(),
                new_normals=[dict(index=i, computational=P.enc(normals[i]),
                                  physical=P.enc(physical(normals[i]))) for i in (55, 56, 57)],
                exact_recipe='A=(k0+k1)/sqrt2,B=(-k0+k1)/sqrt2; u=pu+sqrt(ru/2)*(t*A+B/t); v=pv+(d/sqrt(2ru)+gamma/sqrt2)*t*A+(d/sqrt(2ru)-gamma/sqrt2)*B/t; t is the unique positive root of u.old6=0, gamma^2=d^2/ru-rv.',
                new_contact_checks=checks, old_finite_vertex_products=P.enc(vertex_products),
                center_perturbation_power=power, common_center=P.enc(center),
                physical_common_center=P.enc(physical(center)), center_square=center_square.enc(),
                center_future=center_future.enc(), center_products=P.enc(center_products),
                facet_witnesses=facets, tetrahedron_gram=[P.enc(row) for row in simplex],
                tetrahedron_inverse_pivots=sp, tetrahedron_vertices=vertices,
                proof='Every new product is a prescribed finite Coxeter product or strictly below -2. The displayed perturbation is a strict future center. Each orthogonal projection of that center to one norm2 wall is future timelike and strict for all other walls, so all58 walls are actual facets. The four dual tetrahedron vertices are future timelike, satisfy all58 inequalities, and have PD rank9 incident sets. Their convex hull is the full four-side simplex in the white-complement, hence the entire actual white face is exactly compact.',
                scope='An actual compact H3 face in a coherent 58-wall H9 chamber. This is the H4+I2(8) face, not the original H4+A2 face. No finite-volume assertion and no complete ambient cone census are made here.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
        print('PASS saved 58-wall compact-face replay')
    else:
        out.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'])
