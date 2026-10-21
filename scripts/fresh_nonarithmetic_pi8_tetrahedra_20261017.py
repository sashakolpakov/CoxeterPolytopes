#!/usr/bin/env python3
"""Exact continuous realization screen of all pi5 tetrahedra at the m=8 cap.

Each template has four O(1,1) branches with positive parameter. Necessary
native acuteness/common-interior inequalities are quadratic in that parameter.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_h3_edge730_domain_20261017 as D

P, U, I = D.P, D.U, D.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_tetrahedra_20261017'
CAT = 'fresh_h9_pi5_tetrahedron_templates_20261017'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def run():
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
    records = []
    for template in catalogue['records']:
        labels = template['labels']
        assert labels[0] == 5
        pu = P.combine(U.imv(wedge_inverse, [-D.theta(labels[1]), -D.theta(labels[3])]), wedge)
        pv = P.combine(U.imv(wedge_inverse, [-D.theta(labels[2]), -D.theta(labels[4])]), wedge)
        ru, rv = U.rnd(2-pair(pu, pu)), U.rnd(2-pair(pv, pv))
        d = U.rnd(-D.theta(labels[5])-pair(pu, pv))
        assert ru.l > 0 and rv.l > 0 and d.h < 0
        gamma2 = U.rnd(d*d/ru-rv)
        assert gamma2.l > 0
        gamma, uscale = U.sqrt(gamma2), U.sqrt(ru/2)
        dscale = U.rnd(d/U.sqrt(2*ru))
        branches = []
        for side in (-1, 1):
            aa = [U.rnd((x+side*y)/sq2) for x, y in zip(k0, k1)]
            bb = [U.rnd((-x+side*y)/sq2) for x, y in zip(k0, k1)]
            for sign in (-1, 1):
                ua, ub = [U.rnd(uscale*x) for x in aa], [U.rnd(uscale*x) for x in bb]
                va = [U.rnd((dscale+sign*gamma/sq2)*x) for x in aa]
                vb = [U.rnd((dscale-sign*gamma/sq2)*x) for x in bb]
                quadratics = []
                for name, p, x, y in [('u', pu, ua, ub), ('v', pv, va, vb)]:
                    for wall in list(i for i in range(56) if i not in fixed)+['V']:
                        n = vertex if wall == 'V' else normals[wall]
                        coefficients = [pair(x, n), pair(p, n), pair(y, n)]
                        quadratics.append(dict(new_root=name, old_constraint=wall, values=coefficients))
                lower, upper, positive = [], [], []
                for row in quadratics:
                    x, y, z = row['values']
                    if all(v.l >= 0 for v in (x, y, z)) and any(v.l > 0 for v in (x, y, z)):
                        positive.append(row)
                    if x.h < 0 < z.l:
                        root = U.rnd((-y-U.sqrt(y*y-4*x*z))/(2*x))
                        assert root.l > 0
                        lower.append((root, row))
                    elif x.l > 0 > z.h:
                        root = U.rnd((-y+U.sqrt(y*y-4*x*z))/(2*x))
                        assert root.l > 0
                        upper.append((root, row))
                def encode(row):
                    return dict(new_root=row['new_root'], old_constraint=row['old_constraint'],
                                quadratic=P.enc(row['values']))
                if positive:
                    proof = dict(status='POSITIVE_QUADRATIC', constraint=encode(positive[0]))
                elif lower and upper:
                    lo, lowrow = max(lower, key=lambda x: x[0].l)
                    hi, highrow = min(upper, key=lambda x: x[0].h)
                    if lo.l > hi.h:
                        proof = dict(status='DISJOINT_NECESSARY_PARAMETER_BOUNDS',
                                     lower=lo.enc(), upper=hi.enc(),
                                     lower_constraint=encode(lowrow), upper_constraint=encode(highrow))
                    else:
                        proof = dict(status='UNRESOLVED_BY_TWO_BOUND_TEST', lower=lo.enc(), upper=hi.enc())
                else:
                    proof = dict(status='UNRESOLVED_BY_TWO_BOUND_TEST')
                if proof['status'].startswith('UNRESOLVED'):
                    proof['all_necessary_quadratics'] = [encode(row) for row in quadratics]
                    proof['new_root_parameter_vectors'] = dict(u=[P.enc(z) for z in (pu, ua, ub)],
                                                               v=[P.enc(z) for z in (pv, va, vb)])
                branches.append(dict(kernel_side=side, timelike_sign=sign, **proof))
        records.append(dict(template_index=template['index'], labels=labels,
                            positive_residuals=P.enc([ru, rv]), residual_product=d.enc(),
                            gamma_square=gamma2.enc(), branches=branches))
        print('TEMPLATE', template['index'], [r['status'] for r in branches], flush=True)
    unresolved = sum(r['status'].startswith('UNRESOLVED') for row in records for r in row['branches'])
    paths = [f'scripts/{TAG}.py', f'scripts/{D.TAG}.py', f'results/{D.TAG}.json',
             f'scripts/{CAT}.py', f'results/{CAT}.json', f'scripts/{P.TAG}.py',
             f'results/{P.TAG}.json', P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55]
    return dict(status='PASS_COMPLETE_PI8_TETRAHEDRON_EXCLUSION' if not unresolved else 'PARTIAL_PI8_TETRAHEDRON_SCREEN',
                sources={p: sha(p) for p in paths}, initial_old_walls=55,
                new_endpoint_label_to51=8, new_endpoint_normal=P.enc(cap),
                computational_root_basis=basis_ids, computational_gram=[P.enc(row) for row in exact_gram],
                computational_inverse_pivots=full_inverse_pivots,
                new_H3_support=white, fixed_wedge=[7, 30], full_kernel_support=fixed,
                support_PD_pivots=pd_pivots, inverse_pivots=pivots,
                kernel_future_unit=P.enc(k0), kernel_positive_unit=P.enc(k1),
                kernel_future_seed=coordinate0, kernel_coordinate_seed=coordinate, common_interior_vertex=P.enc(vertex),
                template_count=18, orientation_branches=72, unresolved_branches=unresolved,
                records=records,
                parameterization='With A=(k0+s*k1)/sqrt2,B=(-k0+s*k1)/sqrt2, u=pu+sqrt(ru/2)*(t*A+B/t). Let d=g_uv-pu.pv and gamma=sqrt(d²/ru-rv); v=pv+(d/sqrt(2ru)+e*gamma/sqrt2)*t*A+(d/sqrt(2ru)-e*gamma/sqrt2)*B/t. Both signs s,e and all t>0 are included.',
                obstruction='Each saved native or persistent-vertex inequality is a*t²+b*t+c<=0. Opposite signs of a,c give exactly one positive root, hence the indicated lower or upper bound. Disjoint certified bounds exclude all t, before any additional finite-angle equation.',
                scope='Complete compact Coxeter tetrahedron templates with old7,30 as their pi5 edge, for the particular coherent endpoint cap m=8. No larger polyhedron, other m>=9, or strict-dotted endpoint parameter is excluded.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
        print('PASS saved continuous tetrahedron replay')
    else:
        out.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['unresolved_branches'])
