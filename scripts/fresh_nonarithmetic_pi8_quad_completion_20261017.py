#!/usr/bin/env python3
"""Two coordinated walls close the entire r8/S6 facet as a compact H2 quad."""
import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_pi8_tetrahedron_completion_20261017 as C
import fresh_h9_persistent_cusp_affine_20261016 as AF

P, U, I, D = C.P, C.U, C.I, C.D
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_quad_completion_20261017'
DOMAIN = 'fresh_nonarithmetic_pi8_two_side_domains_20261017'
WHITE = [0, 1, 2, 3, 5, 6, 55]
PRESCRIBED = {
    58: {**{i: 2 for i in WHITE+[49]}, 7: 3},
    59: {**{i: 2 for i in WHITE}, 30: 5, 58: 9},
}


def decode(v):
    return [I(*x) for x in v]


def read(tag):
    data = json.loads((ROOT / f'results/{tag}.json').read_text())
    for path, digest in data['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    return data


def geometry():
    z = C.geometry()
    old, domain = read(C.TAG), read(DOMAIN)
    normals, pair = z['normals'][:], z['pair']
    normals += [decode(row['computational']) for row in old['new_normals'][1:]]
    assert len(normals) == 58
    recipes = []
    for domain_index, profile_index, anchor, label in [(0, 1, 49, 2), (1, 3, 58, 9)]:
        block = domain['domains'][domain_index]
        profile = block['finite_profiles'][profile_index]
        orientation = block['orientations'][0]
        assert orientation['kernel_side'] == -1
        assert profile['labels'] == ([2, 2, 2, 3] if domain_index == 0 else [2, 2, 2, 5])
        projection = decode(profile['projection'])
        rho = I(*profile['positive_kernel_square'])
        scale = U.sqrt(rho/2)
        aa, bb = decode(orientation['A']), decode(orientation['B'])
        qa = U.rnd(scale*pair(aa, normals[anchor]))
        qb = U.rnd(pair(projection, normals[anchor])+D.theta(label))
        qc = U.rnd(scale*pair(bb, normals[anchor]))
        assert qa.h < 0 < qc.l
        discriminant = U.rnd(qb*qb-4*qa*qc)
        assert discriminant.l > 0
        parameter = U.rnd((-qb-U.sqrt(discriminant))/(2*qa))
        assert parameter.l > 0
        normal = [U.rnd(p+scale*(parameter*a+b/parameter)) for p, a, b in zip(projection, aa, bb)]
        endpoint = [U.rnd(parameter*a-b/parameter) for a, b in zip(aa, bb)]
        recipes.append(dict(index=len(normals), domain_index=domain_index,
            profile_index=profile_index, profile_labels=profile['labels'],
            frame=block['full_support'], kernel_side=-1, anchor=anchor, anchor_label=label,
            projection=P.enc(projection), residual_square=rho.enc(), A=P.enc(aa), B=P.enc(bb),
            equation=P.enc([qa, qb, qc]), discriminant=discriminant.enc(),
            parameter=parameter.enc(), endpoint=P.enc(endpoint)))
        normals.append(normal)
    z.update(normals=normals, recipes=recipes, old=old, domain=domain)
    return z


def affine_inventory(old):
    labels = D.E.labels()
    for row in old['new_contact_checks']:
        for line in labels:
            line.append(0)
        labels.append([0]*len(labels[0]))
        assert row['wall'] == len(labels)-1
        for contact in row['contacts']:
            m = contact['label']
            labels[-1][contact['wall']] = labels[contact['wall']][-1] = m if isinstance(m, int) else 0
    for wall, contacts in PRESCRIBED.items():
        assert wall == len(labels)
        for line in labels:
            line.append(0)
        labels.append([0]*(wall+1))
        for j, m in contacts.items():
            labels[wall][j] = labels[j][wall] = m
    queue = [(i,) for i in range(60)]
    seen, counts, affine = set(queue), Counter(), []
    while queue:
        support = queue.pop()
        gram = tuple(tuple(AF.scale(AF.ONE, 2) if i == j else AF.VALUES[labels[i][j]] for j in support) for i in support)
        kind, rank, nullity = AF.classify(gram)
        counts[kind] += 1
        if kind == 'positive_semidefinite':
            affine.append(list(support))
            continue
        if kind != 'positive_definite':
            continue
        assert len(support) <= 9
        for j in range(60):
            if j in support:
                continue
            contacts = [labels[j][i] for i in support]
            if any(m not in AF.VALUES for m in contacts) or all(m == 2 for m in contacts):
                continue
            enlarged = tuple(sorted(support+(j,)))
            if enlarged not in seen:
                seen.add(enlarged)
                queue.append(enlarged)
    assert not affine
    assert dict(counts) == {'positive_definite': 670, 'indefinite': 2}
    return dict(counts=dict(counts), affine_circuits=affine,
        proof='Every irreducible Euclidean Coxeter diagram has only labels2/3/4/6 or a parallel pair. There are no parallel pairs. Every connected affine diagram has a connected positive-definite predecessor. The exact connected-growth traversal therefore covers every possible affine circuit; none occurs. This is only a persistent-cusp screen, not a finite-volume proof.')


def run():
    z = geometry()
    pair, normals = z['pair'], z['normals']
    checks = []
    for i in (58, 59):
        square = pair(normals[i], normals[i])
        assert square.l <= 2 <= square.h
        contacts = []
        for j in range(i):
            value = U.rnd(pair(normals[i], normals[j]))
            if j in PRESCRIBED[i]:
                label = PRESCRIBED[i][j]
                residual = value+D.theta(label)
                assert residual.l <= 0 <= residual.h
                contacts.append(dict(wall=j, label=label, product=value.enc(),
                    equality='prescribed frame projection or displayed quadratic identity'))
            else:
                assert value.h < -2, (i, j, value.enc())
                contacts.append(dict(wall=j, label='dotted', product=value.enc()))
        checks.append(dict(wall=i, norm=square.enc(), contacts=contacts))
    vertex, old_center = z['vertex'], z['h']
    vertex_products = [U.rnd(pair(vertex, normals[i])) for i in (58, 59)]
    assert all(x.h < 0 for x in vertex_products)
    for power in range(64):
        epsilon = I(1)/I(2**power)
        center = [U.rnd(v+epsilon*h) for v, h in zip(vertex, old_center)]
        center_products = [U.rnd(pair(center, n)) for n in normals]
        if all(x.h < 0 for x in center_products):
            break
    else:
        raise AssertionError('common future interior')
    assert pair(center, center).h < 0 and pair(center, old_center).h < 0
    facets = []
    for i, normal in enumerate(normals):
        point = [U.rnd(x-center_products[i]*y/2) for x, y in zip(center, normal)]
        square, future = pair(point, point), pair(point, center)
        assert square.h < 0 and future.h < 0
        products = [U.rnd(pair(point, n)) for n in normals]
        assert products[i].l <= 0 <= products[i].h
        assert all(x.h < 0 for j, x in enumerate(products) if j != i)
        facets.append(dict(wall=i, point=P.enc(point), norm=square.enc(), future=future.enc(), products=P.enc(products)))

    white = [normals[i] for i in WHITE]
    wg = [[pair(a, b) for b in white] for a in white]
    white_pd = P.pd(wg)
    wi, white_inverse_pivots = P.C53.inverse(wg)
    def project(v):
        return P.sub(v, P.combine(U.imv(wi, [pair(v, n) for n in white]), white))
    sides = [7, 30, 59, 58]
    projected = {i: project(normals[i]) for i in sides}
    projected_gram = [[U.rnd(pair(projected[i], projected[j])) for j in sides] for i in sides]
    assert all(x.l <= q <= x.h for x, q in zip([projected_gram[i][i] for i in range(4)], [Fraction(4, 3), 2, 2, 2]))
    basis_ids = [7, 30, 58]
    basis = [projected[i] for i in basis_ids]
    bg = [[pair(a, b) for b in basis] for a in basis]
    bi, basis_inverse_pivots = P.C53.inverse(bg)
    face_center = project(center)
    assert pair(face_center, face_center).h < 0 and pair(face_center, center).h < 0
    assert all(pair(face_center, projected[i]).h < 0 for i in sides)
    records = []
    adjacent = {frozenset((sides[i], sides[(i+1) % 4])) for i in range(4)}
    for a in range(4):
        for b in range(a+1, 4):
            first, second = sides[a], sides[b]
            x = [pair(v, projected[first]) for v in basis]
            y = [pair(v, projected[second]) for v in basis]
            coefficients = [U.rnd(x[1]*y[2]-x[2]*y[1]), U.rnd(x[2]*y[0]-x[0]*y[2]), U.rnd(x[0]*y[1]-x[1]*y[0])]
            point = P.combine(coefficients, basis)
            square = U.rnd(pair(point, point))
            if frozenset((first, second)) in adjacent:
                assert square.h < 0
                if pair(point, center).l > 0:
                    point = [-x for x in point]
                assert pair(point, center).h < 0
                incident = WHITE+[first, second]
                corner = [[pair(normals[i], normals[j]) for j in incident] for i in incident]
                pivots = P.pd(corner)
                products = [U.rnd(pair(point, n)) for n in normals]
                for j, value in enumerate(products):
                    if j in incident:
                        assert value.l <= 0 <= value.h
                    else:
                        assert value.h < 0, (first, second, j, value.enc())
                records.append(dict(pair=[first, second], kind='finite_vertex', point=P.enc(point),
                    norm=square.enc(), future=pair(point, center).enc(), incident=incident,
                    corner_PD_pivots=pivots, products=P.enc(products)))
            else:
                assert square.l > 0
                values = {i: U.rnd(pair(point, projected[i])) for i in sides if i not in (first, second)}
                plus = [i for i, q in values.items() if q.l > 0]
                minus = [i for i, q in values.items() if q.h < 0]
                assert plus and minus, (first, second, {i:q.enc() for i,q in values.items()})
                records.append(dict(pair=[first, second], kind='both_orientations_infeasible',
                    point=P.enc(point), norm=square.enc(), other_side_products={str(i):q.enc() for i,q in values.items()},
                    positive_orientation_blocker=plus[0], negative_orientation_blocker=minus[0]))
    assert sum(r['kind'] == 'finite_vertex' for r in records) == 4
    affine = affine_inventory(z['old'])
    physical = lambda v: P.combine(v, z['basis'])
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json',
        f'scripts/{DOMAIN}.py', f'results/{DOMAIN}.json', f'scripts/{D.TAG}.py',
        f'results/{D.TAG}.json', f'scripts/{D.E.TAG}.py', f'results/{D.E.TAG}.json',
        'scripts/fresh_h9_persistent_cusp_affine_20261016.py', 'results/fresh_h9_persistent_cusp_affine_20261016.json']
    return dict(status='PASS_EXACT_60_WALL_CHAMBER_WITH_COMPACT_R8_QUADRILATERAL_FACE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        old_walls=58, total_walls=60, computational_root_basis=z['basis_ids'],
        computational_gram=[P.enc(row) for row in z['exact_gram']], recipes=z['recipes'],
        exact_recipe='For each saved endpoint domain p,rho,A,B, use r(t)=p+sqrt(rho/2)*(t*A+B/t). First solve r.old49=0 in the pi3-to7 profile; then solve the pi5-to30 profile against the first new wall with product -2cos(pi/9). Each equation has negative leading coefficient and positive constant, hence exactly one positive root.',
        new_normals=[dict(index=i, computational=P.enc(normals[i]), physical=P.enc(physical(normals[i]))) for i in (58, 59)],
        new_contact_checks=checks, old_persistent_vertex_products=P.enc(vertex_products),
        common_center=P.enc(center), physical_common_center=P.enc(physical(center)),
        center_perturbation_power=power, center_square=pair(center, center).enc(),
        center_future=pair(center, old_center).enc(), center_products=P.enc(center_products),
        facet_witnesses=facets, face_support=WHITE, face_cyclic_sides=sides,
        white_PD_pivots=white_pd, white_inverse_pivots=white_inverse_pivots,
        projected_side_normals={str(i):P.enc(v) for i,v in projected.items()},
        projected_side_gram=[P.enc(row) for row in projected_gram],
        face_basis=basis_ids, face_basis_inverse_pivots=basis_inverse_pivots,
        face_center=P.enc(face_center), face_center_square=pair(face_center, face_center).enc(),
        complete_six_pair_intersections=records, affine_inventory=affine,
        proof='All new native contacts are exact finite Coxeter contacts or strictly dotted. A strict future center and one strict future facet witness per wall prove a genuine 60-facet chamber. The seven-root white support is positive definite. In its Lorentz3 orthogonal complement, the four side inequalities have full rank and therefore a pointed cone. All six possible pair-intersection rays are checked: exactly four future timelike rays are feasible; both orientations of each remaining line violate a side inequality. The four feasible rays satisfy every other native wall strictly. Their convex hull is therefore the entire actual compact white H2 face, not merely four locally compatible corners.',
        scope='A coordinated two-wall completion of an actual H2 facet inside the original H4+A2 H3 face. The intrinsic quadrilateral angles need not be Coxeter. No assertion of ambient finite volume or complete ambient cone is made.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    output = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
        print('PASS saved coordinated quadrilateral completion replay')
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
