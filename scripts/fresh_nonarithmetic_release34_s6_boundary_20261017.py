#!/usr/bin/env python3
"""Complete exact intrinsic S6 boundary after deleting34 and replacing52.

Enumerate every triple of positive-square projected sides in a bounded
projective chart in S6-perp.  Causal projections are proved redundant on the
future H3; no ambient65 cone or floating-point acceptance test is used.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path

import fresh_nonarithmetic_s6_release34_exact_20261017 as X

P, U, I, D = X.P, X.U, X.I, X.D
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_release34_s6_boundary_20261017'
S = [0, 1, 2, 3, 5, 6]
FRAME = [8, 30, 44, 52]


def enc_matrix(a):
    return [P.enc(row) for row in a]


def overlap(x, y):
    return x.l <= y.h and y.l <= x.h


def ledger(normals, pair, h, changed):
    normals = dict(enumerate(normals)) if isinstance(normals, list) else normals
    labels = sorted(normals)
    original_labels = D.E.labels()
    prescribed = {
        52: {**{i: 2 for i in S + [8]}, 30: 4},
        63: {**{i: 2 for i in S + [30, 44, 65]}, 64: 27},
        64: {**{i: 2 for i in S + [30, 52, 65]}, 63: 27},
        65: {**{i: 2 for i in S + [8, 63, 64]}, 44: 15, 52: 4},
    } if changed else {}

    def literal(i, j):
        if i == j:
            q = I(2)
        elif j in prescribed.get(i, {}) or i in prescribed.get(j, {}):
            q = -D.theta(prescribed.get(i, {}).get(j, prescribed.get(j, {}).get(i)))
        elif i < 55 and j < 55 and original_labels[i][j]:
            # Replacement52 retains all its old finite contacts; every other
            # relation involving it is newly certified strictly dotted.
            q = -D.theta(original_labels[i][j])
        else:
            return U.rnd(pair(normals[i], normals[j]))
        assert overlap(q, pair(normals[i], normals[j]))
        return q

    white_gram = [[literal(i, j) for j in S] for i in S]
    white_pivots = P.pd(white_gram)
    right_frame = []
    for i in FRAME:
        for j in S:
            if changed and i == 52:
                assert prescribed[i][j] == 2
            else:
                assert original_labels[i][j] == 2
            q = literal(i, j)
            assert q.l == q.h == 0
            right_frame.append([i, j])
    gram = [[literal(i, j) for j in FRAME] for i in FRAME]
    inverse, inverse_pivots = P.C53.inverse(gram)
    h_products = [U.rnd(pair(h, normals[i])) for i in FRAME]
    center = U.imv(inverse, h_products)
    center_square = U.idot(h_products, center)
    assert center_square.h < 0
    native_center_products = {i: U.rnd(pair(h, normals[i])) for i in labels}
    assert all(q.h < 0 for q in native_center_products.values())

    projected, sides, causal = [], {}, []
    for i in labels:
        if i in S:
            continue
        covector = [literal(i, j) for j in FRAME]
        vector = U.imv(inverse, covector)
        square = U.idot(covector, vector)
        future = U.idot(covector, center)
        assert future.h < 0, ('projection future chart', i, future.enc())
        row = dict(wall=i, covector=P.enc(covector), coordinates=P.enc(vector),
                   square=square.enc(), projected_center_product=future.enc())
        if square.l > 0:
            row['kind'] = 'actual_hyperbolic_side'
            sides[i] = covector
        else:
            # This particular geometry has no null or zero projection.  A
            # failed strict test is unresolved, never silently discarded.
            assert square.h < 0, ('unresolved projected norm', i, square.enc())
            row['kind'] = 'future_timelike_redundant_on_H3'
            causal.append(i)
        projected.append(row)
    side_labels = sorted(sides)
    assert all(i in sides for i in FRAME)

    # The sum is a strict dual functional: each side pairing is <=0, and
    # equality of their sum forces all to vanish; FRAME spans S6-perp.
    # Thus its level -1 section of the pointed cone is bounded.  Timelikeness
    # also fixes a common future chart for all finite rays and edge witnesses.
    chart_covector = [U.rnd(sum((v[j] for v in sides.values()), I())) for j in range(4)]
    chart = U.imv(inverse, chart_covector)
    chart_square = U.idot(chart_covector, chart)
    chart_center = U.idot(chart_covector, center)
    assert chart_square.h < 0 and chart_center.h < 0
    normalized_center = [U.rnd(-v / chart_center) for v in center]
    assert all(U.idot(v, normalized_center).h < 0 for v in sides.values())

    triples, rays = [], []
    for support in combinations(side_labels, 3):
        matrix = [sides[i] for i in support] + [chart_covector]
        inv, pivots = P.C53.inverse(matrix)
        vector = U.imv(inv, [I(), I(), I(), I(-1)])
        gates = {i: U.idot(cov, vector) for i, cov in sides.items() if i not in support}
        positive = next((i for i, q in gates.items() if q.l > 0), None)
        record = dict(support=list(support), inverse_pivots=pivots)
        if positive is not None:
            record.update(kind='rejected_by_strict_positive_side', blocker=positive,
                          product=gates[positive].enc())
            triples.append(record)
            continue
        assert all(q.h < 0 for q in gates.values()), ('nonsimple_or_unresolved_vertex', support)
        for i in support:
            q = U.idot(sides[i], vector)
            assert q.l <= 0 <= q.h
        level = U.idot(chart_covector, vector)
        assert level.l <= -1 <= level.h
        square = U.idot(vector, U.imv(gram, vector))
        future = U.idot(h_products, vector)
        if square.h < 0:
            assert future.h < 0
            kind = 'finite_vertex'
        else:
            assert square.l > 0, ('ideal_or_unresolved_ray', support, square.enc())
            kind = 'exterior_ray'
        record.update(kind=kind, ray=len(rays))
        triples.append(record)
        rays.append(dict(index=len(rays), active=list(support), kind=kind,
                         coordinates=P.enc(vector), norm=square.enc(),
                         common_center_product=future.enc(), chart_level=level.enc(),
                         nonincident_side_products=[dict(wall=i, product=q.enc()) for i, q in gates.items()]))
    assert len(triples) == comb(len(sides), 3)
    assert len({tuple(r['active']) for r in rays}) == len(rays)

    pair_index = defaultdict(list)
    for r in rays:
        for support in combinations(r['active'], 2):
            pair_index[support].append(r['index'])
    edges, counts = [], Counter()
    for support, ends in sorted(pair_index.items()):
        assert len(ends) == 2, ('non_edge_or_nonsimple_face', support, ends)
        a, b = [[I(F(lo), F(hi)) for lo, hi in rays[i]['coordinates']] for i in ends]
        aa = U.idot(a, U.imv(gram, a))
        ab = U.idot(a, U.imv(gram, b))
        bb = U.idot(b, U.imv(gram, b))
        determinant = U.rnd(aa * bb - ab * ab)
        finite = sum(rays[i]['kind'] == 'finite_vertex' for i in ends)
        witness = None
        if finite:
            assert determinant.h < 0
            kind = 'compact_geodesic_edge' if finite == 2 else 'one_open_end'
        elif ab.h < 0 and determinant.h < 0:
            ca, cb = U.rnd(bb - ab), U.rnd(aa - ab)
            assert ca.l > 0 and cb.l > 0
            v = [U.rnd(ca * x + cb * y) for x, y in zip(a, b)]
            identity = U.rnd(determinant * (aa + bb - 2 * ab))
            direct = U.idot(v, U.imv(gram, v))
            chart_product = U.idot(v, chart_covector)
            assert identity.h < 0 and chart_product.h < 0 and overlap(identity, direct)
            kind = 'two_open_ends'
            witness = dict(positive_coefficients=P.enc([ca, cb]), coordinates=P.enc(v),
                           norm=identity.enc(), direct_norm=direct.enc(), chart_product=chart_product.enc())
        else:
            assert ab.l > 0 or determinant.l > 0, ('unresolved_edge_cone', support)
            kind = 'outside_hyperbolic_space'
        counts[kind] += 1
        edges.append(dict(sides=list(support), endpoints=ends,
                          endpoint_supports=[rays[i]['active'] for i in ends],
                          kind=kind, square_pairing_bounds=P.enc([aa, ab, bb]),
                          determinant=determinant.enc(), future_timelike_witness=witness))
    assert 3 * len(rays) == 2 * len(edges)
    assert len(rays) - len(edges) + len(sides) == 2
    facets = []
    for i in side_labels:
        vertices = [r['index'] for r in rays if i in r['active']]
        incident = [e for e in edges if i in e['sides']]
        assert len(vertices) == len(incident) >= 3
        compact = all(rays[j]['kind'] == 'finite_vertex' for j in vertices)
        if compact:
            assert all(e['kind'] == 'compact_geodesic_edge' for e in incident)
        facets.append(dict(wall=i, compact=compact, ray_indices=vertices,
                           vertex_supports=[rays[j]['active'] for j in vertices],
                           projective_edges=[e['sides'] for e in incident],
                           edge_counts=dict(Counter(e['kind'] for e in incident))))
    outside = {r['index'] for r in rays if r['kind'] == 'exterior_ray'}
    adjacency = {i: set() for i in outside}
    for e in edges:
        a, b = e['endpoints']
        if a in outside and b in outside:
            adjacency[a].add(b)
            adjacency[b].add(a)
    components, remaining = [], set(outside)
    while remaining:
        start = min(remaining)
        seen, todo = {start}, [start]
        for i in todo:
            for j in sorted(adjacency[i] - seen):
                seen.add(j)
                todo.append(j)
        remaining -= seen
        components.append(dict(ray_indices=sorted(seen),
                               supports=[rays[i]['active'] for i in sorted(seen)],
                               internal_edges=[dict(sides=e['sides'], kind=e['kind']) for e in edges
                                               if set(e['endpoints']) <= seen]))
    return dict(labels=labels, white_support=S, white_positive_pivots=white_pivots,
                residual_basis=FRAME, basis_right_incidence=right_frame,
                residual_gram=enc_matrix(gram), residual_inverse_pivots=inverse_pivots,
                projected_common_center=P.enc(center), projected_center_square=center_square.enc(),
                native_center_products=[dict(wall=i, product=q.enc()) for i, q in native_center_products.items()],
                projections=projected, actual_sides=side_labels, timelike_redundant_walls=causal,
                chart_covector=P.enc(chart_covector), chart_coordinates=P.enc(chart),
                chart_square=chart_square.enc(), chart_center_product=chart_center.enc(),
                normalized_strict_center=P.enc(normalized_center),
                triples=triples, singular_triples=0, nonsimple_rays=0, ideal_rays=0,
                rays=rays, ray_counts=dict(Counter(r['kind'] for r in rays)),
                edges=edges, edge_counts=dict(counts),
                open_ends=counts['one_open_end'] + 2 * counts['two_open_ends'],
                facets=facets, exterior_components=components)


def run():
    primary = json.loads((ROOT / 'results' / (X.TAG + '.json')).read_text())
    assert primary['status'] == 'PASS_EXACT_RELEASE34_65_WALL_CHAMBER_WITH_CLOSED_THREE_CAP_PATCH'
    changed = X.geometry()
    old = X.C.geometry()
    assert changed['basis_ids'] == old['basis_ids']
    assert enc_matrix(changed['exact_gram']) == enc_matrix(old['exact_gram'])
    h = changed['common_center']
    assert primary['wall_labels'] == sorted(changed['normals'])
    assert primary['common_center'] == P.enc(h)
    assert primary['computational_root_basis'] == changed['basis_ids']
    assert primary['computational_gram'] == enc_matrix(changed['exact_gram'])
    for normal in primary['new_normals']:
        assert normal['computational'] == P.enc(changed['normals'][normal['index']])
    before = ledger(old['normals'], old['pair'], h, False)
    print('EXACT_S6_BEFORE', before['ray_counts'], before['edge_counts'], flush=True)
    after = ledger(changed['normals'], changed['pair'], h, True)
    print('EXACT_S6_AFTER', after['ray_counts'], after['edge_counts'], flush=True)
    assert before['ray_counts'] == dict(finite_vertex=19, exterior_ray=9)
    assert after['ray_counts'] == dict(finite_vertex=23, exterior_ray=9)
    assert before['edge_counts'] == dict(compact_geodesic_edge=22, one_open_end=13,
                                       two_open_ends=3, outside_hyperbolic_space=4)
    assert after['edge_counts'] == dict(compact_geodesic_edge=28, one_open_end=13,
                                      two_open_ends=2, outside_hyperbolic_space=5)
    assert before['open_ends'] == 19 and after['open_ends'] == 17
    for i, vertices in [(63, 4), (64, 4), (65, 5)]:
        f = next(r for r in after['facets'] if r['wall'] == i)
        assert f['compact'] and len(f['ray_indices']) == vertices
    actual = lambda d: {tuple(e['sides']): e for e in d['edges'] if e['kind'] != 'outside_hyperbolic_space'}
    b, a = actual(before), actual(after)
    changes = dict(
        removed=[dict(sides=list(k), kind=b[k]['kind']) for k in sorted(b.keys() - a.keys())],
        added=[dict(sides=list(k), kind=a[k]['kind']) for k in sorted(a.keys() - b.keys())],
        changed=[dict(sides=list(k), before=b[k]['kind'], after=a[k]['kind'])
                 for k in sorted(a.keys() & b.keys()) if a[k]['kind'] != b[k]['kind']])
    assert len(changes['removed']) == 5 and len(changes['added']) == 10 and len(changes['changed']) == 6
    paths = [f'scripts/{TAG}.py', f'scripts/{X.TAG}.py', f'results/{X.TAG}.json',
             f'scripts/{X.C.TAG}.py', f'results/{X.C.TAG}.json',
             f'scripts/{P.TAG}.py', f'scripts/{P.C53.TAG}.py', f'scripts/{D.TAG}.py']
    return dict(status='PASS_COMPLETE_EXACT_RELEASE34_S6_GEODESIC_BOUNDARY',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                ambient_dimension=9, intrinsic_dimension=3, fixed_support=S,
                before63=before, after65=after, geodesic_edge_changes=changes,
                projection_scope='Only positive-square projected walls define H3 sides. All other nonwhite projections are strictly future timelike, hence redundant on the entire future H3. They are deliberately omitted from the auxiliary exterior projective cone; its exterior rays are not asserted to equal the full ambient cone slice.',
                completeness_proof='The four right-to-S6 roots8,30,44,52 have an invertible Gram and span the Lorentz4 complement. The projected center is strictly future timelike and strict for every side. The sum of all side covectors is strictly negative on every nonzero feasible vector because the side normals span; its level minus-one affine section is bounded. All triples of side equations with this chart have certified nonsingular4x4 systems. Every triple is either rejected by a strictly positive side or gives a simple vertex with all other sides strict. Therefore the listed28/32 rays exhaust the projective cones, with no singular triple, nonsimple ray or ideal ray. Two common side equations define each projective edge; the complete vertex list gives its two endpoints. Their2x2 Lorentz Gram exactly classifies the positive cone as compact, one-open, two-open or outside H3, with explicit timelike witnesses for every two-open edge.',
                geometric_scope='Exact complete intrinsic boundary comparison for the actual S6 face, not a full ambient65 cone. Label52 denotes different normals before and after the explicit replacement. The three new actual facets are complete compact quadrilaterals63/64 and pentagon65; seventeen geodesic open ends remain. The inherited strict ambient escape survives, so no finite-volume completion is claimed.',
                arithmetic='Shared frozen geometry and outward rational interval engine; new direct4x4 residual enumeration. Numerical boundary files are neither inputs nor acceptance certificates.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], 'open ends', data['before63']['open_ends'], data['after65']['open_ends'], flush=True)
