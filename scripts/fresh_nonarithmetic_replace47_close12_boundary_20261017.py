#!/usr/bin/env python3
"""Complete exact intrinsic S6 boundary after replacing47 and closing12.

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

import fresh_nonarithmetic_replace47_close12_20261017 as X

P, U, I, D = X.P, X.U, X.I, X.D
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_replace47_close12_boundary_20261017'
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
        q=X.literal(dict(normals=normals,pair=pair),i,j)
        assert overlap(q,pair(normals[i],normals[j]))
        return q

    white_gram = [[literal(i, j) for j in S] for i in S]
    white_pivots = P.pd(white_gram)
    white_inverse, _ = P.C53.inverse(white_gram)
    white_products={i:[literal(i,j)for j in S]for i in labels if i not in S}
    white_coefficients={i:U.imv(white_inverse,v)for i,v in white_products.items()}
    def projected_pair(i,j):
        return U.rnd(literal(i,j)-U.idot(white_products[i],white_coefficients[j]))
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
        direct_square = U.idot(vector, U.imv(gram, vector))
        # Same exact ray, evaluated through the Gram of its three active
        # projected normals and the chart. The inverse bottom-right entry
        # is its square because the chart value is -1. This avoids loss of
        # precision for very short finite rays in the residual coordinates.
        active_gram=[[projected_pair(i,j)for j in support]for i in support]
        active_inverse,active_pivots=P.C53.inverse(active_gram)
        chart_products=[U.idot(sides[i],chart)for i in support]
        complement_square=U.rnd(chart_square-U.idot(chart_products,U.imv(active_inverse,chart_products)))
        assert complement_square.l>0 or complement_square.h<0
        gram_square=U.rnd(I(1)/complement_square)
        assert overlap(direct_square,gram_square)
        square=I(max(direct_square.l,gram_square.l),min(direct_square.h,gram_square.h))
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
                         direct_norm=direct_square.enc(),active_gram=enc_matrix(active_gram),
                         active_inverse_pivots=active_pivots,chart_complement_square=complement_square.enc(),
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
        aa = I(*[F(x)for x in rays[ends[0]]['norm']])
        ab = U.idot(a, U.imv(gram, b))
        bb = I(*[F(x)for x in rays[ends[1]]['norm']])
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
    primary_path=ROOT/'results'/(X.TAG+'.json')
    primary=json.loads(primary_path.read_text())
    assert primary['status']=='PASS_EXACT_64_WALL_REPLACEMENT47_CLOSE12'
    before_path=ROOT/'results/fresh_nonarithmetic_release54_replace12_boundary_20261017.json'
    before=json.loads(before_path.read_text())['after64']
    z=X.geometry()
    assert primary['new_normal']==P.enc(z['normals'][47])
    assert primary['common_center']==P.enc(z['common_center'])
    after=ledger(z['normals'],z['pair'],z['common_center'],True)
    print('AFTER',after['ray_counts'],after['edge_counts'],flush=True)
    assert after['ray_counts']==dict(finite_vertex=25,exterior_ray=7)
    assert after['edge_counts']==dict(compact_geodesic_edge=32,one_open_end=11,two_open_ends=2,outside_hyperbolic_space=3)
    assert before['open_ends']==16 and after['open_ends']==15
    assert len(after['triples'])==816 and len(after['actual_sides'])==18
    assert not after['singular_triples'] and not after['nonsimple_rays'] and not after['ideal_rays']
    bc={f['wall']for f in before['facets']if f['compact']}
    ac={f['wall']for f in after['facets']if f['compact']}
    assert bc<=ac and ac-bc=={12}
    for wall,num in [(8,7),(12,4),(63,4),(64,4),(65,5)]:
        f=next(x for x in after['facets']if x['wall']==wall)
        assert f['compact'] and len(f['ray_indices'])==num
    f12=next(x for x in after['facets']if x['wall']==12)
    assert sorted(f12['vertex_supports'])==[[7,8,12],[7,12,47],[8,12,52],[12,47,52]]
    a={tuple(r['sides']):r for r in after['edges']}
    b={tuple(r['sides']):r for r in before['edges']}
    changes=dict(removed=[dict(sides=list(k),kind=b[k]['kind'])for k in sorted(b.keys()-a.keys())],
        added=[dict(sides=list(k),kind=a[k]['kind'])for k in sorted(a.keys()-b.keys())],
        changed=[dict(sides=list(k),before=b[k]['kind'],after=a[k]['kind'])for k in sorted(a.keys()&b.keys())if a[k]['kind']!=b[k]['kind']])
    paths=[Path(__file__),Path(X.__file__),primary_path,before_path,
           ROOT/'scripts/fresh_nonarithmetic_release34_s6_boundary_20261017.py']
    return dict(status='PASS_COMPLETE_EXACT_S6_REPLACEMENT47_BOUNDARY_16_TO_15',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        before64=before,after64=after,changed_edges=changes,
        compact_side_faces_lost=sorted(bc-ac),compact_side_faces_gained=sorted(ac-bc),
        compact12_quadrilateral=f12,
        completeness='The frozen direct4 residual ledger algorithm is copied with all old/new finite products routed through the new literal Gram accessor. All816 triple systems are nonsingular; every candidate is rejected strictly or is a simple ray.32rays/48edges exhaust the bounded projective section. Causal projections are future timelike and redundant on futureH3. Exact two-ray cone Grams distinguish all four geodesic edge types.',
        scope='64actualwalls remain after replacing47. The12face becomes a compact quadrilateral and8heptagon survives. Old compact edge7/30 becomes open, so this is not monotone containment. Complete198 ambient finitevertices and strict escape are separately certified; no full ambient exterior cone is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
