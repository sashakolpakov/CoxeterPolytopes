#!/usr/bin/env python3
"""Exact finite-edge and S6 boundary ledger in the62-wall root coordinates."""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_pi9_prism_cone_20261017 as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi9_prism_edge_ledger_20261017'
Q, U, I, P = C.Q, C.U, C.I, C.P


def run():
    cone = Q.read(C.TAG)
    assert cone['status'] == 'PASS_COMPLETE_62_WALL_CONE_AND_STRICT_ESCAPE'
    rays = cone['final_rays']
    assert all(len(r['active']) == 9 for r in rays if r['kind'] == 'finite_vertices')
    index = defaultdict(list)
    for i, ray in enumerate(rays):
        for support in combinations(ray['active'], 8):
            index[support].append(i)
    finite_edges, counts = [], Counter()
    for support, ends in sorted(index.items()):
        finite = [i for i in ends if rays[i]['kind'] == 'finite_vertices']
        if not finite:
            continue
        assert len(ends) == 2
        kind = 'compact_geodesic_edge' if len(finite) == 2 else 'open_geodesic_edge'
        counts[kind] += 1
        finite_edges.append(dict(support=list(support), endpoints=ends, kind=kind,
                                 endpoint_supports=[rays[i]['active'] for i in ends]))
    assert 2*counts['compact_geodesic_edge']+counts['open_geodesic_edge'] == 9*cone['counts']['finite_vertices']
    q = [Q.vector(row) for row in cone['computational_gram']]
    center = Q.vector(cone['computational_common_center'])
    hc = U.imv(q, center)
    face_edges, face_counts = [], Counter()
    for support, ends in sorted(index.items()):
        if not set(P.S) <= set(support):
            continue
        assert len(ends) == 2 and all(len(rays[i]['active']) == 9 for i in ends)
        a, b = (Q.vector(rays[i]['normalized_computational_coordinate_bounds']) for i in ends)
        aa, ab, bb = U.idot(a, U.imv(q, a)), U.idot(a, U.imv(q, b)), U.idot(b, U.imv(q, b))
        finite = sum(rays[i]['kind'] == 'finite_vertices' for i in ends)
        determinant = U.rnd(aa*bb-ab*ab)
        if finite:
            assert determinant.h < 0
            kind = 'compact_geodesic_edge' if finite == 2 else 'one_open_end'
            witness = None
        elif ab.h < 0 and determinant.h < 0:
            ca, cb = U.rnd(bb-ab), U.rnd(aa-ab)
            assert ca.l > 0 and cb.l > 0
            v = [U.rnd(ca*x+cb*y) for x, y in zip(a, b)]
            direct, future = U.idot(v, U.imv(q, v)), U.idot(v, hc)
            norm = U.rnd(determinant*(aa+bb-2*ab))
            assert norm.l <= direct.h and direct.l <= norm.h
            assert norm.h < 0 and future.h < 0
            witness = dict(positive_coefficients=P.enc([ca, cb]), norm=norm.enc(), future=future.enc(),
                           norm_identity='(aa*bb-ab^2)*(aa+bb-2ab)', direct_coordinate_norm=direct.enc())
            kind = 'two_open_ends'
        else:
            assert ab.l > 0 or determinant.l > 0
            kind, witness = 'outside_hyperbolic_space', None
        face_counts[kind] += 1
        face_edges.append(dict(support=list(support), sides=sorted(set(support)-set(P.S)),
                               endpoints=ends, endpoint_supports=[rays[i]['active'] for i in ends], kind=kind,
                               square_pairing_bounds=P.enc([aa, ab, bb]), determinant=determinant.enc(),
                               future_timelike_witness=witness))
    face_rays = [r for r in rays if set(P.S) <= set(r['active'])]
    assert len(face_rays)*3 == 2*len(face_edges)
    exterior = {i for i, r in enumerate(rays) if set(P.S) <= set(r['active']) and r['kind'] == 'hyperideal_rays'}
    outside_graph = {i: [] for i in exterior}
    for edge in face_edges:
        if set(edge['endpoints']) <= exterior:
            a, b = edge['endpoints']
            outside_graph[a].append(dict(ray=b, kind=edge['kind'], sides=edge['sides']))
            outside_graph[b].append(dict(ray=a, kind=edge['kind'], sides=edge['sides']))
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json',
             f'results/{C.CAP}.json', f'scripts/{Q.TAG}.py']
    return dict(status='PASS_COMPLETE_62_WALL_FINITE_EDGE_AND_S6_BOUNDARY_LEDGER',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                indexed_eight_node_supports=len(index), finite_incident_edge_counts=dict(counts),
                finite_incident_edges=finite_edges, face_support=P.S,
                face_ray_count=len(face_rays), face_projective_edges=len(face_edges),
                face_edge_counts=dict(face_counts), face_edges=face_edges,
                exterior_adjacencies=[dict(ray=i, support=rays[i]['active'], neighbors=outside_graph[i]) for i in sorted(exterior)],
                proof='The complete ambient cone supplies every ray and active set. Finite rays and allS6rays have nine independent active normals. Every indexed eight-node support incident to one of these rays has exactly two endpoints, giving all relevant projective edges. Timelike endpoints give actual geodesic edges. For two spacelike endpoints the exact negative mutual product and negative two-by-two determinant criterion detects whether their positive cone meets hyperbolic space, with a future timelike witness. Computations use the exact root-basis Gram exported by the62cone.',
                scope='Exact extraction from the complete62-wall primary cone, not independent reconstruction of its roots or rays. The outside graph includes both geodesic edges with two open ends and edges wholly outside hyperbolic space; their kinds remain explicit. No finite-volume claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['finite_incident_edge_counts'], result['face_edge_counts'], flush=True)
