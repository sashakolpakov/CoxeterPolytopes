#!/usr/bin/env python3
"""Complete finite-edge incidence and the target three-face boundary at60 walls."""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_pi8_quad_cone_20261017 as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_quad_edge_ledger_20261017'
U, I, P = C.U, C.I, C.P


def run():
    cone = C.read(C.TAG)
    assert cone['status'] == 'PASS_COMPLETE_60_WALL_CONE_AND_STRICT_ESCAPE'
    rays = cone['final_rays']
    assert all(len(r['active']) == 9 for r in rays if r['kind'] == 'finite_vertices')
    index = defaultdict(list)
    for i, r in enumerate(rays):
        for support in combinations(r['active'], 8):
            index[support].append(i)
    # Only supports incident to a simple ray are used below. Exterior rays
    # elsewhere in the cone may have dependent eight-node subsets.
    finite_edges = []
    counts = Counter()
    for support, ends in sorted(index.items()):
        finite = [i for i in ends if rays[i]['kind'] == 'finite_vertices']
        if not finite:
            continue
        assert len(ends) == 2
        kind = 'compact_geodesic_edge' if len(finite) == 2 else 'open_geodesic_edge'
        counts[kind] += 1
        finite_edges.append(dict(support=list(support), endpoints=ends, kind=kind,
            endpoint_supports=[rays[i]['active'] for i in ends],
            finite_endpoint_supports=[rays[i]['active'] for i in finite]))
    assert 2*counts['compact_geodesic_edge']+counts['open_geodesic_edge'] == 9*cone['counts']['finite_vertices']

    _, q, _, _, _, _, _, _ = P.physical_data()
    cap = C.read(C.CAP)
    center = C.vector(cap['physical_common_center'])
    hc = U.imv(q, center)
    face_edges = []
    face_counts = Counter()
    for support, ends in sorted(index.items()):
        if not set(P.S) <= set(support):
            continue
        assert len(ends) == 2 and all(len(rays[i]['active']) == 9 for i in ends)
        a, b = (C.vector(rays[i]['normalized_coordinate_bounds']) for i in ends)
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
            direct_norm, future = U.idot(v, U.imv(q, v)), U.idot(v, hc)
            # The exact identity avoids interval dependency cancellation.
            norm = U.rnd(determinant*(aa+bb-2*ab))
            assert norm.l <= direct_norm.h and direct_norm.l <= norm.h
            assert norm.h < 0 and future.h < 0
            witness = dict(positive_coefficients=P.enc([ca, cb]), norm=norm.enc(), future=future.enc(),
                norm_identity='(aa*bb-ab^2)*(aa+bb-2ab)', direct_coordinate_norm=direct_norm.enc())
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
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json',
             f'results/{C.CAP}.json', f'scripts/{P.TAG}.py']
    return dict(status='PASS_COMPLETE_60_WALL_FINITE_EDGE_AND_S6_BOUNDARY_LEDGER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        indexed_eight_node_supports=len(index), finite_incident_edge_counts=dict(counts),
        finite_incident_edges=finite_edges, face_support=P.S,
        face_ray_count=len(face_rays), face_projective_edges=len(face_edges),
        face_edge_counts=dict(face_counts), face_edges=face_edges,
        proof='Every finite ray and every S6 ray in the inherited complete10-dimensional pointed cone has exactly9 active facets, hence9 independent normals and a simple projective vertex. Other exterior rays may be nonsimple. Removing one active equation at a simple vertex gives a rank8 support and a two-dimensional cone face with exactly two extreme rays. The complete support index checks this literally for every relevant support. At a finite endpoint its orthogonal complement is positive definite, so each such support is spherical and gives an actual geodesic edge. For the S6 face, all39 projective edges are included. Two spacelike endpoints meet hyperbolic space exactly when their mutual product is negative and the2x2 determinant is negative; positive coefficients then give the displayed future timelike witness. A positive product or positive determinant excludes any timelike positive combination.',
        scope='Exact incidence extraction from the complete60-wall primary cone certificate, not independent reconstruction of that cone. It locates the remaining open ends; no finite-volume or compactness claim for the whole chamber or S6 face.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['finite_incident_edge_counts'],data['face_edge_counts'],flush=True)
