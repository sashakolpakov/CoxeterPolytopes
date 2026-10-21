#!/usr/bin/env python3
"""Bounded exact geometry for a next nonarithmetic H9 cap experiment.

This is not a cap search.  Arithmetic is shared with the frozen independent
octic bridge review; no construction-producer arithmetic is imported.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import fresh_h9_joint_corner_20261015_bridge_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_escape_target_20261016'
BASE = 'results/fresh_h9_joint_pair_20261015_two_extension.json'
BRIDGE = 'results/fresh_h9_joint_corner_20261015_bridge.json'
REVIEW = 'results/fresh_h9_joint_corner_20261015_bridge_review.json'
C, K, Z, O = R.C, R.K, R.Z, R.O


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def run():
    base, bridge, review = [json.loads((ROOT / p).read_text())
                            for p in (BASE, BRIDGE, REVIEW)]
    for document in (base, bridge, review):
        for p, digest in document['sources'].items():
            assert sha(p) == digest, p
    assert review['sources'][BRIDGE] == sha(BRIDGE)
    Q, normals, H = R.m4(base['form']), R.m4(base['normals']), R.v4(base['center'])
    normals.append(R.v8(bridge['cap']['normal']))
    E = R.v4(base['strict_escape']['point'])
    pair = lambda x, y: C.pair(x, y, Q)
    covectors = [C.mv(Q, n) for n in normals]
    gram = [[C.dot(n, c) for c in covectors] for n in normals]
    assert len(normals) == 33
    inherited_cycle = 4*gram[11][12]**2/(gram[11][11]*gram[12][12])
    assert inherited_cycle == 4*R.BASIS[1]
    assert R.sg(pair(E, E)) == 1 and R.sg(pair(E, H)) == -1
    assert all(R.sg(pair(E, n)) == -1 for n in normals)
    assert all(R.sg(pair(H, n)) == -1 for n in normals)
    incidence = bridge['finite_vertex_incidence']

    def vertex(index):
        active = incidence[index]
        assert R.pd([[gram[i][j] for j in active] for i in active])[0]
        v = C.null_one([covectors[i] for i in active])
        v = [a / -pair(v, H) for a in v]
        assert R.sg(pair(v, v)) == -1 and pair(v, H) == -O
        values = [pair(v, n) for n in normals]
        assert [R.sg(x) for x in values] == [0 if i in active else -1 for i in range(33)]
        return v

    def projection(v, support):
        G = [[gram[i][j] for j in support] for i in support]
        positive, pivots = R.pd(G)
        assert positive
        inverse = C.inverse(G)
        coefficients = C.mv(inverse, [pair(v, normals[i]) for i in support])
        p = [v[k] - sum((a * normals[i][k] for a, i in zip(coefficients, support)), Z)
             for k in range(10)]
        assert all(not pair(p, normals[i]) for i in support)
        return p, G, pivots, inverse

    # This natural old support is blocked, even for new normals in larger fields.
    S = [0, 1, 2, 5, 6, 10, 29]
    P, GS, pivotsS, _ = projection(E, S)
    assert R.sg(pair(P, P)) == -1 and R.sg(pair(P, H)) == -1
    blocked_vertices = [20, 21, 90]
    vectors = [vertex(i) for i in blocked_vertices]
    assert all(set(S) <= set(incidence[i]) for i in blocked_vertices)
    coordinates = C.eliminate(vectors)[1]
    assert len(coordinates) == 3
    coefficients = C.mv(C.inverse([[v[k] for v in vectors] for k in coordinates]),
                        [P[k] for k in coordinates])
    assert all(R.sg(a) == 1 for a in coefficients)
    assert P == [sum((a * v[k] for a, v in zip(coefficients, vectors)), Z)
                 for k in range(10)]
    print('S: exact positive finite-vertex cone certificate PASS', flush=True)

    # Replace support wall29 by wall30.  This support has a spacelike escape
    # projection and two consecutive actual corners in its two-dimensional face.
    T = [0, 1, 2, 5, 6, 10, 30]
    EP, GT, pivotsT, JT = projection(E, T)
    assert R.sg(pair(EP, EP)) == 1
    phi = (O + R.BASIS[1]) / 2
    expected = [[2*O if i == j else Z for j in range(7)] for i in range(7)]
    for i, j, g in ((0, 1, -phi), (1, 2, -O), (3, 4, -O), (3, 5, -O)):
        expected[i][j] = expected[j][i] = g
    assert GT == expected  # H3 + A3 + A1, in the literal displayed ordering.
    Tvertices = [i for i, ids in enumerate(incidence) if set(T) <= set(ids)]
    assert Tvertices == [21, 91]
    actual = [vertex(i) for i in Tvertices]
    assert [sorted(set(incidence[i]) - set(T)) for i in Tvertices] == [[29, 31], [29, 32]]
    projected_normals = {}
    for i, q in ((29, 3*O/2), (31, O/2), (32, 3*O/2)):
        p, _, _, _ = projection(normals[i], T)
        assert pair(p, p) == q
        projected_normals[str(i)] = dict(point=R.enc(p), norm=R.enc(q))
    # For a norm2 candidate perpendicular to T, a pi/3 contact to wall31
    # has zero Schur determinant; every larger finite label is worse.
    assert (2*O) * (O/2) - O == Z

    # Unique polar line of the two open-side walls: both orientations fail.
    polar_support = T + [31, 32]
    polar = C.null_one([covectors[i] for i in polar_support])
    polar_scale = -pair(polar, H)
    polar = [a / polar_scale for a in polar]
    assert pair(polar, H) == -O
    polar_norm = pair(polar, polar)
    polar_pairs = [pair(polar, n) for n in normals]
    assert R.sg(polar_norm) == 1
    assert all(not polar_pairs[i] for i in polar_support)
    assert R.sg(polar_pairs[29]) == 1
    assert R.sg(polar_pairs[4]) == -1
    assert R.sg(pair(polar, E)) == -1
    print('T: signature, actual corners, escape and both polar orientations PASS', flush=True)

    paths = [f'scripts/{TAG}.py', BASE, BRIDGE, REVIEW,
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    audit_paths = []
    for prefix in ('fresh_sqrt2sqrt5_lanner_catalogue_20261016_review',
                   'fresh_compact21_single_field_20261016',
                   'fresh_h9_witt_twist_20261016_review'):
        audit_paths.extend((f'scripts/{prefix}.py', f'results/{prefix}.json', f'results/{prefix}.md'))
    audit_docs = [json.loads((ROOT / p).read_text()) for p in audit_paths if p.endswith('.json')]
    assert audit_docs[0]['compact_dimension_upper_bound_candidate'] == 20
    assert audit_docs[1]['remaining_fields'] == ['Q(sqrt5)']
    assert audit_docs[1]['dimension21_excluded'] is False
    assert audit_docs[2]['local_Witt_index'] == 3
    paths += audit_paths
    return dict(
        status='PASS_BOUNDED_NONARITHMETIC_ESCAPE_TARGET_GEOMETRY',
        sources={p: sha(p) for p in paths},
        field_encoding='ANP coefficient lists in the frozen octic review primitive element',
        shared_arithmetic='Frozen independent octic bridge review; no producer imports',
        old_chamber=dict(walls=33, finite_vertices=92, finite_edges=253, open_germs=322),
        inherited_nonquasiarithmetic_obstruction=dict(walls=[11, 12],
            Cartan_two_cycle=R.enc(inherited_cycle), expression='4*sqrt(5)',
            negative_conjugate='-4*sqrt(5)',
            scope='Any finite-covolume completion retaining these walls fails the necessary quasi-arithmetic conjugate-PSD condition.'),
        strict_escape=dict(point=R.enc(E), norm=R.enc(pair(E, E))),
        blocked_support=dict(support=S, gram=R.enc(GS), pivots=R.enc(pivotsS),
                             escape_projection=R.enc(P), projection_norm=R.enc(pair(P, P)),
                             finite_vertex_indices=blocked_vertices,
                             finite_vertex_points=R.enc(vectors), positive_coefficients=R.enc(coefficients),
                             conclusion='Every coherent acute additional wall perpendicular to S preserves E.'),
        target_support=dict(support=T, type='H3+A3+A1', gram=R.enc(GT), inverse=R.enc(JT),
                            pivots=R.enc(pivotsT), escape_projection=R.enc(EP),
                            projection_norm=R.enc(pair(EP, EP)),
                            finite_vertex_indices=Tvertices, finite_vertex_points=R.enc(actual),
                            projected_normals=projected_normals,
                            finite_contact_to_31_for_perpendicular_norm2_cap='right only'),
        failed_polar=dict(support=polar_support, normal=R.enc(polar), norm=R.enc(polar_norm),
                          pairings=R.enc(polar_pairs), escape_pairing=R.enc(pair(polar, E)),
                          positive_obstruction_wall=29, opposite_orientation_positive_wall=4),
        proposed_bounded_experiment=dict(
            performed=False, support=T, retained_walls=list(range(33)), cap_norms=[2, 2],
            first_cap_right_wall=31, first_cap_extra_anchor=[3, 7],
            first_anchor_labels=[2, 3, 4, 5], second_cap_endpoint_wall=32,
            second_endpoint_labels=[2, 3, 4, 5], mutual_labels=[2, 3, 4, 5],
            algebraic_recipe_count=128, maximum_isolated_branch_pairs=512,
            required_finite_corners=['T+{31,u}', 'T+{u,v}', 'T+{v,32}'],
            success_conditions=['all 67 new old/mutual Coxeter gates', 'common strict future center',
                                'at least one cap pairs positively with E',
                                'all three displayed full rank9 corners PD and feasible',
                                'both added walls actual facets'],
            scope='A finite proposal family, not an all-label or all-parameter completeness claim. Degenerate conics must be reported separately rather than sampled.'),
        read_only_scope_audit=dict(inputs=audit_paths,
            result='No scope inconsistency found in final F25 review, exact-field composition, or local-Witt setup review. Source/proof read only; their arithmetic replays are inherited.'),
        scope='Read-only reconstruction of two exact supports and one polar line. No new wall or finite-volume chamber is constructed; no unperformed cap search is claimed.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
