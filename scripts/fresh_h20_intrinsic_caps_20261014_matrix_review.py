#!/usr/bin/env python3
"""Independent full matrix replay of the 587-wall construction using FLINT.

Uses coefficient matrices over Q, with phi^2=phi+1, rather than the
constructor's scalar field arithmetic. Checks geometry, not search completeness.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from flint import fmpq, fmpq_mat

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'fresh_h20_intrinsic_caps_20261014'
OUT = ROOT / f'results/{PREFIX}_matrix_review.json'


def sign(a, b):
    a = 2*a+b
    if not b:
        return (a > 0)-(a < 0)
    if not a:
        return (b > 0)-(b < 0)
    if (a > 0) == (b > 0):
        return (a > 0)-(a < 0)
    d = a*a-5*b*b
    assert d
    return ((a > 0)-(a < 0))*((d > 0)-(d < 0))


def split(rows):
    return tuple(fmpq_mat([[fmpq(v[k]) for v in row] for row in rows]) for k in (0, 1))


def pairing(x, y, form):
    a, b = x
    c, d = y
    bd = b*form*d.transpose()
    return a*form*c.transpose()+bd, a*form*d.transpose()+b*form*c.transpose()+bd


def run():
    paths = [f'results/{PREFIX}_{suffix}.json' for suffix in ('golden_seed', 'golden_extension', 'quarter_face_exact')]
    seed, ext, face = [json.loads((ROOT/p).read_text()) for p in paths]
    for record in (seed, ext, face):
        for p, h in record['sources'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h
    raw = seed['walls']+ext['new_walls']+face['new_walls']
    assert len(raw) == 587
    form = fmpq_mat(21, 21)
    for i, row in enumerate(seed['spatial_Gram']):
        for j, value in enumerate(row):
            form[i, j] = fmpq(value)
    form[19, 20] = form[20, 19] = -1
    # Independent rational LDL of the positive spatial block.
    b = [[form[i, j] for j in range(19)] for i in range(19)]
    pivots = []
    for k in range(19):
        pivot = b[k][k]
        assert pivot > 0
        pivots.append(str(pivot))
        for i in range(k+1, 19):
            for j in range(k+1, 19):
                b[i][j] -= b[i][k]*b[k][j]/pivot
    roots = split(raw)
    g0, g1 = pairing(roots, roots, form)
    lengths = [fmpq(1)]*22+[fmpq(2)]+[fmpq(1)]*504+[fmpq(3, 2)]*60
    counts = Counter()
    for i in range(587):
        assert g1[i, i] == 0 and g0[i, i] == lengths[i]**2
        for j in range(i):
            a, b = g0[i, j], g1[i, j]
            den = lengths[i]*lengths[j]
            if a == b == 0:
                label = 'right'
            elif a == -den/2 and b == 0:
                label = 'label3'
            elif a == 0 and b == -den/2:
                label = 'label5'
            else:
                s = sign(a+den, b)
                assert s <= 0, (i, j, str(a), str(b))
                label = 'parallel' if s == 0 else 'dotted'
            counts[label] += 1
    assert dict(counts) == face['all587_pair_counts']

    center_raw = seed['common_strict_center']['point']
    center = split([center_raw])
    c0, c1 = pairing(center, center, form)
    assert c0[0, 0] == -9 and c1[0, 0] == 0
    cr0, cr1 = pairing(center, roots, form)
    assert all(sign(cr0[0, j], cr1[0, j]) < 0 for j in range(587))
    points = [r['point'] for r in seed['finite_endpoints']]
    expected_zeros = [r['incident_wall_indices'] for r in seed['finite_endpoints']]
    norms = [fmpq(-1)]*252
    points.append(seed['central_facet_point']['point'])
    expected_zeros.append([22])
    norms.append(fmpq(-4))
    for block in (ext, face):
        for r in block['new_actual_facet_points']:
            points.append(r['point'])
            expected_zeros.append([r['only_incident_wall']])
            norms.append(fmpq(r['norm'][0]))
            assert fmpq(r['norm'][1]) == 0
    pp = split(points)
    p0, p1 = pairing(pp, roots, form)
    n0, n1 = pairing(pp, pp, form)
    cc0, cc1 = pairing(pp, center, form)
    used = set()
    for i, zeros in enumerate(expected_zeros):
        assert n0[i, i] == norms[i] and n1[i, i] == 0
        assert sign(cc0[i, 0], cc1[i, 0]) < 0
        got = [j for j in range(587) if p0[i, j] == p1[i, j] == 0]
        assert got == sorted(zeros), (i, got, zeros)
        assert all(sign(p0[i, j], p1[i, j]) <= 0 for j in range(587))
        if i < 252:
            # At these all-right cap endpoints, the 20x20 incident block
            # consists of a positive 19x19 seed block and one norm-one root.
            mat = [[g0[u, v] for v in got] for u in got]
            assert all(g1[u, v] == 0 for u in got for v in got)
            for k in range(20):
                pivot = mat[k][k]
                assert pivot > 0
                for a in range(k+1, 20):
                    for b in range(k+1, 20):
                        mat[a][b] -= mat[a][k]*mat[k][b]/pivot
        used.update(got)
    assert used == set(range(587))
    esc = split([face['strict_remaining_escape']['point']])
    e0, e1 = pairing(esc, esc, form)
    assert e0[0, 0] == fmpq(1, 2) and e1[0, 0] == 0
    ec0, ec1 = pairing(esc, center, form)
    assert sign(ec0[0, 0], ec1[0, 0]) < 0
    er0, er1 = pairing(esc, roots, form)
    assert all(sign(er0[0, j], er1[0, j]) < 0 for j in range(587))
    paths.append(str(Path(__file__).resolve().relative_to(ROOT)))
    return {'status': 'PASS_INDEPENDENT_FLINT_587_WALL_MATRIX_REPLAY',
            'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            'signature': [20, 1], 'spatial_LDL_pivots': pivots,
            'wall_count': 587, 'pair_count': sum(counts.values()), 'pair_counts': dict(counts),
            'finite_endpoints': 252, 'facet_witness_points': len(points),
            'all_facets_actual': True, 'escape_norm': '1/2', 'infinite_volume': True,
            'scope': 'Independent full exact Gram, strict center, endpoint rank, facet incidence and strict escape replay. Does not independently reconstruct the selected recipes or enumerate all corners.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == result
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['pair_count'])
