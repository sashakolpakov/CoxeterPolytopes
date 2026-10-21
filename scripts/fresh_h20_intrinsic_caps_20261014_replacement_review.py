#!/usr/bin/env python3
"""Independent rational-matrix verification of the radius-four replacement."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from flint import fmpq, fmpq_mat
from fresh_h20_intrinsic_caps_20261014_matrix_review import split, pairing, sign

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_intrinsic_caps_20261014'
OUT = ROOT / f'results/{TAG}_replacement_review.json'


def run():
    paths = [f'results/{TAG}_{suffix}.json' for suffix in
             ('golden_seed', 'golden_extension', 'quarter_face_exact', 'radius4_replacement', 'matrix_review')]
    seed, half, quarter, replacement, reviewed = [json.loads((ROOT/p).read_text()) for p in paths]
    assert reviewed['status'] == 'PASS_INDEPENDENT_FLINT_587_WALL_MATRIX_REPLAY'
    for rec in (replacement, reviewed):
        for p, h in rec['sources'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h
    raw = seed['walls']+half['new_walls']+quarter['new_walls']
    assert replacement['replacement_index'] == 23
    raw[23] = replacement['new_normal']
    form = fmpq_mat(21, 21)
    for i, row in enumerate(seed['spatial_Gram']):
        for j, v in enumerate(row):
            form[i, j] = fmpq(v)
    form[19, 20] = form[20, 19] = -1
    roots = split(raw)
    new = split([raw[23]])
    a, b = pairing(new, roots, form)
    counts = Counter()
    lengths = [fmpq(1)]*22+[fmpq(2)]+[fmpq(1)]*504+[fmpq(3, 2)]*60
    lengths[23] = 4
    for j in range(587):
        if j == 23:
            assert a[0, j] == 16 and b[0, j] == 0
            continue
        den = 4*lengths[j]
        if a[0, j] == b[0, j] == 0:
            label = 'right'
        elif a[0, j] == -den/2 and b[0, j] == 0:
            label = 'label3'
        else:
            s = sign(a[0, j]+den, b[0, j])
            assert s <= 0
            label = 'parallel' if s == 0 else 'dotted'
        counts[label] += 1
    assert dict(counts) == replacement['new_pair_counts']
    center = split([seed['common_strict_center']['point']])
    c0, c1 = pairing(center, new, form)
    assert sign(c0[0, 0], c1[0, 0]) < 0
    oldpoints = [r['point'] for r in seed['finite_endpoints'][1:]]
    oldpoints += [seed['central_facet_point']['point']]
    oldpoints += [r['point'] for r in half['new_actual_facet_points']+quarter['new_actual_facet_points']]
    p0, p1 = pairing(split(oldpoints), new, form)
    assert all(sign(p0[i, 0], p1[i, 0]) < 0 for i in range(len(oldpoints)))
    ep = replacement['new_actual_finite_endpoint']
    point = split([ep['point']])
    n0, n1 = pairing(point, point, form)
    assert n0[0, 0] == fmpq(-8, 3) and n1[0, 0] == 0
    e0, e1 = pairing(point, roots, form)
    zeros = [j for j in range(587) if e0[0, j] == e1[0, j] == 0]
    assert zeros == ep['incident_indices']
    assert all(sign(e0[0, j], e1[0, j]) <= 0 for j in range(587))
    v0, v1 = pairing(point, center, form)
    assert sign(v0[0, 0], v1[0, 0]) < 0
    inc = split([raw[j] for j in zeros])
    g0, g1 = pairing(inc, inc, form)
    assert all(g1[i, j] == 0 for i in range(20) for j in range(20))
    mat = [[g0[i, j] for j in range(20)] for i in range(20)]
    for k in range(20):
        pivot = mat[k][k]
        assert pivot > 0
        for i in range(k+1, 20):
            for j in range(k+1, 20):
                mat[i][j] -= mat[i][k]*mat[k][j]/pivot
    for key, desired in [('annular_escape_cut', 1), ('strict_remaining_escape', -1)]:
        v = split([replacement[key]['point']])
        t0, t1 = pairing(v, new, form)
        assert sign(t0[0, 0], t1[0, 0]) == desired
    assert replacement['strict_remaining_escape']['point'] == quarter['strict_remaining_escape']['point']
    paths += [f'scripts/{TAG}_matrix_review.py', str(Path(__file__).resolve().relative_to(ROOT))]
    return {'status': 'PASS_INDEPENDENT_RADIUS4_REPLACEMENT_MATRIX_REVIEW',
            'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            'changed_pairs': 586, 'pair_counts': dict(counts), 'new_endpoint_norm': '-8/3',
            'retained_old_facet_witness_points': len(oldpoints), 'retained_old_finite_endpoints': 251,
            'annular_escape_cut': True, 'earlier_eighth_escape_survives': True,
            'scope': 'Rechecks all geometry changed by the one-wall replacement using FLINT rational coefficient matrices; inherits unchanged geometry from the independently replayed 587-wall certificate. No finite-volume claim.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == result
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
