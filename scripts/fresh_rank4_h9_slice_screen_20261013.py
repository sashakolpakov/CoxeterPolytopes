#!/usr/bin/env python3
"""Exact Farkas screening of all two-fundamental-weight H9 slices.

Numerical linear programs propose small supports only. The saved identities
are recomputed over Q(phi), and --verify uses no linear programming.
The surviving geometric direction is not asserted to be reflective.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json

from scripts.fresh_multibasepoint_h9_merge import (
    ROOT, G, read, matrix, inverse, mv, dot, enc,
)

SOURCE = 'results/fresh_rank3_h9_cones_20261013.json'
OUTPUT = ROOT / 'results/fresh_rank4_h9_slice_screen_20261013.json'


def setup():
    source = read(SOURCE)
    gram = matrix(source['gram'])
    walls = matrix(source['cases'][0]['all_walls'])
    v = tuple(map(G.from_json, source['cases'][0]['hyperideal_witness']['vector']))
    einv = inverse(tuple(tuple(row[1:9]) for row in gram[1:9]))
    omega = tuple((G(0), *(einv[k][j] for k in range(8)), G(0)) for j in range(8))
    time = (G(1), *(G(0) for _ in range(9)))
    scalar = (*(G(0) for _ in range(9)), G(1))
    normals = tuple(mv(gram, r) for r in walls)
    target_normal = mv(gram, v)
    return source, gram, walls, v, omega, time, scalar, normals, target_normal


def solve_support(columns, target):
    count = len(columns)
    for rows in combinations(range(len(target)), count):
        square = tuple(tuple(column[i] for column in columns) for i in rows)
        try:
            coeff = mv(inverse(square), tuple(target[i] for i in rows))
        except (ValueError, ZeroDivisionError, AssertionError, StopIteration):
            continue
        assert all(sum((a * col[k] for a, col in zip(coeff, columns)), G(0)) == target[k]
                   for k in range(len(target)))
        assert all(a.sign() >= 0 for a in coeff)
        return coeff
    raise AssertionError('Proposed support has no independent square subsystem')


def check(payload):
    source, gram, walls, v, omega, time, scalar, normals, target_normal = setup()
    assert payload['source'] == SOURCE and payload['source_case'] == 0
    assert payload['wall_count'] == len(walls) == 382
    assert payload['target_escape'] == enc(v)
    assert dot(v, target_normal).sign() > 0
    seen = set()
    for row in payload['slices']:
        i, j = row['fundamental_weight_indices']
        assert 0 <= i < j < 8 and (i, j) not in seen
        seen.add((i, j))
        basis = (time, omega[i], omega[j], scalar)
        target = tuple(dot(s, target_normal) for s in basis)
        if row['status'] == 'EXACT_SLICE_OBSTRUCTION':
            coeff = tuple(map(G.from_json, row['nonnegative_coefficients']))
            indices = row['wall_indices']
            assert len(indices) == len(coeff)
            assert all(a.sign() >= 0 for a in coeff)
            for k, s in enumerate(basis):
                assert sum((a * dot(s, normals[l]) for a, l in zip(coeff, indices)), G(0)) == target[k]
        else:
            assert row['status'] == 'EXACT_GEOMETRIC_CUTTER'
            r = tuple(map(G.from_json, row['normal']))
            b = tuple(map(G.from_json, row['slice_coefficients']))
            assert r == tuple(sum((b[k] * basis[k][l] for k in range(4)), G(0)) for l in range(10))
            assert dot(r, mv(gram, r)) == G.from_json(row['norm'])
            assert dot(r, mv(gram, r)).sign() > 0
            assert all(dot(r, n).sign() <= 0 for n in normals)
            assert dot(r, target_normal) == G.from_json(row['escape_pairing'])
            assert dot(r, target_normal).sign() > 0
            inside = tuple(map(G.from_json, source['cases'][0]['strict_interior']))
            assert dot(r, mv(gram, inside)).sign() < 0
            q = dot(r, mv(gram, r))
            assert any(not (2*x/q).is_integral for x in mv(gram, r))
    assert seen == set(combinations(range(8), 2))
    assert sum(row['status'] == 'EXACT_SLICE_OBSTRUCTION' for row in payload['slices']) == 27
    return {'status': 'PASS_EXACT_28_SLICE_REPLAY', 'closed_slices': 27,
            'open_geometric_slices': 1, 'new_reflective_walls': 0}


def produce():
    import numpy as np
    from scipy.optimize import linprog

    source, gram, walls, v, omega, time, scalar, normals, target_normal = setup()
    real = lambda x: float(x.a) + float(x.b) * (1 + 5**.5) / 2
    rows = []
    for i, j in combinations(range(8), 2):
        basis = (time, omega[i], omega[j], scalar)
        columns = tuple(tuple(dot(s, n) for s in basis) for n in normals)
        target = tuple(dot(s, target_normal) for s in basis)
        prop = linprog(np.zeros(len(walls)), A_eq=np.array([[real(x) for x in col] for col in columns]).T,
                       b_eq=np.array([real(x) for x in target]), bounds=(0, None), method='highs')
        if prop.success:
            ids = tuple(map(int, np.flatnonzero(prop.x > 1e-7)))
            coeff = solve_support(tuple(columns[k] for k in ids), target)
            rows.append({'fundamental_weight_indices': [i, j], 'status': 'EXACT_SLICE_OBSTRUCTION',
                         'wall_indices': list(ids), 'nonnegative_coefficients': enc(coeff)})
        else:
            assert (i, j) == (4, 7)
            coeff = (G(4, 3), G(2), G(0, 2), G(1, 1))
            r = tuple(sum((coeff[k]*basis[k][l] for k in range(4)), G(0)) for l in range(10))
            rows.append({'fundamental_weight_indices': [i, j], 'status': 'EXACT_GEOMETRIC_CUTTER',
                         'slice_coefficients': enc(coeff), 'normal': enc(r),
                         'norm': dot(r, mv(gram, r)).to_json(),
                         'escape_pairing': dot(r, target_normal).to_json(),
                         'reflective': False})
    payload = {'source': SOURCE, 'source_case': 0, 'wall_count': len(walls),
               'target_escape': enc(v), 'slices': rows,
               'scope': 'All 28 spans of time, scalar, and two standard E8 fundamental weights, for this one 382-wall cone and its saved escape. Closed slices exclude any compatible real normal in that slice. The open slice has a geometric cutter, not yet a reflective wall.'}
    payload['replay'] = check(payload)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    print(json.dumps(payload['replay']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(check(json.loads(OUTPUT.read_text()))))
    else:
        produce()
