#!/usr/bin/env python3
"""Exact full-group nonreflectivity witness for the older silver F4^2 host.

Every run also repeats the complete old height30 producer in a temporary
directory and compares all mathematical fields with the archived prefix.
This is a single symmetry certificate, not a complete finite-corner census.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from hashlib import sha256

import fresh_frontier_construction_20261015 as A
import fresh_frontier_construction_20261015_geometry as B
import fresh_frontier_construction_20261015_symmetry as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_silver_full_chamber_20261015'
SOURCE = 'results/fresh_silver_f4_h9_unit_h30_20261013.json'
PRODUCER = 'scripts/fresh_silver_f4_h9_20261013.py'
SRC = [0, 1, 2, 5, 6, 10, 12, 56, 75]
DST = [0, 1, 2, 57, 17, 5, 35, 6, 26]
K = A.K


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def pd(matrix):
    m = [r[:] for r in matrix]
    out = []
    for k in range(len(m)):
        pivot = m[k][k]
        assert pivot.sign() > 0
        out.append(pivot)
        for i in range(k + 1, len(m)):
            for j in range(k + 1, len(m)):
                m[i][j] -= m[i][k] * m[k][j] / pivot
    return out


def run():
    d = json.loads((ROOT / SOURCE).read_text())
    with tempfile.TemporaryDirectory(prefix='silver_full_chamber_replay_') as tmp:
        output = Path(tmp) / 'prefix.json'
        cmd = [sys.executable, str(ROOT / PRODUCER), '--alpha', '1,1',
               '--height', '30', '--output', str(output)]
        env = dict(os.environ, PYTHONPATH=str(ROOT / 'src') + os.pathsep + str(ROOT))
        subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, check=True)
        replay = json.loads(output.read_text())
    timing = {'seconds', 'root_seconds'}
    assert {k: v for k, v in d.items() if k not in timing} == {
        k: v for k, v in replay.items() if k not in timing}
    assert d['completed'] and A.dec(d['completed_height']) == 30
    G = [[A.dec(x) for x in row] for row in d['gram']]
    N = [list(map(A.dec, row)) for row in d['all_walls']]
    assert len(N) == 92
    H = [K(1)] + [K()] * 9
    assert A.pair(H, H, G) == -1 - A.S
    assert C.determinant(G) == -1 - A.S
    positive_pivots = pd([r[1:] for r in G[1:]])
    conjugate_pivots = pd([[x.conjugate() for x in row] for row in G])
    assert list(map(A.dec, d['norm_pool'])) == [K(1), 2 + A.S, K(2)]
    covectors = [A.mv(G, n) for n in N]
    gram = [[A.dot(n, c) for c in covectors] for n in N]
    assert all(gram[i][j].sign() <= 0 for i in range(92) for j in range(i))
    assert all((2 * x / gram[i][i]).integral()
               for i in range(92) for x in covectors[i])
    angles = Counter()
    allowed = {K(0): 2, K(1): 3, K(2): 4, K(3): 6, 2 + A.S: 8}
    for i in range(92):
        for j in range(i):
            q = 4 * gram[i][j] * gram[i][j] / gram[i][i] / gram[j][j]
            if (q - 4).sign() < 0:
                assert q in allowed
                angles['m' + str(allowed[q])] += 1
            else:
                assert (q - 4).sign() > 0
                angles['dotted'] += 1
    center = list(map(A.dec, d['strict_interior']))
    assert A.pair(center, center, G).sign() < 0
    assert A.pair(center, H, G).sign() < 0
    assert all(A.dot(center, c).sign() < 0 for c in covectors)
    facet_points = []
    for i, n in enumerate(N):
        coefficient = A.dot(center, covectors[i]) / gram[i][i]
        v = [x - coefficient * y for x, y in zip(center, n)]
        assert A.pair(v, v, G).sign() < 0 and A.pair(v, H, G).sign() < 0
        values = [A.dot(v, c) for c in covectors]
        assert all(value.sign() == (0 if i == j else -1) for j, value in enumerate(values))
        facet_points.append(v)

    common = [[gram[i][j] for j in SRC] for i in SRC]
    assert common == [[gram[i][j] for j in DST] for i in DST]
    pivots = pd(common)
    inverse = C.inverse(common)
    assert all(x.sign() >= 0 for row in inverse for x in row)
    corners = []
    for ids in (SRC, DST):
        v = B.null([covectors[i] for i in ids])
        assert v
        if A.pair(v, H, G).sign() > 0:
            v = [-x for x in v]
        square = A.pair(v, v, G)
        hv = A.pair(H, v, G)
        assert square.sign() < 0 and hv.sign() < 0
        values = [A.dot(v, c) for c in covectors]
        assert [i for i, x in enumerate(values) if not x] == sorted(ids)
        assert all(x.sign() <= 0 for x in values)
        a = hv / square
        b = A.mv(inverse, [A.dot(H, covectors[i]) for i in ids])
        assert a.sign() > 0 and all(x.sign() <= 0 for x in b)
        assert H == [a * v[k] + sum(b[j] * N[ids[j]][k] for j in range(9)) for k in range(10)]
        corners.append(dict(ordered_indices=ids, point=v, norm=square,
                            basepoint_pairing=hv, wall_pairings=values,
                            basepoint_decomposition=dict(vertex_coefficient=a, normal_coefficients=b)))
    F = C.transpose([N[i] for i in SRC] + [corners[0]['point']])
    D = C.transpose([N[i] for i in DST] + [corners[1]['point']])
    scale = C.determinant(F) / C.determinant(D)
    if scale.sign() < 0:
        scale = -scale
    assert scale.sign() > 0
    assert corners[0]['norm'] == scale * scale * corners[1]['norm']
    for row in D:
        row[-1] *= scale
    matrix = C.mm(D, C.inverse(F))
    inverse_matrix = C.inverse(matrix)
    assert all(x.integral() for row in matrix + inverse_matrix for x in row)
    identity = [[K(i == j) for j in range(10)] for i in range(10)]
    assert C.mm(matrix, inverse_matrix) == C.mm(inverse_matrix, matrix) == identity
    assert C.mm(C.mm(C.transpose(matrix), G), matrix) == G
    assert all(A.mv(matrix, N[i]) == N[j] for i, j in zip(SRC, DST))
    assert A.mv(matrix, corners[0]['point']) == [scale * x for x in corners[1]['point']]
    trace = sum(matrix[i][i] for i in range(10))
    assert trace == 8 + 3 * A.S and (trace - 10).sign() > 0
    sources = [f'scripts/{TAG}.py', SOURCE, PRODUCER,
               'scripts/fresh_frontier_construction_20261015.py',
               'scripts/fresh_frontier_construction_20261015_geometry.py',
               'scripts/fresh_frontier_construction_20261015_symmetry.py',
               'src/coxeter_polytopes/exact/silver.py',
               'src/coxeter_polytopes/search/silver_path_host.py',
               'results/fresh_silver_f4_h9_report_20261013.md',
               'results/fresh_silver_f4_independent_audit_20261013.json',
               'results/fresh_silver_f4_independent_report_20261013.md',
               'results/fresh_silver_f4_identifications_20261013.json',
               'scripts/fresh_silver_f4_identifications_20261013.py']
    return A.enc(dict(
        status='PASS_EXACT_INFINITE_ORDER_FULL_SILVER_F4_SQUARED_CHAMBER_SYMMETRY',
        sources={path: sha(path) for path in sources},
        host='Z[sqrt2]-lattice <-1-sqrt2> + F4 + F4 + <1>',
        complete_prefix_replay=dict(height=30, walls=92, control_levels=len(d['levels']),
                                    timing_fields_excluded=sorted(timing), all_mathematical_fields_equal=True),
        primitive_norm_class_representatives=[K(1), 2 + A.S, K(2)],
        physical_positive_block_pivots=positive_pivots,
        conjugate_form_pivots=conjugate_pivots, form=G,
        finite_angle_counts=dict(sorted(angles.items())), strict_center=center,
        actual_facet_points=facet_points, source_corner=corners[0], target_corner=corners[1],
        common_corner_type='A1 + B4 + I2(4) + I2(4)', common_corner_gram=common,
        common_corner_pivots=pivots, inverse_corner_gram=inverse,
        positive_corner_scale=scale, matrix=matrix, inverse_matrix=inverse_matrix,
        determinant=C.determinant(matrix), trace=trace, trace_minus_ten=trace - 10,
        conclusion='The full integral reflection chamber of this literal lattice has infinite volume; the lattice is nonreflective.',
        scope='Single exact infinite-order stabilizer witness. Full norm alphabet and complete zero sector are inherited from the bound old package, and the entire bounded producer is replayed. Complete outward corners persist under later simple walls by the saved decompositions. No complete finite-corner census, dimension exclusion, or claim about other commensurable lattices. Independent new review pending.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], 'trace=8+3sqrt2')
