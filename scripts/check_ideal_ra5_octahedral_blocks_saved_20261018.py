#!/usr/bin/env python3
"""Independent exact reconstruction of all104 polygons in the four T ridges.

No producer imports, search, optimizer, Lean job, or five-polytope enumeration.
Universal geometric arguments remain in the companion written proof.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import file_digest
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json'
CERTIFICATE = ROOT / 'results/fresh_ideal_ra5_octahedral_blocks_20261018.json'
PINS = {
    'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json':
        'cf5e71eb3720d16070852c19d3274f0a62b052ed38035d3f9d51adce666deb1c',
    'results/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.json':
        '59c12fc327210d80223030ae20428fed127adceee4ebed1ead36bd7efb0e3c2f',
}


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def dot(a, b):
    return sum(((-1 if i == 0 else 1)*x*y for i, (x, y) in enumerate(zip(a, b))), Q())


def rank(rows):
    a = [list(row) for row in rows]
    pivot = 0
    for column in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][column]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        factor = a[pivot][column]
        a[pivot] = [x/factor for x in a[pivot]]
        for i in range(len(a)):
            if i != pivot:
                factor = a[i][column]
                a[i] = [x-factor*y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


data = json.loads(INPUT.read_text())['stages'][-1]
certificate = json.loads(CERTIFICATE.read_text())
assert certificate['sources'] == PINS
for path, expected in PINS.items():
    assert digest(ROOT/path) == expected
assert certificate['source_sha256'] == digest(
    ROOT/'scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py')
ns = [tuple(map(Q, row)) for row in data['walls']]
vs = [tuple(map(Q, row)) for row in data['vertices']]
assert len(ns) == 198 and len(vs) == 258
fv = [{j for j, v in enumerate(vs) if dot(n, v) == 0} for n in ns]
assert [sorted(s) for s in fv] == data['facet_vertices']
large = [i for i, vertices in enumerate(fv) if len(vertices) == 24]
assert large == certificate['ridge_indices'] == [5, 14, 30, 37]
assert rank([ns[i] for i in large]) == 3
assert [[str(dot(ns[i], ns[j])) for j in large] for i in large] == certificate['gram']
for i, j in ((5, 37), (14, 30)):
    summed = tuple(x+y for x, y in zip(ns[i], ns[j]))
    assert summed[0] > 0 and dot(summed, summed) == 0
    assert tuple(x/summed[0] for x in summed) == vs[85]
records = []
for i in large:
    polygons = []
    for j in range(len(ns)):
        if j == i:
            continue
        common = fv[i] & fv[j]
        adjacent = dot(ns[i], ns[j]) == 0
        # Exclude accidental counts from nonadjacent, asymptotic wall pairs.
        assert (len(common) >= 3) == adjacent
        if adjacent:
            assert rank([vs[k] for k in sorted(common)]) == 3
            polygons.append({'neighbor': j, 'polygon_vertices': len(common),
                             'neighbor_ridge_vertices': len(fv[j])})
    assert Counter(row['polygon_vertices'] for row in polygons) == {3: 12, 4: 12, 6: 2}
    assert {row['neighbor'] for row in polygons if row['polygon_vertices'] == 6} == {
        j for j in large if j != i and dot(ns[i], ns[j]) == 0}
    records.append({'ridge': i, 'faces': polygons})
assert records == certificate['ridges']
print(json.dumps({'status': 'PASS_INDEPENDENT_104_POLYGON_RECONSTRUCTION',
                  'ridges': 4, 'polygons': 104, 'all_pair_intersections_checked': 4*197,
                  'geometric_block_proof_machine_verified': False,
                  'dimension_bound_improved': False}))
