#!/usr/bin/env python3
"""Exact local data and finite arithmetic for the restricted RA5 block theorem.

The universal geometric and graph arguments are in the companion report.
This does not enumerate hyperbolic five-polytopes or prove a dimension bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fresh_ideal_ra5_octahedral_blocks_20261018.json'
PINS = {
    'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json':
    'cf5e71eb3720d16070852c19d3274f0a62b052ed38035d3f9d51adce666deb1c',
    'results/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.json':
    '59c12fc327210d80223030ae20428fed127adceee4ebed1ead36bd7efb0e3c2f',
}


def dot(a, b):
    return -a[0]*b[0] + sum(x*y for x, y in zip(a[1:], b[1:]))


def build():
    for path, digest in PINS.items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == digest
    data = json.loads((ROOT/next(iter(PINS))).read_text())['stages'][-1]
    normals = [tuple(map(Q, n)) for n in data['walls']]
    vertices = [tuple(map(Q, v)) for v in data['vertices']]
    faces = [set(i for i, v in enumerate(vertices) if dot(n, v) == 0)
             for n in normals]
    assert [sorted(f) for f in faces] == data['facet_vertices']
    large = [i for i, f in enumerate(faces) if len(f) == 24]
    assert large == [5, 14, 30, 37]
    gram = [[dot(normals[i], normals[j]) for j in large] for i in large]
    assert gram == [[1, 0, 0, -1], [0, 1, -1, 0],
                    [0, -1, 1, 0], [-1, 0, 0, 1]]
    pair_sums = [tuple(x+y for x, y in zip(normals[i], normals[j]))
                 for i, j in ((5, 37), (14, 30))]
    assert all(s[0] > 0 and dot(s, s) == 0 for s in pair_sums)
    assert tuple(x/pair_sums[0][0] for x in pair_sums[0]) == vertices[85]
    assert tuple(x/pair_sums[1][0] for x in pair_sums[1]) == vertices[85]
    ridges = []
    for i in large:
        rows = [{'neighbor': j, 'polygon_vertices': len(faces[i] & faces[j]),
                 'neighbor_ridge_vertices': len(faces[j])}
                for j in range(len(normals)) if dot(normals[i], normals[j]) == 0]
        assert Counter(r['polygon_vertices'] for r in rows) == {3: 12, 4: 12, 6: 2}
        assert all((r['polygon_vertices'] == 6) ==
                   (r['neighbor_ridge_vertices'] == 24) for r in rows)
        ridges.append({'ridge': i, 'faces': rows})

    # Exhaust all proper colorings of an induced K_{2,2,2} by four colors.
    edges = [(i, j) for i, j in combinations(range(6), 2) if i//2 != j//2]
    opposite = ((0, 1), (2, 3), (4, 5))
    color_cases = []
    for sizes in ((6, 8, 8, 8), (7, 7, 8, 8)):
        reasons = Counter()
        survivors = Counter()
        common_bounds = set()
        for colors in product(range(4), repeat=6):
            if any(colors[i] == colors[j] for i, j in edges):
                continue
            used = set(colors)
            assert len(used) in (3, 4)
            if len(used) == 3:
                missing, = set(range(4))-used
                bound = sizes[missing] if sizes[missing] == 8 else (
                    sizes[missing]-sum(8-sizes[c] for c in colors))
                assert bound > 2
                common_bounds.add(bound)
                reasons['too_many_common_neighbors'] += 1
                continue
            split = [(i, j) for i, j in opposite if colors[i] != colors[j]]
            assert len(split) == 1
            i, j = split[0]
            if sizes[colors[i]] == 8 or sizes[colors[j]] == 8:
                reasons['nonedge_incident_to_complete_color'] += 1
                continue
            count = tuple(colors.count(c) for c in range(4))
            assert sizes == (7, 7, 8, 8) and count == (1, 1, 2, 2)
            assert tuple(5*x for x in count) != sizes
            survivors[count] += 1
        color_cases.append({'part_sizes': sizes, 'excluded_colorings': dict(reasons),
                            'three_color_common_neighbor_lower_bounds': sorted(common_bounds),
                            'remaining_block_color_counts': [
                                {'counts': k, 'colorings': v} for k, v in sorted(survivors.items())]})
    assert 22*11 > 8*30  # AES strict threshold.
    assert 3*36-24 == 84
    return {'status': 'RESTRICTED_RA5_OCTAHEDRAL_BLOCK_AND_T30_EXCLUSION_DATA',
            'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
            'sources': PINS, 'ridge_indices': large,
            'gram': [[str(x) for x in row] for row in gram],
            'null_pair_sums': [[str(x) for x in s] for s in pair_sums],
            'ridges': ridges, 'four_color_cases': color_cases,
            'published_external_input': {
                'theorem': 'Andrasfai-Erdos-Sos1974 Theorem1.1, r=5',
                'url': 'https://real.mtak.hu/110445/1/1974-18.pdf',
                'specialization': 'K5-free and minimum degree >8n/11 imply four-colorable'},
            'conclusions_with_companion_geometric_proof': {
                'T_multiple': 6, 'T_minimum': 36, 'double24_minimum': 84,
                'spectrum_excluded': False, 'new_dimension_bound': False,
                'new_geometric_example': False, 'Lean_verified': False},
            'scope': 'Exact inherited T geometry and exhaustive4096-color finite arithmetic; '
                     'universal block, Lorentz, projection and counting proofs are written in companion report.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.verify:
        assert json.loads(OUT.read_text()) == json.loads(json.dumps(result))
    else:
        assert not OUT.exists()
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'verify': args.verify, 'ridges': 4,
                      'polygon_faces': 104, 'colorings_per_partition': 4096,
                      'T_minimum': 36, 'double24_minimum': 84,
                      'new_dimension_bound': False}))
