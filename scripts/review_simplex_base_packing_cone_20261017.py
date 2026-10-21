#!/usr/bin/env python3
"""Independent exact-rank review of the complete <=5-vertex base cones."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'review_simplex_base_packing_cone_20261017'
INPUT = 'fresh_simplex_base_packing_cone_20261017'


def matrix_rank(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    if not a:
        return 0
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [x / scale for x in a[row]]
        for i in range(row + 1, len(a)):
            if a[i][col]:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def set_partitions(n):
    def visit(next_root, blocks):
        if next_root == n:
            yield tuple(blocks)
            return
        bit = 1 << next_root
        for i in range(len(blocks)):
            new = list(blocks); new[i] |= bit
            yield from visit(next_root + 1, new)
        yield from visit(next_root + 1, blocks + [bit])
    return list(visit(0, []))


def run():
    path = ROOT / 'results' / (INPUT + '.json')
    saved = json.loads(path.read_text())
    assert saved['status'] == 'PASS_EXACT_COMPLETE_SIMPLEX_BASE_PACKING_CONES'
    pair_masks = [m for m in range(32) if m.bit_count() == 2]
    disjoint = [(i, j) for i, j in combinations(range(10), 2) if pair_masks[i] & pair_masks[j] == 0]
    ranks = Counter(); exact_vertices = []; feasible = 0
    for p in product((0, 1, 2), repeat=10):
        if any(p[i] + p[j] > 2 for i, j in disjoint):
            continue
        feasible += 1
        equations = []
        for i, v in enumerate(p):
            if v in (0, 2):
                equations.append([int(i == j) for j in range(10)])
        for i, j in disjoint:
            if p[i] + p[j] == 2:
                equations.append([int(k in (i, j)) for k in range(10)])
        rank = matrix_rank(equations); ranks[rank] += 1
        if rank == 10:
            exact_vertices.append(p)
    assert feasible == 1959 and len(exact_vertices) == 294
    histogram = Counter(p.count(1) for p in exact_vertices)
    assert sorted(map(list, histogram.items())) == saved['petersen_vertex_half_coordinate_histogram']
    five = next(c for c in saved['cones'] if c['vertices'] == 5)
    pair_positions = [five['base_masks'].index(m) for m in pair_masks]
    projected_saved = {tuple(row[i] for i in pair_positions) for row in five['doubled_inequality_coefficients']}
    assert set(exact_vertices) == projected_saved
    cone_checks = []
    for cone in saved['cones']:
        n = cone['vertices']; masks = cone['base_masks']; coefficients = cone['doubled_inequality_coefficients']
        assert masks == [m for m in range(1 << n) if m.bit_count() >= 2]
        partitions = set_partitions(n)
        packings = {tuple(sorted(b for b in part if b.bit_count() >= 2)) for part in partitions}
        assert len(packings) == len(partitions)
        assert packings == {tuple(sorted(p)) for p in cone['integer_packings']}
        assert all(not (a & b) for p in packings for a, b in combinations(p, 2))
        incidence = [[int(m in p) for m in masks] for p in sorted(packings)]
        # Origin and all singleton-block vectors give the full ambient dimension.
        assert matrix_rank(incidence) == len(masks)
        facet_ranks = []
        for co in coefficients:
            assert len(co) == len(masks) and all(x in (0, 1, 2) for x in co)
            loads = [sum(x * y for x, y in zip(co, p)) for p in incidence]
            assert max(loads) == cone['doubled_capacity'] == 2
            tight = [p for p, val in zip(incidence, loads) if val == 2]
            origin = tight[0]
            affine_rank = matrix_rank([[x-y for x, y in zip(p, origin)] for p in tight[1:]])
            assert affine_rank == len(masks) - 1
            facet_ranks.append(affine_rank)
        assert facet_ranks == cone['facet_affine_ranks']
        if n == 5:
            for row in coefficients:
                values = dict(zip(masks, row))
                assert all(values[m] + values[31 ^ m] == 2 for m in pair_masks)
                assert all(values[m] == 2 for m in masks if m.bit_count() >= 4)
        cone_checks.append(dict(vertices=n, base_types=len(masks), packings=len(packings),
                               facets=len(coefficients), facet_affine_ranks=facet_ranks))
    witness = saved['witness']; blocks = witness['base_masks']
    assert all(a & b for a, b in combinations(blocks, 2))
    checks = []
    for T in range(1, 32):
        for r in range(1, T.bit_count() + 1):
            load = sum(Fraction((T & B).bit_count() // r, 3) for B in blocks)
            assert load <= T.bit_count() // r
            checks.append([T, r, str(load), T.bit_count() // r])
    assert checks == witness['subset_floor_checks']
    paths = [f'scripts/{TAG}.py', f'scripts/{INPUT}.py', f'results/{INPUT}.json']
    return dict(status='PASS_INDEPENDENT_EXACT_BASE_PACKING_CONE_REVIEW',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                arithmetic='stdlib Fraction only; no producer arithmetic imports',
                feasible_half_grid_points=feasible, vertices_by_active_matrix_rank=len(exact_vertices),
                active_rank_histogram=sorted(map(list, ranks.items())),
                half_coordinate_histogram=sorted(map(list, histogram.items())),
                cone_checks=cone_checks, fractional_witness_all_subset_checks=len(checks),
                scope='Independently reconstructs all feasible half-grid active ranks, all subset packings via set partitions, and every saved facet affine rank. Completeness beyond the grid uses the separately reviewed half-integrality/LP-duality proof; this does not assert geometric realizability or an ambient dimension bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = run(); output = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(result['status'], result['feasible_half_grid_points'], result['vertices_by_active_matrix_rank'], flush=True)
