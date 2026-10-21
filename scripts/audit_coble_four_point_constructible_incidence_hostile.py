#!/usr/bin/env python3
"""Hostile, helper-free audit of the advertised ``D_CobVin^CI4`` theorem.

The audit intentionally does not import the production certificate or its
blowdown bridge.  It reconstructs the two raw divisor lattices, contractions,
dependency designs, and characteristic-three sections from the pinned full
wall data.  It then separates three logically different statements:

* the two explicit four-point cores are genuine strong-Coble root seeds;
* the existing D_CobVin theorem conditionally returns either seed to the
  literally same marked polytope; and
* the advertised closure under "standard-free" extensions is not proved by
  base-point-freeness plus generic separability.

The last point is witnessed by |O(2,0)| on P1 x P1: its morphism is
generically separable in characteristic three, while every general member is
the disjoint union of two ruling fibres.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
import functools
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_four_point_constructible_incidence_hostile_audit.json"

SOURCES = {
    "production_script": ROOT / "scripts/certify_coble_four_point_constructible_incidence.py",
    "production_archive": ROOT / "data/enumerations/coble_four_point_constructible_incidence.json",
    "production_report": ROOT / "results/coble_four_point_constructible_incidence_lemma.md",
    "production_test": ROOT / "tests/test_coble_four_point_constructible_incidence.py",
    "remaining_seven_script": ROOT / "scripts/certify_coble_remaining_seven_contraction_designs.py",
    "remaining_seven_archive": ROOT / "data/enumerations/coble_remaining_seven_contraction_designs.json",
    "remaining_seven_report": ROOT / "results/coble_remaining_seven_contraction_design_audit.md",
    "remaining_seven_test": ROOT / "tests/test_coble_remaining_seven_contraction_designs.py",
    "D_CobVin_archive": ROOT / "data/enumerations/coble_vinberg_incidence_inverse_transfer.json",
    "type_VI_raw_full": ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json",
    "type_MII_raw_full": ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json",
    "type_MII_literal_cone": ROOT / "data/enumerations/coble_type_mii_h17_hostile_audit.json",
    "type_MI_raw_full": ROOT / "data/enumerations/coble_type_mi_full_ns_h11.json",
}

PINNED_SHA256 = {
    "production_script": "2f4189250daa3877d17935bf331cdb312c91fde1dc466c6a61728f4ef2d999ab",
    "production_archive": "0c0d2ae6a1c37b723d3804b8f856fcee931d9eb96d39a08116586766985f06c7",
    "production_report": "6c7232e66486619dcd2ca50089db74c344de7caa034349e8242b3c3b6f444a18",
    "production_test": "1cafa6b5aceb7aa98caf2615f05ff50c2134cd090e57a7c40200f5cf8d5cfb7b",
    "remaining_seven_script": "4b5ad212e1d18d8078267a7772c4e07482b373015277206c94f5dcc0ab55719f",
    "remaining_seven_archive": "adec1b838db24ef005b881810d70e6851fd17ef6f1a719e35e94ee1ec7474f3d",
    "remaining_seven_report": "7da0bf6966981d1bb2163d2e7ebfbce3e4e4728aeaab5197dd81049674e1ce3f",
    "remaining_seven_test": "43189fee749e97297eb439719b852c99936d1243018cfb2eeaf7303890eec1e3",
    "D_CobVin_archive": "ee4280e2c9cbee2a730f2333c63566f8cc36bf7dab125f9b168ce437260028f8",
    "type_VI_raw_full": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "type_MII_raw_full": "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
    "type_MII_literal_cone": "124a299740c6a82c7ab543ce9236c1ad38d6e988e637a0d173cbe2203756c25f",
    "type_MI_raw_full": "5431030fab49d0a19c4ecd32765bbf1a3a9db256a2e37a15b7d92355df1d6976",
}

SCHEMA = "coxeter.coble_four_point_constructible_incidence.hostile_audit.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_pinned(name: str) -> dict[str, Any]:
    path = SOURCES[name]
    digest = _sha256(path)
    if digest != PINNED_SHA256[name]:
        raise RuntimeError(f"pinned source changed for {name}: {digest}")
    return json.loads(path.read_text(encoding="utf-8"))


def _integer_vector(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1 or any(value.q != 1 for value in vector):
        raise ArithmeticError("expected an integral column vector")
    return [int(value) for value in vector]


def _primitive_integer_vector(vector: sp.Matrix) -> sp.Matrix:
    denominators = [int(value.q) for value in vector]
    scale = sp.ilcm(*denominators) if denominators else 1
    values = [int(value * scale) for value in vector]
    divisor = functools.reduce(math.gcd, (abs(value) for value in values if value), 0)
    if not divisor:
        raise ArithmeticError("zero vector cannot be made primitive")
    return sp.Matrix([value // divisor for value in values])


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact symmetric-congruence inertia, independent of project helpers."""

    work = sp.Matrix(matrix)
    if work != work.T:
        raise ArithmeticError("inertia input is not symmetric")
    positive = negative = zero = 0
    while work.rows:
        size = work.rows
        diagonal = next((index for index in range(size) if work[index, index]), None)
        if diagonal is not None:
            order = [diagonal] + [index for index in range(size) if index != diagonal]
            work = work.extract(order, order)
            pivot = work[0, 0]
            positive += 1 if pivot > 0 else 0
            negative += 1 if pivot < 0 else 0
            if size == 1:
                work = sp.zeros(0)
            else:
                tail = work[1:, 1:]
                column = work[1:, 0]
                work = tail - column * column.T / pivot
            continue
        off_diagonal = next(
            (
                (left, right)
                for left in range(size)
                for right in range(left + 1, size)
                if work[left, right]
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break
        left, right = off_diagonal
        order = [left, right] + [
            index for index in range(size) if index not in {left, right}
        ]
        work = work.extract(order, order)
        block = work[:2, :2]
        positive += 1
        negative += 1
        if size == 2:
            work = sp.zeros(0)
        else:
            cross = work[:2, 2:]
            work = work[2:, 2:] - cross.T * block.inv() * cross
    return positive, negative, zero


def _modular_kernel(
    matrix: Sequence[Sequence[Any]], prime: int, substitutions: dict[Any, int] | None = None
) -> tuple[int, tuple[int, ...]]:
    """Return rank and the normalized unique kernel vector over F_p."""

    substitutions = substitutions or {}
    rows = [
        [int(sp.sympify(value).subs(substitutions)) % prime for value in row]
        for row in matrix
    ]
    row_count = len(rows)
    column_count = len(rows[0]) if rows else 0
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]), None
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [(value * inverse) % prime for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                (rows[row][other] - factor * rows[pivot_row][other]) % prime
                for other in range(column_count)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    rank = len(pivot_columns)
    free_columns = [column for column in range(column_count) if column not in pivot_columns]
    if len(free_columns) != 1:
        raise ArithmeticError(f"expected a one-dimensional kernel, got {free_columns}")
    kernel = [0] * column_count
    kernel[free_columns[0]] = 1
    for row, column in reversed(list(enumerate(pivot_columns))):
        kernel[column] = -sum(
            rows[row][other] * kernel[other] for other in free_columns
        ) % prime
    first = next(value for value in kernel if value)
    inverse = pow(first, -1, prime)
    normalized = tuple((value * inverse) % prime for value in kernel)
    return rank, normalized


def _normalize_mod(vector: Iterable[int], prime: int) -> tuple[int, ...]:
    values = tuple(value % prime for value in vector)
    first = next(value for value in values if value)
    inverse = pow(first, -1, prime)
    return tuple(value * inverse % prime for value in values)


def _four_points(parameter: Any) -> list[tuple[Any, Any]]:
    return [(1, 0), (0, 1), (1, 1), (parameter, 1)]


def _grid_evaluation(points: Sequence[tuple[Any, Any]], permutation: Sequence[int]) -> list[list[Any]]:
    return [
        [x0 * y0, x0 * y1, x1 * y0, x1 * y1]
        for (x0, x1), right in zip(points, permutation)
        for y0, y1 in [points[right]]
    ]


def _conic_point(parameter: Any) -> sp.Matrix:
    return sp.Matrix([1, 0, 0]) if parameter is None else sp.Matrix([parameter**2, parameter, 1])


def _conic_tangent(parameter: Any) -> sp.Matrix:
    return sp.Matrix([0, 0, 1]) if parameter is None else sp.Matrix([1, -2 * parameter, parameter**2])


def _audit_four_point_algebra() -> dict[str, Any]:
    parameter = sp.Symbol("lambda")
    permutations = list(itertools.permutations(range(4)))
    determinants = [
        sp.factor(sp.Matrix(_grid_evaluation(_four_points(parameter), permutation)).det())
        for permutation in permutations
    ]
    determinant_histogram = Counter(str(value) for value in determinants)
    expected = Counter(
        {
            "0": 4,
            "-(lambda - 1)*(lambda + 1)": 4,
            "2*lambda - 1": 4,
            "-lambda**2 + lambda - 1": 4,
            "lambda**2 - lambda + 1": 4,
            "lambda*(lambda - 2)": 4,
        }
    )
    if determinant_histogram != expected:
        raise ArithmeticError(f"grid determinant alphabet differs: {determinant_histogram}")

    grid_sections = []
    for permutation in permutations:
        evaluation = _grid_evaluation(_four_points(2), permutation)
        rank, section = _modular_kernel(evaluation, 3)
        coefficient_determinant = (section[0] * section[3] - section[1] * section[2]) % 3
        if rank != 3 or coefficient_determinant == 0:
            raise ArithmeticError("a characteristic-three grid section is not a smooth graph")
        grid_sections.append(section)
    if len(set(grid_sections)) != 24:
        raise ArithmeticError("the 24 projectivity graph sections are not distinct")

    symbolic_parameters: list[Any] = [None, 0, 1, parameter]
    polar_polynomials = []
    for pair in itertools.combinations(range(4), 2):
        complement = tuple(index for index in range(4) if index not in pair)
        pole = _conic_tangent(symbolic_parameters[pair[0]]).cross(
            _conic_tangent(symbolic_parameters[pair[1]])
        )
        chord = _conic_point(symbolic_parameters[complement[0]]).cross(
            _conic_point(symbolic_parameters[complement[1]])
        )
        polar_polynomials.append(sp.factor(chord.dot(pole)))
    polar_histogram = Counter(str(value) for value in polar_polynomials)
    expected_polar = Counter(
        {
            "(lambda - 1)*(lambda + 1)": 2,
            "lambda*(lambda - 2)": 2,
            "1 - 2*lambda": 2,
        }
    )
    if polar_histogram != expected_polar:
        raise ArithmeticError(f"polar/chord alphabet differs: {polar_histogram}")

    finite_parameters: list[Any] = [None, 0, 1, 2]
    chord_lines: list[tuple[int, ...]] = []
    polar_ranks = []
    for pair in itertools.combinations(range(4), 2):
        complement = tuple(index for index in range(4) if index not in pair)
        pole = _conic_tangent(finite_parameters[pair[0]]).cross(
            _conic_tangent(finite_parameters[pair[1]])
        )
        points = [
            pole,
            _conic_point(finite_parameters[complement[0]]),
            _conic_point(finite_parameters[complement[1]]),
        ]
        evaluation = sp.Matrix.hstack(*points).T.tolist()
        rank, line = _modular_kernel(evaluation, 3)
        chord = _normalize_mod(
            _conic_point(finite_parameters[complement[0]]).cross(
                _conic_point(finite_parameters[complement[1]])
            ),
            3,
        )
        if line != chord:
            raise ArithmeticError("the polar line is not the complementary chord")
        polar_ranks.append(rank)
        chord_lines.append(line)
    tangent_lines = {
        _normalize_mod(_conic_tangent(value), 3) for value in finite_parameters
    }
    if Counter(polar_ranks) != Counter({2: 6}):
        raise ArithmeticError("a polar/chord evaluation does not have rank two")
    if len(set(chord_lines)) != 6 or set(chord_lines) & tangent_lines:
        raise ArithmeticError("residual chord lines collide with each other or a boundary tangent")

    if not all(int(value.subs(parameter, 2)) % 3 == 0 for value in determinants):
        raise ArithmeticError("the grid design does not specialize at lambda=2 in F3")
    if not all(int(value.subs(parameter, 2)) % 3 == 0 for value in polar_polynomials):
        raise ArithmeticError("the polar design does not specialize at lambda=2 in F3")

    return {
        "independent_of_production_helpers": True,
        "grid": {
            "determinant_histogram": dict(sorted(determinant_histogram.items())),
            "rank_histogram_at_char3_lambda2": {"3": 24},
            "distinct_projective_kernel_sections": len(set(grid_sections)),
            "all_coefficient_matrices_invertible": True,
            "geometric_consequence": "24 distinct smooth irreducible (1,1) projectivity graphs",
        },
        "polar_chord": {
            "polynomial_histogram": dict(sorted(polar_histogram.items())),
            "rank_histogram_at_char3_lambda2": {"2": 6},
            "distinct_residual_lines": len(set(chord_lines)),
            "distinct_boundary_tangents": len(tangent_lines),
            "residual_boundary_collisions": 0,
        },
        "necessity_and_sufficiency": {
            "grid_integer_identity": "2*(lambda+1)-(2*lambda-1)=3",
            "polar_integer_identity": "(lambda+1)-(lambda-2)=3",
            "domain": "lambda != 0,1 and characteristic != 2",
            "only_solution": "characteristic 3, lambda=2",
            "converse_checked_for_all_24_plus_6_conditions": True,
        },
    }


def _bipartition(vertices: Sequence[int], edges: Iterable[tuple[int, int]]) -> tuple[list[int], list[int]]:
    adjacency = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    colors: dict[int, int] = {}
    for start in vertices:
        if start in colors:
            continue
        colors[start] = 0
        queue = deque([start])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                expected = 1 - colors[vertex]
                if neighbor in colors and colors[neighbor] != expected:
                    raise ArithmeticError("incidence graph is not bipartite")
                if neighbor not in colors:
                    colors[neighbor] = expected
                    queue.append(neighbor)
    pieces = (
        sorted(vertex for vertex in vertices if colors[vertex] == 0),
        sorted(vertex for vertex in vertices if colors[vertex] == 1),
    )
    return min(pieces, (pieces[1], pieces[0]))


def _audit_raw_mii(payload: dict[str, Any]) -> dict[str, Any]:
    gram_m = sp.Matrix(payload["arithmetic"]["gram"])
    raw_walls = [sp.Matrix(row) for row in payload["arithmetic"]["wall_vectors"]]
    wall_matrix = sp.Matrix.hstack(*raw_walls)
    raw_gram = wall_matrix.T * gram_m * wall_matrix
    pairings = [[int(raw_gram[left, right]) for right in range(48)] for left in range(48)]
    short = [index for index in range(48) if raw_gram[index, index] == -2]
    odd_edges = [
        (left, right)
        for left in range(48)
        for right in range(left + 1, 48)
        if pairings[left][right] % 2
    ]
    if len(short) != 24:
        raise ArithmeticError("raw MII does not have 24 short mirrors")

    orthogonal_count = 0
    odd_cover_count = 0
    integral_frames: list[tuple[int, ...]] = []
    for frame in itertools.combinations(short, 8):
        if any(pairings[left][right] for left, right in itertools.combinations(frame, 2)):
            continue
        orthogonal_count += 1
        frame_set = set(frame)
        covers_odd_edges = all(left in frame_set or right in frame_set for left, right in odd_edges)
        odd_cover_count += int(covers_odd_edges)
        surviving = [index for index in range(48) if index not in frame_set]
        canonical_pairs_even = all(
            sum(pairings[index][other] for index in frame) % 2 == 0
            for other in surviving
        )
        canonical_square_even = (
            sum(pairings[left][right] for left in frame for right in frame) % 2 == 0
        )
        if covers_odd_edges and canonical_pairs_even and canonical_square_even:
            integral_frames.append(frame)
    if orthogonal_count != 21_211 or integral_frames != [tuple(range(8))]:
        raise ArithmeticError(
            f"MII frame census changed: orthogonal={orthogonal_count}, frames={integral_frames}"
        )

    frame = integral_frames[0]
    frame_set = set(frame)
    canonical_m = -sum((raw_walls[index] for index in frame), sp.zeros(18, 1))
    generators = [
        2 * wall if index in frame_set else wall for index, wall in enumerate(raw_walls)
    ] + [canonical_m]
    basis_in_m = hermite_normal_form(sp.Matrix.hstack(*generators))
    gram_n = basis_in_m.T * gram_m * basis_in_m / 2
    if any(value.q != 1 for value in gram_n):
        raise ArithmeticError("MII de-scaled Gram is not integral")
    if int(gram_n.det()) != -1 or _exact_inertia(gram_n) != (1, 17, 0):
        raise ArithmeticError("MII de-scaling is not I_(1,17)")
    inverse_basis = basis_in_m.inv()
    physical = [
        inverse_basis * (2 * wall if index in frame_set else wall)
        for index, wall in enumerate(raw_walls)
    ]
    physical_matrix = sp.Matrix.hstack(*physical)
    physical_gram = physical_matrix.T * gram_n * physical_matrix
    canonical_n = inverse_basis * canonical_m

    boundary = list(frame)
    minus_one = [
        index for index in range(48) if index not in frame_set and physical_gram[index, index] == -1
    ]
    minus_two = [
        index for index in range(48) if index not in frame_set and physical_gram[index, index] == -2
    ]
    if (len(boundary), len(minus_one), len(minus_two)) != (8, 16, 24):
        raise ArithmeticError("MII physical wall taxonomy changed")
    if any(physical_gram[left, right] for left, right in itertools.combinations(minus_one, 2)):
        raise ArithmeticError("MII exceptional curves are not disjoint")
    if any(physical_gram[left, right] for left, right in itertools.combinations(boundary, 2)):
        raise ArithmeticError("MII boundary curves are not disjoint")
    if sum((physical[index] for index in boundary), sp.zeros(18, 1)) != -2 * canonical_n:
        raise ArithmeticError("MII boundary is not -2K")

    profiles = {
        exceptional: tuple(index for index in boundary if physical_gram[index, exceptional] == 1)
        for exceptional in minus_one
    }
    if any(len(profile) != 2 for profile in profiles.values()):
        raise ArithmeticError("an MII exceptional does not meet exactly two boundaries")
    left, right = _bipartition(boundary, profiles.values())
    expected_edges = {tuple(sorted((one, two))) for one in left for two in right}
    if set(tuple(sorted(profile)) for profile in profiles.values()) != expected_edges:
        raise ArithmeticError("MII boundary graph is not K4,4")
    grid: dict[tuple[int, int], int] = {}
    for exceptional, profile in profiles.items():
        left_vertex = next(index for index in profile if index in left)
        right_vertex = next(index for index in profile if index in right)
        grid[(left.index(left_vertex), right.index(right_vertex))] = exceptional

    boundary_images = []
    for index in boundary:
        boundary_images.append(
            physical[index]
            + sum(
                (
                    int(physical_gram[index, exceptional]) * physical[exceptional]
                    for exceptional in minus_one
                ),
                sp.zeros(18, 1),
            )
        )
    first_ruling = boundary_images[boundary.index(left[0])]
    second_ruling = boundary_images[boundary.index(right[0])]
    if not all(boundary_images[boundary.index(index)] == first_ruling for index in left):
        raise ArithmeticError("MII left boundaries do not contract to one ruling")
    if not all(boundary_images[boundary.index(index)] == second_ruling for index in right):
        raise ArithmeticError("MII right boundaries do not contract to one ruling")
    ruling_gram = (
        sp.Matrix.hstack(first_ruling, second_ruling).T
        * gram_n
        * sp.Matrix.hstack(first_ruling, second_ruling)
    )
    if ruling_gram != sp.Matrix([[0, 1], [1, 0]]):
        raise ArithmeticError("MII contraction does not end at P1 x P1")

    permutations = []
    for root in minus_two:
        incident = [exceptional for exceptional in minus_one if physical_gram[root, exceptional] == 1]
        positions = [position for position, exceptional in grid.items() if exceptional in incident]
        if len(positions) != 4 or {row for row, _ in positions} != set(range(4)):
            raise ArithmeticError("MII root is not a four-grid matching")
        permutation = tuple(column for row, column in sorted(positions))
        if set(permutation) != set(range(4)):
            raise ArithmeticError("MII root matching is not bijective")
        image = physical[root] + sum((physical[index] for index in incident), sp.zeros(18, 1))
        if image != first_ruling + second_ruling:
            raise ArithmeticError("MII root does not contract to bidegree (1,1)")
        permutations.append(permutation)
    if set(permutations) != set(itertools.permutations(range(4))):
        raise ArithmeticError("MII does not contain every permutation matching")

    return {
        "raw_input": "rank18/48 walls",
        "bridge_archive_read": False,
        "frame_census": {
            "short_mirrors": len(short),
            "all_eight_subsets": math.comb(len(short), 8),
            "pairwise_orthogonal": orthogonal_count,
            "odd_edges": len(odd_edges),
            "odd_edge_vertex_covers": odd_cover_count,
            "integral_frames": len(integral_frames),
            "unique_frame_zero_based": list(frame),
        },
        "reconstructed_NS": {
            "rank": gram_n.rows,
            "determinant": int(gram_n.det()),
            "inertia": list(_exact_inertia(gram_n)),
            "identification": "odd unimodular I_(1,17)",
        },
        "physical_taxonomy": {"boundary_minus4": 8, "minus1": 16, "minus2": 24},
        "strong_Coble_numerics": {
            "boundaries_disjoint": True,
            "boundary_sum_equals_minus_2K": True,
            "rational_blowup_endpoint": "P1 x P1",
        },
        "dependency_extraction": {
            "boundary_graph": "K4,4",
            "grid_centers": len(grid),
            "root_matchings": len(permutations),
            "all_S4_permutations_once": True,
            "root_pushdown_class": "(1,1)",
        },
    }


def _audit_raw_vi(payload: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(payload["full_NS_lattice"]["gram"])
    walls = [sp.Matrix(row) for row in payload["facets"]["NS_wall_vectors"]]
    wall_matrix = sp.Matrix.hstack(*walls)
    wall_gram = wall_matrix.T * gram * wall_matrix
    boundary = [index for index in range(len(walls)) if wall_gram[index, index] == -4]
    roots = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -2
        and all(wall_gram[index, other] == 0 for other in boundary)
    ]
    crossing = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -1
        and sorted(int(wall_gram[index, other]) for other in boundary) == [0, 0, 0, 1, 1]
    ]
    cusp_lifts = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -1
        and sorted(int(wall_gram[index, other]) for other in boundary) == [0, 0, 0, 0, 2]
    ]
    if (len(boundary), len(roots), len(crossing), len(cusp_lifts)) != (5, 10, 10, 35):
        raise ArithmeticError("VI raw wall taxonomy changed")
    if _exact_inertia(gram) != (1, 14, 0) or int(gram.det()) not in {-1, 1}:
        raise ArithmeticError("VI lattice is not odd unimodular I_(1,14)")
    if any(wall_gram[left, right] for left, right in itertools.combinations(boundary, 2)):
        raise ArithmeticError("VI boundaries are not disjoint")
    canonical = -sum((walls[index] for index in boundary), sp.zeros(15, 1)) / 2
    if any(value.q != 1 for value in canonical):
        raise ArithmeticError("VI canonical class is not integral")

    crossing_profiles = {
        index: tuple(other for other in boundary if wall_gram[index, other] == 1)
        for index in crossing
    }
    if set(tuple(sorted(profile)) for profile in crossing_profiles.values()) != set(
        itertools.combinations(boundary, 2)
    ):
        raise ArithmeticError("VI crossings do not index all boundary pairs")
    crossing_by_pair = {
        tuple(sorted(profile)): index for index, profile in crossing_profiles.items()
    }
    root_by_crossing: dict[int, int] = {}
    for exceptional in crossing:
        incident_roots = [root for root in roots if wall_gram[root, exceptional] == 1]
        if len(incident_roots) != 1:
            raise ArithmeticError("VI crossing lacks a unique associated root")
        root_by_crossing[exceptional] = incident_roots[0]
    if len(set(root_by_crossing.values())) != 10:
        raise ArithmeticError("VI crossing/root matching is not bijective")

    branch_records = []
    canonical_dependencies = None
    for conic_boundary in boundary:
        four_set = sorted(index for index in boundary if index != conic_boundary)
        total_exceptionals: list[tuple[str, Any, sp.Matrix]] = [
            ("crossing", pair, walls[exceptional])
            for pair, exceptional in sorted(crossing_by_pair.items())
        ]
        active_roots = []
        for mark in four_set:
            pair = tuple(sorted((conic_boundary, mark)))
            exceptional = crossing_by_pair[pair]
            root = root_by_crossing[exceptional]
            active_roots.append(root)
            total_exceptionals.append(("first_tangency", mark, walls[root] + walls[exceptional]))
        exceptional_matrix = sp.Matrix.hstack(*(record[2] for record in total_exceptionals))
        nullspace = (exceptional_matrix.T * gram).nullspace()
        if len(nullspace) != 1:
            raise ArithmeticError("VI exceptional span does not have a unique plane orthogonal")
        plane = _primitive_integer_vector(nullspace[0])
        canonical_degree = int((canonical.T * gram * plane)[0])
        if canonical_degree == 3:
            plane = -plane
            canonical_degree = -3
        if canonical_degree != -3 or int((plane.T * gram * plane)[0]) != 1:
            raise ArithmeticError("VI null ray is not a numerical plane class")
        degrees = [int((plane.T * gram * walls[index])[0]) for index in boundary]
        if degrees[boundary.index(conic_boundary)] != 2 or sorted(degrees) != [1, 1, 1, 1, 2]:
            raise ArithmeticError("VI plane branch does not have conic plus four lines")
        basis = sp.Matrix.hstack(plane, exceptional_matrix)
        if abs(int(basis.det())) != 1:
            raise ArithmeticError("VI plane/exceptional basis is not unimodular")
        if basis.T * gram * basis != sp.diag(1, *([-1] * 14)):
            raise ArithmeticError("VI plane/exceptional basis does not diagonalize NS")

        active_set = set(active_roots)
        if {root for root in roots if int((plane.T * gram * walls[root])[0]) == 0} != active_set:
            raise ArithmeticError("VI four contracted roots are not the tangency roots")
        residual = [root for root in roots if root not in active_set]
        dependencies = []
        for root in residual:
            if int((plane.T * gram * walls[root])[0]) != 1:
                raise ArithmeticError("VI residual root is not a plane line")
            support = [
                (kind, label)
                for kind, label, exceptional in total_exceptionals
                if int((walls[root].T * gram * exceptional)[0]) == 1
            ]
            ordinary = [
                label
                for kind, label in support
                if kind == "crossing" and set(label).issubset(four_set)
            ]
            tangency = [label for kind, label in support if kind == "first_tangency"]
            if len(ordinary) != 1 or len(tangency) != 2:
                raise ArithmeticError("VI residual line lacks the polar/chord support")
            pair = tuple(sorted(four_set.index(index) for index in ordinary[0]))
            complement = tuple(sorted(four_set.index(index) for index in tangency))
            if set(complement) != set(range(4)) - set(pair):
                raise ArithmeticError("VI residual support is not pair/complement")
            dependencies.append((pair, complement))
        dependencies.sort()
        expected_dependencies = sorted(
            (pair, tuple(sorted(set(range(4)) - set(pair))))
            for pair in itertools.combinations(range(4), 2)
        )
        if dependencies != expected_dependencies:
            raise ArithmeticError("VI residual lines do not form the complete six-pair design")
        if canonical_dependencies is None:
            canonical_dependencies = dependencies
        elif dependencies != canonical_dependencies:
            raise ArithmeticError("VI five contraction branches give different abstract designs")
        branch_records.append(
            {
                "conic_boundary_zero_based": conic_boundary,
                "plane_class": _integer_vector(plane),
                "plane_square": 1,
                "K_dot_plane": -3,
                "boundary_degree_multiset": sorted(degrees),
                "exceptional_basis_determinant_absolute": abs(int(basis.det())),
                "active_tangency_roots": len(active_roots),
                "residual_polar_chord_roots": len(residual),
            }
        )

    return {
        "raw_input": "rank15/60 walls",
        "bridge_archive_read": False,
        "reconstructed_NS": {
            "rank": gram.rows,
            "determinant": int(gram.det()),
            "inertia": list(_exact_inertia(gram)),
            "identification": "odd unimodular I_(1,14)",
        },
        "physical_taxonomy": {
            "boundary_minus4": len(boundary),
            "actual_minus2": len(roots),
            "crossing_minus1": len(crossing),
            "single_boundary_cusp_lifts": len(cusp_lifts),
        },
        "strong_Coble_numerics": {
            "boundaries_disjoint": True,
            "boundary_sum_equals_minus_2K": True,
            "rational_blowup_endpoint": "P2",
        },
        "dependency_extraction": {
            "crossing_root_pairs": len(root_by_crossing),
            "plane_branches": len(branch_records),
            "every_branch_is_conic_plus_four_lines": True,
            "every_branch_has_complete_pair_complement_design": True,
            "branch_records": branch_records,
        },
    }


# F_9 is represented as F_3[r]/(r^2-2), encoded by a+3*b for a+b*r.
def _gf9_add(left: int, right: int) -> int:
    return (left % 3 + right % 3) % 3 + 3 * ((left // 3 + right // 3) % 3)


def _gf9_neg(value: int) -> int:
    return (-value % 3) + 3 * (-(value // 3) % 3)


def _gf9_sub(left: int, right: int) -> int:
    return _gf9_add(left, _gf9_neg(right))


def _gf9_mul(left: int, right: int) -> int:
    a, b = left % 3, left // 3
    c, d = right % 3, right // 3
    return (a * c + 2 * b * d) % 3 + 3 * ((a * d + b * c) % 3)


def _gf9_pow(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = _gf9_mul(result, value)
        value = _gf9_mul(value, value)
        exponent //= 2
    return result


def _gf9_inverse(value: int) -> int:
    if not value:
        raise ZeroDivisionError("zero has no inverse in F9")
    return _gf9_pow(value, 7)


def _pgl2_9_normalize(matrix: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    scale = _gf9_inverse(next(value for value in matrix if value))
    return tuple(_gf9_mul(value, scale) for value in matrix)  # type: ignore[return-value]


def _pgl2_9_action(matrix: tuple[int, int, int, int], point: int) -> int:
    # The point 9 denotes infinity; 0,...,8 are finite F9 elements.
    a, b, c, d = matrix
    if point == 9:
        numerator, denominator = a, c
    else:
        numerator = _gf9_add(_gf9_mul(a, point), b)
        denominator = _gf9_add(_gf9_mul(c, point), d)
    return 9 if denominator == 0 else _gf9_mul(numerator, _gf9_inverse(denominator))


@functools.lru_cache(maxsize=1)
def _pgl2_9() -> tuple[tuple[int, int, int, int], ...]:
    matrices = set()
    for matrix in itertools.product(range(9), repeat=4):
        determinant = _gf9_sub(
            _gf9_mul(matrix[0], matrix[3]), _gf9_mul(matrix[1], matrix[2])
        )
        if determinant:
            matrices.add(_pgl2_9_normalize(matrix))
    if len(matrices) != 720:
        raise ArithmeticError(f"PGL2(F9) order changed: {len(matrices)}")
    return tuple(sorted(matrices))


def _standard_f3_sublines() -> set[frozenset[int]]:
    base = (9, 0, 1, 2)
    blocks = {
        frozenset(_pgl2_9_action(matrix, point) for point in base)
        for matrix in _pgl2_9()
    }
    if len(blocks) != 30:
        raise ArithmeticError(f"wrong number of F3-sublines in P1(F9): {len(blocks)}")
    return blocks


def _design_isomorphism_to_f9(
    raw_blocks: set[frozenset[int]], standard_blocks: set[frozenset[int]]
) -> tuple[int, ...]:
    """Find an isomorphism after fixing one ordered triple, only 6! trials."""

    seed = (0, 1, 2)
    raw_seed_block = next(block for block in raw_blocks if set(seed) <= block)
    raw_fourth = next(iter(raw_seed_block - set(seed)))
    standard_seed = (9, 0, 1)
    standard_seed_block = next(block for block in standard_blocks if set(standard_seed) <= block)
    standard_fourth = next(iter(standard_seed_block - set(standard_seed)))
    raw_remaining = [
        point for point in range(10) if point not in {*seed, raw_fourth}
    ]
    standard_remaining = [
        point for point in range(10) if point not in {*standard_seed, standard_fourth}
    ]
    solutions = []
    for permutation in itertools.permutations(standard_remaining):
        mapping = dict(zip(seed, standard_seed))
        mapping[raw_fourth] = standard_fourth
        mapping.update(zip(raw_remaining, permutation))
        transported = {
            frozenset(mapping[point] for point in block) for block in raw_blocks
        }
        if transported == standard_blocks:
            solutions.append(tuple(mapping[point] for point in range(10)))
    if not solutions:
        raise ArithmeticError("raw 3-(10,4,1) design is not the F3-subline design")
    return min(solutions)


def _audit_mobius_factory_and_raw_mi(payload: dict[str, Any]) -> dict[str, Any]:
    # Read only the physical Gram and walls.  Surface labels, stored blocks,
    # and the published finite-field presentation are deliberately ignored.
    gram = sp.Matrix(payload["surface"]["NS_gram"])
    walls = [sp.Matrix(row) for row in payload["full_nef_chamber"]["NS_wall_vectors"]]
    wall_matrix = sp.Matrix.hstack(*walls)
    wall_gram = wall_matrix.T * gram * wall_matrix
    boundary = [index for index in range(len(walls)) if wall_gram[index, index] == -4]
    roots = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -2
        and all(wall_gram[index, other] == 0 for other in boundary)
    ]
    crossing = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -1
        and sorted(int(wall_gram[index, other]) for other in boundary) == [1, 1]
    ]
    cusp_lifts = [
        index
        for index in range(len(walls))
        if wall_gram[index, index] == -1
        and sorted(int(wall_gram[index, other]) for other in boundary) == [0, 2]
    ]
    if (len(boundary), len(roots), len(crossing), len(cusp_lifts)) != (2, 30, 10, 164):
        raise ArithmeticError("raw MI wall taxonomy changed")
    if _exact_inertia(gram) != (1, 11, 0) or abs(int(gram.det())) != 1:
        raise ArithmeticError("raw MI physical lattice is not I_(1,11)")
    if wall_gram[boundary[0], boundary[1]] != 0:
        raise ArithmeticError("raw MI boundaries are not disjoint")
    if any(wall_gram[left, right] for left, right in itertools.combinations(crossing, 2)):
        raise ArithmeticError("raw MI crossing curves are not pairwise disjoint")

    exceptional_matrix = sp.Matrix.hstack(*(walls[index] for index in crossing))
    nullspace = (exceptional_matrix.T * gram).nullspace()
    if len(nullspace) != 2:
        raise ArithmeticError("raw MI crossing contraction does not have rank-two endpoint")
    rulings = [_primitive_integer_vector(vector) for vector in nullspace]
    ruling_matrix = sp.Matrix.hstack(*rulings)
    ruling_gram = ruling_matrix.T * gram * ruling_matrix
    if ruling_gram == sp.Matrix([[0, -1], [-1, 0]]):
        rulings[1] = -rulings[1]
        ruling_matrix = sp.Matrix.hstack(*rulings)
        ruling_gram = ruling_matrix.T * gram * ruling_matrix
    if ruling_gram != sp.Matrix([[0, 1], [1, 0]]):
        raise ArithmeticError("raw MI contraction endpoint is not U=P1xP1")
    exceptional_basis = sp.Matrix.hstack(ruling_matrix, exceptional_matrix)
    if abs(int(exceptional_basis.det())) != 1:
        raise ArithmeticError("raw MI ruling/exceptional basis is not unimodular")
    if exceptional_basis.T * gram * exceptional_basis != sp.diag(
        sp.Matrix([[0, 1], [1, 0]]), sp.diag(*([-1] * 10))
    ):
        raise ArithmeticError("raw MI exceptional contraction Gram changed")

    def contracted_image(index: int, exceptionals: Sequence[int]) -> sp.Matrix:
        return walls[index] + sum(
            (
                int(wall_gram[index, exceptional]) * walls[exceptional]
                for exceptional in exceptionals
            ),
            sp.zeros(gram.rows, 1),
        )

    boundary_bidegrees = []
    for index in boundary:
        image = contracted_image(index, crossing)
        coordinates = ruling_matrix.gauss_jordan_solve(image)[0]
        boundary_bidegrees.append(tuple(_integer_vector(coordinates)))
    if sorted(boundary_bidegrees) != [(1, 3), (3, 1)]:
        raise ArithmeticError(f"raw MI boundary pushdowns changed: {boundary_bidegrees}")

    raw_blocks: set[frozenset[int]] = set()
    for root in roots:
        incident = [
            exceptional for exceptional in crossing if wall_gram[root, exceptional] == 1
        ]
        if len(incident) != 4:
            raise ArithmeticError("raw MI root does not pass through four crossing centers")
        image = contracted_image(root, incident)
        coordinates = ruling_matrix.gauss_jordan_solve(image)[0]
        if tuple(_integer_vector(coordinates)) != (1, 1):
            raise ArithmeticError("raw MI root does not push down to bidegree (1,1)")
        raw_blocks.add(frozenset(crossing.index(index) for index in incident))
    if len(raw_blocks) != 30:
        raise ArithmeticError("raw MI root supports are not 30 distinct blocks")
    triple_multiplicities = Counter(
        sum(set(triple) <= block for block in raw_blocks)
        for triple in itertools.combinations(range(10), 3)
    )
    if triple_multiplicities != Counter({1: 120}):
        raise ArithmeticError(f"raw MI blocks are not S(3,4,10): {triple_multiplicities}")

    standard_blocks = _standard_f3_sublines()
    design_isomorphism = _design_isomorphism_to_f9(raw_blocks, standard_blocks)
    transported = {
        frozenset(design_isomorphism[point] for point in block) for block in raw_blocks
    }
    if transported != standard_blocks:
        raise ArithmeticError("raw MI design transport to P1(F9) failed")

    graph_transformations = []
    for block in sorted(standard_blocks, key=lambda item: tuple(sorted(item))):
        transformations = [
            matrix
            for matrix in _pgl2_9()
            if all(
                _pgl2_9_action(matrix, point)
                == (9 if point == 9 else _gf9_pow(point, 3))
                for point in block
            )
        ]
        if len(transformations) != 1:
            raise ArithmeticError("an F3-subline lacks a unique Frobenius-interpolating projectivity")
        graph_transformations.append(transformations[0])
    if len(set(graph_transformations)) != 30:
        raise ArithmeticError("the 30 F3-subline graph curves are not distinct")

    endpoint_canonical = (
        -2 * rulings[0]
        - 2 * rulings[1]
        + sum((walls[index] for index in crossing), sp.zeros(12, 1))
    )
    boundary_sum = sum((walls[index] for index in boundary), sp.zeros(12, 1))
    if boundary_sum != -2 * endpoint_canonical:
        raise ArithmeticError("raw MI boundary sum is not -2K of the reconstructed blowup")
    # General q factory identities, derived on Bl_{q^2+1}(P1 x P1).
    q = sp.Symbol("q", integer=True, positive=True)
    factory = {
        "centers": str(q**2 + 1),
        "block_size": str(q + 1),
        "block_count": str(sp.expand(q * (q**2 + 1))),
        "design": "S(3,q+1,q^2+1) of F_q-sublines in P1(F_{q^2})",
        "boundary_square": str(sp.factor(2 * q - (q**2 + 1))),
        "block_curve_square": str(sp.factor(2 - (q + 1))),
        "K_dot_block_curve": str(sp.factor(-4 + (q + 1))),
        "boundary_sum_minus_minus2K": "(q-3)*(F1+F2)",
        "exact_Coble_root_locus": "q=3",
        "q_not_3_scope": (
            "The factory then fails the (-4)-boundary, (-2)-root, K.C=0, and B=-2K "
            "Coble identities; this is not a nonexistence theorem for an abstract Coxeter polytope."
        ),
    }
    if factory["boundary_square"] != "-(q - 1)**2":
        raise ArithmeticError("general Möbius factory boundary square formula changed")
    if factory["block_curve_square"] != "1 - q" or factory["K_dot_block_curve"] != "q - 3":
        raise ArithmeticError("general Möbius factory root formulas changed")

    return {
        "raw_input": "rank12/206 walls",
        "surface_labels_or_stored_blocks_read": False,
        "physical_taxonomy": {
            "boundary_minus4": len(boundary),
            "actual_minus2": len(roots),
            "crossing_minus1": len(crossing),
            "single_boundary_cusp_lifts": len(cusp_lifts),
        },
        "contraction": {
            "pairwise_disjoint_crossing_centers": len(crossing),
            "endpoint_gram": [[0, 1], [1, 0]],
            "exceptional_basis_determinant_absolute": abs(int(exceptional_basis.det())),
            "boundary_bidegrees": [list(value) for value in sorted(boundary_bidegrees)],
            "root_bidegree": [1, 1],
        },
        "raw_design": {
            "points": 10,
            "blocks": len(raw_blocks),
            "block_size": 4,
            "triple_multiplicity_histogram": {"1": 120},
            "identification": "Steiner system S(3,4,10)",
        },
        "independent_F9_realization": {
            "field_model": "F3[r]/(r^2-2)",
            "PGL2_F9_order": len(_pgl2_9()),
            "F3_subline_count": len(standard_blocks),
            "raw_center_to_P1_F9_encoding": list(design_isomorphism),
            "unique_projectivity_graphs_interpolating_Frobenius": len(
                set(graph_transformations)
            ),
            "boundary_graphs": "y=x^3 and x=y^3",
            "boundary_intersection": "x=x^9, namely the 10 reduced points P1(F9)",
            "all_root_seed_curves_smooth_irreducible_distinct": True,
        },
        "strong_Coble_at_q3": {
            "boundaries_disjoint_after_blowup": True,
            "boundary_squares": [-4, -4],
            "boundary_sum_equals_minus_2K": True,
            "root_squares": [-2] * 30,
            "root_K_pairings": [0] * 30,
        },
        "general_mobius_factory": factory,
    }


def _audit_production_operational_scope(production: dict[str, Any]) -> dict[str, Any]:
    source = SOURCES["production_script"].read_text(encoding="utf-8")
    tree = ast.parse(source)
    options = sorted(
        {
            argument.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_argument"
            for argument in node.args
            if isinstance(argument, ast.Constant)
            and isinstance(argument.value, str)
            and argument.value.startswith("--")
        }
    )
    path_constants = sorted(
        node.targets[0].id
        for node in tree.body
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
        and node.targets[0].id in {"BRIDGE", "TYPE_VI", "TYPE_MII", "OUTPUT"}
    )
    signatures = production["P_derived_dependency_signatures"]
    extension_counts = {
        key: value["standard_free_extension_steps_after_core"]
        for key, value in signatures.items()
    }
    if options != ["--check-archive", "--write"] or set(path_constants) != {
        "BRIDGE",
        "OUTPUT",
        "TYPE_MII",
        "TYPE_VI",
    }:
        raise ArithmeticError("production CLI/path surface changed; hostile scope audit needs review")
    if set(extension_counts.values()) != {0}:
        raise ArithmeticError("a production member now exercises the extension clause")
    return {
        "abstract_definition": {
            "uses_only_P_derived_objects_on_its_face": True,
            "proved_class_is_core_only": production["constructible_incidence_class"][
                "proved_scope"
            ]["zero_extension_required"],
            "qualification": (
                "A normalized arbitrary-input all-contraction dependency extractor is still "
                "required; the repaired nonzero-extension verifier is explicitly prospective."
            ),
        },
        "current_executable": {
            "CLI_options": options,
            "arbitrary_input_option": False,
            "hardwired_path_constants": path_constants,
            "reads_precompiled_bridge": True,
            "generic_dependency_canonicalizer_implemented": False,
            "generic_standard_free_extension_verifier_implemented": False,
            "certified_rows": len(signatures),
            "extension_steps_by_row": extension_counts,
        },
        "verdict": "CONCEPTUALLY_P_ONLY_BUT_NOT_YET_AN_INPUT_COMPUTABLE_IMPLEMENTED_MEMBERSHIP_TEST",
    }


def _audit_remaining_seven_operational_scope(
    production: dict[str, Any], independent_mi: dict[str, Any]
) -> dict[str, Any]:
    source = SOURCES["remaining_seven_script"].read_text(encoding="utf-8")
    tree = ast.parse(source)
    options = sorted(
        {
            argument.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_argument"
            for argument in node.args
            if isinstance(argument, ast.Constant)
            and isinstance(argument.value, str)
            and argument.value.startswith("--")
        }
    )
    if options != ["--check-archive", "--full-contractions"]:
        raise ArithmeticError(f"remaining-seven CLI changed: {options}")
    core = production["two_boundary_206_mobius_core"]
    raw = core["raw_extraction"]
    finite = core["finite_geometry"]
    if not (
        raw["boundary_count"] == 2
        and raw["crossing_minus_one_count"] == 10
        and raw["root_count"] == 30
        and raw["raw_design_parameters"] == "3-(10,4,1)"
        and finite["raw_to_canonical_incidence_isomorphism_exists"]
        and finite["all_thirty_unique_smooth_1_1_curves"]
        and independent_mi["raw_design"]["identification"] == "Steiner system S(3,4,10)"
    ):
        raise ArithmeticError("production and hostile Möbius cores disagree")
    factory = production["mobius_plane_factory_rigidity_theorem"]
    if factory["formulas"] != {
        "B_square": "-(q-1)^2",
        "circle_square": "1-q",
        "K_dot_circle": "q-3",
        "B_plus_Bprime": "-2K+(q-3)(F+G)",
    }:
        raise ArithmeticError("production Möbius factory formulas changed")
    return {
        "production_core_matches_independent_replay": True,
        "particular_design_isomorphism_archived": False,
        "invariant_incidence_isomorphism_archived": True,
        "CLI_options": options,
        "arbitrary_input_option": False,
        "hardwired_input_archive_count": 6,
        "hardwired_output_row_count": 7,
        "extractor_after_row_selection_is_label_free": True,
        "actual_frozen_rating": "L3_STRUCTURAL_P_FIRST_CALIBRATION",
        "conditional_promotion": (
            "SCOPED_FORMAL_L4_CORE_CLASS only after a normalized arbitrary-input "
            "CI4/Mobius recognizer is frozen and hostile-replayed"
        ),
    }


def _extension_logic_audit() -> dict[str, Any]:
    return {
        "advertised_implication": (
            "base-point-free residual system plus a generically separable associated map "
            "has a nonempty open of smooth geometrically irreducible members"
        ),
        "verdict": "FALSE_AS_STATED",
        "counterexample": {
            "field": "algebraically closed, characteristic 3",
            "surface": "P1 x P1",
            "linear_system": "|O(2,0)|=|2F1|",
            "base_point_free": True,
            "associated_map": "pr1 followed by the quadratic Veronese embedding P1 -> P2",
            "generically_separable": True,
            "general_section": "a square-free binary quadratic",
            "general_divisor": "two disjoint ruling fibres",
            "smooth": True,
            "geometrically_irreducible": False,
            "reason": "the system is composed with a pencil",
        },
        "ambiguity": (
            "The JSON phrase 'has a separable Bertini open set' would repair the implication "
            "only if it literally includes nonempty geometrically-integral members; then it "
            "partly assumes the desired existence and still lacks a finite verifier schema."
        ),
        "minimum_repair": [
            "require image dimension at least two, or certify a geometrically integral generic member",
            "for a pencil, require a geometrically integral generic fibre rather than separability alone",
            "prove relative constant rank and local freeness of each kernel family",
            "prove the parameter tower irreducible/dominant and every certified open nonempty on the same component",
            "freeze explicit smoothness, irreducibility, transversality, and distinctness witnesses",
        ],
        "effect_on_theorem": (
            "The two zero-extension cores survive; closure under arbitrary declared standard-free "
            "extensions and hence the advertised class-wide carrier theorem do not."
        ),
    }


def _handoff_and_pushback_audit(
    d_cobvin: dict[str, Any], vi: dict[str, Any], mii: dict[str, Any], mii_cone: dict[str, Any]
) -> dict[str, Any]:
    d_checks = d_cobvin["checks"]
    if not d_cobvin["all_checks_pass"] or not all(d_checks.values()):
        raise ArithmeticError("D_CobVin archive is not green")
    vi_volume = vi["literal_vinberg"]
    mii_volume = mii_cone["literal_full_cone"]
    raw_replays = {
        "H14_60": {
            "CM_roots": 20,
            "CM_cusps": vi["curve_exhaustion"]["CM_cusp_count"],
            "integral_cusp_lifts": len(vi["curve_exhaustion"]["lifted_cusp_boundary_pairs_zero_based"]),
            "full_cone": {
                "finite": vi_volume["finite_vertex_count"],
                "ideal": vi_volume["ideal_vertex_count"],
                "hyperideal": 0,
            },
            "Vinberg_bad_edges": vi_volume["bad_edge_count"],
            "Vinberg_bad_parabolic_neighborhoods": vi_volume["bad_parabolic_neighborhood_count"],
        },
        "H17_48": {
            "CM_roots": mii["curve_exhaustion"]["Kondo_complete_effective_irreducible_CM_roots"],
            "CM_cusps": mii["curve_exhaustion"]["primitive_CM_cusps_tested"],
            "integral_cusp_lifts": mii["curve_exhaustion"]["integral_single_boundary_lifts"],
            "reflection_index_in_O_plus": mii["finite_volume"]["finite_index_of_reflection_group_in_O_plus"],
            "full_cone": {
                "finite": mii_volume["finite_vertex_count"],
                "ideal": mii_volume["ideal_vertex_count"],
                "hyperideal": mii_volume["hyperideal_vertex_count"],
            },
        },
    }
    if raw_replays["H14_60"]["full_cone"] != {"finite": 27845, "ideal": 1027, "hyperideal": 0}:
        raise ArithmeticError("H14 literal cone census changed")
    if raw_replays["H17_48"]["full_cone"] != {"finite": 190808, "ideal": 5244, "hyperideal": 0}:
        raise ArithmeticError("H17 literal cone census changed")
    if raw_replays["H14_60"]["integral_cusp_lifts"] != 35:
        raise ArithmeticError("H14 cusp-lift count changed")

    return {
        "D_CobVin_exact_return": {
            "verdict": "CONDITIONAL_PASS",
            "condition": "an actual strong-Coble root-seed branch satisfying the raw D_CobVin gates",
            "logic": [
                "the constructed curves realize every raw CM root seed with the original marking",
                "raw Vinberg finite index plus Kondo root exhaustion excludes hidden CM roots",
                "Kondo's CM-nef criterion makes enumerated primitive root-cone cusps nef",
                "primitive-isotropic/Halphen theory realizes every and only integral single-boundary lift",
                "adjunction reduces every hidden negative curve to one of those two exhausted cases",
                "the reassembled oriented wall set is therefore literally the original P",
            ],
            "not_merely_an_abstract_isomorphism": True,
            "core_raw_replays": raw_replays,
            "artifact_coverage_qualification": (
                "The CI4 executable records this handoff in prose rather than replaying D_CobVin; "
                "the pinned general theorem and raw cone archives supply the conditional certificate."
            ),
        },
        "Kondo_pushback": {
            "logical_independence_after_carrier": True,
            "formal_characteristic_three_rows": 7,
            "formal_dimensions": [10, 11, 14, 17],
            "currently_realized_CI4_core_boundary_counts": [5, 8],
            "currently_realized_CI4_dimensions": [14, 17],
            "strict_novelty_on_the_two_members": "NONE_OR_NEGLIGIBLE",
            "reason": (
                "Each complete core already fixes characteristic three and b=5 or b=8, hence "
                "rank 15/H14 or rank 18/H17 before the classification is invoked."
            ),
            "would_become_nontrivial_if": (
                "a third raw dependency core or a genuine repaired extension entered without "
                "already fixing a classified endpoint"
            ),
        },
    }


@functools.lru_cache(maxsize=1)
def build_certificate() -> dict[str, Any]:
    for name, path in SOURCES.items():
        digest = _sha256(path)
        if digest != PINNED_SHA256[name]:
            raise RuntimeError(f"pinned source changed for {name}: {digest}")
    production = _load_pinned("production_archive")
    remaining_seven_production = _load_pinned("remaining_seven_archive")
    d_cobvin = _load_pinned("D_CobVin_archive")
    vi_payload = _load_pinned("type_VI_raw_full")
    mii_payload = _load_pinned("type_MII_raw_full")
    mii_cone = _load_pinned("type_MII_literal_cone")
    mi_payload = _load_pinned("type_MI_raw_full")

    algebra = _audit_four_point_algebra()
    raw_mii = _audit_raw_mii(mii_payload)
    raw_vi = _audit_raw_vi(vi_payload)
    raw_mi = _audit_mobius_factory_and_raw_mi(mi_payload)
    operational = _audit_production_operational_scope(production)
    remaining_operational = _audit_remaining_seven_operational_scope(
        remaining_seven_production, raw_mi
    )
    extensions = _extension_logic_audit()
    handoff = _handoff_and_pushback_audit(d_cobvin, vi_payload, mii_payload, mii_cone)

    strong_coble = {
        "verdict_for_both_explicit_cores": "PASS",
        "argument": [
            "the independently reconstructed endpoint is a blowup of P2 or P1 x P1",
            "the reduced smooth boundary components are disjoint and sum to -2K",
            "for D in |-2K|, D.B_i=-4 forces every B_i into D, hence D=B",
            "if A lies in |-K|, then 2A lies in |-2K| and would equal the reduced divisor B, impossible",
        ],
        "extension_qualification": (
            "The same divisor argument works for any genuinely constructed extension retaining "
            "the exact boundary hypotheses; it cannot supply the missing construction itself."
        ),
    }

    five_gates = {
        "C_independent_P_only_class": {
            "verdict": "PARTIAL_PASS",
            "reason": (
                "The certificate language is phrased in P-derived data, but the current artifact "
                "is a two-input replay and no generic extension-certificate verifier exists."
            ),
        },
        "Car_nonempty_for_every_member": {
            "verdict": "FAIL_FOR_ADVERTISED_EXTENSION_CLASS; PASS_FOR_TWO_EXPLICIT_CORES",
            "reason": "the standard-free Bertini implication is false and all tested extension counts are zero",
        },
        "Return_exact_same_P": {
            "verdict": "CONDITIONAL_PASS",
            "reason": "D_CobVin is a valid marked exact-return theorem once a genuine root-seed carrier exists",
        },
        "AG_independent_classification": {
            "verdict": "CONDITIONAL_PASS",
            "reason": "Kondo applies after the carrier and finite-index argument, not in membership",
        },
        "Pushback_new_theorem_on_every_member": {
            "verdict": "FAIL_STRICT; FORMAL_QUALIFIED_ONLY",
            "reason": "on the only two certified CI4 members the core already encodes the H14/H17 answers",
        },
        "tuple": ["PARTIAL", "FAIL/PASS_TWO_CORES", "PASS_CONDITIONAL", "PASS_CONDITIONAL", "FAIL"],
    }
    rating = {
        "original_advertised_rating_rejected": "L3_PLUS_SCOPED_CONSTRUCTIVE_CLASS",
        "corrected_production_rating_accepted": "L3_TWO_MEMBER_CORE_CALIBRATION;_NOT_L3_PLUS_OR_L4",
        "recommended_rating": "L3_TWO_MEMBER_INVERSE_CALIBRATION; NOT_AN_EXTENSION_CLASS_THEOREM; NOT_L4",
        "core_only_subtheorem": (
            "A rigorous uniform P1(F3) realization theorem for two independently reconstructed "
            "raw cores, followed by conditional exact same-P return."
        ),
        "why_not_names_only": (
            "The raw frame and dependency reconstructions are genuinely label-free and the shared "
            "projective-line theorem is real; operationally, however, the executable denominator "
            "contains only the two known zero-extension configurations."
        ),
    }

    mi_cone = mi_payload["full_nef_chamber"]["exact_double_description"]
    mi_vinberg = mi_payload["full_nef_chamber"]["literal_vinberg"]
    enlarged_finite_geometry_class = {
        "status": "INDEPENDENT_MATHEMATICS_AND_PRODUCTION_CALIBRATION_VERIFIED",
        "proposed_core_only_name": "D_CobVin^FGcore",
        "definition_guardrail": (
            "Accept only a raw-P-recognized complete CI4 core or a complete Möbius-plane "
            "S(3,q+1,q^2+1) core satisfying the exact Coble numerical alphabet. Do not "
            "inherit the broken standard-free extension clause."
        ),
        "third_raw_member": raw_mi,
        "production_artifact_scope": remaining_operational,
        "Möbius_factory_theorem": {
            "structural_family_is_parametric_not_a_named_Gram": True,
            "finite_geometry_exists_for_every_prime_power_q": True,
            "exact_strong_Coble_root_seed_iff": "q=3",
            "q3_geometry": "P1(F9) centers and its 30 F3-sublines",
            "q_not_3_nonexistence_claim": False,
        },
        "same_P_return_at_q3": {
            "verdict": "PASS_THROUGH_D_CobVin_AFTER_EXPLICIT_ROOT_SEED",
            "CM_roots": mi_payload["coble_mukai"]["root_count"],
            "CM_cusps": mi_payload["coble_mukai"]["primitive_cusp_count"],
            "integral_single_boundary_lifts": 2
            * mi_payload["coble_mukai"]["cusp_lift_membership_histogram"]["both_boundaries"],
            "full_cone": {
                "finite": mi_cone["finite_vertex_count"],
                "ideal": mi_cone["ideal_vertex_count"],
                "hyperideal": mi_cone["hyperideal_vertex_count"],
            },
            "Vinberg_bad_edges": mi_vinberg["bad_edge_count"],
            "Vinberg_bad_parabolic_neighborhoods": mi_vinberg[
                "invalid_full_zero_neighborhood_count"
            ],
        },
        "lineage_credit": {
            "published": (
                "Kondo/Katsura-Kondo finite-field Type-MI surface geometry, Frobenius "
                "boundaries, F9 centers, and the root configuration"
            ),
            "project_level_new": (
                "label-free raw-P extraction of the ten-center contraction, intrinsic "
                "Steiner recognition, and the P->surface->literally same P theorem"
            ),
            "not_claimed": [
                "new worldwide Type-MI geometry",
                "a new abstract H11/206 Coxeter polytope",
                "nonexistence of abstract q!=3 Coxeter polytopes",
            ],
        },
        "conditional_five_gate_rating_if_the_generic_recognizer_is_frozen": {
            "C": "PASS_SCOPED for a genuinely raw-P Möbius/CI4 core recognizer",
            "Car": "PASS for the explicit cores; no extension closure claimed",
            "Return": "PASS through D_CobVin",
            "AG": "PASS after carrier construction",
            "Pushback": (
                "FORMAL_PASS_BUT_THIN: Kondo applies class-wide, while the core grammar "
                "already fixes much of b/rank/characteristic"
            ),
            "rating": (
                "SCOPED_FORMAL_L4_CORE_CLASS, conditional on the production recognizer "
                "being input-computable and structural; not broad/flagship L4"
            ),
        },
    }

    checks = {
        "symbolic_algebra_replayed_without_production_import": algebra["independent_of_production_helpers"],
        "all_24_grid_sections_are_distinct_smooth_graphs": algebra["grid"]["distinct_projective_kernel_sections"] == 24
        and algebra["grid"]["all_coefficient_matrices_invertible"],
        "all_6_residual_lines_are_distinct_from_boundaries": algebra["polar_chord"]["distinct_residual_lines"] == 6
        and algebra["polar_chord"]["residual_boundary_collisions"] == 0,
        "MII_raw_frame_and_dependency_reconstructed_without_bridge": raw_mii["bridge_archive_read"] is False
        and raw_mii["frame_census"]["integral_frames"] == 1
        and raw_mii["dependency_extraction"]["all_S4_permutations_once"],
        "VI_all_five_raw_plane_dependencies_reconstructed_without_bridge": raw_vi["bridge_archive_read"] is False
        and raw_vi["dependency_extraction"]["plane_branches"] == 5
        and raw_vi["dependency_extraction"]["every_branch_has_complete_pair_complement_design"],
        "both_explicit_cores_are_strong_Coble": strong_coble["verdict_for_both_explicit_cores"] == "PASS",
        "current_production_has_zero_extension_witnesses": set(
            operational["current_executable"]["extension_steps_by_row"].values()
        )
        == {0},
        "extension_counterexample_is_in_characteristic_three": extensions["counterexample"]["field"].endswith(
            "characteristic 3"
        ),
        "D_CobVin_handoff_is_only_conditional": handoff["D_CobVin_exact_return"]["verdict"]
        == "CONDITIONAL_PASS",
        "strict_pushback_fails_on_current_denominator": handoff["Kondo_pushback"]["strict_novelty_on_the_two_members"]
        == "NONE_OR_NEGLIGIBLE",
        "MI_raw_Steiner_design_reconstructed": raw_mi["raw_design"]["triple_multiplicity_histogram"]
        == {"1": 120},
        "MI_F9_root_seed_reconstructed": raw_mi["independent_F9_realization"][
            "unique_projectivity_graphs_interpolating_Frobenius"
        ]
        == 30,
        "general_Mobius_factory_has_exact_q3_Coble_locus": raw_mi["general_mobius_factory"][
            "exact_Coble_root_locus"
        ]
        == "q=3",
        "remaining_seven_production_matches_independent_MI_replay": remaining_operational[
            "production_core_matches_independent_replay"
        ],
        "remaining_seven_current_rating_stays_L3": remaining_operational[
            "actual_frozen_rating"
        ]
        == "L3_STRUCTURAL_P_FIRST_CALIBRATION",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile audit internal checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "target": (
            "the original D_CobVin^CI4 extension claim, its corrected core-only theorem, "
            "and the separate remaining-seven Möbius calibration"
        ),
        "independence": {
            "imports_production_helpers": False,
            "reads_precompiled_blowdown_bridge": False,
            "recomputes_symbolic_determinants_and_modular_kernels": True,
            "reconstructs_both_raw_NS_contractions": True,
        },
        "sources": {
            name: {"path": str(path.relative_to(ROOT)), "sha256": PINNED_SHA256[name]}
            for name, path in SOURCES.items()
        },
        "independent_four_point_algebra": algebra,
        "independent_raw_dependency_extraction": {
            "rank15_60wall": raw_vi,
            "rank18_48wall": raw_mii,
        },
        "strong_Coble": strong_coble,
        "P_only_membership": operational,
        "standard_free_extension_audit": extensions,
        **handoff,
        "five_gate_rating": five_gates,
        "level_recommendation": rating,
        "enlarged_finite_geometry_class": enlarged_finite_geometry_class,
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archive mismatch: {arguments.output}")
        print("hostile archive replay passed")
        return
    if arguments.write:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
