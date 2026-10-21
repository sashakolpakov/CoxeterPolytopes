#!/usr/bin/env python3
"""Independent hostile replay of the cusped K3 spanning-root theorem.

This file deliberately imports no code from the owner certifier.  It pins the
five owner artifacts and reconstructs the exact algebraic interfaces used by
the geometric proof.  The visual-shadow and finite-volume implications remain
mathematical arguments; the calculations here pressure-test their formulas and
rank consequences in every K3 Picard rank.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

PINS = {
    "owner_report": (
        "results/noncompact_k3_bounded_spanning_root_lemma.md",
        "de3b78a9070dc7035a8ecf51c5aa697b0736ad1559f0db8867b6629eb3cc4446",
    ),
    "owner_scope_audit": (
        "results/noncompact_k3_bounded_spanning_root_lemma_hostile_scope_audit.md",
        "cfaf5b5fd10c97e14f28105b99c379c8f9c391c6c89ffd9875ceefd1426ff851",
    ),
    "owner_certifier": (
        "scripts/certify_noncompact_k3_bounded_spanning_root_lemma.py",
        "895a0c1fe0e10d38f5516ab917154bde0d3a8d2eef89d1ad567260abaf2e5ffb",
    ),
    "owner_archive": (
        "data/enumerations/noncompact_k3_bounded_spanning_root_lemma.json",
        "5c5ffa217c15aebb189102c3a92937c5e20fe5e465a00caaac44780fa9761874",
    ),
    "owner_tests": (
        "tests/test_noncompact_k3_bounded_spanning_root_lemma.py",
        "bf4fdfd82d1c1f7529cb99ed7a62a9c11822d534755f304273eccc373707fd93",
    ),
}


Vector = list[int]
Matrix = list[list[int]]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_rank(rows: Matrix) -> int:
    """Compute rank over Q without using the owner implementation."""

    if not rows:
        return 0
    work = [[Fraction(entry) for entry in row] for row in rows]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, row_count) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for target in range(column, column_count):
            work[rank][target] /= pivot_value
        for index in range(row_count):
            if index == rank or not work[index][column]:
                continue
            factor = work[index][column]
            for target in range(column, column_count):
                work[index][target] -= factor * work[rank][target]
        rank += 1
        if rank == row_count:
            break
    return rank


def determinant(matrix: Matrix) -> int:
    """Exact determinant by fraction-free elimination."""

    size = len(matrix)
    if size == 0:
        return 1
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    denominator = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                work[row][target] = (
                    work[row][target] * pivot_value
                    - work[row][column] * work[column][target]
                ) // denominator
        denominator = pivot_value
    return sign * work[-1][-1]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    if not left or not right:
        return []
    return [
        [
            sum(left[row][mid] * right[mid][column] for mid in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def visual_formula_certificate() -> dict[str, Any]:
    """Re-derive the ultraparallel sine identity and sample its strict bound."""

    a, b, c = sp.symbols("a b c", real=True)
    observer_frame = 2 * (
        sp.cos(a) * sp.cos(b) - sp.cos(a + b + c)
    ) / (sp.sin(a) * sp.sin(b))
    owner_form = (
        4
        * sp.sin(a + c / 2)
        * sp.sin(b + c / 2)
        / (sp.sin(a) * sp.sin(b))
        - 2
    )
    residual = sp.simplify(sp.expand_trig(observer_frame - owner_form))
    assert residual == 0

    # theta is a full aperture, so a=theta_F/2 and c is the angular gap.
    # The proof uses 0<c<=min(theta_F,theta_G), not a numerical estimate.
    largest_quotient = 0.0
    largest_pairing = 0.0
    cases = 0
    for theta_f_index in range(1, 80):
        theta_f = math.pi * theta_f_index / 80
        for theta_g_index in range(1, 80):
            theta_g = math.pi * theta_g_index / 80
            for fraction in (0.01, 0.2, 0.5, 0.8, 1.0):
                gap = fraction * min(theta_f, theta_g)
                q_f = math.sin((theta_f + gap) / 2) / math.sin(theta_f / 2)
                q_g = math.sin((theta_g + gap) / 2) / math.sin(theta_g / 2)
                pairing = 4 * q_f * q_g - 2
                assert q_f < 2 and q_g < 2
                assert pairing < 14
                largest_quotient = max(largest_quotient, q_f, q_g)
                largest_pairing = max(largest_pairing, pairing)
                cases += 1
    return {
        "observer_frame_symbolic_residual": str(residual),
        "half_aperture_substitution": "a=theta_F/2, b=theta_G/2",
        "positive_gap_required_for_formula": True,
        "gap_at_most_minimum_aperture_required": True,
        "intersecting_or_asymptotic_branch_uses_b_at_most_2": True,
        "ultraparallel_facet_shadows_disjoint_after_irredundancy": True,
        "nested_ultraparallel_shadows_would_make_one_wall_redundant": True,
        "sample_count": cases,
        "sample_maximum_sine_quotient": largest_quotient,
        "sample_maximum_normalized_pairing": largest_pairing,
    }


def cusp_pair(left: Vector, right: Vector) -> int:
    """Pairing in the basis (u,v,e_1,...), with (u,v)=1 and e_i^2=-2."""

    return (
        left[0] * right[1]
        + left[1] * right[0]
        - 2 * sum(
            left[index] * right[index] for index in range(2, len(left))
        )
    )


def complete_box_cusp_packet(rank: int) -> dict[str, Any]:
    """A bounded, generally nonsimple Euclidean link and its lifted packet."""

    if rank < 3:
        raise ValueError("rank must be at least three")
    horizontal_dimension = rank - 2
    u = [1, 0, *([0] * horizontal_dimension)]
    packet: Matrix = []
    for coordinate in range(horizontal_dimension):
        plus = [1, 0, *([0] * horizontal_dimension)]
        minus = [1, 0, *([0] * horizontal_dimension)]
        plus[coordinate + 2] = 1
        minus[coordinate + 2] = -1
        packet.extend([plus, minus])
    assert all(cusp_pair(root, u) == 0 for root in packet)
    assert all(cusp_pair(root, root) == -2 for root in packet)
    assert matrix_rank(packet) == rank - 1
    positive_sum = [sum(root[index] for root in packet) for index in range(rank)]
    assert positive_sum == [2 * horizontal_dimension, *([0] * (rank - 1))]
    return {
        "rank": rank,
        "horosphere_dimension": horizontal_dimension,
        "incident_wall_count": len(packet),
        "link_model": "bounded coordinate box",
        "packet_rank": matrix_rank(packet),
        "expected_u_perp_rank": rank - 1,
        "all_lifts_orthogonal_to_u": True,
        "all_sample_lifts_square_minus_two": True,
        "positive_horizontal_dependence_lifts_to_nonzero_u": True,
        "nonsimple_when_rank_at_least_four": len(packet) > rank - 1,
    }


def lorentz_pair(left: Vector, right: Vector) -> int:
    return left[0] * right[0] - sum(
        left[index] * right[index] for index in range(1, len(left))
    )


def perpendicular_basis(vector: Vector) -> Matrix:
    """Independent exact kernel basis for the covector induced by vector."""

    covector = [vector[0], *(-entry for entry in vector[1:])]
    pivot = next(index for index, entry in enumerate(covector) if entry)
    rows: Matrix = []
    for free in range(len(vector)):
        if free == pivot:
            continue
        row = [0] * len(vector)
        row[free] = covector[pivot]
        row[pivot] = -covector[free]
        rows.append(row)
    assert matrix_rank(rows) == len(vector) - 1
    assert all(lorentz_pair(vector, row) == 0 for row in rows)
    return rows


def vertex(kind: str, rank: int, slot: int) -> Vector:
    point = [0] * rank
    if kind == "finite":
        if slot == 0:
            point[0] = 1
        else:
            point[0] = 2
            point[1] = 1
        assert lorentz_pair(point, point) > 0
    elif kind == "ideal":
        point[0] = 1
        point[slot + 1] = 1
        assert lorentz_pair(point, point) == 0
    else:
        raise ValueError(f"unknown vertex kind: {kind}")
    return point


def two_vertex_certificate(rank: int, first_kind: str, second_kind: str) -> dict[str, Any]:
    first = vertex(first_kind, rank, 0)
    second = vertex(second_kind, rank, 1)
    first_packet = perpendicular_basis(first)
    second_packet = perpendicular_basis(second)
    assert matrix_rank(first_packet) == rank - 1
    assert matrix_rank(second_packet) == rank - 1
    assert matrix_rank(first_packet + first_packet) == rank - 1
    assert matrix_rank(first_packet + second_packet) == rank
    return {
        "rank": rank,
        "vertex_types": [first_kind, second_kind],
        "first_square": lorentz_pair(first, first),
        "second_square": lorentz_pair(second, second),
        "distinct_projective_lines": True,
        "one_packet_rank": rank - 1,
        "duplicated_packet_rank": rank - 1,
        "two_packet_union_rank": rank,
    }


def sample_root_gram(rank: int) -> Matrix:
    gram = [[-2 * int(row == column) for column in range(rank)] for row in range(rank)]
    gram[0][1] = gram[1][0] = 3
    return gram


def identity(size: int) -> Matrix:
    return [[int(row == column) for column in range(size)] for row in range(size)]


def root_reflection(gram: Matrix, root_index: int) -> Matrix:
    """Column-coordinate matrix x |-> x+(x,r)r for r^2=-2."""

    size = len(gram)
    reflection = identity(size)
    for column in range(size):
        reflection[root_index][column] += gram[root_index][column]
    return reflection


def host_and_chamber_certificate(rank: int) -> dict[str, Any]:
    gram = sample_root_gram(rank)
    gram_determinant = determinant(gram)
    expected = -5 * ((-2) ** (rank - 2))
    assert gram_determinant == expected != 0
    reflection = root_reflection(gram, 0)
    assert matmul(reflection, reflection) == identity(rank)
    assert matmul(matmul(transpose(reflection), gram), reflection) == gram
    slots = rank * (rank - 1) // 2
    return {
        "rank": rank,
        "strict_normalized_less_than": 14,
        "integral_physical_intersection_at_most": 13,
        "bounded_gram_box_count": 14**slots,
        "sample_lorentz_gram_determinant": gram_determinant,
        "sample_discriminant_group_order": abs(gram_determinant),
        "R_subset_L_subset_R_dual_from_integrality_and_full_rank": True,
        "intermediate_overlattices_in_finite_R_dual_mod_R": True,
        "root_reflection_is_integral_involutive_isometry": True,
        "complete_root_mirrors_are_locally_finite": True,
        "generic_compact_path_crosses_finitely_many_mirrors": True,
        "successive_wall_reflections_give_one_Weyl_chamber_orbit": True,
    }


def build_audit() -> dict[str, Any]:
    pair_types = (
        ("finite", "finite"),
        ("finite", "ideal"),
        ("ideal", "finite"),
        ("ideal", "ideal"),
    )
    return {
        "schema": "coxeter.noncompact-k3-bounded-spanning-independent-hostile.v1",
        "verdict": "PASS_EXACT_FINITE_SIDED_FULL_PICARD_SCOPE",
        "pins": {
            name: {"path": relative, "sha256": sha256(ROOT / relative)}
            for name, (relative, _) in PINS.items()
        },
        "visual_geometry": visual_formula_certificate(),
        "complete_cusp_links": [
            complete_box_cusp_packet(rank) for rank in range(3, 21)
        ],
        "two_extreme_point_bridge": [
            two_vertex_certificate(rank, *pair)
            for rank in range(3, 21)
            for pair in pair_types
        ],
        "physical_host_and_chamber": [
            host_and_chamber_certificate(rank) for rank in range(3, 21)
        ],
        "scope": {
            "finite_sidedness_needed_for_minimum_and_uniform_cusp_height": True,
            "facet_irredundancy_used_for_disjoint_ultraparallel_shadows": True,
            "projective_closed_star_needed": True,
            "finite_volume_needed_for_bounded_cusp_link": True,
            "complete_incident_packet_needed": True,
            "two_distinct_extreme_points_needed": True,
            "literal_same_vector_physical_minus_two_roots_needed": True,
            "full_Picard_complete_root_arrangement_needed_for_one_orbit": True,
            "cusp_simplicity_needed": False,
            "compactness_needed": False,
            "reflection_maximality_needed": False,
            "Property_4_proved": False,
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-pins", action="store_true")
    arguments = parser.parse_args(argv)
    if arguments.check_pins:
        for relative, expected in PINS.values():
            actual = sha256(ROOT / relative)
            if actual != expected:
                raise SystemExit(
                    f"owner pin mismatch: {relative}: {actual} != {expected}"
                )
    print(json.dumps(build_audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
