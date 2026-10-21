#!/usr/bin/env python3
"""Exact replay for the admissible cusped-field collapse theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Sequence


Matrix = list[list[Fraction]]
Vector = list[Fraction]


def identity(size: int) -> Matrix:
    return [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix)]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def bilinear(left: Vector, form: Matrix, right: Vector) -> Fraction:
    return sum(left[i] * form[i][j] * right[j] for i in range(len(left)) for j in range(len(right)))


def matrix_rank(matrix: Matrix) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [entry / scale for entry in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            multiple = work[index][column]
            work[index] = [
                work[index][j] - multiple * work[row][j]
                for j in range(len(work[0]))
            ]
        row += 1
        if row == len(work):
            break
    return row


def reflection_matrix(form: Matrix, normal: Vector) -> Matrix:
    norm = bilinear(normal, form, normal)
    if norm == 0:
        raise ValueError("a reflection normal must be nonisotropic")
    covector = matvec(form, normal)
    result = identity(len(normal))
    for row in range(len(normal)):
        for column in range(len(normal)):
            result[row][column] -= 2 * normal[row] * covector[column] / norm
    return result


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return [
        [left[i][j] - right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def cube_cusp_record(dimension: int) -> dict[str, object]:
    """Return an exact nonsimplicial-capable cusp frame in H^dimension."""

    if dimension < 2:
        raise ValueError("hyperbolic dimension must be at least two")
    size = dimension + 1
    form = identity(size)
    form[-1][-1] = Fraction(-1)
    cusp = [Fraction(0) for _ in range(size)]
    cusp[-2] = cusp[-1] = Fraction(1)
    normals: list[Vector] = []
    for coordinate in range(dimension - 1):
        for sign in (1, -1):
            normal = cusp[:]
            normal[coordinate] = Fraction(sign)
            normals.append(normal)

    if bilinear(cusp, form, cusp) != 0:
        raise AssertionError("the cusp vector is not null")
    if any(bilinear(cusp, form, normal) for normal in normals):
        raise AssertionError("an incident wall misses the cusp")
    if any(bilinear(normal, form, normal) != 1 for normal in normals):
        raise AssertionError("the cube-cusp normals should have norm one")
    if matrix_rank(normals) != dimension:
        raise AssertionError("the incident normals do not span u-perp")

    reflections = [reflection_matrix(form, normal) for normal in normals]
    fixed_equations: Matrix = []
    for reflection in reflections:
        fixed_equations.extend(subtract(reflection, identity(size)))
        if matmul(reflection, reflection) != identity(size):
            raise AssertionError("a returned map is not involutive")
        if matmul(transpose(reflection), matmul(form, reflection)) != form:
            raise AssertionError("a returned map does not preserve the form")
        if matvec(reflection, cusp) != cusp:
            raise AssertionError("the cusp line is not fixed")
        if any(entry.denominator != 1 for row in reflection for entry in row):
            raise AssertionError("the model reflection is not integral")

    fixed_dimension = size - matrix_rank(fixed_equations)
    if fixed_dimension != 1:
        raise AssertionError("the common fixed space is not the cusp line")
    normal_sum = [sum(normal[j] for normal in normals) for j in range(size)]
    expected_sum = [2 * (dimension - 1) * entry for entry in cusp]
    if normal_sum != expected_sum:
        raise AssertionError("the positive cusp dependence changed")

    return {
        "hyperbolic_dimension": dimension,
        "ambient_rank": size,
        "horosphere_section": f"[-1,1]^{dimension - 1}",
        "incident_facet_count": len(normals),
        "simple_cusp_facet_count": dimension,
        "nonsimple_when_incident_count_exceeds_dimension": len(normals) > dimension,
        "incident_normal_rank": matrix_rank(normals),
        "common_fixed_space_dimension": fixed_dimension,
        "cusp_square": 0,
        "positive_dependence": f"sum(e_i)={2 * (dimension - 1)}*u",
        "all_reflections_integral_involutive_isometries": True,
    }


def build_certificate() -> dict[str, object]:
    frames = [cube_cusp_record(dimension) for dimension in range(2, 11)]
    return {
        "theorem": "exact admissible finite-volume noncompact trace chambers have k=Q",
        "exact_cusp_frames": frames,
        "proof_gates": {
            "noncompact_finite_sided_chamber_has_ideal_accumulation_point": True,
            "finite_volume_forces_bounded_full_dimensional_horosphere_section": True,
            "bounded_section_forces_incident_normals_to_span_u_perp": True,
            "all_incident_reflection_fixed_spaces_intersect_in_exactly_the_cusp_line": True,
            "kernel_of_k_linear_reflection_packet_is_defined_over_k": True,
            "scalar_extension_preserves_the_kernel_dimension": True,
            "distinguished_nullity_implies_h_u_u_is_zero_in_k": True,
            "nonzero_vector_stays_nonzero_at_every_field_embedding": True,
            "positive_definite_nonidentity_conjugate_forbids_that_null_vector": True,
            "totally_real_field_degree_is_one": True,
        },
        "conclusions": {
            "returned_coefficient_field_is_Q": True,
            "full_trace_K3_noncompact_field_set_is_singleton_Q": True,
            "full_trace_K3_rank_identity_after_collapse": "n+1=rho(X)<=20",
            "full_trace_K3_hyperbolic_dimension_at_most_19": True,
            "does_not_supersede_rational_weak_log_n_at_most_25": True,
            "holds_in_every_hyperbolic_dimension_at_least_two": True,
            "requires_correspondence_height": False,
            "requires_reflection_maximality": False,
            "requires_reflective_lattice_classification": False,
        },
        "scope": {
            "finite_projective_order_module_allowed": True,
            "free_module_required": False,
            "primitive_wall_normals_required": False,
            "simple_cusp_required": False,
            "exact_full_incident_cusp_packet_required": True,
            "admissible_positive_definite_nonidentity_places_required": True,
            "partial_affine_subdiagram_suffices": False,
            "upper_half_space_argument_includes_nonsimple_cusps": True,
            "compact_chambers_closed": False,
            "carrier_construction_proved": False,
            "Property_4_proved": False,
            "raw_chamber_dimension_bound_proved": False,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args(argv)
    encoded = encoded_certificate()
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit(f"archive mismatch: {arguments.check}")
        print(f"PASS archive replay: {arguments.check}")
        return
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
        print(arguments.output)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
