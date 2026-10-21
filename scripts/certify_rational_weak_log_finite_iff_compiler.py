#!/usr/bin/env python3
"""Exact route replay for the finite weak-log rational iff compiler.

The companion report proves the universal cluster/Hilbert constructibility
and exact-return theorem.  This script checks the two new numerical inputs:
Fujino's fixed big-branch multiplier and the positive affine-cusp divisor on
the standard rank-ten anticanonical cusp.  It also freezes the logical route
and its deliberately limited scope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_weak_log_finite_iff_compiler.json"
SCHEMA = "rational_weak_log_finite_iff_compiler/v1"


def _strict_nonnegative_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def ceil_fraction(value: Fraction) -> int:
    if value < 0:
        raise ValueError("a must be nonnegative")
    return -(-value.numerator // value.denominator)


def fujino_multiplier(dimension: int, a: Fraction | int) -> int:
    """Return 2^(n+1)(n+1)!(ceil(a)+n), Fujino Remark 1.2."""

    n = _strict_nonnegative_integer(dimension, "dimension")
    parameter = a if isinstance(a, Fraction) else Fraction(a)
    return 2 ** (n + 1) * math.factorial(n + 1) * (ceil_fraction(parameter) + n)


def weak_log_big_cartier_multiplier() -> int:
    """Multiplier of D when L_D=2D and (n,a)=(2,0)."""

    return 2 * fujino_multiplier(2, Fraction(0))


def lorentz_dot(left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != len(right) or not left:
        raise ValueError("Lorentz vectors must have equal positive length")
    return left[0] * right[0] - sum(
        x * y for x, y in zip(left[1:], right[1:])
    )


def _solve_full_column_rank(
    columns: Sequence[Sequence[int]], target: Sequence[int]
) -> tuple[Fraction, ...]:
    """Solve sum x_j columns[j] = target by exact RREF.

    The system may be overdetermined.  A unique consistent solution is
    required, which is the situation for a minimal affine cusp block.
    """

    if not columns:
        raise ValueError("at least one column is required")
    row_count = len(target)
    column_count = len(columns)
    if any(len(column) != row_count for column in columns):
        raise ValueError("column lengths must match the target")
    matrix = [
        [Fraction(columns[column][row]) for column in range(column_count)]
        + [Fraction(target[row])]
        for row in range(row_count)
    ]
    pivot_row = 0
    pivot_columns: list[int] = []
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
    for row in matrix:
        if all(entry == 0 for entry in row[:-1]) and row[-1] != 0:
            raise ValueError("inconsistent cusp relation")
    if len(pivot_columns) != column_count:
        raise ValueError("cusp columns do not have full column rank")
    solution = [Fraction(0) for _ in range(column_count)]
    for row, column in enumerate(pivot_columns):
        solution[column] = matrix[row][-1]
    return tuple(solution)


def positive_cusp_relation(
    target: Sequence[int], normals: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    """Return primitive positive integral coefficients for target.

    The input is a chosen irreducible affine block.  It must consist of
    negative Coxeter normals orthogonal to a nonzero null target.
    """

    if lorentz_dot(target, target) != 0 or not any(target):
        raise ValueError("the cusp target must be nonzero and null")
    for normal in normals:
        if lorentz_dot(normal, normal) >= 0:
            raise ValueError("cusp normals must be negative")
        if lorentz_dot(target, normal) != 0:
            raise ValueError("every selected normal must be incident to the cusp")
    for i, left in enumerate(normals):
        for right in normals[i + 1 :]:
            if lorentz_dot(left, right) < 0:
                raise ValueError("off-diagonal Coxeter intersections must be nonnegative")
    rational = _solve_full_column_rank(normals, target)
    if any(coefficient <= 0 for coefficient in rational):
        raise ValueError("the affine cusp relation is not strictly positive")
    denominator = math.lcm(*(coefficient.denominator for coefficient in rational))
    integral = [int(coefficient * denominator) for coefficient in rational]
    divisor = math.gcd(*integral)
    primitive = tuple(coefficient // divisor for coefficient in integral)
    reconstructed = tuple(
        sum(primitive[j] * normals[j][row] for j in range(len(normals)))
        for row in range(len(target))
    )
    # Primitive coefficients may reconstruct a positive integral multiple or
    # fraction of target.  Collinearity and future orientation are the exact
    # invariant checks.
    nonzero_index = next(index for index, value in enumerate(target) if value)
    ratio = Fraction(reconstructed[nonzero_index], target[nonzero_index])
    if ratio <= 0 or any(
        Fraction(value) != ratio * target[index]
        for index, value in enumerate(reconstructed)
    ):
        raise ArithmeticError("positive coefficients do not reconstruct the cusp ray")
    return primitive


def gram(vectors: Sequence[Sequence[int]]) -> list[list[int]]:
    return [[lorentz_dot(left, right) for right in vectors] for left in vectors]


def physical_square_assignment_bound(facet_count: int) -> int:
    count = _strict_nonnegative_integer(facet_count, "facet_count")
    return 3**count


def sequential_program_bound(facet_count: int, picard_rank: int) -> int:
    facets = _strict_nonnegative_integer(facet_count, "facet_count")
    rank = _strict_nonnegative_integer(picard_rank, "picard_rank")
    depth = min(facets, max(rank - 1, 0))
    return sum(math.perm(facets, step) for step in range(depth + 1))


def determinant_3x3(matrix: Sequence[Sequence[int]]) -> int:
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise ValueError("a 3x3 matrix is required")
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def delta246_scaled_form(x: int, y: int, z: int) -> int:
    """Return -v^t G v / 2 for the pinned (2,4,6) Gram."""

    if any(isinstance(value, bool) or not isinstance(value, int) for value in (x, y, z)):
        raise ValueError("form coordinates must be integers")
    return 6 * x * x - 6 * x * y + 3 * y * y - 3 * y * z + z * z


def delta246_mod3_descent_certificate() -> dict[str, Any]:
    first_zeroes = [
        [x, y, z]
        for x in range(3)
        for y in range(3)
        for z in range(3)
        if delta246_scaled_form(x, y, z) % 3 == 0
    ]
    if any(z != 0 for _, _, z in first_zeroes):
        raise ArithmeticError("the first mod-3 step no longer forces z=0")
    binary_zeroes = [
        [x, y]
        for x in range(3)
        for y in range(3)
        if (2 * x * x + x * y + y * y) % 3 == 0
    ]
    if binary_zeroes != [[0, 0]]:
        raise ArithmeticError("the descended binary form became isotropic mod 3")
    return {
        "first_reduction": "q(x,y,z)=z^2 mod 3, so 3 divides z",
        "first_zero_residues": first_zeroes,
        "divided_reduction": "2x^2+xy+y^2=0 mod 3",
        "divided_binary_zero_residues": binary_zeroes,
        "conclusion": "3 divides x,y,z; primitive isotropic vector impossible",
    }


def _standard_p9_cusp() -> tuple[tuple[int, ...], list[tuple[int, ...]], list[str]]:
    # Coordinates are (H,E_1,...,E_9), with form diag(1,-1^9).
    target = (3,) + (-1,) * 9
    roots: list[tuple[int, ...]] = []
    names: list[str] = []
    for index in range(1, 9):
        vector = [0] * 10
        vector[index] = 1
        vector[index + 1] = -1
        roots.append(tuple(vector))
        names.append(f"d{index}")
    cremona = [0] * 10
    cremona[0] = 1
    cremona[1] = cremona[2] = cremona[3] = -1
    roots.append(tuple(cremona))
    names.append("c")
    return target, roots, names


def sha256_json(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build() -> dict[str, Any]:
    multiplier = fujino_multiplier(2, Fraction(0))
    if multiplier != 96 or weak_log_big_cartier_multiplier() != 192:
        raise ArithmeticError("the fixed effective base-point-free multiplier changed")

    target, roots, names = _standard_p9_cusp()
    coefficients = positive_cusp_relation(target, roots)
    coefficient_map = dict(zip(names, coefficients))
    expected = {
        "d1": 2,
        "d2": 4,
        "d3": 6,
        "d4": 5,
        "d5": 4,
        "d6": 3,
        "d7": 2,
        "d8": 1,
        "c": 3,
    }
    if coefficient_map != expected:
        raise ArithmeticError("the affine E8 cusp divisor changed")
    if any(lorentz_dot(target, root) != 0 for root in roots):
        raise ArithmeticError("a selected P9 root is not incident to the cusp")
    reconstructed = tuple(
        sum(coefficients[j] * roots[j][row] for j in range(len(roots)))
        for row in range(len(target))
    )
    if reconstructed != target:
        raise ArithmeticError("the P9 affine root divisor no longer equals -K")

    p9_gram = gram(roots)
    if any(p9_gram[i][i] != -2 for i in range(len(roots))):
        raise ArithmeticError("the P9 affine block root norm changed")
    if any(
        p9_gram[i][j] < 0
        for i in range(len(roots))
        for j in range(len(roots))
        if i != j
    ):
        raise ArithmeticError("the P9 affine block lost its Coxeter signs")

    delta246_gram = [[-12, 6, 0], [6, -6, 3], [0, 3, -2]]
    if determinant_3x3(delta246_gram) != 36:
        raise ArithmeticError("the (2,4,6) Gram determinant changed")
    delta246_descent = delta246_mod3_descent_certificate()
    i12_null = (1, 1, 0)
    if lorentz_dot(i12_null, i12_null) != 0:
        raise ArithmeticError("the I_(1,2) isotropic control changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": {
            "name": "RLogComp finite iff/EMPTY theorem",
            "base_field": "C",
            "input": (
                "finite-sided finite-volume arithmetic-reflection-maximal "
                "Coxeter chamber"
            ),
            "output": "finite weak-log rational carrier stack or EMPTY",
            "equivalence": (
                "output is nonempty iff a full-Picard complex rational weak-log "
                "carrier with simultaneous prime facets and exact full-nef return exists"
            ),
        },
        "finite_stages": [
            "HOST: 3^f physical square assignments and finite Picard overlattices",
            "MMP: finite facet-supported contraction programs to P2 or F_e",
            "PRIME: simultaneous geometrically integral relative Hilbert loci",
            "NEF: D=0, positive affine cusp divisor, or fixed |192D| gate",
            "RETURN: weak-log reflection saturation and ample orientation",
        ],
        "host_finiteness": {
            "physical_squares": [-1, -2, -4],
            "sample_facet_count": 12,
            "sample_assignment_bound": physical_square_assignment_bound(12),
            "overlattice_interval": "R subset L subset R*",
            "same_rational_quadratic_space_required": True,
            "rational_isotropy_and_local_invariants_preserved": True,
            "canonical_equations": "K.C_i=-2-C_i^2, K characteristic, K^2=10-rho",
        },
        "mmp_finiteness": {
            "sample_facets": 12,
            "sample_picard_rank": 11,
            "coarse_ordered_program_bound": sequential_program_bound(12, 11),
            "boundary_pivots_after_square_change_allowed": True,
            "endpoints": ["P2", "F_e"],
        },
        "nef_routes": {
            "zero": {
                "condition": "D=0",
                "certificate": "tautological nefness",
            },
            "isotropic": {
                "condition": "D nonzero and D^2=0",
                "certificate": "positive effective affine-cusp facet divisor mD",
                "semiampleness_required": False,
            },
            "big": {
                "condition": "D^2>0",
                "cartier_divisor": "L_D=2D",
                "fujino_dimension": 2,
                "fujino_a": 0,
                "fujino_multiplier_of_L_D": multiplier,
                "fixed_basepoint_free_system": "|96(2D)|=|192D|",
                "multiplier_of_D": weak_log_big_cartier_multiplier(),
            },
        },
        "p9_affine_cusp_replay": {
            "lattice": "I_(1,9)",
            "D_minus_K": list(target),
            "D_square": lorentz_dot(target, target),
            "root_order": names,
            "root_coefficients": list(coefficients),
            "coefficient_map": coefficient_map,
            "reconstructed_D": list(reconstructed),
            "root_gram": p9_gram,
            "all_coefficients_positive": all(value > 0 for value in coefficients),
            "all_roots_incident": all(lorentz_dot(target, root) == 0 for root in roots),
        },
        "p10_p11_routing": {
            "decorated_P10": {
                "physical_special_facet": "B=2(-K), reduced rational (-4) requested",
                "D": 0,
                "maximality_forces_prime_effectivity": False,
                "weak_log_compiler_expected_answer_over_C": "EMPTY",
            },
            "uncoloured_P10": {
                "physical_special_facet": "F=-K, primitive genus-one (-1)",
                "weak_log_smooth_rational_species": False,
                "known_exact_flex_carrier": True,
                "hyperbolic_dimension": 10,
            },
            "uncoloured_P11": {
                "physical_special_facet": "F=-K, primitive genus-one (-2)",
                "weak_log_smooth_rational_species": False,
                "known_exact_flex_carrier": True,
                "hyperbolic_dimension": 11,
            },
        },
        "raw_maximal_rational_log_no_go": {
            "coxeter_triangle": [2, 5, 6],
            "compact": True,
            "absolute_reflection_maximality": True,
            "minimal_vinberg_field": "Q(sqrt(5))",
            "literal_rational_full_nef_field": "Q",
            "every_rational_or_rational_log_full_nef_carrier_empty": True,
            "other_carrier_species_excluded": False,
            "companion_report": "results/k3_compact_maximal_256_phys_empty.md",
        },
        "raw_maximal_rational_field_no_go": {
            "coxeter_triangle": [2, 4, 6],
            "compact": True,
            "absolute_reflection_maximality": True,
            "minimal_vinberg_field": "Q",
            "gram": delta246_gram,
            "gram_determinant": determinant_3x3(delta246_gram),
            "scaled_form": "6x^2-6xy+3y^2-3yz+z^2",
            "anisotropy_certificate": delta246_descent,
            "rank3_rational_picard_lattice": "I_(1,2)",
            "rank3_rational_picard_null_vector": list(i12_null),
            "host_gate": "EMPTY: anisotropic raw Q-space cannot be I_(1,2)",
            "geometric_gate": "EMPTY: compact chamber cannot contain the rational ruling cusp",
            "every_rational_or_rational_log_full_nef_carrier_empty": True,
            "companion_report": "results/compact_246_raw_maximal_property4_counterexample.md",
        },
        "consequences_on_nonempty_output": {
            "dimension_bound": 25,
            "minimal_vinberg_field": "Q",
            "boundary_free_dimension_bound": 9,
            "carrier_relative_maximal_chamber_finiteness": True,
        },
        "scope": {
            "universal_property4_claimed": False,
            "universal_nonemptiness_claimed": False,
            "hostile_direct_full_ns_lane_closed": False,
            "genus_one_facet_lane_included": False,
            "generic_groebner_elimination_implemented": False,
            "finite_iff_criterion_proved": True,
            "maximality_used_to_create_sections": False,
            "maximality_used_only_after_prime_and_nef_gates": True,
        },
    }
    payload["content_sha256"] = sha256_json(payload)
    return payload


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = build()
    if args.check_archive:
        if not ARCHIVE.exists():
            raise SystemExit(f"missing archive: {ARCHIVE}")
        frozen = json.loads(ARCHIVE.read_text())
        if frozen != payload:
            raise SystemExit("live weak-log iff compiler certificate differs from archive")
    print(json.dumps(payload, indent=None if args.compact else 2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
