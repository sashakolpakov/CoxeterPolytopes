#!/usr/bin/env python3
"""Independent exact audit for the rational weak-log iff/EMPTY compiler.

This replay deliberately does not import the owner certifier.  It checks the
numerical claims used by the hostile audit, pins the audited owner bytes, and
records the finite-endpoint repair that is easy to miss in the abstract MMP
argument.  The companion markdown report contains the geometric proofs.
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
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_weak_log_finite_iff_compiler_independent_audit.json"
)
REPORT = ROOT / "results/rational_weak_log_finite_iff_compiler.md"
OWNER_SCRIPT = ROOT / "scripts/certify_rational_weak_log_finite_iff_compiler.py"
OWNER_ARCHIVE = ROOT / "data/enumerations/rational_weak_log_finite_iff_compiler.json"
OWNER_TEST = ROOT / "tests/test_rational_weak_log_finite_iff_compiler.py"
SCHEMA = "rational_weak_log_finite_iff_compiler_independent_audit/v1"

EXPECTED_OWNER_SHA256 = {
    "report": "48cf57f67dea15ed7a307047506e1d23ffd6841af6dcb29798d99e9eb5e34ffe",
    "script": "45de835b6a75f2abc0ff9d1bdd02314279ceee09db90f575ff4655a43569f5ee",
    "archive": "87f391c8c82f6ef9739c5b034b53c14f24dc70fba397aa8a52011c05dc31906d",
    "test": "f8f0111c2a955ec58f14ce5786a848bfcabba25c4b6c4ee8f65379395f49720c",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def owner_hashes() -> dict[str, str]:
    return {
        "report": sha256_file(REPORT),
        "script": sha256_file(OWNER_SCRIPT),
        "archive": sha256_file(OWNER_ARCHIVE),
        "test": sha256_file(OWNER_TEST),
    }


def sha256_json(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def rational_rank(rows: Sequence[Sequence[int | Fraction]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("ragged matrix")
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def lorentz_dot(left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != len(right) or not left:
        raise ValueError("Lorentz vectors must have equal positive length")
    return left[0] * right[0] - sum(
        x * y for x, y in zip(left[1:], right[1:])
    )


def standard_p9_affine_block() -> tuple[
    tuple[int, ...], list[tuple[int, ...]], tuple[int, ...]
]:
    target = (3,) + (-1,) * 9
    roots: list[tuple[int, ...]] = []
    for index in range(1, 9):
        root = [0] * 10
        root[index] = 1
        root[index + 1] = -1
        roots.append(tuple(root))
    cremona = [0] * 10
    cremona[0] = 1
    cremona[1] = cremona[2] = cremona[3] = -1
    roots.append(tuple(cremona))
    coefficients = (2, 4, 6, 5, 4, 3, 2, 1, 3)
    return target, roots, coefficients


def cusp_certificate() -> dict[str, Any]:
    target, roots, coefficients = standard_p9_affine_block()
    root_gram = [
        [lorentz_dot(left, right) for right in roots] for left in roots
    ]
    reconstruction = tuple(
        sum(coefficients[column] * roots[column][row] for column in range(9))
        for row in range(10)
    )
    coordinate_rows = [list(row) for row in zip(*roots)]
    checks = {
        "target_is_nonzero_null": any(target) and lorentz_dot(target, target) == 0,
        "root_coordinate_rank": rational_rank(coordinate_rows),
        "root_gram_rank": rational_rank(root_gram),
        "all_roots_have_square_minus_two": all(
            root_gram[index][index] == -2 for index in range(9)
        ),
        "off_diagonal_intersections_nonnegative": all(
            root_gram[i][j] >= 0
            for i in range(9)
            for j in range(9)
            if i != j
        ),
        "all_roots_incident_to_target": all(
            lorentz_dot(target, root) == 0 for root in roots
        ),
        "all_relation_coefficients_positive": all(value > 0 for value in coefficients),
        "positive_relation_reconstructs_target": reconstruction == target,
    }
    if checks != {
        "target_is_nonzero_null": True,
        "root_coordinate_rank": 9,
        "root_gram_rank": 8,
        "all_roots_have_square_minus_two": True,
        "off_diagonal_intersections_nonnegative": True,
        "all_roots_incident_to_target": True,
        "all_relation_coefficients_positive": True,
        "positive_relation_reconstructs_target": True,
    }:
        raise ArithmeticError("the independent affine-cusp certificate failed")
    return {
        "target": list(target),
        "coefficients": list(coefficients),
        "reconstruction": list(reconstruction),
        "checks": checks,
    }


def fujino_multiplier(dimension: int, a: Fraction | int) -> int:
    if isinstance(dimension, bool) or dimension < 0:
        raise ValueError("dimension must be a nonnegative integer")
    parameter = a if isinstance(a, Fraction) else Fraction(a)
    if parameter < 0:
        raise ValueError("a must be nonnegative")
    ceiling = -(-parameter.numerator // parameter.denominator)
    return 2 ** (dimension + 1) * math.factorial(dimension + 1) * (
        ceiling + dimension
    )


def hirzebruch_endpoint_pairs() -> list[list[int]]:
    """Return all (e,t) with e>=2 and -(e+t) in {-1,-2,-4}."""

    return [
        [e, t]
        for e in range(2, 5)
        for t in range(5)
        if e + t in {1, 2, 4}
    ]


def delta246_form(x: int, y: int, z: int) -> int:
    return 6 * x * x - 6 * x * y + 3 * y * y - 3 * y * z + z * z


def delta246_descent() -> dict[str, Any]:
    first = [
        [x, y, z]
        for x in range(3)
        for y in range(3)
        for z in range(3)
        if delta246_form(x, y, z) % 3 == 0
    ]
    second = [
        [x, y]
        for x in range(3)
        for y in range(3)
        if (2 * x * x + x * y + y * y) % 3 == 0
    ]
    if any(z for _, _, z in first) or second != [[0, 0]]:
        raise ArithmeticError("the (2,4,6) mod-3 descent failed")
    return {
        "q_zero_residues_mod_3": first,
        "first_step_forces_z_divisible_by_3": True,
        "descended_binary_zero_residues": second,
        "second_step_forces_x_y_divisible_by_3": True,
        "rational_form_is_anisotropic": True,
    }


def _integer_square_root(value: Fraction) -> int | None:
    if value.denominator != 1 or value < 0:
        return None
    root = math.isqrt(value.numerator)
    return root if root * root == value.numerator else None


def delta246_integral_physical_assignments() -> list[dict[str, int]]:
    """Enumerate physical squares for labels 4 and 6 with integral crossings."""

    survivors: list[dict[str, int]] = []
    for a in (1, 2, 4):
        for b in (1, 2, 4):
            for c in (1, 2, 4):
                intersection_12 = _integer_square_root(Fraction(a * b, 2))
                intersection_23 = _integer_square_root(Fraction(3 * b * c, 4))
                if intersection_12 is not None and intersection_23 is not None:
                    survivors.append(
                        {
                            "minus_C1_square": a,
                            "minus_C2_square": b,
                            "minus_C3_square": c,
                            "C1_C2": intersection_12,
                            "C2_C3": intersection_23,
                        }
                    )
    return survivors


def build() -> dict[str, Any]:
    hashes = owner_hashes()
    if hashes != EXPECTED_OWNER_SHA256:
        raise RuntimeError("audited owner bytes changed")

    endpoint_pairs = hirzebruch_endpoint_pairs()
    if endpoint_pairs != [[2, 0], [2, 2], [3, 1], [4, 0]]:
        raise ArithmeticError("the finite Hirzebruch endpoint list changed")

    multiplier = fujino_multiplier(2, 0)
    if multiplier != 96:
        raise ArithmeticError("Fujino's surface multiplier changed")

    physical_survivors = delta246_integral_physical_assignments()
    if physical_survivors:
        raise ArithmeticError("the (2,4,6) physical HOST gate is no longer empty")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "verdict": "PASS_AFTER_EXPLICIT_F_e_ENDPOINT_REPAIR",
        "audited_owner_sha256": hashes,
        "host_audit": {
            "physical_square_alphabet": [-1, -2, -4],
            "finite_assignment_count_for_f_facets": "3^f",
            "overlattice_interval": "R subset L subset R*",
            "rational_form_test_required": True,
            "real_signature_alone_sufficient": False,
            "delta246_integral_physical_assignments": physical_survivors,
            "delta246_descent": delta246_descent(),
        },
        "mmp_audit": {
            "later_exceptional_strict_transforms_are_original_negative_facets": True,
            "boundary_facets_may_be_later_pivots": True,
            "f1_is_not_minimal": True,
            "allowed_e_t_pairs_for_e_at_least_2": endpoint_pairs,
            "finite_hirzebruch_endpoints": [0, 2, 3, 4],
            "fixed_endpoint_ordered_cluster_schemes_are_finite_type": True,
            "automorphism_quotient_needed_for_finite_presentation": False,
        },
        "prime_audit": {
            "all_facets_parameterized_over_one_cluster_base": True,
            "geometrically_integral_fibre_locus_constructible": True,
            "adjunction_genus_zero_then_implies_smooth_rational": True,
            "zero_intersection_of_distinct_prime_divisors_implies_disjoint": True,
        },
        "nef_audit": {
            "zero_branch": "D=0 is nef",
            "isotropic_branch": cusp_certificate(),
            "big_branch": {
                "pair": "(S,B/2) is klt",
                "L_D": "2D=-2K-B is nef Cartier on an actual return",
                "a": 0,
                "aL_minus_K_plus_Delta": "D is nef and big",
                "fujino_multiplier_of_L_D": multiplier,
                "fixed_system": "|96L_D|=|192D|",
                "global_generation_is_necessary": True,
                "global_generation_is_sufficient_for_D_nef": True,
            },
        },
        "constructibility_audit": {
            "finite_cohomology_base_change_stratification": True,
            "evaluation_cokernel_support_closed_on_each_stratum": True,
            "proper_image_of_support_closed": True,
            "basepoint_free_locus_constructible": True,
            "algebraically_closed_nonemptiness_decidable_by_elimination": True,
            "generic_degree_192_elimination_implemented": False,
        },
        "return_and_scope_audit": {
            "adapted_host_reflection_maximality_used_only_after_prime_nef": True,
            "literal_full_nef_return_then_follows": True,
            "empty_means_no_carrier_in_this_weak_log_species": True,
            "empty_means_no_carrier_in_any_possible_species": False,
            "universal_property4_nonemptiness_proved": False,
            "conditional_dimension_and_field_transfer_only": True,
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
        if not ARCHIVE.is_file():
            raise SystemExit(f"missing archive: {ARCHIVE}")
        frozen = json.loads(ARCHIVE.read_text())
        if frozen != payload:
            raise SystemExit("live independent RLogComp audit differs from archive")
    print(json.dumps(payload, indent=None if args.compact else 2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
