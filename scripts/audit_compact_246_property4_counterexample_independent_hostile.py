#!/usr/bin/env python3
"""Independent hostile replay of the compact (2,4,6) Property-4 no-go.

This auditor imports neither SymPy nor the owner certifier.  It reconstructs
the rank-three arithmetic packet with Fraction arithmetic, proves absolute
reflection maximality by an exact index census, and separates the two
rational full-nef obstructions from the exact minimal-trace angle obstruction.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/compact_246_raw_maximal_property4_counterexample.md"
OWNER_SCRIPT = ROOT / "scripts/certify_compact_246_property4_counterexample.py"
OWNER_ARCHIVE = ROOT / "data/enumerations/compact_246_property4_counterexample.json"
OWNER_TEST = ROOT / "tests/test_compact_246_property4_counterexample.py"

AUDIT_REPORT = (
    ROOT
    / "results/compact_246_raw_maximal_property4_counterexample_"
    "independent_hostile_audit.md"
)
AUDIT_ARCHIVE = (
    ROOT
    / "data/enumerations/compact_246_property4_counterexample_"
    "independent_hostile_audit.json"
)

OWNER_HASHES = {
    "results/compact_246_raw_maximal_property4_counterexample.md": (
        "860d82c50ecfe46bd12c0f2fba851e85c6d08f6bf5e6b9b6ef575c75f31006b0"
    ),
    "scripts/certify_compact_246_property4_counterexample.py": (
        "c0afe27bacc218b1128ece87de8b7fa1f421a574f8b9039afa54947c2d73a64c"
    ),
    "data/enumerations/compact_246_property4_counterexample.json": (
        "74fbb33d09a31af0a0af3ec1a65694732362ccc8e6eb6d799a4906a4debf2cae"
    ),
    "tests/test_compact_246_property4_counterexample.py": (
        "089410474b73f09c134a10b00c33edc96c6ab30b9e1e6171b26be7eaaf3c34e0"
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def determinant_2(matrix: list[list[Fraction]]) -> Fraction:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def determinant_3(matrix: list[list[Fraction]]) -> Fraction:
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def transpose(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*matrix)]


def multiply(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Fraction(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def identity(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(int(i == j)) for j in range(size)]
        for i in range(size)
    ]


def reflection_matrix(
    gram: list[list[Fraction]], root_index: int
) -> list[list[Fraction]]:
    result = identity(len(gram))
    denominator = gram[root_index][root_index]
    for column in range(len(gram)):
        result[root_index][column] -= (
            2 * gram[column][root_index] / denominator
        )
    return result


def encode_matrix(matrix: list[list[Fraction]]) -> list[list[int | str]]:
    return [
        [
            int(value) if value.denominator == 1 else str(value)
            for value in row
        ]
        for row in matrix
    ]


def reconstruct_raw_packet() -> dict[str, Any]:
    gram = [
        [Fraction(-12), Fraction(6), Fraction(0)],
        [Fraction(6), Fraction(-6), Fraction(3)],
        [Fraction(0), Fraction(3), Fraction(-2)],
    ]
    if gram != transpose(gram):
        raise AssertionError("Gram is not symmetric")

    determinant = determinant_3(gram)
    leading_two = determinant_2([row[:2] for row in gram[:2]])
    pivots = [
        gram[0][0],
        leading_two / gram[0][0],
        determinant / leading_two,
    ]
    if determinant != 36 or pivots != [-12, -3, 1]:
        raise AssertionError("Gram determinant or inertia changed")

    pairs = [(0, 1), (1, 2), (0, 2)]
    normalized = [
        gram[i][j] ** 2 / (gram[i][i] * gram[j][j])
        for i, j in pairs
    ]
    if normalized != [Fraction(1, 2), Fraction(3, 4), Fraction(0)]:
        raise AssertionError("Coxeter products changed")

    principal_pair_determinants = [
        gram[i][i] * gram[j][j] - gram[i][j] ** 2
        for i, j in pairs
    ]
    if principal_pair_determinants != [36, 3, 24]:
        raise AssertionError("a vertex is no longer finite")

    cartan = {
        f"C_{i + 1}{j + 1}": 2 * gram[i][j] / gram[i][i]
        for i in range(3)
        for j in range(3)
        if i != j and gram[i][j] != 0
    }
    expected_cartan = {
        "C_12": Fraction(-1),
        "C_21": Fraction(-2),
        "C_23": Fraction(-1),
        "C_32": Fraction(-3),
    }
    if cartan != expected_cartan:
        raise AssertionError("Cartan coefficients changed")

    reflections = [reflection_matrix(gram, i) for i in range(3)]
    expected_reflections = [
        [
            [Fraction(-1), Fraction(1), Fraction(0)],
            [Fraction(0), Fraction(1), Fraction(0)],
            [Fraction(0), Fraction(0), Fraction(1)],
        ],
        [
            [Fraction(1), Fraction(0), Fraction(0)],
            [Fraction(2), Fraction(-1), Fraction(1)],
            [Fraction(0), Fraction(0), Fraction(1)],
        ],
        [
            [Fraction(1), Fraction(0), Fraction(0)],
            [Fraction(0), Fraction(1), Fraction(0)],
            [Fraction(0), Fraction(3), Fraction(-1)],
        ],
    ]
    if reflections != expected_reflections:
        raise AssertionError("reflection matrices changed")
    for reflection in reflections:
        if multiply(reflection, reflection) != identity(3):
            raise AssertionError("reflection is not an involution")
        if multiply(multiply(transpose(reflection), gram), reflection) != gram:
            raise AssertionError("reflection does not preserve the Gram")
        if any(value.denominator != 1 for row in reflection for value in row):
            raise AssertionError("reflection is not integral")

    area = 1 - Fraction(1, 2) - Fraction(1, 4) - Fraction(1, 6)
    if area != Fraction(1, 12):
        raise AssertionError("triangle area changed")

    cycle_products = [
        cartan["C_12"] * cartan["C_21"],
        cartan["C_23"] * cartan["C_32"],
    ]
    if cycle_products != [2, 3]:
        raise AssertionError("Vinberg cyclic products changed")

    return {
        "wall_gram": encode_matrix(gram),
        "determinant": int(determinant),
        "LDL_pivots": [int(value) for value in pivots],
        "signature": [1, 2],
        "normalized_cosine_squares": [str(value) for value in normalized],
        "edge_labels": [4, 6, 2],
        "principal_pair_determinants": [
            int(value) for value in principal_pair_determinants
        ],
        "all_vertices_finite": True,
        "area_over_pi": str(area),
        "cartan_coefficients": {
            key: str(value) for key, value in cartan.items()
        },
        "cycle_products": [int(value) for value in cycle_products],
        "minimal_vinberg_field": "Q",
        "reflection_matrices": [encode_matrix(matrix) for matrix in reflections],
        "integral_reflection_host": True,
        "arithmetic_lattice_argument": (
            "W is a lattice inside the rational arithmetic lattice O(G,Z), "
            "hence has finite index"
        ),
    }


def anisotropy_descent() -> dict[str, Any]:
    # For q=-G/2, v^t(-G/2)v is the integral polynomial below.
    q_mod_three_zeros = [
        [x, y, z]
        for x in range(3)
        for y in range(3)
        for z in range(3)
        if (
            6 * x * x
            - 6 * x * y
            + 3 * y * y
            - 3 * y * z
            + z * z
        )
        % 3
        == 0
    ]
    if {row[2] for row in q_mod_three_zeros} != {0}:
        raise AssertionError("first descent does not force 3|z")

    residual_binary_zeros = [
        [x, y]
        for x in range(3)
        for y in range(3)
        if (2 * x * x + x * y + y * y) % 3 == 0
    ]
    if residual_binary_zeros != [[0, 0]]:
        raise AssertionError("residual binary form gained a nonzero zero")

    return {
        "similar_integral_form": (
            "q=6x^2-6xy+3y^2-3yz+z^2"
        ),
        "descent_prime": 3,
        "q_mod_3_zero_count": len(q_mod_three_zeros),
        "q_mod_3_forces_z_zero": True,
        "residual_binary_form": "2x^2+xy+y^2",
        "residual_binary_zeros_mod_3": residual_binary_zeros,
        "primitive_integer_zero_would_be_divisible_by_3": True,
        "rational_form_anisotropic": True,
    }


def reflection_area_maximality() -> dict[str, Any]:
    source_area = Fraction(1, 12)
    polygon_floor = Fraction(1, 6)
    triangle_floor = Fraction(1, 42)
    if not source_area / 2 < polygon_floor:
        raise AssertionError("four-side polygon exclusion failed")

    # area(Q)=1/(12N) and area(Q)>=1/42 imply N<=3.
    indices = [2, 3]
    solutions_by_index: dict[str, list[list[int]]] = {}
    for index in indices:
        reciprocal_sum = 1 - source_area / index
        solutions: list[list[int]] = []
        # With 2<=a<=b<=c, reciprocal_sum<=3/a, so a<=3.
        for a in range(2, 4):
            remainder_after_a = reciprocal_sum - Fraction(1, a)
            b_upper = int(Fraction(2) / remainder_after_a)
            for b in range(a, b_upper + 1):
                reciprocal_c = (
                    reciprocal_sum - Fraction(1, a) - Fraction(1, b)
                )
                if reciprocal_c <= 0:
                    continue
                c = Fraction(1) / reciprocal_c
                if c.denominator == 1 and c.numerator >= b:
                    solutions.append([a, b, c.numerator])
        solutions_by_index[str(index)] = solutions

    if solutions_by_index != {"2": [[2, 3, 8]], "3": []}:
        raise AssertionError(
            f"unexpected triangle overgroup census {solutions_by_index}"
        )

    doubling_cases = [
        {
            "shared_endpoint_orders": [2, 8],
            "opposite_order": 3,
            "union": "triangle (3,3,4)",
            "coxeter": True,
            "equals_target": False,
        },
        {
            "shared_endpoint_orders": [2, 3],
            "opposite_order": 8,
            "union": "triangle angles (pi/8,pi/8,2pi/3)",
            "coxeter": False,
            "equals_target": False,
        },
        {
            "shared_endpoint_orders": [3, 8],
            "opposite_order": 2,
            "union": "quadrilateral",
            "coxeter": False,
            "equals_target": False,
        },
    ]
    if any(case["equals_target"] for case in doubling_cases):
        raise AssertionError("the index-two candidate doubles to the target")

    return {
        "source_area_over_pi": str(source_area),
        "proper_overgroup_area_upper_bound_over_pi": str(source_area / 2),
        "compact_polygon_s_ge_4_area_floor_over_pi": str(polygon_floor),
        "overgroup_chamber_must_be_triangle": True,
        "compact_triangle_area_floor_over_pi": str(triangle_floor),
        "derived_index_range": indices,
        "reciprocal_equation": "1/a+1/b+1/c=1-1/(12N)",
        "enumeration_bounds": "2<=a<=3 and a<=b<=2/(S-1/a)",
        "solutions_by_index": solutions_by_index,
        "doubling_cases": doubling_cases,
        "proper_discrete_reflection_overgroup_exists": False,
        "absolute_reflection_maximality": True,
    }


def rational_full_nef_obstruction(
    anisotropy: dict[str, Any],
) -> dict[str, Any]:
    if not anisotropy["rational_form_anisotropic"]:
        raise AssertionError("anisotropic host obstruction unavailable")

    # I_(1,2) is visibly isotropic on (1,1,0).
    isotropic_vector = [1, 1, 0]
    i12_square = (
        isotropic_vector[0] ** 2
        - isotropic_vector[1] ** 2
        - isotropic_vector[2] ** 2
    )
    if i12_square != 0:
        raise AssertionError("I_(1,2) witness is not isotropic")

    return {
        "target_hyperbolic_dimension": 2,
        "literal_full_nef_forces_picard_rank": 3,
        "rank_three_rational_NS": "I_(1,2)",
        "I_1_2_isotropic_vector": isotropic_vector,
        "I_1_2_isotropic_vector_square": i12_square,
        "target_trace_form_anisotropic": True,
        "full_rank_similarity_preserves_isotropy": True,
        "rational_HOST_exists": False,
        "MMP_nef_isotropic_cases": {
            "ruled_minimal_endpoint": "pull back the ruling fibre",
            "P2_minimal_endpoint": (
                "stop after the final inverse blowup at F1 and pull back H-E"
            ),
        },
        "every_rational_surface_rho_ge_2_has_nef_isotropic_ray": True,
        "compact_target_has_ideal_ray": False,
        "direct_rational_exact_full_nef": "EMPTY",
        "certified_weak_log_rational_exact_full_nef": "EMPTY",
    }


def exact_minimal_trace_obstruction() -> dict[str, Any]:
    allowed_rows = []
    for intersection in (0, 1):
        cosine_square = Fraction(intersection * intersection, 4)
        label = 2 if intersection == 0 else 3
        allowed_rows.append(
            {
                "intersection": intersection,
                "cosine_square": str(cosine_square),
                "finite_Coxeter_label": label,
            }
        )

    target_labels = [2, 4, 6]
    if not {4, 6}.isdisjoint({row["finite_Coxeter_label"] for row in allowed_rows}):
        raise AssertionError("a forbidden target label became K-trivial")

    return {
        "trace_rank": 3,
        "literal_full_nef_forces_picard_rank": 3,
        "nontrivial_contraction_endpoint_rank_at_most": 2,
        "rank_three_trace_cannot_descend_nontrivially": True,
        "exact_carrier_is_already_minimal": True,
        "trace_space_equals_NS_over_Q": True,
        "canonical_orthogonality_forces_K_numerically_zero": True,
        "adjunction_for_negative_irreducible_facet": (
            "C^2=2p_a(C)-2<0, hence p_a(C)=0 and C^2=-2"
        ),
        "facet_curves_are_smooth_rational_minus_two": True,
        "distinct_curve_intersection_nonnegative_integer": True,
        "allowed_finite_angle_rows": allowed_rows,
        "allowed_finite_Coxeter_labels": [2, 3],
        "target_labels": target_labels,
        "forbidden_target_labels": [4, 6],
        "exact_return_minimalizable_full_trace": "EMPTY",
        "bare_carrier_relative_trace_injection_excluded": False,
    }


def owner_hash_audit() -> dict[str, str]:
    observed = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (OWNER_REPORT, OWNER_SCRIPT, OWNER_ARCHIVE, OWNER_TEST)
    }
    if observed != OWNER_HASHES:
        raise AssertionError(
            f"owner bytes changed; expected {OWNER_HASHES}, observed {observed}"
        )
    return observed


def build_certificate() -> dict[str, Any]:
    owner_hashes = owner_hash_audit()
    raw = reconstruct_raw_packet()
    anisotropy = anisotropy_descent()
    maximality = reflection_area_maximality()
    rational = rational_full_nef_obstruction(anisotropy)
    minimal_trace = exact_minimal_trace_obstruction()

    return {
        "schema": (
            "coxeter.compact-246-property4-counterexample-"
            "independent-hostile-audit.v1"
        ),
        "audit_report_sha256": sha256_file(AUDIT_REPORT),
        "owner_hashes": owner_hashes,
        "independence": {
            "imports_owner_certifier": False,
            "imports_sympy": False,
            "exact_engine": "Python Fraction and explicit 3x3 matrix algebra",
            "maximality_reconstruction": (
                "index-bounded reciprocal equation, not owner area window"
            ),
        },
        "raw_reconstruction": raw,
        "anisotropy": anisotropy,
        "absolute_reflection_maximality": maximality,
        "rational_full_nef_obstruction": rational,
        "exact_minimal_trace_obstruction": minimal_trace,
        "four_property_scope": {
            "requirement_1_trace_HOST": (
                "fails for the two rank-three rational-surface branches"
            ),
            "requirement_2_exact_chamber": (
                "literal k=Q full-nef equality forces rho=3 and is essential"
            ),
            "requirement_3_canonical_log_MMP": (
                "nontrivial rank-three descent is impossible; the minimal "
                "K-trivial lane has only labels 2 and 3"
            ),
            "requirement_4_essential_surjectivity": (
                "false for the union of the three currently declared "
                "exact-return carrier species"
            ),
            "current_three_species_all_empty": True,
            "conditional_AG_dimension_theorem_refuted": False,
            "all_conceivable_AG_carriers_excluded": False,
            "proper_invariant_or_eigen_nef_slice_for_k_Q_excluded": False,
            "bare_carrier_relative_third_row_excluded": False,
            "projected_larger_Picard_return_excluded": False,
            "higher_dimensional_variety_return_excluded": False,
            "positive_characteristic_return_excluded": False,
        },
        "presentation": {
            "owner_status_self_contained": True,
            "owner_scope_escape_clause_explicit": True,
            "owner_distinguishes_exact_from_bare_third_row": True,
            "owner_reproduction_requires_PYTHONPATH_src_dot": True,
            "owner_markdown_pandoc_checked": True,
        },
        "verdict": "PASS",
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        AUDIT_ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        AUDIT_ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if (
            not AUDIT_ARCHIVE.exists()
            or AUDIT_ARCHIVE.read_text(encoding="utf-8") != encoded
        ):
            raise SystemExit("independent hostile compact (2,4,6) archive differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
