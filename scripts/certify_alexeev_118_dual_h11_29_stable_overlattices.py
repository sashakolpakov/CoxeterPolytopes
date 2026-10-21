#!/usr/bin/env python3
"""Stable-overlattice and exact common-lineage certificate for H11/29.

The two 29-wall chambers come from Alexeev's dual 118-wall stocks.  This
script exhausts reflection-stable integral overlattices for both fixed
lattices, keeps global scaling separate from overclosure, and gives an
explicit integral isometry between the normalized 29-wall systems.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _canonical_sha256,
    _matrix,
    _overlattice_basis,
    _payload,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
OTHER_CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
OTHER_STABLE = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
OUTPUT = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json"
)
SCHEMA = "coxeter.alexeev_118_dual_h11_29.stable_overlattices.v1"

FACET_PERMUTATION = [
    1,
    19,
    2,
    18,
    10,
    11,
    15,
    3,
    4,
    12,
    27,
    7,
    6,
    22,
    8,
    24,
    14,
    13,
    17,
    5,
    9,
    25,
    20,
    23,
    29,
    16,
    26,
    28,
    21,
]

FORCED_ROOT_FACTORS = [
    "1/2",
    "1/2",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1/2",
    "1",
    "1",
    "1",
    "1/2",
    "1/2",
    "1",
    "1",
    "1/2",
    "1",
    "1",
    "1/2",
    "1/2",
    "1/2",
    "1",
    "1/2",
    "1",
    "1/2",
    "1/2",
    "1/2",
]

# A^T K_12 A=M_6, where K_12=G_12/2 and M_6 is the unique index-eight
# reflection-stable overclosure of the A16_6 fixed lattice.
COMMON_ISOMETRY_A = sp.Matrix(
    [
        [0, 0, 0, 0, -15, 0, -36, 0, -15, 0, 0, -50],
        [0, 0, 1, 0, 4, 0, 8, 0, 3, 0, 0, 10],
        [0, 0, 0, 0, -3, 0, -8, 1, -3, 0, 0, -10],
        [0, 0, 0, 1, -8, 0, -18, 0, -8, 0, 0, -25],
        [0, 0, 0, 0, 16, 0, 36, 0, 15, 0, 0, 50],
        [1, 0, 0, 0, 10, 0, 24, 0, 10, 0, 0, 33],
        [0, 1, 0, 0, 5, 0, 12, 0, 5, 0, 0, 16],
        [0, 0, 0, 0, -2, 1, -4, 0, -1, 0, 0, -4],
        [0, 0, 0, 0, -12, 0, -29, 0, -12, 0, 1, -41],
        [0, 0, 0, 0, -7, 0, -17, 0, -7, 1, 0, -24],
        [0, 0, 0, 0, 4, 0, 8, 0, 3, 0, 0, 11],
        [0, 0, 0, 0, 2, 0, 5, 0, 2, 0, 0, 7],
    ]
)

# P_6^T(-M_6)P_6 is Scharlau--Kirschmer Gram11 entry 5.
A16_6_TO_SK5 = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 4],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 4],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4],
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 4],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 4],
        [0, 1, 1, 1, 0, 0, 0, 0, 0, -1, -2, -4],
        [0, 0, 1, 1, 0, 0, 0, 0, 0, -1, -2, -4],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, -1, -2, -4],
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 4],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4],
    ]
)

SK_GRAM11_ENTRY5 = sp.diag(
    sp.eye(10), sp.Matrix([[0, -2], [-2, -4]])
)

EXPECTED_MAXIMAL = {
    "A16_6_1_C01": (-4, [1] * 10 + [2] * 2, False),
    "A16_12_1_C01": (-256, [1] * 4 + [2] * 8, False),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _active_roots(case: dict[str, Any]) -> list[sp.Matrix]:
    return [
        sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        for record in case["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]


def _oriented_active_roots(case: dict[str, Any]) -> list[sp.Matrix]:
    return [
        sp.Rational(record["presentation_scale_from_primitive"])
        * sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        for record in case["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]


def _stable_case(case: dict[str, Any]) -> tuple[dict[str, Any], sp.Matrix, sp.Matrix]:
    identifier = case["id"]
    gram = _matrix(case["saturated_fixed_lattice"]["gram"])
    roots = _active_roots(case)
    reflections = [_reflection(gram, root) for root in roots]
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    if audit["stable_integral_subgroup_orders"] != [1, 2, 4, 8]:
        raise ArithmeticError(f"stable integral subgroup census changed for {identifier}")
    if audit["stable_even_subgroup_orders"] != [1, 4]:
        raise ArithmeticError(f"stable even subgroup census changed for {identifier}")
    maximal_positions = audit["maximal_stable_subgroup_positions_zero_based"]
    if len(maximal_positions) != 1:
        raise ArithmeticError(f"maximal stable subgroup not unique for {identifier}")
    position = maximal_positions[0]
    group = groups[position]
    cosets = audit["stable_subgroups"][position]["cosets_in_smith_coordinates"]
    basis = _overlattice_basis(gram, left_inverse, cosets)
    maximal_gram = sp.simplify(basis.T * gram * basis)
    actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
    determinant, smith, even = EXPECTED_MAXIMAL[identifier]
    checks = {
        "unique_maximal_stable_subgroup": len(maximal_positions) == 1,
        "index_eight": abs(sp.Rational(1, basis.det())) == len(group) == 8,
        "determinant_square_law": maximal_gram.det() == gram.det() / 64,
        "expected_determinant": maximal_gram.det() == determinant,
        "expected_smith": _smith(maximal_gram) == smith,
        "expected_parity": (
            all(maximal_gram[i, i] % 2 == 0 for i in range(12)) == even
        ),
        "all_reflections_integral_on_maximal_overlattice": all(
            value.is_Integer for action in actions for value in action
        ),
        "all_reflections_unimodular": all(abs(action.det()) == 1 for action in actions),
        "all_reflections_preserve_maximal_gram": all(
            action.T * maximal_gram * action == maximal_gram for action in actions
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"maximal stable replay failed for {identifier}")
    return (
        {
            "id": identifier,
            "original_fixed_lattice": {
                "rank": gram.rows,
                "determinant": int(gram.det()),
                "smith_invariants": _smith(gram),
                "content": math.gcd(*(abs(int(value)) for value in gram)),
                "discriminant_q_mod_2_histogram": _q_histogram(gram),
            },
            "exhaustive_discriminant_module_audit": audit,
            "unique_maximal_reflection_stable_overlattice": {
                "index_over_original_lattice": 8,
                "basis_columns_in_original_fixed_lattice": _payload(basis),
                "gram": _payload(maximal_gram),
                "determinant": int(maximal_gram.det()),
                "smith_invariants": _smith(maximal_gram),
                "parity": "even" if even else "odd",
                "discriminant_q_mod_2_histogram": _q_histogram(maximal_gram),
                "facet_reflection_actions_sha256": _canonical_sha256(
                    [_payload(action) for action in actions]
                ),
                "checks": checks,
            },
        },
        basis,
        maximal_gram,
    )


def build_certificate() -> dict[str, Any]:
    source = json.loads(CONSTRUCTION.read_text())
    cases = {case["id"]: case for case in source["cases"]}
    records: dict[str, dict[str, Any]] = {}
    bases: dict[str, sp.Matrix] = {}
    maximal_grams: dict[str, sp.Matrix] = {}
    for identifier in ("A16_6_1_C01", "A16_12_1_C01"):
        record, basis, maximal_gram = _stable_case(cases[identifier])
        records[identifier] = record
        bases[identifier] = basis
        maximal_grams[identifier] = maximal_gram

    gram12 = _matrix(cases["A16_12_1_C01"]["saturated_fixed_lattice"]["gram"])
    normalized12 = gram12 / 2
    gram6_maximal = maximal_grams["A16_6_1_C01"]
    permutation = [value - 1 for value in FACET_PERMUTATION]
    factors = [sp.Rational(value) for value in FORCED_ROOT_FACTORS]
    walls6_original = _oriented_active_roots(cases["A16_6_1_C01"])
    walls6 = [bases["A16_6_1_C01"].inv() * wall for wall in walls6_original]
    walls12 = _oriented_active_roots(cases["A16_12_1_C01"])
    reflection6 = [_reflection(gram6_maximal, wall) for wall in walls6]
    reflection12 = [_reflection(normalized12, wall) for wall in walls12]
    wall_checks = {
        "facet_permutation_valid": sorted(permutation) == list(range(29)),
        "forced_factors_positive": all(factor > 0 for factor in factors),
        "all_29_oriented_walls_map_exactly": all(
            COMMON_ISOMETRY_A * (factors[index] * walls6[index])
            == walls12[permutation[index]]
            for index in range(29)
        ),
        "all_29_reflections_conjugate_exactly": all(
            COMMON_ISOMETRY_A
            * reflection6[index]
            * COMMON_ISOMETRY_A.inv()
            == reflection12[permutation[index]]
            for index in range(29)
        ),
    }
    common_checks = {
        "A_integral": all(value.is_Integer for value in COMMON_ISOMETRY_A),
        "A_unimodular": abs(COMMON_ISOMETRY_A.det()) == 1,
        "A_transports_normalized12_to_maximal6": (
            COMMON_ISOMETRY_A.T * normalized12 * COMMON_ISOMETRY_A
            == gram6_maximal
        ),
        **wall_checks,
    }
    if not all(common_checks.values()):
        raise ArithmeticError("common normalized H11/29 reflection system changed")

    p12 = COMMON_ISOMETRY_A * A16_6_TO_SK5
    sk_checks = {
        "P6_unimodular": abs(A16_6_TO_SK5.det()) == 1,
        "P6_congruence": (
            A16_6_TO_SK5.T * (-gram6_maximal) * A16_6_TO_SK5
            == SK_GRAM11_ENTRY5
        ),
        "P12_unimodular": abs(p12.det()) == 1,
        "P12_congruence": p12.T * (-normalized12) * p12 == SK_GRAM11_ENTRY5,
    }
    if not all(sk_checks.values()):
        raise ArithmeticError("Scharlau--Kirschmer entry-5 congruence changed")

    records["A16_12_1_C01"]["global_half_scaling"] = {
        "operation": "divide the complete Gram globally by 2; not an overlattice",
        "reflection_matrices_unchanged": True,
        "normalized_gram": _payload(normalized12),
        "determinant": int(normalized12.det()),
        "smith_invariants": _smith(normalized12),
        "parity": "odd",
        "discriminant_q_mod_2_histogram": _q_histogram(normalized12),
    }

    other_construction = json.loads(OTHER_CONSTRUCTION.read_text())
    c03 = next(case for case in other_construction["cases"] if case["id"] == "A18_4_0_C03")
    c03_gram = _matrix(c03["saturated_fixed_lattice"]["gram"])
    a16_6_gram = _matrix(cases["A16_6_1_C01"]["saturated_fixed_lattice"]["gram"])
    c03_q = _q_histogram(c03_gram)
    a16_q = _q_histogram(a16_6_gram)
    other_stable = json.loads(OTHER_STABLE.read_text())
    c03_stable = next(case for case in other_stable["cases"] if case["id"] == "A18_4_0_C03")
    collision_checks = {
        "same_rank": c03_gram.rows == a16_6_gram.rows == 12,
        "same_determinant": c03_gram.det() == a16_6_gram.det() == -256,
        "same_smith": _smith(c03_gram) == _smith(a16_6_gram) == [1] * 6 + [2] * 4 + [4] * 2,
        "different_discriminant_quadratic_histograms": c03_q != a16_q,
        "different_maximal_stable_determinants": (
            c03_stable["unique_maximal_reflection_stable_overlattice"]["determinant"]
            == -16
            and gram6_maximal.det() == -4
        ),
    }
    if not all(collision_checks.values()):
        raise ArithmeticError("A18-C03/A16_6 collision audit changed")

    return {
        "schema": SCHEMA,
        "source": {
            "construction_path": str(CONSTRUCTION.relative_to(ROOT)),
            "construction_sha256": _sha256(CONSTRUCTION),
            "collision_construction_path": str(OTHER_CONSTRUCTION.relative_to(ROOT)),
            "collision_construction_sha256": _sha256(OTHER_CONSTRUCTION),
            "collision_stable_path": str(OTHER_STABLE.relative_to(ROOT)),
            "collision_stable_sha256": _sha256(OTHER_STABLE),
        },
        "method": {
            "stable_overlattices": (
                "all totally integral discriminant subgroups invariant under every "
                "facet reflection are enumerated exhaustively"
            ),
            "normalization_warning": (
                "actual integral overclosure and global Gram rescaling are distinct "
                "operations and are recorded separately"
            ),
        },
        "cases": [records["A16_6_1_C01"], records["A16_12_1_C01"]],
        "common_normalized_reflection_system": {
            "A16_6_model": "unique index-8 reflection-stable overclosure",
            "A16_12_model": "global half-scaling of the original fixed Gram",
            "basis_isometry_A": _payload(COMMON_ISOMETRY_A),
            "identity": "A^T (G_12/2) A = M_6",
            "source_to_dual_facet_permutation_one_based": FACET_PERMUTATION,
            "source_root_factors_before_A": FORCED_ROOT_FACTORS,
            "checks": common_checks,
            "conclusion": (
                "the two stocks yield the same integral lattice and the same "
                "29-reflection subgroup after the stated normalizations"
            ),
        },
        "maximal_reflection_lineage": {
            "scharlau_kirschmer_list": "List11",
            "entry_one_based": 5,
            "genus_symbol": "I_{11,1}(2^{2})",
            "P6_to_Gram11_entry5": _payload(A16_6_TO_SK5),
            "P12_to_Gram11_entry5": _payload(p12),
            "checks": sk_checks,
            "full_chamber_facets": 14,
            "full_chamber_finite_vertices": 56,
            "full_chamber_cusps": 5,
            "section_to_full_index_status": (
                "finite; numerical value open because both H11 Coxeter Euler "
                "characteristics vanish"
            ),
            "new_arithmetic_commensurability_class": False,
        },
        "same_header_collision_with_A18_4_0_C03": {
            "A18_C03_discriminant_q_mod_2_histogram": c03_q,
            "A16_6_discriminant_q_mod_2_histogram": a16_q,
            "checks": collision_checks,
            "conclusion": (
                "matching determinant and Smith invariants do not give an integral "
                "isometry; the discriminant quadratic forms and stable overclosures differ"
            ),
        },
        "summary": {
            "stable_integral_subgroup_orders_in_each_stock": [1, 2, 4, 8],
            "stable_even_subgroup_orders_in_each_stock": [1, 4],
            "one_common_normalized_integral_reflection_system": True,
            "known_full_reflection_chamber": "Scharlau--Kirschmer List11 entry 5",
            "new_arithmetic_commensurability_class_count": 0,
            "numerical_full_group_index_open": True,
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text()):
            raise ArithmeticError("rebuilt dual H11/29 stable archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
