#!/usr/bin/env python3
"""Certify an equation-level UE8D4 carrier of the A18 C06 H13/31 chamber.

The archived rank-14 wall lattice M has a unique nontrivial even
reflection-stable overclosure of index four.  The explicit matrix below
embeds M with index four in the very-general Neron--Severi lattice

    U + E8(-1) + D4(-1)

of the archived order-six Weierstrass family.  All 31 facet reflections
extend integrally.  Primitive reduction is essential: twelve square -8
presentation normals become primitive square -2 roots.  The resulting
chamber is a mixed divisor-reflection subsystem, not a full-nef claim and
not an automorphism-fixed-section claim on this new carrier.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _matrix,
    _overlattice_basis,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
STABLE = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
VINBERG = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
)
FAMILY = ROOT / "data/enumerations/ue8d4_order6_weierstrass_family.json"
OUTPUT = (
    ROOT / "data/enumerations/a18_c06_ue8d4_divisor_overlattice_carrier.json"
)
SCHEMA = "coxeter.a18_c06.ue8d4_divisor_overlattice_carrier.v1"

EXPECTED_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "stable": "882204e5457e96a3930399201dd061bd6bc35839497693a384d0b931f10924a0",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "family": "0ce8340494520835daf4df1aaaa57239640759111fb45682dc9b6b330b618fff",
}

# Columns are the archived M basis in the explicit UE8D4 family basis.
# This was obtained by splitting a primitive U in the unique even index-four
# overclosure and applying an exact positive-definite qfisom calculation.
# The permanent certificate replays the resulting congruence directly and
# therefore does not depend on PARI/GP.
INCLUSION_M_TO_FAMILY = sp.Matrix(
    [
        [6, 12, -7, -5, 7, 1, -1, 18, -11, -6, -4, 4, -14, -10],
        [1, 2, -1, -1, 1, 0, 0, 3, -2, -1, 0, 0, -1, -2],
        [5, 9, -5, -3, 5, 0, 0, 13, -9, -4, -2, 0, -8, -7],
        [9, 18, -10, -6, 10, 0, 0, 27, -18, -9, -4, 0, -16, -14],
        [13, 27, -15, -9, 14, 0, 0, 41, -27, -13, -6, 0, -25, -22],
        [11, 23, -13, -8, 12, 1, -1, 35, -23, -11, -6, 0, -21, -19],
        [9, 19, -10, -7, 10, 1, -1, 28, -18, -9, -4, 0, -16, -15],
        [6, 13, -7, -5, 7, 0, -1, 19, -12, -6, -2, 0, -11, -10],
        [3, 7, -4, -3, 4, 0, 0, 10, -6, -3, -2, 0, -6, -6],
        [6, 13, -7, -4, 7, 0, 0, 20, -13, -6, -2, 0, -12, -11],
        [-1, -2, 1, 2, -2, 0, 0, -4, 3, 1, 0, -1, 2, 3],
        [-3, -4, 2, 3, -3, 0, 0, -8, 5, 3, 0, -2, 4, 6],
        [-2, -3, 2, 2, -2, 0, 0, -5, 3, 2, 0, -1, 3, 4],
        [-2, -4, 2, 2, -2, 0, 0, -6, 4, 2, 0, -2, 3, 4],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, col]) for col in range(matrix.cols)] for row in range(matrix.rows)]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _root_reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _active_roots(case: dict[str, Any]) -> list[sp.Matrix]:
    return [
        sp.Matrix(row["primitive_integral_same_direction_coordinates"])
        for row in case["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    ]


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "stable": STABLE,
        "vinberg": VINBERG,
        "family": FAMILY,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    stable = json.loads(STABLE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    case = next(row for row in construction["cases"] if row["id"] == "A18_4_0_C06")
    stable_case = next(row for row in stable["cases"] if row["id"] == "A18_4_0_C06")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "A18_4_0_C06")

    gram_m = _matrix(case["saturated_fixed_lattice"]["gram"])
    roots_m = _active_roots(case)
    reflections_m = [_reflection(gram_m, root) for root in roots_m]
    audit = _stable_subgroups(gram_m, reflections_m)
    even_four = [
        row
        for row in audit["stable_subgroups"]
        if row["even_overlattice"] and row["order_and_overlattice_index"] == 4
    ]
    if len(even_four) != 1:
        raise ArithmeticError("the unique even index-four overclosure changed")
    basis_n = _overlattice_basis(
        gram_m,
        audit["left_smith_inverse"],
        even_four[0]["cosets_in_smith_coordinates"],
    )
    gram_n_internal = sp.simplify(basis_n.T * gram_m * basis_n)
    family_gram = sp.Matrix(family["neron_severi"]["gram"])
    inclusion = INCLUSION_M_TO_FAMILY
    internal_to_family = sp.simplify(inclusion * basis_n)

    lattice_checks = {
        "source_rank_signature_determinant": (
            gram_m.rows == 14
            and exact_inertia(gram_m) == (1, 13, 0)
            and gram_m.det() == -64
        ),
        "stable_integral_orders": audit["stable_integral_subgroup_orders"] == [1, 4],
        "stable_even_orders": audit["stable_even_subgroup_orders"] == [1, 4],
        "unique_nontrivial_even_overclosure": len(even_four) == 1,
        "index_four_basis": abs(basis_n.det()) == sp.Rational(1, 4),
        "overclosure_even_determinant_four": (
            gram_n_internal.det() == -4
            and all(gram_n_internal[i, i] % 2 == 0 for i in range(14))
        ),
        "overclosure_smith": _smith(gram_n_internal) == [1] * 12 + [2, 2],
        "overclosure_q_histogram": _q_histogram(gram_n_internal) == {"0": 1, "1": 3},
        "archived_stable_census_agrees": (
            stable_case["exhaustive_discriminant_module_audit"][
                "stable_even_subgroup_orders"
            ]
            == [1, 4]
            and stable_case["unique_maximal_reflection_stable_overlattice"][
                "determinant"
            ]
            == -4
        ),
        "explicit_index_four_embedding": (
            all(value.is_Integer for value in inclusion)
            and abs(inclusion.det()) == 4
            and inclusion.T * family_gram * inclusion == gram_m
        ),
        "explicit_overclosure_isometry": (
            all(value.is_Integer for value in internal_to_family)
            and abs(internal_to_family.det()) == 1
            and internal_to_family.T * family_gram * internal_to_family
            == gram_n_internal
        ),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            f"lattice bridge failed: {[k for k, v in lattice_checks.items() if not v]}"
        )

    wall_records: list[dict[str, Any]] = []
    for index, (source_record, root_m, reflection_m) in enumerate(
        zip(
            [
                row
                for row in case["restricted_mirrors"]["records"]
                if row["hyperbolic_status"] == "active_facet"
            ],
            roots_m,
            reflections_m,
            strict=True,
        )
    ):
        presentation = inclusion * root_m
        content = math.gcd(*(abs(int(value)) for value in presentation))
        primitive = sp.simplify(presentation / content)
        presentation_square = int((presentation.T * family_gram * presentation)[0])
        primitive_square = int((primitive.T * family_gram * primitive)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in family_gram * primitive))
        reflection_n = _root_reflection(family_gram, primitive)
        presentation_reflection = _root_reflection(family_gram, presentation)
        scale = sp.Rational(source_record["presentation_scale_from_primitive"])
        oriented_primitive = sp.sign(scale) * primitive
        checks = {
            "ambient_primitive": math.gcd(*(abs(int(value)) for value in primitive)) == 1,
            "square_law": presentation_square == content**2 * primitive_square,
            "same_hyperplane": presentation == content * primitive,
            "same_reflection_after_reduction": reflection_n == presentation_reflection,
            "reflection_integral": all(value.is_Integer for value in reflection_n),
            "reflection_unimodular_involution": (
                reflection_n.det() == -1 and reflection_n**2 == sp.eye(14)
            ),
            "reflection_preserves_family_NS": (
                reflection_n.T * family_gram * reflection_n == family_gram
            ),
            "reflection_intertwines_M": reflection_n * inclusion
            == inclusion * reflection_m,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"wall {index} failed")
        wall_records.append(
            {
                "facet_zero_based": index,
                "source_orbit_index": source_record["orbit_index"],
                "M_coordinates": [str(value) for value in root_m],
                "M_presentation_square": presentation_square,
                "NS_presentation_coordinates": [str(value) for value in presentation],
                "NS_coordinate_content": content,
                "primitive_NS_coordinates": [str(value) for value in primitive],
                "oriented_primitive_NS_coordinates": [
                    str(value) for value in oriented_primitive
                ],
                "primitive_NS_square": primitive_square,
                "primitive_NS_divisibility": divisibility,
                "checks": checks,
            }
        )

    presentation_histogram = dict(
        sorted(Counter(str(row["M_presentation_square"]) for row in wall_records).items())
    )
    content_histogram = dict(
        sorted(Counter(str(row["NS_coordinate_content"]) for row in wall_records).items())
    )
    primitive_histogram = dict(
        sorted(Counter(str(row["primitive_NS_square"]) for row in wall_records).items())
    )
    primitive_type_histogram = {
        f"{square},{divisibility}": count
        for (square, divisibility), count in sorted(
            Counter(
                (row["primitive_NS_square"], row["primitive_NS_divisibility"])
                for row in wall_records
            ).items()
        )
    }
    normalization_checks = {
        "presentation_split": presentation_histogram == {"-2": 14, "-4": 5, "-8": 12},
        "content_split": content_histogram == {"1": 19, "2": 12},
        "primitive_split": primitive_histogram == {"-2": 26, "-4": 5},
        "primitive_type_split": primitive_type_histogram == {"-2,1": 26, "-4,2": 5},
        "no_primitive_square_minus_eight": all(
            row["primitive_NS_square"] != -8 for row in wall_records
        ),
    }
    if not all(normalization_checks.values()):
        raise ArithmeticError("primitive mirror normalization changed")

    interior_m = sp.Matrix(
        [sp.Rational(value) for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    interior_ns = inclusion * interior_m
    oriented_pairings = [
        (interior_ns.T * family_gram * sp.Matrix(row["oriented_primitive_NS_coordinates"]))[0]
        for row in wall_records
    ]
    effectivity_checks = {
        "strict_interior_positive": (interior_ns.T * family_gram * interior_ns)[0] == 385,
        "all_oriented_pairings_positive": min(oriented_pairings) == 1,
        "minus_two_count": sum(row["primitive_NS_square"] == -2 for row in wall_records)
        == 26,
        "mixed_minus_four_count": sum(row["primitive_NS_square"] == -4 for row in wall_records)
        == 5,
    }
    if not all(effectivity_checks.values()):
        raise ArithmeticError("effectivity/orientation handoff changed")

    family_checks = {
        "family_certificate_green": family["all_checks_pass"] is True,
        "very_general_NS_exact": (
            family["neron_severi"]["lower_bound_lattice"]
            == "U+E8(-1)+D4(-1)"
            and family["neron_severi"]["generic_exactness"][
                "parameter_moduli_dimension"
            ]
            == 3
            and "Hence very general rho=14."
            in family["neron_severi"]["generic_exactness"]["argument"]
        ),
        "NS_gram_exact_target": family_gram == sp.Matrix(family["neron_severi"]["gram"]),
        "equation_recorded": "y^2 = x^3" in family["family"]["affine_equation"],
        "order_six_action_recorded": family["family"]["global_order"] == 6,
    }
    if not all(family_checks.values()):
        raise ArithmeticError("explicit UE8D4 family handoff changed")

    vinberg_checks = {
        "classification": vinberg_case["classification"]
        == "finite_volume_noncompact_Coxeter_polytope",
        "dimension_facets": (vinberg_case["dimension"], vinberg_case["facet_count"])
        == (13, 31),
        "counts": (
            vinberg_case["finite_vertex_count"],
            vinberg_case["ideal_vertex_count"],
            vinberg_case["elliptic_edge_count"],
            vinberg_case["connected_affine_component_count"],
        )
        == (3433, 238, 28766, 529),
        "zero_bad_audits": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["bad_parabolic_neighborhood_count"] == 0
        ),
        "literal_checks_green": all(vinberg_case["checks"].values()),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("literal Vinberg handoff changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "The finite-volume noncompact A18 C06 H13/31 chamber has an "
            "equation-level K3 divisor carrier: its wall lattice embeds with "
            "index four in the very-general U+E8(-1)+D4(-1) Picard lattice."
        ),
        "scope": {
            "proved": "full-rank mixed divisor-reflection subsystem on an explicit K3 family",
            "not_full_nef": "five primitive square -4 mirrors are not K3 curve walls, and extra Picard roots may refine the chamber",
            "not_automorphism_fixed": "no claim is made that the order-six family automorphism preserves this 31-wall chamber",
            "novelty": "new equation-level carrier of an already certified chamber; no new abstract-polytope or arithmetic-class claim",
        },
        "sources": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        },
        "lattice_and_embedding": {
            "wall_lattice_M": {
                "rank": 14,
                "signature": [1, 13, 0],
                "determinant": -64,
                "smith": _smith(gram_m),
                "gram": _payload(gram_m),
            },
            "stable_even_overclosure_census": {
                "stable_integral_orders": audit["stable_integral_subgroup_orders"],
                "stable_even_orders": audit["stable_even_subgroup_orders"],
                "unique_nontrivial_even_index_four": True,
                "cosets_in_smith_coordinates": even_four[0]["cosets_in_smith_coordinates"],
            },
            "overclosure_N": {
                "isometry_class": "U+E8(-1)+D4(-1)",
                "rank": 14,
                "signature": [1, 13, 0],
                "determinant": -4,
                "smith": _smith(gram_n_internal),
                "discriminant_q_mod_2_histogram": _q_histogram(gram_n_internal),
                "basis_columns_in_M_tensor_Q": _payload(basis_n),
                "internal_to_family_isometry": _payload(internal_to_family),
            },
            "explicit_M_into_family_NS": {
                "index": 4,
                "basis_columns": _payload(inclusion),
                "basis_columns_sha256": _payload_sha256(_payload(inclusion)),
                "identity": "I^T Gram(NS) I = Gram(M)",
            },
            "checks": lattice_checks,
        },
        "primitive_mirrors": {
            "count": len(wall_records),
            "M_presentation_square_histogram": presentation_histogram,
            "NS_coordinate_content_histogram": content_histogram,
            "primitive_NS_square_histogram": primitive_histogram,
            "primitive_NS_square_divisibility_histogram": primitive_type_histogram,
            "normalization_checks": normalization_checks,
            "all_31_reflections_extend_integrally": True,
            "records": wall_records,
        },
        "explicit_K3_family": {
            "equation": family["family"]["affine_equation"],
            "parameter_open": family["family"]["parameter_open"],
            "very_general_NS": "U+E8(-1)+D4(-1)",
            "automorphism": family["family"]["automorphism"],
            "checks": family_checks,
        },
        "effectivity_and_scope": {
            "interior_in_family_coordinates": [str(value) for value in interior_ns],
            "minimum_oriented_pairing": str(min(oriented_pairings)),
            "minus_two_effectivity": (
                "after Weyl transport of the marking, K3 Riemann--Roch chooses "
                "the effective signs of all 26 primitive square -2 classes"
            ),
            "mixed_warning": (
                "the five primitive square -4 divisibility-two mirrors are "
                "algebraic divisor reflections, not asserted irreducible curves"
            ),
            "checks": effectivity_checks,
        },
        "vinberg_finite_volume": {
            "classification": vinberg_case["classification"],
            "dimension": 13,
            "facets": 31,
            "finite_vertices": 3433,
            "cusps": 238,
            "elliptic_edges": 28766,
            "connected_affine_components": 529,
            "bad_edges": 0,
            "bad_parabolic_neighborhoods": 0,
            "checks": vinberg_checks,
        },
        "all_checks_pass": True,
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("rebuilt C06 UE8D4 carrier archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
