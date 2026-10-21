#!/usr/bin/env python3
"""Certify a K3 divisor-overlattice carrier of Alexeev's H11/29 chamber.

The 29-wall chamber is carried by the rank-12 even lattice ``M`` extracted
from Alexeev's (16,6,1) chamber.  Its unique nontrivial reflection-stable
even overclosure has index four.  This script constructs that overclosure
and gives an explicit integral isometry

    N = M[index 4]  ~=  U + D8(-1) + A1(-1)^2.

The latter is the exact very-general Neron--Severi lattice of the archived
Brandhorst--Hofmann 0.4.4.2 Weierstrass family.  All 29 reflections extend
integrally to N.  The result is a mixed divisor-subsystem chamber: its
primitive ambient mirrors have square -2 or -4.  Six normals primitive in
``M`` become twice a primitive square-minus-two class in ``N``.  The square
-4 mirrors are divisor cuts, not asserted curve walls, and extra roots of N
may refine the chamber.  In particular this is not advertised as the full
nef cone or as an automorphism-fixed section.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_118_dual_h11_29_stable_overlattices import (
    A16_6_TO_SK5,
    SK_GRAM11_ENTRY5,
    _active_roots,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _matrix,
    _overlattice_basis,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
)
STABLE = ROOT / "data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json"
VINBERG = ROOT / "data/enumerations/alexeev_118_dual_h11_29_vinberg.json"
FAMILY = ROOT / "data/enumerations/ud8a1x2_order4_weierstrass_family.json"
OUTPUT = (
    ROOT
    / "data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json"
)
SCHEMA = "coxeter.alexeev_h11_29.ud8_divisor_overlattice_carrier.v2"

EXPECTED_SHA256 = {
    "construction": "53432322ecc6210e047c9b6e9ca101a2b2f86cd5cfef758532f0975bd53f9a99",
    "stable": "3083308ff2b269b8e50b5466b0433dde70c41ef8df089c62dff486b42d14b45b",
    "vinberg": "f711b1869f9106ae0c4c6d48673ff8b3de84495dc8ac72a03787e0434e9497d9",
    "family": "8494d89aec70059f090795dc2aa757e64a7ac686e878379cb24a021bc2f9bdc3",
}

# PARI/GP qfisom output for C = S^T (D8+A1^2) S.  It is replayed below,
# so PARI is not required by the permanent certificate.
POSITIVE_COMPLEMENT_ISOMETRY = sp.Matrix(
    [
        [-1, 0, 0, 0, 1, 0, 0, 0, 1, 1],
        [-2, 1, 0, 0, 1, 0, 0, 0, 0, 0],
        [-2, 0, 1, 0, 1, 0, 0, 0, -1, -1],
        [-2, 0, 0, 1, 1, 0, 0, 0, -2, -2],
        [-2, 0, 0, 0, 1, 0, -1, -1, 5, 9],
        [-2, 0, 0, 0, 0, 1, -1, -1, 6, 10],
        [-1, 0, 0, 0, 0, 1, -1, -1, 5, 8],
        [-1, 0, 0, 0, 0, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, -2, -3],
        [0, 0, 0, 0, 0, 0, 1, 0, -2, -3],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(row)] for row in range(matrix.rows)]


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _d10_even_basis() -> sp.Matrix:
    """Columns give the standard D10 simple-root basis in Z^10."""

    basis = sp.zeros(10)
    for index in range(9):
        basis[index, index] = 1
        basis[index + 1, index] = -1
    basis[8, 9] = 1
    basis[9, 9] = 1
    return basis


def _root_reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    norm = (root.T * gram * root)[0]
    return sp.simplify(sp.eye(gram.rows) - sp.Rational(2, norm) * root * root.T * gram)


def _overclosure_and_embedding(
    case: dict[str, Any], family: dict[str, Any]
) -> tuple[dict[str, Any], sp.Matrix, sp.Matrix, list[sp.Matrix]]:
    gram_m = _matrix(case["saturated_fixed_lattice"]["gram"])
    roots_m = _active_roots(case)
    reflections_m = [_reflection(gram_m, root) for root in roots_m]
    audit = _stable_subgroups(gram_m, reflections_m)
    left_inverse = audit["left_smith_inverse"]

    even_four = [
        (position, record)
        for position, record in enumerate(audit["stable_subgroups"])
        if record["even_overlattice"]
        and record["order_and_overlattice_index"] == 4
    ]
    if len(even_four) != 1:
        raise ArithmeticError("the unique stable even index-four overclosure changed")
    _, even_record = even_four[0]
    basis_four = _overlattice_basis(
        gram_m, left_inverse, even_record["cosets_in_smith_coordinates"]
    )
    gram_n_internal = sp.simplify(basis_four.T * gram_m * basis_four)

    maximal_position = audit["maximal_stable_subgroup_positions_zero_based"]
    if len(maximal_position) != 1:
        raise ArithmeticError("the maximal stable overclosure stopped being unique")
    maximal_record = audit["stable_subgroups"][maximal_position[0]]
    basis_eight = _overlattice_basis(
        gram_m, left_inverse, maximal_record["cosets_in_smith_coordinates"]
    )
    gram_maximal = sp.simplify(basis_eight.T * gram_m * basis_eight)

    # The odd index-eight overclosure is the pinned SK5 lattice.  Its even
    # sublattice is D10 plus the final rank-two block.
    d10_basis = _d10_even_basis()
    even_basis = sp.diag(d10_basis, sp.eye(2))
    internal_to_d10 = sp.simplify(
        even_basis.inv()
        * A16_6_TO_SK5.inv()
        * basis_eight.inv()
        * basis_four
    )
    d10_u2_gram = sp.simplify(even_basis.T * SK_GRAM11_ENTRY5 * even_basis)

    # Split a primitive U.  In orthonormal D10 coordinates use
    # d=(1^8,0^2), and in U(-2) use p+2q'.  The boundary root e8-e9 pairs
    # once with d, so v=(e8-e9)-u completes a unimodular U pair.
    d = sp.Matrix([1] * 8 + [0] * 2)
    d_coordinates = d10_basis.inv() * d
    u = sp.Matrix.vstack(d_coordinates, sp.Matrix([-1, 2]))
    boundary_root = sp.zeros(12, 1)
    boundary_root[7] = 1
    v = boundary_root - u
    complement_basis = saturated_row_kernel(
        sp.Matrix.hstack(d10_u2_gram * u, d10_u2_gram * v)
    )
    split_basis = sp.Matrix.hstack(u, v, complement_basis.T)
    complement_gram = sp.simplify(
        complement_basis * d10_u2_gram * complement_basis.T
    )

    family_gram = sp.Matrix(family["neron_severi"]["gram"])
    positive_root_gram = -family_gram[2:, 2:]
    complement_isometry = POSITIVE_COMPLEMENT_ISOMETRY
    split_isometry = sp.diag(sp.eye(2), complement_isometry)
    positive_model = sp.diag(sp.Matrix([[0, 1], [1, 0]]), positive_root_gram)
    internal_to_positive_model = sp.simplify(
        split_isometry * split_basis.inv() * internal_to_d10
    )
    hyperbolic_sign_change = sp.diag(-1, *([1] * 11))
    internal_to_family = hyperbolic_sign_change * internal_to_positive_model

    # Columns of inclusion_m_to_family are the original M basis in the
    # explicit family basis.  Its determinant has absolute value four.
    inclusion_m_to_family = sp.simplify(internal_to_family * basis_four.inv())

    construction_checks = {
        "source_M_rank_signature_determinant": (
            gram_m.rows == 12
            and exact_inertia(gram_m) == (1, 11, 0)
            and gram_m.det() == -256
        ),
        "stable_integral_overclosure_orders": (
            audit["stable_integral_subgroup_orders"] == [1, 2, 4, 8]
        ),
        "stable_even_overclosure_orders": audit["stable_even_subgroup_orders"]
        == [1, 4],
        "unique_nontrivial_stable_even_overclosure": len(even_four) == 1,
        "index_four_basis": abs(basis_four.det()) == sp.Rational(1, 4),
        "overclosure_determinant_square_law": gram_n_internal.det() == -16,
        "overclosure_even": all(
            gram_n_internal[index, index] % 2 == 0 for index in range(12)
        ),
        "overclosure_smith": _smith(gram_n_internal) == [1] * 8 + [2] * 4,
        "overclosure_discriminant_q": _q_histogram(gram_n_internal)
        == {"0": 4, "1": 4, "1/2": 2, "3/2": 6},
        "index_eight_SK5_congruence": (
            A16_6_TO_SK5.T * (-gram_maximal) * A16_6_TO_SK5
            == SK_GRAM11_ENTRY5
        ),
        "D10_even_basis_index_two": abs(even_basis.det()) == 2,
        "internal_to_D10_U2_integral_unimodular": (
            all(value.is_Integer for value in internal_to_d10)
            and abs(internal_to_d10.det()) == 1
            and internal_to_d10.T * d10_u2_gram * internal_to_d10
            == -gram_n_internal
        ),
        "primitive_U_split": (
            (u.T * d10_u2_gram * u)[0] == 0
            and (v.T * d10_u2_gram * v)[0] == 0
            and (u.T * d10_u2_gram * v)[0] == 1
            and abs(split_basis.det()) == 1
        ),
        "positive_complement_is_D8_A1x2": (
            abs(complement_isometry.det()) == 1
            and complement_isometry.T
            * positive_root_gram
            * complement_isometry
            == complement_gram
        ),
        "explicit_internal_overclosure_to_family_isometry": (
            all(value.is_Integer for value in internal_to_family)
            and abs(internal_to_family.det()) == 1
            and internal_to_family.T * family_gram * internal_to_family
            == gram_n_internal
        ),
        "explicit_index_four_M_embedding": (
            all(value.is_Integer for value in inclusion_m_to_family)
            and abs(inclusion_m_to_family.det()) == 4
            and inclusion_m_to_family.T
            * family_gram
            * inclusion_m_to_family
            == gram_m
        ),
    }
    if not all(construction_checks.values()):
        failed = [key for key, value in construction_checks.items() if not value]
        raise ArithmeticError(f"H11/29 overclosure construction failed: {failed}")

    wall_records = []
    for index, root_m in enumerate(roots_m):
        presentation_root_n = inclusion_m_to_family * root_m
        coordinate_content = int(
            sp.gcd_list([abs(int(value)) for value in presentation_root_n])
        )
        primitive_root_n = sp.simplify(presentation_root_n / coordinate_content)
        presentation_norm = int(
            (presentation_root_n.T * family_gram * presentation_root_n)[0]
        )
        primitive_norm = int(
            (primitive_root_n.T * family_gram * primitive_root_n)[0]
        )
        primitive_divisibility = int(
            sp.gcd_list(
                [abs(int(value)) for value in family_gram * primitive_root_n]
            )
        )
        presentation_reflection = _root_reflection(
            family_gram, presentation_root_n
        )
        reflection_n = _root_reflection(family_gram, primitive_root_n)
        intertwining = reflection_n * inclusion_m_to_family == (
            inclusion_m_to_family * reflections_m[index]
        )
        checks = {
            "primitive_in_M": sp.gcd_list([abs(int(x)) for x in root_m]) == 1,
            "NS_coordinate_content_is_one_or_two": coordinate_content in (1, 2),
            "primitive_in_NS": (
                sp.gcd_list([abs(int(value)) for value in primitive_root_n]) == 1
            ),
            "presentation_norm_matches_M": presentation_norm
            == int((root_m.T * gram_m * root_m)[0]),
            "primitive_norm_square_law": (
                presentation_norm == coordinate_content**2 * primitive_norm
            ),
            "same_hyperplane": (
                presentation_root_n == coordinate_content * primitive_root_n
            ),
            "presentation_and_primitive_reflections_equal": (
                presentation_reflection == reflection_n
            ),
            "reflection_integral_on_NS": all(value.is_Integer for value in reflection_n),
            "reflection_unimodular": reflection_n.det() == -1,
            "reflection_involutive": reflection_n**2 == sp.eye(12),
            "reflection_preserves_NS": (
                reflection_n.T * family_gram * reflection_n == family_gram
            ),
            "reflection_extends_M_action": intertwining,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"facet reflection {index} failed to extend")
        wall_records.append(
            {
                "facet_zero_based": index,
                "M_presentation_coordinates": [str(value) for value in root_m],
                "M_presentation_square": presentation_norm,
                "NS_presentation_coordinates": [
                    str(value) for value in presentation_root_n
                ],
                "NS_coordinate_content": coordinate_content,
                "primitive_NS_normal_coordinates": [
                    str(value) for value in primitive_root_n
                ],
                "primitive_NS_square": primitive_norm,
                "primitive_NS_divisibility": primitive_divisibility,
                "checks": checks,
            }
        )

    presentation_square_histogram = dict(
        sorted(
            Counter(
                str(row["M_presentation_square"]) for row in wall_records
            ).items()
        )
    )
    content_histogram = dict(
        sorted(
            Counter(
                str(row["NS_coordinate_content"]) for row in wall_records
            ).items()
        )
    )
    primitive_square_histogram = dict(
        sorted(
            Counter(
                str(row["primitive_NS_square"]) for row in wall_records
            ).items()
        )
    )
    primitive_square_divisibility_histogram = {
        f"{key[0]},{key[1]}": value
        for key, value in sorted(
            Counter(
                (
                    row["primitive_NS_square"],
                    row["primitive_NS_divisibility"],
                )
                for row in wall_records
            ).items()
        )
    }
    normalization_checks = {
        "presentation_square_split_unchanged": presentation_square_histogram
        == {"-2": 11, "-4": 12, "-8": 6},
        "NS_coordinate_content_split": content_histogram == {"1": 23, "2": 6},
        "primitive_NS_square_split": primitive_square_histogram
        == {"-2": 17, "-4": 12},
        "primitive_NS_square_divisibility_split": (
            primitive_square_divisibility_histogram
            == {"-2,1": 16, "-2,2": 1, "-4,2": 12}
        ),
        "no_primitive_square_minus_eight_mirrors": all(
            row["primitive_NS_square"] != -8 for row in wall_records
        ),
    }
    if not all(normalization_checks.values()):
        raise ArithmeticError("the primitive ambient-normal split changed")

    return (
        {
            "original_divisor_lattice_M": {
                "rank": 12,
                "signature": [1, 11, 0],
                "determinant": -256,
                "smith": _smith(gram_m),
                "discriminant_q_mod_2_histogram": _q_histogram(gram_m),
                "gram": _payload(gram_m),
            },
            "reflection_stable_even_overclosure_census": {
                "stable_integral_orders": audit["stable_integral_subgroup_orders"],
                "stable_even_orders": audit["stable_even_subgroup_orders"],
                "nontrivial_even_index_four_count": len(even_four),
                "unique": True,
                "cosets_in_smith_coordinates": even_record[
                    "cosets_in_smith_coordinates"
                ],
            },
            "overclosure_N": {
                "isometry_class": "U+D8(-1)+A1(-1)^2",
                "rank": 12,
                "signature": [1, 11, 0],
                "determinant": -16,
                "smith": _smith(gram_n_internal),
                "discriminant_q_mod_2_histogram": _q_histogram(gram_n_internal),
                "basis_columns_in_M_tensor_Q": _payload(basis_four),
                "internal_gram": _payload(gram_n_internal),
                "explicit_internal_to_family_isometry": _payload(
                    internal_to_family
                ),
                "positive_D8_A1x2_complement_isometry": _payload(
                    complement_isometry
                ),
                "checks": construction_checks,
            },
            "explicit_M_into_family_NS": {
                "index": 4,
                "basis_columns": _payload(inclusion_m_to_family),
                "basis_columns_sha256": _payload_sha256(
                    _payload(inclusion_m_to_family)
                ),
                "identity": "I^T Gram(NS) I = Gram(M)",
            },
            "facet_reflection_extensions": {
                "count": len(wall_records),
                "M_presentation_square_histogram": presentation_square_histogram,
                "NS_coordinate_content_histogram": content_histogram,
                "primitive_NS_square_histogram": primitive_square_histogram,
                "primitive_NS_square_divisibility_histogram": (
                    primitive_square_divisibility_histogram
                ),
                "normalization_checks": normalization_checks,
                "all_29_extend_integrally": True,
                "records": wall_records,
            },
        },
        gram_m,
        inclusion_m_to_family,
        roots_m,
    )


def build_certificate() -> dict[str, Any]:
    for key, path in (
        ("construction", CONSTRUCTION),
        ("stable", STABLE),
        ("vinberg", VINBERG),
        ("family", FAMILY),
    ):
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    stable = json.loads(STABLE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    case = next(row for row in construction["cases"] if row["id"] == "A16_6_1_C01")

    lattice, gram_m, inclusion, roots_m = _overclosure_and_embedding(case, family)

    # The explicit order-four family has Picard involution equal to the A1
    # swap.  Under our embedding M is stable, although its 29-simple-wall set
    # is not stable; the latter warning prevents a fixed-section overclaim.
    family_gram = sp.Matrix(family["neron_severi"]["gram"])
    picard_involution = sp.eye(12)
    picard_involution[10, 10] = picard_involution[11, 11] = 0
    picard_involution[10, 11] = picard_involution[11, 10] = 1
    induced_m_action = sp.simplify(inclusion.inv() * picard_involution * inclusion)
    primitive_wall_keys = {
        tuple(int(value) for value in root): index for index, root in enumerate(roots_m)
    }
    preserved_oriented_wall_count = sum(
        tuple(int(value) for value in induced_m_action * root) in primitive_wall_keys
        for root in roots_m
    )
    action_checks = {
        "family_NS_action_is_A1_swap": (
            picard_involution.T * family_gram * picard_involution == family_gram
            and picard_involution**2 == sp.eye(12)
        ),
        "M_is_stable_under_family_Picard_involution": all(
            value.is_Integer for value in induced_m_action
        ),
        "induced_M_action_is_involutive_isometry": (
            induced_m_action**2 == sp.eye(12)
            and induced_m_action.T * gram_m * induced_m_action == gram_m
        ),
        "selected_oriented_simple_wall_set_not_preserved": (
            preserved_oriented_wall_count < 29
        ),
    }
    if not all(action_checks.values()):
        raise ArithmeticError("the scoped family-action handoff changed")

    family_checks = {
        "family_certificate_green": family["all_checks_pass"],
        "very_general_NS_exact": (
            family["neron_severi"]["generic_exactness"]["conclusion"]
            == "very general NS=U+D8(-1)+A1(-1)^2"
        ),
        "NS_gram_matches_explicit_isometry_target": (
            family["neron_severi"]["rank_signature_determinant"]
            == [12, [1, 11, 0], -16]
        ),
        "global_order_four": family["family"]["global_order"] == 4,
        "BH_row_0_4_4_2": family["brandhorst_hofmann_handoff"]["id"]
        == "0.4.4.2",
        "root_stock": family["neron_severi"]["reducible_fiber_roots"]
        == "D8(-1)+A1(-1)^2",
    }
    if not all(family_checks.values()):
        raise ArithmeticError("the explicit K3 family handoff changed")

    vinberg_case = next(row for row in vinberg["audits"] if row["id"] == "A16_6_1_C01")
    vinberg_checks = {
        "literal_finite_volume_noncompact": (
            vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "dimension_facets": (
            vinberg_case["dimension"] == 11 and vinberg_case["facet_count"] == 29
        ),
        "every_elliptic_edge_has_two_endpoints": vinberg_case["checks"][
            "every_elliptic_edge_has_exactly_two_endpoints"
        ],
        "all_affine_neighborhoods_full_rank": vinberg_case["checks"][
            "all_affine_zero_neighborhoods_are_full_rank_parabolic"
        ],
        "zero_bad_audits": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["bad_parabolic_neighborhood_count"] == 0
        ),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("the literal Vinberg handoff changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "The finite-volume noncompact Alexeev H11/29 chamber has a genuine "
            "equation-level K3 divisor carrier: its rank-12 lattice is a "
            "reflection-stable index-four divisor sublattice of the very-general "
            "U+D8(-1)+A1(-1)^2 Neron--Severi lattice in the BH 0.4.4.2 family."
        ),
        "scope": {
            "proved": (
                "full-rank mixed divisor-subsystem chamber with all 29 facet "
                "reflections integral on the ambient Neron--Severi lattice"
            ),
            "not_full_nef": (
                "extra -2 roots of the ambient index-four overclosure may refine "
                "the chamber"
            ),
            "not_all_curve_walls": (
                "the twelve primitive square -4 mirrors are algebraic divisor cuts and are "
                "not asserted to be irreducible curves"
            ),
            "not_automorphism_fixed_section": (
                "the BH Picard involution preserves the divisor sublattice but "
                "not the selected 29-simple-wall set under the explicit embedding"
            ),
            "abstract_novelty": (
                "no new abstract polytope claim; this closes a previously open "
                "algebro-geometric realization of the known project H11/29 row"
            ),
        },
        "sources": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in (
                ("construction", CONSTRUCTION),
                ("stable", STABLE),
                ("vinberg", VINBERG),
                ("family", FAMILY),
            )
        },
        "bounded_primary_source_lineage_audit": {
            "alexeev": {
                "url": "https://arxiv.org/abs/2209.09110",
                "location": "Section 3",
                "records": "the 118-root (16,6,1) chamber and its S6 symmetry",
                "does_not_record": (
                    "the transposition section or its even index-four overclosure"
                ),
            },
            "brandhorst_hofmann": {
                "article": "https://arxiv.org/abs/2112.07715",
                "database_v2_1": "https://zenodo.org/records/10410278",
                "records": "classification row 0.4.4.2",
                "does_not_record": (
                    "this 29-mirror chamber, its root list, or the displayed equation"
                ),
            },
            "roulleau": {
                "url": "https://doi.org/10.46298/epiga.2022.6286",
                "location": "Section 12.2",
                "records": (
                    "U+D8(-1)+A1(-1)^2 through its 19 curve walls"
                ),
                "does_not_record": "this mixed 29-mirror subsystem",
            },
            "conclusion": (
                "the equation-to-overlattice-to-chamber realization appears "
                "unrecorded in these audited primary sources"
            ),
            "scope": "bounded primary-source evidence, not worldwide novelty",
        },
        "lattice_and_embedding": lattice,
        "explicit_K3_family": {
            "equation": family["family"]["affine_equation"],
            "coefficient_polynomials": {
                "P": family["family"]["P"],
                "Q": family["family"]["Q"],
            },
            "parameter_open": family["family"]["parameter_open"],
            "automorphism": family["family"]["automorphism"],
            "very_general_NS": "U+D8(-1)+A1(-1)^2",
            "BH_id": "0.4.4.2",
            "checks": family_checks,
        },
        "family_action_on_divisor_sublattice": {
            "ambient_picard_action": "identity on U+D8(-1), swap the two A1 roots",
            "induced_M_action": _payload(induced_m_action),
            "selected_oriented_simple_walls_preserved": (
                preserved_oriented_wall_count
            ),
            "selected_oriented_simple_wall_count": 29,
            "conclusion": (
                "the selected Coxeter chamber is not preserved, so this is not "
                "an automorphism-fixed-section construction"
            ),
            "checks": action_checks,
        },
        "vinberg_finite_volume": {
            "classification": vinberg_case["classification"],
            "dimension": 11,
            "facets": 29,
            "finite_vertices": vinberg_case["finite_vertex_count"],
            "cusps": vinberg_case["ideal_vertex_count"],
            "elliptic_edges": vinberg_case["elliptic_edge_count"],
            "affine_component_audits": vinberg_case[
                "connected_affine_component_count"
            ],
            "bad_edges": vinberg_case["bad_edge_count"],
            "bad_parabolic_neighborhoods": vinberg_case[
                "bad_parabolic_neighborhood_count"
            ],
            "criterion": vinberg_case["criterion"],
            "checks": vinberg_checks,
        },
        "effectivity_and_divisor_interpretation": {
            "strict_positive_chamber_interior": case["strict_fixed_interior"],
            "generic_ample_choice": (
                "choose a generic point in the strict chamber away from every "
                "ambient -2 hyperplane, then Weyl re-mark the K3 so it is ample"
            ),
            "minus_two_wall_effectivity": (
                "K3 Riemann--Roch chooses the effective sign of each of the seventeen "
                "square -2 wall classes against that ample point"
            ),
            "mixed_mirror_warning": (
                "the remaining twelve square -4 mirrors are divisor reflections, not "
                "claimed effective irreducible curve walls"
            ),
        },
        "all_checks_pass": True,
    }
    payload["certificate_payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H11/29 K3 carrier archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
