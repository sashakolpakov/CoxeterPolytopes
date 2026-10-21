#!/usr/bin/env python3
"""Certify nef framing and mirror effectivity for the complex A16 H11/29 carrier.

The existing carrier embeds the A16_6 rank-12 divisor lattice with index four
in the very-general Picard lattice ``U+D8(-1)+A1(-1)^2`` of the explicit
Brandhorst--Hofmann 0.4.4.2 family.  Its coordinates use the standard
``U+D8+A1^2`` basis.  This certificate first applies an exact unimodular
change to Roulleau's basis ``A1,...,A12``, fixes the positive-cone sign using
the published ample class D136, and then Weyl-reduces one strict carrier
interior with the complete stock of nineteen irreducible (-2)-curves.

All twenty-nine mirrors are transported simultaneously.  The seventeen
square-minus-two orientations are decomposed exactly in the complete curve
stock (one is itself A16 and sixteen are reducible); the twelve square-minus-
four divisibility-two vectors are retained separately.  Finally, H+D136 is
checked to be a primitive ample class strict on every imported wall.  None of
this promotes the imported 29-wall chamber to the full nef chamber.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CARRIER = (
    ROOT / "data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json"
)
ROULLEAU = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
INVENTORY = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
FAMILY = ROOT / "data/enumerations/ud8a1x2_order4_weierstrass_family.json"
OUTPUT = ROOT / "data/enumerations/a16_h11_29_ud8a1x2_nef_effectivity.json"

SCHEMA = "coxeter.a16_h11_29.ud8a1x2_nef_effectivity.v1"
EXPECTED_SHA256 = {
    "carrier": "5d4c2299ecf2b8b853063c142472c46c61f8554197168b06d144924cdc6d07fd",
    "roulleau": "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    "inventory": "c3e9c0fba30c76c3e90bdb34fcf1d7d8b25574d1cd67ce1937b88e6905dbb317",
    "family": "8494d89aec70059f090795dc2aa757e64a7ac686e878379cb24a021bc2f9bdc3",
}
EXPECTED_CARRIER_PAYLOAD_SHA256 = (
    "24a73754ddfa9f1254f4a6f7a83575dedc19f50e54e1f96d67576b5c52d8634a"
)
EXPECTED_ROULLEAU_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)

# Rows are a standard U+D8+A1^2 basis written in Roulleau's A1,...,A12
# basis.  Thus its transpose sends standard column coordinates to Roulleau
# column coordinates.  The permanent replay verifies B*G_R*B^t=G_standard.
STANDARD_BASIS_ROWS_IN_ROULLEAU = sp.Matrix(
    [
        [1, 3, 0, 0, -3, -6, -9, -12, -3, -6, -6, -1],
        [0, 1, 0, 0, -1, -2, -3, -4, -1, -2, -2, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 3, 0, 0, -3, -6, -9, -12, -3, -6, -6, -1],
        [2, 3, 2, 1, 0, -2, -4, -6, -2, -3, -3, 1],
    ]
)

# Roulleau's published ample divisor in the A1,...,A12 basis.
D136 = sp.Matrix([2, 5, 9, 4, 10, 12, 15, 19, 9, 9, 6, 0])

# Deterministic first-negative-wall descent.  The literal word is frozen so
# the certificate catches changes in either the marking or wall order.
WEYL_WORD_ONE_BASED = [
    1, 4, 3, 5, 3, 6, 5, 3, 7, 8, 7, 11, 8, 7, 6, 5, 3, 4, 10,
    8, 7, 6, 5, 11, 8, 7, 12, 1, 13, 1, 12, 15, 18, 13, 1, 2, 18,
    13, 1, 2, 3, 5, 6, 7, 8, 11, 12, 1, 2, 3, 5, 6, 7, 8, 13, 1,
    2, 3, 5, 6, 7, 12, 1, 2, 3, 5, 6, 13, 1, 2, 3, 5, 12, 1, 2,
    3, 4, 13, 1, 2, 3, 5, 6, 12, 1, 2, 13,
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("nonintegral matrix entered an integral certificate")
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _vector_payload(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1 or not all(value.is_Integer for value in vector):
        raise ArithmeticError("nonintegral column entered an integral certificate")
    return [int(value) for value in vector]


def _dot(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> int:
    value = sp.Rational((left.T * gram * right)[0])
    if value.q != 1:
        raise ArithmeticError("nonintegral pairing entered the certificate")
    return int(value)


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _root_reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _dot(gram, root, root)
    return sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram


def _support_is_connected(coefficients: list[int], wall_gram: sp.Matrix) -> bool:
    support = [index for index, coefficient in enumerate(coefficients) if coefficient]
    if not support:
        return False
    reached = {support[0]}
    frontier = [support[0]]
    while frontier:
        vertex = frontier.pop()
        for neighbor in support:
            if neighbor not in reached and wall_gram[vertex, neighbor] > 0:
                reached.add(neighbor)
                frontier.append(neighbor)
    return reached == set(support)


def _decompose_minus_two_root(
    root: sp.Matrix,
    gram: sp.Matrix,
    curves: list[sp.Matrix],
    ample: sp.Matrix,
) -> tuple[list[int], list[dict[str, int | str]], int]:
    """Descend an effective root to one member of the complete curve stock."""

    coefficients = [0] * len(curves)
    current = root
    trace: list[dict[str, int | str]] = []
    seen: set[tuple[int, ...]] = set()
    for _ in range(1000):
        key = tuple(_vector_payload(current))
        if key in seen:
            raise ArithmeticError("minus-two descent cycled")
        seen.add(key)

        terminal = next(
            (index for index, curve in enumerate(curves) if current == curve), None
        )
        if terminal is not None:
            coefficients[terminal] += 1
            return coefficients, trace, terminal + 1

        pairings = [_dot(gram, current, curve) for curve in curves]
        negative = [index for index, pairing in enumerate(pairings) if pairing < 0]
        if not negative:
            raise ArithmeticError("root descent missed the complete curve stock")
        index = negative[0]
        pairing = pairings[index]
        multiplicity = -pairing
        height_before = _dot(gram, current, ample)
        following = current - multiplicity * curves[index]
        height_after = _dot(gram, following, ample)
        if not (
            multiplicity > 0
            and _dot(gram, following, following) == -2
            and 0 < height_after < height_before
        ):
            raise ArithmeticError("invalid strict ample-height root descent")
        coefficients[index] += multiplicity
        trace.append(
            {
                "curve": f"A{index + 1}",
                "curve_one_based": index + 1,
                "negative_pairing": pairing,
                "multiplicity_added": multiplicity,
                "D136_degree_before": height_before,
                "D136_degree_after": height_after,
            }
        )
        current = following
    raise ArithmeticError("root descent exceeded its hard finite bound")


def build_certificate() -> dict[str, Any]:
    paths = {
        "carrier": CARRIER,
        "roulleau": ROULLEAU,
        "inventory": INVENTORY,
        "family": FAMILY,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    roulleau = json.loads(ROULLEAU.read_text(encoding="utf-8"))
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    stock = next(
        row for row in roulleau["stocks"] if row["stock_key"] == "R12_UD8A1x2"
    )
    source_stock = next(
        row for row in inventory["stocks"] if row["key"] == "R12_UD8A1x2"
    )

    carrier_without_digest = dict(carrier)
    carrier_digest = carrier_without_digest.pop("certificate_payload_sha256")
    source_checks = {
        "all_four_input_archives_are_green": (
            carrier["all_checks_pass"]
            and roulleau["all_checks_pass"]
            and inventory["all_checks_pass"]
            and family["all_checks_pass"]
        ),
        "carrier_payload_digest_replays": (
            carrier_digest == EXPECTED_CARRIER_PAYLOAD_SHA256
            and _payload_sha256(carrier_without_digest) == carrier_digest
        ),
        "correct_finite_volume_carrier": (
            carrier["vinberg_finite_volume"]["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and carrier["vinberg_finite_volume"]["dimension"] == 11
            and carrier["vinberg_finite_volume"]["facets"] == 29
        ),
        "correct_exact_family": (
            family["neron_severi"]["generic_exactness"]["conclusion"]
            == "very general NS=U+D8(-1)+A1(-1)^2"
            and family["brandhorst_hofmann_handoff"]["id"] == "0.4.4.2"
        ),
        "published_19_curve_stock_is_complete": (
            stock["wall_count"] == len(stock["wall_coordinate_rows"]) == 19
            and source_stock["wall_count"] == 19
            and source_stock["source_completeness_method"]
            == "direct_curve_enumeration_and_completeness_test"
            and source_stock["reconstruction_status"] == "already_exactly_archived"
        ),
        "primary_source_is_pinned": (
            roulleau["source_tex_sha256"] == EXPECTED_ROULLEAU_TEX_SHA256
            and inventory["source"]["sha256"]["tex"]
            == EXPECTED_ROULLEAU_TEX_SHA256
            and source_stock["source"]["section_lines"] == [8648, 8803]
        ),
    }
    if not all(source_checks.values()):
        raise ArithmeticError("source handoff for H11/29 nef framing failed")

    gram = sp.Matrix(stock["basis_gram"])
    wall_gram = sp.Matrix(stock["wall_gram"])
    curves = [sp.Matrix(row) for row in stock["wall_coordinate_rows"]]
    curve_matrix = sp.Matrix.hstack(*curves)
    standard_gram = sp.Matrix(family["neron_severi"]["gram"])
    basis_rows = STANDARD_BASIS_ROWS_IN_ROULLEAU
    standard_to_roulleau = basis_rows.T
    ample_degrees = [_dot(gram, D136, curve) for curve in curves]
    expected_ample_degrees = [1] * 10 + [7, 2, 2, 12, 18, 18, 12, 18, 18]
    stock_checks = {
        "A1_through_A12_are_the_NS_basis": curve_matrix[:, :12] == sp.eye(12),
        "all_19_curves_have_square_minus_two": all(
            _dot(gram, curve, curve) == -2 for curve in curves
        ),
        "archived_complete_wall_gram_recomputes": (
            curve_matrix.T * gram * curve_matrix == wall_gram
        ),
        "D136_has_published_square": _dot(gram, D136, D136) == 136,
        "D136_has_published_positive_degrees": (
            ample_degrees == expected_ample_degrees and min(ample_degrees) == 1
        ),
        "standard_to_Roulleau_bridge_is_unimodular": basis_rows.det() == -1,
        "standard_to_Roulleau_bridge_is_exact_isometry": (
            basis_rows * gram * basis_rows.T == standard_gram
            and standard_to_roulleau.T * gram * standard_to_roulleau
            == standard_gram
        ),
    }
    if not all(stock_checks.values()):
        raise ArithmeticError("Roulleau curve/ample/marking stock changed")

    inclusion = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in carrier["lattice_and_embedding"]
            ["explicit_M_into_family_NS"]["basis_columns"]
        ]
    )
    gram_m = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in carrier["lattice_and_embedding"]
            ["original_divisor_lattice_M"]["gram"]
        ]
    )
    interior_m = sp.Matrix(
        [
            sp.Rational(value)
            for value in carrier["effectivity_and_divisor_interpretation"]
            ["strict_positive_chamber_interior"]["fixed_lattice_coordinates"]
        ]
    )
    standard_interior = inclusion * interior_m
    initial = standard_to_roulleau * standard_interior
    opposite = -initial
    marking_checks = {
        "index_four_carrier_embedding_replays": (
            inclusion.T * standard_gram * inclusion == gram_m
            and abs(int(inclusion.det())) == 4
        ),
        "strict_interior_has_archived_square": (
            _dot(gram_m, interior_m, interior_m) == 37418
            and _dot(gram, initial, initial) == 37418
        ),
        "plus_sign_is_the_published_ample_component": (
            _dot(gram, initial, D136) == 581903
            and _dot(gram, opposite, D136) == -581903
        ),
        "combined_initial_embedding_has_index_four": (
            abs(int((standard_to_roulleau * inclusion).det())) == 4
        ),
    }
    if not all(marking_checks.values()):
        raise ArithmeticError("positive-cone sign or carrier marking changed")

    current = initial
    weyl = sp.eye(12)
    firing_trace: list[dict[str, int | str]] = []
    computed_word: list[int] = []
    for _ in range(1000):
        pairings = [_dot(gram, current, curve) for curve in curves]
        negative = [index for index, pairing in enumerate(pairings) if pairing < 0]
        if not negative:
            break
        index = negative[0]
        label = index + 1
        pairing = pairings[index]
        height_before = _dot(gram, current, D136)
        reflection = sp.eye(12) + curves[index] * curves[index].T * gram
        following = reflection * current
        height_after = _dot(gram, following, D136)
        if not (pairing < 0 and 0 < height_after < height_before):
            raise ArithmeticError("Weyl firing is not strict negative descent")
        computed_word.append(label)
        firing_trace.append(
            {
                "step_one_based": len(computed_word),
                "curve": f"A{label}",
                "curve_one_based": label,
                "negative_pairing": pairing,
                "D136_height_before": height_before,
                "D136_height_after": height_after,
            }
        )
        current = following
        weyl = reflection * weyl
    else:
        raise ArithmeticError("Weyl reduction exceeded its hard finite bound")

    nef_interior = current
    nef_pairings = [_dot(gram, nef_interior, curve) for curve in curves]
    expected_nef_interior = [178, 267, 418, 209, 360, 302, 271, 240, 58, 31, 120, 89]
    expected_nef_pairings = [
        0, 62, 0, 0, 0, 27, 0, 0, 124, 178, 0, 0, 178, 240, 294,
        62, 418, 116, 240,
    ]
    standard_to_nef = weyl * standard_to_roulleau
    m_to_nef = standard_to_nef * inclusion
    weyl_checks = {
        "deterministic_word_is_the_frozen_87_word": (
            computed_word == WEYL_WORD_ONE_BASED and len(firing_trace) == 87
        ),
        "all_firings_strictly_lower_positive_D136_height": all(
            row["negative_pairing"] < 0
            and 0 < row["D136_height_after"] < row["D136_height_before"]
            for row in firing_trace
        ),
        "W_is_integral_determinant_minus_one_Gram_isometry": (
            all(value.is_Integer for value in weyl)
            and weyl.det() == -1
            and weyl.T * gram * weyl == gram
        ),
        "standard_to_nef_is_unimodular_Gram_isometry": (
            standard_to_nef.det() == 1
            and standard_to_nef.T * gram * standard_to_nef == standard_gram
        ),
        "combined_M_embedding_is_exact_index_four": (
            abs(int(m_to_nef.det())) == 4
            and m_to_nef.T * gram * m_to_nef == gram_m
        ),
        "final_nef_vector_is_frozen": (
            _vector_payload(nef_interior) == expected_nef_interior
            and nef_pairings == expected_nef_pairings
        ),
        "final_vector_is_square_37418_closed_nef": (
            _dot(gram, nef_interior, nef_interior) == 37418
            and _dot(gram, nef_interior, D136) == 3352
            and min(nef_pairings) == 0
        ),
    }
    if not all(weyl_checks.values()):
        raise ArithmeticError("explicit H11/29 Weyl reduction changed")

    minus_two_records: list[dict[str, Any]] = []
    minus_four_records: list[dict[str, Any]] = []
    transported_roots: list[sp.Matrix] = []
    original_pairings: list[int] = []
    transported_pairings: list[int] = []
    transport_checks_flat: list[bool] = []

    for source_record in carrier["lattice_and_embedding"][
        "facet_reflection_extensions"
    ]["records"]:
        facet = int(source_record["facet_zero_based"])
        standard_root = sp.Matrix(
            [sp.Rational(value) for value in source_record["primitive_NS_normal_coordinates"]]
        )
        initial_root = standard_to_roulleau * standard_root
        transported = weyl * initial_root
        transported_roots.append(transported)
        square = _dot(gram, transported, transported)
        divisibility = _content(gram * transported)
        initial_pairing = _dot(gram, initial, initial_root)
        nef_pairing = _dot(gram, nef_interior, transported)
        ample_degree = _dot(gram, D136, transported)
        reflection = _root_reflection(gram, transported)
        checks = {
            "combined_transport_agrees": transported
            == standard_to_nef * standard_root,
            "square_and_divisibility_preserved": (
                square == source_record["primitive_NS_square"]
                and divisibility == source_record["primitive_NS_divisibility"]
            ),
            "strict_oriented_halfspace_pairing_preserved": (
                initial_pairing == nef_pairing and nef_pairing > 0
            ),
            "reflection_integral_unimodular_involutive": (
                all(value.is_Integer for value in reflection)
                and reflection.det() == -1
                and reflection**2 == sp.eye(12)
            ),
            "reflection_preserves_Roulleau_NS_Gram": (
                reflection.T * gram * reflection == gram
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"transported facet {facet} failed")
        transport_checks_flat.extend(checks.values())
        original_pairings.append(initial_pairing)
        transported_pairings.append(nef_pairing)
        common: dict[str, Any] = {
            "facet_zero_based": facet,
            "standard_family_coordinates": _vector_payload(standard_root),
            "initial_Roulleau_coordinates": _vector_payload(initial_root),
            "nef_framed_coordinates": _vector_payload(transported),
            "square": square,
            "divisibility": divisibility,
            "nef_carrier_pairing": nef_pairing,
            "D136_degree": ample_degree,
            "reflection_matrix": _matrix_payload(reflection),
            "transport_checks": checks,
        }

        if square == -2:
            if ample_degree <= 0:
                raise ArithmeticError("oriented minus-two root is not effective")
            coefficients, descent_trace, terminal = _decompose_minus_two_root(
                transported, gram, curves, D136
            )
            reconstructed = sum(
                (
                    coefficient * curve
                    for coefficient, curve in zip(coefficients, curves)
                ),
                sp.zeros(12, 1),
            )
            sparse = [
                {"curve": f"A{index + 1}", "coefficient": coefficient}
                for index, coefficient in enumerate(coefficients)
                if coefficient
            ]
            is_published_curve = any(transported == curve for curve in curves)
            decomposition_checks = {
                "nonnegative_integral_coefficients": all(
                    isinstance(value, int) and value >= 0 for value in coefficients
                ),
                "exact_coordinate_equality": reconstructed == transported,
                "D136_degree_adds_exactly": ample_degree
                == sum(
                    coefficient * degree
                    for coefficient, degree in zip(coefficients, ample_degrees)
                ),
                "support_connected_in_complete_curve_graph": (
                    _support_is_connected(coefficients, wall_gram)
                ),
                "finite_strict_D136_descent": all(
                    row["negative_pairing"] < 0
                    and 0 < row["D136_degree_after"] < row["D136_degree_before"]
                    for row in descent_trace
                ),
                "irreducibility_decided_by_complete_stock": (
                    is_published_curve == (sum(coefficients) == 1)
                ),
            }
            if not all(decomposition_checks.values()):
                raise ArithmeticError(f"facet {facet} effective decomposition failed")
            minus_two_records.append(
                {
                    **common,
                    "effectivity": (
                        "effective divisor by the displayed nonnegative sum of "
                        "the complete published irreducible (-2)-curve stock"
                    ),
                    "decomposition_coefficients_A1_through_A19": coefficients,
                    "decomposition_sparse": sparse,
                    "total_component_multiplicity": sum(coefficients),
                    "maximum_component_coefficient": max(coefficients),
                    "descent_terminal_curve": f"A{terminal}",
                    "descent_trace": descent_trace,
                    "is_itself_an_irreducible_curve": is_published_curve,
                    "irreducible_curve_label": (
                        f"A{terminal}" if is_published_curve else None
                    ),
                    "decomposition_checks": decomposition_checks,
                }
            )
        elif square == -4:
            minus_four_records.append(
                {
                    **common,
                    "classification": (
                        "square-minus-four divisibility-two divisor reflection vector"
                    ),
                    "effectivity_status": (
                        "not asserted; no square-minus-two Riemann--Roch or "
                        "irreducible-curve inference is used"
                    ),
                }
            )
        else:
            raise ArithmeticError(f"unexpected transported square at facet {facet}")

    effective_checks = {
        "exactly_17_effective_square_minus_two_records": len(minus_two_records) == 17,
        "all_17_have_exact_complete_stock_decompositions": all(
            all(row["decomposition_checks"].values()) for row in minus_two_records
        ),
        "exactly_one_is_irreducible_A16": (
            [
                (row["facet_zero_based"], row["irreducible_curve_label"])
                for row in minus_two_records
                if row["is_itself_an_irreducible_curve"]
            ]
            == [(23, "A16")]
        ),
        "remaining_16_are_reducible": sum(
            not row["is_itself_an_irreducible_curve"] for row in minus_two_records
        )
        == 16,
        "component_multiplicity_range_is_1_through_18": (
            min(row["total_component_multiplicity"] for row in minus_two_records),
            max(row["total_component_multiplicity"] for row in minus_two_records),
        )
        == (1, 18),
        "maximum_component_coefficient_is_4": max(
            row["maximum_component_coefficient"] for row in minus_two_records
        )
        == 4,
        "all_minus_two_D136_degrees_are_positive": min(
            row["D136_degree"] for row in minus_two_records
        )
        == 2,
    }
    minus_four_checks = {
        "exactly_12_square_minus_four_records": len(minus_four_records) == 12,
        "all_have_divisibility_two": all(
            (row["square"], row["divisibility"]) == (-4, 2)
            for row in minus_four_records
        ),
        "exact_facet_partition": [
            row["facet_zero_based"] for row in minus_four_records
        ]
        == [0, 1, 6, 7, 8, 9, 12, 18, 20, 24, 27, 28],
        "no_effectivity_claim_attached": all(
            row["effectivity_status"].startswith("not asserted")
            for row in minus_four_records
        ),
    }
    mirror_checks = {
        "all_29_transports_green": (
            len(transport_checks_flat) == 29 * 5 and all(transport_checks_flat)
        ),
        "all_strict_carrier_pairings_preserved": (
            original_pairings == transported_pairings
            and min(transported_pairings) == 27
        ),
        "17_plus_12_is_complete_partition": (
            len(minus_two_records) + len(minus_four_records) == 29
            and sorted(
                [row["facet_zero_based"] for row in minus_two_records]
                + [row["facet_zero_based"] for row in minus_four_records]
            )
            == list(range(29))
        ),
        "effective_decomposition_audit_green": all(effective_checks.values()),
        "square_minus_four_separation_green": all(minus_four_checks.values()),
    }
    if not all(mirror_checks.values()):
        raise ArithmeticError("simultaneous 29-mirror effectivity audit failed")

    D136_mirror_pairings = [
        _dot(gram, D136, root) for root in transported_roots
    ]
    simultaneous_ample = nef_interior + D136
    simultaneous_curve_degrees = [
        _dot(gram, simultaneous_ample, curve) for curve in curves
    ]
    simultaneous_mirror_pairings = [
        _dot(gram, simultaneous_ample, root) for root in transported_roots
    ]
    ample_checks = {
        "D136_alone_has_three_nonpositive_imported_facets": [
            (index, pairing)
            for index, pairing in enumerate(D136_mirror_pairings)
            if pairing <= 0
        ]
        == [(0, 0), (7, 0), (8, -6)],
        "H_plus_D136_has_expected_primitive_coordinates": (
            _vector_payload(simultaneous_ample)
            == [180, 272, 427, 213, 370, 314, 286, 259, 67, 40, 126, 89]
            and _content(simultaneous_ample) == 1
        ),
        "H_plus_D136_has_square_44258": _dot(
            gram, simultaneous_ample, simultaneous_ample
        )
        == 44258,
        "H_plus_D136_is_strictly_inside_complete_ample_chamber": (
            min(simultaneous_curve_degrees) == 1
            and all(degree > 0 for degree in simultaneous_curve_degrees)
        ),
        "H_plus_D136_is_strict_on_all_29_imported_walls": (
            len(simultaneous_mirror_pairings) == 29
            and min(simultaneous_mirror_pairings) == 30
            and [
                index
                for index, pairing in enumerate(simultaneous_mirror_pairings)
                if pairing == 30
            ]
            == [5, 11]
        ),
        "k_equals_one_is_first_nonnegative_strict_member": (
            min(D136_mirror_pairings) < 0
            and min(simultaneous_mirror_pairings) > 0
        ),
    }
    if not all(ample_checks.values()):
        raise ArithmeticError("simultaneous ample/Coxeter interior audit failed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "After the exact positive marking and an 87-reflection Weyl transport, "
            "the complex A16 H11/29 carrier lies in Roulleau's complete closed "
            "nef chamber on NS=U+D8(-1)+A1(-1)^2.  All seventeen transported "
            "(-2) orientations are effective: facet 23 is the irreducible curve "
            "A16 and the other sixteen are explicit reducible divisors.  The "
            "twelve (-4,div2) vectors remain separate.  H+D136 is a primitive "
            "square-44258 ample class strict on all twenty-nine imported walls."
        ),
        "scope": {
            "positive_gain": (
                "an explicit nef-framed, effective-mirror, ample-interior "
                "realization on the equation-level BH 0.4.4.2 K3 family"
            ),
            "correct_marking_sign": (
                "after the standard-to-Roulleau bridge, the plus sign has positive "
                "D136-degree 581903; the negative sign is in the opposite component"
            ),
            "irreducibility_decision": (
                "the complete 19-curve stock identifies exactly one transported "
                "(-2) mirror with an irreducible curve, namely facet 23=A16"
            ),
            "not_full_nef": (
                "the imported chamber is not the full nef chamber: sixteen (-2) "
                "mirrors are reducible divisors and twelve primitive facets have "
                "square -4 rather than being K3 root walls"
            ),
            "square_minus_four": (
                "the twelve (-4,div2) reflection vectors are archived without an "
                "effectivity or irreducible-curve claim"
            ),
            "abstract_novelty": (
                "this strengthens a known finite-volume H11/29 polytope carrier; "
                "it does not claim a new abstract Coxeter polytope"
            ),
        },
        "sources": {
            "carrier": {
                "path": str(CARRIER.relative_to(ROOT)),
                "file_sha256": EXPECTED_SHA256["carrier"],
                "payload_sha256": EXPECTED_CARRIER_PAYLOAD_SHA256,
            },
            "Roulleau_complete_nef_stock": {
                "path": str(ROULLEAU.relative_to(ROOT)),
                "file_sha256": EXPECTED_SHA256["roulleau"],
                "inventory_path": str(INVENTORY.relative_to(ROOT)),
                "inventory_sha256": EXPECTED_SHA256["inventory"],
                "primary_reference": (
                    "Xavier Roulleau, An atlas of K3 surfaces with finite "
                    "automorphism group, EPIGA 6 (2022), article 19"
                ),
                "doi": "10.46298/epiga.2022.6286",
                "arxiv": "2003.08985",
                "pinned_primary_TeX_sha256": EXPECTED_ROULLEAU_TEX_SHA256,
                "pinned_subsection_lines": [8648, 8803],
            },
            "equation_level_family": {
                "path": str(FAMILY.relative_to(ROOT)),
                "file_sha256": EXPECTED_SHA256["family"],
                "BH_id": "0.4.4.2",
            },
            "checks": source_checks,
        },
        "published_nef_chamber_stock": {
            "lattice": "U+D8(-1)+A1(-1)^2",
            "coordinate_basis": "A1,...,A12",
            "gram": _matrix_payload(gram),
            "complete_irreducible_minus_two_curve_count": 19,
            "curve_coordinates_A1_through_A19": [
                _vector_payload(curve) for curve in curves
            ],
            "curve_intersection_gram": _matrix_payload(wall_gram),
            "ample_class": {
                "label": "D136",
                "coordinates": _vector_payload(D136),
                "square": 136,
                "degrees_on_A1_through_A19": ample_degrees,
            },
            "checks": stock_checks,
        },
        "exact_marking_and_sign": {
            "standard_family_gram": _matrix_payload(standard_gram),
            "standard_basis_rows_in_Roulleau_basis": _matrix_payload(basis_rows),
            "standard_to_Roulleau_matrix": _matrix_payload(standard_to_roulleau),
            "carrier_M_to_standard_family_matrix": _matrix_payload(inclusion),
            "chosen_sign": "+",
            "initial_Roulleau_interior_coordinates": _vector_payload(initial),
            "initial_D136_degree": _dot(gram, initial, D136),
            "opposite_sign_D136_degree": _dot(gram, opposite, D136),
            "checks": marking_checks,
        },
        "explicit_weyl_transport": {
            "reflection_convention": (
                "for A_i^2=-2, s_i(v)=v+(v.A_i)A_i; the word acts left-to-right"
            ),
            "deterministic_rule": "fire the first A_i with negative current pairing",
            "word_curve_labels_one_based": computed_word,
            "word_length": len(computed_word),
            "firing_trace": firing_trace,
            "W_matrix": _matrix_payload(weyl),
            "W_determinant": int(weyl.det()),
            "standard_family_to_nef_matrix": _matrix_payload(standard_to_nef),
            "carrier_M_to_nef_matrix": _matrix_payload(m_to_nef),
            "final_nef_carrier_coordinates": _vector_payload(nef_interior),
            "final_nef_carrier_square": _dot(gram, nef_interior, nef_interior),
            "final_nef_carrier_D136_degree": _dot(gram, nef_interior, D136),
            "final_pairings_with_A1_through_A19": nef_pairings,
            "zero_wall_curve_labels": [
                f"A{index + 1}"
                for index, pairing in enumerate(nef_pairings)
                if pairing == 0
            ],
            "checks": weyl_checks,
        },
        "explicit_simultaneous_ample_interior": {
            "family": "B_k=kH+D136 for nonnegative integers k",
            "H_coordinates": _vector_payload(nef_interior),
            "D136_coordinates": _vector_payload(D136),
            "D136_pairings_with_transported_facets_0_through_28": (
                D136_mirror_pairings
            ),
            "B1_equals_H_plus_D136": {
                "coordinates": _vector_payload(simultaneous_ample),
                "primitive": True,
                "square": _dot(gram, simultaneous_ample, simultaneous_ample),
                "D136_degree": _dot(gram, simultaneous_ample, D136),
                "degrees_on_A1_through_A19": simultaneous_curve_degrees,
                "minimum_complete_curve_degree": min(simultaneous_curve_degrees),
                "pairings_with_transported_facets_0_through_28": (
                    simultaneous_mirror_pairings
                ),
                "minimum_imported_facet_pairing": min(
                    simultaneous_mirror_pairings
                ),
                "minimum_facet_indices_zero_based": [5, 11],
                "classification": (
                    "primitive ample class strictly inside the transported "
                    "29-wall Coxeter chamber"
                ),
            },
            "minimality_within_integral_family": (
                "k=0 has pairings 0,0,-6 on facets 0,7,8; k=1 is strict on "
                "all twenty-nine imported facets"
            ),
            "checks": ample_checks,
        },
        "transported_mirrors": {
            "count": 29,
            "minimum_strict_carrier_pairing": min(transported_pairings),
            "square_minus_two_effective_count": len(minus_two_records),
            "square_minus_two_irreducible_count": sum(
                row["is_itself_an_irreducible_curve"] for row in minus_two_records
            ),
            "square_minus_two_reducible_count": sum(
                not row["is_itself_an_irreducible_curve"]
                for row in minus_two_records
            ),
            "square_minus_four_divisibility_two_count": len(minus_four_records),
            "square_minus_two_effective_decompositions": minus_two_records,
            "square_minus_four_divisor_classes": minus_four_records,
            "effective_checks": effective_checks,
            "square_minus_four_checks": minus_four_checks,
            "checks": mirror_checks,
        },
        "hostile_audit": {
            "status": "green",
            "exact_replays": [
                "the standard-family to Roulleau integral marking and its sign",
                "all 87 strict negative wall firings and ample-height decreases",
                "the Weyl Gram isometry and the index-four carrier embedding",
                "all 29 transported roots, halfspaces, and reflection matrices",
                "all 17 effective decompositions and connected supports",
                "the exact 1 irreducible plus 16 reducible decision",
                "the separate inventory of all twelve (-4,divisibility-two) vectors",
                "the primitive ample class H+D136 strict on all 19 curve and 29 Coxeter walls",
            ],
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
        if archived != payload:
            raise ArithmeticError("rebuilt H11/29 nef-effectivity archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
