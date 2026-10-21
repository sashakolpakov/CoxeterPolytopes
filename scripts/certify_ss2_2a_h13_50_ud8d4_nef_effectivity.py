#!/usr/bin/env python3
"""Replay a simultaneous nef-chamber effectivity theorem for SS2-2A H13/50.

The input carrier uses an isometry X from the doubled dual N=M^vee(2) to
the Roulleau basis of NS=U+D8(-1)+D4(-1).  Its archived interior is in the
positive-cone component opposite Roulleau's published ample class D_154.
We therefore choose the equally valid signed marking -X and then replay an
explicit word in reflections in the complete published stock A_1,...,A_20
of irreducible (-2)-curves.

The resulting Weyl isometry puts the carrier interior in the closed nef
chamber.  All fifty carrier mirrors are transported simultaneously.  Each
of the forty square-minus-two orientations is written exactly as a
nonnegative sum of the A_i; the remaining ten square-minus-four reflection
vectors are archived separately without an effectivity claim.
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
    ROOT / "data/enumerations/ss2_2a_h13_50_ud8d4_divisor_carrier.json"
)
ROULLEAU = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
OUTPUT = ROOT / "data/enumerations/ss2_2a_h13_50_ud8d4_nef_effectivity.json"

SCHEMA = "coxeter.ss2_2a_h13_50.ud8d4_nef_effectivity.v1"
EXPECTED_SHA256 = {
    "carrier": "b683aa65e929872747ff0cec7a4f64ee8395d95b3614aea6e8f97a5c883f8412",
    "roulleau": "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
}
EXPECTED_CARRIER_PAYLOAD_SHA256 = (
    "8d5c2e14661654b8ab4b9e97dc9d75233b2c604aa40418193c0aef3831c2faf4"
)
EXPECTED_ROULLEAU_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)
EXPECTED_SOURCE_INVENTORY_SHA256 = (
    "c3e9c0fba30c76c3e90bdb34fcf1d7d8b25574d1cd67ce1937b88e6905dbb317"
)

# Coordinates in Roulleau's basis A_1,...,A_14.  This is the ample class
# printed in the U+D8+D4 subsection of the pinned primary source.
D154 = sp.Matrix([1, 3, 6, 10, 7, 15, 14, 14, 15, 17, 8, 8, 0, 4])

# One-based labels in the published 20-curve stock.  If v.A_i<0, the step
# is the root reflection s_i(v)=v+(v.A_i)A_i.  The replay verifies every
# negative firing and every strict D154-height decrease rather than trusting
# this literal word.
WEYL_WORD_ONE_BASED = [
    1, 2, 1, 3, 2, 4, 3, 2, 1, 5,
    6, 4, 3, 2, 7, 6, 4, 3, 2, 1,
    8, 7, 6, 4, 3, 2, 1, 9, 8, 7,
    6, 4, 3, 2, 11, 13, 10, 12, 14, 10,
    9, 12, 13, 10, 14, 16, 1, 20, 1, 2,
    3, 4, 6, 7, 8, 16, 1, 2, 3, 20,
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("nonintegral matrix entered an integral certificate")
    return [[int(value) for value in matrix.row(row)] for row in range(matrix.rows)]


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
    support = [index for index, value in enumerate(coefficients) if value > 0]
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
) -> tuple[list[int], list[dict[str, Any]], int]:
    """Deterministically descend a positive-degree root to a published curve.

    At a negative wall A_i, put m=-(q.A_i)>0 and replace q by
    q-m A_i=s_i(q).  Thus the original class is the new root plus mA_i,
    and its ample degree strictly decreases by m(D154.A_i).  The stored
    finite trace is itself the termination certificate.
    """

    coefficients = [0] * len(curves)
    current = root
    trace: list[dict[str, Any]] = []
    seen: set[tuple[int, ...]] = set()
    for _ in range(1000):
        key = tuple(_vector_payload(current))
        if key in seen:
            raise ArithmeticError("root descent cycled")
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
            raise ArithmeticError("root descent did not reach the complete curve stock")
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
            raise ArithmeticError("invalid negative-wall root descent")
        coefficients[index] += multiplicity
        trace.append(
            {
                "curve": f"A{index + 1}",
                "curve_one_based": index + 1,
                "negative_pairing": pairing,
                "multiplicity_added": multiplicity,
                "D154_degree_before": height_before,
                "D154_degree_after": height_after,
            }
        )
        current = following
    raise ArithmeticError("root descent exceeded its hard finite bound")


def build_certificate() -> dict[str, Any]:
    paths = {"carrier": CARRIER, "roulleau": ROULLEAU}
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    roulleau_archive = json.loads(ROULLEAU.read_text(encoding="utf-8"))
    source = roulleau_archive["source"]

    carrier_without_digest = dict(carrier)
    embedded_carrier_digest = carrier_without_digest.pop("payload_sha256")
    source_checks = {
        "carrier_file_checksum_pinned": _sha256(CARRIER)
        == EXPECTED_SHA256["carrier"],
        "roulleau_file_checksum_pinned": _sha256(ROULLEAU)
        == EXPECTED_SHA256["roulleau"],
        "carrier_payload_digest_replays": (
            embedded_carrier_digest == EXPECTED_CARRIER_PAYLOAD_SHA256
            and _payload_sha256(carrier_without_digest) == embedded_carrier_digest
        ),
        "both_input_certificates_green": (
            carrier["all_checks_pass"] is True
            and roulleau_archive["all_checks_pass"] is True
        ),
        "correct_carrier_and_picard_lattice": (
            carrier["vinberg_finite_volume"]["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and carrier["vinberg_finite_volume"]["dimension"] == 13
            and carrier["primitive_mirrors"]["count"] == 50
            and source["lattice"] == "U+D8(-1)+D4(-1)"
            and source["picard_rank"] == 14
        ),
        "published_curve_stock_is_complete": (
            source["complete_minus_two_curve_count"] == 20
            and source["completeness"]
            == "Roulleau_published_complete_curve_configuration"
            and len(source["wall_coordinate_rows"]) == 20
        ),
        "primary_source_and_inventory_pinned": (
            source["roulleau_tex_sha256"] == EXPECTED_ROULLEAU_TEX_SHA256
            and source["source_inventory_sha256"]
            == EXPECTED_SOURCE_INVENTORY_SHA256
            and source["roulleau_tex_section_lines"] == [9669, 9827]
        ),
    }
    if not all(source_checks.values()):
        raise ArithmeticError(
            f"source handoff failed: {[key for key, ok in source_checks.items() if not ok]}"
        )

    gram = sp.Matrix(source["basis_gram"])
    wall_gram = sp.Matrix(source["wall_gram"])
    curves = [sp.Matrix(row) for row in source["wall_coordinate_rows"]]
    curve_matrix = sp.Matrix.hstack(*curves)
    computed_wall_gram = curve_matrix.T * gram * curve_matrix
    ample_degrees = [_dot(gram, D154, curve) for curve in curves]
    expected_ample_degrees = [1] * 12 + [17, 9, 17, 17, 17, 1, 9, 1]
    curve_stock_checks = {
        "first_14_curves_are_the_published_NS_basis": curve_matrix[:, :14]
        == sp.eye(14),
        "all_20_curves_have_square_minus_two": all(
            _dot(gram, curve, curve) == -2 for curve in curves
        ),
        "archived_wall_gram_recomputes": computed_wall_gram == wall_gram,
        "D154_has_published_square": _dot(gram, D154, D154) == 154,
        "D154_has_published_positive_curve_degrees": (
            ample_degrees == expected_ample_degrees and min(ample_degrees) == 1
        ),
    }
    if not all(curve_stock_checks.values()):
        raise ArithmeticError("published curve/ample stock changed")

    gram_n = sp.Matrix(
        carrier["scaled_dual_bridge"]["doubled_dual_N"]["gram"]
    )
    old_transport = sp.Matrix(
        carrier["scaled_dual_bridge"]["explicit_N_to_K3_isometry"]["matrix"]
    )
    interior_n = sp.Matrix(carrier["strict_carrier_interior"]["N_coordinates"])
    old_interior = sp.Matrix(
        carrier["strict_carrier_interior"]["K3_coordinates"]
    )
    signed_transport = -old_transport
    initial = signed_transport * interior_n
    marking_checks = {
        "old_X_reproduces_archived_interior": old_transport * interior_n
        == old_interior,
        "signed_marking_is_an_isometry": signed_transport.T * gram * signed_transport
        == gram_n,
        "signed_marking_is_unimodular": abs(int(signed_transport.det())) == 1,
        "old_interior_is_opposite_published_ample_component": _dot(
            gram, old_interior, D154
        )
        == -10008,
        "signed_interior_is_in_published_ample_component": (
            initial == -old_interior and _dot(gram, initial, D154) == 10008
        ),
        "signed_interior_keeps_square_2012": _dot(gram, initial, initial) == 2012,
    }
    if not all(marking_checks.values()):
        raise ArithmeticError("signed marking/component check changed")

    weyl = sp.eye(14)
    current = initial
    firing_trace: list[dict[str, int | str]] = []
    for step, label in enumerate(WEYL_WORD_ONE_BASED, start=1):
        curve = curves[label - 1]
        pairing = _dot(gram, current, curve)
        height_before = _dot(gram, current, D154)
        reflection = sp.eye(14) + curve * curve.T * gram
        following = reflection * current
        height_after = _dot(gram, following, D154)
        if not (pairing < 0 and 0 < height_after < height_before):
            raise ArithmeticError(f"Weyl firing {step} is not strict negative descent")
        firing_trace.append(
            {
                "step_one_based": step,
                "curve": f"A{label}",
                "curve_one_based": label,
                "negative_pairing": pairing,
                "D154_height_before": height_before,
                "D154_height_after": height_after,
            }
        )
        current = following
        weyl = reflection * weyl

    nef_interior = current
    nef_pairings = [_dot(gram, nef_interior, curve) for curve in curves]
    expected_nef_interior = [26, 52, 78, 120, 81, 162, 123, 84, 61, 38, -13, 3, 19, 6]
    expected_nef_pairings = [0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 64, 32, 0, 26, 38, 0, 32, 64, 38, 26]
    combined_n_to_nef = weyl * signed_transport
    signed_old_k3_to_nef = -weyl
    weyl_checks = {
        "word_has_60_strict_negative_firings": (
            len(WEYL_WORD_ONE_BASED) == len(firing_trace) == 60
            and all(row["negative_pairing"] < 0 for row in firing_trace)
        ),
        "D154_height_strictly_descends_and_stays_positive": all(
            0 < row["D154_height_after"] < row["D154_height_before"]
            for row in firing_trace
        ),
        "W_is_integral_unimodular_Gram_isometry": (
            all(value.is_Integer for value in weyl)
            and weyl.det() == 1
            and weyl.T * gram * weyl == gram
        ),
        "combined_signed_map_is_explicit_N_to_NS_isometry": (
            combined_n_to_nef == signed_old_k3_to_nef * old_transport
            and combined_n_to_nef.T * gram * combined_n_to_nef == gram_n
            and abs(int(combined_n_to_nef.det())) == 1
        ),
        "combined_map_takes_archived_N_interior_to_nef_vector": (
            combined_n_to_nef * interior_n == nef_interior
            and _vector_payload(nef_interior) == expected_nef_interior
        ),
        "final_vector_has_square_2012_and_positive_ample_degree": (
            _dot(gram, nef_interior, nef_interior) == 2012
            and _dot(gram, nef_interior, D154) == 1192
        ),
        "final_vector_is_in_closed_published_nef_chamber": (
            nef_pairings == expected_nef_pairings and min(nef_pairings) == 0
        ),
    }
    if not all(weyl_checks.values()):
        raise ArithmeticError("explicit Weyl transport changed")

    minus_two_records: list[dict[str, Any]] = []
    minus_four_records: list[dict[str, Any]] = []
    transported_roots: list[sp.Matrix] = []
    all_transport_checks: list[bool] = []
    original_pairings: list[int] = []
    transported_pairings: list[int] = []

    for source_record in carrier["primitive_mirrors"]["records"]:
        facet = int(source_record["facet_zero_based"])
        root_n = sp.Matrix(source_record["primitive_N_coordinates"])
        old_root = sp.Matrix(source_record["primitive_K3_coordinates"])
        signed_root = -old_root
        transported = weyl * signed_root
        transported_roots.append(transported)
        transported_from_n = combined_n_to_nef * root_n
        square = _dot(gram, transported, transported)
        divisibility = _content(gram * transported)
        old_pairing = _dot(gram, initial, signed_root)
        nef_pairing = _dot(gram, nef_interior, transported)
        ample_degree = _dot(gram, D154, transported)
        reflection = _root_reflection(gram, transported)
        common_checks = {
            "signed_input_matches_minus_old_X_image": signed_root
            == signed_transport * root_n,
            "combined_transport_agrees": transported == transported_from_n,
            "square_and_divisibility_preserved": (
                square == source_record["primitive_K3_square"]
                and divisibility == source_record["primitive_K3_divisibility"]
            ),
            "carrier_halfspace_pairing_preserved": (
                old_pairing == nef_pairing and nef_pairing > 0
            ),
            "reflection_remains_integral_Gram_preserving": (
                all(value.is_Integer for value in reflection)
                and reflection.det() == -1
                and reflection.T * gram * reflection == gram
            ),
        }
        all_transport_checks.extend(common_checks.values())
        original_pairings.append(old_pairing)
        transported_pairings.append(nef_pairing)

        common: dict[str, Any] = {
            "facet_zero_based": facet,
            "source_orbit_index_one_based": source_record[
                "source_orbit_index_one_based"
            ],
            "signed_initial_K3_coordinates": _vector_payload(signed_root),
            "nef_chamber_coordinates": _vector_payload(transported),
            "square": square,
            "divisibility": divisibility,
            "nef_carrier_pairing": nef_pairing,
            "D154_degree": ample_degree,
            "transport_checks": common_checks,
        }
        if square == -2:
            coefficients, descent_trace, terminal = _decompose_minus_two_root(
                transported, gram, curves, D154
            )
            reconstructed = sum(
                (coefficient * curve for coefficient, curve in zip(coefficients, curves)),
                sp.zeros(14, 1),
            )
            support = [
                index + 1 for index, coefficient in enumerate(coefficients)
                if coefficient > 0
            ]
            sparse = [
                {"curve": f"A{index + 1}", "coefficient": coefficient}
                for index, coefficient in enumerate(coefficients)
                if coefficient > 0
            ]
            is_published_curve = any(transported == curve for curve in curves)
            decomposition_checks = {
                "all_coefficients_nonnegative_integers": all(
                    isinstance(value, int) and value >= 0 for value in coefficients
                ),
                "exact_coordinate_equality": reconstructed == transported,
                "ample_degree_adds_exactly": ample_degree
                == sum(
                    coefficient * degree
                    for coefficient, degree in zip(coefficients, ample_degrees)
                ),
                "support_is_connected_in_published_curve_graph": (
                    _support_is_connected(coefficients, wall_gram)
                ),
                "descent_is_finite_strict_D154_descent": all(
                    0 < row["D154_degree_after"] < row["D154_degree_before"]
                    and row["negative_pairing"] < 0
                    for row in descent_trace
                ),
                "decomposition_is_nontrivial": sum(coefficients) >= 2,
                "not_one_of_complete_20_irreducible_minus_two_curves": (
                    not is_published_curve
                    and source["completeness"]
                    == "Roulleau_published_complete_curve_configuration"
                ),
            }
            if not all(decomposition_checks.values()):
                raise ArithmeticError(f"facet {facet} decomposition failed")
            minus_two_records.append(
                {
                    **common,
                    "effectivity": (
                        "effective reducible divisor, by the displayed nonnegative "
                        "sum of published irreducible (-2)-curves"
                    ),
                    "decomposition_coefficients_A1_through_A20": coefficients,
                    "decomposition_sparse": sparse,
                    "support_curve_labels_one_based": support,
                    "total_component_multiplicity": sum(coefficients),
                    "maximum_component_coefficient": max(coefficients),
                    "descent_terminal_curve": f"A{terminal}",
                    "descent_trace": descent_trace,
                    "is_itself_an_irreducible_curve": False,
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
                        "not asserted; no (-2)-root Riemann--Roch inference is used"
                    ),
                }
            )
        else:
            raise ArithmeticError(f"unexpected transported square at facet {facet}")

    effective_checks = {
        "exactly_40_square_minus_two_records": len(minus_two_records) == 40,
        "every_minus_two_class_has_replayed_effective_decomposition": all(
            all(row["decomposition_checks"].values()) for row in minus_two_records
        ),
        "none_of_the_40_is_itself_irreducible": all(
            row["is_itself_an_irreducible_curve"] is False
            for row in minus_two_records
        ),
        "all_40_have_positive_D154_degree": min(
            row["D154_degree"] for row in minus_two_records
        )
        > 0,
        "all_40_decompositions_have_connected_support": all(
            row["decomposition_checks"][
                "support_is_connected_in_published_curve_graph"
            ]
            for row in minus_two_records
        ),
        "component_multiplicity_range_is_2_through_28": (
            min(row["total_component_multiplicity"] for row in minus_two_records),
            max(row["total_component_multiplicity"] for row in minus_two_records),
        )
        == (2, 28),
        "maximum_coefficient_is_6": max(
            row["maximum_component_coefficient"] for row in minus_two_records
        )
        == 6,
    }
    minus_four_checks = {
        "exactly_10_square_minus_four_records": len(minus_four_records) == 10,
        "all_have_divisibility_two": all(
            (row["square"], row["divisibility"]) == (-4, 2)
            for row in minus_four_records
        ),
        "exact_facet_partition": [
            row["facet_zero_based"] for row in minus_four_records
        ]
        == [3, 6, 8, 12, 14, 16, 19, 21, 22, 23],
        "no_effectivity_claim_attached": all(
            row["effectivity_status"].startswith("not asserted")
            for row in minus_four_records
        ),
    }
    mirror_checks = {
        "all_50_transports_green": len(all_transport_checks) == 50 * 5
        and all(all_transport_checks),
        "all_strict_carrier_pairings_preserved": (
            original_pairings == transported_pairings
            and min(transported_pairings) == 6
        ),
        "40_plus_10_is_complete_partition": (
            len(minus_two_records) + len(minus_four_records) == 50
            and sorted(
                [row["facet_zero_based"] for row in minus_two_records]
                + [row["facet_zero_based"] for row in minus_four_records]
            )
            == list(range(50))
        ),
        "effective_decomposition_audit_green": all(effective_checks.values()),
        "square_minus_four_separation_green": all(minus_four_checks.values()),
    }
    if not all(mirror_checks.values()):
        raise ArithmeticError("simultaneous 50-mirror transport audit failed")

    # H itself is nef and strict for the imported Coxeter chamber, but lies on
    # ten genuine K3 curve walls.  Adding the published ample class moves into
    # the ample cone.  The first integral translate H+D154 still crosses one
    # imported wall; the next translate 2H+D154 is simultaneously ample and
    # strict for all fifty imported walls.
    first_ample_translate = nef_interior + D154
    simultaneous_ample = 2 * nef_interior + D154
    D154_mirror_pairings = [
        _dot(gram, D154, root) for root in transported_roots
    ]
    first_translate_curve_degrees = [
        _dot(gram, first_ample_translate, curve) for curve in curves
    ]
    first_translate_mirror_pairings = [
        _dot(gram, first_ample_translate, root) for root in transported_roots
    ]
    simultaneous_curve_degrees = [
        _dot(gram, simultaneous_ample, curve) for curve in curves
    ]
    simultaneous_mirror_pairings = [
        _dot(gram, simultaneous_ample, root) for root in transported_roots
    ]
    ample_interior_checks = {
        "D154_alone_fails_exactly_facets_8_14_22": [
            index for index, pairing in enumerate(D154_mirror_pairings)
            if pairing <= 0
        ]
        == [8, 14, 22],
        "H_plus_D154_is_ample": (
            _dot(gram, first_ample_translate, first_ample_translate) == 4550
            and min(first_translate_curve_degrees) == 1
        ),
        "H_plus_D154_fails_only_imported_facet_22_by_minus_two": (
            [
                (index, pairing)
                for index, pairing in enumerate(first_translate_mirror_pairings)
                if pairing <= 0
            ]
            == [(22, -2)]
        ),
        "2H_plus_D154_has_expected_primitive_coordinates": (
            _vector_payload(simultaneous_ample)
            == [53, 107, 162, 250, 169, 339, 260, 182, 137, 93, -18, 14, 38, 16]
            and _content(simultaneous_ample) == 1
        ),
        "2H_plus_D154_has_square_12970": _dot(
            gram, simultaneous_ample, simultaneous_ample
        )
        == 12970,
        "2H_plus_D154_is_strictly_inside_complete_ample_chamber": (
            min(simultaneous_curve_degrees) == 1
            and all(degree > 0 for degree in simultaneous_curve_degrees)
        ),
        "2H_plus_D154_is_strict_on_all_50_imported_walls": (
            len(simultaneous_mirror_pairings) == 50
            and min(simultaneous_mirror_pairings) == 4
            and simultaneous_mirror_pairings[22] == 4
        ),
        "k_equals_2_is_first_nonnegative_integral_translate_that_is_Coxeter_strict": (
            min(D154_mirror_pairings) < 0
            and min(first_translate_mirror_pairings) < 0
            and min(simultaneous_mirror_pairings) > 0
            and all(pairing > 0 for pairing in transported_pairings)
        ),
    }
    if not all(ample_interior_checks.values()):
        raise ArithmeticError("simultaneous ample/Coxeter interior audit failed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "After the explicit signed marking and 60-letter Weyl transport, "
            "the SS2-2A H13/50 carrier interior lies in the published closed nef "
            "chamber of NS=U+D8(-1)+D4(-1); all forty transported (-2) classes "
            "are explicit effective reducible divisors in Roulleau's complete "
            "20-curve stock; moreover 2H+D154 is an explicit primitive ample "
            "class strictly inside all 50 imported walls, while the ten "
            "(-4,div2) vectors remain separate algebraic divisor reflections "
            "without an effectivity claim."
        ),
        "scope": {
            "signed_marking_is_essential": (
                "the old X-image has D154-degree -10008, so it is in the opposite "
                "positive-cone component; -X is chosen before applying Weyl reflections"
            ),
            "nef_means_closed_chamber": (
                "the transported square-2012 vector has nonnegative intersection "
                "with all 20 published irreducible (-2)-curves and lies on ten walls"
            ),
            "decompositions": (
                "one exact nonnegative decomposition is certified for each of the "
                "40 classes; uniqueness of decompositions is not claimed"
            ),
            "irreducibility": (
                "none of the 40 classes is among the complete published stock of "
                "20 irreducible (-2)-curves, so all 40 are reducible"
            ),
            "simultaneous_ample_interior": (
                "2H+D154 is primitive, ample, and strict on all 50 imported "
                "walls; H+D154 fails only facet 22, so coefficient 2 is minimal "
                "among the nonnegative integral family kH+D154"
            ),
            "square_minus_four": (
                "the 10 (-4,div2) reflection vectors are recorded separately; no "
                "effectivity or irreducible-curve assertion is made for them"
            ),
            "abstract_polytope_novelty": (
                "this strengthens the divisor realization of the already certified "
                "H13/50 chamber; it is not a new abstract Coxeter polytope"
            ),
        },
        "sources": {
            "carrier": {
                "path": str(CARRIER.relative_to(ROOT)),
                "file_sha256": EXPECTED_SHA256["carrier"],
                "payload_sha256": EXPECTED_CARRIER_PAYLOAD_SHA256,
            },
            "published_curve_stock": {
                "path": str(ROULLEAU.relative_to(ROOT)),
                "file_sha256": EXPECTED_SHA256["roulleau"],
                "primary_reference": (
                    "Xavier Roulleau, An atlas of K3 surfaces with finite "
                    "automorphism group, EPIGA 6 (2022), article 19"
                ),
                "doi": "10.46298/epiga.2022.6286",
                "arxiv": "2003.08985",
                "pinned_primary_TeX_sha256": EXPECTED_ROULLEAU_TEX_SHA256,
                "pinned_subsection_lines": source["roulleau_tex_section_lines"],
                "pinned_source_inventory_sha256": EXPECTED_SOURCE_INVENTORY_SHA256,
            },
            "checks": source_checks,
        },
        "published_nef_chamber_stock": {
            "lattice": "U+D8(-1)+D4(-1)",
            "coordinate_basis": "A1,...,A14",
            "gram": _matrix_payload(gram),
            "complete_irreducible_minus_two_curve_count": 20,
            "curve_coordinates_A1_through_A20": [
                _vector_payload(curve) for curve in curves
            ],
            "curve_intersection_gram": _matrix_payload(wall_gram),
            "ample_class": {
                "label": "D154",
                "coordinates": _vector_payload(D154),
                "square": 154,
                "degrees_on_A1_through_A20": ample_degrees,
            },
            "checks": curve_stock_checks,
        },
        "signed_marking": {
            "old_X_interior_coordinates": _vector_payload(old_interior),
            "old_X_interior_D154_degree": -10008,
            "chosen_marking": "-X",
            "signed_initial_interior_coordinates": _vector_payload(initial),
            "signed_initial_D154_degree": 10008,
            "signed_N_to_K3_matrix": _matrix_payload(signed_transport),
            "checks": marking_checks,
        },
        "explicit_weyl_transport": {
            "reflection_convention": (
                "for A_i^2=-2, s_i(v)=v+(v.A_i)A_i; listed word acts left-to-right"
            ),
            "word_curve_labels_one_based": WEYL_WORD_ONE_BASED,
            "word_length": len(WEYL_WORD_ONE_BASED),
            "firing_trace": firing_trace,
            "W_matrix": _matrix_payload(weyl),
            "W_determinant": int(weyl.det()),
            "signed_old_K3_to_nef_matrix_minus_W": _matrix_payload(
                signed_old_k3_to_nef
            ),
            "combined_N_to_nef_matrix_W_times_minus_X": _matrix_payload(
                combined_n_to_nef
            ),
            "final_nef_carrier_coordinates": _vector_payload(nef_interior),
            "final_nef_carrier_square": _dot(gram, nef_interior, nef_interior),
            "final_nef_carrier_D154_degree": _dot(gram, nef_interior, D154),
            "final_pairings_with_A1_through_A20": nef_pairings,
            "zero_wall_curve_labels": [
                f"A{index + 1}" for index, pairing in enumerate(nef_pairings)
                if pairing == 0
            ],
            "checks": weyl_checks,
        },
        "explicit_simultaneous_ample_interior": {
            "family": "B_k = k H + D154 for nonnegative integers k",
            "H_coordinates": _vector_payload(nef_interior),
            "D154_coordinates": _vector_payload(D154),
            "D154_pairings_with_transported_facets_0_through_49": (
                D154_mirror_pairings
            ),
            "B1_equals_H_plus_D154": {
                "coordinates": _vector_payload(first_ample_translate),
                "square": _dot(
                    gram, first_ample_translate, first_ample_translate
                ),
                "degrees_on_A1_through_A20": first_translate_curve_degrees,
                "pairings_with_transported_facets_0_through_49": (
                    first_translate_mirror_pairings
                ),
                "unique_nonpositive_imported_facet": {
                    "facet_zero_based": 22,
                    "pairing": -2,
                },
            },
            "B2_equals_2H_plus_D154": {
                "coordinates": _vector_payload(simultaneous_ample),
                "primitive": True,
                "square": _dot(gram, simultaneous_ample, simultaneous_ample),
                "D154_degree": _dot(gram, simultaneous_ample, D154),
                "degrees_on_A1_through_A20": simultaneous_curve_degrees,
                "minimum_complete_curve_degree": min(simultaneous_curve_degrees),
                "pairings_with_transported_facets_0_through_49": (
                    simultaneous_mirror_pairings
                ),
                "minimum_imported_facet_pairing": min(
                    simultaneous_mirror_pairings
                ),
                "minimum_facet_zero_based": simultaneous_mirror_pairings.index(
                    min(simultaneous_mirror_pairings)
                ),
                "classification": (
                    "primitive ample class strictly inside the transported "
                    "50-wall Coxeter chamber"
                ),
            },
            "minimality_within_integral_family": (
                "k=0 fails facets 8,14,22; k=1 fails only facet 22 by -2; "
                "k=2 is strict on every imported facet"
            ),
            "checks": ample_interior_checks,
        },
        "transported_mirrors": {
            "count": 50,
            "minimum_strict_carrier_pairing": min(transported_pairings),
            "square_minus_two_effective_count": len(minus_two_records),
            "square_minus_two_irreducible_count": sum(
                row["is_itself_an_irreducible_curve"] for row in minus_two_records
            ),
            "square_minus_two_reducible_count": len(minus_two_records),
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
                "positive-cone component sign against the published ample D154",
                "all 60 negative wall firings and strict ample-height decreases",
                "the integral unimodular Gram isometry W and combined map W(-X)",
                "all 50 transported roots, strict carrier pairings, and reflections",
                "all 40 nonnegative coordinate decompositions and connected supports",
                "exclusion from the complete 20-curve irreducible stock",
                "separate inventory of all ten (-4,divisibility-two) vectors",
                "primitive ample class 2H+D154 strict on all 20 curve and 50 Coxeter walls",
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
            raise ArithmeticError("rebuilt nef-effectivity archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
