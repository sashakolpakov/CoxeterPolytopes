#!/usr/bin/env python3
"""Certify nef effectivity and an ample interior for A18-C06 H13/31.

The equation-level carrier is marked in the standard fiber basis of
U+E8(-1)+D4(-1).  Roulleau's complete 15-curve chamber supplies a second,
geometric basis of the same lattice.  The archived carrier marking initially
lands in the positive-cone component opposite Roulleau's ample class D506, so
we negate the marking and then replay a deterministic word in the 15 curve
reflections.

The endpoint is in the closed nef chamber.  All 31 carrier normals are moved
simultaneously.  Every square-minus-two orientation is decomposed into the
complete curve stock, which also decides irreducibility; the five square-minus-
four reflection vectors remain separate.  Finally, with H twice the rational
carrier interior, D506 + 7H is a primitive ample class strict on all 31 walls,
and seven is minimal in the integral family D506 + kH, k >= 0.
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
CARRIER = ROOT / "data/enumerations/a18_c06_ue8d4_divisor_overlattice_carrier.json"
FAMILY = ROOT / "data/enumerations/ue8d4_order6_weierstrass_family.json"
ROULLEAU = ROOT / "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json"
INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
OUTPUT = ROOT / "data/enumerations/a18_c06_ue8d4_nef_effectivity.json"

SCHEMA = "coxeter.a18_c06.ue8d4_nef_effectivity.v1"
EXPECTED_SHA256 = {
    "carrier": "c4a5ee745791ac5a52e182d21d32bc0c4ee77cf7753ff5b2763eb6eebf5770bc",
    "family": "0ce8340494520835daf4df1aaaa57239640759111fb45682dc9b6b330b618fff",
    "roulleau": "de41f13f9b47fa2c1e82bbd9255da74667694658e5acac2461054bb5d2996f04",
    "inventory": "c3e9c0fba30c76c3e90bdb34fcf1d7d8b25574d1cd67ce1937b88e6905dbb317",
}
EXPECTED_CARRIER_PAYLOAD_SHA256 = (
    "d34d7a9b36dd9d8476db53e253530635d2e51d0d4199e12e147c46b385eebf05"
)
EXPECTED_ROULLEAU_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)

# Roulleau, rank-14 U+E8+D4 subsection, in the basis A1,...,A14.
D506 = sp.Matrix([46, 93, 141, 70, 120, 100, 81, 63, 46, 30, 15, 1, -12, 0])
EXPECTED_D506_DEGREES = [1] * 12 + [25, 1, 1]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _jsonify(value: Any) -> Any:
    """Normalize exact SymPy values and tuples before archival."""

    if isinstance(value, dict):
        return {key: _jsonify(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(item) for item in value]
    if value is sp.true:
        return True
    if value is sp.false:
        return False
    if isinstance(value, sp.Integer):
        return int(value)
    if isinstance(value, sp.Rational):
        return int(value) if value.q == 1 else str(value)
    return value


def _scalar_payload(value: sp.Expr) -> int | str:
    value = sp.Rational(value)
    return int(value) if value.q == 1 else str(value)


def _vector_payload(vector: sp.Matrix) -> list[int | str]:
    if vector.cols != 1:
        raise ArithmeticError("expected a column vector")
    return [_scalar_payload(value) for value in vector]


def _integer_vector_payload(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1 or not all(value.is_Integer for value in vector):
        raise ArithmeticError("nonintegral vector entered an integral field")
    return [int(value) for value in vector]


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("nonintegral matrix entered an integral field")
    return [[int(value) for value in matrix.row(row)] for row in range(matrix.rows)]


def _dot(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _integer_dot(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> int:
    value = _dot(gram, left, right)
    if value.q != 1:
        raise ArithmeticError("nonintegral divisor pairing")
    return int(value)


def _content(vector: sp.Matrix) -> int:
    if not all(value.is_Integer for value in vector):
        raise ArithmeticError("content requested for a nonintegral vector")
    return math.gcd(*(abs(int(value)) for value in vector))


def _root_reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _dot(gram, root, root)
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


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


def _standard_to_roulleau(curves: list[sp.Matrix]) -> sp.Matrix:
    """Columns give [F,F+O,E8-roots,D4-roots] in the A_i basis.

    The I0* fiber is A11+2A12+A13+A14+A15 and O=A10.  The E8 roots are
    A1,A2,A3,A5,A6,A7,A8,A4; the D4 roots are A13,A12,A14,A15.
    """

    fiber = curves[10] + 2 * curves[11] + curves[12] + curves[13] + curves[14]
    return sp.Matrix.hstack(
        fiber,
        fiber + curves[9],
        *(curves[index] for index in [0, 1, 2, 4, 5, 6, 7, 3]),
        *(curves[index] for index in [12, 11, 13, 14]),
    )


def _reduce_to_nef(
    initial: sp.Matrix,
    gram: sp.Matrix,
    curves: list[sp.Matrix],
    ample: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, list[int], list[dict[str, int | str]]]:
    """Fire the most negative curve wall, breaking ties by its label."""

    current = initial
    action = sp.eye(gram.rows)
    word: list[int] = []
    trace: list[dict[str, int | str]] = []
    for step in range(1, 10_001):
        pairings = [_dot(gram, current, curve) for curve in curves]
        negative = [
            (pairing, index)
            for index, pairing in enumerate(pairings)
            if pairing < 0
        ]
        if not negative:
            return action, current, word, trace
        pairing, index = min(negative)
        if pairing.q != 1:
            raise ArithmeticError("Weyl firing pairing became nonintegral")
        height_before = _integer_dot(gram, current, ample)
        reflection = sp.eye(gram.rows) + curves[index] * curves[index].T * gram
        following = reflection * current
        height_after = _integer_dot(gram, following, ample)
        if not (pairing < 0 and 0 < height_after < height_before):
            raise ArithmeticError(f"Weyl firing {step} is not strict descent")
        word.append(index + 1)
        trace.append(
            {
                "step_one_based": step,
                "curve": f"A{index + 1}",
                "curve_one_based": index + 1,
                "negative_pairing": int(pairing),
                "D506_height_before": height_before,
                "D506_height_after": height_after,
            }
        )
        current = following
        action = reflection * action
    raise ArithmeticError("nef reduction exceeded its hard finite bound")


def _decompose_minus_two_root(
    root: sp.Matrix,
    gram: sp.Matrix,
    curves: list[sp.Matrix],
    ample: sp.Matrix,
) -> tuple[list[int], list[dict[str, int | str]], int]:
    """Descend an effective root to the complete irreducible curve stock."""

    coefficients = [0] * len(curves)
    current = root
    trace: list[dict[str, int | str]] = []
    seen: set[tuple[int, ...]] = set()
    for _ in range(1000):
        key = tuple(_integer_vector_payload(current))
        if key in seen:
            raise ArithmeticError("minus-two root descent cycled")
        seen.add(key)

        terminal = next(
            (index for index, curve in enumerate(curves) if current == curve), None
        )
        if terminal is not None:
            coefficients[terminal] += 1
            return coefficients, trace, terminal + 1

        pairings = [_integer_dot(gram, current, curve) for curve in curves]
        negative = [index for index, pairing in enumerate(pairings) if pairing < 0]
        if not negative:
            raise ArithmeticError("root descent did not reach the complete stock")
        index = negative[0]
        pairing = pairings[index]
        multiplicity = -pairing
        height_before = _integer_dot(gram, current, ample)
        following = current - multiplicity * curves[index]
        height_after = _integer_dot(gram, following, ample)
        if not (
            multiplicity > 0
            and _integer_dot(gram, following, following) == -2
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
                "D506_degree_before": height_before,
                "D506_degree_after": height_after,
            }
        )
        current = following
    raise ArithmeticError("root descent exceeded its hard finite bound")


def build_certificate() -> dict[str, Any]:
    paths = {
        "carrier": CARRIER,
        "family": FAMILY,
        "roulleau": ROULLEAU,
        "inventory": INVENTORY,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    roulleau = json.loads(ROULLEAU.read_text(encoding="utf-8"))
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    source = roulleau["source"]
    stock = next(row for row in inventory["stocks"] if row["key"] == "R14_UE8D4")

    carrier_without_digest = dict(carrier)
    embedded_carrier_digest = carrier_without_digest.pop("payload_sha256")
    source_checks = {
        "four_input_archives_green": (
            carrier["all_checks_pass"] is True
            and family["all_checks_pass"] is True
            and roulleau["all_checks_pass"] is True
            and inventory["all_checks_pass"] is True
        ),
        "carrier_payload_digest_replays": (
            embedded_carrier_digest == EXPECTED_CARRIER_PAYLOAD_SHA256
            and _payload_sha256(carrier_without_digest) == embedded_carrier_digest
        ),
        "correct_equation_level_carrier": (
            carrier["vinberg_finite_volume"]["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and carrier["vinberg_finite_volume"]["dimension"] == 13
            and carrier["primitive_mirrors"]["count"] == 31
            and family["neron_severi"]["lower_bound_lattice"]
            == "U+E8(-1)+D4(-1)"
        ),
        "published_stock_is_complete": (
            stock["picard_rank"] == 14
            and stock["wall_count"] == 15
            and stock["source_completeness_method"]
            == "direct_curve_enumeration_and_completeness_test"
            and stock["reconstruction_status"] == "already_exactly_archived"
        ),
        "primary_source_is_pinned": (
            source["tex_sha256"] == EXPECTED_ROULLEAU_TEX_SHA256
            and source["tex_section_lines"] == [9586, 9668]
            and stock["source"]["section_lines"] == [9586, 9668]
        ),
    }
    if not all(source_checks.values()):
        raise ArithmeticError(
            f"source gate failed: {[key for key, ok in source_checks.items() if not ok]}"
        )

    gram = sp.Matrix(source["basis_gram"])
    wall_gram = sp.Matrix(source["wall_gram"])
    curves = [sp.Matrix(row) for row in source["wall_coordinate_rows"]]
    curve_matrix = sp.Matrix.hstack(*curves)
    standard_gram = sp.Matrix(family["neron_severi"]["gram"])
    standard_to_roulleau = _standard_to_roulleau(curves)
    ample_degrees = [_integer_dot(gram, D506, curve) for curve in curves]
    fiber = standard_to_roulleau[:, 0]
    curve_stock_checks = {
        "first_14_curves_are_the_published_NS_basis": (
            curve_matrix[:, :14] == sp.eye(14)
        ),
        "all_15_curves_have_square_minus_two": all(
            _integer_dot(gram, curve, curve) == -2 for curve in curves
        ),
        "archived_wall_gram_recomputes": (
            curve_matrix.T * gram * curve_matrix == wall_gram
        ),
        "D506_has_published_square_and_degrees": (
            _integer_dot(gram, D506, D506) == 506
            and ample_degrees == EXPECTED_D506_DEGREES
        ),
        "geometric_fiber_is_I0_star_and_meets_A10_once": (
            fiber
            == curves[10] + 2 * curves[11] + curves[12] + curves[13] + curves[14]
            and _integer_dot(gram, fiber, fiber) == 0
            and _integer_dot(gram, fiber, curves[9]) == 1
        ),
        "fiber_basis_bridge_is_integral_unimodular": (
            all(value.is_Integer for value in standard_to_roulleau)
            and standard_to_roulleau.det() == -1
        ),
        "fiber_basis_bridge_isometry": (
            standard_to_roulleau.T * gram * standard_to_roulleau
            == standard_gram
        ),
    }
    if not all(curve_stock_checks.values()):
        raise ArithmeticError("Roulleau curve/fiber-basis bridge changed")

    old_standard_interior = sp.Matrix(
        [
            sp.Rational(value)
            for value in carrier["effectivity_and_scope"][
                "interior_in_family_coordinates"
            ]
        ]
    )
    old_roulleau_interior = standard_to_roulleau * old_standard_interior
    signed_standard_to_roulleau = -standard_to_roulleau
    signed_initial = -old_roulleau_interior
    marking_checks = {
        "old_interior_has_square_385": (
            _dot(standard_gram, old_standard_interior, old_standard_interior) == 385
            and _dot(gram, old_roulleau_interior, old_roulleau_interior) == 385
        ),
        "old_marking_is_opposite_published_ample_component": (
            _integer_dot(gram, old_roulleau_interior, D506) == -15889
        ),
        "signed_marking_enters_published_ample_component": (
            _integer_dot(gram, signed_initial, D506) == 15889
        ),
        "signed_marking_is_integral_unimodular_isometry": (
            all(value.is_Integer for value in signed_standard_to_roulleau)
            and signed_standard_to_roulleau.det() == -1
            and signed_standard_to_roulleau.T
            * gram
            * signed_standard_to_roulleau
            == standard_gram
        ),
    }
    if not all(marking_checks.values()):
        raise ArithmeticError("positive-cone component correction changed")

    weyl, nef_interior, word, firing_trace = _reduce_to_nef(
        signed_initial, gram, curves, D506
    )
    nef_pairings = [_integer_dot(gram, nef_interior, curve) for curve in curves]
    expected_nef_interior: list[int | str] = [
        38, 77, 116, 58, 98, 80, 63, 46, 34, 22, 11, 0, "-9/2", "-11/2"
    ]
    expected_nef_pairings = [1, 0, 1, 0, 0, 1, 0, 5, 0, 1, 0, 1, 9, 11, 0]
    combined_standard_to_nef = weyl * signed_standard_to_roulleau
    H = 2 * nef_interior
    H_curve_degrees = [_integer_dot(gram, H, curve) for curve in curves]
    weyl_checks = {
        "word_has_235_strict_negative_firings": (
            len(word) == len(firing_trace) == 235
            and all(row["negative_pairing"] < 0 for row in firing_trace)
        ),
        "D506_height_strictly_descends_and_stays_positive": all(
            0 < row["D506_height_after"] < row["D506_height_before"]
            for row in firing_trace
        ),
        "W_is_integral_Gram_isometry_of_determinant_minus_one": (
            all(value.is_Integer for value in weyl)
            and weyl.det() == -1
            and weyl.T * gram * weyl == gram
        ),
        "combined_map_is_integral_unimodular_isometry": (
            all(value.is_Integer for value in combined_standard_to_nef)
            and combined_standard_to_nef.det() == 1
            and combined_standard_to_nef.T * gram * combined_standard_to_nef
            == standard_gram
        ),
        "exact_nef_endpoint": (
            _vector_payload(nef_interior) == expected_nef_interior
            and nef_pairings == expected_nef_pairings
        ),
        "endpoint_is_square_385_closed_nef": (
            _dot(gram, nef_interior, nef_interior) == 385
            and _integer_dot(gram, nef_interior, D506) == 525
            and min(nef_pairings) == 0
        ),
        "twice_endpoint_is_primitive_integral_big_nef_H": (
            _integer_vector_payload(H)
            == [76, 154, 232, 116, 196, 160, 126, 92, 68, 44, 22, 0, -9, -11]
            and _content(H) == 1
            and _integer_dot(gram, H, H) == 1540
            and _integer_dot(gram, H, D506) == 1050
            and H_curve_degrees == [2, 0, 2, 0, 0, 2, 0, 10, 0, 2, 0, 2, 18, 22, 0]
        ),
    }
    if not all(weyl_checks.values()):
        raise ArithmeticError(
            f"Weyl gate failed: {[key for key, ok in weyl_checks.items() if not ok]}"
        )

    minus_two_records: list[dict[str, Any]] = []
    minus_four_records: list[dict[str, Any]] = []
    transported_roots: list[sp.Matrix] = []
    all_transport_checks: list[bool] = []

    for source_record in carrier["primitive_mirrors"]["records"]:
        facet = int(source_record["facet_zero_based"])
        root_standard = sp.Matrix(
            [sp.Rational(value) for value in source_record["oriented_primitive_NS_coordinates"]]
        )
        signed_root = signed_standard_to_roulleau * root_standard
        transported = combined_standard_to_nef * root_standard
        transported_roots.append(transported)
        square = _integer_dot(gram, transported, transported)
        divisibility = _content(gram * transported)
        original_pairing = _dot(
            standard_gram, old_standard_interior, root_standard
        )
        transported_pairing = _dot(gram, nef_interior, transported)
        D506_degree = _integer_dot(gram, D506, transported)
        reflection = _root_reflection(gram, transported)
        common_checks = {
            "signed_input_transport_replays": (
                signed_root == -standard_to_roulleau * root_standard
            ),
            "combined_Weyl_transport_replays": transported == weyl * signed_root,
            "primitive_integral_transport": (
                all(value.is_Integer for value in transported)
                and _content(transported) == 1
            ),
            "square_and_divisibility_preserved": (
                square == source_record["primitive_NS_square"]
                and divisibility == source_record["primitive_NS_divisibility"]
            ),
            "strict_carrier_halfspace_pairing_preserved": (
                original_pairing == transported_pairing and transported_pairing > 0
            ),
            "reflection_is_integral_unimodular_Gram_isometry": (
                all(value.is_Integer for value in reflection)
                and reflection.det() == -1
                and reflection**2 == sp.eye(14)
                and reflection.T * gram * reflection == gram
            ),
        }
        all_transport_checks.extend(common_checks.values())
        common: dict[str, Any] = {
            "facet_zero_based": facet,
            "source_orbit_index": source_record["source_orbit_index"],
            "old_equation_fiber_basis_coordinates": _integer_vector_payload(root_standard),
            "signed_initial_Roulleau_coordinates": _integer_vector_payload(signed_root),
            "nef_chamber_coordinates": _integer_vector_payload(transported),
            "square": square,
            "divisibility": divisibility,
            "nef_carrier_pairing": _scalar_payload(transported_pairing),
            "primitive_H_pairing": _integer_dot(gram, H, transported),
            "D506_degree": D506_degree,
            "transport_checks": common_checks,
        }

        if square == -2:
            if D506_degree <= 0:
                raise ArithmeticError(f"facet {facet} has the wrong effective sign")
            coefficients, descent_trace, terminal = _decompose_minus_two_root(
                transported, gram, curves, D506
            )
            reconstructed = sum(
                (
                    coefficient * curve
                    for coefficient, curve in zip(coefficients, curves, strict=True)
                ),
                sp.zeros(14, 1),
            )
            is_published_curve = any(transported == curve for curve in curves)
            sparse = [
                {"curve": f"A{index + 1}", "coefficient": coefficient}
                for index, coefficient in enumerate(coefficients)
                if coefficient > 0
            ]
            decomposition_checks = {
                "all_coefficients_nonnegative_integers": all(
                    isinstance(value, int) and value >= 0 for value in coefficients
                ),
                "exact_coordinate_equality": reconstructed == transported,
                "ample_degree_adds_exactly": D506_degree
                == sum(
                    coefficient * degree
                    for coefficient, degree in zip(
                        coefficients, ample_degrees, strict=True
                    )
                ),
                "support_is_connected_in_complete_curve_graph": (
                    _support_is_connected(coefficients, wall_gram)
                ),
                "finite_strict_D506_descent": all(
                    row["negative_pairing"] < 0
                    and 0 < row["D506_degree_after"] < row["D506_degree_before"]
                    for row in descent_trace
                ),
                "irreducibility_decided_by_complete_stock": (
                    (is_published_curve and sum(coefficients) == 1)
                    or (not is_published_curve and sum(coefficients) >= 2)
                ),
            }
            if not all(decomposition_checks.values()):
                raise ArithmeticError(f"facet {facet} decomposition failed")
            minus_two_records.append(
                {
                    **common,
                    "effectivity": (
                        "effective by the displayed nonnegative decomposition into "
                        "Roulleau's complete irreducible (-2)-curve stock"
                    ),
                    "decomposition_coefficients_A1_through_A15": coefficients,
                    "decomposition_sparse": sparse,
                    "total_component_multiplicity": sum(coefficients),
                    "maximum_component_coefficient": max(coefficients),
                    "descent_terminal_curve": f"A{terminal}",
                    "descent_trace": descent_trace,
                    "is_itself_an_irreducible_curve": is_published_curve,
                    "published_curve_label": (
                        f"A{terminal}" if is_published_curve else None
                    ),
                    "classification": (
                        "published irreducible (-2)-curve"
                        if is_published_curve
                        else "effective reducible (-2)-divisor"
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
                        "not asserted; no (-2)-root Riemann--Roch inference is used"
                    ),
                }
            )
        else:
            raise ArithmeticError(f"unexpected square at facet {facet}")

    irreducible_records = [
        row for row in minus_two_records if row["is_itself_an_irreducible_curve"]
    ]
    reducible_records = [
        row for row in minus_two_records if not row["is_itself_an_irreducible_curve"]
    ]
    effective_checks = {
        "exactly_26_square_minus_two_effective_classes": len(minus_two_records) == 26,
        "all_26_decompositions_replay": all(
            all(row["decomposition_checks"].values()) for row in minus_two_records
        ),
        "exact_irreducible_facet_and_curve_list": [
            (row["facet_zero_based"], row["published_curve_label"])
            for row in irreducible_records
        ]
        == [(0, "A1"), (3, "A12"), (8, "A3"), (10, "A10"), (14, "A6")],
        "exactly_5_irreducible_and_21_reducible": (
            len(irreducible_records) == 5 and len(reducible_records) == 21
        ),
        "reducible_multiplicity_range_3_through_55": (
            min(row["total_component_multiplicity"] for row in reducible_records),
            max(row["total_component_multiplicity"] for row in reducible_records),
        )
        == (3, 55),
        "maximum_component_coefficient_seven": max(
            row["maximum_component_coefficient"] for row in minus_two_records
        )
        == 7,
    }
    minus_four_checks = {
        "exactly_five_minus_four_divisibility_two_vectors": (
            len(minus_four_records) == 5
            and all(
                (row["square"], row["divisibility"]) == (-4, 2)
                for row in minus_four_records
            )
        ),
        "exact_minus_four_facet_partition": [
            row["facet_zero_based"] for row in minus_four_records
        ]
        == [4, 5, 6, 12, 13],
        "no_minus_four_effectivity_claim": all(
            row["effectivity_status"].startswith("not asserted")
            for row in minus_four_records
        ),
    }
    mirror_checks = {
        "all_31_simultaneous_transports_green": (
            len(all_transport_checks) == 31 * 6 and all(all_transport_checks)
        ),
        "complete_26_plus_5_partition": (
            len(minus_two_records) + len(minus_four_records) == 31
            and sorted(
                [row["facet_zero_based"] for row in minus_two_records]
                + [row["facet_zero_based"] for row in minus_four_records]
            )
            == list(range(31))
        ),
        "carrier_stays_strict_after_nef_transport": min(
            _dot(gram, nef_interior, root) for root in transported_roots
        )
        == 1,
        "effective_decomposition_and_irreducibility_gate": all(
            effective_checks.values()
        ),
        "minus_four_separation_gate": all(minus_four_checks.values()),
    }
    if not all(mirror_checks.values()):
        raise ArithmeticError("simultaneous 31-mirror gate failed")

    D506_mirror_pairings = [
        _integer_dot(gram, D506, root) for root in transported_roots
    ]
    H_mirror_pairings = [_integer_dot(gram, H, root) for root in transported_roots]
    first_strict_k = next(
        k
        for k in range(101)
        if min(
            _integer_dot(gram, D506 + k * H, root)
            for root in transported_roots
        )
        > 0
    )
    B6 = D506 + 6 * H
    B7 = D506 + 7 * H
    B6_mirror_pairings = [
        _integer_dot(gram, B6, root) for root in transported_roots
    ]
    B7_curve_degrees = [_integer_dot(gram, B7, curve) for curve in curves]
    B7_mirror_pairings = [
        _integer_dot(gram, B7, root) for root in transported_roots
    ]
    minimum_failures = [
        [
            (index, _integer_dot(gram, D506 + k * H, root))
            for index, root in enumerate(transported_roots)
            if _integer_dot(gram, D506 + k * H, root) <= 0
        ]
        for k in range(7)
    ]
    ample_checks = {
        "D506_is_published_ample": (
            _integer_dot(gram, D506, D506) == 506 and min(ample_degrees) == 1
        ),
        "H_is_primitive_big_nef_and_imported_strict": (
            _content(H) == 1
            and _integer_dot(gram, H, H) == 1540
            and min(H_curve_degrees) == 0
            and min(H_mirror_pairings) == 2
        ),
        "facet_five_controls_the_exact_boundary": (
            D506_mirror_pairings[5] == -24
            and H_mirror_pairings[5] == 4
            and [failures for failures in minimum_failures]
            == [[(5, -24)], [(5, -20)], [(5, -16)], [(5, -12)], [(5, -8)], [(5, -4)], [(5, 0)]]
        ),
        "k_seven_is_first_strict_nonnegative_integral_translate": (
            first_strict_k == 7
            and B6_mirror_pairings[5] == 0
            and min(B7_mirror_pairings) == 4
        ),
        "B7_has_expected_primitive_coordinates_and_square": (
            _integer_vector_payload(B7)
            == [578, 1171, 1765, 882, 1492, 1220, 963, 707, 522, 338, 169, 1, -75, -77]
            and _content(B7) == 1
            and _integer_dot(gram, B7, B7) == 90666
            and _integer_dot(gram, B7, D506) == 7856
        ),
        "B7_is_strictly_inside_complete_ample_chamber": (
            B7_curve_degrees
            == [15, 1, 15, 1, 1, 15, 1, 71, 1, 15, 1, 15, 151, 155, 1]
            and min(B7_curve_degrees) == 1
        ),
        "B7_is_strict_on_all_31_imported_walls": (
            B7_mirror_pairings
            == [
                15, 337, 17, 15, 152, 4, 358, 17, 15, 339, 15, 17, 150,
                360, 15, 17, 17, 17, 442, 442, 225, 225, 443, 650, 650,
                441, 225, 753, 651, 649, 755,
            ]
            and min(B7_mirror_pairings) == 4
        ),
    }
    if not all(ample_checks.values()):
        raise ArithmeticError(
            f"ample gate failed: {[key for key, ok in ample_checks.items() if not ok]}"
        )

    vinberg = carrier["vinberg_finite_volume"]
    vinberg_checks = {
        "classification": vinberg["classification"]
        == "finite_volume_noncompact_Coxeter_polytope",
        "dimension_facets": (vinberg["dimension"], vinberg["facets"]) == (13, 31),
        "literal_edge_gate": (
            vinberg["elliptic_edges"] == 28766 and vinberg["bad_edges"] == 0
        ),
        "literal_parabolic_gate": (
            vinberg["connected_affine_components"] == 529
            and vinberg["bad_parabolic_neighborhoods"] == 0
        ),
        "all_archived_Vinberg_checks": all(vinberg["checks"].values()),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("inherited literal Vinberg gate changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "After the necessary signed marking and a 235-reflection Weyl "
            "transport, the equation-level A18-C06 H13/31 carrier meets "
            "Roulleau's closed nef chamber on NS=U+E8(-1)+D4(-1).  Its 26 "
            "square-minus-two mirrors are effective: exactly five are published "
            "irreducible curves and 21 are explicit reducible divisors.  The "
            "primitive ample class D506+7H has square 90666 and is strict on "
            "all 31 imported walls; seven is minimal in this integral family."
        ),
        "scope": {
            "positive_gain": (
                "the known finite-volume H13/31 chamber now has an explicit "
                "ample-interior equation-level K3 realization and a complete "
                "effectivity/irreducibility decision for all (-2) mirrors"
            ),
            "not_full_nef": (
                "the imported 31-wall chamber is not asserted to equal the K3 "
                "nef chamber; five mirrors have square -4 and the complete K3 "
                "nef chamber has 15 curve walls"
            ),
            "minus_four": (
                "the five (-4,div2) vectors remain divisor reflections without "
                "effectivity or irreducible-curve claims"
            ),
            "minimality": (
                "k=7 is minimal only in the chosen family D506+kH for integral k>=0"
            ),
            "no_surface_automorphism_claim": True,
            "not_a_new_abstract_polytope": True,
            "worldwide_priority_claim": "none",
        },
        "sources": {
            key: {"path": str(path.relative_to(ROOT)), "file_sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        }
        | {
            "carrier_payload_sha256": EXPECTED_CARRIER_PAYLOAD_SHA256,
            "primary_reference": (
                "Xavier Roulleau, An atlas of K3 surfaces with finite "
                "automorphism group, EPIGA 6 (2022), article 19"
            ),
            "primary_doi": "10.46298/epiga.2022.6286",
            "primary_arxiv": "2003.08985",
            "pinned_primary_TeX_sha256": EXPECTED_ROULLEAU_TEX_SHA256,
            "checks": source_checks,
        },
        "published_nef_chamber_stock": {
            "lattice": "U+E8(-1)+D4(-1)",
            "coordinate_basis": "A1,...,A14",
            "gram": _matrix_payload(gram),
            "complete_irreducible_minus_two_curve_count": 15,
            "curve_coordinates_A1_through_A15": [
                _integer_vector_payload(curve) for curve in curves
            ],
            "curve_intersection_gram": _matrix_payload(wall_gram),
            "standard_fiber_basis_columns_in_Roulleau_basis": _matrix_payload(
                standard_to_roulleau
            ),
            "standard_fiber_basis_labels": [
                "F", "F+O", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8",
                "D1", "D2", "D3", "D4",
            ],
            "ample_class": {
                "label": "D506",
                "coordinates": _integer_vector_payload(D506),
                "square": 506,
                "degrees_on_A1_through_A15": ample_degrees,
            },
            "checks": curve_stock_checks,
        },
        "signed_marking": {
            "old_equation_fiber_interior_coordinates": _vector_payload(
                old_standard_interior
            ),
            "old_Roulleau_interior_coordinates": _vector_payload(
                old_roulleau_interior
            ),
            "old_Roulleau_interior_D506_degree": -15889,
            "chosen_marking": "negative of the geometric fiber-basis bridge",
            "signed_initial_coordinates": _vector_payload(signed_initial),
            "signed_initial_D506_degree": 15889,
            "signed_standard_to_Roulleau_matrix": _matrix_payload(
                signed_standard_to_roulleau
            ),
            "checks": marking_checks,
        },
        "explicit_weyl_transport": {
            "deterministic_rule": (
                "fire the most negative A_i pairing, breaking ties by the "
                "smallest one-based curve label"
            ),
            "word_length": len(word),
            "word_one_based": word,
            "firing_trace": firing_trace,
            "W_matrix": _matrix_payload(weyl),
            "W_determinant": int(weyl.det()),
            "combined_equation_fiber_to_nef_matrix": _matrix_payload(
                combined_standard_to_nef
            ),
            "combined_determinant": int(combined_standard_to_nef.det()),
            "final_nef_carrier_coordinates": _vector_payload(nef_interior),
            "final_nef_carrier_square": 385,
            "final_nef_carrier_D506_degree": 525,
            "final_pairings_with_A1_through_A15": nef_pairings,
            "primitive_integral_H_equals_2v_coordinates": _integer_vector_payload(H),
            "H_square": 1540,
            "H_D506_degree": 1050,
            "H_degrees_on_A1_through_A15": H_curve_degrees,
            "checks": weyl_checks,
        },
        "transported_mirrors": {
            "count": 31,
            "square_minus_two_effective_count": len(minus_two_records),
            "square_minus_two_irreducible_count": len(irreducible_records),
            "square_minus_two_reducible_count": len(reducible_records),
            "square_minus_two_effective_decompositions": minus_two_records,
            "square_minus_four_divisibility_two_count": len(minus_four_records),
            "square_minus_four_divisor_classes": minus_four_records,
            "effective_checks": effective_checks,
            "square_minus_four_checks": minus_four_checks,
            "checks": mirror_checks,
        },
        "explicit_simultaneous_ample_interior": {
            "family": "B_k = D506 + kH for nonnegative integral k",
            "ampleness_reason": (
                "D506 is ample and H is nef, hence every B_k is ample; the "
                "complete curve degrees are replayed explicitly for B_7"
            ),
            "D506_pairings_with_transported_facets_0_through_30": (
                D506_mirror_pairings
            ),
            "H_pairings_with_transported_facets_0_through_30": H_mirror_pairings,
            "nonpositive_imported_facets_for_k_0_through_6": minimum_failures,
            "B6_equals_D506_plus_6H": {
                "coordinates": _integer_vector_payload(B6),
                "square": _integer_dot(gram, B6, B6),
                "pairings_with_transported_facets_0_through_30": B6_mirror_pairings,
                "unique_nonpositive_imported_facet": {
                    "facet_zero_based": 5,
                    "pairing": 0,
                },
            },
            "B7_equals_D506_plus_7H": {
                "coordinates": _integer_vector_payload(B7),
                "primitive": _content(B7) == 1,
                "square": 90666,
                "D506_degree": 7856,
                "degrees_on_A1_through_A15": B7_curve_degrees,
                "pairings_with_transported_facets_0_through_30": B7_mirror_pairings,
                "minimum_complete_curve_degree": min(B7_curve_degrees),
                "minimum_imported_facet_pairing": min(B7_mirror_pairings),
                "minimum_facet_zero_based": B7_mirror_pairings.index(
                    min(B7_mirror_pairings)
                ),
            },
            "minimal_strict_k_in_chosen_family": first_strict_k,
            "checks": ample_checks,
        },
        "literal_Vinberg_gate": {
            "classification": vinberg["classification"],
            "dimension": vinberg["dimension"],
            "facet_count": vinberg["facets"],
            "finite_vertex_count": vinberg["finite_vertices"],
            "ideal_vertex_count": vinberg["cusps"],
            "elliptic_edge_count": vinberg["elliptic_edges"],
            "connected_affine_component_count": vinberg[
                "connected_affine_components"
            ],
            "bad_edge_count": vinberg["bad_edges"],
            "bad_parabolic_neighborhood_count": vinberg[
                "bad_parabolic_neighborhoods"
            ],
            "checks": vinberg_checks,
        },
        "all_checks_pass": True,
    }
    payload = _jsonify(payload)
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
            raise ArithmeticError("rebuilt A18-C06 nef/effectivity archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
