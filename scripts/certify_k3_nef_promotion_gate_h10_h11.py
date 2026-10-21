#!/usr/bin/env python3
"""Exact K3-nef promotion gate for the H10/16, H10/45, and H11/29 carriers.

The three archived chambers are genuine effective divisor/reflection
subsystems on K3 surfaces.  This certificate asks the stronger question:
can the displayed chamber itself be the full nef cone, or the full invariant
section of a nef cone, in the finite-index Picard carriers already found?

For H10/16 and every one of the ten admissible H10/45 overextensions, an
explicit ambient square-minus-two root cuts the strict interior of the
displayed fixed chamber.  For the full-rank H11/29 chamber, twelve genuine
facets have primitive ambient normal of square -4, so they cannot be K3 nef
walls.  The script also records the exact general root-hyperplane criterion
used by these obstructions.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
)
from scripts.certify_r12_h10_45_bh_divisor_sublattice_carrier import (
    _closure,
    _target_pair,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
H10_16_FOLD = (
    ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
)
H10_16_CARRIER = (
    ROOT / "data/enumerations/r12_h10_16_order4_divisor_overlattice_carrier.json"
)
H10_45_FOLD = (
    ROOT
    / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
)
H10_45_CARRIER = (
    ROOT / "data/enumerations/r12_h10_45_bh_divisor_sublattice_carrier.json"
)
H11_29_CARRIER = (
    ROOT / "data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json"
)
OUTPUT = ROOT / "data/enumerations/k3_nef_promotion_gate_h10_h11.json"

EXPECTED_SHA256 = {
    "h10_16_fold": "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    "h10_16_carrier": "5dbf9fabab666864baf89d9b44cb89addaca5beabe80422ad696f04743e32e8b",
    "h10_45_fold": "448bcc47937a2763bbd4ea97466fbc48924dc862e88ba8719e5ab98b0cc4c916",
    "h10_45_carrier": "90ef1da5c34185143bed0a9341dfa565248d55ba1df8d4a74a12bbfc4886d9a8",
    "h11_29_carrier": "5d4c2299ecf2b8b853063c142472c46c61f8554197168b06d144924cdc6d07fd",
}

H10_16_EXTRA_ROOT = [-1, -1, -1, -1, -1, 0, 0, 1, 1, 0, 1, -1]

# One short exact witness in each overextension, aligned with the canonical
# ordering (subgroup order, then lexicographically sorted subgroup elements).
H10_45_EXTRA_ROOTS = [
    [0, 0, 0, -1, -1, -1, 1, -1, -1, -1, 0, -1],
    [-1, 0, 0, 0, -1, 0, 1, 1, 0, 1, 1, -1],
    [0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, -1, 0, -1, 0, 1],
    [0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 0, 0, -1, 1, 0, 1, 0, 0, 1],
    [0, 1, 1, 1, -1, 0, 0, 0, 1, 0, 0, 1],
    [0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, -1],
    [0, -1, -1, 0, -1, 0, 1, 1, 0, 1, 1, -1],
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atom(value: Any) -> int | str:
    rational = sp.Rational(value)
    return int(rational) if rational.q == 1 else str(rational)


def _vector_payload(vector: sp.MatrixBase | Iterable[Any]) -> list[int | str]:
    return [_atom(value) for value in vector]


def _matrix_payload(matrix: sp.MatrixBase) -> list[list[int | str]]:
    return [
        [_atom(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    raw = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _primitive_integer(vector: sp.MatrixBase) -> bool:
    values = [abs(int(value)) for value in vector]
    return all(sp.Rational(value).q == 1 for value in vector) and math.gcd(*values) == 1


def _column_orbit(
    vector: sp.Matrix, action_column: sp.Matrix, maximum: int = 16
) -> list[sp.Matrix]:
    orbit: list[sp.Matrix] = []
    current = vector
    while not any(current == old for old in orbit):
        orbit.append(current)
        current = action_column * current
        if len(orbit) > maximum:
            raise ArithmeticError("root orbit exceeded the certified action order")
    if current != orbit[0]:
        raise ArithmeticError("root orbit did not close at its starting vector")
    return orbit


def _orbit_sum_rows(
    wall_rows: sp.Matrix, cycles_one_based: Sequence[Sequence[int]]
) -> sp.Matrix:
    return sp.Matrix.vstack(
        *(
            sum(
                (wall_rows.row(index - 1) for index in cycle),
                sp.zeros(1, wall_rows.cols),
            )
            for cycle in cycles_one_based
        )
    )


def _point_in_fixed_basis(
    coordinates: sp.Matrix,
    basis_rows: sp.Matrix,
    gram: sp.Matrix,
    wall_rows: sp.Matrix,
    root: sp.Matrix,
) -> dict[str, Any]:
    source = basis_rows.T * coordinates
    return {
        "coordinates": _vector_payload(coordinates),
        "square": _atom((source.T * gram * source)[0]),
        "extra_root_pairing": _atom((source.T * gram * root)[0]),
        "minimum_old_wall_pairing": _atom(min(wall_rows * gram * source)),
        "strict_for_all_old_walls": all(
            value > 0 for value in wall_rows * gram * source
        ),
        "timelike": bool((source.T * gram * source)[0] > 0),
    }


def _h10_16_gate(fold: dict[str, Any], carrier: dict[str, Any]) -> dict[str, Any]:
    stock = next(row for row in fold["stocks"] if row["stock_key"] == "R12_UD8A1x2")
    survivor = next(row for row in stock["survivors"] if row["class_index"] == 3)
    class_record = next(row for row in stock["classes"] if row["class_index"] == 3)
    geometry = survivor["geometry"]
    gram = sp.Matrix(stock["basis_gram"])
    wall_rows = sp.Matrix(stock["wall_coordinate_rows"])
    saturation_basis = sp.Matrix(
        carrier["unique_saturation"]["basis_rows_in_archived_M_coordinates"]
    )
    ambient_gram = sp.simplify(saturation_basis * gram * saturation_basis.T)

    permutation = class_record["representative_generators_one_based"][0]
    action_column = sp.zeros(12)
    for source_index in range(12):
        action_column[:, source_index] = wall_rows[
            permutation[source_index] - 1, :
        ].T

    cycles = geometry["orbit_cycles"]
    orbit_sums = _orbit_sum_rows(wall_rows, cycles)
    basis_indices = [
        cycles.index(orbit) for orbit in geometry["vertex_basis_orbits"]
    ]
    fixed_basis = orbit_sums[basis_indices, :]
    fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
    archived_orbit_gram = sp.Matrix(geometry["orbit_gram"])

    ambient_coordinates = sp.Matrix(H10_16_EXTRA_ROOT)
    root = saturation_basis.T * ambient_coordinates
    root_image = action_column * root
    orbit = _column_orbit(root, action_column)
    vertex_vectors = [
        sp.Matrix(row["basis_coordinates"]) for row in geometry["vertices"]
    ]
    root_values = [
        (vertex.T * fixed_basis * gram * root)[0] for vertex in vertex_vectors
    ]
    interior = sum(vertex_vectors, sp.zeros(fixed_basis.rows, 1))
    positive_ray_index = 44
    negative_ray_index = 17
    positive_point = interior + vertex_vectors[positive_ray_index]
    negative_point = interior + vertex_vectors[negative_ray_index]
    zero_point = positive_point + negative_point

    points = {
        "archived_extreme_ray_sum": _point_in_fixed_basis(
            interior, fixed_basis, gram, wall_rows, root
        ),
        "strict_positive_point": _point_in_fixed_basis(
            positive_point, fixed_basis, gram, wall_rows, root
        ),
        "strict_negative_point": _point_in_fixed_basis(
            negative_point, fixed_basis, gram, wall_rows, root
        ),
        "strict_zero_point": _point_in_fixed_basis(
            zero_point, fixed_basis, gram, wall_rows, root
        ),
    }
    expected_histogram = {
        "-7": 1,
        "-5": 3,
        "-4": 2,
        "-3": 4,
        "-2": 5,
        "-1": 8,
        "0": 30,
        "1": 8,
        "2": 5,
        "3": 4,
        "4": 2,
        "5": 3,
        "7": 1,
    }
    histogram = {
        str(value): count
        for value, count in sorted(Counter(root_values).items())
    }
    checks = {
        "upstream_certificates_pass": fold["all_checks_pass"] and carrier["all_checks_pass"],
        "unique_admissible_nontrivial_saturation": carrier[
            "inverse_overextension_census"
        ]["nontrivial_action_stable_even_overgroup_count"]
        == 1,
        "ambient_is_even_lorentzian": (
            exact_inertia(ambient_gram) == (1, 11, 0)
            and all(ambient_gram[index, index] % 2 == 0 for index in range(12))
        ),
        "ambient_root_is_primitive_minus_two": (
            _primitive_integer(ambient_coordinates)
            and (ambient_coordinates.T * ambient_gram * ambient_coordinates)[0]
            == -2
        ),
        "root_is_new_outside_M": any(value.q != 1 for value in root),
        "action_orbit_is_two_orthogonal_roots": (
            len(orbit) == 2 and (root.T * gram * root_image)[0] == 0
        ),
        "orbit_restrictions_are_identical": (
            fixed_basis * gram * (root_image - root)
            == sp.zeros(fixed_basis.rows, 1)
        ),
        "fixed_basis_gram_replays": fixed_gram
        == archived_orbit_gram.extract(basis_indices, basis_indices),
        "all_76_extreme_rays_replayed": len(vertex_vectors) == 76,
        "extreme_ray_pairing_histogram_replays": histogram == expected_histogram,
        "chosen_boundary_rays_have_opposite_pairing": (
            root_values[positive_ray_index] == 1
            and root_values[negative_ray_index] == -1
        ),
        "strict_points_are_timelike_and_inside_old_chamber": all(
            row["timelike"] and row["strict_for_all_old_walls"]
            for row in points.values()
        ),
        "strict_zero_point_lies_on_extra_root_wall": points[
            "strict_zero_point"
        ]["extra_root_pairing"]
        == 0,
        "opposite_strict_sides_certified": (
            points["strict_positive_point"]["extra_root_pairing"] > 0
            and points["strict_negative_point"]["extra_root_pairing"] < 0
        ),
        "witness_points_are_fixed": all(
            action_column * fixed_basis.T * vector == fixed_basis.T * vector
            for vector in (interior, positive_point, negative_point, zero_point)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H10/16 promotion gate changed: {checks}")

    return {
        "case_id": "R12-D8A1x2-C2b-H10-16",
        "dimension_facets": [10, 16],
        "ambient_NS": "U+E8(-1)+A1(-1)^2",
        "ambient_index_over_divisor_lattice": 2,
        "admissible_nontrivial_saturation_count": 1,
        "extra_ambient_root": {
            "coordinates_in_saturation_basis": H10_16_EXTRA_ROOT,
            "coordinates_in_M_tensor_Q": _vector_payload(root),
            "square": -2,
            "primitive": True,
            "lies_in_M": False,
            "action_orbit_size": len(orbit),
            "inner_product_with_nontrivial_image": _atom(
                (root.T * gram * root_image)[0]
            ),
            "same_restriction_on_fixed_space": True,
        },
        "extreme_ray_test": {
            "ray_count": len(vertex_vectors),
            "pairing_histogram": histogram,
            "positive_ray": {
                "zero_based_index": positive_ray_index,
                "one_based_index": positive_ray_index + 1,
                "coordinates": _vector_payload(vertex_vectors[positive_ray_index]),
                "pairing": _atom(root_values[positive_ray_index]),
            },
            "negative_ray": {
                "zero_based_index": negative_ray_index,
                "one_based_index": negative_ray_index + 1,
                "coordinates": _vector_payload(vertex_vectors[negative_ray_index]),
                "pairing": _atom(root_values[negative_ray_index]),
            },
        },
        "strict_cut_certificate": points,
        "decision": "not_full_invariant_nef_section_in_the_unique_rank12_saturation",
        "reason": (
            "the displayed ambient (-2)-root hyperplane contains a strict "
            "timelike point of the H10/16 chamber and has strict chamber "
            "points on both sides"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def _h10_45_subgroups(target: dict[str, Any]) -> list[frozenset[tuple[int, ...]]]:
    gram = target["fixed_gram"]
    moduli, _, discriminant_gram = _discriminant_data(gram)
    elements = list(itertools.product(*(range(modulus) for modulus in moduli)))
    gamma = next(
        fixed_class
        for fixed_class, coinvariant_class in target["glue"]
        if any(coinvariant_class)
    )

    def q(element: Sequence[int]) -> sp.Rational:
        vector = sp.Matrix(element)
        return sp.Rational((vector.T * discriminant_gram * vector)[0]) % 2

    def bilinear(left: Sequence[int], right: Sequence[int]) -> sp.Rational:
        return sp.Rational(
            (sp.Matrix(left).T * discriminant_gram * sp.Matrix(right))[0]
        ) % 1

    isotropic_perpendicular = [
        element
        for element in elements
        if q(element) == 0 and bilinear(element, gamma) == 0
    ]
    subgroups: set[frozenset[tuple[int, ...]]] = set()
    for element in isotropic_perpendicular:
        subgroup = _closure([element], moduli)
        if len(subgroup) in (2, 4) and all(q(value) == 0 for value in subgroup):
            subgroups.add(subgroup)
    for left, right in itertools.combinations(isotropic_perpendicular, 2):
        subgroup = _closure([left, right], moduli)
        if (
            len(subgroup) == 4
            and all(q(value) == 0 for value in subgroup)
            and all(bilinear(x, y) == 0 for x in subgroup for y in subgroup)
        ):
            subgroups.add(subgroup)
    return sorted(subgroups, key=lambda value: (len(value), sorted(value)))


def _extension_basis(
    subgroup: frozenset[tuple[int, ...]], target: dict[str, Any]
) -> sp.Matrix:
    fixed_gram = target["fixed_gram"]
    _, left_inverse, _ = _discriminant_data(fixed_gram)
    rows = [sp.eye(12).row(index) for index in range(12)]
    for element in subgroup:
        if not any(element):
            continue
        fixed_coordinates = fixed_gram.inv() * left_inverse * sp.Matrix(element)
        rows.append(fixed_coordinates.T * target["fixed"])
    denominator = math.lcm(
        *(sp.Rational(value).q for row in rows for value in row)
    )
    integral_generators = sp.Matrix.vstack(
        *(denominator * row for row in rows)
    ).applyfunc(int)
    return hermite_normal_form(integral_generators.T).T / denominator


def _source_point(
    point: sp.Matrix,
    gram: sp.Matrix,
    walls: sp.Matrix,
    root: sp.Matrix,
    action: sp.Matrix,
) -> dict[str, Any]:
    wall_pairings = point * gram * walls.T
    return {
        "coordinates_in_M_tensor_Q": _vector_payload(point),
        "square": _atom((point * gram * point.T)[0]),
        "extra_root_pairing": _atom((point * gram * root)[0]),
        "minimum_old_wall_pairing": _atom(min(wall_pairings)),
        "strict_for_all_59_old_walls": all(value > 0 for value in wall_pairings),
        "timelike": bool((point * gram * point.T)[0] > 0),
        "fixed_by_action": point * action == point,
    }


def _h10_45_gate(fold: dict[str, Any], carrier: dict[str, Any]) -> dict[str, Any]:
    target = _target_pair(fold)
    gram = target["gram"]
    action = target["action"]
    wall_rows = sp.Matrix(fold["coordinate_reconstruction"]["wall_coordinate_rows"])
    case = next(
        row
        for row in fold["cases"]
        if row["id"] == "R12-D4x2A1x2-S6transposition-H10-45"
    )
    wall_orbits = case["subgroup"]["wall_orbits_one_based"]
    orbit_roots = _orbit_sum_rows(wall_rows, wall_orbits)
    interior = sp.Matrix([[3, 3, 3, 0, 0, 2, 2, 6, 6, 3, 1, 7]])

    facet_witnesses: list[sp.Matrix] = []
    for orbit_root in orbit_roots.tolist():
        root_row = sp.Matrix([orbit_root])
        projection = sp.simplify(
            interior
            - (interior * gram * root_row.T)[0]
            / (root_row * gram * root_row.T)[0]
            * root_row
        )
        facet_witnesses.append(projection)
    archived_witnesses = case["strict_facet_witnesses"]["records"]
    # The source routine averages over the order-two subgroup, hence its
    # archived representative is exactly twice our projectively equivalent
    # witness based on the already invariant vector ``interior``.
    witness_replay = all(
        str(((2 * point) * gram * (2 * point).T)[0]) == archived["square"]
        and str(min(
            value
            for index, value in enumerate((2 * point) * gram * orbit_roots.T)
            if index != position
        ))
        == archived["minimum_other_pairing"]
        and (point * gram * orbit_roots.row(position).T)[0] == 0
        for position, (point, archived) in enumerate(
            zip(facet_witnesses, archived_witnesses)
        )
    )

    subgroups = _h10_45_subgroups(target)
    records: list[dict[str, Any]] = []
    for extension_index, (subgroup, root_coordinates) in enumerate(
        zip(subgroups, H10_45_EXTRA_ROOTS)
    ):
        basis = _extension_basis(subgroup, target)
        ambient_gram = sp.simplify(basis * gram * basis.T)
        ambient_action = sp.simplify(basis * action * basis.inv())
        ambient_root = sp.Matrix(root_coordinates)
        root = basis.T * ambient_root
        root_image = action.T * root
        orbit = _column_orbit(root, action.T)
        facet_values = [(point * gram * root)[0] for point in facet_witnesses]
        positive_index = next(
            index for index, value in enumerate(facet_values) if value > 0
        )
        negative_index = next(
            index for index, value in enumerate(facet_values) if value < 0
        )
        interior_pairing = (interior * gram * root)[0]
        positive_scale = next(
            scale
            for scale in range(1, 100)
            if interior_pairing + scale * facet_values[positive_index] > 0
        )
        negative_scale = next(
            scale
            for scale in range(1, 100)
            if interior_pairing + scale * facet_values[negative_index] < 0
        )
        positive_point = interior + positive_scale * facet_witnesses[positive_index]
        negative_point = interior + negative_scale * facet_witnesses[negative_index]
        positive_pairing = (positive_point * gram * root)[0]
        negative_pairing = (negative_point * gram * root)[0]
        zero_point = sp.simplify(
            -negative_pairing * positive_point + positive_pairing * negative_point
        )
        points = {
            "strict_positive_point": _source_point(
                positive_point, gram, wall_rows, root, action
            ),
            "strict_negative_point": _source_point(
                negative_point, gram, wall_rows, root, action
            ),
            "strict_zero_point": _source_point(
                zero_point, gram, wall_rows, root, action
            ),
        }
        expected_determinant = -16 if len(subgroup) == 2 else -4
        checks = {
            "extension_basis_has_expected_index": abs(basis.det())
            == sp.Rational(1, len(subgroup)),
            "ambient_lattice_is_even_lorentzian": (
                ambient_gram.det() == expected_determinant
                and exact_inertia(ambient_gram) == (1, 11, 0)
                and all(value.q == 1 for value in ambient_gram)
                and all(ambient_gram[index, index] % 2 == 0 for index in range(12))
            ),
            "action_extends_integrally": (
                all(value.q == 1 for value in ambient_action)
                and ambient_action * ambient_gram * ambient_action.T == ambient_gram
                and ambient_action**2 == sp.eye(12)
            ),
            "ambient_root_is_primitive_minus_two": (
                _primitive_integer(ambient_root)
                and (ambient_root.T * ambient_gram * ambient_root)[0] == -2
            ),
            "root_is_new_outside_M": any(value.q != 1 for value in root),
            "root_orbit_has_size_one_or_two": len(orbit) in (1, 2),
            "orbit_restrictions_are_identical": (
                target["fixed"] * gram * (root_image - root)
                == sp.zeros(target["fixed"].rows, 1)
            ),
            "facet_witness_pairings_have_both_signs": (
                bool(min(facet_values) < 0 < max(facet_values))
            ),
            "strict_points_are_fixed_timelike_and_inside": all(
                row["fixed_by_action"]
                and row["timelike"]
                and row["strict_for_all_59_old_walls"]
                for row in points.values()
            ),
            "opposite_strict_sides_certified": (
                bool(positive_pairing > 0 and negative_pairing < 0)
            ),
            "strict_zero_point_lies_on_extra_root_wall": (
                bool((zero_point * gram * root)[0] == 0)
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"H10/45 extension {extension_index} gate changed: {checks}"
            )
        records.append(
            {
                "extension_index_zero_based": extension_index,
                "subgroup_order_and_ambient_index": len(subgroup),
                "subgroup_elements_in_fixed_smith_coordinates": [
                    list(element) for element in sorted(subgroup)
                ],
                "ambient_NS_isometry_class": (
                    "U+D8(-1)+A1(-1)^2"
                    if len(subgroup) == 2
                    else "U+E8(-1)+A1(-1)^2"
                ),
                "ambient_determinant": expected_determinant,
                "extra_ambient_root": {
                    "coordinates_in_extension_basis": root_coordinates,
                    "coordinates_in_M_tensor_Q": _vector_payload(root),
                    "square": -2,
                    "primitive": True,
                    "lies_in_M": False,
                    "action_orbit_size": len(orbit),
                    "inner_product_with_action_image": _atom(
                        (root.T * gram * root_image)[0]
                    ),
                    "same_restriction_on_fixed_space": True,
                },
                "facet_witness_pairing_range": [
                    _atom(min(facet_values)),
                    _atom(max(facet_values)),
                ],
                "chosen_facet_witnesses": {
                    "positive_one_based": positive_index + 1,
                    "negative_one_based": negative_index + 1,
                    "positive_scale_from_interior": positive_scale,
                    "negative_scale_from_interior": negative_scale,
                },
                "strict_cut_certificate": points,
                "decision": "extra_root_strictly_refines_H10_45",
                "checks": checks,
                "all_checks_pass": True,
            }
        )

    counts = Counter(len(subgroup) for subgroup in subgroups)
    checks = {
        "upstream_certificates_pass": fold["all_checks_pass"] and carrier["all_checks_pass"],
        "archived_45_strict_facet_witnesses_replayed": witness_replay,
        "complete_admissible_extension_count": len(subgroups) == 10,
        "seven_index_two_extensions": counts == {2: 7, 4: 3},
        "one_cutting_root_for_every_extension": (
            len(records) == 10 and all(row["all_checks_pass"] for row in records)
        ),
        "upstream_census_agrees": (
            carrier["finite_index_carrier_criterion"][
                "admissible_order_two_subgroups"
            ]
            == 7
            and carrier["finite_index_carrier_criterion"][
                "admissible_order_four_subgroups"
            ]
            == 3
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H10/45 complete promotion gate changed: {checks}")
    return {
        "case_id": "R12-D4x2A1x2-S6transposition-H10-45",
        "dimension_facets": [10, 45],
        "admissible_rank12_action_compatible_overextensions": 10,
        "index_histogram": {str(order): count for order, count in sorted(counts.items())},
        "criterion_scope": (
            "the complete order-two/order-four even isotropic fixed-discriminant "
            "overextension census compatible with the eigen-glue"
        ),
        "extensions": records,
        "decision": "no_full_invariant_nef_section_in_any_admissible_rank12_overextension",
        "reason": (
            "each of the ten extensions contains an explicit new (-2)-root "
            "whose restricted hyperplane cuts the strict H10/45 interior"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def _h11_29_gate(carrier: dict[str, Any]) -> dict[str, Any]:
    facets = carrier["lattice_and_embedding"]["facet_reflection_extensions"]
    records = facets["records"]
    minus_four = [row for row in records if row["primitive_NS_square"] == -4]
    indices = [row["facet_zero_based"] + 1 for row in minus_four]
    checks = {
        "upstream_certificate_passes": carrier["all_checks_pass"],
        "full_rank_lorentzian_chamber": (
            carrier["vinberg_finite_volume"]["dimension"] == 11
            and carrier["lattice_and_embedding"]["overclosure_N"]["rank"] == 12
        ),
        "all_29_facets_are_present": facets["count"] == 29,
        "primitive_square_histogram_replays": facets[
            "primitive_NS_square_histogram"
        ]
        == {"-2": 17, "-4": 12},
        "twelve_obstructing_facets_found": len(minus_four) == 12,
        "all_obstructing_normals_are_primitive_divisibility_two": all(
            row["primitive_NS_divisibility"] == 2
            and row["checks"]["primitive_in_NS"]
            and row["checks"]["same_hyperplane"]
            for row in minus_four
        ),
        "listed_indices_replay": indices
        == [1, 2, 7, 8, 9, 10, 13, 19, 21, 25, 28, 29],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H11/29 facet-normal gate changed: {checks}")
    return {
        "case_id": "A16_6_1_C01-H11-29-UD8-overclosure",
        "dimension_facets": [11, 29],
        "ambient_NS": "U+D8(-1)+A1(-1)^2",
        "full_rank_in_ambient_positive_cone": True,
        "primitive_facet_normal_square_histogram": {"-2": 17, "-4": 12},
        "obstructing_square_minus_four_facets_one_based": indices,
        "hyperplane_argument": (
            "A codimension-one facet of a K3 nef cone must be supported by "
            "delta-perp for an integral root delta with delta^2=-2.  If such "
            "a hyperplane equalled n-perp for one of these primitive n with "
            "n^2=-4, then delta=c*n for rational c and c^2=1/2, impossible."
        ),
        "decision": "not_a_full_nef_cone_for_the_rank12_ambient_NS",
        "marking_independent": True,
        "reason": (
            "twelve genuine facets have primitive ambient square -4 normals "
            "and therefore are not K3 Weyl-root walls"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "h10_16_fold": H10_16_FOLD,
        "h10_16_carrier": H10_16_CARRIER,
        "h10_45_fold": H10_45_FOLD,
        "h10_45_carrier": H10_45_CARRIER,
        "h11_29_carrier": H11_29_CARRIER,
    }
    hashes = {key: _sha256(path) for key, path in paths.items()}
    if hashes != EXPECTED_SHA256:
        raise ArithmeticError(f"a pinned promotion source changed: {hashes}")
    archives = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in paths.items()
    }

    h10_16 = _h10_16_gate(
        archives["h10_16_fold"], archives["h10_16_carrier"]
    )
    h10_45 = _h10_45_gate(
        archives["h10_45_fold"], archives["h10_45_carrier"]
    )
    h11_29 = _h11_29_gate(archives["h11_29_carrier"])
    checks = {
        "all_three_case_gates_pass": all(
            row["all_checks_pass"] for row in (h10_16, h10_45, h11_29)
        ),
        "h10_16_unique_saturation_blocked": h10_16["decision"].startswith("not_full"),
        "h10_45_all_ten_overextensions_blocked": (
            h10_45["admissible_rank12_action_compatible_overextensions"] == 10
            and h10_45["decision"].startswith("no_full")
        ),
        "h11_29_full_rank_facet_norm_obstruction": h11_29["decision"].startswith(
            "not_a_full"
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"top-level promotion gate changed: {checks}")

    certificate: dict[str, Any] = {
        "schema": "k3_nef_promotion_gate_h10_h11_v1",
        "theorem": (
            "Let M be a finite-index g-stable divisor lattice in a K3 Picard "
            "lattice N, and let C be the full-dimensional pointed fixed "
            "section of a complete M-root chamber.  Assuming the displayed "
            "old facets are genuine N-roots "
            "and no N-root is orthogonal to the fixed space, C is the fixed "
            "section of an N-Weyl chamber if and only if no N-root hyperplane "
            "meets the relative interior of C.  For a finitely generated "
            "polyhedral C this is equivalent to every N-root restriction "
            "having one weak sign, but not vanishing identically, on all "
            "extreme rays.  The unique admissible H10/16 saturation and all ten "
            "admissible H10/45 overextensions fail by explicit (-2)-root cuts; "
            "the full-rank H11/29 carrier fails because twelve genuine facets "
            "have primitive ambient square -4 normals."
        ),
        "source_sha256": hashes,
        "promotion_criterion": {
            "setting": (
                "M is a finite-index g-stable even hyperbolic sublattice of "
                "the K3 Picard lattice N; F=Fix(g)_R; C is a full-dimensional "
                "pointed polyhedral cone in one component of the positive cone of F; "
                "each old facet is a genuine restriction of an N (-2)-root wall"
            ),
            "root_orthogonality_gate": (
                "there is no delta in N with delta^2=-2 and delta perpendicular "
                "to all of F; equivalently an invariant ample class is not "
                "forced onto an ambient root wall"
            ),
            "exact_if_and_only_if": (
                "C equals the closure of an ambient N-Weyl chamber intersected "
                "with F if and only if delta-perp misses the relative interior "
                "of C for every ambient root delta"
            ),
            "extreme_ray_form": (
                "if v_1,...,v_k are all extreme rays, then for every ambient "
                "root delta the values (v_i,delta) are all >=0 or all <=0, "
                "and are not all zero"
            ),
            "finite_failure_witness": (
                "two strict timelike points of C with opposite delta-pairings, "
                "or equivalently one strict timelike point on delta-perp plus "
                "strict points on both sides, proves strict refinement"
            ),
            "full_rank_facet_corollary": (
                "when F=N_R, every genuine facet hyperplane must have a "
                "primitive normal proportional over Q to a (-2)-root; a "
                "primitive square -4 normal can never qualify"
            ),
            "sufficiency_reason": (
                "choose a generic point of C and its ambient Weyl chamber.  "
                "The no-cut condition keeps all of C in its closure, while "
                "the old N-root facet inequalities prevent that chamber's "
                "fixed section from extending beyond C"
            ),
        },
        "cases": {
            "H10_16": h10_16,
            "H10_45": h10_45,
            "H11_29": h11_29,
        },
        "positive_context_retained": {
            "divisor_subsystem_realizations_remain_valid": True,
            "finite_volume_certificates_remain_valid": True,
            "ambient_complete_nef_replacements": {
                "U+E8(-1)+A1(-1)^2": "the already certified H10/12 fixed complete-nef section",
                "U+D8(-1)+A1(-1)^2": "the already certified H10/15 fixed complete-nef section",
            },
            "interpretation": (
                "the result narrows the positive claim rather than retracting "
                "it: H10/16, H10/45, and H11/29 remain effective algebraic "
                "Coxeter divisor/reflection carriers, but the displayed coarse "
                "chambers are not the full nef objects in these rank-12 carriers"
            ),
        },
        "scope": {
            "proved": (
                "no promotion for the unique admissible rank-12 H10/16 "
                "saturation, any of the ten admissible rank-12 H10/45 finite "
                "overextensions, or the existing rank-12 full-rank H11/29 carrier"
            ),
            "not_proved": (
                "a worldwide impossibility for every unrelated embedding into "
                "a higher-Picard K3 lattice, or nonexistence of the abstract "
                "hyperbolic Coxeter polytopes"
            ),
            "next_positive_design_rule": (
                "search embeddings whose primitive facet normals are all (-2) "
                "and then run the ambient-root one-sign test before investing "
                "in a geometric family"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }
    certificate["payload_sha256"] = _payload_digest(certificate)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if certificate != archived:
            raise SystemExit("recomputed K3 nef-promotion gate differs from archive")
        print(f"verified {args.output}")
        return
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.output}")
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
