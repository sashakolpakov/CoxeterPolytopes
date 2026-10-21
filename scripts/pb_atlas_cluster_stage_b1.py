#!/usr/bin/env python3
"""Compile one PB-Atlas Stage-A endpoint program into numerical Stage B1.

The compiler consumes only the arrays and integer indices in one frozen
Stage-A ``endpointProgram``.  It reverses the exact contraction witness,
recovers the Enriques proximity DAG from strict-versus-total exceptional
coordinates, and tabulates every supplied wall as a divisor on ``P2`` or
``F_e``.  It also emits Riemann--Roch data and fully indexed numerical
templates for the ``-K`` and ``-2K`` evaluation matrices.

No surface/type name, archived geometric point, incidence equation, or
known classification label enters this module.  Matrix entries depending
on universal point/blowup-chart coordinates deliberately belong to Stage
B2; Stage B1 fixes their section columns and jet rows but does not claim a
rank, a nonempty incidence stratum, irreducibility, a bad-factor list,
same-P return, or L4.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


INPUT_SCHEMA = "coxeter.pb_atlas.numerical_contraction_stage_a.v1#endpointProgram"
OUTPUT_SCHEMA = "coxeter.pb_atlas.numerical_cluster_stage_b1.v1"

_PROGRAM_KEYS = {
    "endpoint_kind",
    "hirzebruch_e",
    "hirzebruch_minimal",
    "contracted_wall_indices_zero_based",
    "sequential_contraction_order_zero_based",
    "strict_contracted_wall_gram",
    "model_basis_in_input_picard_coordinates",
    "total_exceptional_basis_in_input_picard_coordinates",
    "endpoint_blowup_basis_gram",
    "canonical_coordinates_in_endpoint_blowup_basis",
    "nef_rays_in_model_basis",
    "all_wall_coordinates_in_endpoint_blowup_basis",
    "all_wall_endpoint_image_coordinates",
    "boundary",
    "contraction_trace",
    "permutation_invariant_profile_sha256",
    "stage_B_payload_complete",
}


def _json_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _integer_vector(value: Any, name: str) -> list[int]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be an integer array")
    return [_integer(entry, name) for entry in value]


def _integer_matrix(value: Any, name: str) -> list[list[int]]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be an integer matrix")
    rows = [_integer_vector(row, name) for row in value]
    if rows and any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{name} is ragged")
    return rows


def _model_gram(kind: str, e: int | None) -> list[list[int]]:
    if kind == "P2":
        return [[1]]
    assert e is not None
    return [[-e, 1], [1, 0]]


def _block_gram(kind: str, e: int | None, exceptional_count: int) -> list[list[int]]:
    model = _model_gram(kind, e)
    size = len(model) + exceptional_count
    result = [[0] * size for _ in range(size)]
    for left, row in enumerate(model):
        for right, value in enumerate(row):
            result[left][right] = value
    for index in range(exceptional_count):
        result[len(model) + index][len(model) + index] = -1
    return result


def _pair(
    left: Sequence[int], right: Sequence[int], kind: str, e: int | None
) -> int:
    model_rank = 1 if kind == "P2" else 2
    if len(left) != len(right):
        raise ValueError("intersection vectors have different lengths")
    if kind == "P2":
        answer = int(left[0]) * int(right[0])
    else:
        assert e is not None
        answer = (
            -e * int(left[0]) * int(right[0])
            + int(left[0]) * int(right[1])
            + int(left[1]) * int(right[0])
        )
    return answer - sum(
        int(left[index]) * int(right[index])
        for index in range(model_rank, len(left))
    )


def _transpose(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    if not matrix:
        return []
    return [list(column) for column in zip(*matrix)]


def _matmul(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> list[list[int]]:
    if not left:
        return []
    if not right:
        return [[] for _ in left]
    if len(left[0]) != len(right):
        raise ValueError("matrix dimensions do not match")
    return [
        [
            sum(int(left[i][k]) * int(right[k][j]) for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _compress_colours(values: Sequence[Any]) -> list[int]:
    encoded = [json.dumps(value, sort_keys=True, separators=(",", ":")) for value in values]
    table = {value: index for index, value in enumerate(sorted(set(encoded)))}
    return [table[value] for value in encoded]


def _refined_unlabelled_profile(
    *,
    endpoint_kind: str,
    hirzebruch_e: int | None,
    centers: Sequence[Mapping[str, Any]],
    divisor_rows: Sequence[Mapping[str, Any]],
    swap_f0_rulings: bool = False,
) -> dict[str, Any]:
    """Return a relabelling-invariant colour-refined cluster profile.

    The profile is intentionally an invariant, not a claim of a complete
    graph-isomorphism canonical form.  Exact program data and a chronological
    witness remain in the packet; mutation tests separately construct the
    explicit relabelling isomorphism.
    """

    center_count = len(centers)
    wall_count = len(divisor_rows)
    parents = [tuple(int(value) for value in row["parent_center_positions_zero_based"]) for row in centers]
    children = [tuple(int(value) for value in row["child_center_positions_zero_based"]) for row in centers]
    depths = [int(row["dependency_depth"]) for row in centers]
    coefficients = [
        tuple(int(value) for value in row["exceptional_coefficients_chronological"])
        for row in divisor_rows
    ]

    wall_bases = []
    for row in divisor_rows:
        endpoint = list(row["endpoint_class_coordinates"])
        if swap_f0_rulings:
            endpoint[0], endpoint[1] = endpoint[1], endpoint[0]
        wall_bases.append(
            (
                tuple(endpoint),
                int(row["anti_bicanonical_boundary_weight"]),
                bool(row["is_contracted_wall"]),
                int(row["self_intersection"]),
                int(row["canonical_intersection"]),
            )
        )

    center_colours = _compress_colours(
        [
            (
                len(parents[index]),
                len(children[index]),
                depths[index],
                str(centers[index]["center_kind"]),
            )
            for index in range(center_count)
        ]
    )
    wall_colours = _compress_colours(wall_bases)
    for _iteration in range(center_count + wall_count + 2):
        center_signatures = []
        for center in range(center_count):
            incident = sorted(
                (coefficients[wall][center], wall_colours[wall])
                for wall in range(wall_count)
                if coefficients[wall][center]
            )
            center_signatures.append(
                (
                    center_colours[center],
                    tuple(sorted(center_colours[parent] for parent in parents[center])),
                    tuple(sorted(center_colours[child] for child in children[center])),
                    tuple(incident),
                )
            )
        wall_signatures = [
            (
                wall_colours[wall],
                tuple(
                    sorted(
                        (coefficients[wall][center], center_colours[center])
                        for center in range(center_count)
                        if coefficients[wall][center]
                    )
                ),
            )
            for wall in range(wall_count)
        ]
        new_center_colours = _compress_colours(center_signatures)
        new_wall_colours = _compress_colours(wall_signatures)
        if new_center_colours == center_colours and new_wall_colours == wall_colours:
            break
        center_colours = new_center_colours
        wall_colours = new_wall_colours

    center_histogram = Counter(center_colours)
    edge_colours = sorted(
        (center_colours[child], center_colours[parent])
        for child in range(center_count)
        for parent in parents[child]
    )
    wall_records = []
    for wall in range(wall_count):
        wall_records.append(
            (
                wall_bases[wall],
                tuple(
                    sorted(
                        (coefficients[wall][center], center_colours[center])
                        for center in range(center_count)
                        if coefficients[wall][center]
                    )
                ),
            )
        )
    layer_records = []
    for depth in sorted(set(depths)):
        layer_records.append(
            (
                depth,
                tuple(sorted(center_colours[index] for index in range(center_count) if depths[index] == depth)),
            )
        )
    return {
        "endpoint_kind": endpoint_kind,
        "hirzebruch_e": hirzebruch_e,
        "center_colour_histogram": [
            [int(colour), int(count)] for colour, count in sorted(center_histogram.items())
        ],
        "proximity_edge_colour_pairs": [list(value) for value in edge_colours],
        "dependency_layer_colour_multisets": [
            [depth, list(colours)] for depth, colours in layer_records
        ],
        "wall_colour_records": [
            [
                [list(value[0][0]), value[0][1], value[0][2], value[0][3], value[0][4]],
                [list(pair) for pair in value[1]],
            ]
            for value in sorted(wall_records)
        ],
        "colour_refinement_is_an_invariant_not_a_completeness_claim": True,
    }


def _section_basis(
    kind: str, e: int | None, endpoint_class: Sequence[int]
) -> list[list[int]] | None:
    """Index the standard endpoint section basis without choosing points."""

    if kind == "P2":
        degree = int(endpoint_class[0])
        if degree < 0:
            return None
        return [
            [x_power, y_power, degree - x_power - y_power]
            for x_power in range(degree + 1)
            for y_power in range(degree - x_power + 1)
        ]
    assert e is not None
    section_coefficient, fiber_coefficient = map(int, endpoint_class)
    if section_coefficient < 0:
        return None
    basis: list[list[int]] = []
    for section_power in range(section_coefficient + 1):
        base_degree = fiber_coefficient - section_power * e
        for left_base_power in range(max(-1, base_degree) + 1):
            basis.append(
                [section_power, left_base_power, base_degree - left_base_power]
            )
    return basis


def _rr_record(
    *,
    self_intersection: int,
    canonical_intersection: int,
    endpoint_section_count: int | None,
    imposed_condition_count: int | None,
) -> dict[str, Any]:
    numerator = self_intersection - canonical_intersection
    genus_numerator = self_intersection + canonical_intersection
    if numerator % 2 or genus_numerator % 2:
        raise ArithmeticError("Riemann--Roch parity failed")
    chi = 1 + numerator // 2
    virtual_projective_dimension = numerator // 2
    arithmetic_genus = 1 + genus_numerator // 2
    return {
        "arithmetic_genus": arithmetic_genus,
        "euler_characteristic_chi_O_D": chi,
        "RR_virtual_projective_dimension": virtual_projective_dimension,
        "endpoint_section_count": endpoint_section_count,
        "ordinary_cluster_condition_count": imposed_condition_count,
        "endpoint_sections_minus_conditions": (
            None
            if endpoint_section_count is None or imposed_condition_count is None
            else endpoint_section_count - imposed_condition_count
        ),
        "expected_projective_dimension_if_conditions_are_independent": (
            None
            if endpoint_section_count is None or imposed_condition_count is None
            else endpoint_section_count - imposed_condition_count - 1
        ),
        "no_vanishing_or_independence_claim": True,
    }


def _proximity_inequalities(
    multiplicities: Sequence[int], children: Sequence[Sequence[int]]
) -> list[dict[str, Any]]:
    answer = []
    for parent, child_list in enumerate(children):
        right = sum(int(multiplicities[child]) for child in child_list)
        left = int(multiplicities[parent])
        answer.append(
            {
                "parent_center_position_zero_based": parent,
                "child_center_positions_zero_based": list(child_list),
                "parent_multiplicity": left,
                "sum_proximate_child_multiplicities": right,
                "proximity_excess": left - right,
                "inequality_holds": left >= right,
            }
        )
    return answer


def _linear_system_packet(
    *,
    name: str,
    multiple: int,
    kind: str,
    e: int | None,
    canonical: Sequence[int],
    center_count: int,
    children: Sequence[Sequence[int]],
) -> dict[str, Any]:
    coordinates = [-multiple * int(value) for value in canonical]
    model_rank = 1 if kind == "P2" else 2
    endpoint_class = coordinates[:model_rank]
    multiplicities = [-value for value in coordinates[model_rank:]]
    if multiplicities != [multiple] * center_count:
        raise ArithmeticError(f"{name} does not have the expected uniform multiplicities")
    section_basis = _section_basis(kind, e, endpoint_class)
    row_descriptors = []
    for center in range(center_count):
        for total_order in range(multiple):
            for first_order in range(total_order + 1):
                row_descriptors.append(
                    {
                        "center_position_zero_based": center,
                        "local_jet_total_order": total_order,
                        "local_jet_multiindex": [first_order, total_order - first_order],
                    }
                )
    self_intersection = _pair(coordinates, coordinates, kind, e)
    canonical_intersection = _pair(canonical, coordinates, kind, e)
    rr = _rr_record(
        self_intersection=self_intersection,
        canonical_intersection=canonical_intersection,
        endpoint_section_count=None if section_basis is None else len(section_basis),
        imposed_condition_count=len(row_descriptors),
    )
    column_kind = "P2_homogeneous_exponents" if kind == "P2" else "F_e_pushforward_summand_and_base_exponents"
    inequalities = _proximity_inequalities(multiplicities, children)
    return {
        "system": name,
        "surface_class_coordinates": coordinates,
        "endpoint_class_coordinates": endpoint_class,
        "center_multiplicities_chronological": multiplicities,
        "self_intersection": self_intersection,
        "canonical_intersection": canonical_intersection,
        "riemann_roch": rr,
        "proximity_inequalities": inequalities,
        "all_proximity_inequalities_hold": all(row["inequality_holds"] for row in inequalities),
        "evaluation_matrix_template": {
            "matrix_dimensions_rows_columns": [
                len(row_descriptors),
                0 if section_basis is None else len(section_basis),
            ],
            "section_column_index_kind": column_kind,
            "section_column_indices": section_basis,
            "jet_row_descriptors": row_descriptors,
            "jet_order_count_matrix_by_center": [
                [order + 1 for order in range(multiple)] for _ in range(center_count)
            ],
            "formal_entry": "jet(center,local_multiindex)(endpoint_section_column)",
            "entries_require_stage_B2_universal_center_and_blowup_chart_coordinates": True,
            "matrix_rank_not_asserted": True,
        },
    }


def compile_endpoint_program(program: Mapping[str, Any]) -> dict[str, Any]:
    """Compile one arrays-only Stage-A endpoint program."""

    unknown = set(program) - _PROGRAM_KEYS
    missing = _PROGRAM_KEYS - set(program)
    if unknown:
        raise ValueError(
            "Stage B1 accepts the frozen arrays-only endpoint program; unknown fields: "
            + ", ".join(sorted(unknown))
        )
    if missing:
        raise ValueError("endpoint program is missing fields: " + ", ".join(sorted(missing)))
    if program["stage_B_payload_complete"] is not True:
        raise ValueError("Stage-A endpoint program is not marked complete")

    kind = program["endpoint_kind"]
    if kind not in ("P2", "F_e"):
        raise ValueError("endpoint_kind must be P2 or F_e")
    if kind == "P2":
        if program["hirzebruch_e"] is not None or program["hirzebruch_minimal"] is not None:
            raise ValueError("a P2 endpoint cannot have Hirzebruch data")
        e: int | None = None
        model_rank = 1
    else:
        e = _integer(program["hirzebruch_e"], "hirzebruch_e")
        if e < 0 or not isinstance(program["hirzebruch_minimal"], bool):
            raise ValueError("F_e requires e >= 0 and a minimality boolean")
        if bool(program["hirzebruch_minimal"]) != (e != 1):
            raise ValueError("hirzebruch_minimal is inconsistent with e")
        model_rank = 2

    contracted = _integer_vector(
        program["contracted_wall_indices_zero_based"], "contracted_wall_indices"
    )
    order = _integer_vector(
        program["sequential_contraction_order_zero_based"], "sequential_contraction_order"
    )
    if len(set(contracted)) != len(contracted) or len(set(order)) != len(order):
        raise ValueError("contracted indices and contraction order must be sets")
    if set(contracted) != set(order):
        raise ValueError("contraction order does not enumerate the contracted wall set")
    center_count = len(order)
    if center_count == 0:
        raise ValueError("Stage B1 requires at least one contraction")

    input_picard_rank = model_rank + center_count
    model_basis_input = _integer_matrix(
        program["model_basis_in_input_picard_coordinates"],
        "model_basis_in_input_picard_coordinates",
    )
    total_exceptional_input = _integer_matrix(
        program["total_exceptional_basis_in_input_picard_coordinates"],
        "total_exceptional_basis_in_input_picard_coordinates",
    )
    if (
        len(model_basis_input) != model_rank
        or len(total_exceptional_input) != center_count
        or any(len(row) != input_picard_rank for row in [*model_basis_input, *total_exceptional_input])
    ):
        raise ValueError("the input-Picard model/exceptional basis arrays have inconsistent dimensions")
    nef_rays = _integer_matrix(program["nef_rays_in_model_basis"], "nef_rays_in_model_basis")
    expected_nef_ray_count = 1 if kind == "P2" else 2
    if len(nef_rays) != expected_nef_ray_count or any(len(row) != model_rank for row in nef_rays):
        raise ValueError("the endpoint nef-ray array has inconsistent dimensions")
    profile_hash = program["permutation_invariant_profile_sha256"]
    if (
        not isinstance(profile_hash, str)
        or len(profile_hash) != 64
        or any(character not in "0123456789abcdef" for character in profile_hash)
    ):
        raise ValueError("the Stage-A permutation-invariant profile hash is malformed")
    trace = program["contraction_trace"]
    if not isinstance(trace, list) or len(trace) != center_count:
        raise ValueError("the Stage-A contraction trace has the wrong length")
    trace_keys = {
        "step",
        "source_wall_index_zero_based",
        "current_total_exceptional_in_input_basis",
        "current_square",
        "current_canonical_pairing",
    }
    for position, record in enumerate(trace):
        if not isinstance(record, Mapping) or set(record) != trace_keys:
            raise ValueError("a contraction-trace entry is not the frozen arrays-only record")
        if (
            _integer(record["step"], "contraction trace step") != position + 1
            or _integer(record["source_wall_index_zero_based"], "trace wall index") != order[position]
            or _integer_vector(
                record["current_total_exceptional_in_input_basis"],
                "trace total exceptional",
            )
            != total_exceptional_input[position]
            or _integer(record["current_square"], "trace square") != -1
            or _integer(record["current_canonical_pairing"], "trace canonical pairing") != -1
        ):
            raise ValueError("the contraction trace disagrees with the endpoint program arrays")

    coordinates = _integer_matrix(
        program["all_wall_coordinates_in_endpoint_blowup_basis"],
        "all_wall_coordinates_in_endpoint_blowup_basis",
    )
    if not coordinates:
        raise ValueError("the Stage-A wall coordinate array is empty")
    width = model_rank + center_count
    if any(len(row) != width for row in coordinates):
        raise ValueError("a wall coordinate row has the wrong endpoint/blowup width")
    wall_count = len(coordinates)
    if any(index < 0 or index >= wall_count for index in contracted):
        raise ValueError("a contracted wall index is outside the wall coordinate array")

    endpoint_rows = _integer_matrix(
        program["all_wall_endpoint_image_coordinates"],
        "all_wall_endpoint_image_coordinates",
    )
    if endpoint_rows != [row[:model_rank] for row in coordinates]:
        raise ValueError("endpoint image rows are not the model-coordinate prefixes")

    canonical = _integer_vector(
        program["canonical_coordinates_in_endpoint_blowup_basis"],
        "canonical_coordinates_in_endpoint_blowup_basis",
    )
    expected_canonical = (
        [-3] if kind == "P2" else [-2, -(int(e) + 2)]
    ) + [1] * center_count
    if canonical != expected_canonical:
        raise ValueError("the endpoint/blowup canonical coordinates are inconsistent")
    basis_gram = _integer_matrix(
        program["endpoint_blowup_basis_gram"], "endpoint_blowup_basis_gram"
    )
    if basis_gram != _block_gram(kind, e, center_count):
        raise ValueError("the endpoint/blowup basis Gram is inconsistent")

    strict_gram = _integer_matrix(
        program["strict_contracted_wall_gram"], "strict_contracted_wall_gram"
    )
    sorted_contracted = sorted(contracted)
    reconstructed_strict_gram = [
        [
            _pair(coordinates[left], coordinates[right], kind, e)
            for right in sorted_contracted
        ]
        for left in sorted_contracted
    ]
    if strict_gram != reconstructed_strict_gram:
        raise ValueError("strict contracted-wall Gram does not match the wall coordinates")

    strict_rows_blowdown = [coordinates[index] for index in order]
    strict_exceptional_blowdown = []
    for step, row in enumerate(strict_rows_blowdown):
        if any(row[:model_rank]):
            raise ValueError("a contracted strict exceptional has nonzero endpoint image")
        exceptional = row[model_rank:]
        if exceptional[step] != 1:
            raise ValueError("strict/total exceptional matrix does not have diagonal one")
        if any(exceptional[index] for index in range(step + 1, center_count)):
            raise ValueError("strict/total exceptional matrix is not blowdown-triangular")
        if any(value not in (0, -1) for value in exceptional[:step]):
            raise ValueError("a strict exceptional coefficient is not 0 or -1")
        strict_exceptional_blowdown.append(exceptional)

    # Reverse both rows and columns: rows/columns now run from the oldest
    # blowup centre to the newest.  Thus P[parent,child] = -1 records that
    # the child centre lies on the parent's strict exceptional transform.
    strict_to_total = [
        [strict_exceptional_blowdown[center_count - 1 - row][center_count - 1 - column]
         for column in range(center_count)]
        for row in range(center_count)
    ]
    proximity = [[0] * center_count for _ in range(center_count)]
    parents: list[list[int]] = [[] for _ in range(center_count)]
    children: list[list[int]] = [[] for _ in range(center_count)]
    for parent in range(center_count):
        for child in range(center_count):
            if parent == child:
                if strict_to_total[parent][child] != 1:
                    raise ArithmeticError("chronological proximity matrix lost its diagonal")
                continue
            value = strict_to_total[parent][child]
            if child < parent and value:
                raise ValueError("a centre is proximate to a later blowup centre")
            if value not in (0, -1):
                raise ValueError("chronological strict/total matrix is not a proximity matrix")
            if value == -1:
                proximity[child][parent] = 1
                parents[child].append(parent)
                children[parent].append(child)
    if any(len(row) > 2 for row in parents):
        raise ValueError("a blowup centre is proximate to more than two older centres")

    satellite_parent_intersections: list[dict[str, Any]] = []
    for child, parent_list in enumerate(parents):
        if len(parent_list) != 2:
            continue
        first, second = parent_list
        # Immediately before centre ``child`` is blown up, the strict
        # transforms of its two claimed exceptional parents must intersect
        # once.  Build their prefix strict rows using only earlier centres;
        # later blowups cannot be used to manufacture a satellite relation.
        first_row = [
            int(first == column) - proximity[column][first]
            for column in range(child)
        ]
        second_row = [
            int(second == column) - proximity[column][second]
            for column in range(child)
        ]
        intersection_before_blowup = -sum(
            first_row[index] * second_row[index] for index in range(child)
        )
        satellite_parent_intersections.append(
            {
                "satellite_center_position_zero_based": child,
                "parent_center_positions_zero_based": parent_list,
                "parent_strict_intersection_before_blowup": intersection_before_blowup,
            }
        )
        if intersection_before_blowup != 1:
            raise ValueError(
                "the two claimed satellite parents do not intersect immediately before the blowup"
            )

    depths: list[int] = []
    for center, parent_list in enumerate(parents):
        if any(parent >= center for parent in parent_list):
            raise ValueError("the proximity graph is not chronological")
        depths.append(0 if not parent_list else 1 + max(depths[parent] for parent in parent_list))
    centers = []
    for center in range(center_count):
        blowdown_position = center_count - 1 - center
        parent_count = len(parents[center])
        center_kind = "proper" if parent_count == 0 else "free" if parent_count == 1 else "satellite"
        centers.append(
            {
                "center_position_zero_based": center,
                "center_id": f"p{center + 1}",
                "stage_A_total_exceptional_column_zero_based": blowdown_position,
                "source_contracted_wall_index_zero_based": order[blowdown_position],
                "parent_center_positions_zero_based": parents[center],
                "child_center_positions_zero_based": children[center],
                "center_kind": center_kind,
                "dependency_depth": depths[center],
                "strict_exceptional_coordinates_in_chronological_total_basis": strict_to_total[center],
            }
        )
    layers = [
        [center for center, value in enumerate(depths) if value == depth]
        for depth in sorted(set(depths))
    ]

    boundary_weight = [0] * wall_count
    boundary_records = program["boundary"]
    if not isinstance(boundary_records, list) or not boundary_records:
        raise ValueError("the endpoint program has no marked boundary")
    for record in boundary_records:
        if not isinstance(record, Mapping):
            raise ValueError("a boundary entry is not an object")
        allowed_boundary = {
            "wall_index_zero_based",
            "anti_bicanonical_weight",
            "endpoint_blowup_coordinates",
            "endpoint_image_coordinates",
        }
        if set(record) != allowed_boundary:
            raise ValueError("a boundary entry is not the frozen arrays-only record")
        index = _integer(record["wall_index_zero_based"], "boundary wall index")
        weight = _integer(record["anti_bicanonical_weight"], "boundary weight")
        if index < 0 or index >= wall_count or weight <= 0 or boundary_weight[index]:
            raise ValueError("boundary indices/weights are invalid")
        if _integer_vector(record["endpoint_blowup_coordinates"], "boundary coordinates") != coordinates[index]:
            raise ValueError("boundary coordinates disagree with the wall table")
        if _integer_vector(record["endpoint_image_coordinates"], "boundary endpoint image") != coordinates[index][:model_rank]:
            raise ValueError("boundary endpoint image disagrees with the wall table")
        boundary_weight[index] = weight
    weighted_boundary = [
        sum(boundary_weight[wall] * coordinates[wall][column] for wall in range(wall_count))
        for column in range(width)
    ]
    minus_two_K = [-2 * value for value in canonical]
    if weighted_boundary != minus_two_K:
        raise ValueError("the weighted boundary coordinates do not sum to -2K")

    divisor_rows = []
    for wall_index, raw_row in enumerate(coordinates):
        chronological_exceptional = list(reversed(raw_row[model_rank:]))
        multiplicities = [-value for value in chronological_exceptional]
        endpoint_class = raw_row[:model_rank]
        self_intersection = _pair(raw_row, raw_row, kind, e)
        canonical_intersection = _pair(canonical, raw_row, kind, e)
        section_basis = _section_basis(kind, e, endpoint_class)
        nonnegative_multiplicities = all(value >= 0 for value in multiplicities)
        condition_count = (
            sum(value * (value + 1) // 2 for value in multiplicities)
            if nonnegative_multiplicities
            else None
        )
        inequalities = _proximity_inequalities(multiplicities, children)
        rr = _rr_record(
            self_intersection=self_intersection,
            canonical_intersection=canonical_intersection,
            endpoint_section_count=None if section_basis is None else len(section_basis),
            imposed_condition_count=condition_count,
        )
        row: dict[str, Any] = {
            "wall_index_zero_based": wall_index,
            "is_contracted_wall": wall_index in set(contracted),
            "anti_bicanonical_boundary_weight": boundary_weight[wall_index],
            "class_coordinates_chronological": [*endpoint_class, *chronological_exceptional],
            "endpoint_class_coordinates": endpoint_class,
            "exceptional_coefficients_chronological": chronological_exceptional,
            "multiplicities_chronological": multiplicities,
            "self_intersection": self_intersection,
            "canonical_intersection": canonical_intersection,
            "riemann_roch": rr,
            "all_multiplicities_nonnegative": nonnegative_multiplicities,
            "proximity_inequalities": inequalities,
            "all_proximity_inequalities_hold": all(value["inequality_holds"] for value in inequalities),
        }
        if kind == "P2":
            row["plane_degree"] = endpoint_class[0]
        else:
            row["hirzebruch_section_fiber_bidegree"] = endpoint_class
        divisor_rows.append(row)

    # P is the strict-versus-total exceptional matrix; Q has one in
    # (child,parent) precisely for a proximity edge, so P = I - Q^T.
    identity_minus_q_transpose = [
        [int(row == column) - proximity[column][row] for column in range(center_count)]
        for row in range(center_count)
    ]
    if identity_minus_q_transpose != strict_to_total:
        raise ArithmeticError("P = I - Q^T failed")
    strict_exceptional_gram = _matmul(
        _matmul(strict_to_total, [[-int(i == j) for j in range(center_count)] for i in range(center_count)]),
        _transpose(strict_to_total),
    )

    column_labels = (
        ["H"] if kind == "P2" else ["S_negative", "F_fiber"]
    ) + [f"E_p{index + 1}_total" for index in range(center_count)]
    canonical_chronological = [*canonical[:model_rank], *reversed(canonical[model_rank:])]
    class_rows = {
        "K": canonical_chronological,
        "minus_K": [-value for value in canonical_chronological],
        "minus_two_K": [-2 * value for value in canonical_chronological],
        "weighted_boundary_sum": [*weighted_boundary[:model_rank], *reversed(weighted_boundary[model_rank:])],
    }
    if class_rows["minus_two_K"] != class_rows["weighted_boundary_sum"]:
        raise ArithmeticError("chronological boundary sum changed")

    anti_systems = {
        "minus_K": _linear_system_packet(
            name="-K",
            multiple=1,
            kind=kind,
            e=e,
            canonical=canonical_chronological,
            center_count=center_count,
            children=children,
        ),
        "minus_two_K": _linear_system_packet(
            name="-2K",
            multiple=2,
            kind=kind,
            e=e,
            canonical=canonical_chronological,
            center_count=center_count,
            children=children,
        ),
    }

    center_kind_histogram = Counter(str(row["center_kind"]) for row in centers)
    invariant_candidates = [
        _refined_unlabelled_profile(
            endpoint_kind=kind,
            hirzebruch_e=e,
            centers=centers,
            divisor_rows=divisor_rows,
            swap_f0_rulings=False,
        )
    ]
    if kind == "F_e" and e == 0:
        invariant_candidates.append(
            _refined_unlabelled_profile(
                endpoint_kind=kind,
                hirzebruch_e=e,
                centers=centers,
                divisor_rows=divisor_rows,
                swap_f0_rulings=True,
            )
        )
    invariant_profile = min(
        invariant_candidates,
        key=lambda value: json.dumps(value, sort_keys=True, separators=(",", ":")),
    )

    result = {
        "schema": OUTPUT_SCHEMA,
        "scope": (
            "PB-Atlas numerical Stage B1 only: reverse blowup/proximity, divisor/RR tables, "
            "and indexed -K/-2K evaluation-matrix templates; no coordinate incidence entries, "
            "rank/nonemptiness, bad factors, irreducibility, same-P return, or L4"
        ),
        "source_contract": {
            "schema": INPUT_SCHEMA,
            "arrays_and_integer_indices_only": True,
            "stage_A_program_sha256": _json_hash(dict(program)),
            "no_source_or_surface_labels_entered": True,
        },
        "endpoint": {
            "kind": kind,
            "hirzebruch_e": e,
            "model_rank": model_rank,
            "model_basis_column_labels": column_labels[:model_rank],
            "F0_ruling_swap_is_an_unoriented_symmetry": kind == "F_e" and e == 0,
        },
        "reverse_blowup_program": {
            "chronological_order_is_reverse_of_stage_A_sequential_contraction_order": True,
            "stage_A_contraction_order_zero_based": order,
            "reverse_source_wall_order_zero_based": list(reversed(order)),
            "center_records": centers,
            "dependency_layers_center_positions_zero_based": layers,
            "proximity_edges_child_to_parent_zero_based": [
                [child, parent]
                for child in range(center_count)
                for parent in parents[child]
            ],
            "proximity_incidence_matrix_Q_child_parent": proximity,
            "strict_exceptional_to_total_matrix_P_chronological": strict_to_total,
            "identity_P_equals_I_minus_Q_transpose": True,
            "strict_exceptional_gram_minus_P_P_transpose": strict_exceptional_gram,
            "satellite_parent_intersection_audit": satellite_parent_intersections,
            "center_kind_histogram": dict(sorted(center_kind_histogram.items())),
            "invariant_object_is_DAG_with_commuting_orders_not_the_raw_wall_indices": True,
        },
        "divisor_table": divisor_rows,
        "class_matrices": {
            "basis_column_labels": column_labels,
            "intersection_gram": _block_gram(kind, e, center_count),
            "row_order": ["K", "minus_K", "minus_two_K", "weighted_boundary_sum"],
            "rows": class_rows,
            "minus_two_K_equals_weighted_boundary_sum": True,
        },
        "anti_canonical_evaluation_matrix_templates": anti_systems,
        "invariant_cluster_profile": invariant_profile,
        "invariant_cluster_profile_sha256": _json_hash(invariant_profile),
        "checks": {
            "strict_contracted_gram_reconstructed": True,
            "reverse_blowup_order_is_chronological": True,
            "proximity_DAG_is_acyclic": True,
            "every_center_has_at_most_two_proximity_parents": True,
            "every_satellite_parent_pair_intersects_before_blowup": True,
            "P_equals_I_minus_Q_transpose": True,
            "weighted_boundary_sum_is_minus_two_K": True,
            "all_RR_parities_integral": True,
            "minus_K_and_minus_two_K_matrix_rows_and_columns_indexed": True,
        },
        "stage_B2_remaining": [
            "universal endpoint point coordinates and iterated blowup-chart coordinates",
            "formal evaluation-matrix entries and exact incidence/jet ideals",
            "collision, Jacobian, irreducibility, Hilbert and avoidance bad factors",
            "rank and nonempty regular-open certificates",
            "hidden-wall exclusion and literal same-P return",
        ],
        "all_checks_pass": True,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="one endpointProgram JSON object")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    program = json.loads(arguments.input.read_text(encoding="utf-8"))
    certificate = compile_endpoint_program(program)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(rendered, end="")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
