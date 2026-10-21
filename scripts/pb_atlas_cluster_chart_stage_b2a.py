#!/usr/bin/env python3
"""PB-Atlas Stage B2a: exact P2 cluster charts and jet matrices.

The compiler consumes one frozen numerical Stage-B1 packet.  It uses only
the endpoint kind, chronological centre records, proximity parents, and the
already indexed P2 section/jet templates.  It constructs one explicit chart
of the ordered cluster space and replaces every formal ``-K``/``-2K`` jet
token by an exact rational-polynomial entry.

The chart is row agnostic.  Proper centres use affine P2 patches; free
centres use one of the two standard blowup charts; satellites use the chart
forced by their two exceptional axes.  A deterministic four-point PGL3
gauge is used on its declared open set when at least four proper centres are
available.  Alternate gauge, affine, and free-direction charts are recorded
as a finite-cover boundary rather than silently asserted unnecessary.

This layer makes no rank, nonemptiness, irreducibility, hidden-wall,
full-nef, or same-P assertion.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib
import itertools
import json
import math
from typing import Any, Mapping, Sequence

import sympy as sp


INPUT_SCHEMA = "coxeter.pb_atlas.numerical_cluster_stage_b1.v1"
OUTPUT_SCHEMA = "coxeter.pb_atlas.cluster_chart_stage_b2a.v1"


def _json_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _rational(value: Any) -> sp.Rational:
    if isinstance(value, bool):
        raise ValueError("a chart-frame entry cannot be boolean")
    if isinstance(value, (float, sp.Float)):
        raise ValueError("a chart-frame entry must be an exact rational, not an approximate number")
    try:
        result = sp.Rational(value)
    except (TypeError, ValueError) as exception:
        raise ValueError("a chart-frame entry must be an exact rational") from exception
    return result


def _strict_integer(value: Any, name: str, *, nonnegative: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if nonnegative and value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _option_index(value: Any, name: str) -> int:
    """Read an integer option index, allowing JSON object keys like ``"3"``."""

    if isinstance(value, str):
        if (
            not value.isascii()
            or not value.isdecimal()
            or (value != "0" and value.startswith("0"))
        ):
            raise ValueError(f"{name} must be a nonnegative integer index")
        return int(value)
    return _strict_integer(value, name, nonnegative=True)


def _indexed_option_items(
    mapping: Mapping[Any, Any], name: str
) -> list[tuple[int, Any]]:
    result: list[tuple[int, Any]] = []
    seen: set[int] = set()
    for raw_key, value in mapping.items():
        key = _option_index(raw_key, name)
        if key in seen:
            raise ValueError(f"{name} occurs more than once after index normalization")
        seen.add(key)
        result.append((key, value))
    return result


def _expression_record(expression: sp.Expr) -> dict[str, str]:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    return {
        "numerator": sp.sstr(sp.expand(numerator)),
        "denominator": sp.sstr(sp.expand(denominator)),
    }


def _expression_string(expression: sp.Expr) -> str:
    record = _expression_record(expression)
    if record["denominator"] == "1":
        return record["numerator"]
    return f"({record['numerator']})/({record['denominator']})"


def _validate_packet_shape(packet: Mapping[str, Any]) -> None:
    if not isinstance(packet, Mapping):
        raise ValueError("the Stage-B1 packet must be a mapping")
    if packet.get("schema") != INPUT_SCHEMA:
        raise ValueError("Stage B2a requires the frozen Stage-B1 packet schema")
    if packet.get("all_checks_pass") is not True:
        raise ValueError("the Stage-B1 packet is not certified")
    endpoint = packet.get("endpoint")
    if not isinstance(endpoint, Mapping) or endpoint.get("kind") != "P2":
        raise ValueError("Stage B2a currently accepts P2 endpoint packets only")
    centers = packet.get("reverse_blowup_program", {}).get("center_records")
    if not isinstance(centers, list) or not centers:
        raise ValueError("the Stage-B1 packet has no chronological centre records")
    for position, center in enumerate(centers):
        if not isinstance(center, Mapping):
            raise ValueError("a Stage-B1 centre record is malformed")
        center_position = center.get("center_position_zero_based")
        if (
            isinstance(center_position, bool)
            or not isinstance(center_position, int)
            or center_position != position
        ):
            raise ValueError("Stage-B1 centre positions are not chronological")
        parents = center.get("parent_center_positions_zero_based")
        kind = center.get("center_kind")
        if not isinstance(parents, list) or any(
            isinstance(parent, bool)
            or not isinstance(parent, int)
            or parent < 0
            or parent >= position
            for parent in parents
        ):
            raise ValueError("a Stage-B1 proximity parent is malformed")
        if len(set(parents)) != len(parents):
            raise ValueError("a Stage-B1 centre has duplicate proximity parents")
        expected = {0: "proper", 1: "free", 2: "satellite"}.get(len(parents))
        if expected is None or kind != expected:
            raise ValueError("centre kind and proximity-parent count disagree")


def _normalized_options(
    centers: Sequence[Mapping[str, Any]], options: Mapping[str, Any] | None
) -> dict[str, Any]:
    if options is None:
        options = {}
    elif not isinstance(options, Mapping):
        raise ValueError("Stage-B2a chart options must be a mapping")
    else:
        options = dict(options)
    allowed = {
        "pgl_gauge_center_positions_zero_based",
        "proper_affine_patch_by_center",
        "proper_frame_by_center",
        "free_chart_mode_by_center",
    }
    unknown = set(options) - allowed
    if unknown:
        raise ValueError(
            "unknown Stage-B2a chart options: "
            + ", ".join(sorted(repr(value) for value in unknown))
        )

    proper = [
        index for index, center in enumerate(centers) if center["center_kind"] == "proper"
    ]
    raw_gauge = options.get("pgl_gauge_center_positions_zero_based", "automatic")
    if raw_gauge == "automatic":
        gauge = proper[:4] if len(proper) >= 4 else []
        gauge_rule = (
            "first four chronological proper centres"
            if gauge
            else "no gauge because fewer than four proper centres exist"
        )
    else:
        if not isinstance(raw_gauge, (list, tuple)):
            raise ValueError("the PGL3 gauge must be automatic or a centre array")
        gauge = [_option_index(value, "a PGL3 gauge centre") for value in raw_gauge]
        gauge_rule = (
            "caller-supplied ordered proper-centre quadruple"
            if gauge
            else "gauge explicitly disabled by caller"
        )
    if gauge and (len(gauge) != 4 or len(set(gauge)) != 4):
        raise ValueError("a full PGL3 gauge requires four distinct proper centres")
    if any(position not in proper for position in gauge):
        raise ValueError("the PGL3 gauge contains a nonproper centre")

    patch_input = options.get("proper_affine_patch_by_center", {})
    frame_input = options.get("proper_frame_by_center", {})
    mode_input = options.get("free_chart_mode_by_center", {})
    if not all(isinstance(value, Mapping) for value in (patch_input, frame_input, mode_input)):
        raise ValueError("chart-option maps must be mappings")

    patches = {
        key: _strict_integer(
            value, "a proper affine patch", nonnegative=True
        )
        for key, value in _indexed_option_items(
            patch_input, "a proper-patch centre"
        )
    }
    if any(position not in proper or patch not in (0, 1, 2) for position, patch in patches.items()):
        raise ValueError("a proper affine patch choice is invalid")
    if set(patches).intersection(gauge):
        raise ValueError("a proper affine patch choice cannot target a PGL3-gauged centre")
    frames: dict[int, list[list[sp.Rational]]] = {}
    for position, value in _indexed_option_items(
        frame_input, "a proper-frame centre"
    ):
        if position not in proper or not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError("a proper local frame is malformed")
        if any(not isinstance(row, (list, tuple)) for row in value):
            raise ValueError("a proper local frame is malformed")
        rows = [[_rational(entry) for entry in row] for row in value]
        if any(len(row) != 2 for row in rows):
            raise ValueError("a proper local frame must be 2 by 2")
        if sp.det(sp.Matrix(rows)) == 0:
            raise ValueError("a proper local frame must be invertible")
        frames[position] = rows

    modes = {
        key: value
        for key, value in _indexed_option_items(mode_input, "a free-chart centre")
    }
    free = [index for index, center in enumerate(centers) if center["center_kind"] == "free"]
    if any(position not in free or mode not in ("A", "B") for position, mode in modes.items()):
        raise ValueError("a free blowup chart choice is invalid")
    return {
        "gauge": gauge,
        "gauge_rule": gauge_rule,
        "proper": proper,
        "patches": patches,
        "frames": frames,
        "modes": modes,
    }


def _proper_point_and_local_map(
    *,
    position: int,
    gauge: Sequence[int],
    patch_choice: int | None,
    frame: Sequence[Sequence[sp.Rational]] | None,
    u: sp.Symbol,
    v: sp.Symbol,
) -> tuple[list[sp.Expr], int, list[sp.Symbol], dict[sp.Symbol, sp.Expr], list[str]]:
    fixed_targets = (
        (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(1), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(0), sp.Integer(1)),
        (sp.Integer(1), sp.Integer(1), sp.Integer(1)),
    )
    variables: list[sp.Symbol] = []
    open_conditions: list[str] = []
    if position in gauge:
        target = list(fixed_targets[list(gauge).index(position)])
        patch = next(index for index, value in enumerate(target) if value)
        target = [sp.cancel(value / target[patch]) for value in target]
    else:
        patch = 2 if patch_choice is None else patch_choice
        q0, q1 = sp.symbols(f"q{position}_0 q{position}_1")
        variables.extend((q0, q1))
        target = [sp.Integer(0)] * 3
        target[patch] = sp.Integer(1)
        free = [index for index in range(3) if index != patch]
        target[free[0]], target[free[1]] = q0, q1
        open_conditions.append(f"original homogeneous coordinate {patch} of p{position + 1} != 0")

    free = [index for index in range(3) if index != patch]
    frame_matrix = (
        [[sp.Integer(1), sp.Integer(0)], [sp.Integer(0), sp.Integer(1)]]
        if frame is None
        else [[sp.Rational(value) for value in row] for row in frame]
    )
    local = list(target)
    local[free[0]] = target[free[0]] + frame_matrix[0][0] * u + frame_matrix[0][1] * v
    local[free[1]] = target[free[1]] + frame_matrix[1][0] * u + frame_matrix[1][1] * v
    substitution = {sp.Symbol("X"): local[0], sp.Symbol("Y"): local[1], sp.Symbol("Z"): local[2]}
    return target, patch, variables, substitution, open_conditions


def _build_chart(
    packet: Mapping[str, Any], options: Mapping[str, Any] | None
) -> tuple[dict[str, Any], dict[str, Any]]:
    _validate_packet_shape(packet)
    centers = packet["reverse_blowup_program"]["center_records"]
    normalized = _normalized_options(centers, options)
    gauge = normalized["gauge"]
    proper = normalized["proper"]

    local_variables = [sp.symbols(f"u{index} v{index}") for index in range(len(centers))]
    global_variables: list[sp.Symbol] = []
    proper_substitutions: dict[int, dict[sp.Symbol, sp.Expr]] = {}
    chart_records: list[dict[str, Any]] = []
    axes: list[dict[int, str]] = [dict() for _ in centers]
    primary_parent: list[int | None] = [None for _ in centers]
    parent_to_child_substitution: dict[int, dict[sp.Symbol, sp.Expr]] = {}
    chart_mode: list[str] = ["proper" for _ in centers]
    direction_parameter: dict[int, sp.Expr] = {}
    open_conditions: list[str] = []

    for position, center in enumerate(centers):
        u, v = local_variables[position]
        kind = center["center_kind"]
        parents = list(center["parent_center_positions_zero_based"])
        if kind == "proper":
            target, patch, variables, substitution, local_open = _proper_point_and_local_map(
                position=position,
                gauge=gauge,
                patch_choice=normalized["patches"].get(position),
                frame=normalized["frames"].get(position),
                u=u,
                v=v,
            )
            global_variables.extend(variables)
            open_conditions.extend(local_open)
            proper_substitutions[position] = substitution
            chart_records.append(
                {
                    "center_position_zero_based": position,
                    "center_kind": kind,
                    "parent_center_positions_zero_based": [],
                    "primary_parent_center_position_zero_based": None,
                    "local_variables": [str(u), str(v)],
                    "homogeneous_center_coordinates": [_expression_string(value) for value in target],
                    "proper_affine_patch_coordinate_zero_based": patch,
                    "proper_local_frame": [
                        [str(value) for value in row]
                        for row in normalized["frames"].get(
                            position,
                            [[sp.Integer(1), sp.Integer(0)], [sp.Integer(0), sp.Integer(1)]],
                        )
                    ],
                    "exceptional_axis_by_parent_center": {},
                    "chart_mode": "proper_affine",
                    "direction_parameter": None,
                    "parent_local_substitution": None,
                }
            )
            continue

        if kind == "free":
            parent = parents[0]
            mode = normalized["modes"].get(position, "A")
            parameter = sp.Symbol(f"t{position}")
            global_variables.append(parameter)
            primary_parent[position] = parent
            chart_mode[position] = mode
            direction_parameter[position] = parameter
            up, vp = local_variables[parent]
            if mode == "A":
                substitution = {up: u, vp: u * (parameter + v)}
                axes[position] = {parent: "u"}
                if "v" in axes[parent].values():
                    open_conditions.append(f"t{position} != 0")
            else:
                substitution = {up: v * (parameter + u), vp: v}
                axes[position] = {parent: "v"}
                if "u" in axes[parent].values():
                    open_conditions.append(f"t{position} != 0")
            parent_to_child_substitution[position] = substitution
        else:
            chosen: tuple[int, int, str] | None = None
            for candidate in sorted(parents, reverse=True):
                other = next(parent for parent in parents if parent != candidate)
                if other in axes[candidate]:
                    chosen = (candidate, other, axes[candidate][other])
                    break
            if chosen is None:
                raise ValueError(
                    "satellite parents are not visible as exceptional axes in the chosen chart"
                )
            parent, other, other_axis = chosen
            primary_parent[position] = parent
            up, vp = local_variables[parent]
            if other_axis == "v":
                mode = "A"
                substitution = {up: u, vp: u * v}
                axes[position] = {parent: "u", other: "v"}
            else:
                mode = "B"
                substitution = {up: v * u, vp: v}
                axes[position] = {other: "u", parent: "v"}
            chart_mode[position] = mode
            direction_parameter[position] = sp.Integer(0)
            parent_to_child_substitution[position] = substitution

        chart_records.append(
            {
                "center_position_zero_based": position,
                "center_kind": kind,
                "parent_center_positions_zero_based": parents,
                "primary_parent_center_position_zero_based": primary_parent[position],
                "local_variables": [str(u), str(v)],
                "homogeneous_center_coordinates": None,
                "proper_affine_patch_coordinate_zero_based": None,
                "proper_local_frame": None,
                "exceptional_axis_by_parent_center": {
                    str(parent): axis for parent, axis in sorted(axes[position].items())
                },
                "chart_mode": f"blowup_{chart_mode[position]}",
                "direction_parameter": _expression_string(direction_parameter[position]),
                "parent_local_substitution": {
                    str(variable): _expression_string(value)
                    for variable, value in parent_to_child_substitution[position].items()
                },
            }
        )

    # In the primary chart, finite-slope siblings on the same exceptional
    # must have distinct parameters.  Other collision/avoidance factors are
    # explicitly deferred.
    free_by_parent: dict[int, list[int]] = defaultdict(list)
    for position, center in enumerate(centers):
        if center["center_kind"] == "free" and chart_mode[position] == "A":
            free_by_parent[int(primary_parent[position])].append(position)
    for siblings in free_by_parent.values():
        for left, right in itertools.combinations(siblings, 2):
            open_conditions.append(f"t{left} - t{right} != 0")

    gauge_candidates = [list(value) for value in itertools.combinations(proper, 4)]
    gauge_open = []
    if gauge:
        for omitted in range(4):
            triple = [position for index, position in enumerate(gauge) if index != omitted]
            gauge_open.append("Delta_" + "_".join(f"p{position + 1}" for position in triple))
        open_conditions.extend(f"{factor} != 0" for factor in gauge_open)

    free_count = sum(center["center_kind"] == "free" for center in centers)
    ungauged_proper_count = len(proper) - len(gauge)
    public = {
        "pgl3_gauge": {
            "full_gauge_used": bool(gauge),
            "selected_proper_center_positions_zero_based": gauge,
            "fixed_target_points_in_selected_order": (
                [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]] if gauge else []
            ),
            "required_nonzero_triple_determinants": gauge_open,
            "deterministic_rule": normalized["gauge_rule"],
            "gauge_is_only_valid_on_the_declared_open_set": True,
        },
        "coordinate_variables": sorted({str(variable) for variable in global_variables}),
        "center_charts": chart_records,
        "chart_open_conditions_not_asserted_nonempty": sorted(set(open_conditions)),
        "finite_chart_cover_boundary": {
            "proper_center_count": len(proper),
            "possible_four_point_gauge_subsets": len(gauge_candidates),
            "four_point_gauge_subsets_zero_based": gauge_candidates,
            "ungauged_proper_centers_need_three_affine_patches_each": ungauged_proper_count,
            "free_centers_need_two_blowup_charts_each": free_count,
            "coarse_cartesian_chart_count_for_this_gauge": 3**ungauged_proper_count
            * 2**free_count,
            "configurations_with_no_general_position_quadruple_require_the_ungauged_mode": True,
            "this_packet_contains_one_selected_chart_not_a_nonemptiness_proof": True,
        },
    }
    internal = {
        "centers": centers,
        "local_variables": local_variables,
        "proper_substitutions": proper_substitutions,
        "primary_parent": primary_parent,
        "parent_to_child_substitution": parent_to_child_substitution,
        "chart_mode": chart_mode,
        "direction_parameter": direction_parameter,
        "global_variables": global_variables,
    }
    return public, internal


def _split_low_high(
    expression: sp.Expr, u: sp.Symbol, v: sp.Symbol, multiplicity: int
) -> tuple[dict[tuple[int, int], sp.Expr], sp.Expr]:
    parsed = sp.Poly(sp.expand(expression), u, v)
    low: dict[tuple[int, int], sp.Expr] = {}
    high = sp.Integer(0)
    for (left, right), coefficient in parsed.terms():
        if left + right < multiplicity:
            low[(left, right)] = sp.expand(coefficient)
        else:
            high += coefficient * u**left * v**right
    return low, sp.expand(high)


def _condition_matrix(
    *,
    internal: Mapping[str, Any],
    degree: int,
    multiplicities: Sequence[int],
    section_columns: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    centers = internal["centers"]
    degree = _strict_integer(degree, "the plane-section degree", nonnegative=True)
    if not isinstance(multiplicities, (list, tuple)):
        raise ValueError("cluster multiplicities must be an integer array")
    if len(multiplicities) != len(centers) or any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in multiplicities
    ):
        raise ValueError("cluster multiplicities must be nonnegative integers, one per centre")
    if section_columns is None:
        section_columns = [
            [left, middle, degree - left - middle]
            for left in range(degree + 1)
            for middle in range(degree - left + 1)
        ]
    elif not isinstance(section_columns, (list, tuple)):
        raise ValueError("P2 section columns must be an array of exponent rows")
    columns: list[list[int]] = []
    for row in section_columns:
        if not isinstance(row, (list, tuple)) or len(row) != 3:
            raise ValueError("a P2 section column is not a homogeneous degree row")
        parsed = [
            _strict_integer(value, "a P2 section exponent", nonnegative=True)
            for value in row
        ]
        if sum(parsed) != degree:
            raise ValueError("a P2 section column is not a homogeneous degree row")
        columns.append(parsed)
    if not columns:
        raise ValueError("P2 section columns must contain at least one monomial index")
    if len({tuple(row) for row in columns}) != len(columns):
        raise ValueError("P2 section columns must be distinct monomial indices")

    X, Y, Z = sp.symbols("X Y Z")
    monomials = [X**row[0] * Y**row[1] * Z**row[2] for row in columns]
    current: dict[int, list[sp.Expr]] = {}
    for position, center in enumerate(centers):
        if center["center_kind"] == "proper":
            substitution = internal["proper_substitutions"][position]
            current[position] = [sp.expand(monomial.subs(substitution)) for monomial in monomials]

    rows: list[list[sp.Expr]] = []
    descriptors: list[dict[str, Any]] = []
    children_by_primary: dict[int, list[int]] = defaultdict(list)
    for child, parent in enumerate(internal["primary_parent"]):
        if parent is not None:
            children_by_primary[int(parent)].append(child)

    for position, center in enumerate(centers):
        if position not in current:
            raise ValueError("the chosen chart did not construct a local polynomial at a centre")
        u, v = internal["local_variables"][position]
        multiplicity = multiplicities[position]
        splits = [_split_low_high(expression, u, v, multiplicity) for expression in current[position]]
        for total in range(multiplicity):
            for first in range(total + 1):
                multiindex = (first, total - first)
                rows.append([low.get(multiindex, sp.Integer(0)) for low, _high in splits])
                descriptors.append(
                    {
                        "center_position_zero_based": position,
                        "local_jet_total_order": total,
                        "local_jet_multiindex": list(multiindex),
                    }
                )

        for child in children_by_primary.get(position, []):
            uc, vc = internal["local_variables"][child]
            substitution = internal["parent_to_child_substitution"][child]
            divisor = uc if internal["chart_mode"][child] == "A" else vc
            transformed = []
            for _low, high in splits:
                value = sp.cancel(high.subs(substitution, simultaneous=True) / divisor**multiplicity)
                numerator, denominator = sp.fraction(value)
                if denominator.has(uc) or denominator.has(vc):
                    raise ArithmeticError("a blowup transform retained a local Laurent denominator")
                transformed.append(sp.expand(value))
            current[child] = transformed

    serialized = [[_expression_record(entry) for entry in row] for row in rows]
    result = {
        "degree": degree,
        "center_multiplicities_chronological": list(multiplicities),
        "matrix_dimensions_rows_columns": [len(rows), len(columns)],
        "section_column_indices": columns,
        "jet_row_descriptors": descriptors,
        "entries": serialized,
        "entry_domain": "exact rational functions over Q in the selected chart variables",
        "formal_jet_tokens_remaining": False,
        "matrix_rank_not_asserted": True,
    }
    result["entries_sha256"] = _json_hash(serialized)
    return result


def compile_p2_divisor_evaluation_matrix(
    packet: Mapping[str, Any],
    *,
    degree: int,
    multiplicities: Sequence[int],
    chart_options: Mapping[str, Any] | None = None,
    section_columns: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    """Compile any nonnegative P2 cluster multiplicity row on one chart."""

    _public, internal = _build_chart(packet, chart_options)
    return _condition_matrix(
        internal=internal,
        degree=degree,
        multiplicities=multiplicities,
        section_columns=section_columns,
    )


def compile_p2_cluster_chart(
    packet: Mapping[str, Any], chart_options: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Compile the selected universal chart and exact B1 ``-K``/``-2K`` entries."""

    public_chart, internal = _build_chart(packet, chart_options)
    systems = {}
    template_parent = packet["anti_canonical_evaluation_matrix_templates"]
    for key in ("minus_K", "minus_two_K"):
        source = template_parent[key]
        template = source["evaluation_matrix_template"]
        degree = _strict_integer(
            source["endpoint_class_coordinates"][0],
            f"the frozen B1 {key} plane-section degree",
            nonnegative=True,
        )
        multiplicities = [
            _strict_integer(
                value,
                f"a frozen B1 {key} centre multiplicity",
                nonnegative=True,
            )
            for value in source["center_multiplicities_chronological"]
        ]
        matrix = _condition_matrix(
            internal=internal,
            degree=degree,
            multiplicities=multiplicities,
            section_columns=template["section_column_indices"],
        )
        if matrix["matrix_dimensions_rows_columns"] != template["matrix_dimensions_rows_columns"]:
            raise ArithmeticError("the Stage-B2a matrix dimensions changed from Stage B1")
        if matrix["jet_row_descriptors"] != template["jet_row_descriptors"]:
            raise ArithmeticError("the Stage-B2a row descriptors changed from Stage B1")
        if matrix["section_column_indices"] != template["section_column_indices"]:
            raise ArithmeticError("the Stage-B2a section columns changed from Stage B1")
        matrix.update(
            {
                "system": source["system"],
                "source_B1_formal_entry": template["formal_entry"],
                "source_B1_template_sha256": _json_hash(template),
                "coefficient_normalization": "Taylor coefficient, not factorial-scaled derivative",
            }
        )
        systems[key] = matrix

    checks = {
        "input_is_frozen_B1_P2_packet": True,
        "every_center_has_one_exact_local_chart": len(public_chart["center_charts"])
        == len(packet["reverse_blowup_program"]["center_records"]),
        "all_B1_formal_minus_K_tokens_replaced": not systems["minus_K"][
            "formal_jet_tokens_remaining"
        ],
        "all_B1_formal_minus_two_K_tokens_replaced": not systems["minus_two_K"][
            "formal_jet_tokens_remaining"
        ],
        "B1_row_and_column_indices_preserved": True,
        "no_local_Laurent_denominator_survived": True,
    }
    result = {
        "schema": OUTPUT_SCHEMA,
        "scope": "one exact universal-coordinate chart for one frozen P2 Stage-B1 packet",
        "source_contract": {
            "schema": INPUT_SCHEMA,
            "stage_B1_packet_sha256": _json_hash(packet),
            "reads_only_B1_cluster_and_evaluation_arrays": True,
            "no_surface_type_or_calibration_name_entered": True,
        },
        "chart": public_chart,
        "evaluation_matrices": systems,
        "checks": checks,
        "claims_not_made": [
            "chart-open nonemptiness",
            "matrix rank or condition independence",
            "collision/Jacobian/irreducibility/Hilbert bad-factor completeness",
            "strong-Coble boundary geometry",
            "hidden-wall or full-Mori completeness",
            "same-P return",
            "L4",
        ],
        "stage_B2b_remaining": [
            "cover every required alternate chart or prove the selected chart meets the incidence locus",
            "derive and remove collision, Jacobian, reducibility, and Hilbert bad factors",
            "prove a nonempty strong-Coble branch and literal same-P return",
        ],
        "all_checks_pass": all(checks.values()),
    }
    return result


def main() -> None:
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="JSON file containing one Stage-B1 packet")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    packet = json.loads(args.input.read_text(encoding="utf-8"))
    result = compile_p2_cluster_chart(packet)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
