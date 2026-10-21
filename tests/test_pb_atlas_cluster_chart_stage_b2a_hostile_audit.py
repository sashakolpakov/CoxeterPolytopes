"""Independent hostile checks for the PB-Atlas P2 Stage-B2a chart compiler.

The coefficient oracles below do not import the production splitter,
strict-transform routine, or either legacy plane-matrix helper.  They pull
back monomials through the displayed blowup maps and read the shifted Taylor
coefficients directly.
"""

from __future__ import annotations

import copy
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import pytest
import sympy as sp

from scripts.pb_atlas_cluster_chart_stage_b2a import (
    compile_p2_cluster_chart,
    compile_p2_divisor_evaluation_matrix,
)
from scripts.pb_atlas_cluster_stage_b1 import compile_endpoint_program


ROOT = Path(__file__).resolve().parents[1]
STAGE_A = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
STAGE_B1 = ROOT / "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json"
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
B2A_ARCHIVE = ROOT / "data/enumerations/pb_atlas_cluster_chart_stage_b2a.json"
OLD_H10 = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json"
OLD_H11 = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse.json"
CORE = ROOT / "scripts/pb_atlas_cluster_chart_stage_b2a.py"
CERTIFIER = ROOT / "scripts/certify_pb_atlas_cluster_chart_stage_b2a.py"
OWNER_TEST = ROOT / "tests/test_pb_atlas_cluster_chart_stage_b2a.py"
REPORT = ROOT / "results/pb_atlas_cluster_chart_stage_b2a.md"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture_packet(fixture_id: str) -> dict[str, Any]:
    source = next(
        row
        for row in _load(FIXTURES)["valid_cases"]
        if row["fixture_id"] == fixture_id
    )
    return compile_endpoint_program(source["program"])


def _forest_sizes(packet: Mapping[str, Any]) -> list[int]:
    centers = packet["reverse_blowup_program"]["center_records"]
    graph = [set() for _ in centers]
    for child, center in enumerate(centers):
        for parent in center["parent_center_positions_zero_based"]:
            graph[child].add(parent)
            graph[parent].add(child)
    unseen = set(range(len(centers)))
    sizes: list[int] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        size = 0
        while stack:
            vertex = stack.pop()
            size += 1
            for neighbor in graph[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        sizes.append(size)
    return sorted(sizes, reverse=True)


@lru_cache(maxsize=1)
def _p2_sources() -> tuple[dict[str, Any], ...]:
    records: list[dict[str, Any]] = []
    for case in _load(STAGE_A)["cases"]:
        for orbit in case["certificate"]["P2_endpoint_orbits"]:
            records.append(
                {
                    "representative_mask": orbit["representative_mask"],
                    "orbit_size": orbit["orbit_size"],
                    "packet": compile_endpoint_program(orbit["program"]),
                }
            )
    return tuple(records)


def _packet_by_forest(sizes: Sequence[int]) -> dict[str, Any]:
    matches = [row["packet"] for row in _p2_sources() if _forest_sizes(row["packet"]) == list(sizes)]
    assert len(matches) == 1
    return matches[0]


def _expression(entry: Mapping[str, str]) -> sp.Expr:
    return sp.cancel(
        sp.sympify(entry["numerator"]) / sp.sympify(entry["denominator"])
    )


def _matrix(record: Mapping[str, Any]) -> list[list[sp.Expr]]:
    return [[_expression(entry) for entry in row] for row in record["entries"]]


def _coefficient(
    expression: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
    first_degree: int,
    second_degree: int,
) -> sp.Expr:
    return sp.expand(
        sp.Poly(sp.expand(expression), first, second).coeff_monomial(
            first**first_degree * second**second_degree
        )
    )


def _monomials(columns: Sequence[Sequence[int]]) -> list[sp.Expr]:
    X, Y, Z = sp.symbols("X Y Z")
    return [X**left * Y**middle * Z**right for left, middle, right in columns]


def _specialize_matrix(
    matrix: Sequence[Sequence[sp.Expr]], substitutions: Mapping[sp.Symbol, sp.Expr]
) -> list[list[sp.Expr]]:
    return [
        [sp.cancel(value.subs(substitutions, simultaneous=True)) for value in row]
        for row in matrix
    ]


def _assert_exact_rows(
    record: Mapping[str, Any],
    expected_by_descriptor: Mapping[tuple[int, int, int], Sequence[sp.Expr]],
    substitutions: Mapping[sp.Symbol, sp.Expr] | None = None,
) -> None:
    matrix = _matrix(record)
    if substitutions:
        matrix = _specialize_matrix(matrix, substitutions)
    assert len(matrix) == len(expected_by_descriptor)
    for descriptor, row in zip(record["jet_row_descriptors"], matrix):
        key = (
            descriptor["center_position_zero_based"],
            descriptor["local_jet_multiindex"][0],
            descriptor["local_jet_multiindex"][1],
        )
        expected = expected_by_descriptor[key]
        assert len(row) == len(expected)
        assert all(sp.cancel(actual - wanted) == 0 for actual, wanted in zip(row, expected))


def test_direct_proper_free_satellite_oracle_in_both_blowup_charts() -> None:
    """Read shifted coefficients from the composite maps, not the compiler."""

    packet = _fixture_packet("satellite_center_two_parents")
    frame = [[2, 1], [1, 1]]
    multiplicities = [2, 2, 2]
    q0, q1, t = sp.symbols("q0_0 q0_1 t1")
    u0, v0, u1, v1, u2, v2 = sp.symbols("u0 v0 u1 v1 u2 v2")

    for free_mode in ("A", "B"):
        options = {
            "pgl_gauge_center_positions_zero_based": [],
            "proper_affine_patch_by_center": {0: 2},
            "proper_frame_by_center": {0: frame},
            "free_chart_mode_by_center": {1: free_mode},
        }
        result = compile_p2_divisor_evaluation_matrix(
            packet,
            degree=4,
            multiplicities=multiplicities,
            chart_options=options,
        )
        chart = compile_p2_cluster_chart(packet, options)["chart"]["center_charts"]
        assert chart[0]["proper_local_frame"] == [["2", "1"], ["1", "1"]]

        X0 = q0 + 2 * u0 + v0
        Y0 = q1 + u0 + v0
        Z0 = sp.Integer(1)
        local_monomials = [
            sp.expand(monomial.subs({sp.Symbol("X"): X0, sp.Symbol("Y"): Y0, sp.Symbol("Z"): Z0}))
            for monomial in _monomials(result["section_column_indices"])
        ]
        expected: dict[tuple[int, int, int], list[sp.Expr]] = {}

        for total in range(multiplicities[0]):
            for left in range(total + 1):
                right = total - left
                expected[(0, left, right)] = [
                    _coefficient(value, u0, v0, left, right)
                    for value in local_monomials
                ]

        if free_mode == "A":
            assert chart[1]["parent_local_substitution"] == {
                "u0": "u1",
                "v0": "t1*u1 + u1*v1",
            }
            assert chart[2]["chart_mode"] == "blowup_B"
            at_free = [
                sp.expand(value.subs({u0: u1, v0: u1 * (t + v1)}, simultaneous=True))
                for value in local_monomials
            ]
            at_satellite = [
                sp.expand(
                    value.subs(
                        {u0: u2 * v2, v0: u2 * v2 * (t + v2)},
                        simultaneous=True,
                    )
                )
                for value in local_monomials
            ]
            free_shift = lambda left, right: (multiplicities[0] + left, right)
            satellite_shift = lambda left, right: (
                multiplicities[0] + left,
                multiplicities[0] + multiplicities[1] + right,
            )
        else:
            assert chart[1]["parent_local_substitution"] == {
                "u0": "t1*v1 + u1*v1",
                "v0": "v1",
            }
            assert chart[2]["chart_mode"] == "blowup_A"
            at_free = [
                sp.expand(value.subs({u0: v1 * (t + u1), v0: v1}, simultaneous=True))
                for value in local_monomials
            ]
            at_satellite = [
                sp.expand(
                    value.subs(
                        {u0: u2 * v2 * (t + u2), v0: u2 * v2},
                        simultaneous=True,
                    )
                )
                for value in local_monomials
            ]
            free_shift = lambda left, right: (left, multiplicities[0] + right)
            satellite_shift = lambda left, right: (
                multiplicities[0] + multiplicities[1] + left,
                multiplicities[0] + right,
            )

        for total in range(multiplicities[1]):
            for left in range(total + 1):
                right = total - left
                shifted = free_shift(left, right)
                expected[(1, left, right)] = [
                    _coefficient(value, u1, v1, *shifted) for value in at_free
                ]
        for total in range(multiplicities[2]):
            for left in range(total + 1):
                right = total - left
                shifted = satellite_shift(left, right)
                expected[(2, left, right)] = [
                    _coefficient(value, u2, v2, *shifted) for value in at_satellite
                ]
        _assert_exact_rows(result, expected)

        local_symbols = {u0, v0, u1, v1, u2, v2}
        for row in result["entries"]:
            for entry in row:
                assert _expression(entry).free_symbols.isdisjoint(local_symbols)


def test_commuting_centers_are_equivariant_after_transporting_all_indices() -> None:
    first = _fixture_packet("free_forest_two_commuting_roots")
    second = _fixture_packet("same_free_forest_commuted_order")
    center_bijection = {0: 1, 1: 0, 2: 2, 3: 3}
    first_multiplicities = [2, 1, 2, 1]
    second_multiplicities = [0] * 4
    for source, multiplicity in enumerate(first_multiplicities):
        second_multiplicities[center_bijection[source]] = multiplicity
    options = {"pgl_gauge_center_positions_zero_based": []}
    left = compile_p2_divisor_evaluation_matrix(
        first, degree=4, multiplicities=first_multiplicities, chart_options=options
    )
    right = compile_p2_divisor_evaluation_matrix(
        second, degree=4, multiplicities=second_multiplicities, chart_options=options
    )
    rename = {
        sp.Symbol("q0_0"): sp.Symbol("q1_0"),
        sp.Symbol("q0_1"): sp.Symbol("q1_1"),
        sp.Symbol("q1_0"): sp.Symbol("q0_0"),
        sp.Symbol("q1_1"): sp.Symbol("q0_1"),
    }
    right_rows = {
        (
            descriptor["center_position_zero_based"],
            *descriptor["local_jet_multiindex"],
        ): row
        for descriptor, row in zip(right["jet_row_descriptors"], _matrix(right))
    }
    for descriptor, row in zip(left["jet_row_descriptors"], _matrix(left)):
        key = (
            center_bijection[descriptor["center_position_zero_based"]],
            *descriptor["local_jet_multiindex"],
        )
        transported = [value.subs(rename, simultaneous=True) for value in row]
        assert all(
            sp.cancel(actual - wanted) == 0
            for actual, wanted in zip(transported, right_rows[key])
        )


def _h10_label_map(packet: Mapping[str, Any]) -> dict[int, str]:
    old = _load(OLD_H10)
    walls = {
        center["source_contracted_wall_index_zero_based"]
        for center in packet["reverse_blowup_program"]["center_records"]
    }
    branches = [
        branch
        for branch in old["raw_polytope_extraction"]["branches"]
        if set(branch["contracted_wall_set_zero_based"]) == walls
    ]
    assert len(branches) == 1
    branch = branches[0]
    label_by_wall = dict(
        zip(
            branch["lexicographic_valid_contraction_order_zero_based"],
            branch["selected_point_labels_in_E_order"],
        )
    )
    return {
        center["center_position_zero_based"]: label_by_wall[
            center["source_contracted_wall_index_zero_based"]
        ]
        for center in packet["reverse_blowup_program"]["center_records"]
    }


def test_H10_calibration_rebuilt_by_direct_monomial_coefficients() -> None:
    """No call to either the new or legacy private matrix helper occurs."""

    packet = _packet_by_forest([2, 2, 2, 1, 1, 1, 1])
    labels = _h10_label_map(packet)
    proper_points = {
        "A": (1, 0, 0),
        "D": (0, 1, 0),
        "G": (0, 0, 1),
        "B": (1, 1, 1),
        "C": (1, 1, -1),
        "E": (-1, 1, 1),
        "F": (-1, 1, -1),
    }
    near_modes = {"a": "B", "d": "A", "g": "B"}
    label_to_center = {label: center for center, label in labels.items()}
    patches: dict[int, int] = {}
    modes: dict[int, str] = {}
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    for center, label in labels.items():
        if label in proper_points:
            point = proper_points[label]
            patch = next(index for index, value in enumerate(point) if value)
            normalized = [sp.Rational(value, point[patch]) for value in point]
            free = [index for index in range(3) if index != patch]
            patches[center] = patch
            substitutions[sp.Symbol(f"q{center}_0")] = normalized[free[0]]
            substitutions[sp.Symbol(f"q{center}_1")] = normalized[free[1]]
        else:
            modes[center] = near_modes[label]
            substitutions[sp.Symbol(f"t{center}")] = sp.Integer(0)
    compiled = compile_p2_cluster_chart(
        packet,
        {
            "pgl_gauge_center_positions_zero_based": [],
            "proper_affine_patch_by_center": patches,
            "free_chart_mode_by_center": modes,
        },
    )

    u, v, U, V = sp.symbols("u v U V")
    for system_key in ("minus_K", "minus_two_K"):
        record = compiled["evaluation_matrices"][system_key]
        columns = record["section_column_indices"]
        monomials = _monomials(columns)
        multiplicities = record["center_multiplicities_chronological"]
        expected: dict[tuple[int, int, int], list[sp.Expr]] = {}
        for center, label in labels.items():
            parent_label = label.upper()
            point = proper_points[parent_label]
            patch = next(index for index, value in enumerate(point) if value)
            normalized = [sp.Rational(value, point[patch]) for value in point]
            free = [index for index in range(3) if index != patch]
            coordinates = list(normalized)
            coordinates[free[0]] += u
            coordinates[free[1]] += v
            local = [
                sp.expand(
                    monomial.subs(
                        dict(zip((sp.Symbol("X"), sp.Symbol("Y"), sp.Symbol("Z")), coordinates))
                    )
                )
                for monomial in monomials
            ]
            if label.islower():
                parent = label_to_center[parent_label]
                parent_multiplicity = multiplicities[parent]
                if near_modes[label] == "A":
                    pulled = [
                        value.subs({u: U, v: U * V}, simultaneous=True)
                        for value in local
                    ]
                    shift = lambda left, right: (parent_multiplicity + left, right)
                else:
                    pulled = [
                        value.subs({u: U * V, v: V}, simultaneous=True)
                        for value in local
                    ]
                    shift = lambda left, right: (left, parent_multiplicity + right)
            else:
                pulled = local
                shift = lambda left, right: (left, right)
            for total in range(multiplicities[center]):
                for left in range(total + 1):
                    right = total - left
                    first_degree, second_degree = shift(left, right)
                    variables = (U, V) if label.islower() else (u, v)
                    expected[(center, left, right)] = [
                        _coefficient(value, *variables, first_degree, second_degree)
                        for value in pulled
                    ]
        _assert_exact_rows(record, expected, substitutions)


def _h11_label_map(packet: Mapping[str, Any]) -> dict[int, str]:
    old = _load(OLD_H11)
    walls = {
        center["source_contracted_wall_index_zero_based"]
        for center in packet["reverse_blowup_program"]["center_records"]
    }
    branches = [
        branch
        for branch in old["raw_extraction"]["branches"]
        if set(branch["target_raw_wall_indices_zero_based"]) == walls
    ]
    assert len(branches) == 1
    branch = branches[0]
    order = branch["lexicographic_valid_contraction_order_zero_based"]
    label_to_position = branch["raw_derived_canonical_program"][
        "canonical_label_to_exceptional_position_zero_based"
    ]
    label_by_wall = {order[position]: label for label, position in label_to_position.items()}
    return {
        center["center_position_zero_based"]: label_by_wall[
            center["source_contracted_wall_index_zero_based"]
        ]
        for center in packet["reverse_blowup_program"]["center_records"]
    }


def test_H11_quintic_calibration_rebuilt_without_legacy_matrix_helper() -> None:
    packet = _packet_by_forest([3, 3, 3, 2])
    labels = _h11_label_map(packet)
    label_to_center = {label: center for center, label in labels.items()}
    frames = {
        "A": [[1, 0], [0, 1]],
        "B": [[1, 0], [0, 1]],
        "C": [[0, 1], [1, 0]],
        "D": [[1, 0], [-2, 1]],
    }
    points = {"A": (0, 0), "B": (1, 0), "C": (0, 1), "D": (-1, 2)}
    second_slopes = {"A": sp.Integer(0), "C": sp.Symbol("a"), "D": sp.Symbol("b")}
    q_and_t: dict[sp.Symbol, sp.Expr] = {}
    frame_options: dict[int, Sequence[Sequence[int]]] = {}
    for center, label in labels.items():
        letter, level = label[0], int(label[1:])
        if level == 0:
            q_and_t[sp.Symbol(f"q{center}_0")] = points[letter][0]
            q_and_t[sp.Symbol(f"q{center}_1")] = points[letter][1]
            frame_options[center] = frames[letter]
        else:
            slope = sp.Integer(0) if level == 1 else second_slopes[letter]
            q_and_t[sp.Symbol(f"t{center}")] = slope
    quintic_rows = [
        row
        for row in packet["divisor_table"]
        if row["anti_bicanonical_boundary_weight"] > 0 and row["plane_degree"] == 5
    ]
    assert len(quintic_rows) == 1
    multiplicities = quintic_rows[0]["multiplicities_chronological"]
    record = compile_p2_divisor_evaluation_matrix(
        packet,
        degree=5,
        multiplicities=multiplicities,
        chart_options={
            "pgl_gauge_center_positions_zero_based": [],
            "proper_affine_patch_by_center": {
                center: 2 for center, label in labels.items() if label.endswith("0")
            },
            "proper_frame_by_center": frame_options,
        },
    )

    monomials = _monomials(record["section_column_indices"])
    expected: dict[tuple[int, int, int], list[sp.Expr]] = {}
    for center, label in labels.items():
        root_u, root_v = sp.symbols("U V")
        letter, level = label[0], int(label[1:])
        root_center = label_to_center[f"{letter}0"]
        point = points[letter]
        frame = frames[letter]
        root_x = point[0] + frame[0][0] * root_u + frame[0][1] * root_v
        root_y = point[1] + frame[1][0] * root_u + frame[1][1] * root_v
        root_local = [
            sp.expand(
                monomial.subs(
                    {
                        sp.Symbol("X"): root_x,
                        sp.Symbol("Y"): root_y,
                        sp.Symbol("Z"): 1,
                    }
                )
            )
            for monomial in monomials
        ]
        if level == 0:
            pulled = root_local
            pulled_u, pulled_v = root_u, root_v
            divisor_shift = 0
        elif level == 1:
            u1, v1 = sp.symbols("u1o v1o")
            slope1 = q_and_t[sp.Symbol(f"t{center}")]
            pulled = [
                value.subs(
                    {root_u: u1, root_v: u1 * (slope1 + v1)},
                    simultaneous=True,
                )
                for value in root_local
            ]
            pulled_u, pulled_v = u1, v1
            divisor_shift = multiplicities[root_center]
        else:
            level1_center = label_to_center[f"{letter}1"]
            slope1 = q_and_t[sp.Symbol(f"t{level1_center}")]
            slope2 = q_and_t[sp.Symbol(f"t{center}")]
            u2, v2 = sp.symbols("u2o v2o")
            pulled = [
                value.subs(
                    {
                        root_u: u2,
                        root_v: u2 * (slope1 + u2 * (slope2 + v2)),
                    },
                    simultaneous=True,
                )
                for value in root_local
            ]
            pulled_u, pulled_v = u2, v2
            divisor_shift = multiplicities[root_center] + multiplicities[level1_center]
        for total in range(multiplicities[center]):
            for left in range(total + 1):
                right = total - left
                expected[(center, left, right)] = [
                    _coefficient(
                        value,
                        pulled_u,
                        pulled_v,
                        divisor_shift + left,
                        right,
                    )
                    for value in pulled
                ]
    _assert_exact_rows(record, expected, q_and_t)


def test_PGL_gauge_and_finite_cover_are_scoped_not_promoted_to_existence() -> None:
    packet = _packet_by_forest([2, 2, 2, 1, 1, 1, 1])
    automatic = compile_p2_cluster_chart(packet)
    proper = [
        index
        for index, center in enumerate(packet["reverse_blowup_program"]["center_records"])
        if center["center_kind"] == "proper"
    ]
    selected = proper[:4]
    gauge = automatic["chart"]["pgl3_gauge"]
    assert gauge["selected_proper_center_positions_zero_based"] == selected
    targets = [sp.Matrix(row) for row in gauge["fixed_target_points_in_selected_order"]]
    assert all(
        sp.det(sp.Matrix.hstack(*(targets[index] for index in triple))) != 0
        for triple in itertools.combinations(range(4), 3)
    )
    expected_factors = {
        "Delta_" + "_".join(f"p{selected[index] + 1}" for index in triple)
        for triple in itertools.combinations(range(4), 3)
    }
    assert set(gauge["required_nonzero_triple_determinants"]) == expected_factors
    with pytest.raises(ValueError, match="cannot target a PGL3-gauged centre"):
        compile_p2_cluster_chart(
            packet, {"proper_affine_patch_by_center": {selected[0]: 2}}
        )

    boundary = automatic["chart"]["finite_chart_cover_boundary"]
    free_count = sum(
        center["center_kind"] == "free"
        for center in packet["reverse_blowup_program"]["center_records"]
    )
    assert boundary["possible_four_point_gauge_subsets"] == math.comb(len(proper), 4)
    assert boundary["coarse_cartesian_chart_count_for_this_gauge"] == (
        3 ** (len(proper) - 4) * 2**free_count
    )
    ungauged = compile_p2_cluster_chart(
        packet, {"pgl_gauge_center_positions_zero_based": []}
    )
    assert ungauged["chart"]["pgl3_gauge"]["full_gauge_used"] is False
    assert len(ungauged["chart"]["coordinate_variables"]) == 2 * len(proper) + free_count
    assert "chart-open nonemptiness" in ungauged["claims_not_made"]
    assert "matrix rank or condition independence" in ungauged["claims_not_made"]


def test_malformed_public_inputs_are_rejected_without_numeric_coercion() -> None:
    packet = _fixture_packet("satellite_center_two_parents")
    for degree in (True, 3.75, -1):
        with pytest.raises(ValueError):
            compile_p2_divisor_evaluation_matrix(
                packet, degree=degree, multiplicities=[1, 1, 1]
            )
    for columns in (
        [[1.9, 1.1, 1.0]],
        [[True, 1, 1]],
        [[1, 1, 1], [1, 1, 1]],
        [],
    ):
        with pytest.raises(ValueError):
            compile_p2_divisor_evaluation_matrix(
                packet,
                degree=3,
                multiplicities=[1, 1, 1],
                section_columns=columns,
            )
    for options in (
        {"proper_affine_patch_by_center": {"0": 2, "00": 1}},
        {"proper_frame_by_center": {0: [1, 2]}},
        {7: "forbidden"},
        {"proper_frame_by_center": {0: [[sp.Float("0.1"), 0], [0, 1]]}},
    ):
        with pytest.raises(ValueError):
            compile_p2_cluster_chart(packet, options)
    with pytest.raises(ValueError):
        compile_p2_cluster_chart([])  # type: ignore[arg-type]

    duplicate_parent = copy.deepcopy(packet)
    duplicate_parent["reverse_blowup_program"]["center_records"][2][
        "parent_center_positions_zero_based"
    ] = [1, 1]
    with pytest.raises(ValueError, match="duplicate proximity parents"):
        compile_p2_cluster_chart(duplicate_parent)

    late_parent = copy.deepcopy(packet)
    late_parent["reverse_blowup_program"]["center_records"][1][
        "parent_center_positions_zero_based"
    ] = [2]
    with pytest.raises(ValueError, match="proximity parent is malformed"):
        compile_p2_cluster_chart(late_parent)

    malformed_template = copy.deepcopy(packet)
    malformed_template["anti_canonical_evaluation_matrix_templates"]["minus_K"][
        "endpoint_class_coordinates"
    ][0] = 3.5
    with pytest.raises(ValueError, match="frozen B1.*degree"):
        compile_p2_cluster_chart(malformed_template)


def test_frozen_archive_pins_stable_B1_and_declared_claim_boundary() -> None:
    archive = _load(B2A_ARCHIVE)
    assert archive["all_checks_pass"] is True
    assert archive["aggregate"]["P2_orbit_representative_count"] == 42
    assert archive["aggregate"]["P2_labelled_endpoint_count"] == 630
    assert archive["aggregate"]["satellite_center_occurrence_count"] == 0
    assert archive["aggregate"]["rank_or_incidence_nonemptiness_asserted"] is False
    assert archive["sources"]["stage_B1_archive"]["sha256"] == _sha256(STAGE_B1)
    assert archive["sources"]["stage_B2a_core"]["sha256"] == _sha256(CORE)
    assert archive["sources"]["stage_B2a_certifier"]["sha256"] == _sha256(CERTIFIER)
    assert archive["sources"]["stage_B2a_tests"]["sha256"] == _sha256(OWNER_TEST)
    assert archive["sources"]["stage_B2a_report"]["sha256"] == _sha256(REPORT)
    assert archive["exact_boundary"]["no_incidence_nonemptiness_or_rank_claim"] is True
    assert archive["exact_boundary"][
        "no_irreducibility_hidden_wall_full_nef_or_same_P_claim"
    ] is True
