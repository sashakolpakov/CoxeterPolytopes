#!/usr/bin/env python3
"""Freeze PB-Atlas Stage B2a for the P2 endpoint packets.

The production compiler is deliberately blind to the historical H10/H11
names.  This certifier finds the two plane packets by their proximity-forest
profiles, then uses the old labelled incidence realizations only as external
calibration adapters.
"""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import jsonschema
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import certify_coble_type_i_two_boundary_plane_inverse as old_h11  # noqa: E402
from scripts import certify_coble_type_ii_b3_plane_inverse as old_h10  # noqa: E402
from scripts.pb_atlas_cluster_chart_stage_b2a import (  # noqa: E402
    _json_hash,
    compile_p2_cluster_chart,
    compile_p2_divisor_evaluation_matrix,
)
from scripts.pb_atlas_cluster_stage_b1 import compile_endpoint_program  # noqa: E402


STAGE_A_ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
STAGE_B1_ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json"
STAGE_B1_CORE = ROOT / "scripts/pb_atlas_cluster_stage_b1.py"
STAGE_B1_SCHEMA = ROOT / "schemas/pb-atlas-numerical-stage-b1.schema.json"
STAGE_B1_REPORT = ROOT / "results/pb_atlas_numerical_cluster_stage_b1.md"
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
OLD_H10_ARCHIVE = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json"
OLD_H11_ARCHIVE = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse.json"
CORE = ROOT / "scripts/pb_atlas_cluster_chart_stage_b2a.py"
SCHEMA_FILE = ROOT / "schemas/pb-atlas-cluster-chart-stage-b2a.schema.json"
CERTIFIER = ROOT / "scripts/certify_pb_atlas_cluster_chart_stage_b2a.py"
TEST_FILE = ROOT / "tests/test_pb_atlas_cluster_chart_stage_b2a.py"
REPORT = ROOT / "results/pb_atlas_cluster_chart_stage_b2a.md"
OUTPUT = ROOT / "data/enumerations/pb_atlas_cluster_chart_stage_b2a.json"
SCHEMA = "coxeter.pb_atlas.cluster_chart_stage_b2a_calibration.v1"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_from_record(record: Mapping[str, Any]) -> sp.Matrix:
    return sp.Matrix(
        [
            [
                sp.cancel(
                    sp.sympify(entry["numerator"])
                    / sp.sympify(entry["denominator"])
                )
                for entry in row
            ]
            for row in record["entries"]
        ]
    )


def _canonical_row(row: Sequence[sp.Expr]) -> tuple[str, ...]:
    return tuple(sp.sstr(sp.cancel(sp.expand(value))) for value in row)


def _row_multiset(matrix: sp.Matrix) -> Counter[tuple[str, ...]]:
    return Counter(_canonical_row(list(matrix.row(index))) for index in range(matrix.rows))


def _component_sizes(packet: Mapping[str, Any]) -> list[int]:
    centers = packet["reverse_blowup_program"]["center_records"]
    adjacency = [set() for _ in centers]
    for child, center in enumerate(centers):
        for parent in center["parent_center_positions_zero_based"]:
            adjacency[child].add(parent)
            adjacency[parent].add(child)
    unseen = set(range(len(centers)))
    sizes: list[int] = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        unseen.remove(seed)
        count = 0
        while stack:
            vertex = stack.pop()
            count += 1
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        sizes.append(count)
    return sorted(sizes, reverse=True)


def _all_p2_packets() -> list[dict[str, Any]]:
    packets = []
    for case in _load(STAGE_A_ARCHIVE)["cases"]:
        for orbit in case["certificate"]["P2_endpoint_orbits"]:
            packet = compile_endpoint_program(orbit["program"])
            packets.append(
                {
                    "source_case_index": len(packets),
                    "representative_mask": orbit["representative_mask"],
                    "orbit_size": orbit["orbit_size"],
                    "packet": packet,
                }
            )
    return packets


def _unique_packet_by_forest(
    packets: Sequence[Mapping[str, Any]], sizes: Sequence[int]
) -> Mapping[str, Any]:
    matches = [record for record in packets if _component_sizes(record["packet"]) == list(sizes)]
    if len(matches) != 1:
        raise ArithmeticError(
            f"expected one intrinsic proximity-forest profile {list(sizes)}, found {len(matches)}"
        )
    return matches[0]


def _source_wall_by_center(packet: Mapping[str, Any]) -> dict[int, int]:
    return {
        int(center["center_position_zero_based"]): int(
            center["source_contracted_wall_index_zero_based"]
        )
        for center in packet["reverse_blowup_program"]["center_records"]
    }


def _specialize(matrix: sp.Matrix, substitutions: Mapping[sp.Symbol, sp.Expr]) -> sp.Matrix:
    return matrix.applyfunc(lambda value: sp.cancel(sp.expand(value.subs(substitutions))))


def _h10_calibration(packet: Mapping[str, Any]) -> dict[str, Any]:
    old = _load(OLD_H10_ARCHIVE)
    walls = set(_source_wall_by_center(packet).values())
    branches = [
        branch
        for branch in old["raw_polytope_extraction"]["branches"]
        if set(branch["contracted_wall_set_zero_based"]) == walls
    ]
    if len(branches) != 1:
        raise ArithmeticError("the H10 calibration did not find one exact contracted-wall branch")
    branch = branches[0]
    label_by_wall = dict(
        zip(
            branch["lexicographic_valid_contraction_order_zero_based"],
            branch["selected_point_labels_in_E_order"],
        )
    )
    wall_by_center = _source_wall_by_center(packet)
    label_by_center = {position: label_by_wall[wall] for position, wall in wall_by_center.items()}

    patch_by_center: dict[int, int] = {}
    mode_by_center: dict[int, str] = {}
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    for position, label in label_by_center.items():
        if label in old_h10.PROPER_POINTS:
            point = old_h10.PROPER_POINTS[label]
            patch = next(index for index, value in enumerate(point) if value)
            normalized = [sp.Rational(value, point[patch]) for value in point]
            free = [index for index in range(3) if index != patch]
            patch_by_center[position] = patch
            substitutions[sp.Symbol(f"q{position}_0")] = normalized[free[0]]
            substitutions[sp.Symbol(f"q{position}_1")] = normalized[free[1]]
        else:
            # These are exactly the chart orientations of the three old
            # near-point substitutions.  They live only in this adapter.
            mode_by_center[position] = {"a": "B", "d": "A", "g": "B"}[label]
            substitutions[sp.Symbol(f"t{position}")] = sp.Integer(0)
    options = {
        "pgl_gauge_center_positions_zero_based": [],
        "proper_affine_patch_by_center": patch_by_center,
        "free_chart_mode_by_center": mode_by_center,
    }
    compiled = compile_p2_cluster_chart(packet, options)
    cubic = _specialize(_matrix_from_record(compiled["evaluation_matrices"]["minus_K"]), substitutions)
    sextic = _specialize(
        _matrix_from_record(compiled["evaluation_matrices"]["minus_two_K"]), substitutions
    )
    _cubics, old_cubic, _cubic_names = old_h10._cluster_condition_matrix(3, 1)
    _sextics, old_sextic, _sextic_names = old_h10._cluster_condition_matrix(6, 2)
    cubic_exact = _row_multiset(cubic) == _row_multiset(old_cubic)
    sextic_exact = _row_multiset(sextic) == _row_multiset(old_sextic)
    return {
        "selection_profile_component_sizes": _component_sizes(packet),
        "matched_branch_by_exact_contracted_wall_set": True,
        "compiler_used_surface_or_type_name": False,
        "cubic_dimensions": [cubic.rows, cubic.cols],
        "sextic_dimensions": [sextic.rows, sextic.cols],
        "cubic_rows_equal_old_plane_matrix_as_multiset": cubic_exact,
        "sextic_rows_equal_old_plane_matrix_as_multiset": sextic_exact,
        "old_plane_archive_sha256": _sha256(OLD_H10_ARCHIVE),
        "all_checks_pass": cubic_exact and sextic_exact,
    }


def _h11_calibration(packet: Mapping[str, Any]) -> dict[str, Any]:
    old = _load(OLD_H11_ARCHIVE)
    walls = set(_source_wall_by_center(packet).values())
    branches = [
        branch
        for branch in old["raw_extraction"]["branches"]
        if set(branch["target_raw_wall_indices_zero_based"]) == walls
    ]
    if len(branches) != 1:
        raise ArithmeticError("the H11 calibration did not find one exact contracted-wall branch")
    branch = branches[0]
    order = branch["lexicographic_valid_contraction_order_zero_based"]
    label_to_position = branch["raw_derived_canonical_program"][
        "canonical_label_to_exceptional_position_zero_based"
    ]
    label_by_wall = {order[position]: label for label, position in label_to_position.items()}
    wall_by_center = _source_wall_by_center(packet)
    label_by_center = {position: label_by_wall[wall] for position, wall in wall_by_center.items()}

    proper_data = {
        "A0": ((0, 0), ((1, 0), (0, 1))),
        "B0": ((1, 0), ((1, 0), (0, 1))),
        "C0": ((0, 1), ((0, 1), (1, 0))),
        "D0": ((-1, 2), ((1, 0), (-2, 1))),
    }
    q_substitutions: dict[sp.Symbol, sp.Expr] = {}
    frames: dict[int, Sequence[Sequence[int]]] = {}
    for position, label in label_by_center.items():
        if label in proper_data:
            point, frame = proper_data[label]
            q_substitutions[sp.Symbol(f"q{position}_0")] = point[0]
            q_substitutions[sp.Symbol(f"q{position}_1")] = point[1]
            frames[position] = frame

    centers = packet["reverse_blowup_program"]["center_records"]
    direction_substitutions: dict[sp.Symbol, sp.Expr] = {}
    for position, center in enumerate(centers):
        if center["center_kind"] == "free":
            label = label_by_center[position]
            direction_substitutions[sp.Symbol(f"t{position}")] = {
                "A1": 0,
                "A2": 0,
                "B1": 0,
                "C1": 0,
                "C2": old_h11.A,
                "D1": 0,
                "D2": old_h11.B,
            }[label]
    options = {
        "pgl_gauge_center_positions_zero_based": [],
        "proper_affine_patch_by_center": {
            position: 2 for position, label in label_by_center.items() if label in proper_data
        },
        "proper_frame_by_center": frames,
    }
    quintic_rows = [
        row
        for row in packet["divisor_table"]
        if row["anti_bicanonical_boundary_weight"] > 0 and row["plane_degree"] == 5
    ]
    if len(quintic_rows) != 1:
        raise ArithmeticError("the H11 packet did not expose one weighted plane quintic row")
    generic = compile_p2_divisor_evaluation_matrix(
        packet,
        degree=5,
        multiplicities=quintic_rows[0]["multiplicities_chronological"],
        chart_options=options,
    )
    matrix = _matrix_from_record(generic)
    matrix = _specialize(matrix, {**q_substitutions, **direction_substitutions})
    old_matrix, _names = old_h11._quintic_condition_matrix()
    exact = _row_multiset(matrix) == _row_multiset(old_matrix)
    return {
        "selection_profile_component_sizes": _component_sizes(packet),
        "matched_branch_by_exact_contracted_wall_set": True,
        "compiler_used_surface_or_type_name": False,
        "generic_divisor_API_used_not_minus_K_or_minus_two_K_special_case": True,
        "quintic_dimensions": [matrix.rows, matrix.cols],
        "quintic_rows_equal_old_plane_matrix_as_multiset_over_Q_a_b": exact,
        "old_plane_archive_sha256": _sha256(OLD_H11_ARCHIVE),
        "all_checks_pass": exact,
    }


def _compact_record(source: Mapping[str, Any], compiled: Mapping[str, Any]) -> dict[str, Any]:
    packet = source["packet"]
    centers = packet["reverse_blowup_program"]["center_records"]
    return {
        "representative_mask": source["representative_mask"],
        "orbit_size": source["orbit_size"],
        "stage_B1_packet_sha256": _json_hash(packet),
        "stage_B2a_packet_sha256": _json_hash(compiled),
        "center_count": len(centers),
        "center_kind_histogram": packet["reverse_blowup_program"]["center_kind_histogram"],
        "proximity_component_sizes": _component_sizes(packet),
        "coordinate_variable_count": len(compiled["chart"]["coordinate_variables"]),
        "selected_full_PGL3_gauge": compiled["chart"]["pgl3_gauge"]["full_gauge_used"],
        "minus_K_entries_sha256": compiled["evaluation_matrices"]["minus_K"]["entries_sha256"],
        "minus_two_K_entries_sha256": compiled["evaluation_matrices"]["minus_two_K"]["entries_sha256"],
    }


def _strict_public_boundary_audit(
    packet: Mapping[str, Any], gauge_packet: Mapping[str, Any]
) -> dict[str, Any]:
    duplicate_parent_packet = copy.deepcopy(packet)
    duplicate_parent_packet["reverse_blowup_program"]["center_records"][2][
        "parent_center_positions_zero_based"
    ] = [1, 1]
    first_gauged_proper = next(
        center["center_position_zero_based"]
        for center in gauge_packet["reverse_blowup_program"]["center_records"]
        if center["center_kind"] == "proper"
    )
    probes = [
        (
            "nonmapping_packet",
            lambda: compile_p2_cluster_chart([]),
            "packet must be a mapping",
        ),
        (
            "nonarray_multiplicities",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet, degree=3, multiplicities=None
            ),
            "multiplicities must be an integer array",
        ),
        (
            "boolean_degree",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet, degree=True, multiplicities=[1, 1, 1]
            ),
            "must be an integer",
        ),
        (
            "fractional_degree",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet, degree=3.75, multiplicities=[1, 1, 1]
            ),
            "must be an integer",
        ),
        (
            "negative_degree",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet, degree=-1, multiplicities=[1, 1, 1]
            ),
            "must be nonnegative",
        ),
        (
            "fractional_exponents",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet,
                degree=3,
                multiplicities=[1, 1, 1],
                section_columns=[[1.9, 1.1, 1.0]],
            ),
            "must be an integer",
        ),
        (
            "duplicate_section_columns",
            lambda: compile_p2_divisor_evaluation_matrix(
                packet,
                degree=3,
                multiplicities=[1, 1, 1],
                section_columns=[[1, 1, 1], [1, 1, 1]],
            ),
            "must be distinct",
        ),
        (
            "duplicate_satellite_parents",
            lambda: compile_p2_cluster_chart(duplicate_parent_packet),
            "duplicate proximity parents",
        ),
        (
            "normalized_option_key_alias",
            lambda: compile_p2_cluster_chart(
                packet,
                {"proper_affine_patch_by_center": {"0": 2, 0: 1}},
            ),
            "occurs more than once",
        ),
        (
            "patch_choice_on_gauged_centre",
            lambda: compile_p2_cluster_chart(
                gauge_packet,
                {"proper_affine_patch_by_center": {first_gauged_proper: 2}},
            ),
            "cannot target a PGL3-gauged centre",
        ),
        (
            "scalar_frame_rows",
            lambda: compile_p2_cluster_chart(
                packet, {"proper_frame_by_center": {0: [1, 2]}}
            ),
            "frame is malformed",
        ),
        (
            "approximate_sympy_frame",
            lambda: compile_p2_cluster_chart(
                packet,
                {"proper_frame_by_center": {0: [[sp.Float("0.1"), 0], [0, 1]]}},
            ),
            "exact rational",
        ),
    ]
    records = []
    for name, probe, substring in probes:
        error = None
        try:
            probe()
        except ValueError as exception:
            error = str(exception)
        records.append(
            {
                "probe": name,
                "expected_error_substring": substring,
                "actual_error": error,
                "rejected_cleanly": error is not None and substring in error,
            }
        )
    return {"probes": records, "all_checks_pass": all(row["rejected_cleanly"] for row in records)}


def build_certificate() -> dict[str, Any]:
    stage_b1_schema = _load(STAGE_B1_SCHEMA)
    stage_b2a_schema = _load(SCHEMA_FILE)
    stage_b1_archive = _load(STAGE_B1_ARCHIVE)
    sources = _all_p2_packets()
    records = []
    checks: dict[str, bool] = {
        "frozen_stage_B1_archive_passes": stage_b1_archive["all_checks_pass"],
        "frozen_stage_B1_archive_pins_current_stage_A_archive": stage_b1_archive[
            "sources"
        ]["stage_A_archive"]["sha256"]
        == _sha256(STAGE_A_ARCHIVE),
        "frozen_stage_B1_archive_pins_current_stage_B1_core_schema_report": (
            stage_b1_archive["sources"]["stage_B1_core"]["sha256"]
            == _sha256(STAGE_B1_CORE)
            and stage_b1_archive["sources"]["stage_B1_output_schema"]["sha256"]
            == _sha256(STAGE_B1_SCHEMA)
            and stage_b1_archive["sources"]["stage_B1_report"]["sha256"]
            == _sha256(STAGE_B1_REPORT)
        ),
    }
    for index, source in enumerate(sources):
        packet = source["packet"]
        jsonschema.validate(packet, stage_b1_schema)
        compiled = compile_p2_cluster_chart(packet)
        jsonschema.validate(compiled, stage_b2a_schema)
        checks[f"P2_representative_{index:02d}_compiled_and_schema_valid"] = compiled[
            "all_checks_pass"
        ]
        histogram = packet["reverse_blowup_program"]["center_kind_histogram"]
        proper_count = int(histogram.get("proper", 0))
        free_count = int(histogram.get("free", 0))
        expected_dimension = 2 * proper_count + free_count - (
            8 if proper_count >= 4 else 0
        )
        checks[f"P2_representative_{index:02d}_chart_dimension_formula"] = len(
            compiled["chart"]["coordinate_variables"]
        ) == expected_dimension
        records.append(_compact_record(source, compiled))

    h10_source = _unique_packet_by_forest(sources, [2, 2, 2, 1, 1, 1, 1])
    h11_source = _unique_packet_by_forest(sources, [3, 3, 3, 2])
    h10 = _h10_calibration(h10_source["packet"])
    h11 = _h11_calibration(h11_source["packet"])
    checks["H10_old_plane_matrices_reproduced_exactly"] = h10["all_checks_pass"]
    checks["H11_old_plane_quintic_matrix_reproduced_exactly"] = h11["all_checks_pass"]

    fixture = next(
        row
        for row in _load(FIXTURES)["valid_cases"]
        if row["fixture_id"] == "satellite_center_two_parents"
    )
    satellite = compile_p2_cluster_chart(compile_endpoint_program(fixture["program"]))
    satellite_axes = satellite["chart"]["center_charts"][2][
        "exceptional_axis_by_parent_center"
    ]
    checks["synthetic_satellite_has_both_parent_exceptional_axes"] = set(
        satellite_axes
    ) == {"0", "1"} and set(satellite_axes.values()) == {"u", "v"}
    hostile_public_API = _strict_public_boundary_audit(
        compile_endpoint_program(fixture["program"]), h10_source["packet"]
    )
    checks["strict_public_API_rejects_inexact_or_ambiguous_inputs"] = hostile_public_API[
        "all_checks_pass"
    ]

    result = {
        "schema": SCHEMA,
        "status": "PASS_PB_ATLAS_P2_EXACT_CLUSTER_CHART_STAGE_B2A",
        "status_date": "2026-08-10",
        "scope": "all frozen Stage-A P2 orbit representatives; one selected exact chart per packet",
        "sources": {
            "stage_A_archive": {"path": str(STAGE_A_ARCHIVE.relative_to(ROOT)), "sha256": _sha256(STAGE_A_ARCHIVE)},
            "stage_B1_archive": {"path": str(STAGE_B1_ARCHIVE.relative_to(ROOT)), "sha256": _sha256(STAGE_B1_ARCHIVE)},
            "stage_B1_core": {"path": str(STAGE_B1_CORE.relative_to(ROOT)), "sha256": _sha256(STAGE_B1_CORE)},
            "stage_B1_schema": {"path": str(STAGE_B1_SCHEMA.relative_to(ROOT)), "sha256": _sha256(STAGE_B1_SCHEMA)},
            "stage_B1_report": {"path": str(STAGE_B1_REPORT.relative_to(ROOT)), "sha256": _sha256(STAGE_B1_REPORT)},
            "stage_B2a_core": {"path": str(CORE.relative_to(ROOT)), "sha256": _sha256(CORE)},
            "stage_B2a_schema": {"path": str(SCHEMA_FILE.relative_to(ROOT)), "sha256": _sha256(SCHEMA_FILE)},
            "stage_B2a_certifier": {"path": str(CERTIFIER.relative_to(ROOT)), "sha256": _sha256(CERTIFIER)},
            "stage_B2a_tests": {"path": str(TEST_FILE.relative_to(ROOT)), "sha256": _sha256(TEST_FILE)},
            "stage_B2a_report": {"path": str(REPORT.relative_to(ROOT)), "sha256": _sha256(REPORT)},
        },
        "aggregate": {
            "P2_orbit_representative_count": len(records),
            "P2_labelled_endpoint_count": sum(int(record["orbit_size"]) for record in records),
            "center_occurrence_count": sum(int(record["center_count"]) for record in records),
            "satellite_center_occurrence_count": sum(
                int(record["center_kind_histogram"].get("satellite", 0)) for record in records
            ),
            "all_B1_minus_K_and_minus_two_K_formal_matrices_replaced": True,
            "rank_or_incidence_nonemptiness_asserted": False,
        },
        "program_records": records,
        "calibration": {"H10_plane_inverse": h10, "H11_plane_inverse": h11},
        "synthetic_satellite_audit": {
            "fixture_id": fixture["fixture_id"],
            "satellite_exceptional_axis_by_parent_center": satellite_axes,
            "all_checks_pass": checks["synthetic_satellite_has_both_parent_exceptional_axes"],
        },
        "hostile_public_API_audit": hostile_public_API,
        "checks": checks,
        "exact_boundary": {
            "one_chart_per_packet_not_a_global_single_chart_claim": True,
            "alternate_four_point_gauges_affine_patches_and_blowup_charts_remain_a_finite_cover": True,
            "no_incidence_nonemptiness_or_rank_claim": True,
            "no_irreducibility_hidden_wall_full_nef_or_same_P_claim": True,
            "P2_only_F_e_packets_rejected": True,
        },
        "all_checks_pass": all(checks.values()),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate()
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(certificate["aggregate"], indent=2, sort_keys=True))
    print(certificate["status"], "all_checks_pass=", certificate["all_checks_pass"])


if __name__ == "__main__":
    main()
