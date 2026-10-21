#!/usr/bin/env python3
"""Freeze the arbitrary-program PB-Atlas numerical cluster Stage B1."""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.freeze_pb_atlas_stage_b1_fixtures import build_fixtures  # noqa: E402
from scripts.pb_atlas_cluster_stage_b1 import (  # noqa: E402
    _json_hash,
    _pair,
    compile_endpoint_program,
)


STAGE_A_ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
STAGE_A_CORE = ROOT / "scripts/pb_atlas_contraction_core.py"
STAGE_A_SCHEMA = ROOT / "schemas/pb-atlas-numerical-stage-a.schema.json"
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
FIXTURE_SCRIPT = ROOT / "scripts/freeze_pb_atlas_stage_b1_fixtures.py"
CORE = ROOT / "scripts/pb_atlas_cluster_stage_b1.py"
SCHEMA_FILE = ROOT / "schemas/pb-atlas-numerical-stage-b1.schema.json"
CERTIFIER = ROOT / "scripts/certify_pb_atlas_numerical_stage_b1.py"
TEST_FILE = ROOT / "tests/test_pb_atlas_numerical_cluster_stage_b1.py"
REPORT = ROOT / "results/pb_atlas_numerical_cluster_stage_b1.md"
OUTPUT = ROOT / "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json"
SCHEMA = "coxeter.pb_atlas.numerical_cluster_stage_b1_calibration.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _endpoint_key(packet: Mapping[str, Any]) -> str:
    endpoint = packet["endpoint"]
    return "P2" if endpoint["kind"] == "P2" else f"F_{endpoint['hirzebruch_e']}"


def _matrix_shape(packet: Mapping[str, Any], system: str) -> str:
    rows, columns = packet["anti_canonical_evaluation_matrix_templates"][system][
        "evaluation_matrix_template"
    ]["matrix_dimensions_rows_columns"]
    return f"{rows}x{columns}"


def _compact_record(
    packet: Mapping[str, Any], *, family: str, mask: str, orbit_size: int
) -> dict[str, Any]:
    reverse = packet["reverse_blowup_program"]
    boundary_rows = [
        row for row in packet["divisor_table"] if row["anti_bicanonical_boundary_weight"]
    ]
    return {
        "family": family,
        "representative_mask": mask,
        "orbit_size": orbit_size,
        "endpoint": _endpoint_key(packet),
        "center_count": len(reverse["center_records"]),
        "wall_divisor_count": len(packet["divisor_table"]),
        "center_kind_histogram": reverse["center_kind_histogram"],
        "proximity_edge_count": len(reverse["proximity_edges_child_to_parent_zero_based"]),
        "dependency_depth_maximum": max(
            row["dependency_depth"] for row in reverse["center_records"]
        ),
        "boundary_component_count": len(boundary_rows),
        "boundary_components_with_nonnegative_proximity_consistent_multiplicities": sum(
            row["all_multiplicities_nonnegative"]
            and row["all_proximity_inequalities_hold"]
            for row in boundary_rows
        ),
        "minus_K_evaluation_shape": _matrix_shape(packet, "minus_K"),
        "minus_two_K_evaluation_shape": _matrix_shape(packet, "minus_two_K"),
        "divisor_table_sha256": _json_hash(packet["divisor_table"]),
        "invariant_cluster_profile_sha256": packet["invariant_cluster_profile_sha256"],
        "packet_sha256": _json_hash(packet),
    }


def _aggregate(records: list[dict[str, Any]]) -> dict[str, Any]:
    endpoint_programs = Counter(record["endpoint"] for record in records)
    endpoint_labelled = Counter()
    representative_centers = Counter()
    labelled_centers = Counter()
    minus_k_shapes = Counter()
    minus_two_k_shapes = Counter()
    for record in records:
        endpoint_labelled[record["endpoint"]] += int(record["orbit_size"])
        minus_k_shapes[record["minus_K_evaluation_shape"]] += 1
        minus_two_k_shapes[record["minus_two_K_evaluation_shape"]] += 1
        for kind, count in record["center_kind_histogram"].items():
            representative_centers[kind] += int(count)
            labelled_centers[kind] += int(count) * int(record["orbit_size"])
    return {
        "endpoint_program_count": len(records),
        "labelled_endpoint_count": sum(int(row["orbit_size"]) for row in records),
        "endpoint_program_histogram": dict(sorted(endpoint_programs.items())),
        "endpoint_labelled_histogram": dict(sorted(endpoint_labelled.items())),
        "representative_center_occurrence_histogram": dict(sorted(representative_centers.items())),
        "labelled_center_occurrence_histogram": dict(sorted(labelled_centers.items())),
        "wall_divisor_record_count": sum(int(row["wall_divisor_count"]) for row in records),
        "proximity_edge_occurrence_count": sum(int(row["proximity_edge_count"]) for row in records),
        "unique_invariant_cluster_profile_count": len(
            {row["invariant_cluster_profile_sha256"] for row in records}
        ),
        "minus_K_evaluation_shape_histogram": dict(sorted(minus_k_shapes.items())),
        "minus_two_K_evaluation_shape_histogram": dict(sorted(minus_two_k_shapes.items())),
    }


def _reverse_wall_array(program: Mapping[str, Any]) -> dict[str, Any]:
    """Transport a program through a hostile reversal of its source wall array."""

    result = copy.deepcopy(dict(program))
    coordinates = program["all_wall_coordinates_in_endpoint_blowup_basis"]
    endpoint_rows = program["all_wall_endpoint_image_coordinates"]
    wall_count = len(coordinates)
    old_to_new = [wall_count - 1 - index for index in range(wall_count)]
    result["all_wall_coordinates_in_endpoint_blowup_basis"] = list(reversed(coordinates))
    result["all_wall_endpoint_image_coordinates"] = list(reversed(endpoint_rows))
    result["contracted_wall_indices_zero_based"] = sorted(
        old_to_new[index] for index in program["contracted_wall_indices_zero_based"]
    )
    result["sequential_contraction_order_zero_based"] = [
        old_to_new[index] for index in program["sequential_contraction_order_zero_based"]
    ]
    for row in result["contraction_trace"]:
        row["source_wall_index_zero_based"] = old_to_new[
            row["source_wall_index_zero_based"]
        ]
    transported_boundary = []
    for row in result["boundary"]:
        row["wall_index_zero_based"] = old_to_new[row["wall_index_zero_based"]]
        transported_boundary.append(row)
    result["boundary"] = sorted(
        transported_boundary, key=lambda row: row["wall_index_zero_based"]
    )
    kind = result["endpoint_kind"]
    e = result["hirzebruch_e"]
    sorted_contracted = result["contracted_wall_indices_zero_based"]
    new_coordinates = result["all_wall_coordinates_in_endpoint_blowup_basis"]
    result["strict_contracted_wall_gram"] = [
        [_pair(new_coordinates[left], new_coordinates[right], kind, e) for right in sorted_contracted]
        for left in sorted_contracted
    ]
    return result


def _swap_f0_rulings(program: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(dict(program))
    if result["endpoint_kind"] != "F_e" or result["hirzebruch_e"] != 0:
        raise ValueError("the F0 mutation requires an F0 endpoint")
    result["model_basis_in_input_picard_coordinates"][0], result[
        "model_basis_in_input_picard_coordinates"
    ][1] = (
        result["model_basis_in_input_picard_coordinates"][1],
        result["model_basis_in_input_picard_coordinates"][0],
    )
    for row in result["all_wall_coordinates_in_endpoint_blowup_basis"]:
        row[0], row[1] = row[1], row[0]
    for row in result["all_wall_endpoint_image_coordinates"]:
        row[0], row[1] = row[1], row[0]
    for row in result["boundary"]:
        coordinates = row["endpoint_blowup_coordinates"]
        coordinates[0], coordinates[1] = coordinates[1], coordinates[0]
        endpoint = row["endpoint_image_coordinates"]
        endpoint[0], endpoint[1] = endpoint[1], endpoint[0]
    for ray in result["nef_rays_in_model_basis"]:
        ray[0], ray[1] = ray[1], ray[0]
    return result


def _synthetic_audit(fixtures: Mapping[str, Any]) -> dict[str, Any]:
    packets: dict[str, dict[str, Any]] = {}
    checks: dict[str, bool] = {}
    expected_by_id = {
        row["fixture_id"]: row["expected_center_kind_histogram"]
        for row in fixtures["valid_cases"]
    }
    for row in fixtures["valid_cases"]:
        packet = compile_endpoint_program(row["program"])
        packets[row["fixture_id"]] = packet
        checks[f"{row['fixture_id']}_histogram"] = (
            packet["reverse_blowup_program"]["center_kind_histogram"]
            == expected_by_id[row["fixture_id"]]
        )
    invalid = []
    for row in fixtures["invalid_cases"]:
        error = None
        try:
            compile_endpoint_program(row["program"])
        except ValueError as exception:
            error = str(exception)
        passed = error is not None and row["expected_error_substring"] in error
        checks[f"{row['fixture_id']}_rejected"] = passed
        invalid.append(
            {
                "fixture_id": row["fixture_id"],
                "expected_error_substring": row["expected_error_substring"],
                "actual_error": error,
                "rejected_as_expected": passed,
            }
        )
    left = packets["free_forest_two_commuting_roots"]
    right = packets["same_free_forest_commuted_order"]
    commuting_equal = (
        left["invariant_cluster_profile_sha256"]
        == right["invariant_cluster_profile_sha256"]
    )
    checks["commuting_reverse_blowup_orders_have_equal_invariant_profile"] = commuting_equal
    center_old_to_new = [1, 0, 2, 3]
    wall_old_to_new = [0, 1, 3, 2, 4, 5, 6]
    left_q = left["reverse_blowup_program"]["proximity_incidence_matrix_Q_child_parent"]
    right_q = right["reverse_blowup_program"]["proximity_incidence_matrix_Q_child_parent"]
    proximity_transports = all(
        left_q[child][parent]
        == right_q[center_old_to_new[child]][center_old_to_new[parent]]
        for child in range(4)
        for parent in range(4)
    )
    divisor_transport_checks = []
    for old_wall, new_wall in enumerate(wall_old_to_new):
        left_row = left["divisor_table"][old_wall]
        right_row = right["divisor_table"][new_wall]
        transported_multiplicities = [0] * 4
        for old_center, multiplicity in enumerate(left_row["multiplicities_chronological"]):
            transported_multiplicities[center_old_to_new[old_center]] = multiplicity
        divisor_transport_checks.append(
            left_row["endpoint_class_coordinates"] == right_row["endpoint_class_coordinates"]
            and transported_multiplicities == right_row["multiplicities_chronological"]
            and left_row["anti_bicanonical_boundary_weight"]
            == right_row["anti_bicanonical_boundary_weight"]
            and left_row["is_contracted_wall"] == right_row["is_contracted_wall"]
        )
    explicit_transport = proximity_transports and all(divisor_transport_checks)
    checks["commuting_order_explicit_DAG_and_divisor_transport_isomorphism"] = explicit_transport
    satellite = packets["satellite_center_two_parents"]
    checks["satellite_fixture_has_one_two_parent_center"] = (
        satellite["reverse_blowup_program"]["center_kind_histogram"].get("satellite")
        == 1
        and len(
            satellite["reverse_blowup_program"]["satellite_parent_intersection_audit"]
        )
        == 1
    )
    return {
        "valid_packets": packets,
        "invalid_controls": invalid,
        "mutation": {
            "commuting_order_left_profile_sha256": left["invariant_cluster_profile_sha256"],
            "commuting_order_right_profile_sha256": right["invariant_cluster_profile_sha256"],
            "equal": commuting_equal,
            "center_old_to_new_zero_based": center_old_to_new,
            "wall_old_to_new_zero_based": wall_old_to_new,
            "explicit_DAG_and_divisor_transport_isomorphism": explicit_transport,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    stage_a = _load(STAGE_A_ARCHIVE)
    fixtures = build_fixtures()
    frozen_fixtures = _load(FIXTURES)
    if fixtures != frozen_fixtures:
        raise ArithmeticError("the Stage-B1 fixture generator does not replay the frozen fixture")
    if _sha256(STAGE_A_CORE) != stage_a["sources"]["core"]["sha256"]:
        raise ArithmeticError("the frozen Stage-A core changed before Stage-B1 certification")

    cases = []
    all_records: list[dict[str, Any]] = []
    mutation_checks: dict[str, bool] = {}
    first_f0_program: Mapping[str, Any] | None = None
    for source_case in stage_a["cases"]:
        records = []
        first_program: Mapping[str, Any] | None = None
        for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
            for orbit in source_case["certificate"][family]:
                program = orbit["program"]
                if first_program is None:
                    first_program = program
                if first_f0_program is None and program["endpoint_kind"] == "F_e" and program["hirzebruch_e"] == 0:
                    first_f0_program = program
                packet = compile_endpoint_program(program)
                compact = _compact_record(
                    packet,
                    family=family,
                    mask=orbit["representative_mask"],
                    orbit_size=int(orbit["orbit_size"]),
                )
                records.append(compact)
                all_records.append(compact)
        if first_program is None:
            raise ArithmeticError("a Stage-A calibration case has no terminal endpoint program")
        original_packet = compile_endpoint_program(first_program)
        reversed_packet = compile_endpoint_program(_reverse_wall_array(first_program))
        wall_mutation_passes = (
            original_packet["invariant_cluster_profile_sha256"]
            == reversed_packet["invariant_cluster_profile_sha256"]
        )
        mutation_checks[f"{source_case['fixture_id']}_source_wall_reversal"] = wall_mutation_passes
        cases.append(
            {
                "calibration_alias_outside_compiler": source_case["fixture_id"],
                "source_stage_A_array_payload_sha256": source_case["certificate"]["input"][
                    "array_payload_sha256"
                ],
                "program_records": records,
                "aggregate": _aggregate(records),
                "source_wall_reversal_invariant_profile_check": wall_mutation_passes,
                "all_checks_pass": wall_mutation_passes,
            }
        )
    if first_f0_program is None:
        raise ArithmeticError("the calibrations supplied no F0 ruling-swap test")
    f0_original = compile_endpoint_program(first_f0_program)
    f0_swapped = compile_endpoint_program(_swap_f0_rulings(first_f0_program))
    mutation_checks["F0_unoriented_ruling_swap"] = (
        f0_original["invariant_cluster_profile_sha256"]
        == f0_swapped["invariant_cluster_profile_sha256"]
    )
    labeled = copy.deepcopy(dict(first_f0_program))
    labeled["known_surface_type"] = "forbidden routing label"
    rejected_label = False
    try:
        compile_endpoint_program(labeled)
    except ValueError as exception:
        rejected_label = "unknown fields" in str(exception)
    mutation_checks["compiler_rejects_surface_type_labels"] = rejected_label

    synthetic = _synthetic_audit(fixtures)
    checks = {
        "frozen_stage_A_core_hash_matches": True,
        "every_one_of_282_endpoint_orbit_programs_compiled": len(all_records) == 282,
        "all_four_stage_A_calibrations_pass": len(cases) == 4
        and all(row["all_checks_pass"] for row in cases),
        "all_program_level_mutations_pass": all(mutation_checks.values()),
        "synthetic_free_satellite_and_invalid_controls_pass": synthetic["all_checks_pass"],
        "no_satellite_was_silently_inferred_in_free_only_calibrations": all(
            record["center_kind_histogram"].get("satellite", 0) == 0
            for record in all_records
        ),
        "every_packet_hash_and_divisor_table_hash_frozen": all(
            len(record["packet_sha256"]) == 64
            and len(record["divisor_table_sha256"]) == 64
            for record in all_records
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Stage-B1 certificate checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "status": "PASS_PB_ATLAS_NUMERICAL_CLUSTER_STAGE_B1",
        "status_date": "2026-08-10",
        "scope": (
            "numerical reverse-blowup/proximity and divisor/evaluation-template Stage B1 only; "
            "no universal coordinate entries, incidence nonemptiness, bad factors, "
            "irreducibility, hidden-wall/full-nef proof, same-P return, or L4"
        ),
        "sources": {
            "stage_A_archive": {
                "path": str(STAGE_A_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha256(STAGE_A_ARCHIVE),
            },
            "stage_A_core_frozen_unchanged": {
                "path": str(STAGE_A_CORE.relative_to(ROOT)),
                "sha256": _sha256(STAGE_A_CORE),
            },
            "stage_A_output_schema": {
                "path": str(STAGE_A_SCHEMA.relative_to(ROOT)),
                "sha256": _sha256(STAGE_A_SCHEMA),
            },
            "stage_B1_core": {
                "path": str(CORE.relative_to(ROOT)),
                "sha256": _sha256(CORE),
            },
            "stage_B1_fixture": {
                "path": str(FIXTURES.relative_to(ROOT)),
                "sha256": _sha256(FIXTURES),
            },
            "stage_B1_fixture_generator": {
                "path": str(FIXTURE_SCRIPT.relative_to(ROOT)),
                "sha256": _sha256(FIXTURE_SCRIPT),
            },
            "stage_B1_output_schema": {
                "path": str(SCHEMA_FILE.relative_to(ROOT)),
                "sha256": _sha256(SCHEMA_FILE),
            },
            "stage_B1_certifier": {
                "path": str(CERTIFIER.relative_to(ROOT)),
                "sha256": _sha256(CERTIFIER),
            },
            "stage_B1_tests": {
                "path": str(TEST_FILE.relative_to(ROOT)),
                "sha256": _sha256(TEST_FILE),
            },
            "stage_B1_report": {
                "path": str(REPORT.relative_to(ROOT)),
                "sha256": _sha256(REPORT),
            },
        },
        "theorem": {
            "input": "one arbitrary arrays-only frozen Stage-A endpointProgram",
            "output": (
                "chronological reverse blowups, exact Enriques proximity matrix/DAG, every-wall "
                "endpoint degree or bidegree and multiplicity row, proximity excesses, exact RR "
                "data, and fully indexed -K/-2K evaluation matrix templates"
            ),
            "exact_identities": [
                "P = I - Q^T",
                "Gram(strict exceptionals) = -P P^T",
                "chi(O(D)) = 1 + (D^2-K.D)/2",
                "weighted boundary sum = -2K",
            ],
        },
        "cases": cases,
        "aggregate": _aggregate(all_records),
        "synthetic_hostile_audit": synthetic,
        "mutation_checks": mutation_checks,
        "checks": checks,
        "stage_B2_exact_remaining": [
            "choose universal P2/F_e gauges and generate proper/free/satellite chart coordinates",
            "replace every formal jet token by its polynomial evaluation entry",
            "generate incidence plus collision/Jacobian/irreducibility/Hilbert avoidance factors",
            "prove one surviving regular open yields the literal strong-Coble/full-nef carrier",
            "run hidden-wall exclusion and exact same-P replay",
        ],
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
