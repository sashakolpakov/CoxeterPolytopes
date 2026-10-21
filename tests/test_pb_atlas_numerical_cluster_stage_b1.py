from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from scripts.certify_pb_atlas_numerical_stage_b1 import (
    _reverse_wall_array,
    _swap_f0_rulings,
    build_certificate,
)
from scripts.freeze_pb_atlas_stage_b1_fixtures import build_fixtures
from scripts.pb_atlas_cluster_stage_b1 import compile_endpoint_program


ROOT = Path(__file__).resolve().parents[1]
STAGE_A_CORE = ROOT / "scripts/pb_atlas_contraction_core.py"
STAGE_A_ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json"
SCHEMA = ROOT / "schemas/pb-atlas-numerical-stage-b1.schema.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _fixture_by_id(fixture_id: str) -> dict:
    return next(row for row in _load(FIXTURES)["valid_cases"] if row["fixture_id"] == fixture_id)


def test_fixture_generator_is_exact_and_programs_validate() -> None:
    assert build_fixtures() == _load(FIXTURES)
    schema = _load(SCHEMA)
    for row in _load(FIXTURES)["valid_cases"]:
        packet = compile_endpoint_program(row["program"])
        jsonschema.validate(packet, schema)
        assert packet["all_checks_pass"] is True
        assert packet["reverse_blowup_program"]["center_kind_histogram"] == row[
            "expected_center_kind_histogram"
        ]


def test_free_and_satellite_proximity_matrices_are_exact() -> None:
    satellite = compile_endpoint_program(
        _fixture_by_id("satellite_center_two_parents")["program"]
    )
    reverse = satellite["reverse_blowup_program"]
    assert reverse["reverse_source_wall_order_zero_based"] == [2, 1, 0]
    assert reverse["proximity_incidence_matrix_Q_child_parent"] == [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
    ]
    assert reverse["strict_exceptional_to_total_matrix_P_chronological"] == [
        [1, -1, -1],
        [0, 1, -1],
        [0, 0, 1],
    ]
    assert reverse["strict_exceptional_gram_minus_P_P_transpose"] == [
        [-3, 0, 1],
        [0, -2, 1],
        [1, 1, -1],
    ]
    assert reverse["satellite_parent_intersection_audit"] == [
        {
            "satellite_center_position_zero_based": 2,
            "parent_center_positions_zero_based": [0, 1],
            "parent_strict_intersection_before_blowup": 1,
        }
    ]
    assert [row["center_kind"] for row in reverse["center_records"]] == [
        "proper",
        "free",
        "satellite",
    ]


def test_hostile_proximity_controls_and_label_routing_are_rejected() -> None:
    fixtures = _load(FIXTURES)
    for row in fixtures["invalid_cases"]:
        with pytest.raises(ValueError, match=row["expected_error_substring"]):
            compile_endpoint_program(row["program"])

    labeled = copy.deepcopy(fixtures["valid_cases"][0]["program"])
    labeled["known_surface_or_type_name"] = "forbidden"
    with pytest.raises(ValueError, match="unknown fields"):
        compile_endpoint_program(labeled)

    nontriangular = copy.deepcopy(fixtures["valid_cases"][0]["program"])
    first = nontriangular["sequential_contraction_order_zero_based"][0]
    model_rank = 1
    nontriangular["all_wall_coordinates_in_endpoint_blowup_basis"][first][
        model_rank + 1
    ] = -1
    with pytest.raises(
        ValueError,
        match="strict contracted-wall Gram does not match|not blowdown-triangular",
    ):
        compile_endpoint_program(nontriangular)


def test_divisor_proximity_RR_and_anticanonical_matrix_templates() -> None:
    packet = compile_endpoint_program(
        _fixture_by_id("satellite_center_two_parents")["program"]
    )
    consistent = packet["divisor_table"][4]
    hostile = packet["divisor_table"][5]
    assert consistent["plane_degree"] == 3
    assert consistent["multiplicities_chronological"] == [2, 1, 1]
    assert [row["proximity_excess"] for row in consistent["proximity_inequalities"]] == [
        0,
        0,
        1,
    ]
    assert consistent["all_proximity_inequalities_hold"] is True
    assert hostile["all_proximity_inequalities_hold"] is False

    classes = packet["class_matrices"]
    assert classes["basis_column_labels"] == ["H", "E_p1_total", "E_p2_total", "E_p3_total"]
    assert classes["rows"] == {
        "K": [-3, 1, 1, 1],
        "minus_K": [3, -1, -1, -1],
        "minus_two_K": [6, -2, -2, -2],
        "weighted_boundary_sum": [6, -2, -2, -2],
    }
    minus_k = packet["anti_canonical_evaluation_matrix_templates"]["minus_K"]
    minus_two_k = packet["anti_canonical_evaluation_matrix_templates"]["minus_two_K"]
    assert minus_k["evaluation_matrix_template"]["matrix_dimensions_rows_columns"] == [3, 10]
    assert minus_two_k["evaluation_matrix_template"]["matrix_dimensions_rows_columns"] == [9, 28]
    assert minus_k["riemann_roch"]["RR_virtual_projective_dimension"] == 6
    assert minus_two_k["riemann_roch"]["RR_virtual_projective_dimension"] == 18
    assert minus_k["evaluation_matrix_template"]["matrix_rank_not_asserted"] is True

    ruled = compile_endpoint_program(_fixture_by_id("ruled_free_center")["program"])
    ruled_minus_k = ruled["anti_canonical_evaluation_matrix_templates"]["minus_K"]
    ruled_minus_two_k = ruled["anti_canonical_evaluation_matrix_templates"]["minus_two_K"]
    assert ruled_minus_k["endpoint_class_coordinates"] == [2, 4]
    assert ruled_minus_k["riemann_roch"]["endpoint_section_count"] == 9
    assert ruled_minus_two_k["endpoint_class_coordinates"] == [4, 8]
    assert ruled_minus_two_k["riemann_roch"]["endpoint_section_count"] == 25


def test_commuting_orders_wall_permutations_and_F0_swap_are_invariant() -> None:
    left = compile_endpoint_program(
        _fixture_by_id("free_forest_two_commuting_roots")["program"]
    )
    right = compile_endpoint_program(
        _fixture_by_id("same_free_forest_commuted_order")["program"]
    )
    assert left["reverse_blowup_program"]["stage_A_contraction_order_zero_based"] != right[
        "reverse_blowup_program"
    ]["reverse_source_wall_order_zero_based"]
    assert left["invariant_cluster_profile_sha256"] == right[
        "invariant_cluster_profile_sha256"
    ]

    stage_a = _load(STAGE_A_ARCHIVE)
    first_program = stage_a["cases"][0]["certificate"]["P2_endpoint_orbits"][0]["program"]
    original = compile_endpoint_program(first_program)
    reversed_walls = compile_endpoint_program(_reverse_wall_array(first_program))
    assert original["invariant_cluster_profile_sha256"] == reversed_walls[
        "invariant_cluster_profile_sha256"
    ]

    f0 = next(
        orbit["program"]
        for case in stage_a["cases"]
        for orbit in case["certificate"]["hirzebruch_endpoint_orbits"]
        if orbit["program"]["hirzebruch_e"] == 0
    )
    assert compile_endpoint_program(f0)["invariant_cluster_profile_sha256"] == compile_endpoint_program(
        _swap_f0_rulings(f0)
    )["invariant_cluster_profile_sha256"]


def test_all_282_stage_A_programs_replay_to_frozen_hashes_and_schema() -> None:
    stage_a = _load(STAGE_A_ARCHIVE)
    frozen = _load(ARCHIVE)
    schema = _load(SCHEMA)
    frozen_records = {
        (case["calibration_alias_outside_compiler"], row["family"], row["representative_mask"]): row
        for case in frozen["cases"]
        for row in case["program_records"]
    }
    count = 0
    for case in stage_a["cases"]:
        for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
            for orbit in case["certificate"][family]:
                packet = compile_endpoint_program(orbit["program"])
                jsonschema.validate(packet, schema)
                record = frozen_records[(case["fixture_id"], family, orbit["representative_mask"])]
                assert record["packet_sha256"] == hashlib.sha256(
                    json.dumps(packet, sort_keys=True, separators=(",", ":")).encode("utf-8")
                ).hexdigest()
                assert record["divisor_table_sha256"] == hashlib.sha256(
                    json.dumps(packet["divisor_table"], sort_keys=True, separators=(",", ":")).encode(
                        "utf-8"
                    )
                ).hexdigest()
                count += 1
    assert count == 282


def test_frozen_stage_B1_archive_counts_and_stage_A_core_hash() -> None:
    archive = _load(ARCHIVE)
    stage_a = _load(STAGE_A_ARCHIVE)
    assert archive["status"] == "PASS_PB_ATLAS_NUMERICAL_CLUSTER_STAGE_B1"
    assert archive["all_checks_pass"] is True
    assert all(archive["checks"].values())
    assert all(archive["mutation_checks"].values())
    assert archive["aggregate"]["endpoint_program_count"] == 282
    assert archive["aggregate"]["endpoint_labelled_histogram"] == {
        "F_0": 1029,
        "F_1": 2586,
        "F_2": 564,
        "F_3": 44,
        "F_4": 12,
        "P2": 630,
    }
    assert archive["aggregate"]["representative_center_occurrence_histogram"] == {
        "free": 1600,
        "proper": 1116,
    }
    assert hashlib.sha256(STAGE_A_CORE.read_bytes()).hexdigest() == stage_a["sources"]["core"][
        "sha256"
    ]


def test_full_frozen_certificate_replay() -> None:
    assert build_certificate() == _load(ARCHIVE)
