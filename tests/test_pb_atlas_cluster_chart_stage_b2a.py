from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import jsonschema
import pytest
import sympy as sp

import scripts.certify_pb_atlas_cluster_chart_stage_b2a as certifier
from scripts.pb_atlas_cluster_chart_stage_b2a import (
    compile_p2_cluster_chart,
    compile_p2_divisor_evaluation_matrix,
)
from scripts.pb_atlas_cluster_stage_b1 import compile_endpoint_program


ROOT = Path(__file__).resolve().parents[1]
STAGE_A = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
ARCHIVE = ROOT / "data/enumerations/pb_atlas_cluster_chart_stage_b2a.json"
SCHEMA = ROOT / "schemas/pb-atlas-cluster-chart-stage-b2a.schema.json"
CORE = ROOT / "scripts/pb_atlas_cluster_chart_stage_b2a.py"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _valid_fixture(fixture_id: str) -> dict:
    return next(
        row for row in _load(FIXTURES)["valid_cases"] if row["fixture_id"] == fixture_id
    )


def _packet_by_component_sizes(sizes: list[int]) -> dict:
    packets = certifier._all_p2_packets()
    return certifier._unique_packet_by_forest(packets, sizes)["packet"]


def test_one_selected_chart_replaces_every_B1_formal_anticanonical_entry() -> None:
    packet = _packet_by_component_sizes([2, 2, 2, 1, 1, 1, 1])
    compiled = compile_p2_cluster_chart(packet)
    jsonschema.validate(compiled, _load(SCHEMA))
    assert compiled["all_checks_pass"] is True
    assert compiled["source_contract"]["no_surface_type_or_calibration_name_entered"] is True
    for key in ("minus_K", "minus_two_K"):
        source = packet["anti_canonical_evaluation_matrix_templates"][key][
            "evaluation_matrix_template"
        ]
        target = compiled["evaluation_matrices"][key]
        assert target["matrix_dimensions_rows_columns"] == source[
            "matrix_dimensions_rows_columns"
        ]
        assert target["section_column_indices"] == source["section_column_indices"]
        assert target["jet_row_descriptors"] == source["jet_row_descriptors"]
        assert target["formal_jet_tokens_remaining"] is False
        assert all(
            set(entry) == {"numerator", "denominator"}
            for row in target["entries"]
            for entry in row
        )


def test_deterministic_gauge_and_explicit_ungauged_chart_are_both_available() -> None:
    packet = _packet_by_component_sizes([2, 2, 2, 1, 1, 1, 1])
    automatic = compile_p2_cluster_chart(packet)
    proper = [
        center["center_position_zero_based"]
        for center in packet["reverse_blowup_program"]["center_records"]
        if center["center_kind"] == "proper"
    ]
    assert automatic["chart"]["pgl3_gauge"][
        "selected_proper_center_positions_zero_based"
    ] == proper[:4]
    assert automatic["chart"]["pgl3_gauge"]["fixed_target_points_in_selected_order"] == [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 1],
    ]
    assert len(
        automatic["chart"]["pgl3_gauge"]["required_nonzero_triple_determinants"]
    ) == 4
    with pytest.raises(ValueError, match="cannot target a PGL3-gauged centre"):
        compile_p2_cluster_chart(
            packet, {"proper_affine_patch_by_center": {proper[0]: 2}}
        )

    explicit_gauge = [proper[1], proper[2], proper[3], proper[0]]
    explicit = compile_p2_cluster_chart(
        packet, {"pgl_gauge_center_positions_zero_based": explicit_gauge}
    )
    assert explicit["chart"]["pgl3_gauge"][
        "selected_proper_center_positions_zero_based"
    ] == explicit_gauge
    assert explicit["chart"]["pgl3_gauge"]["deterministic_rule"] == (
        "caller-supplied ordered proper-centre quadruple"
    )

    ungauged = compile_p2_cluster_chart(
        packet, {"pgl_gauge_center_positions_zero_based": []}
    )
    assert ungauged["chart"]["pgl3_gauge"]["full_gauge_used"] is False
    assert ungauged["chart"]["pgl3_gauge"]["deterministic_rule"] == (
        "gauge explicitly disabled by caller"
    )
    assert len(ungauged["chart"]["coordinate_variables"]) == 2 * len(proper) + 3
    assert ungauged["chart"]["finite_chart_cover_boundary"][
        "configurations_with_no_general_position_quadruple_require_the_ungauged_mode"
    ] is True


def test_free_and_satellite_substitutions_transport_both_exceptional_axes() -> None:
    packet = compile_endpoint_program(_valid_fixture("satellite_center_two_parents")["program"])
    compiled = compile_p2_cluster_chart(packet)
    charts = compiled["chart"]["center_charts"]
    assert [row["center_kind"] for row in charts] == ["proper", "free", "satellite"]
    assert charts[1]["chart_mode"] == "blowup_A"
    assert charts[1]["parent_local_substitution"] == {
        "u0": "u1",
        "v0": "t1*u1 + u1*v1",
    }
    assert charts[2]["chart_mode"] == "blowup_B"
    assert charts[2]["exceptional_axis_by_parent_center"] == {"0": "u", "1": "v"}
    assert charts[2]["parent_local_substitution"] == {"u1": "u2*v2", "v1": "v2"}

    generic = compile_p2_divisor_evaluation_matrix(
        packet, degree=4, multiplicities=[2, 1, 1]
    )
    assert generic["matrix_dimensions_rows_columns"] == [5, 15]
    assert generic["formal_jet_tokens_remaining"] is False
    local_symbols = {sp.Symbol(f"u{index}") for index in range(3)} | {
        sp.Symbol(f"v{index}") for index in range(3)
    }
    for row in generic["entries"]:
        for entry in row:
            denominator = sp.sympify(entry["denominator"])
            assert denominator.free_symbols.isdisjoint(local_symbols)

    # Independent monomial oracle through A followed by the forced satellite B
    # chart: X^2 Y^2 -> u1^2 v1^2 -> u2^2 v2^3 after divisions by 2 and 1.
    oracle = compile_p2_divisor_evaluation_matrix(
        packet,
        degree=4,
        multiplicities=[2, 1, 6],
        section_columns=[[2, 2, 0]],
    )
    substitutions = {
        sp.Symbol("q0_0"): 0,
        sp.Symbol("q0_1"): 0,
        sp.Symbol("t1"): 0,
    }
    satellite_entries = []
    for descriptor, row in zip(oracle["jet_row_descriptors"], oracle["entries"]):
        if descriptor["center_position_zero_based"] != 2:
            continue
        expression = sp.cancel(
            sp.sympify(row[0]["numerator"]) / sp.sympify(row[0]["denominator"])
        ).subs(substitutions)
        satellite_entries.append((descriptor["local_jet_multiindex"], sp.expand(expression)))
    assert dict((tuple(index), value) for index, value in satellite_entries)[(2, 3)] == 1
    assert sum(value != 0 for _index, value in satellite_entries) == 1


def test_generic_API_is_row_agnostic_and_validates_multiplicity_and_columns() -> None:
    packet = compile_endpoint_program(_valid_fixture("satellite_center_two_parents")["program"])
    matrix = compile_p2_divisor_evaluation_matrix(
        packet, degree=3, multiplicities=[0, 2, 3]
    )
    assert matrix["matrix_dimensions_rows_columns"] == [9, 10]
    assert [row["center_position_zero_based"] for row in matrix["jet_row_descriptors"]] == [
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
    ]
    with pytest.raises(ValueError, match="one per centre"):
        compile_p2_divisor_evaluation_matrix(packet, degree=3, multiplicities=[1])
    with pytest.raises(ValueError, match="homogeneous degree"):
        compile_p2_divisor_evaluation_matrix(
            packet,
            degree=3,
            multiplicities=[1, 1, 1],
            section_columns=[[1, 1, 0]],
        )
    for hostile_degree in (True, 3.75, -1):
        with pytest.raises(ValueError, match="degree.*must be"):
            compile_p2_divisor_evaluation_matrix(
                packet, degree=hostile_degree, multiplicities=[1, 1, 1]
            )
    for hostile_columns in (
        [[1.9, 1.1, 1.0]],
        [[True, 1, 1]],
        [[1, 1, 1], [1, 1, 1]],
    ):
        with pytest.raises(ValueError, match="exponent|distinct"):
            compile_p2_divisor_evaluation_matrix(
                packet,
                degree=3,
                multiplicities=[1, 1, 1],
                section_columns=hostile_columns,
            )


def test_duplicate_satellite_parents_are_rejected_at_the_public_boundary() -> None:
    packet = compile_endpoint_program(_valid_fixture("satellite_center_two_parents")["program"])
    hostile = copy.deepcopy(packet)
    hostile["reverse_blowup_program"]["center_records"][2][
        "parent_center_positions_zero_based"
    ] = [1, 1]
    with pytest.raises(ValueError, match="duplicate proximity parents"):
        compile_p2_cluster_chart(hostile)

    hostile_degree = copy.deepcopy(packet)
    hostile_degree["anti_canonical_evaluation_matrix_templates"]["minus_K"][
        "endpoint_class_coordinates"
    ][0] = 3.75
    with pytest.raises(ValueError, match="frozen B1 minus_K.*degree.*integer"):
        compile_p2_cluster_chart(hostile_degree)

    hostile_multiplicity = copy.deepcopy(packet)
    hostile_multiplicity["anti_canonical_evaluation_matrix_templates"]["minus_K"][
        "center_multiplicities_chronological"
    ][0] = True
    with pytest.raises(ValueError, match="frozen B1 minus_K.*multiplicity.*integer"):
        compile_p2_cluster_chart(hostile_multiplicity)

    with pytest.raises(ValueError, match="packet must be a mapping"):
        compile_p2_cluster_chart([])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="multiplicities must be an integer array"):
        compile_p2_divisor_evaluation_matrix(
            packet, degree=3, multiplicities=None  # type: ignore[arg-type]
        )


def test_F_e_packet_and_label_routing_options_are_rejected() -> None:
    ruled = compile_endpoint_program(_valid_fixture("ruled_free_center")["program"])
    with pytest.raises(ValueError, match="P2 endpoint packets only"):
        compile_p2_cluster_chart(ruled)

    plane = compile_endpoint_program(_valid_fixture("satellite_center_two_parents")["program"])
    with pytest.raises(ValueError, match="unknown Stage-B2a chart options"):
        compile_p2_cluster_chart(plane, {"surface_type_name": "forbidden"})
    with pytest.raises(ValueError, match="gauge centre.*integer"):
        compile_p2_cluster_chart(
            plane, {"pgl_gauge_center_positions_zero_based": [0.0, 1, 2, 3]}
        )
    with pytest.raises(ValueError, match="exact rational, not an approximate number"):
        compile_p2_cluster_chart(
            plane, {"proper_frame_by_center": {0: [[1.0, 0], [0, 1]]}}
        )
    with pytest.raises(ValueError, match="proper local frame is malformed"):
        compile_p2_cluster_chart(
            plane, {"proper_frame_by_center": {0: [1, 2]}}
        )
    with pytest.raises(ValueError, match="unknown Stage-B2a chart options"):
        compile_p2_cluster_chart(plane, {7: "forbidden"})
    with pytest.raises(ValueError, match="occurs more than once"):
        compile_p2_cluster_chart(
            plane, {"proper_affine_patch_by_center": {"0": 2, 0: 1}}
        )
    with pytest.raises(ValueError, match="nonnegative integer index"):
        compile_p2_cluster_chart(
            plane, {"proper_affine_patch_by_center": {"00": 2}}
        )
    with pytest.raises(ValueError, match="approximate number"):
        compile_p2_cluster_chart(
            plane,
            {"proper_frame_by_center": {0: [[sp.Float("0.1"), 0], [0, 1]]}},
        )


def test_H10_and_H11_old_plane_matrices_are_exact_entry_calibrations_only() -> None:
    packets = certifier._all_p2_packets()
    h10 = certifier._unique_packet_by_forest(packets, [2, 2, 2, 1, 1, 1, 1])
    h11 = certifier._unique_packet_by_forest(packets, [3, 3, 3, 2])
    h10_result = certifier._h10_calibration(h10["packet"])
    h11_result = certifier._h11_calibration(h11["packet"])
    assert h10_result["cubic_rows_equal_old_plane_matrix_as_multiset"] is True
    assert h10_result["sextic_rows_equal_old_plane_matrix_as_multiset"] is True
    assert h11_result[
        "quintic_rows_equal_old_plane_matrix_as_multiset_over_Q_a_b"
    ] is True
    assert h11_result["generic_divisor_API_used_not_minus_K_or_minus_two_K_special_case"] is True
    assert h10_result["compiler_used_surface_or_type_name"] is False
    assert h11_result["compiler_used_surface_or_type_name"] is False


def test_frozen_archive_scope_counts_hashes_and_claim_boundary() -> None:
    archive = _load(ARCHIVE)
    assert archive["status"] == "PASS_PB_ATLAS_P2_EXACT_CLUSTER_CHART_STAGE_B2A"
    assert archive["all_checks_pass"] is True
    assert all(archive["checks"].values())
    assert archive["aggregate"] == {
        "P2_labelled_endpoint_count": 630,
        "P2_orbit_representative_count": 42,
        "all_B1_minus_K_and_minus_two_K_formal_matrices_replaced": True,
        "center_occurrence_count": 433,
        "rank_or_incidence_nonemptiness_asserted": False,
        "satellite_center_occurrence_count": 0,
    }
    assert archive["exact_boundary"]["no_incidence_nonemptiness_or_rank_claim"] is True
    assert archive["exact_boundary"][
        "no_irreducibility_hidden_wall_full_nef_or_same_P_claim"
    ] is True
    assert archive["hostile_public_API_audit"]["all_checks_pass"] is True
    assert all(
        row["rejected_cleanly"] for row in archive["hostile_public_API_audit"]["probes"]
    )
    assert hashlib.sha256(CORE.read_bytes()).hexdigest() == archive["sources"][
        "stage_B2a_core"
    ]["sha256"]
