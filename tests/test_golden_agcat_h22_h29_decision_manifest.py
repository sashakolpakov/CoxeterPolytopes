import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from coxeter_polytopes.classification.golden_agcat_manifest import (
    carrier_logic_records,
    ceil_integral_golden,
    cycle_box_record,
    floor_integral_golden,
    golden_closed_cauchy_edge_box,
    golden_coxeter_frame_pair_alphabet,
    golden_finite_angle_triangle_packets,
    golden_finite_pair_values,
    golden_frame_edge_box,
    rank_manifest_rows,
)
from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.certify_golden_agcat_h22_h29_decision_manifest import (
    OUTPUT,
    build_manifest,
    verify_manifest,
)


def test_exact_integral_golden_floor_and_ceiling_at_both_places() -> None:
    for a in range(-8, 9):
        for b in range(-8, 9):
            value = Golden(a, b)
            for embedding in (0, 1):
                floor = floor_integral_golden(value, embedding)
                ceiling = ceil_integral_golden(value, embedding)
                assert (value - floor).sign(embedding) >= 0
                assert (value - (floor + 1)).sign(embedding) < 0
                assert (value - ceiling).sign(embedding) <= 0
                assert (value - (ceiling - 1)).sign(embedding) > 0


def test_complete_golden_frame_pair_box_and_coxeter_filter() -> None:
    closed = golden_closed_cauchy_edge_box()
    box = golden_frame_edge_box()
    assert len(closed) == 6_878 == len(set(closed))
    assert len(box) == 6_877 == len(set(box))
    assert set(closed) - set(box) == {Golden(4)}
    assert all(
        value.sign(0) >= 0
        and (value - 62**2).sign(0) < 0
        and value.sign(1) >= 0
        and (value - 4).sign(1) < 0
        for value in box
    )

    finite = {row.product for row in golden_finite_pair_values()}
    assert finite == {Golden(0), Golden(1), Golden(2), 1 + PHI, Golden(3), 2 + PHI}
    alphabet = golden_coxeter_frame_pair_alphabet()
    assert len(alphabet) == 6_875
    assert sum((value - 4).sign(0) > 0 for value in alphabet) == 6_869
    assert Golden(4) not in alphabet
    assert {
        value for value in box if (value - 4).sign(0) < 0 and value not in alphabet
    } == {2 - PHI, 3 - PHI}


def test_cycle_rectangles_and_eight_rank_rows_are_explicit() -> None:
    records = tuple(cycle_box_record(length) for length in range(2, 31))
    assert len(records) == 29
    assert records[0]["selected_absolute_bound_strict"] == 62**2
    assert records[-1]["selected_absolute_bound_strict"] == 62**30
    assert records[0]["box_materialized"] is True
    assert all(record["rectangle_is_a_proved_superset"] for record in records)

    rows = rank_manifest_rows()
    assert [row["hyperbolic_dimension"] for row in rows] == list(range(22, 30))
    assert [row["native_rank"] for row in rows] == list(range(23, 31))
    assert sum(row["required_cycle_box_count"] for row in rows) == 204
    assert all(row["pipeline_status"]["empty"] == "not_proved" for row in rows)


def test_complete_finite_angle_triangle_cycle_tranche() -> None:
    packets = golden_finite_angle_triangle_packets()
    assert len(packets) == 10
    assert [row.labels for row in packets if row.admissible_lorentzian] == [
        (3, 3, 5),
        (3, 5, 5),
        (4, 4, 5),
        (5, 5, 5),
        (5, 10, 10),
    ]
    assert all(
        row.directed_three_cycle**2
        == row.edge_products[0] * row.edge_products[1] * row.edge_products[2]
        for row in packets
    )
    assert all(
        row.gram_determinant
        == 8 - 2 * sum(row.edge_products, start=Golden(0)) + 2 * row.directed_three_cycle
        for row in packets
    )
    assert all(
        row.selected_determinant_sign < 0 and row.conjugate_determinant_sign > 0
        for row in packets
        if row.admissible_lorentzian
    )


def test_AGCat_and_coherent_K3_logics_are_not_conflated() -> None:
    arrangement, k3 = carrier_logic_records()
    assert arrangement["exact_return"] is True
    assert arrangement["independent_AG_dimension_bound"] is False
    assert "what_it_does_not_prove" in arrangement
    assert k3["conditional_quadratic_dimension_bound"] == 8
    assert k3["essential_surjectivity_for_Qsqrt5_dimensions_22_29"] is False
    assert k3["can_be_used_to_exclude_raw_dimensions_22_29"] is False


def test_frozen_manifest_rebuilds_exactly_and_says_open() -> None:
    payload = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    verify_manifest(payload)
    assert payload == build_manifest()
    assert payload["verdict"]["existential_finiteness_of_reflection_maximal_slice_proved"]
    assert payload["verdict"]["enumerated_emptiness"] is False
    assert payload["verdict"]["nonexistence_in_dimensions_22_29_over_Qsqrt5"] == "OPEN"
    assert payload["finite_enumeration_debt"]["rank_cycle_box_applications"] == 204
    assert payload["finite_enumeration_debt"]["complete_full_host_Vinberg_transcripts_completed"] == 0
    bridge = payload["decision_target"]["reduction_to_all_compact_chambers"]
    assert bridge["maximal_overgroup_exists"]
    assert bridge["compactness_preserved"] and bridge["minimal_Vinberg_field_preserved"]
    scope = payload["decision_target"]["primary_enumeration_scope"]
    assert scope["every_candidate_full_host_must_be_completely_classified"] is True
    assert scope["direct_terminal_empty_set_must_be_empty"] is True
    assert scope["all_full_native_hosts_themselves_must_be_absent"] is False
    assert payload["finite_angle_triangle_cycle_tranche"]["raw_unordered_label_multisets"] == 35
    assert payload["finite_angle_triangle_cycle_tranche"]["admissible_lorentzian_packet_count"] == 5
    full_triangles = payload["full_triangle_cycle_tranche"]
    assert full_triangles["integral_golden_square_product_packet_count"] == 4_768_141
    assert full_triangles["conjugate_positive_packet_count"] == 2_633_885
    assert full_triangles["raw_selected_degenerate_packet_count"] == 1
    assert full_triangles["retained_selected_degenerate_packet_count"] == 0
    four_nodes = payload["four_node_at_most_one_dotted_tranche"]
    assert four_nodes["finite_only_cycle_compatible_orbits"] == 327
    assert four_nodes["one_dotted_cycle_compatible_orbits"] == 532_734
    assert four_nodes["finite_only_two_place_retained_orbits"] == 119
    assert four_nodes["one_dotted_two_place_retained_orbits"] == 180_402
    assert four_nodes["global_one_dotted_label3_path_family"][
        "admissible_dotted_value_counts"
    ] == [3_590, 3_584, 3_577, 3_571, 3_566, 3_561, 3_556, 3_552]
    assert payload["finite_enumeration_debt"][
        "local_four_node_at_most_one_dotted_tranches_completed"
    ] == 1
    two_dotted_trees = payload["four_node_exactly_two_dotted_tree_tranche"]
    assert two_dotted_trees["canonical_tree_packets_before_gram_gate"] == 471_865_955
    assert two_dotted_trees["two_place_retained_packets"] == 97_172_716
    assert two_dotted_trees["retained_by_topology"] == {
        "star": 25_023_281,
        "disjoint_dotted_path": 81,
        "adjacent_dotted_path": 72_149_354,
    }
    assert two_dotted_trees["independent_hostile_audit_status"].startswith("PASS")
    screens = payload["existing_normalized_diagonal_and_root_block_screens"]
    assert any(row["family"] == "<-a> plus S23 odd-Leech section" for row in screens)
    assert "golden_odd_Leech_section_H23" in payload["source_artifacts"]
    assert "golden_four_node_at_most_one_dotted_packets" in payload["source_artifacts"]
    assert payload["four_node_at_most_one_dotted_tranche"][
        "independent_hostile_audit_status"
    ] == "PASS"
    assert (
        "golden_four_node_at_most_one_dotted_independent_hostile_audit"
        in payload["source_artifacts"]
    )
    assert "golden_two_dotted_four_node_trees" in payload["source_artifacts"]
    assert (
        "golden_two_dotted_four_node_trees_independent_hostile_audit"
        in payload["source_artifacts"]
    )
