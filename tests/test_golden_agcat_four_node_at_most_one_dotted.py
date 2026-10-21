from __future__ import annotations

from coxeter_polytopes.classification.golden_four_node_packets import (
    finite_only_canonical_patterns,
    golden_four_node_determinant,
    golden_four_node_packet_census,
    golden_one_dotted_path_rank_packets,
    golden_path_leading_determinant,
    golden_triangle_determinant,
    one_dotted_canonical_patterns,
)


def test_four_node_orbit_quotients_are_exact() -> None:
    finite = finite_only_canonical_patterns()
    dotted = one_dotted_canonical_patterns()
    assert len(finite) == 2_380
    assert len(set(finite)) == len(finite)
    assert len(dotted) == 2_080
    assert len(set(dotted)) == len(dotted)
    assert all(pattern[0] == 6 for pattern in dotted)


def test_determinant_expansion_on_tree_special_cases() -> None:
    # Star with products 1,2,3: det=16-4*(1+2+3)=-8.
    star = ((1, 0), (2, 0), (3, 0), (0, 0), (0, 0), (0, 0))
    assert golden_four_node_determinant(star) == (-8, 0)
    # Path products 1,2,3: det=16-4*(1+2+3)+1*3=-5.
    path = ((1, 0), (0, 0), (0, 0), (2, 0), (0, 0), (3, 0))
    assert golden_four_node_determinant(path) == (-5, 0)
    # The leading 012 restriction of the path has products 1,0,2.
    assert golden_triangle_determinant(path, (0, 1, 3)) == (2, 0)


def test_complete_at_most_one_dotted_census_counts_and_hashes() -> None:
    census = golden_four_node_packet_census()
    assert census.finite_value_count == 6
    assert census.dotted_value_count == 6_869
    assert census.finite_only_connected_orbit_count == 2_380
    assert census.one_dotted_complement_orbit_count == 2_080
    assert census.one_dotted_raw_orbit_count == 14_287_520
    assert census.finite_only_cycle_compatible_count == 327
    assert census.one_dotted_cycle_compatible_count == 532_734
    assert census.finite_only_conjugate_positive_count == 119
    assert census.one_dotted_conjugate_positive_count == 180_402
    assert census.finite_only_retained_count == 119
    assert census.finite_only_retained_selected_positive_count == 5
    assert census.finite_only_retained_selected_lorentzian_count == 114
    assert census.one_dotted_retained_count == 180_402
    assert census.one_dotted_retained_selected_lorentzian_count == 180_402
    assert census.one_dotted_selected_nonlorentzian_count == 0
    assert census.finite_only_selected_leading_triple_degenerate_count == 21
    assert census.finite_only_selected_full_degenerate_count == 3
    assert census.one_dotted_selected_leading_triple_degenerate_count == 0
    assert census.one_dotted_selected_full_degenerate_count == 0
    assert census.finite_only_conjugate_leading_triple_degenerate_count == 21
    assert census.finite_only_conjugate_full_degenerate_count == 3
    assert census.one_dotted_conjugate_leading_triple_degenerate_count == 0
    assert census.one_dotted_conjugate_full_degenerate_count == 0
    assert sum(row.count for row in census.strata) == 327 + 532_734
    assert census.finite_only_compatible_stream_sha256 == (
        "687762835339ac1f56fed84bcc41975e51da471abba4d97a7ab4cafc64590e7a"
    )
    assert census.finite_only_retained_stream_sha256 == (
        "a8a04a852ce054ebfa1668156200998263939d7417d7c437662130b2320c50e8"
    )
    assert census.one_dotted_compatible_stream_sha256 == (
        "0c988e10b93d804457dcfb9f238a832c546fd749b0a66751e67c875462153cc5"
    )
    assert census.one_dotted_retained_stream_sha256 == (
        "7de03d51631e744ba34c248a23a933da855e841256f7f27a33c1cf09b8a0b0a1"
    )


def test_overlapping_cycle_topologies_really_occur_and_are_retained() -> None:
    census = golden_four_node_packet_census()
    retained = {}
    for row in census.strata:
        conjugate_positive = (
            row.conjugate_leading_triple_sign > 0
            and row.conjugate_full_determinant_sign > 0
        )
        selected_allowed = (
            row.selected_leading_triple_sign < 0
            and row.selected_full_determinant_sign < 0
            if row.tranche == "exactly_one_dotted"
            else row.selected_full_determinant_sign < 0
            or (
                row.selected_leading_triple_sign > 0
                and row.selected_full_determinant_sign > 0
            )
        )
        if conjugate_positive and selected_allowed:
            retained[(row.tranche, row.topology)] = (
                retained.get((row.tranche, row.topology), 0) + row.count
            )
    assert retained[("exactly_one_dotted", "diamond_K4_minus_edge")] == 3_277
    assert retained[("exactly_one_dotted", "complete_K4")] == 1_168
    assert retained[("finite_only", "diamond_K4_minus_edge")] == 19
    assert retained[("finite_only", "complete_K4")] == 9


def test_global_rank_23_30_one_dotted_path_family() -> None:
    rows = golden_one_dotted_path_rank_packets()
    assert [row.native_rank for row in rows] == list(range(23, 31))
    assert [row.hyperbolic_dimension for row in rows] == list(range(22, 30))
    assert [row.admissible_dotted_value_count for row in rows] == [
        3_590, 3_584, 3_577, 3_571, 3_566, 3_561, 3_556, 3_552
    ]
    for row in rows:
        assert row.sample_dotted_product == (2, 2)
        assert row.sample_full_determinant == (
            2,
            -2 * (row.native_rank - 1),
        )
        assert row.sample_full_determinant == golden_path_leading_determinant(
            row.native_rank, row.sample_dotted_product
        )
