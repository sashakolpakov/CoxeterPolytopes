from __future__ import annotations

from scripts.verify_h9_14_record31_a1_orbit4_partitioned_transcript import verify


def test_source_orbit4_partitioned_incidence_no_go() -> None:
    result = verify()
    assert result["completed_original_root_class_ids"] == [1, 2, 3, 4, 5, 6]
    assert result["first_fragment_ridge_valid_target_orbit_count"] == 0
    assert result["complementary_fragment_ridge_valid_target_orbit_count"] == 0
    assert result["canonical_target_orbit_union_count"] == 0
    assert result["source_orbit_incidence_impossible"] is True
    assert result["fragment_internal_counts_combined"] is False
