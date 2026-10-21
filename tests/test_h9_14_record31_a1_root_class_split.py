from __future__ import annotations

from scripts.verify_h9_14_record31_a1_root_class_split import verify


def test_root_class_split_preserves_full_bounded_target_orbit_union() -> None:
    result = verify()
    assert result["unrestricted_transcript_core_identical"] is True
    assert result["bounded_full_search_node_count"] == 231_773
    assert result["bounded_split_search_node_count"] == 266_434
    assert result["bounded_target_orbit_union_identical"] is True
    assert result["split_internal_counts_nonadditive"] is True
    assert result["prepared_source4_root_class_range"] == [4, 6]
