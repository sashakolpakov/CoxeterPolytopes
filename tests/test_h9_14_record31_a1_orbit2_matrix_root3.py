from scripts.verify_h9_14_record31_a1_orbit2_matrix_root3 import verify


def test_source2_matrix_root3_partitioned_completion() -> None:
    result = verify()
    assert result["root_class_id"] == 3
    assert result["depth_two_class_count"] == 19
    assert result["terminal_fragment_count"] == 19
    assert result["raw_fragment_target_orbit_count"] == 5
    assert result["canonical_target_union_count"] == 5
    assert result["topology_survivor_count"] == 5
    assert result["compact_label_survivor_count"] == 0
