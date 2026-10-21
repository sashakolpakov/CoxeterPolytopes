from scripts.verify_h9_14_record31_a1_orbit2_matrix_root2 import verify


def test_source2_matrix_root2_partitioned_completion() -> None:
    result = verify()
    assert result["root_class_id"] == 2
    assert result["depth_two_class_count"] == 26
    assert result["terminal_fragment_count"] == 13
    assert result["raw_fragment_target_orbit_count"] == 7
    assert result["canonical_target_union_count"] == 3
    assert result["topology_survivor_count"] == 3
    assert result["compact_label_survivor_count"] == 0
