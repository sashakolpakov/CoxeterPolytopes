from scripts.verify_h9_14_record31_a1_orbit2_matrix_root1 import verify


def test_source2_matrix_root1_terminal_tranche() -> None:
    result = verify()
    assert result == {
        "root_class_id": 1,
        "root_class_count": 6,
        "incidence_target_orbit_count": 8,
        "topology_survivor_count": 8,
        "compact_label_survivor_count": 0,
        "generated_cpp_sha256": (
            "0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc"
        ),
        "scope": "one of six source-2 root classes only",
    }
