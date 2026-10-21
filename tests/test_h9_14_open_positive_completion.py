from scripts.explore_h9_14_open_positive_a4 import mask
from scripts.explore_h9_14_open_positive_completion import (
    audit_compact_labels,
    build_problem,
)
from scripts.independent_h9_14_record1_a1_target0_label_audit import (
    replay as independent_target0_replay,
)


FIRST_RECORD1_ORBIT0_TARGET = (
    (1, 3, 4, 7, 8),
    (2, 3, 4, 7, 8),
    (3, 4, 5, 7, 8),
    (1, 3, 4, 8, 9),
    (0, 4, 5, 9, 10),
    (1, 4, 5, 9, 10),
    (0, 5, 6, 9, 10),
    (1, 2, 7, 8, 11),
    (1, 2, 8, 9, 11),
    (1, 4, 8, 9, 11),
    (1, 4, 9, 10, 11),
    (1, 6, 9, 10, 11),
    (0, 3, 5, 7, 12),
    (0, 2, 6, 7, 12),
    (2, 3, 7, 8, 12),
    (0, 5, 6, 10, 12),
    (2, 6, 7, 11, 12),
    (2, 6, 10, 11, 12),
    (6, 9, 10, 11, 12),
)


def test_record1_a1_source_orbit0_problem_shape() -> None:
    problem = build_problem(1, 1, 0)
    assert problem["fours"] == (mask((0, 3, 4, 5)),)
    assert len(problem["shadow"]) == 9
    assert problem["addition_count"] == 19
    assert problem["scheduled_facets"] == 158
    assert len(problem["candidates"]) == 1278
    assert len(problem["small"]) == 259
    assert len(problem["fives"]) == 634


def test_first_record1_a1_target_has_unrestricted_compact_label_no_go() -> None:
    problem = build_problem(1, 1, 0)
    target = tuple(mask(face) for face in FIRST_RECORD1_ORBIT0_TARGET)
    audit = audit_compact_labels(problem["edge_rows"], problem["fours"], target)
    assert audit["rank_nine_actual_vertex_count"] == 158
    assert [stage["passing_state_count"] for stage in audit["stages"]] == [
        108,
        84,
        8,
        0,
    ]
    assert audit["terminal_sparse_label_state_count"] == 0
    assert audit["unrestricted_compact_label_no_go"]


def test_first_record1_a1_target_independent_replay() -> None:
    audit = independent_target0_replay()
    assert audit["lanner_table_counts"] == {"4": 108, "5": 420}
    assert [stage["passing_state_count"] for stage in audit["stages"]] == [
        108,
        84,
        8,
        0,
    ]
    assert [stage["passing_state_sha256"] for stage in audit["stages"]] == [
        "857c6268544a6e22216a4068a1642a2e48730ff5288ae810cdaa398a9517f730",
        "de2dcf81429b02fe4970c82cc168428caddf019335da066bed4c228047e8b2e5",
        "a3cd5a07b59bf0c8e0b299bd334049f235dc243abded9be5f31d1ad6f011f270",
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    ]
    assert audit["unrestricted_compact_label_no_go"]
