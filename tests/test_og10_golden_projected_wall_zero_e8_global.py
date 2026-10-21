import json

from scripts.certify_og10_golden_projected_wall_zero_e8_global import (
    CERTIFICATE_PATH,
    build_certificate,
    verify_certificate,
)
from scripts.independent_og10_golden_projected_wall_zero_e8_global import run


def test_global_zero_e8_certificate_replays() -> None:
    archived = json.loads(CERTIFICATE_PATH.read_text())
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert archived["signature_admissible_pair_count"] == 26
    assert archived["sign_histogram"] == {"-1,-1": 13, "1,1": 13}
    assert archived["cutting_covector_count"] == 0


def test_independent_direct_covector_superset_has_no_cut() -> None:
    audit = run()
    assert audit["literal_superset_pair_count"] == 525
    assert audit["signature_admissible_count"] == 30
    assert audit["sign_histogram"] == {"-1,-1": 15, "1,1": 15}
    assert audit["cutting_records"] == []


def test_q_rectangle_and_square_unit_classes_are_exact() -> None:
    payload = build_certificate()
    assert payload["totally_positive_q_rectangle"] == [
        [1, 0],
        [1, 1],
        [2, -1],
        [2, 0],
        [2, 1],
        [3, -1],
        [3, 0],
    ]
    assert [len(group) for group in payload["projective_square_unit_classes"]] == [
        3,
        1,
        2,
        1,
    ]
