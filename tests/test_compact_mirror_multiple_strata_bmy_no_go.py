from __future__ import annotations

import json

from scripts.certify_compact_mirror_multiple_strata_bmy_no_go import (
    OUTPUT,
    build_certificate,
    half_weight_record,
    reduced_blowup_record,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_live_certificate_matches_archive() -> None:
    assert build_certificate() == _archive()


def test_reduced_multiple_points_have_the_wrong_sign_for_feedback() -> None:
    for b in range(2, 30):
        record = reduced_blowup_record(b)
        assert record["log_c2_change"] == 0
        assert record["log_c1_square_change"] == -(b - 2) ** 2
        assert record["bmy_defect_change"] == (b - 2) ** 2
        assert record["nonnegative"] is True


def test_half_weight_dihedral_values_and_resolution_signs() -> None:
    records = [half_weight_record(m) for m in (2, 3, 4)]
    assert [record["langer_local_eorb"] for record in records] == [
        "1/4",
        "1/16",
        0,
    ]
    assert [record["global_bmy_point_term"] for record in records] == [
        "1/4",
        "3/16",
        0,
    ]
    assert [
        record["reduced_exceptional_resolution"]["bmy_defect_change"]
        for record in records
    ] == ["1/4", "1/16", 0]
    assert [
        record["crepant_resolution"]["bmy_defect_change"]
        for record in records
    ] == ["9/4", "9/16", 0]
    assert all(record["all_relevant_terms_nonnegative"] for record in records)


def test_scope_keeps_the_uniform_determinant_degree_gate_open() -> None:
    payload = _archive()
    assert "class-uniform noncircular bound for U_det" in payload["remaining_gate"]
    assert "do not provide a uniformly negative" in payload["conclusion"]
    assert payload["all_checks_pass"] is True

