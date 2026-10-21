from __future__ import annotations

import json

import pytest

from scripts.certify_raw_to_pb_atlas_stage0_a_bridge import (
    ARCHIVE,
    PB_FIXTURES,
    RAW_FIXTURES,
    build_certificate,
    candidate_to_pb_input,
)
from scripts.certify_raw_adjunction_log_stage0 import solve_raw_input


def _archive() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_frozen_bridge_is_exact_rebuild() -> None:
    assert build_certificate() == _archive()


def test_four_coble_payloads_are_derived_exactly() -> None:
    result = _archive()
    derived = [
        row for row in result["cases"] if row["routing"] == "PB_ATLAS_NUMERICAL_STAGE_A"
    ]
    expected = {
        case["calibration_alias"]: case["input"]
        for case in json.loads(PB_FIXTURES.read_text(encoding="utf-8"))["cases"]
    }
    assert len(derived) == 4
    for row in derived:
        assert row["derived_pb_atlas_input"] == expected[row["calibration_alias"]]
        assert all(row["checks"].values())


def test_raw_external_interface_contains_no_stage_a_metadata() -> None:
    fixtures = json.loads(RAW_FIXTURES.read_text(encoding="utf-8"))
    for case in fixtures["cases"]:
        assert set(case["raw_input"]) == {"schema", "picard_gram", "wall_vectors"}


def test_weak_and_maximal_negative_controls_do_not_enter_pb_atlas() -> None:
    rows = {row["case_id"]: row for row in _archive()["cases"]}
    assert rows["toric_h6_9"]["routing"] == "DIRECT_WEAK_ANTICANONICAL_BRANCH"
    assert (
        rows["blowup_F4_maximal_triangle_negative_control"]["routing"]
        == "REJECTED_AT_LOG_EQUATION"
    )


def test_converter_rejects_non_coble_candidates() -> None:
    fixtures = json.loads(RAW_FIXTURES.read_text(encoding="utf-8"))
    weak = next(case for case in fixtures["cases"] if case["case_id"] == "toric_h6_9")
    certificate = solve_raw_input(weak["raw_input"])
    with pytest.raises(ValueError, match="only the index-two Coble branch"):
        candidate_to_pb_input(weak["raw_input"], certificate["accepted_candidates"][0])
