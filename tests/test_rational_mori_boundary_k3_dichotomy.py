import json
from pathlib import Path

import pytest

from scripts.certify_rational_mori_boundary_k3_dichotomy import (
    ARCHIVE,
    adjusted_intersection,
    build_certificate,
    canonical_intersection,
    certify,
    dimension_ceiling,
    picard_ceiling,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_mori_boundary_k3_dimension_dichotomy.md"


def test_symbolic_mori_alphabets():
    assert canonical_intersection(-1) == -1
    assert canonical_intersection(-2) == 0
    assert canonical_intersection(-4) == 2
    assert adjusted_intersection(-4, -4) == 0
    assert adjusted_intersection(-2, 0) == 0
    assert adjusted_intersection(-1, 2) == 0


def test_exact_ceilings_and_char_two_guard():
    assert picard_ceiling("weak_anticanonical", 0) == 10
    assert dimension_ceiling("weak_anticanonical", 0) == 9
    assert picard_ceiling("index_two_boundary", 0) == 20
    assert dimension_ceiling("index_two_boundary", 0) == 19
    assert picard_ceiling("index_two_boundary", 3) == 22
    assert dimension_ceiling("index_two_boundary", 3) == 21
    with pytest.raises(ValueError):
        picard_ceiling("index_two_boundary", 2)


def test_certificate_and_controls():
    payload = build_certificate()
    certify(payload)
    assert payload["all_checks_pass"]
    assert payload["controls"]["toric"]["passes"]
    assert payload["controls"]["toric"]["alphabet_replay"]["all_checks_pass"]
    assert payload["controls"]["toric"]["alphabet_replay"]["square_histogram"] == {
        "-2": 6,
        "-1": 3,
    }
    assert (
        payload["controls"]["toric"]["alphabet_replay"]["wall_gram_rank"] == 7
    )
    assert len(payload["controls"]["coble"]) == 4
    assert all(row["passes"] for row in payload["controls"]["coble"])
    assert all(
        row["alphabet_replay"]["all_checks_pass"]
        and row["alphabet_replay"]["checks"]["physical_walls_span_N1"]
        and row["alphabet_replay"]["checks"][
            "source_certifies_complete_mori_generation"
        ]
        for row in payload["controls"]["coble"]
    )
    assert [
        row["alphabet_replay"]["other_generator_pair_histogram"]
        for row in payload["controls"]["coble"]
    ] == [
        {"(-2,0)": 12, "(-1,2)": 17},
        {"(-2,0)": 10, "(-1,2)": 14},
        {"(-2,0)": 10, "(-1,2)": 45},
        {"(-2,0)": 24, "(-1,2)": 16},
    ]
    assert [
        row["alphabet_replay"]["physical_wall_gram_rank"]
        for row in payload["controls"]["coble"]
    ] == [11, 12, 15, 18]
    assert payload["controls"]["coble"][0]["alphabet_replay"][
        "archived_boundary_scale_to_physical_curve"
    ] == 2
    assert all(
        row["alphabet_replay"]["archived_boundary_scale_to_physical_curve"] == 1
        for row in payload["controls"]["coble"][1:]
    )


def test_frozen_archive_matches_builder():
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == build_certificate()


def test_report_scope_and_clean_bytes():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    assert "This does **not** prove that every reflection-maximal" in text
    assert "B=-2K_S" in text
    assert "PB-Atlas" in text
    assert "generators automatically span" in text
    assert "integral divisor \\(2K_S+B\\) is numerically trivial" in text
    assert "Locally along \\(B\\) it has equation" in text
    assert "The cover is connected" in text
    assert "No K3-cover conclusion is asserted in characteristic" in text
    assert "entire closed Mori cone" in text
    assert "Completeness of the Mori generators is imported" in text
    assert "stored first normal is \\(b=-K_S=B/2\\)" in text
    assert "seven source pins, and five control rows" in text
