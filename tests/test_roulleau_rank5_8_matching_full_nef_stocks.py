from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.certify_roulleau_rank5_8_matching_full_nef_stocks import (
    OUTPUT,
    analyze,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/roulleau_rank5_8_matching_full_nef_stocks.md"
CANONICAL_SHA256 = "e05af6a969ad9b27fc3929c8fe499cf4975f1f066409bfbff46df58dc8d8a201"


@pytest.fixture(scope="module")
def archived() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_digest_and_live_exact_replay(archived: dict) -> None:
    body = dict(archived)
    claimed = body.pop("canonical_sha256")
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(encoded).hexdigest() == claimed == CANONICAL_SHA256
    assert analyze() == archived


def test_exact_target_header_denominator(archived: dict) -> None:
    assert archived["summary"] == {
        "stock_count": 17,
        "header_stock_counts": {
            "5,5": 2,
            "6,6": 3,
            "7,7": 2,
            "8,8": 3,
            "8,9": 5,
            "8,10": 2,
        },
        "decision_counts": {"finite_volume_noncompact": 17},
        "all_complete_full_nef": True,
        "all_literal_vinberg_finite_volume": True,
    }
    assert archived["primary_source"]["source_tar_sha256"] == (
        "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9"
    )
    assert archived["primary_source"]["source_tex_sha256"] == (
        "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
    )


def test_every_stock_is_complete_strict_and_literal_finite_volume(archived: dict) -> None:
    for case in archived["cases"]:
        matrix = case["matrix_certificate"]
        vinberg = case["literal_vinberg"]
        assert case["all_checks_pass"]
        assert "full_nef_K3" in case["classification"]
        assert matrix["all_checks_pass"]
        assert all(matrix["checks"].values())
        assert min(matrix["ample_wall_pairings"]) > 0
        assert vinberg["all_checks_pass"]
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["ideal_vertex_count"] > 0
        assert not vinberg["bad_edges"]
        assert not vinberg["bad_parabolic_neighborhoods"]


def test_report_freezes_source_and_claim_boundaries() -> None:
    report = REPORT.read_text(encoding="utf-8")
    assert "17 entries" in report
    assert "finite-volume noncompact" in report
    assert "complete irreducible (-2)-curve count" in report
    assert "not a claim to" in report
    assert CANONICAL_SHA256 in report
