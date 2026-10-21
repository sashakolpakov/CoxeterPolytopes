from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.certify_bredon_adeg24_remaining_k3_stock_promotions import (
    OUTPUT,
    analyze,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/bredon_adeg24_remaining_k3_stock_promotions.md"
LEDGER = ROOT / "results/bredon_adeg24_remaining_k3_stock_promotion_ledger.md"
CANONICAL_SHA256 = "91a4867d958bb05af27ef08fbbc752779c65657ee7e94350fb1cb87db40302fa"
PROMOTIONS = {
    (
        "ADEG-N2-01-path-66::orientation-1",
        "R5-P4D4-NEF-H4-5",
    ): (12, 3, 198),
    (
        "ADEG-N3-01-path-633::orientation-0",
        "R6-UD4-NEF-H5-6",
    ): (6, 4, 630),
    (
        "ADEG-N3-04-star-6::orientation-0",
        "R6-U2D4-NEF-H5-6",
    ): (120, 4, 198),
    (
        "ADEG-N3-05-triangle-leaf-6::orientation-1",
        "R8-UA2x3-NEF-H7-10",
    ): (48, 4, 162),
    (
        "ADEG-N5-01-path-633336::orientation-0",
        "R10-UD4D4-NEF-H9-11",
    ): (72, 6, 504),
    (
        "ADEG-N5-02-path-triangle::orientation-0",
        "R8-UA2D4-NEF-H7-9",
    ): (12, 6, 504),
    (
        "ADEG-N6-01-diamond-tail::orientation-1",
        "R9-UD4A3-NEF-H8-10",
    ): (12, 7, 504),
    (
        "ADEG-N7-02-path-star::orientation-0",
        "R10-UD4D4-NEF-H9-11",
    ): (72, 8, 504),
}


@pytest.fixture(scope="module")
def archived() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_digest_and_live_replay(archived: dict) -> None:
    body = dict(archived)
    claimed = body.pop("canonical_sha256")
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(encoded).hexdigest() == claimed == CANONICAL_SHA256
    assert analyze() == archived


def test_complete_local_denominator_and_scope(archived: dict) -> None:
    summary = archived["summary"]
    assert summary["remaining_stable_packages"] == 30
    assert summary["local_complete_k3_nef_stocks"] == 55
    assert summary["local_stock_rank_distribution"] == {
        "5": 2,
        "6": 3,
        "7": 2,
        "8": 10,
        "9": 10,
        "10": 9,
        "11": 4,
        "12": 4,
        "13": 3,
        "14": 3,
        "15": 1,
        "16": 1,
        "17": 1,
        "18": 1,
        "19": 1,
    }
    assert len(archived["comparisons"]) == 30
    assert len(archived["stock_denominator"]) == 55
    assert len({stock["stock_id"] for stock in archived["stock_denominator"]}) == 55
    assert "not a worldwide or all-rank" in archived["scope"]["stock_denominator"]


def test_eight_exact_full_nef_lattice_returns(archived: dict) -> None:
    found = {}
    for comparison in archived["comparisons"]:
        for promotion in comparison["promotions"]:
            key = (comparison["package"], promotion["stock_id"])
            found[key] = (
                promotion["exact_weighted_gram_isomorphism_count"],
                promotion["triality_fixed_rank"],
                promotion["invariant_ample_square"],
            )
            assert promotion["integral_wall_lattice_isometry"]
            assert abs(promotion["wall_lattice_basis_change_determinant"]) == 1
            assert promotion["triality_order"] == 3
            assert promotion["triality_preserves_full_nef_wall_set"]
            assert promotion["invariant_ample_minimum_wall_pairing"] == 3
            assert promotion["coinvariant_minus_two_rootless_by_invariant_ample"]
            assert promotion["same_P_return"] == {
                "source_orbit_sum_gram_equals_original_folded_gram": True,
                "fixed_nef_section_has_exactly_original_facet_orbits": True,
                "hidden_nef_wall_count": 0,
                "level": "complete K3 nef cone plus divisor-lattice triality fixed section",
            }
            assert not promotion["triality_action_on_discriminant_is_trivial"]
            assert not promotion["geometric_order_three_K3_automorphism"]["certified"]
            assert promotion["complete_full_nef_source"]["literal_vinberg_decision"] \
                == "finite_volume_noncompact"
    assert found == PROMOTIONS
    assert archived["summary"]["exact_complete_stock_promotions"] == 8
    assert archived["summary"]["distinct_promoted_stocks"] == 7
    assert archived["summary"]["geometric_order_three_K3_automorphism_promotions"] == 0


def test_every_unpromoted_package_has_an_explicit_relative_gap(archived: dict) -> None:
    no_header = [
        row
        for row in archived["comparisons"]
        if row["status"] == "no_matching_stock_in_local_complete_denominator"
    ]
    rejected = [
        row
        for row in archived["comparisons"]
        if row["status"]
        == "matching_headers_all_rejected_by_lattice_or_exact_gram"
    ]
    assert len(no_header) == 7
    assert len(rejected) == 15
    assert all(not row["matching_header_stock_ids"] and row["gap"] for row in no_header)
    assert {row["package"] for row in no_header} == {
        "ADEG-N5-01-path-633336::orientation-1",
        "ADEG-N5-03-octagon-alternating::orientation-0",
        "ADEG-N5-03-octagon-alternating::orientation-1",
        "ADEG-N5-04-octagon-clustered::orientation-0",
        "ADEG-N5-04-octagon-clustered::orientation-1",
        "ADEG-N7-03-three-armed-path::orientation-0",
        "ADEG-N7-05-two-triangles-branch::orientation-1",
    }
    assert all(
        row["matching_header_stock_ids"]
        and not row["matching_lattice_stock_ids"]
        and not row["exact_weighted_gram_stock_ids"]
        and row["gap"]
        for row in rejected
    )
    assert archived["summary"]["packages_without_matching_local_header"] == 7
    assert archived["summary"]["packages_rejected_after_header_match"] == 15


def test_report_freezes_claim_boundary() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "TOOT: eight ADEG polytopes" in report
    assert "no new abstract polytope" in report.lower()
    assert "not a worldwide or all-rank K3 classification" in report
    assert "positive inverse calibrations" in report
    assert "not yet the requested flagship" in report
    assert "class-wide" in report and "pushback theorem" in report
    assert "zero**\norder-three automorphisms" in report
    assert "divisor-lattice fixed sections" in report
    assert "P = 8,    R = 15,    O = 7" in ledger
    assert ledger.count("| `ADEG-") == 30
    assert "`O` is not a nonexistence result" in ledger
