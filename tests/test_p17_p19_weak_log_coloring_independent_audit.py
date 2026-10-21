from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_p17_p19_weak_log_coloring_independent_audit import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/p17_p19_three_branch_property4_escape_independent_hostile_audit.md"
EXPECTED_CERTIFICATE_SHA256 = (
    "d9a4cc82dc1c4acaaff9568da4a928e07e014d10d081f0ba96cf12a82115215f"
)


def test_live_audit_passes_and_is_frozen() -> None:
    payload = build_certificate()
    assert payload["all_checks_pass"]
    assert payload["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["checks"]["producer_quartet_digests_pinned"]
    assert set(payload["producer_artifacts"]) == {
        "report",
        "certificate_generator",
        "archive",
        "tests",
    }


def test_every_physical_coloring_is_rejected() -> None:
    payload = build_certificate()
    assert [row["n"] for row in payload["rows"]] == [17, 18, 19]
    assert [row["sign_assignment_count"] for row in payload["rows"]] == [4, 32768, 32]
    assert [row["survivor_count"] for row in payload["rows"]] == [0, 0, 0]


def test_norm_two_facets_span_each_full_host() -> None:
    expected = {
        17: (18, 20),
        18: (19, 22),
        19: (20, 45),
    }
    for row in build_certificate()["rows"]:
        rank, count = expected[row["n"]]
        assert row["rank"] == rank
        assert row["norm_two_facet_count"] == count
        assert row["norm_two_span_rank"] == rank


def test_audit_does_not_reuse_the_old_canonical_orbits() -> None:
    payload = build_certificate()
    rendered = json.dumps(payload, sort_keys=True)
    assert "canonical_classes" not in rendered
    assert "MMP" in payload["scope"]["repair"]


def test_scope_keeps_maximality_and_host_quantifiers_open() -> None:
    not_proved = build_certificate()["scope"]["not_proved"]
    assert "absolute reflection maximality of P17, P18, or P19" in not_proved
    assert "exhaustion of different commensurable integral wall hosts or similarities" in not_proved


def test_report_records_repaired_pass_and_k3_scope() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for token in (
        "PASS in the exact standard-host/relative-maximal scope",
        "historical defect removed",
        "the norm-two walls already span",
        "zero survivors",
        "same-line/common-similarity contract",
        "absolute or maximal-arithmetic reflection maximality",
        "compact \\(\\Delta(2,4,6)\\) theorem",
    ):
        assert token in text
