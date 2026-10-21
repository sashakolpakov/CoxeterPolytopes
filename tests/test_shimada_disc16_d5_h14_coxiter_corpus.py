from __future__ import annotations

import json

from scripts.audit_shimada_disc16_d5_h14_coxiter_corpus import OUTPUT


def test_pinned_official_coxiter_corpus_has_no_literal_candidate_match() -> None:
    audit = json.loads(OUTPUT.read_text(encoding="utf-8"))
    snapshot = audit["corpus_snapshot"]
    scan = audit["scan"]
    decision = audit["decision"]

    assert audit["candidate"] == {
        "hyperbolic_dimension": 14,
        "facet_count": 60,
        "intrinsic_lattice": "U(2) + E8(-1) + A1(-1)^5",
    }
    assert snapshot["commit"] == "b800be48d8240fbabd2d1c9542e61054584ad7fe"
    assert snapshot["archive_sha256"] == (
        "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
    )
    assert scan["coxiter_file_count"] == 1214
    assert scan["parse_error_count"] == 0
    assert scan["global_maximum_facet_count"] == 50
    assert scan["h14_record_count"] == 3
    assert scan["h14_facet_histogram"] == {"17": 2, "20": 1}
    assert scan["candidate_dimension_and_facet_count_matches"] == []
    assert decision["exact_60_facet_h14_chamber_absent_from_snapshot"] is True
    assert decision["does_not_claim_underlying_lattice_absent"] is True
    assert decision["does_not_claim_global_literature_novelty"] is True
