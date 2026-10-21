from __future__ import annotations

import json

from scripts.certify_coxeter_congruence_dessin import (
    DEFAULT_OUTPUT,
    ROOT,
    RSS_CEILING_BYTES,
    build,
)
from scripts.independent_coxeter_congruence_dessin_hostile_audit import (
    DEFAULT_OUTPUT as HOSTILE_OUTPUT,
    build as hostile_build,
)


def test_primary_archive_and_e10_level_five() -> None:
    certificate = build()
    assert json.loads(DEFAULT_OUTPUT.read_text()) == certificate
    assert certificate["all_checks_pass"]
    assert certificate["finite_order_level"]["rational_prime"] == 5
    assert certificate["finite_order_level"]["marked_involution_count"] == 10
    assert certificate["finite_order_level"]["order_2_pair_count"] == 36
    assert len(certificate["finite_order_level"]["order_3_pairs"]) == 9
    assert RSS_CEILING_BYTES == 1 << 29


def test_finite_nonrelation_tower_is_literal() -> None:
    tower = build()["finite_nonrelation_tower"]
    assert tower["all_prescribed_nonrelations_separated"]
    assert tower["prescribed_powers"] == list(range(1, 25))
    assert tower["difference_content_histogram"] == {"1": 24}
    assert tower["tower_exponent"] == 1


def test_height_reconstruction_is_separate_and_strict() -> None:
    height = build()["height_one_level_reconstruction"]
    assert height["logical_role"] == "separate from finite-order preservation"
    assert height["prime_ideal_norm"] == 5
    assert height["strict_norm_threshold"] == 4
    assert height["strict_inequality"]
    assert height["exact_gram_reconstructed"]


def test_independent_hostile_archive_and_failure_modes() -> None:
    audit = hostile_build()
    assert json.loads(HOSTILE_OUTPUT.read_text()) == audit
    assert audit["all_checks_pass"]
    assert audit["e10_level_5_replay"]["finite_pair_order_histogram"] == {"2": 36, "3": 9}
    assert set(audit["hostile_counterexamples"]) == {
        "relation_without_nonidentity_witness",
        "shallow_layer_without_valuation_bound",
        "single_layer_infinite_order_claim",
        "nonstrict_height_threshold",
        "unmarked_reflection_orientation",
    }


def test_report_keeps_stack_and_infinite_order_scope_exact() -> None:
    report = (ROOT / "results/coxeter_congruence_dessin_theorem.md").read_text()
    lowered = report.lower()
    assert "a neat cover has no inertia" in lowered
    assert "finite deck group" in lowered
    assert "no single congruence dessin certifies" in lowered
    assert "n\\mathfrak p>(2b)^d" in lowered
    assert "ideal " + "faces" not in lowered


def test_report_has_no_corrupt_control_bytes() -> None:
    payload = (ROOT / "results/coxeter_congruence_dessin_theorem.md").read_bytes()
    forbidden = [
        (offset, byte)
        for offset, byte in enumerate(payload)
        if byte < 32 and byte not in (9, 10)
    ]
    assert forbidden == []
    assert b"\r" not in payload
    assert b"\x08" not in payload
