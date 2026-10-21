from __future__ import annotations

from pathlib import Path

import sympy as sp

from scripts.certify_projected_k3_universality_rank_saturation_no_go import (
    certify,
    direct_block,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/projected_k3_universality_rank_saturation_no_go.md"


def test_frozen_certificate_passes() -> None:
    payload = certify()
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())


def test_h19_saturates_the_k3_picard_rank() -> None:
    obstruction = certify()["h19_obstruction"]
    assert obstruction["input"]["hyperbolic_dimension"] == 19
    assert obstruction["input"]["trace_rank"] == 20
    assert obstruction["rank_saturation"]["forced_projector_kernel_rank"] == 0


def test_h19_has_two_incompatible_facet_square_classes() -> None:
    obstruction = certify()["h19_obstruction"]
    assert obstruction["input"]["facet_norm_counts"] == {
        "minus_two": 24,
        "minus_six": 12,
    }
    assert not obstruction["square_class_equations"]["three_is_a_rational_square"]
    assert obstruction["square_class_equations"]["three_adic_valuation"] % 2 == 1


def test_general_block_projector_formula() -> None:
    gram = sp.Matrix([[2, 1], [1, -2]])
    facets = sp.Matrix([[1, 0, 1], [0, 1, 1]])
    ambient, projector, schur = direct_block(gram, facets, 10)
    assert ambient.shape == (5, 5)
    assert projector.shape == (5, 5)
    assert schur.shape == (3, 3)
    assert projector * projector == projector
    assert projector.T * ambient == ambient * projector
    assert projector[:2, 2:] == facets / 10


def test_direct_block_root_identity() -> None:
    payload = certify()["direct_block_calibration"]
    assert payload["checks"]["root_identity_checked_on_729_pairs"]
    assert payload["checks"]["scale_exceeds_exact_row_sum_shell_bound"]
    assert payload["checks"]["stable_k3_embedding_inequality"]


def test_weak_mixed_divisor_scaled_trace_is_root_free_and_embeddable() -> None:
    payload = certify()["weak_mixed_divisor_calibration"]
    assert payload["scale"] > 2
    assert payload["rank"] <= 10
    assert "cannot equal -2" in payload["root_free_reason"]
    assert payload["all_checks_pass"]


def test_scope_does_not_promote_host_completeness_to_absolute_maximality() -> None:
    scope = certify()["scope"]
    assert "absolute maximality" in scope["counterexample_maximality"]
    assert "not asserted" in scope["counterexample_maximality"]
    assert "maximal-arithmetic-only" in scope["not_disproved"]


def test_report_preserves_theorem_and_scope_needles() -> None:
    report = REPORT.read_text()
    assert "Compact low-total-rank theorem" in report
    assert "Residue-shell completeness `RSC(s)`" in report
    assert "(a/b)^2=3" in report
    assert "does not by itself disprove a conjecture quantified only over" in report
    assert "absolutely maximal arithmetic reflection groups" in report
    assert "mixed-divisor K3 models" in report
    assert "Weak mixed-divisor theorem" in report
