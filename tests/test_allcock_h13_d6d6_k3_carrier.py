from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_allcock_h13_d6d6_k3_carrier import build


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/enumerations/allcock_h13_d6d6_k3_carrier.json"
REPORT = ROOT / "results/allcock_h13_d6d6_k3_carrier.md"


def _payload() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_archived_certificate_rebuilds_exactly() -> None:
    assert build() == _payload()


def test_corrected_figure3_and_full_wall_lattice_are_frozen() -> None:
    payload = _payload()
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    forcing = payload["rank_forcing_for_dashed_distances"]
    assert forcing["groebner_basis"] == ["x - 4", "y - 4"]
    assert forcing["unique_solution"] == {"x": 4, "y": 4}
    walls = payload["wall_lattice"]
    assert walls["determinant"] == -16
    assert walls["smith_diagonal"] == [1] * 10 + [2] * 4
    assert walls["inertia"] == [13, 1, 0]
    assert walls["wall_type_histogram"] == {
        "norm_2_div_1": 15,
        "norm_2_div_2": 4,
        "norm_4_div_2": 2,
    }
    assert walls["F_relation_histogram"] == {
        "0": 12, "2": 176, "3": 14, "4": 8
    }
    assert walls["Q_relation_histogram"] == {
        "0": 39, "2": 233, "3": 16, "4": 12
    }


def test_F_Q_and_every_finite_subtree_have_the_same_full_picard_carrier() -> None:
    payload = _payload()
    walls = payload["wall_lattice"]
    assert set(walls["F_delete_wall_span_indices"].values()) == {1}
    assert set(walls["Q_delete_wall_span_indices"].values()) == {1}
    k3 = payload["primitive_complement_and_K3_embedding"]
    assert k3["lattice_identification"]["sign_reversal"] == (
        "N=-K is isometric to U+D6(-1)+D6(-1)"
    )
    assert k3["primitive_K3_embedding"]["exact_Picard_moduli_dimension"] == 6
    assert "leaf copy" in walls["finite_subtree_span_argument"]


def test_ample_effective_and_equation_level_strengthenings_are_exact() -> None:
    payload = _payload()
    ample = payload["ample_intersection_framing"]
    assert ample["h_square_in_NS_minus_K"] == 1252
    assert min(ample["NS_intersections_with_oriented_F_wall_divisors"]) > 0
    assert ample["short_vector_enumeration"]["minimum"] == 4
    assert ample["short_vector_enumeration"]["nonzero_vectors_square_at_most_2"] == 0
    family = payload["weierstrass_family"]
    assert family["family"]["generic_fibers"] == "2 I_2^* + 8 I_1"
    assert family["family"]["generic_NS"] == "U+D6(-1)+D6(-1)"
    assert family["parameter_count"]["coefficient_image_rank"] == 8
    assert family["parameter_count"]["frozen_8_by_8_minor"] == 108
    assert family["specialization"]["two_torsion_groebner_basis"] == ["1"]


def test_novelty_and_nonclaims_are_explicit() -> None:
    payload = _payload()
    report = REPORT.read_text(encoding="utf-8")
    assert "worldwide publication priority not asserted" in payload["outcome"]["novelty_scope"]
    assert "the polytope is the full nef cone" in payload["excluded_claims"]
    assert "surface automorphism" in report
    assert "not a claim that either endpoint" in report
    assert "previously unknown" in report
    assert "full-Picard" in report
    assert "finite-subtree" in report
