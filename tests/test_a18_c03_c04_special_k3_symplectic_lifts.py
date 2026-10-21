from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data"
    / "enumerations"
    / "a18_c03_c04_special_k3_symplectic_lifts.json"
)
REPORT = ROOT / "results" / "a18_c03_c04_special_k3_symplectic_lifts.md"
LEDGER = ROOT / "results" / "higher_dimensional_polytope_lineage_ledger.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_both_special_period_lifts_replay_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["case_count"] == 2
    assert archived["certificate_payload_sha256"] == (
        "44fbce0d83222fa8f1e868b6e03516f726a1b3bcf9592effa07b01c3183fac0f"
    )
    assert [case["id"] for case in archived["cases"]] == [
        "A18_4_0_C03",
        "A18_4_0_C04",
    ]


def test_nikulin_lattice_and_fixed_picard_lattice_are_exact() -> None:
    archived = _archived()
    for case in archived["cases"]:
        action = case["k3_lattice_action"]
        assert action["action_order"] == 2
        assert action["action_characteristic_polynomial"] == (
            "(x-1)^14*(x+1)^8"
        )
        assert action["fixed_lattice_rank_signature_determinant"] == [
            14,
            [3, 11],
            -256,
        ]
        assert action["anti_invariant_lattice"] == "E8(-2)"
        assert action["anti_invariant_rank"] == 8
        assert action["anti_invariant_determinant"] == 256
        assert action["anti_invariant_has_no_minus_two_vectors"]

        period = case["period_and_picard_lattice"]
        assert period["period_eigenvalue"] == "1"
        assert period["period_square"] == "0"
        assert period["ambient_N_is_contained"]
        assert period["fixed_NS_is_exactly_archived_N_fixed_lattice"]
        assert "symplectic Nikulin involution" in case["torelli_and_chamber"][
            "global_torelli_conclusion"
        ]


def test_special_picard_ranks_and_vinberg_handoffs() -> None:
    cases = {case["id"]: case for case in _archived()["cases"]}
    c03 = cases["A18_4_0_C03"]
    assert c03["period_and_picard_lattice"]["picard_rank"] == 20
    assert c03["period_and_picard_lattice"]["picard_determinant"] == -16
    assert c03["period_and_picard_lattice"]["fixed_NS_rank"] == 12
    assert c03["hyperbolic_fixed_section"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "compact": False,
        "criterion": (
            "Vinberg edge-endpoint criterion plus every irreducible affine "
            "component's full zero-neighborhood"
        ),
        "dimension": 11,
        "elliptic_edge_count": 8851,
        "facet_count": 28,
        "finite_vertex_count": 1252,
        "finite_volume": True,
        "ideal_vertex_count": 64,
    }

    c04 = cases["A18_4_0_C04"]
    assert c04["period_and_picard_lattice"]["picard_rank"] == 19
    assert c04["period_and_picard_lattice"]["picard_determinant"] == 16
    assert c04["period_and_picard_lattice"]["fixed_NS_rank"] == 11
    chamber = c04["hyperbolic_fixed_section"]
    assert (chamber["dimension"], chamber["facet_count"]) == (10, 28)
    assert (chamber["finite_vertex_count"], chamber["ideal_vertex_count"]) == (
        1087,
        74,
    )
    assert chamber["elliptic_edge_count"] == 6956
    assert chamber["finite_volume"] and not chamber["compact"]
    assert chamber["bad_edge_count"] == 0
    assert chamber["bad_parabolic_neighborhood_count"] == 0


def test_report_and_ledger_keep_special_period_and_nef_scope_explicit() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "Picard rank 20" in report
    assert "Picard rank 19" in report
    assert "symplectic Nikulin involution" in report
    assert "not the literal nef cone" in report
    assert "not new abstract" in report
    assert "A18-C2a-H11-28" in ledger and "Picard-rank-20" in ledger
    assert "A18-C2b-H10-28" in ledger and "Picard-rank-19" in ledger

