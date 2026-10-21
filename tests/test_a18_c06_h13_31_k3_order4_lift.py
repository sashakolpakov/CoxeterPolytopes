from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a18_c06_h13_31_k3_order4_lift import build_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "a18_c06_h13_31_k3_order4_lift.json"
)
REPORT = ROOT / "results" / "a18_c06_h13_31_k3_order4_lift.md"
LEDGER = ROOT / "results" / "higher_dimensional_polytope_lineage_ledger.md"


def test_exact_order_four_k3_glue_replays_frozen_certificate() -> None:
    rebuilt = build_certificate()
    archived = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert archived["all_checks_pass"]
    assert archived["certificate_payload_sha256"] == (
        "c6dd715b4667021a7fa7a01189ed37b68ec20617f1ca58cf0cf141432278845e"
    )

    picard = archived["ambient_picard_lattice"]
    assert picard["invariants"] == [18, 4, 0]
    assert picard["selected_action_order"] == 2
    assert picard["selected_action_fixed_rank"] == 14

    glue = archived["discriminant_gluing"]
    assert glue["transcendental_lattice"] == "U(2)+U(2)"
    assert glue["action_square_is_minus_identity"]
    assert glue["actions_equivariant"]
    assert len(glue["quadratic_form_table"]) == 16
    assert glue["picard_discriminant_action"] == [
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
    ]

    k3 = archived["k3_lattice"]
    assert (k3["rank"], k3["signature"]) == (22, [3, 19])
    assert k3["even"] and k3["unimodular_determinant"] == -1
    assert k3["action_order"] == 4
    assert k3["action_characteristic_polynomial"] == (
        "(x-1)^14*(x+1)^4*(x^2+1)^2"
    )


def test_period_kahler_and_vinberg_handoffs_are_exact() -> None:
    archived = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    period = archived["period_and_torelli"]
    assert period["period_in_T_tensor_C"] == ["1", "sqrt(2)", "i", "i*sqrt(2)"]
    assert period["period_square"] == "0"
    assert period["period_positive_norm"] == "8*sqrt(2)"
    assert period["resulting_picard_lattice_is_exactly_N"]
    assert period["source_48_wall_chamber_is_preserved"]
    assert period["primitive_integral_polarization_scale"] == 2
    assert period["primitive_integral_polarization_square"] == "1540"
    assert "genuine non-symplectic automorphism of order four" in period[
        "global_torelli_conclusion"
    ]

    chamber = archived["hyperbolic_fixed_section"]
    assert (chamber["dimension"], chamber["facet_count"]) == (13, 31)
    assert chamber["finite_volume"] and not chamber["compact"]
    assert (chamber["finite_vertex_count"], chamber["ideal_vertex_count"]) == (
        3433,
        238,
    )
    assert chamber["elliptic_edge_count"] == 28766
    assert chamber["bad_edge_count"] == 0
    assert chamber["bad_parabolic_neighborhood_count"] == 0


def test_report_and_ledger_state_the_positive_result_without_overclaiming() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "genuine order-four K3 lift" in report
    assert "not a claim that" in report
    assert "does not turn the 48-wall mixed-root chamber" in report
    assert "worldwide priority" in report
    assert "A18-C2c-H13-31" in ledger
    assert "genuine K3 realization" in ledger
    assert "mixed-root source is not asserted" in ledger
