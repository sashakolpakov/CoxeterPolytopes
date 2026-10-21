from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_alexeev_parked_h9_k3_action_carriers import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def _archive() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_builder_matches_frozen_archive() -> None:
    assert build_certificate() == _archive()


def test_four_original_actions_have_nikulin_k3_carriers() -> None:
    payload = _archive()
    rows = payload["strong_original_action_promotions"]
    assert [row["id"] for row in rows] == [
        "A14_8_0_C01",
        "A14_8_1_C02",
        "A14_10_1_C01",
        "A14_10_1_C02",
    ]
    for row in rows:
        assert row["original_ambient_action"]["action_intertwining_exact"] is True
        assert row["final_Picard_lattice"]["rank_signature_determinant"] == [
            18,
            [1, 17],
            -16,
        ]
        assert row["final_Picard_lattice"]["smith_diagonal_absolute"] == [
            *([1] * 16),
            4,
            4,
        ]
        coinvariant = row["final_Picard_lattice"]["coinvariant_lattice"]
        assert coinvariant["isometry_type"] == "E8(-2)"
        assert coinvariant["minus_two_rootless"] is True
        assert row["transcendental_and_K3_glue"][
            "Picard_primitive_and_exact_for_period"
        ] is True
        assert row["torelli_and_geometry"]["genuine_K3_automorphism"] is True
        assert row["torelli_and_geometry"]["automorphism_type"] == (
            "symplectic Nikulin involution"
        )
        assert row["torelli_and_geometry"]["effective_wall_claim"] is False
        assert row["torelli_and_geometry"]["imported_cone_equals_nef_claim"] is False


def test_fifth_obstruction_and_all_five_alternative_carriers_are_scoped() -> None:
    payload = _archive()
    obstruction = payload["fifth_route_obstruction"]
    assert obstruction["id"] == "A16_14_1_C01"
    assert obstruction["half_coinvariant_even"] is False
    assert obstruction["odd_norm_witness"]["norm"] == 3
    assert obstruction["symplectic_completion_to_E8_minus_2_impossible"] is True
    assert "not a no-surface-action theorem" in obstruction["scope"]

    alternatives = payload["pointwise_invariant_alternative_carriers"]
    assert len(alternatives) == 5
    assert all(
        row["all_Coxeter_reflections_extend_integrally"] is True
        and row["involution_fixes_full_Picard_pointwise"] is True
        and row["realizes_original_ambient_action"] is False
        and row["effective_wall_claim"] is False
        and row["imported_cone_equals_nef_claim"] is False
        for row in alternatives
    )
    phi = payload["surjectivity_Phi_C"]
    assert phi["original_action_genuine_invariant_divisor_carrier_count"] == 4
    assert phi["alternative_genuine_pointwise_invariant_divisor_carrier_count"] == 5
    assert phi["effective_or_nef_promotion_count"] == 0
    assert phi["surface_free_Coxeter_diagram_forces_carrier_type"] is False
