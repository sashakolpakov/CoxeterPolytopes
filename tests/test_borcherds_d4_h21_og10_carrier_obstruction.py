from __future__ import annotations

import json

from scripts.certify_borcherds_d4_h21_og10_carrier_obstruction import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_polytope_first_scope_and_exact_replay_checks() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    accounting = payload["directional_accounting"]
    assert accounting["direction"] == (
        "POLYTOPE -> OG10 VARIETY -> POLYTOPE-CONE CRITERION"
    )
    assert "210-by-210 Coxeter Gram" in accounting["first_input"]
    assert "does not obstruct" in accounting["scope_warning"]


def test_unique_primitive_og10_embedding_orbit() -> None:
    payload = _payload()
    embedding = payload["primitive_OG10_embedding"]
    assert embedding["ambient_signature"] == [3, 21]
    assert embedding["ambient_determinant"] == -3
    assert embedding["ambient_smith_diagonal"] == [1] * 23 + [3]
    assert embedding["ambient_discriminant_q_mod_2"] == "4/3"
    assert embedding["L_embedding_smith_diagonal"] == [1] * 22
    assert embedding["ambient_wall_type_counts"] == {
        "square_-2_div_1": 42,
        "square_-4_div_1": 168,
    }

    exhaustion = payload["embedding_orbit_exhaustion"]
    assert exhaustion["compatible_complement_count"] == 1
    assert exhaustion["A2_2_integral_isometry_group_order"] == 12
    assert exhaustion["induced_two_torsion_permutation_group_order"] == 6
    assert "one primitive embedding orbit" in exhaustion["conclusion"]


def test_hidden_og10_wall_obstructs_complete_nef_promotion() -> None:
    payload = _payload()
    witness = payload["hidden_wall_certificate"]
    assert witness["hidden_root_in_L_basis"] == ["-1", "1"] + ["0"] * 20
    assert witness["hidden_root_square"] == -4
    assert witness["hidden_root_ambient_divisibility"] == 1
    assert witness["hidden_root_absent_from_210_facets"] is True
    assert witness["interior_square"] == 581
    assert witness["input_facet_pairing_minimum"] == "1"
    assert witness["input_facet_pairing_maximum"] == "421"
    assert witness["hidden_root_pairing"] == "0"
    assert "not a complete Kahler/nef cone" in witness["conclusion"]


def test_claim_is_carrier_specific_not_a_polytope_no_go() -> None:
    payload = _payload()
    accounting = payload["directional_accounting"]
    assert accounting["positive_result"] == (
        "the full-Picard OG10 period-carrier moduli problem is nonempty"
    )
    assert accounting["negative_result"] == (
        "the moduli problem requiring these 210 facets to be the complete "
        "Kahler/nef chamber is empty"
    )
    gate = payload["OG10_wall_divisor_gate"]
    assert gate["all_210_input_facets_are_wall_divisors"] is True
    assert gate["no_divisibility_three_argument"]["gram_invertible_mod_3"] is True
