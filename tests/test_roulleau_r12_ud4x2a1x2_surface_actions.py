import json

from scripts.certify_roulleau_r12_ud4x2a1x2_surface_actions import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_archive_recomputes_exactly() -> None:
    assert build_certificate() == _archive()


def test_two_actual_deck_actions_are_identified_exactly() -> None:
    payload = _archive()
    first = payload["first_double_plane_deck_action"]
    second = payload["second_F4_deck_action"]
    assert first["fixed_rank"] == 4
    assert len(first["nontrivial_cycles_by_source_name"]) == 29
    assert first["source_forced_pairs"]["B1_12"] == "C1_12"
    assert first["source_forced_pairs"]["B2_12"] == "C2_12"
    assert second["fixed_rank"] == 12
    assert second["permutation_zero_based"] == list(range(59))
    assert payload["general_surface"]["NS_image_order"] == 2


def test_p67_is_a_non_geometric_central_diagram_action() -> None:
    payload = _archive()
    p67 = payload["p67_central_diagram_action"]
    assert p67["fixed_rank"] == 11
    assert p67["central_in_full_weighted_group"]
    assert p67["status"].endswith("not_realized_by_Aut_of_the_general_K3")
    assert "not a quotient-nef chamber of the general K3" in p67[
        "strongest_geometric_wording"
    ]
    actual = {
        tuple(permutation)
        for permutation in payload["general_surface"]["NS_image_permutations_zero_based"]
    }
    assert tuple(p67["permutation_zero_based"]) not in actual


def test_all_source_and_lattice_checks_pass() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert all(payload["source"]["checks"].values())
    assert all(payload["checks"].values())
    assert all(payload["first_double_plane_deck_action"]["checks"].values())
    assert all(payload["p67_central_diagram_action"]["checks"].values())
