import json

from scripts.certify_h9_14_record0_a2_source29_triple123_cut import (
    EXPECTED_MATHEMATICAL_DIGEST,
    OUTPUT,
    _exact_structure,
    _near_seed_localization,
    verify_document,
)


def test_frozen_triple123_refinement() -> None:
    document = json.loads(OUTPUT.read_text())
    assert verify_document(document) == {
        "status": "PASS",
        "source29_status": "OPEN",
        "combined_refined_configuration_branch_count": 144,
        "strict_link_seed_found": False,
    }
    assert document["mathematical_digest"] == EXPECTED_MATHEMATICAL_DIGEST


def test_triple123_exact_fibre_budget() -> None:
    structure = _exact_structure()
    assert structure["base_link_g"] == [2, 3, 4]
    assert structure["allowed_six_fibre_count"] == 84
    assert structure["T_containing_candidate_A5_count"] == 36
    assert structure["noncontaining_killer_candidate_A5_count"] == 252
    assert structure["outside_fibre_budget"] == {
        "d0": 4,
        "d1": 2,
        "d2_intersecting_pairs": 1,
        "d2_disjoint_pairs": 0,
    }
    assert structure["zero_support_configuration_orbits"]["total"] == 8
    assert structure["active_support_configuration_orbits"]["total"] == 30
    assert structure["combined_refined_configuration_branch_count"] == 144


def test_near_seed_defect_is_exactly_two_outside_fibres() -> None:
    row = _near_seed_localization()
    assert row["T_containing_masks"] == [2094, 302]
    assert row["base_missing_six_fibres"] == 13
    assert [item["mask"] for item in row["selected_outside_killer_fibres"]] == [342, 1578]
    assert row["actual_missing_six_fibres"] == 15
    assert row["strict_max_missing_six_fibres"] == 14
    assert row["strict_defect"] == 1
