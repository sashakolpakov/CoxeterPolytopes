from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a16_h11_29_ud8a1x2_nef_effectivity import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_exact_marking_uses_the_positive_Roulleau_component(
    certificate: dict,
) -> None:
    marking = certificate["exact_marking_and_sign"]
    assert marking["chosen_sign"] == "+"
    assert marking["initial_D136_degree"] == 581903
    assert marking["opposite_sign_D136_degree"] == -581903
    assert all(marking["checks"].values())
    assert all(certificate["published_nef_chamber_stock"]["checks"].values())


def test_explicit_87_reflection_word_reaches_the_closed_nef_chamber(
    certificate: dict,
) -> None:
    transport = certificate["explicit_weyl_transport"]
    assert transport["word_length"] == len(transport["firing_trace"]) == 87
    assert transport["W_determinant"] == -1
    assert transport["final_nef_carrier_coordinates"] == [
        178, 267, 418, 209, 360, 302, 271, 240, 58, 31, 120, 89
    ]
    assert transport["final_nef_carrier_square"] == 37418
    assert transport["final_nef_carrier_D136_degree"] == 3352
    assert transport["final_pairings_with_A1_through_A19"] == [
        0, 62, 0, 0, 0, 27, 0, 0, 124, 178, 0, 0, 178, 240, 294,
        62, 418, 116, 240,
    ]
    assert transport["zero_wall_curve_labels"] == [
        "A1", "A3", "A4", "A5", "A7", "A8", "A11", "A12"
    ]
    assert all(transport["checks"].values())


def test_seventeen_minus_two_mirrors_have_exact_effectivity_decisions(
    certificate: dict,
) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_two_effective_decompositions"]
    assert mirrors["square_minus_two_effective_count"] == len(records) == 17
    assert mirrors["square_minus_two_irreducible_count"] == 1
    assert mirrors["square_minus_two_reducible_count"] == 16
    assert [
        (row["facet_zero_based"], row["irreducible_curve_label"])
        for row in records
        if row["is_itself_an_irreducible_curve"]
    ] == [(23, "A16")]
    assert min(row["total_component_multiplicity"] for row in records) == 1
    assert max(row["total_component_multiplicity"] for row in records) == 18
    assert max(row["maximum_component_coefficient"] for row in records) == 4
    assert all(all(row["decomposition_checks"].values()) for row in records)
    assert all(mirrors["effective_checks"].values())


def test_H_plus_D136_is_primitive_ample_and_coxeter_strict(
    certificate: dict,
) -> None:
    simultaneous = certificate["explicit_simultaneous_ample_interior"]
    ample = simultaneous["B1_equals_H_plus_D136"]
    assert ample["coordinates"] == [
        180, 272, 427, 213, 370, 314, 286, 259, 67, 40, 126, 89
    ]
    assert ample["primitive"] is True
    assert ample["square"] == 44258
    assert ample["D136_degree"] == 3488
    assert ample["minimum_complete_curve_degree"] == 1
    assert ample["minimum_imported_facet_pairing"] == 30
    assert ample["minimum_facet_indices_zero_based"] == [5, 11]
    assert all(value > 0 for value in ample["degrees_on_A1_through_A19"])
    assert all(
        value > 0
        for value in ample["pairings_with_transported_facets_0_through_28"]
    )
    assert all(simultaneous["checks"].values())


def test_twelve_minus_four_vectors_remain_separate(certificate: dict) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_four_divisor_classes"]
    assert mirrors["square_minus_four_divisibility_two_count"] == len(records) == 12
    assert [row["facet_zero_based"] for row in records] == [
        0, 1, 6, 7, 8, 9, 12, 18, 20, 24, 27, 28
    ]
    assert all((row["square"], row["divisibility"]) == (-4, 2) for row in records)
    assert all(
        row["effectivity_status"].startswith("not asserted") for row in records
    )
    assert all(mirrors["square_minus_four_checks"].values())
    assert all(mirrors["checks"].values())


def test_report_freezes_the_positive_claim_and_boundary() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a16_h11_29_ud8a1x2_nef_effectivity.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "87-reflection" in normalized
    assert "17/17" in normalized
    assert "facet `23` is exactly `a16`" in normalized
    assert "sixteen are reducible" in normalized
    assert "square `44258`" in normalized
    assert "minimum imported-wall pairing is `30`" in normalized
    assert "not the full nef chamber" in normalized
    assert "no effectivity claim" in normalized
