from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c06_ue8d4_nef_effectivity import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_signed_fiber_bridge_corrects_the_positive_cone_component(
    certificate: dict,
) -> None:
    marking = certificate["signed_marking"]
    stock = certificate["published_nef_chamber_stock"]
    assert marking["old_Roulleau_interior_D506_degree"] == -15889
    assert marking["signed_initial_D506_degree"] == 15889
    assert stock["complete_irreducible_minus_two_curve_count"] == 15
    assert stock["ample_class"]["square"] == 506
    assert stock["ample_class"]["degrees_on_A1_through_A15"] == [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 25, 1, 1
    ]
    assert all(marking["checks"].values())
    assert all(stock["checks"].values())


def test_explicit_235_reflection_word_reaches_the_closed_nef_chamber(
    certificate: dict,
) -> None:
    transport = certificate["explicit_weyl_transport"]
    assert transport["word_length"] == len(transport["firing_trace"]) == 235
    assert transport["W_determinant"] == -1
    assert transport["combined_determinant"] == 1
    assert transport["final_nef_carrier_coordinates"] == [
        38, 77, 116, 58, 98, 80, 63, 46, 34, 22, 11, 0, "-9/2", "-11/2"
    ]
    assert transport["final_nef_carrier_square"] == 385
    assert transport["final_pairings_with_A1_through_A15"] == [
        1, 0, 1, 0, 0, 1, 0, 5, 0, 1, 0, 1, 9, 11, 0
    ]
    assert transport["primitive_integral_H_equals_2v_coordinates"] == [
        76, 154, 232, 116, 196, 160, 126, 92, 68, 44, 22, 0, -9, -11
    ]
    assert transport["H_square"] == 1540
    assert all(transport["checks"].values())


def test_all_minus_two_mirrors_have_exact_effectivity_and_irreducibility(
    certificate: dict,
) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_two_effective_decompositions"]
    irreducible = [
        (row["facet_zero_based"], row["published_curve_label"])
        for row in records
        if row["is_itself_an_irreducible_curve"]
    ]
    reducible = [row for row in records if not row["is_itself_an_irreducible_curve"]]
    assert mirrors["square_minus_two_effective_count"] == len(records) == 26
    assert mirrors["square_minus_two_irreducible_count"] == 5
    assert mirrors["square_minus_two_reducible_count"] == len(reducible) == 21
    assert irreducible == [
        (0, "A1"), (3, "A12"), (8, "A3"), (10, "A10"), (14, "A6")
    ]
    assert min(row["total_component_multiplicity"] for row in reducible) == 3
    assert max(row["total_component_multiplicity"] for row in reducible) == 55
    assert max(row["maximum_component_coefficient"] for row in records) == 7
    assert all(all(row["decomposition_checks"].values()) for row in records)
    assert all(mirrors["effective_checks"].values())


def test_D506_plus_7H_is_primitive_ample_and_minimal_in_the_family(
    certificate: dict,
) -> None:
    interior = certificate["explicit_simultaneous_ample_interior"]
    B6 = interior["B6_equals_D506_plus_6H"]
    B7 = interior["B7_equals_D506_plus_7H"]
    assert interior["minimal_strict_k_in_chosen_family"] == 7
    assert interior["nonpositive_imported_facets_for_k_0_through_6"] == [
        [[5, -24]], [[5, -20]], [[5, -16]], [[5, -12]],
        [[5, -8]], [[5, -4]], [[5, 0]],
    ]
    assert B6["unique_nonpositive_imported_facet"] == {
        "facet_zero_based": 5,
        "pairing": 0,
    }
    assert B7["coordinates"] == [
        578, 1171, 1765, 882, 1492, 1220, 963, 707, 522, 338, 169, 1, -75, -77
    ]
    assert B7["primitive"] is True
    assert B7["square"] == 90666
    assert B7["D506_degree"] == 7856
    assert B7["minimum_complete_curve_degree"] == 1
    assert B7["minimum_imported_facet_pairing"] == 4
    assert B7["minimum_facet_zero_based"] == 5
    assert all(interior["checks"].values())


def test_five_minus_four_vectors_are_kept_separate(certificate: dict) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_four_divisor_classes"]
    assert mirrors["square_minus_four_divisibility_two_count"] == len(records) == 5
    assert [row["facet_zero_based"] for row in records] == [4, 5, 6, 12, 13]
    assert all((row["square"], row["divisibility"]) == (-4, 2) for row in records)
    assert all(
        row["effectivity_status"].startswith("not asserted") for row in records
    )
    assert all(mirrors["square_minus_four_checks"].values())
    assert all(mirrors["checks"].values())


def test_standalone_report_freezes_the_claim_boundary() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c06_ue8d4_nef_effectivity.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "235-reflection" in normalized
    assert "exactly five" in normalized
    assert "21 are reducible" in normalized
    assert "b_7^2=90666" in normalized
    assert "minimal in the chosen family" in normalized
    assert "not the full nef cone" in normalized
    assert "no effectivity claim" in normalized
