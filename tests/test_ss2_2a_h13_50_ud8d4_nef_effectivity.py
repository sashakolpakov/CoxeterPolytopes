from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ss2_2a_h13_50_ud8d4_nef_effectivity import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_signed_marking_resolves_the_positive_cone_component(certificate: dict) -> None:
    marking = certificate["signed_marking"]
    assert marking["chosen_marking"] == "-X"
    assert marking["old_X_interior_D154_degree"] == -10008
    assert marking["signed_initial_D154_degree"] == 10008
    assert all(marking["checks"].values())


def test_explicit_60_reflection_word_reaches_the_closed_nef_chamber(
    certificate: dict,
) -> None:
    transport = certificate["explicit_weyl_transport"]
    assert transport["word_length"] == len(transport["firing_trace"]) == 60
    assert transport["W_determinant"] == 1
    assert transport["final_nef_carrier_square"] == 2012
    assert transport["final_nef_carrier_D154_degree"] == 1192
    assert transport["final_nef_carrier_coordinates"] == [
        26, 52, 78, 120, 81, 162, 123, 84, 61, 38, -13, 3, 19, 6
    ]
    assert transport["final_pairings_with_A1_through_A20"] == [
        0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 64, 32, 0, 26, 38, 0, 32, 64, 38, 26
    ]
    assert all(transport["checks"].values())


def test_all_40_minus_two_classes_have_exact_effective_decompositions(
    certificate: dict,
) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_two_effective_decompositions"]
    assert mirrors["square_minus_two_effective_count"] == len(records) == 40
    assert mirrors["square_minus_two_irreducible_count"] == 0
    assert mirrors["square_minus_two_reducible_count"] == 40
    assert min(row["total_component_multiplicity"] for row in records) == 2
    assert max(row["total_component_multiplicity"] for row in records) == 28
    assert max(row["maximum_component_coefficient"] for row in records) == 6
    assert all(row["square"] == -2 and row["divisibility"] == 1 for row in records)
    assert all(row["is_itself_an_irreducible_curve"] is False for row in records)
    assert all(all(row["decomposition_checks"].values()) for row in records)
    assert all(mirrors["effective_checks"].values())


def test_2H_plus_D154_is_simultaneously_ample_and_coxeter_strict(
    certificate: dict,
) -> None:
    interior = certificate["explicit_simultaneous_ample_interior"]
    first = interior["B1_equals_H_plus_D154"]
    second = interior["B2_equals_2H_plus_D154"]
    assert first["square"] == 4550
    assert first["unique_nonpositive_imported_facet"] == {
        "facet_zero_based": 22,
        "pairing": -2,
    }
    assert [
        (index, value)
        for index, value in enumerate(
            first["pairings_with_transported_facets_0_through_49"]
        )
        if value <= 0
    ] == [(22, -2)]
    assert second["coordinates"] == [
        53, 107, 162, 250, 169, 339, 260, 182, 137, 93, -18, 14, 38, 16
    ]
    assert second["primitive"] is True
    assert second["square"] == 12970
    assert second["D154_degree"] == 2538
    assert second["degrees_on_A1_through_A20"] == [
        1, 1, 33, 1, 1, 1, 1, 33, 1, 1, 129, 65, 17, 61, 93, 17, 81, 129, 85, 53
    ]
    assert second["pairings_with_transported_facets_0_through_49"] == [
        36, 36, 34, 74, 38, 35, 78, 38, 36, 38,
        52, 39, 80, 35, 44, 54, 72, 36, 56, 70,
        58, 28, 4, 82, 41, 38, 143, 161, 168, 181,
        182, 144, 163, 165, 149, 169, 162, 167, 171, 150,
        166, 187, 151, 183, 165, 183, 151, 168, 188, 167,
    ]
    assert second["minimum_complete_curve_degree"] == 1
    assert second["minimum_imported_facet_pairing"] == 4
    assert second["minimum_facet_zero_based"] == 22
    assert all(interior["checks"].values())


def test_ten_minus_four_vectors_are_separate_without_effectivity_claim(
    certificate: dict,
) -> None:
    mirrors = certificate["transported_mirrors"]
    records = mirrors["square_minus_four_divisor_classes"]
    assert mirrors["square_minus_four_divisibility_two_count"] == len(records) == 10
    assert [row["facet_zero_based"] for row in records] == [
        3, 6, 8, 12, 14, 16, 19, 21, 22, 23
    ]
    assert all((row["square"], row["divisibility"]) == (-4, 2) for row in records)
    assert all(
        row["effectivity_status"].startswith("not asserted") for row in records
    )
    assert all(mirrors["square_minus_four_checks"].values())
    assert all(mirrors["checks"].values())


def test_standalone_report_freezes_the_exact_claims() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "ss2_2a_h13_50_ud8d4_nef_effectivity.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "60-letter weyl word" in normalized
    assert "40/40" in normalized
    assert "0/40" in normalized
    assert "all forty are reducible" in normalized
    assert "no effectivity claim" in normalized
    assert "closed nef chamber" in normalized
    assert "simultaneous ample and coxeter-strict interior" in normalized
    assert "b_2^2 = 12970" in normalized
    assert "minimum is `4`, at facet `22`" in normalized
