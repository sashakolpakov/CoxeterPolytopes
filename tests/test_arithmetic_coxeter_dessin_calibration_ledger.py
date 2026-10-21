from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_arithmetic_coxeter_dessin_calibration_ledger import (
    DEFAULT_OUTPUT,
    PINNED_SHA256,
    ROOT,
    RSS_CEILING_BYTES,
    build,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    return build()


def test_archive_is_exact_builder_output(certificate: dict) -> None:
    assert json.loads(DEFAULT_OUTPUT.read_text()) == certificate
    assert certificate["all_checks_pass"]
    assert RSS_CEILING_BYTES == 1 << 29


def test_every_source_is_sha_pinned(certificate: dict) -> None:
    used = {
        pin["path"]
        for group in (certificate["rows"], certificate["determinant_line_controls"])
        for section in group
        for pin in section["source_pins"]
    }
    assert used <= set(PINNED_SHA256)
    assert len(used) >= 20
    assert all(len(PINNED_SHA256[path]) == 64 for path in used)


def test_exact_required_polytope_matrix(certificate: dict) -> None:
    got = {
        row["id"]: (
            row["polytope"]["dimension"],
            row["polytope"]["facets"],
            row["polytope"]["vertices"]["finite"],
            row["polytope"]["vertices"]["ideal_vertices_or_cusps"],
            row["polytope"]["vertices"]["hyperideal"],
        )
        for row in certificate["rows"]
    }
    assert got == {
        "BUGAENKO-H8-11": (8, 11, 41, 0, 0),
        "TORIC-H6-9": (6, 9, 18, 3, 0),
        "COBLE-II-H10-30": (10, 30, 1333, 207, 0),
        "COBLE-I-H11-26": (11, 26, 508, 42, 0),
        "COBLE-VI3-H14-60": (14, 60, 27845, 1027, 0),
        "COBLE-MII3-H17-48": (17, 48, 190808, 5244, 0),
    }


def test_field_lane_is_not_surface_characteristic(certificate: dict) -> None:
    rows = {row["id"]: row for row in certificate["rows"]}
    assert rows["BUGAENKO-H8-11"]["polytope"]["arithmetic_definition_field"] == "Q(sqrt(5))"
    for key in rows.keys() - {"BUGAENKO-H8-11"}:
        assert rows[key]["polytope"]["arithmetic_definition_field"] == "Q"
    assert "characteristic 3" in rows["COBLE-VI3-H14-60"]["surface_carrier"]["ground_field_or_characteristic"]
    assert "characteristic 3" in rows["COBLE-MII3-H17-48"]["surface_carrier"]["ground_field_or_characteristic"]


def test_directional_synthesis_is_honest(certificate: dict) -> None:
    synthesis = certificate["synthesis"]["general_compiler_support"]
    assert synthesis["direct_determinant_and_heegner_calibration"] == ["BUGAENKO-H8-11"]
    assert synthesis["direct_cusp_calibration"] == ["TORIC-H6-9"]
    assert set(synthesis["native_exact_P_to_surface_experiments"]) == {
        "TORIC-H6-9",
        "COBLE-II-H10-30",
        "COBLE-I-H11-26",
        "COBLE-VI3-H14-60",
        "COBLE-MII3-H17-48",
    }
    assert set(synthesis["surface_to_polytope_only_controls"]) == {
        "E10-H9-10-POSITIVE",
        "UE8D4-H13-15-NEGATIVE",
    }
    assert "L4" in certificate["rows"][2]["surface_carrier"]["exact_return_level"]


def test_automorphic_levels_are_not_conflated(certificate: dict) -> None:
    controls = {row["id"]: row for row in certificate["determinant_line_controls"]}
    e10 = controls["E10-H9-10-POSITIVE"]["automorphic_status"]
    assert e10["pure_borcherds_product"] == "proved"
    assert "weight 165" in e10["exact_global_reflective_product"]
    negative = controls["UE8D4-H13-15-NEGATIVE"]["automorphic_status"]
    assert negative["arrangement"].startswith("available")
    assert negative["exact_global_reflective_product"].startswith("impossible")
    bug = controls["BUGAENKO-H8-11-OBSTRUCTED"]["automorphic_status"]
    assert bug["pure_borcherds_product"].startswith("obstructed")
    assert bug["first_complete_compensator_shell"] == "obstructed"
    assert "not explicit" in bug["compensated_slice_exact_product"]


def test_cusp_encoding_and_terminology(certificate: dict) -> None:
    text = json.dumps(certificate, sort_keys=True).lower()
    report = (ROOT / "results/arithmetic_coxeter_dessin_calibration_ledger.md").read_text().lower()
    producer = Path(__file__).resolve().parents[1].joinpath(
        "scripts/certify_arithmetic_coxeter_dessin_calibration_ledger.py"
    ).read_text().lower()
    forbidden_boundary_phrase = "ideal " + "faces"
    assert forbidden_boundary_phrase not in text
    assert forbidden_boundary_phrase not in report
    assert forbidden_boundary_phrase not in producer
    for row in certificate["rows"]:
        if row["polytope"]["vertices"]["ideal_vertices_or_cusps"]:
            encoding = row["cusp_transfer"]["surface_boundary_encoding"]
            assert "cluster" in encoding
            assert "dual-complex" in encoding
    assert "no single-point encoding" in certificate["rows"][1]["cusp_transfer"]["surface_boundary_encoding"]


def test_single_missing_invariant_is_explicit(certificate: dict) -> None:
    invariant = certificate["synthesis"]["single_missing_classification_invariant"]
    assert invariant["name"] == "minimal faithful Coxeter-dessin Hilbert degree T_dess(P)"
    assert "exceptional boundary clusters" in invariant["definition"]
    assert "Hilbert/Hurwitz stacks" in invariant["why_it_matters"]
    assert invariant["current_status"].startswith("open")


def test_uncomputed_row_specific_automorphic_data_stay_open(certificate: dict) -> None:
    for row in certificate["rows"][1:]:
        assert row["determinant_line_and_heegner"]["status"] == "open_not_row_specifically_compiled"
        status = row["automorphic_status"]
        assert status["status"] == "open"
        assert status["pure_borcherds_product"] is None
        assert status["compensated_slice_exact_product"] is None
