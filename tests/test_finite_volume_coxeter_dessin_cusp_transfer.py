from __future__ import annotations

import json

from scripts.certify_finite_volume_coxeter_dessin_cusp_transfer import (
    OUTPUT,
    build_certificate,
)


def test_archive_matches_exact_builder() -> None:
    archived = json.loads(OUTPUT.read_text())
    assert archived == build_certificate()
    assert archived["all_checks_pass"]


def test_all_three_cusps_have_literal_seven_plus_two_split() -> None:
    cusps = build_certificate()["calibration"]["cusps"]
    assert len(cusps) == 3
    assert all(cusp["vertical_wall_count"] == 7 for cusp in cusps)
    assert all(cusp["bounded_penetration_wall_count"] == 2 for cusp in cusps)


def test_vertical_and_bounded_records_are_distinguished() -> None:
    for cusp in build_certificate()["calibration"]["cusps"]:
        assert len(cusp["wall_records"]) == 9
        for wall in cusp["wall_records"]:
            if wall["vertical_through_cusp"]:
                assert wall["height_cap_squared"] is None
                assert wall["positive_K_norm_w"] == wall["rho_minus_root_square"]
            else:
                assert wall["height_cap_squared"] is not None


def test_finiteness_hypotheses_are_literal() -> None:
    certificate = build_certificate()
    hypotheses = " ".join(certificate["theorem"]["hypotheses"])
    assert "fixed-denominator" in hypotheses
    assert "parabolic rank n-1" in hypotheses
    assert "compact horosphere cell" in hypotheses
    assert "fixed rational Witt partner" in hypotheses
    assert certificate["field_split"]["scope"].startswith("admissible quadratic")
