import json
from pathlib import Path

from scripts.certify_h9_14_s16_strong_transfer_corollary import (
    OUTPUT,
    build_certificate,
)


def test_synthesis_archive_is_exact() -> None:
    result = build_certificate()
    assert json.loads(Path(OUTPUT).read_text()) == result
    assert result["all_checks_pass"] is True
    assert all(result["checks"].values())


def test_oriented_s16_closes_all_four_obligations() -> None:
    result = build_certificate()
    strong = result["strong_marked_S16_class"]
    assert strong["name"] == "D16_plus"
    assert strong["finite_A_D"]["diagram_orbits"] == 1
    assert strong["finite_A_D"]["stable_even_fixed_overlattice_packages"] == 16
    assert strong["finite_A_D"]["strong_nonempty_packages"] == 1
    assert strong["marked_equivalence_proved"] is True
    assert set(strong["obligations"]) == {"U", "R", "C", "V"}
    assert all(strong["checks"].values())
    assert strong["unmarked_discrete_Gram_image_count"] == 1


def test_bare_scaled_dual_warning_is_explicit() -> None:
    result = build_certificate()
    warning = result["bare_S16_corollary"]
    assert "unmarked" in warning["statement"]
    assert "formal marked U is not asserted" in warning["marked_warning"]
    assert "complementary" in warning["marked_warning"]


def test_s17_is_an_exact_carrier_countercalibration() -> None:
    result = build_certificate()
    counter = result["S17_countercalibration"]
    assert counter["finite_A_D"] == {
        "diagram_orbits": 1,
        "stable_even_fixed_overlattice_packages": 5,
        "weak_K3_packages": 1,
        "strong_nonempty_packages": 0,
    }
    assert all(counter["checks"].values())
    assert counter["minimal_carrier_equivalence_false"] is True
    assert "other K3 lift alphabets" in counter["scope"]


def test_milestone_assessment_stays_scoped() -> None:
    result = build_certificate()
    status = result["milestone_assessment"]
    assert status["M3b"].startswith("closed for the sharply scoped")
    assert "broad milestone should remain open" in status["M3c"]
    assert status["broad_surface_softening"].startswith("open")
