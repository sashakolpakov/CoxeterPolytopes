from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a18_c01_c07_k3_nef_promotions import build_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "a18_c01_c07_k3_nef_promotions.json"
)
REPORT = ROOT / "results" / "a18_c01_c07_k3_nef_promotions.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_both_exact_k3_promotions_replay() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["case_count"] == 2
    assert [case["id"] for case in archived["cases"]] == [
        "A18_4_0_C01",
        "A18_4_0_C07",
    ]


def test_discriminant_glue_and_rootlessness_are_exact() -> None:
    cases = {case["id"]: case for case in _archived()["cases"]}
    c01 = cases["A18_4_0_C01"]
    assert c01["ambient_NS"]["discriminant_action_mod_two"] == [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    assert c01["ambient_NS"]["coinvariant_lattice"]["isometry_type"] == "E8(-2)"
    assert c01["transcendental_and_glue"]["glued_action_characteristic_polynomial"] == (
        "(lambda - 1)**14*(lambda + 1)**8"
    )

    c07 = cases["A18_4_0_C07"]
    assert c07["ambient_NS"]["discriminant_action_mod_two"] == [
        [0, 1, 1, 0],
        [0, 1, 1, 1],
        [1, 0, 0, 1],
        [1, 1, 0, 1],
    ]
    coin = c07["ambient_NS"]["coinvariant_lattice"]
    assert coin["isometry_type"] == "-(A2 tensor D4)"
    assert coin["minimum_positive_norm"] == 4
    assert coin["minus_two_rootless"]
    assert c07["transcendental_and_glue"]["glued_action_characteristic_polynomial"] == (
        "(lambda - 1)**10*(lambda**2 + lambda + 1)**6"
    )


def test_torelli_nef_and_normal_quotient_scope() -> None:
    for case in _archived()["cases"]:
        period = case["period_and_torelli"]
        assert period["picard_lattice_is_exactly_N"]
        assert period["NS_coinvariant_has_no_minus_two_vectors"]
        assert period["invariant_kahler_chamber_exists"]
        assert "genuine" in period["global_torelli_conclusion"]

        nef = case["nef_and_quotient_cone"]
        assert nef["archived_dimension_facets"] == [9, 14]
        assert nef["active_source_wall_norm_histogram"] == {"-2": 14, "-4": 0}
        assert nef["all_restricted_long_walls_redundant"]
        assert "Nef(X) intersect Fix(g)" in nef["invariant_nef_section_equality"]
        assert "normal finite quotient" in nef["quotient_nef_conclusion"]
        assert "resolution" in nef["quotient_nef_conclusion"]

    report = REPORT.read_text(encoding="utf-8")
    assert "symplectic Nikulin involution" in report
    assert "purely non-symplectic" in report
    assert "full invariant nef cone" in report
    assert "not the nef cone of a resolution" in report
