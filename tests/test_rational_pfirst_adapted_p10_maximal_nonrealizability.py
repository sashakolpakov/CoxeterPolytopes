from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_pfirst_adapted_p10_maximal_nonrealizability import (
    ARCHIVE,
    build,
    dot,
    wall_packet,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_pfirst_adapted_p10_maximal_nonrealizability.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_frozen_archive() -> None:
    assert build() == frozen()


def test_physical_colours_and_both_canonical_escapes() -> None:
    payload = frozen()
    assert payload["walls"]["physical_square_split"] == {
        "minus_1": 1,
        "minus_2": 10,
        "minus_4": 1,
    }
    assert payload["walls"]["canonical_pairings_physical"] == [
        -1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
    ]
    log = payload["canonical_escape_and_log"]
    assert log["K_fails_E10_by"] == -1
    assert log["F_fails_boundary_half_by"] == -1
    assert log["D_equals_F_minus_B_over_2"] == [0] * 11
    assert log["boundary_incidence"]["weak_H1_H2_H3_gate"].startswith("PASS")


def test_exact_cone_has_no_hyperideal_rays_and_all_facets_survive() -> None:
    cone = frozen()["cone"]
    assert (
        cone["finite_ray_count"],
        cone["ideal_ray_count"],
        cone["hyperideal_ray_count"],
    ) == (17, 2, 0)
    assert cone["extreme_ray_count"] == 19
    assert cone["all_facets_irredundant"] is True
    assert cone["strict_interior_square"] == 3940
    assert all(z > 0 for z in cone["strict_interior_slacks"])
    assert all(row["square"] >= 0 and row["vector"][0] > 0 for row in cone["rays"])


def test_coxeter_support_is_the_twelve_wall_p10_support() -> None:
    labels, roots = wall_packet()
    assert labels == [
        "E10",
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
        "D6",
        "D7",
        "D8",
        "D9",
        "C",
        "F_boundary_half",
    ]
    assert [dot(v, v) for v in roots] == [-1] + [-2] * 10 + [-1]
    edges = {tuple(row) for row in frozen()["walls"]["Coxeter_edges_order_0_means_parallel"]}
    assert ("E10", "F_boundary_half", 0) in edges
    assert ("E10", "D9", 4) in edges
    assert ("D3", "C", 3) in edges
    assert len(edges) == 11


def test_classification_free_saturation_regression() -> None:
    saturation = frozen()["classification_free_reflection_saturation"]
    assert saturation["primitive_integral_reflector_norms"] == [-1, -2]
    assert "every integral hyperbolic reflection" in saturation["conclusion"]
    rows = saturation["bounded_regression"]["rows"]
    assert len(rows) == saturation["bounded_regression"]["reduced_root_count"]
    assert any(row["simple_F"] for row in rows)
    for row in rows:
        if not row["simple_F"]:
            assert abs(row["reflected_a"]) < row["a"]


def test_complex_no_go_and_h11_nonisomorphism_are_scoped() -> None:
    lineage = frozen()["lineage_and_AG_no_go"]
    assert lineage["identified_raw_polytope"].startswith("classical P10")
    assert lineage["complex_physical_full_nef_carrier"] is False
    assert "twelve effective irreducible roots" in lineage["frozen_AG_reason"]
    assert lineage["H11_26_comparison"]["same_polytope"] is False
    scope = frozen()["scope"]
    assert scope["base_field_for_nonrealizability"] == "C"
    assert "positive-characteristic nonrealizability" in scope["not_claimed"]
    assert "nonexistence of rational carriers with a different geometric colouring" in scope[
        "not_claimed"
    ]


def test_report_states_uncoloured_flex_carrier_and_missing_effectivity_gate() -> None:
    text = REPORT.read_text()
    assert "uncoloured chamber is classical `P10`" in text
    assert "irreducible genus-one `(-1)` curve" in text
    assert "not a reduced smooth rational curve" in text
    assert "P-native simultaneous reduced-irreducible effectivity" in text
    assert "No Nikulin--Vinberg" in text
