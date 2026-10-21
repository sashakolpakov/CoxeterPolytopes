from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_pfirst_anticanonical_genus_one_facet_dimension_bound import (
    ARCHIVE,
    build,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_archive() -> None:
    assert build() == frozen()


def test_intrinsic_bound_has_exact_hypotheses_and_no_classification_input() -> None:
    theorem = frozen()["theorem"]
    assert theorem["bound"] == "n<=11"
    assert theorem["sharp"] is True
    assert theorem["classification_inputs"] == []
    assert "F^2=K_S^2 is -1 or -2" in theorem["deduction"]
    assert "rho(S)=11 or 12" in theorem["deduction"]


def test_both_norms_are_realized_and_p11_is_sharp() -> None:
    rows = {row["n"]: row for row in frozen()["exact_returns"]}
    assert set(rows) == {10, 11}
    assert rows[10]["F_square"] == -1
    assert rows[11]["F_square"] == -2
    assert rows[10]["rho"] == 11
    assert rows[11]["rho"] == 12
    assert rows[10]["arithmetic_genus"] == rows[11]["arithmetic_genus"] == 1
    assert rows[10]["reflection_integral"] and rows[11]["reflection_integral"]
    assert rows[10]["complete_nef_cone_certified"]
    assert rows[11]["complete_nef_cone_certified"]
    assert rows[11]["known_chamber"].startswith("Vinberg P11")


def test_endpoint_volume_data_are_exact() -> None:
    rows = {row["n"]: row for row in frozen()["exact_returns"]}
    assert (rows[10]["finite_vertices"], rows[10]["ideal_vertices"], rows[10]["hyperideal_rays"]) == (
        17,
        2,
        0,
    )
    assert (rows[11]["finite_vertices"], rows[11]["ideal_vertices"], rows[11]["hyperideal_rays"]) == (
        26,
        2,
        0,
    )


def test_scope_does_not_conflate_anticanonical_and_boundary_branches() -> None:
    scope = frozen()["scope"]
    assert scope["facet_is_literal_curve_class_not_only_same_hyperplane"] is True
    assert "a bound for adapted (-4) boundary facets B=2b" in scope["not_claimed"]
    text = REPORT.read_text()
    assert "`P11`, not `P10`" in text
    assert "attains the upper endpoint" in text
    assert "irreducible genus-one" in text
    assert "reduced rational boundary" in text
