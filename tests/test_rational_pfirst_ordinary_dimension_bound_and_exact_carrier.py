"""Hostile tests for the ordinary P-first n<=9 exact-carrier theorem."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.certify_rational_pfirst_ordinary_dimension_bound_and_exact_carrier import (
    DEPENDENCIES,
    ROOT,
    payload,
)


REPORT = (
    ROOT
    / "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md"
)
AUDIT = ROOT / (
    "results/"
    "rational_pfirst_ordinary_dimension_bound_and_exact_carrier_hostile_audit.md"
)
TIMELIKE = (
    ROOT / "results/rational_pfirst_timelike_anticanonical_vertex_closure.md"
)
ARCHIVE = ROOT / (
    "data/enumerations/"
    "rational_pfirst_ordinary_dimension_bound_and_exact_carrier.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_exact_branch_partition_and_dimension_bound() -> None:
    data = payload()
    assert data["sharp_dimension_bound"] == {
        "hyperbolic_dimension_at_most": 9,
        "rank_at_most": 10,
        "n_at_least_10": "IMPOSSIBLE_BY_FINITE_VOLUME",
    }

    rows = data["records"]
    assert [row["hyperbolic_dimension"] for row in rows[:7]] == list(
        range(3, 10)
    )
    assert [row["branch"] for row in rows[:6]] == [
        "timelike_exact_carrier"
    ] * 6
    assert rows[6]["branch"] == "characteristic_cusp_exact_carrier"
    assert all(
        row["branch"] == "spacelike_wall_cone_infinite_volume"
        for row in rows[7:]
    )
    assert all(row["exact_rational_carrier_proved"] for row in rows[:7])
    assert all(
        not row["finite_volume_possible_in_ordinary_lane"]
        for row in rows[7:]
    )


def test_anticanonical_color_and_escape_polynomial_are_exact() -> None:
    data = payload()
    for row in data["records"]:
        n = row["hyperbolic_dimension"]
        assert row["canonical_square"] == 9 - n
        assert [
            species["anticanonical_wall_slack"]
            for species in row["ordinary_facet_species"]
        ] == [1, 0]
        assert row["escape_polynomial"]["quadratic_coefficient"] == 9 - n

    # For q(t)=a+2bt+ct^2 with a>0 and c<0, the root product a/c is
    # negative.  Hence the two real roots have opposite signs and the
    # positive one is the strict light-cone crossing used in the proof.
    for a in (1, 2, 11):
        for b in (-7, 0, 5):
            for c in (-1, -2, -9):
                discriminant = 4 * (b * b - a * c)
                root_product_numerator = a
                root_product_denominator = c
                assert discriminant > 0
                assert root_product_numerator * root_product_denominator < 0


def test_positive_branch_dependencies_are_pinned_and_replayed() -> None:
    data = payload()
    for item in DEPENDENCIES.values():
        assert sha256(ROOT / item["path"]) == item["sha256"]
    assert data["positive_branch_composition"] == {
        "dependency_hashes_verified": True,
        "timelike_exceptional_counts": [3, 4, 5, 6, 7, 8],
        "timelike_all_exact_full_nef_returns": True,
        "rank_ten_characteristic_cusp_exact_full_nef_return": True,
    }


def test_orientation_repair_is_explicit_and_non_circular() -> None:
    combined = REPORT.read_text(encoding="utf-8")
    audit = AUDIT.read_text(encoding="utf-8")
    timelike = TIMELIKE.read_text(encoding="utf-8")

    assert "This does **not** by itself choose the positive sheet" in timelike
    assert "an oriented reflection isometry \\(w\\) carries \\(P\\)" in timelike
    assert "Only now, because the standard \\(F\\) is the finite vertex" in (
        timelike
    )
    assert "makes it future-directed and its\nray an ideal vertex" in combined
    assert "canonical escape is a conclusion of the raw finite-volume" in (
        combined
    )
    assert "K\\notin P" in combined
    assert "The timelike orientation was repaired" in audit
    assert "No assertion that \\(w\\) is nef, future-directed" in audit
    assert "Twice-primitive square-\\(-4\\) facets are not ordinary" in combined


def test_frozen_archive_is_canonical_replay() -> None:
    expected = json.dumps(payload(), indent=2, sort_keys=True) + "\n"
    assert ARCHIVE.read_text(encoding="utf-8") == expected
