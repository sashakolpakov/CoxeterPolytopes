from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_pfirst_corrected_genus_one_pn_boundary import (
    ARCHIVE,
    build,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_pfirst_corrected_genus_one_pn_boundary.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_archive() -> None:
    assert build() == frozen()


def test_local_genus_one_packet_is_identical_across_all_ranks() -> None:
    rows = frozen()["rows"]
    assert [row["n"] for row in rows] == list(range(12, 20))
    assert {
        (row["C_square"], row["K_dot_C"], row["arithmetic_genus"], row["C_reflection_integral"])
        for row in rows
    } == {(-2, 2, 1, True)}
    assert frozen()["local_NS_theorem"]["conclusion"].endswith("once C != -K")


def test_corrected_anticanonical_formula_tracks_all_rank_dependence() -> None:
    for row in frozen()["rows"]:
        n = row["n"]
        assert row["C"] == [3] + [-1] * 11 + [0] * (n - 11)
        assert row["Z_equals_C_plus_K"] == [0] * 12 + [1] * (n - 11)
        assert row["C_dot_Z"] == 0
        assert row["Z_square"] == 11 - n


def test_scoped_Pn_carrier_boundary_is_exactly_sixteen() -> None:
    theorem = frozen()["scoped_Pn_AG_theorem"]
    assert theorem["complex_rational_full_nef_carrier_iff"] == "n<=16"
    assert theorem["positive_dimensions"] == [12, 13, 14, 15, 16]
    assert theorem["negative_dimensions"] == [17, 18, 19]
    assert "bare reflection maximality" in theorem["not_from"]
    rows = {row["n"]: row for row in frozen()["rows"]}
    assert all(rows[n]["complex_rational_full_nef_carrier"] for n in range(12, 17))
    assert all(not rows[n]["complex_rational_full_nef_carrier"] for n in range(17, 20))
    assert theorem["obstructing_quartic_classes"] == {
        "17": [4] + [-1] * 17,
        "18": [4] + [-1] * 18,
        "19": [4] + [-1] * 18 + [0],
    }


def test_negative_cases_have_marking_independent_MMP_Abel_certificates() -> None:
    rows = {row["n"]: row for row in frozen()["rows"]}
    assert rows[17]["marking_independent_MMP_Abel_obstruction"]["canonical_class_count"] == 8
    assert rows[18]["marking_independent_MMP_Abel_obstruction"]["canonical_class_count"] == 24
    assert rows[19]["marking_independent_MMP_Abel_obstruction"]["canonical_class_count"] == 120
    dependencies = frozen()["upstream_MMP_archives"]
    assert [item["dimensions"] for item in dependencies] == [[17], [18, 19]]
    assert dependencies[0]["certificate_payload_sha256"] == (
        "87586f25471616fe9e5d59fcd6e3d3ab6a1c38deab1b5e6cdb12cc19b98a873c"
    )
    assert dependencies[1]["certificate_payload_sha256"] == (
        "451072212c9e2def74f9f8f3d996e29e2f8816ec1a2dd21e1cc4b9ff3cce0bfa"
    )


def test_scope_keeps_local_and_family_theorems_separate() -> None:
    scope = frozen()["scope"]
    assert "a universal n<=16 theorem for all reflection-maximal chambers with a genus-one facet" in scope[
        "not_claimed"
    ]
    text = REPORT.read_text()
    assert "local Neron--Severi/canonical packet" in text
    assert "**rank-blind**" in text
    assert "P-native contraction/incidence compatibility" in text
    assert "arbitrary" in text
    assert "unique-genus-one maximal-chamber problem remains open" in text
    assert "integral plane cubic" in text
    assert "global smoothness of the putative cubic is not assumed" in " ".join(text.split())
    assert "a_n=\\min(n,18)" in text
    assert "literal quartic class changes" in text
