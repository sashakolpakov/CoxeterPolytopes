from __future__ import annotations

import json

from scripts.certify_roulleau_rank10_remaining_full_nef import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_certificate_replays_exactly() -> None:
    assert build_certificate() == _archive()


def test_all_seven_are_literal_finite_volume_full_nef_chambers() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert [row["facet_count"] for row in payload["cases"]] == [
        11,
        11,
        14,
        18,
        27,
        145,
        11,
    ]
    assert [
        (
            row["literal_vinberg"]["finite_vertex_count"],
            row["literal_vinberg"]["ideal_vertex_count"],
            row["literal_vinberg"]["elliptic_edge_count"],
        )
        for row in payload["cases"]
    ] == [
        (15, 2, 80),
        (16, 10, 125),
        (37, 7, 231),
        (186, 71, 1840),
        (313, 45, 2245),
        (13824, 1249, 123520),
        (17, 5, 105),
    ]
    for row in payload["cases"]:
        assert row["complete_full_nef_reason"]
        assert row["literal_vinberg"]["decision"] == "finite_volume_noncompact"
        assert row["literal_vinberg"]["bad_edge_count"] == 0
        assert row["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0


def test_lineage_and_source_erratum_are_bounded_precisely() -> None:
    payload = _archive()
    exact = [
        row["id"]
        for row in payload["cases"]
        if row["lineage"].get("literal_relation_graph_match")
    ]
    assert exact == [
        "R10-UE7A1-NEF-H9-11",
        "R10-UD4D4-NEF-H9-11",
        "R10-UA2E6-NEF-H9-11",
    ]
    errata = [row for row in payload["cases"] if row["source_erratum"]]
    assert [row["id"] for row in errata] == ["R10-UD4A1x4-NEF-H9-27"]
    assert errata[0]["matrix_certificate"]["ample_wall_pairings"][5] == 2
    assert errata[0]["matrix_certificate"]["ample_wall_pairings"][4] == 1
