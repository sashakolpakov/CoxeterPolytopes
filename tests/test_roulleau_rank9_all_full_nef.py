from __future__ import annotations

import json

from scripts.certify_roulleau_rank9_all_full_nef import OUTPUT, build_certificate


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_certificate_replays_exactly() -> None:
    assert build_certificate() == _archive()


def test_all_ten_are_literal_finite_volume_full_nef_chambers() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert [row["facet_count"] for row in payload["cases"]] == [
        9,
        10,
        15,
        37,
        240,
        9,
        10,
        10,
        9,
        10,
    ]
    assert [
        (
            row["literal_vinberg"]["finite_vertex_count"],
            row["literal_vinberg"]["ideal_vertex_count"],
            row["literal_vinberg"]["elliptic_edge_count"],
        )
        for row in payload["cases"]
    ] == [
        (8, 1, 36),
        (12, 3, 63),
        (53, 11, 315),
        (569, 78, 4188),
        (17280, 2160, 207360),
        (7, 2, 36),
        (17, 4, 90),
        (18, 1, 81),
        (8, 1, 36),
        (14, 1, 63),
    ]
    for row in payload["cases"]:
        assert row["complete_full_nef_reason"]
        assert row["literal_vinberg"]["decision"] == "finite_volume_noncompact"
        assert row["literal_vinberg"]["bad_edge_count"] == 0
        assert row["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0


def test_lineage_large_del_pezzo_row_and_erratum_are_exact() -> None:
    payload = _archive()
    exact = [
        row["id"]
        for row in payload["cases"]
        if row["lineage"].get("literal_relation_graph_match")
    ]
    assert exact == [
        "R9-UE7-NEF-H8-9",
        "R9-UD6A1-NEF-H8-10",
        "R9-UA7-NEF-H8-9",
        "R9-UD4A3-NEF-H8-10",
        "R9-UD5A2-NEF-H8-10",
        "R9-UD7-NEF-H8-9",
        "R9-UA1E6-NEF-H8-10",
    ]
    del_pezzo = next(
        row for row in payload["cases"] if row["id"] == "R9-U2A1x7-NEF-H8-240"
    )
    assert del_pezzo["matrix_certificate"]["ample_square"] == 2
    assert set(del_pezzo["matrix_certificate"]["ample_wall_pairings"]) == {2}
    assert del_pezzo["literal_vinberg"]["edge_endpoint_type_histogram"] == {
        "0,2": 69120,
        "1,1": 138240,
    }
    erratum = next(row for row in payload["cases"] if row["source_erratum"])
    assert erratum["id"] == "R9-UD7-NEF-H8-9"
    assert erratum["matrix_certificate"]["ample_wall_pairings"] == [
        1,
        1,
        1,
        2,
        1,
        1,
        1,
        1,
        1,
    ]
