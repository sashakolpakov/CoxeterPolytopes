import json

from scripts.certify_roulleau_rank10_two_full_nef_simplices import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_two_rank_ten_full_nef_simplices_are_frozen() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert [
        (
            row["id"],
            row["lattice"],
            row["matrix_certificate"]["determinant"],
            row["literal_vinberg"]["finite_vertex_count"],
            row["literal_vinberg"]["ideal_vertex_count"],
            row["literal_vinberg"]["elliptic_edge_count"],
            row["matrix_certificate"]["reference"]["name"],
        )
        for row in payload["cases"]
    ] == [
        ("R10-UE8-NEF-H9-10", "U+E8(-1)", -1, 9, 1, 45, "JKRT T9"),
        ("R10-UD8-NEF-H9-10", "U+D8(-1)", -4, 7, 3, 45, "JKRT Q9"),
    ]
    assert all(
        row["literal_vinberg"]["decision"] == "finite_volume_noncompact"
        for row in payload["cases"]
    )


def test_exact_vertices_facets_and_ample_orientation() -> None:
    for row in _archive()["cases"]:
        matrix = row["matrix_certificate"]
        assert matrix["inertia"] == [1, 9, 0]
        assert matrix["ample_wall_pairings"] == [1] * 10
        assert all(value >= 0 for value in matrix["opposite_vertex_norms"])
        assert all(witness["square"] > 0 for witness in matrix["facet_witnesses"])
        assert all(
            witness["pairings"].count(0) == 1
            and witness["pairings"][witness["wall"] - 1] == 0
            for witness in matrix["facet_witnesses"]
        )
        assert matrix["checks"]["literal_jkrt_graph_isomorphism"]


def test_full_certificate_replays() -> None:
    assert build_certificate() == _archive()
