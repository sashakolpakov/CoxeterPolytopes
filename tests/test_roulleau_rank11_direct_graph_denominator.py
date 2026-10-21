import json

from scripts.certify_roulleau_rank11_direct_graph_denominator import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_archive_recomputes_exactly() -> None:
    assert build_certificate() == _archive()


def test_weighted_groups_and_subgroup_denominator() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert [row["weighted_wall_group"]["order"] for row in payload["cases"]] == [1, 2]
    assert [
        row["subgroup_denominator"]["subgroup_count"] for row in payload["cases"]
    ] == [1, 2]
    assert [
        record["fixed_rank"]
        for row in payload["cases"]
        for record in row["subgroup_denominator"]["records"]
        if record["order"] > 1
    ] == [10]
    assert payload["summary"]["proper_fixed_rank_at_least_11_count"] == 0


def test_published_source_chambers_pass_literal_vinberg_gate() -> None:
    observed = []
    for row in _archive()["cases"]:
        geometry = row["published_full_nef_chamber"]["geometry"]
        observed.append(
            (
                geometry["facet_count"],
                geometry["finite_vertex_count"],
                geometry["ideal_vertex_count"],
                geometry["vinberg_volume_criterion"]["elliptic_edge_count"],
            )
        )
        assert geometry["finite_volume_certified"]
        assert geometry["compact"] is False
        assert geometry["all_checks_pass"]
    assert observed == [(12, 18, 1, 99), (14, 35, 4, 218)]


def test_lineage_separates_known_source_from_new_fold_search() -> None:
    cases = {row["id"]: row for row in _archive()["cases"]}
    assert cases["R11_UE8A1"]["published_full_nef_chamber"]["lineage"][
        "reflection_subgroup_index"
    ] == 1
    d8 = cases["R11_UD8A1"]["published_full_nef_chamber"]["lineage"]
    assert [row["parity"] for row in d8["maximal_stable_lattices"]] == [
        "even",
        "even",
        "odd",
    ]
    assert [
        row["containing_full_reflection_index"]
        for row in d8["maximal_stable_lattices"]
    ] == [1023, 1023, None]
    assert _archive()["summary"]["new_abstract_polytope_count"] == 0

