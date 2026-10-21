import json
from pathlib import Path

from scripts.certify_roulleau_ue8a1x3_all_subgroup_folds import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(Path(OUTPUT).read_text())


def test_archive_recomputes_exactly() -> None:
    assert build_certificate() == _archive()


def test_s3_denominator_and_vinberg_counts() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert [row["order"] for row in payload["subgroup_census"]["classes"]] == [
        1,
        2,
        3,
        6,
    ]
    types = {row["type_id"]: row for row in payload["types"]}
    c2 = types["R13_UE8A1x3_C2"]
    c3 = types["R13_UE8A1x3_C3_S3"]
    assert (
        c2["geometry"]["hyperbolic_dimension"],
        c2["geometry"]["facet_count"],
        c2["geometry"]["finite_vertex_count"],
        c2["geometry"]["ideal_vertex_count"],
    ) == (11, 15, 54, 3)
    assert (
        c3["geometry"]["hyperbolic_dimension"],
        c3["geometry"]["facet_count"],
        c3["geometry"]["finite_vertex_count"],
        c3["geometry"]["ideal_vertex_count"],
    ) == (10, 13, 31, 2)


def test_full_reflection_lineage_is_intrinsic_index_one() -> None:
    for record in _archive()["types"]:
        lineage = record["maximal_reflection_lineage"]
        assert lineage["full_reflection_group"]["reflection_subgroup_index"] == 1
        assert all(lineage["intrinsic_all_root_vinberg"]["checks"].values())
        assert lineage["all_integral_overlattice_census"][
            "unique_maximal_integral_overlattice"
        ]
    types = {row["type_id"]: row for row in _archive()["types"]}
    assert types["R13_UE8A1x3_C2"]["maximal_reflection_lineage"][
        "all_integral_overlattice_census"
    ]["integral_overlattice_indices"] == [1, 2]
    assert types["R13_UE8A1x3_C3_S3"]["maximal_reflection_lineage"][
        "all_integral_overlattice_census"
    ]["integral_overlattice_indices"] == [1]


def test_priority_and_surface_action_are_not_conflated() -> None:
    for record in _archive()["types"]:
        status = record["status_separation"]
        assert status == {
            "source_coordinate_validity": "proved_exactly",
            "geometric_surface_action": "not_realized_on_the_general_K3",
            "abstract_hyperbolic_chamber": "proved_finite_volume_noncompact",
            "full_integral_reflection_group": "proved_exact_index_one",
            "worldwide_abstract_priority": "known_not_new",
            "algebro_geometric_realization_priority": "open_literature_check",
        }

