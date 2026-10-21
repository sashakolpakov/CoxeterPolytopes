from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_alexeev_parked_h9_doubling_products import (
    OUTPUT,
    build,
)


@lru_cache(maxsize=1)
def _rebuilt():
    return build()


def test_product_family_archive_and_binary_artifacts_rebuild_exactly() -> None:
    payload, artifacts = _rebuilt()
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload == archived
    for path, content in artifacts.items():
        assert path.read_bytes() == content


def test_complete_cube_and_infinite_family_formulas() -> None:
    payload, _artifacts = _rebuilt()
    cube = payload["finite_C2_cubed_census"]
    assert cube["complete_nonempty_subset_count"] == 7
    assert cube["labelled_isomorphism_type_count"] == 3
    assert [row["facet_count"] for row in cube["types"]] == [14, 15, 16]

    families = payload["infinite_families"]
    assert families["A14_8_1_C02"]["facet_formula"] == "4*N+13"
    assert families["A14_8_1_C02"]["index_and_volume_multiple"] == "2*N"
    assert families["A14_10_1_C02"]["facet_formula"] == "3*N+14"
    assert families["A14_10_1_C02"]["index_and_volume_multiple"] == "4*N"
    rectangle = families["A16_14_1_C01"]
    assert rectangle["facet_formula"] == "m*n+m+n+12"
    assert rectangle["index_and_volume_multiple"] == "m*n"
    assert rectangle["group_factor"] == "D_infinity x D_infinity"


def test_six_direct_volume_and_vinberg_replays_are_green() -> None:
    payload, _artifacts = _rebuilt()
    theorem = payload["theorem"]
    assert theorem == {
        "all_direct_bad_edge_count": 0,
        "all_direct_bad_parabolic_neighborhood_count": 0,
        "direct_GMP_Vinberg_representative_count": 6,
        "finite_new_project_diagram_count": 3,
        "infinite_family_count": 3,
    }
    expected = {
        "A14_10_1_C01_D2_W3": (14, 94, 8),
        "A14_10_1_C01_D4_W3_5": (15, 121, 10),
        "A14_10_1_C01_D8_W3_5_11": (16, 148, 14),
        "A14_8_1_C02_P2_X_C2": (21, 334, 60),
        "A14_10_1_C02_P2_X_C2SQ": (20, 312, 49),
        "A16_14_1_C01_RECT_2X2": (20, 210, 24),
    }
    for row in payload["direct_exact_audits"]:
        exact = row["exact_cone"]["exact_double_description"]
        assert (
            row["facet_count"],
            exact["finite_vertex_count"],
            exact["ideal_vertex_count"],
        ) == expected[row["id"]]
        assert exact["hyperideal_vertex_count"] == 0
        assert row["time_orientation"]["all_extreme_rays_same_time_cone"]
        assert row["literal_vinberg"]["bad_edge_count"] == 0
        assert row["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0
    assert not any(payload["local_prior_H9_diagram_matches"].values())
