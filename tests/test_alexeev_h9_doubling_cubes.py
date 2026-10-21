from __future__ import annotations

import json

from scripts.certify_alexeev_h9_doubling_cubes import OUTPUT, build


def test_alexeev_h9_doubling_cube_archive_rebuilds() -> None:
    payload, artifacts = build()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert all(path.read_bytes() == content for path, content in artifacts.items())
    assert payload["theorem"] == {
        "seed_count": 2,
        "nonempty_doubling_subset_count_per_seed": 15,
        "new_project_ledger_isomorphism_type_count": 9,
        "facet_count_range": [15, 18],
        "copy_count_range": [2, 16],
        "all_types_finite_volume_noncompact": True,
        "all_types_pairwise_nonisomorphic": True,
    }
    assert len(payload["certified_types"]) == 9
    assert payload["bounded_lineage_audit"][
        "H9_facet_15_through_18_header_histogram"
    ] == {"15": 8}
    assert payload["bounded_lineage_audit"][
        "total_exact_labelled_graph_match_count"
    ] == 0
    assert all(
        row["exact_cone"]["decision"]["classification"]
        == "finite_volume_noncompact"
        and row["strict_timelike_facet_witness_count"] == row["facet_count"]
        and row["time_orientation"]["all_extreme_rays_same_time_cone"] is True
        and row["literal_vinberg"]["literal_vinberg_passes"] is True
        and row["literal_vinberg"]["bad_edge_count"] == 0
        and row["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0
        for row in payload["certified_types"]
    )
