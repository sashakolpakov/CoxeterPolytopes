from __future__ import annotations

import json

from scripts.certify_alexeev_k3_h9_14_factory import OUTPUT, build


def test_five_construction_three_diagram_factory_rebuilds() -> None:
    payload, artifacts = build()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert all(path.read_bytes() == content for path, content in artifacts.items())

    assert payload["summary"] == {
        "raw_marked_subgroup_row_count": 5,
        "ambient_weighted_diagram_section_cone_orbit_count": 3,
        "distinct_abstract_relation_diagram_count": 3,
        "finite_volume_noncompact_count": 5,
        "compact_count": 0,
        "total_bad_edge_count": 0,
        "total_bad_parabolic_neighborhood_count": 0,
    }
    assert payload["diagram_types"] == {
        "T1": ["A14_8_1_C01"],
        "T2": ["A18_4_0_C01", "A18_4_0_C14"],
        "T3": ["A18_4_0_C07", "A18_4_0_C35"],
    }
    assert payload["scope"]["complete_nef_section_carriers_certified"] is False
    assert payload["scope"]["worldwide_novelty_claimed"] is False
    assert payload["decision"]["worldwide_novelty_status"] == (
        "mixed_T1_T2_open_T3_known"
    )
    assert payload["lineage"]["T3"]["abstract_diagram_status"] == (
        "known_Bredon_P_star_and_Scharlau_Kirschmer_II_9_1_3_power_4"
    )
    assert payload["lineage"]["worldwide_new_diagram_count_certified"] == 0
    assert payload["ambient_section_cone_orbit_audit"][
        "ambient_section_cone_orbit_count"
    ] == 3
    assert {
        key: row["conjugator_count"]
        for key, row in payload["ambient_section_cone_orbit_audit"][
            "duplicate_pair_audits"
        ].items()
    } == {
        "A18_4_0_C01__A18_4_0_C14": 192,
        "A18_4_0_C07__A18_4_0_C35": 72,
    }
    assert all(
        row["exact_cone"]["decision"]["classification"]
        == "finite_volume_noncompact"
        and row["literal_vinberg"]["literal_vinberg_passes"] is True
        and row["time_orientation"]["all_extreme_rays_same_time_cone"] is True
        for row in payload["constructions"]
    )
    assert {
        row["diagram_type"]: row["source_wall_audit"][
            "active_facet_source_norm_histogram"
        ]
        for row in payload["constructions"]
    } == {
        "T1": {"-4": 8, "-2": 6},
        "T2": {"-4": 0, "-2": 14},
        "T3": {"-4": 0, "-2": 14},
    }
    assert {
        row["diagram_type"]: row["source_wall_audit"][
            "folded_primitive_facet_norm_histogram"
        ]
        for row in payload["constructions"]
    } == {
        "T1": {"-8": 2, "-4": 10, "-2": 2},
        "T2": {"-4": 10, "-2": 4},
        "T3": {"-6": 5, "-2": 9},
    }
    assert {
        row["id"]: (
            row["time_orientation"]["minimum_strict_interior_pairing"],
            row["time_orientation"]["maximum_strict_interior_pairing"],
        )
        for row in payload["constructions"]
    } == {
        "A14_8_1_C01": ("21", "226"),
        "A18_4_0_C01": ("22", "136"),
        "A18_4_0_C07": ("26", "276"),
        "A18_4_0_C14": ("22", "129"),
        "A18_4_0_C35": ("26", "269"),
    }
