from __future__ import annotations

import json

from scripts.certify_shimada_disc16_d4_duplicate_survivors import OUTPUT


def _load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_all_three_survivors_are_exact_duplicate_chambers() -> None:
    certificate = _load()
    records = {record["source"]: record for record in certificate["identifications"]}
    assert set(records) == {"D0004/C001", "D0004/C003", "D0004/C004"}
    assert records["D0004/C001"]["target"] == "D5-H14-60"
    assert records["D0004/C001"]["hyperbolic_dimension"] == 14
    for source in ("D0004/C003", "D0004/C004"):
        assert records[source]["target"] == "D5-H10-30"
        assert records[source]["hyperbolic_dimension"] == 10
    for record in records.values():
        comparison = record["comparison"]
        assert sorted(comparison["source_to_target_facet_permutation_one_based"]) == list(
            range(1, comparison["facet_count"] + 1)
        )
        assert comparison["complete_edge_coloured_graph_isomorphic"] is True
        assert comparison["all_normalized_gram_entries_equal_after_permutation"] is True
        assert comparison["same_hyperbolic_chamber_up_to_isometry"] is True
        assert record["finite_volume_noncompact_inherited_from_exact_target_isometry"] is True


def test_target_volume_facts_and_denominator_decision() -> None:
    certificate = _load()
    assert certificate["target_volume_facts"] == {
        "D5-H10-30": {
            "cusps": 207,
            "finite_vertices": 1333,
            "hyperideal_rays": 0,
            "vinberg_edges": 10030,
        },
        "D5-H14-60": {
            "cusps": 1027,
            "finite_vertices": 27845,
            "hyperideal_rays": 0,
        },
    }
    assert all(certificate["checks"].values())
    assert certificate["decision"] == {
        "D0004_C001_is_D5_H14_60": True,
        "D0004_C003_is_D5_H10_30": True,
        "D0004_C004_is_D5_H10_30": True,
        "D0004_C007_is_the_only_new_literal_type_among_the_four_D0004_survivors": True,
        "three_cast_instances_add_zero_abstract_chamber_types": True,
    }
