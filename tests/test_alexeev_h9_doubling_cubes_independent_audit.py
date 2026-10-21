from __future__ import annotations

import json

from scripts.audit_alexeev_h9_doubling_cubes_independent import OUTPUT, audit


def test_independent_alexeev_h9_doubling_cube_audit_rebuilds() -> None:
    payload = audit()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["independence"] == {
        "imports_primary_certifier": False,
        "calls_primary_build": False,
        "reconstructs_boundaries_from_factory_roots": True,
        "recomputes_relation_graphs": True,
    }
    assert payload["theorem"] == {
        "seed_count": 2,
        "total_nonempty_subset_count": 30,
        "all_subgroups_have_order_two_to_k": True,
        "all_boundaries_subgroup_invariant": True,
        "all_glued_unions_pass_exact_orthogonal_convexity_lemma": True,
        "all_facets_have_strict_timelike_witnesses": True,
        "all_subset_diagram_classes_equal_seed_automorphism_orbits": True,
        "primary_representative_count": 9,
        "primary_representatives_pairwise_nonisomorphic": True,
        "primary_representative_collision_count": 0,
    }
    assert all(
        seed["nonempty_subset_count"] == 15
        and all(
            row["copy_count"] == row["reflection_subgroup_order"]
            == row["sign_sector_count"]
            and row["facet_count"] == row["strict_facet_witness_count"]
            and row["boundary_subgroup_invariant"] is True
            and row["no_opposite_oriented_boundary_duplicates"] is True
            for row in seed["subsets"]
        )
        for seed in payload["seeds"].values()
    )
