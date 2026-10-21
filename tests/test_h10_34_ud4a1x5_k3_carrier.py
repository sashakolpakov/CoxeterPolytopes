from __future__ import annotations

import json
from functools import lru_cache

from scripts.certify_h10_34_ud4a1x5_k3_carrier import (
    OUTPUT,
    _payload_sha256,
    build_certificate,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return json.loads(OUTPUT.read_text())


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == _payload()
    assert _payload_sha256(_payload()) == (
        "55b990d3002e787f7199ca6577e2e39a05ccf713ca05365e08fafd8212395e88"
    )


def test_both_stable_overlattice_denominators_are_exhaustive() -> None:
    census = _payload()["stable_closure_and_scaled_dual_census"]
    assert census["A15"]["stable_integral_subgroup_orders"] == [1, 4]
    assert census["A15"]["stable_even_subgroup_orders"] == [1, 4]
    assert census["SS2"]["stable_integral_subgroup_orders"] == [1, 2, 4, 8]
    assert census["SS2"]["stable_even_subgroup_orders"] == [1, 4]
    assert census["SS2"]["closures"][3]["scaled_duals"]["2"]["determinant"] == 128
    assert census["SS2"]["closures"][3]["scaled_duals"]["2"]["parity"] == "even"
    assert all(census["checks"].values())


def test_all_34_oriented_roots_share_one_marked_k3_lattice() -> None:
    bridge = _payload()["common_lattice_bridge"]
    assert bridge["discriminant_q_mod_2_histogram"] == {
        "0": 36,
        "1": 28,
        "1/2": 28,
        "3/2": 36,
    }
    assert sorted(bridge["SS2_to_A15_facet_permutation_zero_based"]) == list(
        range(34)
    )
    assert all(bridge["checks"].values())


def test_primitive_k3_embedding_and_ample_framing() -> None:
    payload = _payload()
    glue = payload["actual_K3_carrier"]["primitive_K3_lattice_embedding"]
    assert all(glue["checks"].values())
    ample = payload["ample_framing"]
    assert ample["Weyl_sequence_one_based"] == [20, 4, 14, 16, 3]
    assert ample["big_nef_sum_interior_square"] == 1648
    assert ample["big_nef_minimum_90_wall_pairing"] == 0
    assert ample["ample_square"] == 2178
    assert ample["ample_minimum_90_wall_pairing"] == 1
    assert ample["ample_minimum_imported_pairing"] == 4
    assert ample["generation_index"] == 1
    assert all(ample["checks"].values())


def test_every_target_reflection_is_frozen_and_integral() -> None:
    facets = _payload()["facet_records"]
    assert len(facets) == 34
    assert {(row["square"], row["divisibility"]) for row in facets} == {
        (-2, 1),
        (-2, 2),
        (-4, 2),
    }
    assert sum(row["square"] == -2 for row in facets) == 24
    assert sum(row["square"] == -4 for row in facets) == 10
    assert all(len(row["reflection_matrix"]) == 11 for row in facets)
    assert all(all(row["checks"].values()) for row in facets)


def test_literal_vinberg_counts_and_scope_are_explicit() -> None:
    payload = _payload()
    audit = payload["literal_Vinberg_audit"]
    assert (
        audit["finite_vertices"],
        audit["ideal_vertices"],
        audit["elliptic_rank9_edges"],
        audit["connected_affine_components"],
    ) == (1774, 146, 12792, 426)
    assert audit["bad_edges"] == audit["bad_parabolic_neighborhoods"] == 0
    assert all(audit["checks"].values())
    scope = payload["lineage_and_scope"]
    assert scope["not_a_new_abstract_polytope"] is True
    assert scope["not_the_complete_90_wall_nef_chamber"] is True
    assert scope["reflections_not_asserted_surface_automorphisms"] is True
    assert payload["all_checks_pass"] is True
