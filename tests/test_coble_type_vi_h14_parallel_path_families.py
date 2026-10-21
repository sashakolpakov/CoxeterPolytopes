from __future__ import annotations

import json

import pytest

from scripts.certify_coble_type_vi_h14_parallel_path_families import (
    OUTPUT,
    _facet_formula,
    certify,
)


@pytest.fixture(scope="module")
def live() -> dict:
    # N=4 keeps the redundant wall-union replay tiny.  The collision proof
    # inside certify() is still exhaustive over every integer exponent.
    return certify(direct_max_n=4)


def test_archived_certificate_is_green() -> None:
    payload = json.loads(OUTPUT.read_text())
    assert payload["schema"] == "coxeter.coble_type_vi_h14_parallel_path_families.v1"
    assert payload["all_checks_pass"] is True
    assert payload["method_scope"] == {
        "all_N_proof_is_independent_of_direct_range": True,
        "cdd_used": False,
        "coxiter_used": False,
        "direct_reconstruction_range": "1<=N<=12",
        "extreme_rays_loaded_or_enumerated": False,
        "memory_profile": (
            "15x15 exact matrices, O(M^2) scalar collision equations for "
            "M<=101, and bounded O(NM) wall sets only"
        ),
    }


@pytest.mark.parametrize(
    ("parent", "expected"),
    [
        ("H14/100", [100, 180, 261, 341]),
        ("H14/101", [101, 182, 264, 345]),
    ],
)
def test_live_all_integer_collision_proof_and_formula(
    live: dict, parent: str, expected: list[int]
) -> None:
    family = live["families"][parent]
    collision = family["all_integer_collision_audit"]
    assert collision["bounded_exponent_scan_used"] is False
    assert collision["cross_seed_wall_collision_count"] == 0
    assert collision["negative_root_collision_count"] == 0
    assert collision["endpoint_to_side_collision_count"] == 0
    assert collision["identities_are_exactly_the_14_common_orthogonals"] is True
    assert [
        row["facet_count"] for row in family["direct_exact_wall_union_checks"]
    ] == expected
    assert [_facet_formula(parent, n) for n in range(1, 5)] == expected


def test_exact_second_generation_anchors(live: dict) -> None:
    anchors = {
        parent: family["N_equals_2_frozen_lineage_anchor"]
        for parent, family in live["families"].items()
    }
    assert anchors["H14/100"]["frozen_second_generation_id"] == (
        "COBLE-VI3-D4A02-H14-180"
    )
    assert anchors["H14/101"]["frozen_second_generation_id"] == (
        "COBLE-VI3-D4B04-H14-182"
    )
    assert {row["facet_count"] for row in anchors.values()} == {180, 182}
    assert all(row["exact_coordinate_list_match"] for row in anchors.values())


def test_full_divisor_lattice_and_support_for_every_n(live: dict) -> None:
    for family in live["families"].values():
        carrier = family["full_lattice_for_every_N"]
        assert carrier["arithmetic_seed_walls_excluding_b_smith_diagonal"] == [1] * 15
        assert carrier["physical_seed_walls_excluding_b_smith_diagonal"] == [1] * 15
        assert family["facet_support"][
            "all_parent_facets_have_strict_timelike_witnesses"
        ] is True
        assert family["all_facets_irredundant_for_every_N"] is True
        assert family["finite_volume"] is True
        assert family["compact"] is False


def test_two_rays_are_pairwise_separated(live: dict) -> None:
    theorem = live["theorem"]
    assert theorem["indices_in_H14_60"] == "both have index 2*N"
    assert theorem["two_polytopes_at_each_N_are_nonisometric"] == (
        "their facet counts differ by exactly N"
    )
    assert all(
        _facet_formula("H14/101", n) - _facet_formula("H14/100", n) == n
        for n in range(1, 100)
    )
