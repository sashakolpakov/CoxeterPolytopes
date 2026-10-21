from __future__ import annotations

import ast
from functools import lru_cache
import json

from scripts.audit_jacobian_kummer_h12_four_vertex_subtree_hostile import (
    OUTPUT,
    ROOT,
    audit,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return audit()


def test_hostile_replay_matches_frozen_archive() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["result"] == {
        "decision": "PASS_FREEZE",
        "P4_type_count": 1964,
        "K1_3_type_count": 680,
        "pairwise_nonisomorphic_ordinary_Coxeter_type_count": 2644,
        "facet_histogram": {"607": 23, "616": 170, "625": 844, "634": 1607},
        "reflection_index_and_volume_multiplier": 4,
        "finite_volume": True,
        "compact": False,
        "worldwide_priority": "open",
    }


def test_independent_orbit_denominators_and_boundary_census() -> None:
    payload = _payload()
    seed = payload["independent_seed_and_action"]
    assert seed["folded_action_order"] == 384
    assert seed["stock_color_orbit_sizes"] == [16, 96]
    assert seed["fixed_core_Smith_diagonal"] == [1] * 13
    assert all(seed["checks"].values())

    denominator = payload["coloring_denominator"]
    assert denominator["P4"]["ordered_words"] == 1_379_952
    assert denominator["P4"]["unrooted_words"] == 696_192
    assert denominator["P4"]["orbit_count"] == 1_964
    assert denominator["K1_3"]["unrooted_words"] == 227_920
    assert denominator["K1_3"]["orbit_count"] == 680

    census = payload["facet_and_relation_census"]
    assert census["type_histogram"] == {
        "607": 23,
        "616": 170,
        "625": 844,
        "634": 1607,
    }
    assert census["shape_type_histograms"] == {
        "P4": {"607": 17, "616": 127, "625": 626, "634": 1194},
        "K1,3": {"607": 6, "616": 43, "625": 218, "634": 413},
    }
    assert census["raw_unrooted_histogram"] == {
        "607": 2480,
        "616": 47616,
        "625": 293376,
        "634": 580640,
    }


def test_weighted_two_walk_nonisomorphism_and_production_agreement() -> None:
    payload = _payload()
    ordinary = payload["ordinary_nonisomorphism"]
    assert ordinary["initial_degree_profile_count"] == 2535
    assert ordinary["collision_group_count"] == 102
    assert ordinary["collision_member_count"] == 211
    assert ordinary["maximum_collision_group_size"] == 4
    assert ordinary["independently_refined_member_count"] == 211
    assert ordinary["weighted_two_walk_maximum_entry"] == 30089
    assert ordinary["weighted_two_walk_term_bound_at_634_facets"] == 31066
    assert ordinary["weighted_two_walk_permutation_invariance_checks"] == 211
    assert ordinary["all_collision_groups_separated"] is True
    assert ordinary["all_2644_types_distinct"] is True

    comparison = payload["production_comparison"]
    assert comparison["exact_row_match_count"] == 2644
    assert len(payload["representatives"]) == 2644
    hashes = [
        row["independent_ordinary_type_sha256"]
        for row in payload["representatives"]
    ]
    assert len(set(hashes)) == 2644


def test_index_volume_carrier_scope_and_no_production_import() -> None:
    payload = _payload()
    carrier = payload["carrier_and_volume"]
    assert carrier["Allcock_finite_subtree_theorem_applies"] is True
    assert carrier["reflection_index"] == carrier["volume_multiplier"] == 4
    assert carrier["finite_volume_noncompact"] is True
    assert carrier["fixed_index_one_Picard_core"] is True
    assert carrier["same_seven_dimensional_full_Picard_K3_family"] is True
    assert carrier["weak_divisor_reflection_carrier_only"] is True
    assert carrier["full_nef_effectivity_or_surface_automorphism_claim"] is False

    script = (
        ROOT
        / "scripts/audit_jacobian_kummer_h12_four_vertex_subtree_hostile.py"
    )
    tree = ast.parse(script.read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert not any(name == "scripts" or name.startswith("scripts.") for name in imported)
    assert not any(
        name == "coxeter_polytopes" or name.startswith("coxeter_polytopes.")
        for name in imported
    )
