from __future__ import annotations

import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.inventory_shimada_disc16_weighted_wall_symmetries import (
    SOURCE,
    _weighted_wall_graph_automorphisms_fast,
)
from scripts.search_apery_fermi_chamber_folds import (
    _weighted_wall_graph_automorphisms,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_weighted_wall_symmetry_inventory.json"
)


def _archive() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_full_weighted_denominator_and_new_high_rank_queue() -> None:
    payload = _archive()
    assert payload["schema"] == (
        "coxeter.shimada_disc16_weighted_wall_symmetry_inventory.v1"
    )
    assert payload["scope"] == {
        "all_graph_symmetries_lifted_to_rational_lorentzian_isometries": True,
        "all_published_chambers": True,
        "all_subgroups_of_nontrivial_graph_groups": True,
        "exact_weighted_wall_graph_automorphisms": True,
        "minimum_fixed_rank": 11,
        "volume_gate_performed": False,
        "wall_angle_gate_performed": False,
    }
    assert payload["summary"] == {
        "archived_and_weighted_groups_equal_count": 3964,
        "chamber_count": 4539,
        "nontrivial_weighted_group_chamber_count": 673,
        "proper_high_dimension_histogram": {
            "10": 4,
            "11": 101,
            "13": 1,
            "14": 2,
            "18": 574,
        },
        "proper_high_dimensional_subgroup_class_count": 682,
        "total_subgroup_conjugacy_class_count": 1407,
        "total_subgroup_count_in_nontrivial_weighted_groups": 1496,
        "weighted_group_order_histogram": {
            "1": 3866,
            "2": 659,
            "3": 4,
            "4": 2,
            "6": 5,
            "8": 1,
            "20": 1,
            "48": 1,
        },
        "weighted_group_strict_extension_count": 575,
    }

    strict = [
        row
        for row in payload["nontrivial_weighted_group_chambers"]
        if not row["groups_equal"]
    ]
    assert len(strict) == 575
    assert {row["archived_stabilizer_order"] for row in strict} == {1}
    assert {row["weighted_wall_group_order"] for row in strict} == {2}
    assert {
        group_element["primal_action_integral"]
        for row in strict
        for group_element in row["weighted_group_elements"]
        if group_element["order"] > 1
    } == {False}

    extra = [
        class_row
        for row in strict
        for class_row in row["classes"]
        if class_row["order"] > 1 and class_row["above_h9"]
    ]
    assert len(extra) == 575
    assert sum(row["hyperbolic_dimension"] == 18 for row in extra) == 574
    assert sum(row["hyperbolic_dimension"] == 13 for row in extra) == 1


@pytest.mark.skipif(not SOURCE.is_file(), reason="pinned disc16.txt is unavailable")
def test_refined_enumerator_matches_generic_exact_graph_matcher() -> None:
    atlas = ShimadaDisc11Atlas(SOURCE)
    dual_gram = sp.Matrix(atlas.gram).inv()
    scale = sp.ilcm(*[sp.Rational(entry).q for entry in dual_gram])
    integral_scaled_dual_gram = scale * dual_gram

    # D57 has a small archived order-six group; D451 is one of the newly
    # detected rational C2 extensions of an archived-trivial chamber.
    for chamber_index, expected_order in ((57, 6), (451, 2)):
        chamber = atlas.chamber(chamber_index)
        walls = tuple(wall for orbit in chamber.wall_orbits for wall in orbit)
        generic = _weighted_wall_graph_automorphisms(dual_gram, walls)
        refined = _weighted_wall_graph_automorphisms_fast(
            integral_scaled_dual_gram, walls
        )
        assert len(generic) == len(refined) == expected_order
        assert set(generic) == set(refined)
