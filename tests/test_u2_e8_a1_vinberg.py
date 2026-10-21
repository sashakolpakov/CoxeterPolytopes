from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

from coxeter_polytopes.search.u2_e8_a1_vinberg import (
    divisibility,
    height_zero_simple_roots,
    pair,
    run_search,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)


def test_u2_e8_a1_vinberg_stops_at_list10_entry11() -> None:
    search = run_search(16)
    assert len(height_zero_simple_roots()) == 10
    assert len(search.roots) == 13
    assert search.stopping.fired
    assert search.stopping.finite_vertex_count == 29
    assert search.stopping.cusp_count == 3
    assert search.stopping.hyperideal_ray_count == 0
    assert search.stopping.opposite_sheet_ray_count == 0
    assert Counter((pair(root, root), divisibility(root)) for root in search.roots) == Counter(
        {(-2, 1): 10, (-2, 2): 2, (-4, 2): 1}
    )


def test_u2_e8_a1_full_chamber_solomon_sum() -> None:
    search = run_search(16)
    positive_gram = [
        [str(-pair(left, right)) for right in search.roots]
        for left in search.roots
    ]
    relations = _relation_matrix(positive_gram)
    known_f2 = json.loads(
        Path("data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json").read_text(
            encoding="utf-8"
        )
    )["vinberg_maximal_integral_reflection_chamber"][
        "normalized_coxeter_relation_matrix"
    ]
    assert relations == known_f2
    assert Counter(
        relations[left][right]
        for left in range(len(relations))
        for right in range(left + 1, len(relations))
    ) == Counter({"2": 65, "3": 9, "4": 2, "0": 2})
    connected = _connected_masks_direct(relations)
    euler = _euler_from_connected_components(relations, connected)
    assert euler["finite_special_subset_count_by_rank"] == [
        1,
        13,
        76,
        264,
        606,
        966,
        1092,
        876,
        488,
        175,
        29,
    ]
    assert Fraction(euler["rational_euler_characteristic"]) == Fraction(
        -31, 11_147_673_600
    )
    assert Fraction(-527, 1_451_520) / Fraction(
        euler["rational_euler_characteristic"]
    ) == 130_560
