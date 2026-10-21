from __future__ import annotations

from collections import Counter
from fractions import Fraction

from coxeter_polytopes.search.i9_u2_vinberg import (
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


def test_i9_u2_vinberg_stops_at_the_scharlau_kirschmer_chamber() -> None:
    search = run_search(16)
    assert len(height_zero_simple_roots()) == 10
    assert len(search.roots) == 13
    assert search.stopping.fired
    assert search.stopping.finite_vertex_count == 34
    assert search.stopping.cusp_count == 4
    assert search.stopping.hyperideal_ray_count == 0
    assert search.stopping.opposite_sheet_ray_count == 0
    assert Counter((pair(root, root), divisibility(root)) for root in search.roots) == Counter(
        {(-2, 1): 9, (-4, 2): 2, (-1, 1): 2}
    )


def test_i9_u2_full_chamber_solomon_sum() -> None:
    search = run_search(16)
    positive_gram = [
        [str(-pair(left, right)) for right in search.roots]
        for left in search.roots
    ]
    relations = _relation_matrix(positive_gram)
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
        265,
        615,
        1001,
        1168,
        975,
        565,
        208,
        34,
    ]
    assert Fraction(euler["rational_euler_characteristic"]) == Fraction(
        -527, 33_443_020_800
    )
    assert Fraction(-527, 4_976_640) / Fraction(
        euler["rational_euler_characteristic"]
    ) == 6_720
