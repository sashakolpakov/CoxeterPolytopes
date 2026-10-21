from __future__ import annotations

from collections import Counter
from fractions import Fraction

from coxeter_polytopes.search.u_e8_a1_vinberg import (
    divisibility,
    height_zero_simple_roots,
    pair,
    run_search,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import _relation_matrix


def test_u_e8_a1_vinberg_stops_on_the_full_twelve_wall_chamber() -> None:
    search = run_search(40)
    assert len(height_zero_simple_roots()) == 10
    assert len(search.roots) == 12
    assert search.levels[-1].height == Fraction(1, 2)
    assert search.stopping.fired
    assert search.stopping.finite_vertex_count == 18
    assert search.stopping.cusp_count == 1
    assert search.stopping.hyperideal_ray_count == 0
    assert search.stopping.opposite_sheet_ray_count == 0
    assert Counter((pair(root, root), divisibility(root)) for root in search.roots) == Counter(
        {(-2, 1): 11, (-2, 2): 1}
    )


def test_u_e8_a1_full_chamber_has_the_expected_coxeter_diagram() -> None:
    search = run_search(40)
    positive_gram = [
        [str(-pair(left, right)) for right in search.roots]
        for left in search.roots
    ]
    relations = _relation_matrix(positive_gram)
    assert Counter(
        relations[left][right]
        for left in range(len(relations))
        for right in range(left + 1, len(relations))
    ) == Counter({"2": 55, "3": 10, "0": 1})

