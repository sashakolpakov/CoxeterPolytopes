"""Exact tests for silver arithmetic tree-host cone triage."""

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver, silver_determinant
from coxeter_polytopes.search.silver_arithmetic_tree_hosts import (
    ARITHMETIC_EDGE_PRODUCT,
    GENTLE_DOTTED_PRODUCT,
    adjacent_two_dotted_rank23_host,
    analyze_tree_host,
    attach_arithmetic_leaf,
    genuine_dotted_hyperideal_lower_bound,
    rank23_tree_host_census,
)


def test_arithmetic_edge_and_integral_symmetric_grams() -> None:
    assert ARITHMETIC_EDGE_PRODUCT == 2 + SQRT2
    assert ARITHMETIC_EDGE_PRODUCT.sign(0) > 0
    assert (ARITHMETIC_EDGE_PRODUCT - 4).sign(0) < 0
    assert ARITHMETIC_EDGE_PRODUCT.conjugate().sign(0) > 0
    for kind in ("A", "D", "B"):
        for attachment in range(22):
            gram = attach_arithmetic_leaf(kind, 22, attachment)
            assert all(entry.is_integral for row in gram for entry in row)
            assert all(gram[i][j] == gram[j][i] for i in range(23) for j in range(23))
            assert silver_determinant(gram)


def test_complete_rank23_attachment_census_and_best_cones() -> None:
    rows = rank23_tree_host_census()
    assert len(rows) == 66
    arithmetic = {kind: [row for row in rows if row.kind == kind and row.arithmetic]
                  for kind in ("A", "D", "B")}
    assert {kind: [row.attachment for row in values]
            for kind, values in arithmetic.items()} == {
        "A": [0, 1, 2, 3, 18, 19, 20, 21],
        "D": [19, 20, 21],
        "B": [19, 20, 21],
    }
    assert {kind: max(row.finite_count for row in values)
            for kind, values in arithmetic.items()} == {"A": 3, "D": 3, "B": 3}
    assert all(row.ideal_count == 0 for values in arithmetic.values() for row in values)
    for kind in ("A", "D", "B"):
        best = analyze_tree_host(kind, 22, 0 if kind == "A" else 21)
        assert best.selected_inertia == Inertia(22, 1, 0)
        assert best.conjugate_inertia == Inertia(23, 0, 0)
        assert (best.finite_count, best.ideal_count, best.hyperideal_count) == (3, 0, 20)


def test_terminal_d22_and_b22_inverse_diagonal_is_exact() -> None:
    d = analyze_tree_host("D", 22, 21)
    assert d.inverse_diagonal[:4] == (
        Silver(19, -1),
        Silver(21, -1) / 4,
        Silver(21, -1) / 4,
        Silver(18, -1),
    )
    assert d.inverse_diagonal[-3:] == (
        Silver(1, -1),
        -SQRT2,
        Silver(1, -1) / 2,
    )
    b = analyze_tree_host("B", 22, 21)
    assert b.inverse_diagonal[:2] == (
        Silver(21, -1) / 2,
        Silver(10) - SQRT2 / 2,
    )
    assert b.inverse_diagonal[-3:] == (
        Silver(1, -1) / 2,
        -SQRT2 / 2,
        Silver(1, -1) / 4,
    )


def test_two_genuine_dotted_edges_force_almost_every_ray_hyperideal() -> None:
    assert genuine_dotted_hyperideal_lower_bound(23, ((0, 1), (0, 2))) == 22
    assert genuine_dotted_hyperideal_lower_bound(23, ((0, 1), (2, 3))) == 23
    assert genuine_dotted_hyperideal_lower_bound(30, ((4, 5), (4, 6), (4, 7))) == 29

    gram = adjacent_two_dotted_rank23_host()
    assert GENTLE_DOTTED_PRODUCT == Silver(3, 2)
    assert (GENTLE_DOTTED_PRODUCT - 4).sign(0) > 0
    assert GENTLE_DOTTED_PRODUCT.conjugate().sign(0) > 0
    assert (GENTLE_DOTTED_PRODUCT.conjugate() - 4).sign(0) < 0
    from coxeter_polytopes.search.silver_arithmetic_tree_hosts import inverse_diagonal, ray_kind
    from coxeter_polytopes.exact.silver import silver_inertia
    diagonal = inverse_diagonal(gram)
    assert silver_inertia(gram, embedding=0) == Inertia(22, 1, 0)
    assert silver_inertia(gram, embedding=1) == Inertia(23, 0, 0)
    assert tuple(map(ray_kind, diagonal)).count("finite") == 1
    assert tuple(map(ray_kind, diagonal)).count("ideal") == 0
    assert tuple(map(ray_kind, diagonal)).count("hyperideal") == 22
    # The Schur complement at the common dotted endpoint is exact.
    assert diagonal[0] == 1 / (Silver(22) / 21 - GENTLE_DOTTED_PRODUCT)
