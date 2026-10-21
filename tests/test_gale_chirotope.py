from itertools import combinations, product

import pytest

from coxeter_polytopes.combinatorics.gale_chirotope import (
    EXPECTED_CHIROTOPE_SIGN_COUNT,
    EXPECTED_GRASSMANN_PLUECKER_COUNT,
    RANK_FOUR_BASES,
    audit_uniform_rank_four_gale,
    solve_uniform_rank_four_gale,
)
from coxeter_polytopes.exact.matrix import determinant


FULL_MASK = (1 << 14) - 1


# A small integral, uniform Gale realization of the dual boundary of
# (Delta_1)^3 x (Delta_3)^2.  This fixture pins the cofactor-sign convention
# independently of the Boolean solver.
EXACT_MATCHING_PRODUCT_GALE = (
    (6, 0, -1, -1),
    (5, -1, 0, -1),
    (1, 4, -1, 0),
    (0, 5, 0, 1),
    (1, -1, 6, 1),
    (-1, -1, 5, 0),
    (-1, 1, 0, 6),
    (-1, 1, -1, 5),
    (1, 1, 1, 4),
    (0, -1, 1, 5),
    (-5, -6, -6, -6),
    (-5, -6, -4, -6),
    (-6, -5, -5, -4),
    (-6, -6, -5, -5),
)


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _matching_product_facets() -> tuple[int, ...]:
    simplex_blocks = ((6, 7, 8, 9), (10, 11, 12, 13))
    complements = []
    for binary_choices in product((0, 1), (2, 3), (4, 5)):
        for omitted_first in simplex_blocks[0]:
            for omitted_second in simplex_blocks[1]:
                facet = set(binary_choices)
                facet.update(set(simplex_blocks[0]) - {omitted_first})
                facet.update(set(simplex_blocks[1]) - {omitted_second})
                complements.append(FULL_MASK ^ _mask(facet))
    return tuple(sorted(complements))


def _pentagon_product_facets() -> tuple[int, ...]:
    pentagon_edges = {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}
    simplex_blocks = ((5, 6, 7, 8), (9, 10, 11, 12, 13))
    complements = []
    for pentagon_edge in pentagon_edges:
        for omitted_first in simplex_blocks[0]:
            for omitted_second in simplex_blocks[1]:
                facet = set(pentagon_edge)
                facet.update(set(simplex_blocks[0]) - {omitted_first})
                facet.update(set(simplex_blocks[1]) - {omitted_second})
                complements.append(FULL_MASK ^ _mask(facet))
    return tuple(sorted(complements))


def test_uniform_rank_four_counts() -> None:
    assert EXPECTED_CHIROTOPE_SIGN_COUNT == len(tuple(combinations(range(14), 4))) == 1001
    assert EXPECTED_GRASSMANN_PLUECKER_COUNT == 45_045


def test_exact_integral_gale_matrix_pins_circuit_sign_convention() -> None:
    minors = tuple(
        determinant([EXACT_MATCHING_PRODUCT_GALE[index] for index in basis])
        for basis in RANK_FOUR_BASES
    )
    assert all(minor != 0 for minor in minors)
    signs = tuple(1 if minor > 0 else -1 for minor in minors)
    audit = audit_uniform_rank_four_gale(_matching_product_facets(), signs)
    assert audit.passes is True
    assert audit.positive_circuit_masks == _matching_product_facets()


@pytest.mark.parametrize(
    "facets",
    (_matching_product_facets(), _pentagon_product_facets()),
)
def test_polytopal_products_have_exact_gale_certificates(facets) -> None:
    pytest.importorskip("z3")
    result = solve_uniform_rank_four_gale(facets, timeout_ms=30_000)
    assert result.status == "sat"
    audit = audit_uniform_rank_four_gale(facets, result.signs)
    assert audit.passes is True
    assert audit.positive_circuit_masks == facets
    assert audit.every_element_in_positive_circuit is True
    assert audit.every_element_avoided_by_positive_circuit is True

    # The sign certificate itself is independent of the claimed facet list.
    tampered = facets[:-1]
    assert audit_uniform_rank_four_gale(tampered, result.signs).passes is False
