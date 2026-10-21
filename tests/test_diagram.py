import pytest

from coxeter_polytopes.exact.diagram import (
    connected_components,
    induced_subgram,
    is_affine_parabolic_component,
    is_elliptic,
    is_parabolic,
    parabolic_rank,
)
from coxeter_polytopes.exact.matrix import rank
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram


def cartan_from_edges(size: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    gram = [[2 if i == j else 0 for j in range(size)] for i in range(size)]
    for i, j in edges:
        gram[i][j] = gram[j][i] = -1
    return gram


def type_a(size: int) -> list[list[int]]:
    return cartan_from_edges(size, [(i, i + 1) for i in range(size - 1)])


def type_d(size: int) -> list[list[int]]:
    assert size >= 4
    return cartan_from_edges(
        size,
        [(i, i + 1) for i in range(size - 3)]
        + [(size - 3, size - 2), (size - 3, size - 1)],
    )


def test_induced_subgram_and_connected_components() -> None:
    gram = block_diagonal(type_a(3), type_a(2))
    assert connected_components(gram) == ((0, 1, 2), (3, 4))
    assert induced_subgram(gram, [4, 3]) == [[2, -1], [-1, 2]]
    assert connected_components([]) == ()


@pytest.mark.parametrize("gram", [type_a(1), type_a(8), type_d(4), type_d(8)])
def test_finite_a_and_d_are_elliptic(gram) -> None:
    assert is_elliptic(gram)
    assert not is_parabolic(gram)


def test_finite_e6_e7_e8_are_elliptic() -> None:
    e8 = e8_cartan_gram()
    e7 = induced_subgram(e8, [0, 1, 2, 3, 4, 5, 7])
    e6 = induced_subgram(e8, [0, 1, 2, 3, 4, 7])
    assert is_elliptic(e6)
    assert is_elliptic(e7)
    assert is_elliptic(e8)
    assert is_elliptic(block_diagonal(e6, type_a(2)))


def test_connected_affine_a2_and_d4() -> None:
    affine_a2 = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
    affine_d4 = cartan_from_edges(5, [(0, 1), (0, 2), (0, 3), (0, 4)])
    for gram, expected_rank in [(affine_a2, 2), (affine_d4, 4)]:
        assert is_affine_parabolic_component(gram)
        assert is_parabolic(gram)
        assert parabolic_rank(gram) == expected_rank
        assert rank(gram) == expected_rank
        assert not is_elliptic(gram)


def test_reducible_parabolic_rank_is_sum_of_component_ranks() -> None:
    affine_a2 = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
    affine_a1 = [[2, -2], [-2, 2]]
    gram = block_diagonal(affine_a2, affine_a1)
    assert connected_components(gram) == ((0, 1, 2), (3, 4))
    assert not is_affine_parabolic_component(gram)
    assert is_parabolic(gram)
    assert parabolic_rank(gram) == 3
    assert rank(gram) == 3


@pytest.mark.parametrize(
    "gram",
    [
        [[2, -3], [-3, 2]],  # indefinite
        [[2, 2], [2, 2]],  # corank one, but wrong Coxeter signs
        [[2, -2, 0], [-2, 2, 0], [0, 0, 2]],  # affine plus finite
        [[2, -1, -1, 0], [-1, 2, -1, 0], [-1, -1, 2, 0], [0, 0, 0, 0]],
    ],
)
def test_nonparabolic_matrices_are_rejected(gram) -> None:
    assert not is_parabolic(gram)
    with pytest.raises(ValueError, match="defined only for a parabolic"):
        parabolic_rank(gram)


def test_malformed_subdiagram_requests_are_rejected() -> None:
    with pytest.raises(ValueError, match="distinct"):
        induced_subgram(type_a(3), [0, 0])
    with pytest.raises(IndexError, match="out of range"):
        induced_subgram(type_a(3), [3])
    with pytest.raises(ValueError, match="symmetric"):
        is_elliptic([[2, -1], [0, 2]])
