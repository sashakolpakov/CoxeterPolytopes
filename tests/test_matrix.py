import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia, rank, rectangular_pairing


def test_determinant_and_rectangular_rank() -> None:
    assert determinant([]) == Golden(1)
    assert determinant([[2, -1], [-1, 2]]) == Golden(3)
    assert determinant([[1, 2], [2, 4]]) == Golden(0)
    assert rank([[1, 2, 3], [2, 4, 6], [0, 1, 0]]) == 2


def test_inertia_positive_affine_and_lorentzian() -> None:
    assert inertia([[2, -1], [-1, 2]]) == Inertia(2, 0, 0)
    affine_a2 = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
    assert inertia(affine_a2) == Inertia(2, 0, 1)
    assert inertia([[2, -3], [-3, 2]]) == Inertia(1, 1, 0)


def test_zero_diagonal_two_by_two_pivot() -> None:
    assert inertia([[0, 1], [1, 0]]) == Inertia(1, 1, 0)
    assert inertia([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) == Inertia(1, 1, 1)


def test_rectangular_pairing() -> None:
    gram = [[2, 0], [0, -1]]
    left = [[1], [1]]
    right = [[1, 0], [0, 1]]
    assert rectangular_pairing(left, gram, right) == [[2, -1]]


def test_golden_matrix_and_conjugate_embedding() -> None:
    matrix = [[2, -PHI], [-PHI, 2]]
    assert determinant(matrix) == 3 - PHI
    assert inertia(matrix) == Inertia(2, 0, 0)
    assert inertia([[PHI]], embedding=1) == Inertia(0, 1, 0)


def test_matrix_validation() -> None:
    with pytest.raises(ValueError):
        rank([[1], [1, 2]])
    with pytest.raises(ValueError):
        determinant([[1, 2]])
    with pytest.raises(ValueError):
        inertia([[1, 2], [3, 4]])
