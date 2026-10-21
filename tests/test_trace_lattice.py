import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    SQRT5,
    integral_gram,
    is_codifferent,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import e8_cartan_gram


def test_codifferent_and_trace_pairing() -> None:
    assert SQRT5 * SQRT5 == 5
    assert is_codifferent(INV_SQRT5)
    assert is_codifferent(-PHI * INV_SQRT5)
    assert trace_gram([[1]], INV_SQRT5) == [[0, 1], [1, 1]]


def test_integrality_and_evenness_checks() -> None:
    assert integral_gram([[0, 1], [1, 2]]) == [[0, 1], [1, 2]]
    assert is_even_integral([[0, 1], [1, 2]])
    assert not is_even_integral([[1]])
    assert not is_even_integral([[Golden(1).a / 2]])
    with pytest.raises(ValueError):
        integral_gram([[Golden(1).a / 2]])


def test_e8_standard_gram() -> None:
    e8 = e8_cartan_gram()
    assert determinant(e8) == 1
    assert inertia(e8) == Inertia(8, 0, 0)

