import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.k3 import (
    K3_WALL_NORMS,
    field_inner_product,
    is_integral_reflection_root,
    is_k3_wall_norm,
    k3_trace_square,
    reflection_coefficients,
)


def test_five_k3_wall_norm_shells_are_exact() -> None:
    expected = (
        5 - 3 * PHI,
        4 - 2 * PHI,
        3 - PHI,
        Golden(2),
        1 + PHI,
    )
    assert K3_WALL_NORMS == expected
    assert all(k3_trace_square(value) == -2 for value in expected)
    assert all(value.sign(0) > 0 and value.sign(1) > 0 for value in expected)

    # The trace equation gives a=2-b.  Exact positivity leaves only -3<=b<=1.
    surviving = tuple(
        Golden(2 - b, b)
        for b in range(-100, 101)
        if is_k3_wall_norm(Golden(2 - b, b))
    )
    assert surviving == expected

    # In an even O_F-lattice both golden coefficients of a represented norm
    # are even.  The two surviving shells form one unit-square orbit.
    even_surviving = tuple(
        value
        for value in expected
        if value.a.denominator == value.b.denominator == 1
        and value.a.numerator % 2 == value.b.numerator % 2 == 0
    )
    assert even_surviving == (4 - 2 * PHI, Golden(2))
    assert (PHI - 1) ** 2 * 2 == 4 - 2 * PHI


def test_k3_wall_norm_rejections() -> None:
    assert not is_k3_wall_norm(4 - 3 * PHI)  # wrong trace square
    assert not is_k3_wall_norm(2 - 4 * PHI)  # negative at the identity place
    assert not is_k3_wall_norm(1 + 2 * PHI)  # negative at the conjugate place
    assert not is_k3_wall_norm(Golden(1, 1) / 2)  # not in Z[phi]


def test_field_pairing_and_reflection_integrality() -> None:
    gram = [[2, 1], [1, 2]]
    root = [1, 0]
    assert field_inner_product(gram, root, root) == 2
    assert reflection_coefficients(gram, root) == (Golden(2), Golden(1))
    assert is_integral_reflection_root(gram, root)

    nonintegral = [[4, 1], [1, 4]]
    assert reflection_coefficients(nonintegral, root) == (Golden(2), Golden(1) / 2)
    assert not is_integral_reflection_root(nonintegral, root)
    assert not is_integral_reflection_root([[2]], [Golden(1) / 2])


def test_reflection_input_validation() -> None:
    with pytest.raises(ValueError, match="symmetric"):
        reflection_coefficients([[1, 2], [0, 1]], [1, 0])
    with pytest.raises(ValueError, match="nonzero norm"):
        reflection_coefficients([[0]], [1])
    with pytest.raises(ValueError, match="vector lengths"):
        field_inner_product([[1]], [1, 0], [1])
