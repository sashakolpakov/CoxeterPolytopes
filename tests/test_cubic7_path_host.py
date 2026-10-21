from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_path_host import (
    COEFFICIENT,
    DETERMINANT,
    admissible_coefficients_through_selected_coefficient,
    exact_host_checks,
    pair,
    seed_simplex_vertex_types,
    small_admissible_coefficients,
    weyl_vector,
)


def test_cubic7_path_host_arithmetic() -> None:
    exact_host_checks()
    assert COEFFICIENT.norm() == 1
    assert DETERMINANT.norm() == -12488
    assert seed_simplex_vertex_types() == ((0, 1, 2), (), tuple(range(3, 23)))
    assert len(small_admissible_coefficients(10)) == 14
    assert admissible_coefficients_through_selected_coefficient() == (COEFFICIENT,)


def test_cubic7_path_weyl_vector() -> None:
    rho = weyl_vector()
    assert pair(rho, rho) == DETERMINANT / 23
    for index in range(1, 23):
        basis = tuple(Cubic7(1) if j == index else Cubic7(0) for j in range(23))
        assert pair(rho, basis) == 0
