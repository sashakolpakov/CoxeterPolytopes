from fractions import Fraction

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO
from coxeter_polytopes.known_examples.bugaenko import (
    bugaenko_basis_gram,
    bugaenko_interior_point,
    bugaenko_simple_roots,
)
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse, _mat_vec
from coxeter_polytopes.search.h11_vinberg import (
    _standard_e8_gram,
    h11_basis_gram,
    run_h11_vinberg_search,
)
from coxeter_polytopes.verify.chamber import (
    _enumerate_extreme_rays,
    _enumerate_extreme_rays_incremental,
)


def _compare(roots, gram, interior, dimension: int) -> None:
    assert _enumerate_extreme_rays_incremental(
        roots, gram, interior, dimension
    ) == _enumerate_extreme_rays(roots, gram, interior, dimension)


def test_incremental_matches_exhaustive_on_triangle_fixtures() -> None:
    roots = [
        [Golden(1), ZERO, ZERO],
        [ZERO, Golden(1), ZERO],
        [ZERO, ZERO, Golden(1)],
    ]
    _compare(
        roots,
        [[Golden(2), Golden(-1), -PHI], [-1, 2, -1], [-PHI, -1, 2]],
        [Golden(1), Golden(Fraction(1, 2)), Golden(1)],
        2,
    )
    _compare(
        roots,
        [[Golden(2), Golden(-1), Golden(-2)], [-1, 2, -1], [-2, -1, 2]],
        [Golden(1), Golden(Fraction(1, 2)), Golden(1)],
        2,
    )


def test_incremental_matches_exhaustive_on_bugaenko_h8() -> None:
    _compare(
        bugaenko_simple_roots(),
        bugaenko_basis_gram(),
        bugaenko_interior_point(),
        8,
    )


def test_incremental_matches_exhaustive_on_h11_bound_five() -> None:
    search = run_h11_vinberg_search(5)
    inverse_e8 = _fraction_inverse(_standard_e8_gram())
    rho = _mat_vec(inverse_e8, [Golden(1)] * 8)
    inward = (ZERO, Golden(1), *rho)
    interior = tuple(
        search.control_vector[i] + inward[i] / 10 for i in range(10)
    )
    _compare(search.roots, h11_basis_gram(), interior, 9)
