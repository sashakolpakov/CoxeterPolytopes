#!/usr/bin/env python3
"""Certify the cusp pencils over the symbolic rational-carrier parameter.

This module complements ``certify_rational_i1n_cusp_fibrations``.  Its
basic coefficient ring is ``ZZ[a]`` (or exact integer fibres during the
interpolation replay), and the selected sixteenth-centre branch is encoded
by the exact universal jet

    z = x^3 + a*x^7 + x^9 + (2*a^2-1)*x^11
        + (5*a+1)*x^13 + (a-1)*(5*a^2+5*a-1)*x^15 mod x^16.

The final certificate is deliberately separate from the two-fibre archive:
it records a dense-open theorem and the explicit rank minors/second fibres
which make that theorem effective.
"""

from __future__ import annotations

import argparse
import gc
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, gcd, lcm
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from flint import fmpz_mat, nmod_mat
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/rational_i1n_symbolic_cusp_family.json"
CERTIFIER = Path(__file__).resolve()
TWO_FIBRE_CERTIFICATE = (
    ROOT / "data/enumerations/rational_i1n_cusp_fibrations.json"
)
FAMILY_CERTIFICATE = (
    ROOT / "data/enumerations/rational_i1n_j_moduli_carriers.json"
)
SCHEMA = "coxeter.rational_i1n.symbolic_cusp_family.v1"
PIVOT_PRIME = 1_000_003


IntegerPolynomial = tuple[int, ...]


def _polynomial(values: Iterable[int]) -> IntegerPolynomial:
    answer = list(values)
    while answer and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def _polynomial_add(
    left: IntegerPolynomial, right: IntegerPolynomial
) -> IntegerPolynomial:
    return _polynomial(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    )


def _polynomial_scale(
    value: IntegerPolynomial, scalar: int
) -> IntegerPolynomial:
    return _polynomial(scalar * coefficient for coefficient in value)


def _polynomial_multiply(
    left: IntegerPolynomial, right: IntegerPolynomial
) -> IntegerPolynomial:
    if not left or not right:
        return ()
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if not x:
            continue
        for j, y in enumerate(right):
            if y:
                answer[i + j] += x * y
    return _polynomial(answer)


def _polynomial_evaluate(value: IntegerPolynomial, parameter: int) -> int:
    answer = 0
    for coefficient in reversed(value):
        answer = answer * parameter + coefficient
    return answer


def _polynomial_content(values: Iterable[IntegerPolynomial]) -> int:
    answer = 0
    for value in values:
        for coefficient in value:
            answer = gcd(answer, abs(coefficient))
    return answer


def _polynomial_gcd(values: Iterable[IntegerPolynomial]) -> IntegerPolynomial:
    a = sp.symbols("a")
    answer: sp.Poly | None = None
    for value in values:
        if not value:
            continue
        current = sp.Poly(
            sum(coefficient * a**index for index, coefficient in enumerate(value)),
            a,
            domain=sp.ZZ,
        )
        answer = current if answer is None else sp.gcd(answer, current)
        if answer.degree() == 0:
            break
    if answer is None:
        return ()
    primitive = answer.primitive()[1]
    if primitive.LC() < 0:
        primitive = -primitive
    return _polynomial(reversed(tuple(map(int, primitive.all_coeffs()))))


def _polynomial_to_string(value: IntegerPolynomial) -> str:
    a = sp.symbols("a")
    return str(sp.factor(sum(
        coefficient * a**index for index, coefficient in enumerate(value)
    )))


def _sympy_to_integer_polynomial(expression: sp.Expr) -> IntegerPolynomial:
    a = sp.symbols("a")
    polynomial = sp.Poly(sp.expand(expression), a, domain=sp.ZZ)
    return _polynomial(reversed(tuple(map(int, polynomial.all_coeffs()))))


def _polynomial_product(values: Iterable[IntegerPolynomial]) -> IntegerPolynomial:
    answer: IntegerPolynomial = (1,)
    for value in values:
        answer = _polynomial_multiply(answer, value)
    return answer


def selected_s_jet_coefficients(parameter: int) -> dict[int, int]:
    """Return the exact selected-S branch through the sixteenth centre."""

    a = parameter
    return {
        3: 1,
        7: a,
        9: 1,
        11: 2 * a * a - 1,
        13: 5 * a + 1,
        15: (a - 1) * (5 * a * a + 5 * a - 1),
    }


def symbolic_selected_s_jet_coefficients() -> dict[int, IntegerPolynomial]:
    """Return the selected branch in ``ZZ[a][[x]]/(x^16)``."""

    return {
        3: (1,),
        7: (0, 1),
        9: (1,),
        11: (-1, 0, 2),
        13: (1, 5),
        15: (1, -6, 0, 5),
    }


def verify_symbolic_cluster_jet() -> dict[str, object]:
    """Derive the universal jet by fifteen exact strict transforms."""

    a, x, z, u = sp.symbols("a x z u")
    cubic = z - x**3 - a * x * z**2 - z**3
    quartic = sp.expand(z**4 + (x + z) * cubic)
    sextic = sp.expand(cubic**2 + z**2 * quartic)
    q_transform = quartic
    s_transform = sextic
    q_multiplicities = (2, *((1,) * 14))
    s_multiplicities = (*((2,) * 7), *((1,) * 8))
    expected_slopes = (
        0,
        0,
        1,
        0,
        0,
        0,
        a,
        0,
        1,
        0,
        2 * a**2 - 1,
        0,
        5 * a + 1,
        0,
        5 * a**3 - 6 * a + 2,
    )
    replay: list[dict[str, object]] = []

    def minimum_total_degree(expression: sp.Expr) -> int:
        polynomial = sp.Poly(expression, x, z, domain=sp.ZZ.poly_ring(a))
        return min(sum(monomial) for monomial, coefficient in polynomial.terms()
                   if coefficient)

    def blow_up(expression: sp.Expr, slope: sp.Expr, multiplicity: int) -> sp.Expr:
        numerator = sp.Poly(
            sp.expand(expression.subs(z, x * (u + slope))),
            x,
            u,
            domain=sp.ZZ.poly_ring(a),
        )
        terms = numerator.terms()
        if min(monomial[0] for monomial, coefficient in terms if coefficient) < multiplicity:
            raise ArithmeticError("a symbolic strict transform is not divisible")
        return sp.expand(sum(
            coefficient * x ** (monomial[0] - multiplicity) * z ** monomial[1]
            for monomial, coefficient in terms
        ))

    for stage in range(1, 16):
        q_multiplicity = minimum_total_degree(q_transform)
        s_multiplicity = minimum_total_degree(s_transform)
        if (
            q_multiplicity != q_multiplicities[stage - 1]
            or s_multiplicity != s_multiplicities[stage - 1]
        ):
            raise ArithmeticError("the symbolic multiplicity sequence changed")
        if stage == 1:
            slope = sp.Integer(0)
        else:
            polynomial = sp.Poly(
                q_transform, x, z, domain=sp.ZZ.poly_ring(a)
            )
            linear_x = polynomial.coeff_monomial(x)
            linear_z = polynomial.coeff_monomial(z)
            if linear_z not in (1, -1):
                raise ArithmeticError("the selected symbolic branch lost unit slope")
            slope = sp.expand(-linear_x / linear_z)
        if sp.expand(slope - expected_slopes[stage - 1]) != 0:
            raise ArithmeticError(f"the symbolic slope at p{stage} changed")
        replay.append({
            "centre": stage,
            "Q_multiplicity": q_multiplicity,
            "S_multiplicity": s_multiplicity,
            "selected_Q_slope": str(sp.factor(slope)),
        })
        q_transform = blow_up(
            q_transform, slope, q_multiplicities[stage - 1]
        )
        s_transform = blow_up(
            s_transform, slope, s_multiplicities[stage - 1]
        )

    q_exceptional = sp.expand(q_transform.subs(x, 0))
    s_exceptional = sp.expand(s_transform.subs(x, 0))
    if q_exceptional != z or s_exceptional != z + 1:
        raise ArithmeticError("the symbolic p16 exceptional restrictions changed")
    derived_s = {
        exponent: _sympy_to_integer_polynomial(
            slope - (1 if exponent == 15 else 0)
        )
        for exponent, slope in enumerate(expected_slopes, start=1)
        if slope - (1 if exponent == 15 else 0) != 0
    }
    if derived_s != symbolic_selected_s_jet_coefficients():
        raise ArithmeticError("the derived S jet differs from the matrix input")
    return {
        "derivation": "fifteen coefficientwise strict transforms over ZZ[a]",
        "Q_multiplicity_sequence": list(q_multiplicities),
        "S_multiplicity_sequence": list((*s_multiplicities, 1)),
        "post_p15_exceptional_restrictions": {"Q": "z", "S": "z+1"},
        "selected_S_jet_coefficients": {
            str(exponent): list(value) for exponent, value in derived_s.items()
        },
        "replay": replay,
    }


def affine_monomials(degree: int) -> tuple[tuple[int, int], ...]:
    """Use x^i z^j for X^i Z^j Y^(degree-i-j)."""

    return tuple(
        (i, j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    )


def _multiply_univariate(
    left: Sequence[int], right: Sequence[int]
) -> list[int]:
    if not left or not right:
        return []
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if not x:
            continue
        for j, y in enumerate(right):
            if y:
                answer[i + j] += x * y
    while answer and not answer[-1]:
        answer.pop()
    return answer


def _phi_powers(
    coefficients: dict[int, int], centre: int, maximum_power: int
) -> list[list[int]]:
    phi = [0] * centre
    for exponent, value in coefficients.items():
        if exponent < centre:
            phi[exponent] = value
    while phi and not phi[-1]:
        phi.pop()
    powers = [[1]]
    for _ in range(maximum_power):
        powers.append(_multiply_univariate(powers[-1], phi))
    return powers


def integer_jet_condition_matrix(
    degree: int,
    multiplicities: Sequence[int],
    parameter: int,
) -> tuple[
    tuple[tuple[int, int], ...],
    list[list[int]],
    list[tuple[int, int, int, tuple[tuple[int, int], ...]]],
]:
    """Build the same complete-ideal matrix as the v3 certifier, quickly.

    The implementation uses elementary coefficient convolution rather than
    repeatedly expanding SymPy expressions.  It is exact over ``ZZ`` and is
    independently compared with the original builder in the test suite.
    """

    coefficients = selected_s_jet_coefficients(parameter)
    monomials = affine_monomials(degree)
    conditions: dict[tuple[int, int, int], dict[int, int]] = {}
    assigned_before = 0

    for centre, required in enumerate(multiplicities, start=1):
        powers = _phi_powers(coefficients, centre, degree)
        for column, (x_power, z_power) in enumerate(monomials):
            for u_power in range(z_power + 1):
                scalar = comb(z_power, u_power)
                for phi_power, coefficient in enumerate(
                    powers[z_power - u_power]
                ):
                    if not coefficient:
                        continue
                    residual_x_power = (
                        x_power
                        + (centre - 1) * u_power
                        + phi_power
                        - assigned_before
                    )
                    if residual_x_power < 0 or (
                        residual_x_power + u_power < required
                    ):
                        key = (centre, residual_x_power, u_power)
                        row = conditions.setdefault(key, {})
                        row[column] = row.get(column, 0) + scalar * coefficient
        assigned_before += required

    dense: list[list[int]] = []
    sparse: list[tuple[int, int, int, tuple[tuple[int, int], ...]]] = []
    for key in sorted(conditions):
        entries = tuple(
            (column, value)
            for column, value in sorted(conditions[key].items())
            if value
        )
        if not entries:
            continue
        row = [0] * len(monomials)
        for column, value in entries:
            row[column] = value
        dense.append(row)
        sparse.append((*key, entries))
    return monomials, dense, sparse


def _multiply_univariate_over_parameter(
    left: Sequence[IntegerPolynomial],
    right: Sequence[IntegerPolynomial],
) -> list[IntegerPolynomial]:
    if not left or not right:
        return []
    answer: list[IntegerPolynomial] = [()] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if not x:
            continue
        for j, y in enumerate(right):
            if y:
                answer[i + j] = _polynomial_add(
                    answer[i + j], _polynomial_multiply(x, y)
                )
    while answer and not answer[-1]:
        answer.pop()
    return answer


def _symbolic_phi_powers(
    centre: int, maximum_power: int
) -> list[list[IntegerPolynomial]]:
    phi: list[IntegerPolynomial] = [()] * centre
    for exponent, value in symbolic_selected_s_jet_coefficients().items():
        if exponent < centre:
            phi[exponent] = value
    while phi and not phi[-1]:
        phi.pop()
    powers: list[list[IntegerPolynomial]] = [[(1,)]]
    for _ in range(maximum_power):
        powers.append(_multiply_univariate_over_parameter(powers[-1], phi))
    return powers


def symbolic_jet_condition_matrix(
    degree: int,
    multiplicities: Sequence[int],
) -> tuple[
    tuple[tuple[int, int], ...],
    list[list[IntegerPolynomial]],
    list[
        tuple[
            int,
            int,
            int,
            tuple[tuple[int, IntegerPolynomial], ...],
        ]
    ],
]:
    """Build the complete-ideal matrix coefficientwise in ``ZZ[a]``."""

    monomials = affine_monomials(degree)
    conditions: dict[
        tuple[int, int, int], dict[int, IntegerPolynomial]
    ] = {}
    assigned_before = 0

    for centre, required in enumerate(multiplicities, start=1):
        powers = _symbolic_phi_powers(centre, degree)
        for column, (x_power, z_power) in enumerate(monomials):
            for u_power in range(z_power + 1):
                scalar = comb(z_power, u_power)
                for phi_power, coefficient in enumerate(
                    powers[z_power - u_power]
                ):
                    if not coefficient:
                        continue
                    residual_x_power = (
                        x_power
                        + (centre - 1) * u_power
                        + phi_power
                        - assigned_before
                    )
                    if residual_x_power < 0 or (
                        residual_x_power + u_power < required
                    ):
                        key = (centre, residual_x_power, u_power)
                        row = conditions.setdefault(key, {})
                        row[column] = _polynomial_add(
                            row.get(column, ()),
                            _polynomial_scale(coefficient, scalar),
                        )
        assigned_before += required

    dense: list[list[IntegerPolynomial]] = []
    sparse: list[
        tuple[
            int,
            int,
            int,
            tuple[tuple[int, IntegerPolynomial], ...],
        ]
    ] = []
    for key in sorted(conditions):
        entries = tuple(
            (column, value)
            for column, value in sorted(conditions[key].items())
            if value
        )
        if not entries:
            continue
        row: list[IntegerPolynomial] = [()] * len(monomials)
        for column, value in entries:
            row[column] = value
        dense.append(row)
        sparse.append((*key, entries))
    return monomials, dense, sparse


def polynomial_values_from_sparse_terms(
    terms: Iterable[Sequence[int]], monomials: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    """Expand ``[i,j,c]`` sparse terms into the monomial coefficient vector."""

    lookup = {(i, j): index for index, (i, j) in enumerate(monomials)}
    values = [0] * len(monomials)
    for i, j, coefficient in terms:
        values[lookup[(int(i), int(j))]] = int(coefficient)
    return tuple(values)


CASES: dict[str, dict[str, object]] = {
    "genus_5": {
        "degree": 11,
        "multiplicities": (*((3,) * 13), 2),
        "expected_rank": 76,
        "wall_factors": {"C": 1, "Q": 2},
        "rank_minor": (-3, -1),  # -(a+3)
        "interpolation_scale": (3, 1),
        "expected_section_max_degree": 4,
        "arithmetic_genus": 5,
    },
    "genus_10": {
        "degree": 21,
        "multiplicities": (*((7,) * 5), 6, *((4,) * 10)),
        "expected_rank": 251,
        "wall_factors": {"L": 1, "Q": 2, "S": 2},
        "rank_minor": (1,),
        "interpolation_scale": (1,),
        "expected_section_max_degree": 7,
        "arithmetic_genus": 10,
    },
    "genus_11": {
        "degree": 23,
        "multiplicities": (*((7,) * 8), 5, *((4,) * 7)),
        "expected_rank": 298,
        "wall_factors": {"C": 1, "Q": 2, "S": 2},
        "rank_minor": (-1, 1),  # a-1
        "interpolation_scale": (-1, 1),
        "expected_section_max_degree": 8,
        "arithmetic_genus": 11,
    },
    "fourth_multiple": {
        "degree": 16,
        "multiplicities": (4,) * 16,
        "expected_rank": 151,
        "wall_factors": {"Q": 4},
        "rank_minor": (-12, -12),  # -12(a+1)
        "interpolation_scale": (1, 1),
        "expected_section_max_degree": 6,
        "arithmetic_genus": 9,
    },
}


ALTERNATE_CHARTS: dict[str, dict[str, object]] = {
    "genus_5": {
        "base_parameter": -3,
        "rank_minor": (-3,),
        "interpolation_scale": (1,),
        "expected_section_max_degree": 5,
    },
    "genus_11": {
        "base_parameter": 1,
        "rank_minor": (-7,),
        "interpolation_scale": (1,),
        "expected_section_max_degree": 9,
    },
    "fourth_multiple": {
        "base_parameter": -1,
        "rank_minor": (-26, -20, 2),  # 2(a^2-10a-13)
        "interpolation_scale": (-26, -20, 2),
        "expected_section_max_degree": 7,
    },
}


RIGID_MULTIPLES: dict[str, dict[str, object]] = {
    "F": {
        "degree": 4,
        "multiplicities": (1,) * 16,
        "expected_rank": 14,
        "rank_minor": (1,),
    },
    "2F": {
        "degree": 8,
        "multiplicities": (2,) * 16,
        "expected_rank": 44,
        "rank_minor": (-2,),
    },
    "3F": {
        "degree": 12,
        "multiplicities": (3,) * 16,
        "expected_rank": 90,
        "rank_minor": (2,),
    },
}


def _pivot_minor_indices(
    rows: Sequence[Sequence[int]], expected_rank: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Select a deterministic full-rank square minor modulo a fixed prime."""

    matrix = nmod_mat(rows, PIVOT_PRIME)
    reduced_transpose, rank = matrix.transpose().rref()
    if rank != expected_rank:
        raise ArithmeticError(f"modular rank changed: {rank} != {expected_rank}")
    selected_rows: list[int] = []
    for row in range(rank):
        pivot = next(
            column
            for column in range(reduced_transpose.ncols())
            if int(reduced_transpose[row, column])
        )
        selected_rows.append(pivot)

    restricted = nmod_mat([rows[index] for index in selected_rows], PIVOT_PRIME)
    reduced, second_rank = restricted.rref()
    if second_rank != expected_rank:
        raise ArithmeticError("the deterministic row restriction lost rank")
    selected_columns: list[int] = []
    for row in range(second_rank):
        pivot = next(
            column
            for column in range(reduced.ncols())
            if int(reduced[row, column])
        )
        selected_columns.append(pivot)
    return tuple(selected_rows), tuple(selected_columns)


def _minor_determinant(
    symbolic_rows: Sequence[Sequence[IntegerPolynomial]],
    selected_rows: Sequence[int],
    selected_columns: Sequence[int],
) -> IntegerPolynomial:
    a = sp.symbols("a")
    entries = [
        [
            sum(
                coefficient * a**exponent
                for exponent, coefficient in enumerate(
                    symbolic_rows[row][column]
                )
            )
            for column in selected_columns
        ]
        for row in selected_rows
    ]
    matrix = sp.polys.matrices.DomainMatrix.from_list_sympy(
        len(selected_rows), len(selected_columns), entries
    )
    return _sympy_to_integer_polynomial(matrix.det().as_expr())


def _normalized_kernel_vector(
    rows: Sequence[Sequence[int]], free_columns: Sequence[int]
) -> tuple[Fraction, ...]:
    """Normalize the second kernel vector to free coordinates ``(0,1)``."""

    nullspace, nullity = fmpz_mat(rows).nullspace()
    if nullity != 2:
        raise ArithmeticError(f"expected a two-dimensional kernel, got {nullity}")
    first, second = map(int, free_columns)
    b00 = int(nullspace[first, 0])
    b01 = int(nullspace[first, 1])
    b10 = int(nullspace[second, 0])
    b11 = int(nullspace[second, 1])
    determinant = b00 * b11 - b01 * b10
    if determinant == 0:
        raise ArithmeticError("the chosen free-coordinate normalization failed")
    coefficient_0 = Fraction(-b01, determinant)
    coefficient_1 = Fraction(b00, determinant)
    answer = tuple(
        Fraction(int(nullspace[row, 0])) * coefficient_0
        + Fraction(int(nullspace[row, 1])) * coefficient_1
        for row in range(nullspace.nrows())
    )
    if (answer[first], answer[second]) != (0, 1):
        raise ArithmeticError("the canonical kernel normalization changed")
    return answer


def _interpolate_second_section(
    degree: int,
    multiplicities: Sequence[int],
    free_columns: Sequence[int],
    scale: IntegerPolynomial,
    expected_max_degree: int,
) -> tuple[IntegerPolynomial, ...]:
    """Reconstruct and primitively normalize a section in ``ZZ[a][x,z]``."""

    a = sp.symbols("a")
    scale_roots = {
        parameter
        for parameter in range(-9, 10)
        if _polynomial_evaluate(scale, parameter) == 0
    }
    training_parameters = tuple(
        parameter
        for parameter in range(-6, 7)
        if parameter not in scale_roots
    )
    if len(training_parameters) <= expected_max_degree:
        raise ArithmeticError("too few interpolation fibres")

    values: list[tuple[Fraction, ...]] = []
    monomials: tuple[tuple[int, int], ...] | None = None
    for parameter in training_parameters:
        current_monomials, rows, _sparse = integer_jet_condition_matrix(
            degree, multiplicities, parameter
        )
        if monomials is None:
            monomials = current_monomials
        elif current_monomials != monomials:
            raise ArithmeticError("the monomial basis changed")
        vector = _normalized_kernel_vector(rows, free_columns)
        scalar = _polynomial_evaluate(scale, parameter)
        values.append(tuple(scalar * entry for entry in vector))

    assert monomials is not None
    rational_polynomials: list[sp.Poly] = []
    for column in range(len(monomials)):
        points = [
            (
                parameter,
                sp.Rational(
                    values[index][column].numerator,
                    values[index][column].denominator,
                ),
            )
            for index, parameter in enumerate(training_parameters)
        ]
        polynomial = sp.Poly(sp.interpolate(points, a), a, domain=sp.QQ)
        if polynomial.degree() > expected_max_degree:
            raise ArithmeticError(
                f"interpolated section degree {polynomial.degree()} is too large"
            )
        rational_polynomials.append(polynomial)

    common_denominator = 1
    for polynomial in rational_polynomials:
        for coefficient in polynomial.all_coeffs():
            common_denominator = lcm(common_denominator, int(coefficient.q))
    integer_polynomials = tuple(
        _polynomial(reversed(tuple(
            int(coefficient * common_denominator)
            for coefficient in polynomial.all_coeffs()
        )))
        for polynomial in rational_polynomials
    )
    content = _polynomial_content(integer_polynomials)
    if content == 0:
        raise ArithmeticError("the interpolated section vanished")
    integer_polynomials = tuple(
        _polynomial(coefficient // content for coefficient in polynomial)
        for polynomial in integer_polynomials
    )
    first_coefficient = next(
        coefficient
        for polynomial in integer_polynomials
        for coefficient in reversed(polynomial)
        if coefficient
    )
    if first_coefficient < 0:
        integer_polynomials = tuple(
            _polynomial_scale(polynomial, -1)
            for polynomial in integer_polynomials
        )
    return integer_polynomials


def _verify_symbolic_kernel(
    symbolic_rows: Sequence[Sequence[IntegerPolynomial]],
    section: Sequence[IntegerPolynomial],
) -> None:
    for row in symbolic_rows:
        value: IntegerPolynomial = ()
        for entry, coefficient in zip(row, section):
            if entry and coefficient:
                value = _polynomial_add(
                    value, _polynomial_multiply(entry, coefficient)
                )
        if value:
            raise ArithmeticError("an interpolated section left the ZZ[a] kernel")


def _carrier_forms() -> tuple[sp.Symbol, sp.Symbol, sp.Symbol, dict[str, sp.Expr]]:
    a, x, z = sp.symbols("a x z")
    cubic = z - x**3 - a * x * z**2 - z**3
    quartic = sp.expand(z**4 + (x + z) * cubic)
    sextic = sp.expand(cubic**2 + z**2 * quartic)
    return a, x, z, {"L": z, "C": cubic, "Q": quartic, "S": sextic}


def _expression_coefficient_vector(
    expression: sp.Expr, monomials: Sequence[tuple[int, int]]
) -> tuple[IntegerPolynomial, ...]:
    a, x, z = sp.symbols("a x z")
    polynomial = sp.Poly(sp.expand(expression), x, z, domain=sp.ZZ.poly_ring(a))
    return tuple(
        _sympy_to_integer_polynomial(polynomial.coeff_monomial(x**i * z**j))
        for i, j in monomials
    )


def _wall_data(
    factors: Mapping[str, int], monomials: Sequence[tuple[int, int]]
) -> tuple[sp.Expr, tuple[IntegerPolynomial, ...]]:
    _a, _x, _z, forms = _carrier_forms()
    expression = sp.expand(sp.prod(
        forms[label] ** exponent for label, exponent in factors.items()
    ))
    return expression, _expression_coefficient_vector(expression, monomials)


def _section_expression(
    section: Sequence[IntegerPolynomial],
    monomials: Sequence[tuple[int, int]],
) -> sp.Expr:
    a, x, z = sp.symbols("a x z")
    return sp.expand(sum(
        coefficient * a**exponent * x**i * z**j
        for polynomial, (i, j) in zip(section, monomials)
        for exponent, coefficient in enumerate(polynomial)
    ))


def _actual_multiplicity_gate(
    section: Sequence[IntegerPolynomial],
    monomials: Sequence[tuple[int, int]],
    multiplicities: Sequence[int],
) -> tuple[tuple[int, ...], tuple[IntegerPolynomial, ...]]:
    maximum_z_power = max(j for _i, j in monomials)
    assigned_before = 0
    actual: list[int] = []
    gates: list[IntegerPolynomial] = []
    for centre, required in enumerate(multiplicities, start=1):
        powers = _symbolic_phi_powers(centre, maximum_z_power)
        terms: dict[tuple[int, int], IntegerPolynomial] = {}
        for scalar, (x_power, z_power) in zip(section, monomials):
            if not scalar:
                continue
            for u_power in range(z_power + 1):
                binomial = comb(z_power, u_power)
                for phi_power, coefficient in enumerate(
                    powers[z_power - u_power]
                ):
                    if not coefficient:
                        continue
                    key = (
                        x_power + (centre - 1) * u_power,
                        u_power,
                    )
                    # phi_power is its x exponent.
                    key = (key[0] + phi_power, key[1])
                    terms[key] = _polynomial_add(
                        terms.get(key, ()),
                        _polynomial_scale(
                            _polynomial_multiply(scalar, coefficient), binomial
                        ),
                    )
        nonzero_terms = {key: value for key, value in terms.items() if value}
        if not nonzero_terms:
            raise ArithmeticError("the symbolic plane section vanished")
        minimum = min(
            x_power - assigned_before + u_power
            for x_power, u_power in nonzero_terms
        )
        if minimum != required:
            raise ArithmeticError(
                f"generic multiplicity {minimum} != assigned {required} at p{centre}"
            )
        leading = [
            value
            for (x_power, u_power), value in nonzero_terms.items()
            if x_power - assigned_before + u_power == minimum
        ]
        gate = _polynomial_gcd(leading)
        if not gate:
            raise ArithmeticError("the leading multiplicity gate vanished")
        actual.append(minimum)
        gates.append(gate)
        assigned_before += minimum
    return tuple(actual), tuple(gates)


def _coprimality_gates(
    section_expression: sp.Expr, factor_labels: Iterable[str]
) -> dict[str, IntegerPolynomial]:
    a, x, z, forms = _carrier_forms()
    coefficient_ring = sp.ZZ.poly_ring(a, z)
    answer: dict[str, IntegerPolynomial] = {}
    for label in sorted(set(factor_labels)):
        factor = forms[label]
        if label == "L":
            remainder = sp.expand(section_expression.subs(z, 0))
        else:
            dividend = sp.Poly(section_expression, x, domain=coefficient_ring)
            divisor = sp.Poly(factor, x, domain=coefficient_ring)
            if abs(divisor.LC()) != 1:
                raise ArithmeticError(f"{label} is not monic in x")
            remainder = dividend.rem(divisor).as_expr()
        coefficient_polynomial = sp.Poly(
            sp.expand(remainder), x, z, domain=sp.ZZ.poly_ring(a)
        )
        coefficients = [
            _sympy_to_integer_polynomial(coefficient)
            for (_monomial, coefficient) in coefficient_polynomial.terms()
        ]
        gate = _polynomial_gcd(coefficients)
        if not gate:
            raise ArithmeticError(f"the section is generically divisible by {label}")
        answer[label] = gate
    return answer


def _gate_lcm(values: Iterable[IntegerPolynomial]) -> IntegerPolynomial:
    a = sp.symbols("a")
    answer = sp.Poly(1, a, domain=sp.QQ)
    for value in values:
        if not value:
            raise ArithmeticError("a chart gate is the zero polynomial")
        current = sp.Poly(sum(
            coefficient * a**index
            for index, coefficient in enumerate(value)
        ), a, domain=sp.QQ)
        answer = sp.lcm(answer, current.sqf_part())
    answer = answer.monic()
    denominator = 1
    for coefficient in answer.all_coeffs():
        denominator = lcm(denominator, int(coefficient.q))
    integral = sp.Poly(answer.as_expr() * denominator, a, domain=sp.ZZ)
    integral = integral.primitive()[1]
    if integral.LC() < 0:
        integral = -integral
    return _polynomial(reversed(tuple(map(int, integral.all_coeffs()))))


def _sparse_section_terms(
    section: Sequence[IntegerPolynomial],
    monomials: Sequence[tuple[int, int]],
) -> list[list[object]]:
    return [
        [i, j, list(polynomial)]
        for polynomial, (i, j) in zip(section, monomials)
        if polynomial
    ]


def _digest_json(value: object) -> str:
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def _chart_record(
    name: str,
    case: Mapping[str, object],
    chart_name: str,
    chart: Mapping[str, object],
    monomials: Sequence[tuple[int, int]],
    symbolic_rows: Sequence[Sequence[IntegerPolynomial]],
    symbolic_sparse: Sequence[
        tuple[
            int,
            int,
            int,
            tuple[tuple[int, IntegerPolynomial], ...],
        ]
    ],
    wall_vector: Sequence[IntegerPolynomial],
    wall_expression: sp.Expr,
) -> dict[str, object]:
    degree = int(case["degree"])
    multiplicities = tuple(map(int, case["multiplicities"]))
    expected_rank = int(case["expected_rank"])
    base_parameter = int(chart.get("base_parameter", 0))
    _base_monomials, base_rows, _base_sparse = integer_jet_condition_matrix(
        degree, multiplicities, base_parameter
    )
    selected_rows, selected_columns = _pivot_minor_indices(
        base_rows, expected_rank
    )
    free_columns = tuple(sorted(
        set(range(len(monomials))) - set(selected_columns)
    ))
    if len(free_columns) != 2:
        raise ArithmeticError("a moving-system minor did not leave two columns")
    determinant = _minor_determinant(
        symbolic_rows, selected_rows, selected_columns
    )
    expected_determinant = tuple(map(int, chart["rank_minor"]))
    if determinant != expected_determinant:
        raise ArithmeticError(
            f"{name}/{chart_name} minor {determinant} != {expected_determinant}"
        )

    section = _interpolate_second_section(
        degree,
        multiplicities,
        free_columns,
        tuple(map(int, chart["interpolation_scale"])),
        int(chart["expected_section_max_degree"]),
    )
    _verify_symbolic_kernel(symbolic_rows, section)
    independence = _polynomial_add(
        _polynomial_multiply(
            wall_vector[free_columns[0]], section[free_columns[1]]
        ),
        _polynomial_scale(
            _polynomial_multiply(
                wall_vector[free_columns[1]], section[free_columns[0]]
            ),
            -1,
        ),
    )
    if not independence:
        raise ArithmeticError("the symbolic section is the wall member")
    actual, multiplicity_gates = _actual_multiplicity_gate(
        section, monomials, multiplicities
    )
    if actual != multiplicities:
        raise ArithmeticError("the symbolic section multiplicities changed")
    section_expression = _section_expression(section, monomials)
    coprimality = _coprimality_gates(
        section_expression, case["wall_factors"].keys()
    )
    chart_gate = _gate_lcm((
        determinant,
        independence,
        *multiplicity_gates,
        *coprimality.values(),
    ))
    section_terms = _sparse_section_terms(section, monomials)
    selected_row_keys = [
        list(symbolic_sparse[index][:3]) for index in selected_rows
    ]
    selected_column_monomials = [
        list(monomials[index]) for index in selected_columns
    ]
    return {
        "chart": chart_name,
        "pivot_parameter": base_parameter,
        "pivot_prime": PIVOT_PRIME,
        "rank_minor_size": expected_rank,
        "rank_minor_determinant_coefficients_ascending": list(determinant),
        "rank_minor_determinant_factorization": _polynomial_to_string(determinant),
        "selected_row_keys": selected_row_keys,
        "selected_row_keys_sha256": _digest_json(selected_row_keys),
        "selected_column_monomials": selected_column_monomials,
        "selected_column_monomials_sha256": _digest_json(
            selected_column_monomials
        ),
        "free_column_indices": list(free_columns),
        "free_column_monomials": [list(monomials[index]) for index in free_columns],
        "second_section_sparse_terms": section_terms,
        "second_section_sparse_terms_sha256": _digest_json(section_terms),
        "second_section_term_count": len(section_terms),
        "second_section_parameter_degree": max(
            len(polynomial) - 1 for polynomial in section if polynomial
        ),
        "second_section_in_ZZ_a_jet_kernel": True,
        "wall_second_section_independence_determinant": (
            _polynomial_to_string(independence)
        ),
        "actual_multiplicity_sequence_over_Q_a": list(actual),
        "actual_multiplicity_leading_gcds": [
            _polynomial_to_string(value) for value in multiplicity_gates
        ],
        "coprimality_remainder_coefficient_gcds": {
            label: _polynomial_to_string(value)
            for label, value in coprimality.items()
        },
        "chart_gate_coefficients_ascending": list(chart_gate),
        "chart_gate_factorization": _polynomial_to_string(chart_gate),
        "chart_conclusion": (
            "where the chart gate is nonzero, the maximal minor gives the "
            "rank lower bound; the wall member and polynomial second section "
            "give the matching kernel upper bound and a coprime exact-"
            "multiplicity pair, hence a base-point-free pencil"
        ),
    }


def _moving_case_record(name: str, case: Mapping[str, object]) -> dict[str, object]:
    degree = int(case["degree"])
    multiplicities = tuple(map(int, case["multiplicities"]))
    monomials, symbolic_rows, symbolic_sparse = symbolic_jet_condition_matrix(
        degree, multiplicities
    )
    wall_expression, wall_vector = _wall_data(case["wall_factors"], monomials)
    _verify_symbolic_kernel(symbolic_rows, wall_vector)
    primary = {
        "base_parameter": 0,
        "rank_minor": case["rank_minor"],
        "interpolation_scale": case["interpolation_scale"],
        "expected_section_max_degree": case["expected_section_max_degree"],
    }
    chart_inputs = [("primary", primary)]
    if name in ALTERNATE_CHARTS:
        chart_inputs.append(("alternate", ALTERNATE_CHARTS[name]))
    charts = [
        _chart_record(
            name,
            case,
            chart_name,
            chart,
            monomials,
            symbolic_rows,
            symbolic_sparse,
            wall_vector,
            wall_expression,
        )
        for chart_name, chart in chart_inputs
    ]
    chart_gates = [
        tuple(map(int, chart["chart_gate_coefficients_ascending"]))
        for chart in charts
    ]
    cover_gcd = _polynomial_gcd(chart_gates)
    if cover_gcd != (1,):
        raise ArithmeticError(f"{name} chart gates do not cover A1: {cover_gcd}")
    return {
        "case": name,
        "degree": degree,
        "assigned_multiplicities": list(multiplicities),
        "arithmetic_genus": int(case["arithmetic_genus"]),
        "monomial_count": len(monomials),
        "symbolic_condition_row_count": len(symbolic_rows),
        "expected_rank": int(case["expected_rank"]),
        "expected_h0": 2,
        "wall_factorization": "*".join(
            label if exponent == 1 else f"{label}^{exponent}"
            for label, exponent in case["wall_factors"].items()
        ),
        "wall_member_in_ZZ_a_jet_kernel": True,
        "charts": charts,
        "chart_gate_gcd": _polynomial_to_string(cover_gcd),
        "all_parameter_values_covered": True,
        "base_point_free_for_every_parameter_where_the_cluster_is_a_carrier": True,
    }


def _rigid_multiple_record(name: str, case: Mapping[str, object]) -> dict[str, object]:
    degree = int(case["degree"])
    multiplicities = tuple(map(int, case["multiplicities"]))
    expected_rank = int(case["expected_rank"])
    monomials, base_rows, _base_sparse = integer_jet_condition_matrix(
        degree, multiplicities, 0
    )
    selected_rows, selected_columns = _pivot_minor_indices(
        base_rows, expected_rank
    )
    symbolic_monomials, symbolic_rows, symbolic_sparse = (
        symbolic_jet_condition_matrix(degree, multiplicities)
    )
    if symbolic_monomials != monomials:
        raise ArithmeticError("the rigid monomial basis changed")
    determinant = _minor_determinant(
        symbolic_rows, selected_rows, selected_columns
    )
    expected_determinant = tuple(map(int, case["rank_minor"]))
    if determinant != expected_determinant:
        raise ArithmeticError(f"the {name} rank minor changed")
    _a, _x, _z, forms = _carrier_forms()
    multiple = degree // 4
    wall_expression = sp.expand(forms["Q"] ** multiple)
    wall_vector = _expression_coefficient_vector(wall_expression, monomials)
    _verify_symbolic_kernel(symbolic_rows, wall_vector)
    selected_row_keys = [
        list(symbolic_sparse[index][:3]) for index in selected_rows
    ]
    selected_column_monomials = [
        list(monomials[index]) for index in selected_columns
    ]
    return {
        "multiple": name,
        "degree": degree,
        "assigned_multiplicities": list(multiplicities),
        "monomial_count": len(monomials),
        "symbolic_condition_row_count": len(symbolic_rows),
        "rank_minor_size": expected_rank,
        "rank_minor_determinant": _polynomial_to_string(determinant),
        "selected_row_keys": selected_row_keys,
        "selected_row_keys_sha256": _digest_json(selected_row_keys),
        "selected_column_monomials": selected_column_monomials,
        "selected_column_monomials_sha256": _digest_json(
            selected_column_monomials
        ),
        "wall_member": f"Q^{multiple}",
        "wall_member_in_ZZ_a_jet_kernel": True,
        "h0_for_every_parameter": 1,
        "moves": False,
    }


def _file_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _payload_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_payload_sha256", None)
    return _digest_json(payload)


def build_certificate() -> dict[str, object]:
    cluster_replay = verify_symbolic_cluster_jet()
    moving: dict[str, object] = {}
    for name, case in CASES.items():
        moving[name] = _moving_case_record(name, case)
        gc.collect()
    rigid: dict[str, object] = {}
    for name, case in RIGID_MULTIPLES.items():
        rigid[name] = _rigid_multiple_record(name, case)
        gc.collect()

    if not all(
        record["all_parameter_values_covered"] for record in moving.values()
    ):
        raise ArithmeticError("a moving system lost parameter coverage")
    if [rigid[name]["h0_for_every_parameter"] for name in ("F", "2F", "3F")] != [1, 1, 1]:
        raise ArithmeticError("the exact-index lower bound changed")
    fourth = moving["fourth_multiple"]
    if fourth["expected_h0"] != 2:
        raise ArithmeticError("the fourth multiple stopped moving")

    document: dict[str, object] = {
        "schema": SCHEMA,
        "claim": (
            "On every parameter of the pre-existing nonempty carrier-open "
            "locus, all 28 P12--P16 cusp classes are semiample. Twenty-seven "
            "move primitively and F=4H-E1-...-E16 has exact semiample index "
            "four. There is no additional cusp-specific exceptional "
            "parameter polynomial: two explicit ZZ[a] rank/section charts "
            "have coprime gates in every nontrivial case."
        ),
        "scope": {
            "base_locus": (
                "the carrier-open subset from the non-isotrivial family theorem, "
                "where C_a,Q_a,S_a are geometrically integral and the complete "
                "nef chambers P12--P16 persist"
            ),
            "weierstrass_smoothness_factor": "4*a^3+27",
            "additional_cusp_exceptional_polynomial": "1",
            "new_abstract_hyperbolic_polytope": False,
            "worldwide_priority": "open",
        },
        "symbolic_selected_branch_jet": {
            "modulus": "x^16",
            "expression": (
                "x^3+a*x^7+x^9+(2*a^2-1)*x^11+(5*a+1)*x^13+"
                "(a-1)*(5*a^2+5*a-1)*x^15"
            ),
            "coefficients": {
                str(exponent): list(value)
                for exponent, value in symbolic_selected_s_jet_coefficients().items()
            },
            "coefficient_convention": (
                "lists are coefficients in ascending powers of a"
            ),
            "strict_transform_certificate": cluster_replay,
        },
        "moving_systems": moving,
        "rigid_lower_multiples": rigid,
        "totals": {
            "cusps": 28,
            "primitive_base_point_free_pencils": 27,
            "exact_index_four_pencils": 1,
            "additional_parameter_exceptions": 0,
        },
        "proof_logic": {
            "old_22": (
                "the 16 multi-affine-component pencils, five line pencils, "
                "and P13 quartic pencil use parameter-independent wall cycles "
                "and apply throughout the carrier-open locus"
            ),
            "new_five": (
                "coefficientwise ZZ[a] jet kernels provide two exact charts; "
                "their gate gcd is one, and in each valid chart the second "
                "member has exact cluster multiplicities and is coprime to "
                "every irreducible wall factor"
            ),
            "last_cusp": (
                "constant nonzero maximal minors give h0(F)=h0(2F)=h0(3F)=1 "
                "for every parameter; two coprime fourth-multiple charts give "
                "h0(4F)=2 and a base-point-free pencil for every parameter"
            ),
            "connected_genus_nine_fibre": (
                "on a rational surface the Stein base has genus zero; h0(4F)=2 "
                "forces the finite Stein degree to be one, and adjunction gives "
                "arithmetic genus nine"
            ),
        },
        "source_provenance": {
            "certifier": {
                "path": str(CERTIFIER.relative_to(ROOT)),
                "sha256": _file_digest(CERTIFIER),
            },
            "two_fibre_certificate": {
                "path": str(TWO_FIBRE_CERTIFICATE.relative_to(ROOT)),
                "sha256": _file_digest(TWO_FIBRE_CERTIFICATE),
            },
            "family_certificate": {
                "path": str(FAMILY_CERTIFICATE.relative_to(ROOT)),
                "sha256": _file_digest(FAMILY_CERTIFICATE),
            },
        },
    }
    document["certificate_payload_sha256"] = _payload_digest(document)
    return document


def validate_archive(path: Path = OUTPUT) -> dict[str, object]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != SCHEMA:
        raise ArithmeticError("wrong symbolic-cusp-family schema")
    if document.get("certificate_payload_sha256") != _payload_digest(document):
        raise ArithmeticError("symbolic-cusp-family payload digest changed")
    if document.get("totals") != {
        "additional_parameter_exceptions": 0,
        "cusps": 28,
        "exact_index_four_pencils": 1,
        "primitive_base_point_free_pencils": 27,
    }:
        raise ArithmeticError("symbolic-cusp-family totals changed")
    if document["scope"]["additional_cusp_exceptional_polynomial"] != "1":
        raise ArithmeticError("the cusp parameter cover is no longer complete")
    for name, record in document["moving_systems"].items():
        if not (
            record["chart_gate_gcd"] == "1"
            and record["all_parameter_values_covered"]
            and record[
                "base_point_free_for_every_parameter_where_the_cluster_is_a_carrier"
            ]
        ):
            raise ArithmeticError(f"the {name} two-chart cover changed")
    provenance = document["source_provenance"]
    expected = {
        "certifier": CERTIFIER,
        "two_fibre_certificate": TWO_FIBRE_CERTIFICATE,
        "family_certificate": FAMILY_CERTIFICATE,
    }
    for label, source in expected.items():
        if provenance[label] != {
            "path": str(source.relative_to(ROOT)),
            "sha256": _file_digest(source),
        }:
            raise ArithmeticError(f"stale symbolic-family source: {label}")
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    if args.check_archive:
        document = validate_archive(args.output)
    else:
        document = build_certificate()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({
        "payload_sha256": document["certificate_payload_sha256"],
        "additional_cusp_exceptional_polynomial": document["scope"][
            "additional_cusp_exceptional_polynomial"
        ],
        "totals": document["totals"],
    }, sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
