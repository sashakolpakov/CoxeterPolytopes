#!/usr/bin/env python3
"""Certify the isotropic cusp cycles and fibrations of the P12--P16 surfaces.

The complete Mori cones of the rational surfaces X_12,...,X_16 are already
certified in ``certify_rational_i1n_inverse_design.py``.  This replay works on
their dual nef cones.  It enumerates every extreme ray exactly, extracts the
isotropic rays, reconstructs every connected affine wall component and its
primitive null cycle, and records the cases where two disjoint effective
cycles prove semiampleness.

The Lorentzian convention is the one used by the parent certificate:
``diag(-1,+1,...,+1)``.  A divisor ``dH-sum(m_i E_i)`` is stored as
``(d,-m_1,...,-m_n)``.  Surface intersection is the negative of the stored
Lorentzian pairing.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import gc
from hashlib import sha256
import json
from math import comb, gcd
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp

from scripts.certify_rational_i1n_inverse_design import (
    VINBERG_VOLUME_DATA,
    _extreme_rays_via_cdd,
    _primitive_integer_vector,
    canonical_pairing,
    enumerate_simple_roots,
    lorentz_pair,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/rational_i1n_cusp_fibrations.json"
A1_WORKER_OUTPUT = ROOT / "data/enumerations/rational_i1n_a1_cusp_jet_workers.json"
PARENT_JSON = ROOT / "data/enumerations/rational_i1n_inverse_design.json"
PARENT_SCRIPT = ROOT / "scripts/certify_rational_i1n_inverse_design.py"
FAMILY_JSON = ROOT / "data/enumerations/rational_i1n_j_moduli_carriers.json"
FAMILY_SCRIPT = ROOT / "scripts/certify_rational_i1n_j_moduli_carriers.py"
CERTIFIER = Path(__file__).resolve()
SCHEMA = "coxeter.rational_i1n.cusp_fibrations.v3"
A1_WORKER_SCHEMA = "coxeter.rational_i1n.a1_cusp_jet_workers.v1"
DIMENSIONS = tuple(range(12, 17))


def _plane_forms(
    parameter: int,
) -> tuple[sp.Symbol, sp.Symbol, sp.Expr, sp.Expr, sp.Expr]:
    """Return the affine carrier forms in the chart ``Y=1``."""

    x, z = sp.symbols("x z")
    cubic = z - x**3 - parameter * x * z**2 - z**3
    quartic = sp.expand(z**4 + (x + z) * cubic)
    sextic = sp.expand(cubic**2 + z**2 * quartic)
    return x, z, cubic, quartic, sextic


def _a0_plane_forms() -> tuple[sp.Symbol, sp.Symbol, sp.Expr, sp.Expr, sp.Expr]:
    """Compatibility wrapper for the original ``a=0`` carrier."""

    return _plane_forms(0)


def _minimum_total_degree(polynomial: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> int:
    terms = sp.Poly(polynomial, x, z).terms()
    if not terms:
        raise ArithmeticError("a strict transform vanished")
    return min(sum(monomial) for monomial, coefficient in terms if coefficient)


@lru_cache(maxsize=2)
def _cluster_jet_chart(parameter: int) -> dict[str, object]:
    """Derive the curvilinear cluster coordinates by exact blow-up replay.

    The first fifteen centres follow the contact-eleven branch of ``Q``.
    After the fifteenth blow-up, the selected branch of ``S`` is at coordinate
    ``-1`` while the branch of ``Q`` is at zero.  Thus the coefficients below
    are not guessed Taylor data: they are recovered from the defining plane
    equations by fifteen exact strict-transform calculations.
    """

    if parameter not in (0, 1):
        raise ValueError("the frozen explicit jet replay covers a=0 and a=1")
    x, z, _cubic, quartic, sextic = _plane_forms(parameter)
    u = sp.symbols("u")
    q_transform = quartic
    s_transform = sextic
    q_multiplicities = (2, *((1,) * 14))
    s_multiplicities = (*((2,) * 7), *((1,) * 8))
    slopes: list[sp.Rational] = []
    replay: list[dict[str, object]] = []

    for stage in range(1, 16):
        q_multiplicity = _minimum_total_degree(q_transform, x, z)
        s_multiplicity = _minimum_total_degree(s_transform, x, z)
        if q_multiplicity != q_multiplicities[stage - 1]:
            raise ArithmeticError("the Q branch multiplicity sequence changed")
        if s_multiplicity != s_multiplicities[stage - 1]:
            raise ArithmeticError("the S branch multiplicity sequence changed")

        if stage == 1:
            # Q has tangent cone z(x+z); the contact-eleven branch is z=0.
            slope = sp.Rational(0)
        else:
            linear_x = sp.Poly(q_transform, x, z).coeff_monomial(x)
            linear_z = sp.Poly(q_transform, x, z).coeff_monomial(z)
            if linear_z == 0:
                raise ArithmeticError("the selected Q branch became vertical")
            slope = sp.Rational(-linear_x, linear_z)
        slopes.append(slope)
        replay.append({
            "centre": stage,
            "Q_multiplicity": q_multiplicity,
            "S_multiplicity": s_multiplicity,
            "selected_Q_tangent_slope": str(slope),
        })

        q_numerator = sp.expand(q_transform.subs(z, x * (u + slope)))
        s_numerator = sp.expand(s_transform.subs(z, x * (u + slope)))
        q_transform = sp.expand(sp.cancel(
            q_numerator / x**q_multiplicities[stage - 1]
        ).subs(u, z))
        s_transform = sp.expand(sp.cancel(
            s_numerator / x**s_multiplicities[stage - 1]
        ).subs(u, z))
        if not (sp.Poly(q_transform, x, z).is_multivariate
                and sp.Poly(s_transform, x, z).is_multivariate):
            raise ArithmeticError("a strict transform stopped being polynomial")

    q_exceptional = sp.Poly(q_transform.subs(x, 0), z)
    s_exceptional = sp.Poly(s_transform.subs(x, 0), z)
    if q_exceptional.as_expr() != z or s_exceptional.as_expr() != z + 1:
        raise ArithmeticError(
            f"the p16 separation changed: Q={q_exceptional}, S={s_exceptional}"
        )
    s_p16_coordinate = sp.Rational(-1)

    # After r blow-ups, the old coordinate is
    # z=sum_(j<=r) slope_j*x^j+x^r*z_r.  At p16 the last local coordinate
    # equals -1, so it modifies the coefficient of x^15.
    q_coefficients = {
        index: slope for index, slope in enumerate(slopes, start=1) if slope
    }
    s_coefficients = dict(q_coefficients)
    s_coefficients[15] = s_coefficients.get(15, sp.Rational(0)) + s_p16_coordinate
    s_coefficients = {
        exponent: value for exponent, value in s_coefficients.items() if value
    }
    expected = {
        0: (
            {3: 1, 9: 1, 11: -1, 13: 1, 15: 2},
            {3: 1, 9: 1, 11: -1, 13: 1, 15: 1},
        ),
        1: (
            {3: 1, 7: 1, 9: 1, 11: 1, 13: 6, 15: 1},
            {3: 1, 7: 1, 9: 1, 11: 1, 13: 6},
        ),
    }
    expected_q, expected_s = expected[parameter]
    if q_coefficients != expected_q or s_coefficients != expected_s:
        raise ArithmeticError("the exact selected-branch coefficients changed")

    return {
        "carrier_parameter": parameter,
        "chart": "Y=1, x=X/Y, z=Z/Y",
        "derivation": (
            "fifteen exact strict transforms of Q and S; no numerical "
            "Taylor approximation"
        ),
        "Q_selected_branch_coefficients_through_x15": {
            str(key): str(value) for key, value in sorted(q_coefficients.items())
        },
        "S_selected_branch_coefficients_through_x15": {
            str(key): str(value) for key, value in sorted(s_coefficients.items())
        },
        "post_p15_exceptional_restrictions": {"Q": "z", "S": "z+1"},
        "p16_S_coordinate_in_Q_chart": "-1",
        "blowup_replay": replay,
    }


def _a0_cluster_jet_chart() -> dict[str, object]:
    """Compatibility wrapper used by the original replay tests."""

    return _cluster_jet_chart(0)


def _selected_s_coefficients(parameter: int = 0) -> dict[int, int]:
    record = _cluster_jet_chart(parameter)
    return {
        int(key): int(value)
        for key, value in record[
            "S_selected_branch_coefficients_through_x15"
        ].items()
    }


def _affine_monomials(degree: int) -> tuple[tuple[int, int], ...]:
    # x^i z^j represents X^i Z^j Y^(degree-i-j).
    return tuple(
        (i, j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    )


def _primitive_row(values: Sequence[object]) -> tuple[int, ...]:
    row = _primitive_integer_vector(values)
    first = next(value for value in row if value)
    return row if first > 0 else tuple(-value for value in row)


def _jet_condition_matrix(
    degree: int,
    multiplicities: Sequence[int],
    parameter: int = 0,
) -> tuple[
    tuple[tuple[int, int], ...],
    list[list[int]],
    list[list[object]],
]:
    """Build the exact complete-ideal matrix for a curvilinear cluster.

    At the k-th centre put ``A_(k-1)=m_1+...+m_(k-1)`` and

        phi_(k-1)(x)=sum_(r<k) c_r x^r.

    A degree-d form is a section of ``dH-sum m_i E_i`` precisely when

        F(x,phi_(k-1)(x)+x^(k-1)u) / x^A_(k-1)

    is regular and has multiplicity at least ``m_k`` for every k.  The rows
    below are exactly the forbidden Laurent/low-total-degree coefficients.
    This automatically allows exceptional components; replacing it by the
    stronger coordinatewise plane-multiplicity test would be wrong.
    """

    x = sp.symbols("x")
    coefficients = _selected_s_coefficients(parameter)
    monomials = _affine_monomials(degree)
    condition_rows: dict[tuple[int, int, int], dict[int, int]] = {}
    assigned_before = 0

    for centre, required in enumerate(multiplicities, start=1):
        phi = sum(
            value * x**exponent
            for exponent, value in coefficients.items()
            if exponent < centre
        )
        phi_powers = [sp.Poly(sp.expand(phi**power), x)
                      for power in range(degree + 1)]
        for column, (x_power, z_power) in enumerate(monomials):
            for u_power in range(z_power + 1):
                scalar = comb(z_power, u_power)
                for (phi_power,), coefficient in phi_powers[
                    z_power - u_power
                ].terms():
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
                        row = condition_rows.setdefault(key, {})
                        row[column] = row.get(column, 0) + scalar * int(coefficient)
        assigned_before += required

    dense_rows: list[list[int]] = []
    sparse_rows: list[list[object]] = []
    for key in sorted(condition_rows):
        entries = sorted(
            (column, value)
            for column, value in condition_rows[key].items()
            if value
        )
        if not entries:
            continue
        dense = [0] * len(monomials)
        for column, value in entries:
            dense[column] = value
        dense_rows.append(dense)
        sparse_rows.append([*key, [[column, value] for column, value in entries]])
    return monomials, dense_rows, sparse_rows


def _coefficient_vector(
    polynomial: sp.Expr,
    monomials: Sequence[tuple[int, int]],
    x: sp.Symbol,
    z: sp.Symbol,
) -> tuple[int, ...]:
    poly = sp.Poly(polynomial, x, z)
    return tuple(
        int(poly.coeff_monomial(x**i * z**j)) for i, j in monomials
    )


def _polynomial_from_vector(
    values: Sequence[int],
    monomials: Sequence[tuple[int, int]],
    x: sp.Symbol,
    z: sp.Symbol,
) -> sp.Expr:
    return sp.expand(sum(
        value * x**i * z**j
        for value, (i, j) in zip(values, monomials)
    ))


def _actual_cluster_multiplicities(
    values: Sequence[int],
    monomials: Sequence[tuple[int, int]],
    centre_count: int,
    parameter: int = 0,
) -> tuple[int, ...]:
    """Recover the strict plane curve's actual multiplicity sequence."""

    x = sp.symbols("x")
    coefficients = _selected_s_coefficients(parameter)
    maximum_z_power = max(j for _i, j in monomials)
    assigned_before = 0
    answer: list[int] = []
    for centre in range(1, centre_count + 1):
        phi = sum(
            value * x**exponent
            for exponent, value in coefficients.items()
            if exponent < centre
        )
        phi_powers = [sp.Poly(sp.expand(phi**power), x)
                      for power in range(maximum_z_power + 1)]
        terms: dict[tuple[int, int], int] = {}
        for scalar, (x_power, z_power) in zip(values, monomials):
            if not scalar:
                continue
            for u_power in range(z_power + 1):
                binomial = comb(z_power, u_power)
                for (phi_power,), coefficient in phi_powers[
                    z_power - u_power
                ].terms():
                    key = (
                        x_power + (centre - 1) * u_power + phi_power,
                        u_power,
                    )
                    terms[key] = (
                        terms.get(key, 0) + scalar * binomial * int(coefficient)
                    )
        nonzero_terms = {key: value for key, value in terms.items() if value}
        if not nonzero_terms:
            raise ArithmeticError("the candidate plane section vanished")
        actual = min(
            x_power - assigned_before + u_power
            for x_power, u_power in nonzero_terms
        )
        if actual < 0:
            raise ArithmeticError("an actual strict transform is not regular")
        answer.append(actual)
        assigned_before += actual
    return tuple(answer)


def _jet_case(
    n: int, cusp: Sequence[int]
) -> tuple[str, int, tuple[int, ...], dict[str, int], int] | None:
    vector = tuple(cusp)
    genus_five = (11, *((-3,) * 13), -2, *((0,) * (n - 14)))
    if n in (14, 15, 16) and vector == genus_five:
        return "genus_5", 11, (*((3,) * 13), 2, *((0,) * (n - 14))), {
            "C": 1,
            "Q": 2,
        }, 5
    if n == 16 and vector == (4, *((-1,) * 16)):
        return "genus_3_rigid", 4, (1,) * 16, {"Q": 1}, 3
    if n == 16 and vector == (
        21, *((-7,) * 5), -6, *((-4,) * 10)
    ):
        return "genus_10", 21, (*((7,) * 5), 6, *((4,) * 10)), {
            "L": 1,
            "Q": 2,
            "S": 2,
        }, 10
    if n == 16 and vector == (
        23, *((-7,) * 8), -5, *((-4,) * 7)
    ):
        return "genus_11", 23, (*((7,) * 8), 5, *((4,) * 7)), {
            "C": 1,
            "Q": 2,
            "S": 2,
        }, 11
    return None


@lru_cache(maxsize=None)
def _plane_jet_linear_system_record(
    n: int, cusp: tuple[int, ...], parameter: int = 0
) -> dict[str, object]:
    case = _jet_case(n, cusp)
    if case is None:
        raise ValueError("not one of the six a=0 single-component jet cases")
    name, degree, multiplicities, wall_factors, genus = case
    if parameter not in (0, 1):
        raise ValueError("the frozen explicit plane-jet systems use a=0 or a=1")
    x, z, cubic, quartic, sextic = _plane_forms(parameter)
    forms = {"L": z, "C": cubic, "Q": quartic, "S": sextic}
    wall_polynomial = sp.expand(sp.prod(
        forms[label] ** exponent for label, exponent in wall_factors.items()
    ))

    monomials, dense_rows, sparse_rows = _jet_condition_matrix(
        degree, multiplicities, parameter
    )
    matrix = sp.polys.matrices.DomainMatrix.from_list(dense_rows, sp.ZZ).to_field()
    rank = matrix.rank()
    nullity = len(monomials) - rank
    nullspace = matrix.nullspace().to_Matrix()
    if nullspace.rows != nullity:
        raise ArithmeticError("the exact jet kernel dimension changed")

    # Row-reduce the kernel itself so the stored sparse plane section is
    # canonical and independent of DomainMatrix's denominator clearing.
    canonical_kernel, _pivots = nullspace.rref()
    kernel_basis = tuple(
        _primitive_row(tuple(canonical_kernel.row(row)))
        for row in range(canonical_kernel.rows)
    )
    wall_vector = _coefficient_vector(wall_polynomial, monomials, x, z)
    if any(
        sum(row[column] * wall_vector[column]
            for column in range(len(monomials)))
        for row in dense_rows
    ):
        raise ArithmeticError("the affine wall cycle left the jet kernel")

    condition_payload_sha256 = sha256(json.dumps(
        sparse_rows, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    common = {
        "case": name,
        "carrier_parameter": f"a={parameter}",
        "degree": degree,
        "assigned_multiplicities": list(multiplicities),
        "homogeneous_monomial_count": len(monomials),
        "forbidden_jet_condition_count": len(dense_rows),
        "exact_condition_matrix_rank": rank,
        "complete_linear_system_h0": nullity,
        "complete_projective_dimension": nullity - 1,
        "condition_matrix_sparse_payload_sha256": condition_payload_sha256,
        "wall_cycle_plane_factorization": "*".join(
            label if exponent == 1 else f"{label}^{exponent}"
            for label, exponent in wall_factors.items()
        ),
        "wall_cycle_in_exact_kernel": True,
        "arithmetic_genus": genus,
        "jet_rule": (
            "vanishing of every negative or total-degree<m_k coefficient in "
            "F(x,phi_(k-1)+x^(k-1)u)/x^(m_1+...+m_(k-1))"
        ),
    }

    if nullity == 1:
        if _primitive_row(wall_vector) != kernel_basis[0]:
            raise ArithmeticError("the rigid kernel is not the affine wall cycle")
        return {
            **common,
            "primitive_linear_system_rigid": True,
            "unique_section_is_affine_wall_cycle": True,
            "semiample_certified": False,
            "semiampleness_boundary": (
                "h0(F)=1 proves that the primitive class does not move; "
                "this primitive calculation alone does not decide movement "
                "of a multiple mF, which is audited separately"
            ),
        }
    if nullity != 2:
        raise ArithmeticError(f"unexpected single-component cusp h0={nullity}")

    candidates: list[tuple[tuple[object, ...], tuple[int, ...], sp.Expr]] = []
    for vector in kernel_basis:
        polynomial = _polynomial_from_vector(vector, monomials, x, z)
        polynomial_gcd = sp.gcd(
            sp.Poly(polynomial, x, z), sp.Poly(wall_polynomial, x, z)
        )
        if polynomial_gcd.total_degree() != 0:
            continue
        actual = _actual_cluster_multiplicities(vector, monomials, n, parameter)
        if actual != multiplicities:
            continue
        key = (
            sum(value != 0 for value in vector),
            max(abs(value) for value in vector),
            vector,
        )
        candidates.append((key, vector, polynomial))
    if not candidates:
        raise ArithmeticError("no coprime exact plane member was found")
    _key, second_vector, second_polynomial = min(candidates, key=lambda item: item[0])
    second_terms = [
        [i, j, int(coefficient)]
        for coefficient, (i, j) in zip(second_vector, monomials)
        if coefficient
    ]
    second_digest = sha256(json.dumps(
        second_terms, separators=(",", ":")
    ).encode()).hexdigest()
    actual_multiplicities = _actual_cluster_multiplicities(
        second_vector, monomials, n, parameter
    )
    if actual_multiplicities != multiplicities:
        raise ArithmeticError("the second plane member gained an exceptional component")
    if sp.gcd(
        sp.Poly(second_polynomial, x, z), sp.Poly(wall_polynomial, x, z)
    ).total_degree() != 0:
        raise ArithmeticError("the two pencil members acquired a plane component")

    return {
        **common,
        "primitive_linear_system_rigid": False,
        "second_plane_section_sparse_terms": second_terms,
        "sparse_term_convention": (
            "[i,j,c] means c*X^i*Z^j*Y^(degree-i-j)"
        ),
        "second_plane_section_sparse_sha256": second_digest,
        "second_plane_section_term_count": len(second_terms),
        "second_section_actual_multiplicities": list(actual_multiplicities),
        "second_section_has_no_exceptional_component": True,
        "second_section_coprime_to_wall_plane_product": True,
        "two_effective_members_have_no_common_curve": True,
        "base_point_free_pencil": True,
        "semiample_certified": True,
        "semiample_reason": (
            "the affine wall cycle and the exact coprime plane-jet section are "
            "distinct effective divisors in the primitive square-zero class; "
            "they have no common component, hence intersection zero makes them "
            "disjoint and their complete pencil base-point-free"
        ),
    }


@lru_cache(maxsize=2)
def _genus_three_multiple_audit(parameter: int = 0) -> dict[str, object]:
    """Find the exact first moving multiple of the primitive genus-three cusp."""

    if parameter not in (0, 1):
        raise ValueError("the frozen multiple audit uses a=0 or a=1")
    x, z, _cubic, quartic, _sextic = _plane_forms(parameter)
    records: list[dict[str, object]] = []
    for multiple in range(1, 5):
        degree = 4 * multiple
        multiplicities = (multiple,) * 16
        monomials, dense_rows, sparse_rows = _jet_condition_matrix(
            degree, multiplicities, parameter
        )
        matrix = sp.polys.matrices.DomainMatrix.from_list(
            dense_rows, sp.ZZ
        ).to_field()
        rank = matrix.rank()
        nullity = len(monomials) - rank
        nullspace = matrix.nullspace().to_Matrix()
        canonical_kernel, _pivots = nullspace.rref()
        kernel_basis = tuple(
            _primitive_row(tuple(canonical_kernel.row(row)))
            for row in range(canonical_kernel.rows)
        )
        wall_polynomial = sp.expand(quartic**multiple)
        wall_vector = _coefficient_vector(
            wall_polynomial, monomials, x, z
        )
        if any(
            sum(row[column] * wall_vector[column]
                for column in range(len(monomials)))
            for row in dense_rows
        ):
            raise ArithmeticError("a Q^m wall member left the multiple kernel")
        record: dict[str, object] = {
            "multiple": multiple,
            "divisor_class": f"{multiple}F",
            "degree": degree,
            "assigned_multiplicities": list(multiplicities),
            "homogeneous_monomial_count": len(monomials),
            "forbidden_jet_condition_count": len(dense_rows),
            "exact_condition_matrix_rank": rank,
            "complete_linear_system_h0": nullity,
            "complete_projective_dimension": nullity - 1,
            "condition_matrix_sparse_payload_sha256": sha256(json.dumps(
                sparse_rows, sort_keys=True, separators=(",", ":")
            ).encode()).hexdigest(),
            "wall_member": f"Q^{multiple}",
            "wall_member_in_exact_kernel": True,
        }
        if multiple < 4:
            if nullity != 1 or _primitive_row(wall_vector) != kernel_basis[0]:
                raise ArithmeticError(
                    f"the unique Q^{multiple} kernel statement changed"
                )
            record.update({
                "unique_section_is_wall_member": True,
                "moves": False,
            })
        else:
            if nullity != 2:
                raise ArithmeticError("the fourth multiple stopped being a pencil")
            candidates: list[
                tuple[tuple[object, ...], tuple[int, ...], sp.Expr]
            ] = []
            for vector in kernel_basis:
                polynomial = _polynomial_from_vector(
                    vector, monomials, x, z
                )
                if sp.gcd(
                    sp.Poly(polynomial, x, z),
                    sp.Poly(wall_polynomial, x, z),
                ).total_degree() != 0:
                    continue
                actual = _actual_cluster_multiplicities(
                    vector, monomials, 16, parameter
                )
                if actual != multiplicities:
                    continue
                key = (
                    sum(value != 0 for value in vector),
                    max(abs(value) for value in vector),
                    vector,
                )
                candidates.append((key, vector, polynomial))
            if not candidates:
                raise ArithmeticError("the fourth multiple lost its coprime fibre")
            _key, second_vector, second_polynomial = min(
                candidates, key=lambda item: item[0]
            )
            second_terms = [
                [i, j, int(coefficient)]
                for coefficient, (i, j) in zip(second_vector, monomials)
                if coefficient
            ]
            actual_multiplicities = _actual_cluster_multiplicities(
                second_vector, monomials, 16, parameter
            )
            if actual_multiplicities != multiplicities:
                raise ArithmeticError("the fourth-multiple fibre gained exceptionals")
            if sp.gcd(
                sp.Poly(second_polynomial, x, z),
                sp.Poly(wall_polynomial, x, z),
            ).total_degree() != 0:
                raise ArithmeticError("the fourth-multiple fibres share a component")
            record.update({
                "moves": True,
                "second_plane_section_sparse_terms": second_terms,
                "sparse_term_convention": (
                    "[i,j,c] means c*X^i*Z^j*Y^(degree-i-j)"
                ),
                "second_plane_section_sparse_sha256": sha256(json.dumps(
                    second_terms, separators=(",", ":")
                ).encode()).hexdigest(),
                "second_plane_section_term_count": len(second_terms),
                "second_section_actual_multiplicities": list(
                    actual_multiplicities
                ),
                "second_section_has_no_exceptional_component": True,
                "second_section_coprime_to_Q4": True,
                "two_effective_members_have_no_common_curve": True,
                "base_point_free_pencil": True,
                "fibre_arithmetic_genus": 9,
            })
        records.append(record)

    expected = [(15, 14, 1), (45, 44, 1), (91, 90, 1), (153, 151, 2)]
    actual = [
        (
            record["homogeneous_monomial_count"],
            record["exact_condition_matrix_rank"],
            record["complete_linear_system_h0"],
        )
        for record in records
    ]
    if actual != expected:
        raise ArithmeticError(f"the semiample-index audit changed: {actual}")
    return {
        "carrier_parameter": f"a={parameter}",
        "primitive_cusp_class": [4, *((-1,) * 16)],
        "primitive_arithmetic_genus": 3,
        "semiample_index": 4,
        "index_is_exact_because_h0_mF_is_one_for_m_1_2_3": True,
        "fourth_multiple_base_point_free": True,
        "fourth_multiple_fibre_arithmetic_genus": 9,
        "terminology": (
            "Halphen-type semiample index 4; not an assertion that this is a "
            "classical anticanonical Halphen surface"
        ),
        "multiples": records,
    }


def _file_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _payload_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_payload_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _nef_rows(roots: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    # -<D,r> >= 0 is precisely D.C_r >= 0 on the surface.
    return tuple((root[0], *(-entry for entry in root[1:])) for root in roots)


def _connected_components(
    active: Sequence[int], roots: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    unseen = set(active)
    components: list[tuple[int, ...]] = []
    while unseen:
        stack = [unseen.pop()]
        component: list[int] = []
        while stack:
            left = stack.pop()
            component.append(left)
            neighbours = [
                right
                for right in tuple(unseen)
                if lorentz_pair(roots[left], roots[right]) != 0
            ]
            for right in neighbours:
                unseen.remove(right)
                stack.append(right)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def _positive_primitive(values: Iterable[object]) -> tuple[int, ...]:
    result = _primitive_integer_vector(tuple(values))
    if all(value < 0 for value in result):
        result = tuple(-value for value in result)
    if not all(value > 0 for value in result):
        raise ArithmeticError("an affine null vector is not strictly positive")
    common = 0
    for value in result:
        common = gcd(common, value)
    if common != 1:
        raise ArithmeticError("an affine null vector is not primitive")
    return result


def _component_record(
    component: Sequence[int],
    roots: Sequence[Sequence[int]],
    cusp: Sequence[int],
) -> dict[str, object]:
    gram = sp.Matrix(
        [
            [lorentz_pair(roots[left], roots[right]) for right in component]
            for left in component
        ]
    )
    nullspace = gram.nullspace()
    if len(nullspace) != 1:
        raise ArithmeticError("an active component is not affine of corank one")
    weights = _positive_primitive(nullspace[0])
    cycle = tuple(
        sum(weights[position] * roots[label][coordinate]
            for position, label in enumerate(component))
        for coordinate in range(len(cusp))
    )
    if cycle[0] < 0:
        cycle = tuple(-entry for entry in cycle)
        weights = tuple(-entry for entry in weights)
    if cycle != tuple(cusp):
        raise ArithmeticError("a primitive affine wall cycle does not equal its cusp")
    inertia = exact_inertia(gram)
    if inertia != (len(component) - 1, 0, 1):
        raise ArithmeticError(f"active component inertia changed: {inertia}")
    return {
        "facet_labels_one_based": [label + 1 for label in component],
        "affine_rank": len(component) - 1,
        "gram_inertia": list(inertia),
        "primitive_null_weights": list(weights),
        "wall_cycle_vector": list(cycle),
        "wall_cycle_equals_primitive_cusp": True,
    }


def _semiample_reason(
    n: int,
    cusp: Sequence[int],
    component_count: int,
) -> str | None:
    if component_count >= 2:
        return (
            "two connected affine components give distinct disjoint effective "
            "divisors, each linearly equivalent to the primitive cusp; their "
            "pencil is base-point-free"
        )
    line_pencil = (1, -1, *((0,) * (n - 1)))
    if tuple(cusp) == line_pencil:
        return (
            "pullback of the base-point-free pencil of plane lines through the "
            "first blow-up centre"
        )
    q13 = (4, -2, *((-1,) * 12))
    if n == 13 and tuple(cusp) == q13:
        return (
            "the absolutely irreducible quartic Q and the affine wall cycle are "
            "distinct effective divisors in the same primitive square-zero class"
        )
    return None


def _dimension_record(n: int) -> dict[str, object]:
    roots = enumerate_simple_roots(n)
    rays = _extreme_rays_via_cdd(_nef_rows(roots))
    if rays is None:
        raise RuntimeError("cddexec_gmp is required for the independent cusp census")
    expected_total = int(VINBERG_VOLUME_DATA[n]["finite_vertices"]) + int(
        VINBERG_VOLUME_DATA[n]["ideal_vertices"]
    )
    if len(rays) != expected_total:
        raise ArithmeticError(f"P{n} nef extreme-ray count changed")

    cusps: list[dict[str, object]] = []
    positive_square_count = 0
    for raw in rays:
        ray = tuple(raw if raw[0] > 0 else tuple(-entry for entry in raw))
        surface_square = -lorentz_pair(ray, ray)
        if surface_square > 0:
            positive_square_count += 1
            continue
        if surface_square < 0:
            raise ArithmeticError("a nef extreme ray has negative surface square")

        active = tuple(
            label for label, root in enumerate(roots)
            if lorentz_pair(ray, root) == 0
        )
        components = _connected_components(active, roots)
        component_records = [
            _component_record(component, roots, ray) for component in components
        ]
        if len(active) != n + len(components) - 1:
            raise ArithmeticError("Vinberg full zero-neighbourhood identity failed")
        if sum(record["affine_rank"] for record in component_records) != n - 1:
            raise ArithmeticError("total parabolic rank is not n-1")
        for left_index, left in enumerate(components):
            for right in components[left_index + 1 :]:
                if any(
                    lorentz_pair(roots[i], roots[j]) != 0
                    for i in left for j in right
                ):
                    raise ArithmeticError("distinct affine components are not orthogonal")

        k_pairing = canonical_pairing(ray)
        if k_pairing % 2:
            raise ArithmeticError("isotropic adjunction parity failed")
        genus = 1 + k_pairing // 2
        reason = _semiample_reason(n, ray, len(components))
        primitive_linear_system: dict[str, object] | None = None
        semiample_multiple_system: dict[str, object] | None = None
        if reason is None and len(components) == 1:
            primitive_linear_system = _plane_jet_linear_system_record(n, ray)
            if primitive_linear_system["semiample_certified"]:
                reason = str(primitive_linear_system["semiample_reason"])
            elif primitive_linear_system["case"] == "genus_3_rigid":
                semiample_multiple_system = _genus_three_multiple_audit(0)
                if semiample_multiple_system["fourth_multiple_base_point_free"]:
                    reason = (
                        "the primitive system is rigid, but the exact first "
                        "moving multiple is |4F|; its Q^4 wall member and a "
                        "coprime degree-16 plane-jet member are disjoint, so "
                        "|4F| is base-point-free"
                    )
        cusps.append({
            "primitive_nef_vector": list(ray),
            "surface_square": 0,
            "canonical_pairing": k_pairing,
            "arithmetic_genus": genus,
            "active_facet_count": len(active),
            "affine_component_count": len(components),
            "full_zero_neighbourhood_identity": f"{len(active)}={n}+{len(components)}-1",
            "components": component_records,
            "primitive_linear_system": primitive_linear_system,
            "semiample_multiple_system": semiample_multiple_system,
            "semiample_certified": reason is not None,
            "semiample_reason": reason,
        })

    cusps.sort(key=lambda row: tuple(row["primitive_nef_vector"]))
    expected_ideal = int(VINBERG_VOLUME_DATA[n]["ideal_vertices"])
    if len(cusps) != expected_ideal or positive_square_count != int(
        VINBERG_VOLUME_DATA[n]["finite_vertices"]
    ):
        raise ArithmeticError(f"P{n} finite/ideal split changed")
    return {
        "hyperbolic_dimension": n,
        "facet_count": len(roots),
        "nef_extreme_ray_count": len(rays),
        "positive_square_ray_count": positive_square_count,
        "isotropic_cusp_count": len(cusps),
        "semiample_cusp_count": sum(row["semiample_certified"] for row in cusps),
        "cusp_genus_histogram": {
            str(key): value
            for key, value in sorted(Counter(row["arithmetic_genus"] for row in cusps).items())
        },
        "cusps": cusps,
    }


def _a1_nonisomorphic_carrier_audit() -> dict[str, object]:
    """Replay the six primitive systems on the non-isomorphic a=1 carrier."""

    cases = [
        (16, (23, *((-7,) * 8), -5, *((-4,) * 7))),
        (16, (21, *((-7,) * 5), -6, *((-4,) * 10))),
        (14, (11, *((-3,) * 13), -2)),
        (15, (11, *((-3,) * 13), -2, 0)),
        (16, (11, *((-3,) * 13), -2, 0, 0)),
        (16, (4, *((-1,) * 16))),
    ]

    # The exact matrices are independent. Replay them in isolated Python
    # processes before building the a=0 census, so the two 300-column exact
    # nullspace calculations never coexist in one GMP/FLINT heap. Each child
    # returns a freshly generated record, not a trusted precomputed answer.
    def replay_case(item: tuple[int, tuple[int, ...]]) -> dict[str, object]:
        dimension, vector = item
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "scripts.certify_rational_i1n_cusp_fibrations",
                "--emit-plane-jet-record",
                "--dimension",
                str(dimension),
                "--parameter",
                "1",
                "--cusp-vector",
                ",".join(str(value) for value in vector),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        a1_record = json.loads(completed.stdout)
        return {
            "hyperbolic_dimension": dimension,
            "primitive_nef_vector": list(vector),
            **a1_record,
        }

    def replay_multiple() -> dict[str, object]:
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "scripts.certify_rational_i1n_cusp_fibrations",
                "--emit-genus-three-multiple-audit",
                "--parameter",
                "1",
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        return json.loads(completed.stdout)

    # Deliberately serial: every heavy GMP/FLINT heap belongs to one fresh
    # worker and is released before the next exact nullspace begins.
    records = [replay_case(item) for item in cases]
    multiple_audit = replay_multiple()

    moving = [record for record in records if record["semiample_certified"]]
    rigid = [record for record in records if not record["semiample_certified"]]
    if not (len(records) == 6 and len(moving) == 5 and len(rigid) == 1):
        raise ArithmeticError("the a=1 five-moving/one-rigid split changed")
    a1_fourth = multiple_audit["multiples"][3]
    if not (
        multiple_audit["semiample_index"] == 4
        and a1_fourth["base_point_free_pencil"]
    ):
        raise ArithmeticError("the a=1 index-four audit changed")
    return {
        "carrier": (
            "the explicit a=1 surface with unique genus-one Mori facet "
            "j=6912/31"
        ),
        "nonisomorphic_to_a0_j0_carrier": True,
        "same_complete_P14_P15_P16_nef_chambers": True,
        "audited_primitive_system_count": len(records),
        "base_point_free_pencil_count": len(moving),
        "primitive_rigid_count": len(rigid),
        "semiample_cusp_count": 28,
        "primitive_base_point_free_pencil_count": 27,
        "index_four_base_point_free_pencil_count": 1,
        "genus_three_multiple_audit": multiple_audit,
        "records": records,
    }


def _validate_a1_worker_certificate(document: Mapping[str, object]) -> None:
    if document.get("schema") != A1_WORKER_SCHEMA:
        raise ArithmeticError("wrong a=1 worker-certificate schema")
    if document.get("certificate_payload_sha256") != _payload_digest(document):
        raise ArithmeticError("a=1 worker-certificate payload changed")
    audit = document.get("audit")
    if not isinstance(audit, Mapping):
        raise ArithmeticError("the a=1 worker audit is malformed")
    if (
        audit.get("audited_primitive_system_count"),
        audit.get("base_point_free_pencil_count"),
        audit.get("primitive_rigid_count"),
        audit.get("semiample_cusp_count"),
        audit.get("index_four_base_point_free_pencil_count"),
    ) != (6, 5, 1, 28, 1):
        raise ArithmeticError("the a=1 worker summary changed")
    multiple = audit.get("genus_three_multiple_audit")
    if not isinstance(multiple, Mapping) or not (
        multiple.get("semiample_index") == 4
        and multiple["multiples"][3][
            "second_plane_section_sparse_sha256"
        ] == "f8b43b63e6e0c560303c398d983884f737623cc732b61ead6c2d604b18866b0b"
    ):
        raise ArithmeticError("the a=1 worker index-four audit changed")
    provenance = document.get("source_provenance")
    if not isinstance(provenance, Mapping):
        raise ArithmeticError("the a=1 worker provenance is missing")
    expected_sources = {
        "certifier": CERTIFIER,
        "nonisotrivial_family_script": FAMILY_SCRIPT,
        "nonisotrivial_family_certificate": FAMILY_JSON,
    }
    for label, path in expected_sources.items():
        record = provenance.get(label)
        if not isinstance(record, Mapping) or record != {
            "path": str(path.relative_to(ROOT)),
            "sha256": _file_digest(path),
        }:
            raise ArithmeticError(f"the a=1 worker source {label} is stale")


def _build_a1_worker_certificate() -> dict[str, object]:
    audit = _a1_nonisomorphic_carrier_audit()
    document: dict[str, object] = {
        "schema": A1_WORKER_SCHEMA,
        "claim": (
            "Six primitive cusp systems and the first four multiples of the "
            "rigid genus-three system are regenerated in serial fresh "
            "processes on the explicit a=1 carrier."
        ),
        "audit": audit,
        "source_provenance": {
            "certifier": {
                "path": str(CERTIFIER.relative_to(ROOT)),
                "sha256": _file_digest(CERTIFIER),
            },
            "nonisotrivial_family_script": {
                "path": str(FAMILY_SCRIPT.relative_to(ROOT)),
                "sha256": _file_digest(FAMILY_SCRIPT),
            },
            "nonisotrivial_family_certificate": {
                "path": str(FAMILY_JSON.relative_to(ROOT)),
                "sha256": _file_digest(FAMILY_JSON),
            },
        },
    }
    document["certificate_payload_sha256"] = _payload_digest(document)
    _validate_a1_worker_certificate(document)
    return document


def _load_a1_worker_certificate() -> dict[str, object]:
    document = json.loads(A1_WORKER_OUTPUT.read_text(encoding="utf-8"))
    _validate_a1_worker_certificate(document)
    return document


def _attach_a1_a0_comparison(
    a1_audit: dict[str, object],
    dimensions: Mapping[str, Mapping[str, object]],
) -> None:
    """Compare the independently generated a=1 and a=0 records exactly."""

    a0_records = {
        (int(dimension), tuple(cusp["primitive_nef_vector"])):
            cusp["primitive_linear_system"]
        for dimension, row in dimensions.items()
        for cusp in row["cusps"]
        if cusp.get("primitive_linear_system") is not None
    }
    for record in a1_audit["records"]:
        key = (
            int(record["hyperbolic_dimension"]),
            tuple(record["primitive_nef_vector"]),
        )
        a0_record = a0_records.get(key)
        if a0_record is None:
            raise ArithmeticError("an a=1 primitive record has no a=0 mate")
        if record["complete_linear_system_h0"] != a0_record[
            "complete_linear_system_h0"
        ]:
            raise ArithmeticError("a=0 and a=1 jet nullities diverged")
        if record["semiample_certified"]:
            if record["second_plane_section_sparse_sha256"] == a0_record[
                "second_plane_section_sparse_sha256"
            ]:
                raise ArithmeticError(
                    "an a=1 second plane section unexpectedly equals a=0"
                )
            record["second_plane_section_differs_from_a0"] = True
        else:
            # The unique plane factor is Q_a. Its x*z^2 coefficient changes
            # with a even though both primitive kernels have h0=1.
            record["unique_Q_section_differs_from_a0"] = True

    a0_multiple = next(
        cusp["semiample_multiple_system"]
        for row in dimensions.values()
        for cusp in row["cusps"]
        if cusp.get("semiample_multiple_system") is not None
    )
    a1_fourth = a1_audit["genus_three_multiple_audit"]["multiples"][3]
    a0_fourth = a0_multiple["multiples"][3]
    if a1_fourth["second_plane_section_sparse_sha256"] == a0_fourth[
        "second_plane_section_sparse_sha256"
    ]:
        raise ArithmeticError("the a=1 index-four section equals a=0")
    a1_audit["all_six_nullities_equal_a0"] = True
    a1_audit["all_five_second_plane_sections_differ_from_a0"] = True
    a1_audit["index_four_second_plane_section_differs_from_a0"] = True


def build_certificate() -> dict[str, object]:
    a1_worker_document = _load_a1_worker_certificate()
    a1_audit = json.loads(json.dumps(a1_worker_document["audit"]))
    dimensions: dict[str, Mapping[str, object]] = {}
    for n in DIMENSIONS:
        dimensions[str(n)] = _dimension_record(n)
        sp.core.cache.clear_cache()
        gc.collect()
    _attach_a1_a0_comparison(a1_audit, dimensions)
    total_cusps = sum(row["isotropic_cusp_count"] for row in dimensions.values())
    semiample = sum(row["semiample_cusp_count"] for row in dimensions.values())
    document: dict[str, object] = {
        "schema": SCHEMA,
        "claim": (
            "Every ideal vertex of the complete P12--P16 rational-surface nef "
            "chambers is a primitive effective isotropic divisor reconstructed "
            "exactly from each connected affine wall component.  All twenty-"
            "eight cusp classes on each of the non-isomorphic explicit a=0 "
            "and a=1 carriers are certified semiample: twenty-seven move in "
            "their primitive class, while the remaining primitive system is "
            "rigid and has exact first moving multiple |4F|."
        ),
        "scope": {
            "complete_cusp_census": True,
            "finite_volume_input": "the already certified complete P12--P16 wall sets",
            "semiampleness_claimed_for_every_cusp": True,
            "remaining_single_component_cusps": total_cusps - semiample,
            "all_six_formerly_open_primitive_linear_systems_audited": True,
            "genus_three_semiample_index_four_audited_on_a0_and_a1": True,
            "new_abstract_hyperbolic_polytope": False,
        },
        "algebro_geometric_vinberg_criterion": (
            "For a complete rational-surface Coxeter nef cone, Vinberg's ideal "
            "endpoint condition is witnessed inside Pic(X) by a full-rank "
            "orthogonal union of affine configurations of negative curves; each "
            "primitive affine null cycle is the nef isotropic cusp class."
        ),
        "algebro_geometric_mobility_criterion": (
            "For either explicit curvilinear a=0 or a=1 cluster, "
            "H0(dH-sum(m_iE_i)) is the exact "
            "kernel of the forbidden Laurent/low-jet coefficient matrix after "
            "successive blow-up substitution.  A second kernel vector whose "
            "plane form is coprime to the affine wall product and has the exact "
            "assigned multiplicity sequence gives a disjoint fibre and hence a "
            "base-point-free cusp pencil."
        ),
        "explicit_cluster_jet_charts": {
            "a=0": _cluster_jet_chart(0),
            "a=1": _cluster_jet_chart(1),
        },
        "a1_nonisomorphic_carrier_primitive_jet_audit": a1_audit,
        "totals": {
            "isotropic_cusps": total_cusps,
            "semiample_cusps": semiample,
            "primitive_base_point_free_pencils": 27,
            "index_four_base_point_free_pencils": 1,
            "unresolved_single_component_cusps": total_cusps - semiample,
        },
        "dimensions": dimensions,
        "source_provenance": {
            "parent_script": {
                "path": str(PARENT_SCRIPT.relative_to(ROOT)),
                "sha256": _file_digest(PARENT_SCRIPT),
            },
            "parent_certificate": {
                "path": str(PARENT_JSON.relative_to(ROOT)),
                "sha256": _file_digest(PARENT_JSON),
            },
            "certifier": {
                "path": str(CERTIFIER.relative_to(ROOT)),
                "sha256": _file_digest(CERTIFIER),
            },
            "nonisotrivial_family_script": {
                "path": str(FAMILY_SCRIPT.relative_to(ROOT)),
                "sha256": _file_digest(FAMILY_SCRIPT),
            },
            "nonisotrivial_family_certificate": {
                "path": str(FAMILY_JSON.relative_to(ROOT)),
                "sha256": _file_digest(FAMILY_JSON),
            },
            "a1_worker_certificate": {
                "path": str(A1_WORKER_OUTPUT.relative_to(ROOT)),
                "sha256": _file_digest(A1_WORKER_OUTPUT),
                "payload_sha256": a1_worker_document[
                    "certificate_payload_sha256"
                ],
            },
        },
    }
    document["certificate_payload_sha256"] = _payload_digest(document)
    validate(document)
    return document


def validate(document: Mapping[str, object]) -> None:
    if document.get("schema") != SCHEMA:
        raise ArithmeticError("wrong cusp-fibration schema")
    if document.get("certificate_payload_sha256") != _payload_digest(document):
        raise ArithmeticError("cusp-fibration payload digest changed")
    worker_document = _load_a1_worker_certificate()
    provenance = document.get("source_provenance")
    if not isinstance(provenance, Mapping) or provenance.get(
        "a1_worker_certificate"
    ) != {
        "path": str(A1_WORKER_OUTPUT.relative_to(ROOT)),
        "sha256": _file_digest(A1_WORKER_OUTPUT),
        "payload_sha256": worker_document["certificate_payload_sha256"],
    }:
        raise ArithmeticError("the aggregate a=1 worker provenance is stale")
    totals = document.get("totals")
    if totals != {
        "isotropic_cusps": 28,
        "semiample_cusps": 28,
        "primitive_base_point_free_pencils": 27,
        "index_four_base_point_free_pencils": 1,
        "unresolved_single_component_cusps": 0,
    }:
        raise ArithmeticError(f"cusp-fibration totals changed: {totals}")
    dimensions = document.get("dimensions")
    if not isinstance(dimensions, Mapping):
        raise ArithmeticError("cusp-fibration dimensions are malformed")
    expected = {
        "12": (37, 2, 2),
        "13": (46, 3, 3),
        "14": (94, 5, 5),
        "15": (145, 6, 6),
        "16": (325, 12, 12),
    }
    for key, triple in expected.items():
        row = dimensions.get(key)
        if not isinstance(row, Mapping):
            raise ArithmeticError(f"missing P{key} cusp record")
        actual = (
            row.get("nef_extreme_ray_count"),
            row.get("isotropic_cusp_count"),
            row.get("semiample_cusp_count"),
        )
        if actual != triple:
            raise ArithmeticError(f"P{key} cusp record changed: {actual}")

    jet_records = [
        cusp["primitive_linear_system"]
        for row in dimensions.values()
        for cusp in row["cusps"]
        if cusp.get("primitive_linear_system") is not None
    ]
    if len(jet_records) != 6:
        raise ArithmeticError("the six-case primitive jet audit changed")
    if Counter(record["complete_linear_system_h0"] for record in jet_records) != {
        1: 1,
        2: 5,
    }:
        raise ArithmeticError("the primitive jet kernel dimensions changed")
    moving = [record for record in jet_records if record["semiample_certified"]]
    rigid = [record for record in jet_records if not record["semiample_certified"]]
    if not (
        len(moving) == 5
        and all(record["base_point_free_pencil"] for record in moving)
        and len(rigid) == 1
        and rigid[0]["case"] == "genus_3_rigid"
        and rigid[0]["unique_section_is_affine_wall_cycle"]
    ):
        raise ArithmeticError("the five-moving/one-rigid cusp partition changed")
    multiple_records = [
        cusp["semiample_multiple_system"]
        for row in dimensions.values()
        for cusp in row["cusps"]
        if cusp.get("semiample_multiple_system") is not None
    ]
    if len(multiple_records) != 1:
        raise ArithmeticError("the unique multiple-system audit changed")
    a0_multiple = multiple_records[0]
    if not (
        a0_multiple["semiample_index"] == 4
        and a0_multiple["fourth_multiple_base_point_free"] is True
        and [
            record["complete_linear_system_h0"]
            for record in a0_multiple["multiples"]
        ] == [1, 1, 1, 2]
        and a0_multiple["multiples"][3][
            "second_plane_section_sparse_sha256"
        ] == "1393b70fc23154a593c95781a8990caab364716c3cbcc308296f91883657ec12"
    ):
        raise ArithmeticError("the a=0 exact semiample-index-four audit changed")

    a1 = document.get("a1_nonisomorphic_carrier_primitive_jet_audit")
    if not isinstance(a1, Mapping):
        raise ArithmeticError("the a=1 primitive jet audit is missing")
    if (
        a1.get("audited_primitive_system_count"),
        a1.get("base_point_free_pencil_count"),
        a1.get("primitive_rigid_count"),
        a1.get("semiample_cusp_count"),
        a1.get("index_four_base_point_free_pencil_count"),
    ) != (6, 5, 1, 28, 1):
        raise ArithmeticError("the a=1 primitive jet summary changed")
    if not (
        a1.get("all_six_nullities_equal_a0") is True
        and a1.get("all_five_second_plane_sections_differ_from_a0") is True
        and a1.get("index_four_second_plane_section_differs_from_a0") is True
        and a1["genus_three_multiple_audit"]["semiample_index"] == 4
        and a1["genus_three_multiple_audit"]["multiples"][3][
            "second_plane_section_sparse_sha256"
        ] == "f8b43b63e6e0c560303c398d983884f737623cc732b61ead6c2d604b18866b0b"
        and all(
            record.get("second_plane_section_differs_from_a0") is True
            for record in a1["records"]
            if record["semiample_certified"]
        )
    ):
        raise ArithmeticError("the a=1/a=0 section comparison changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--rebuild-a1-workers", action="store_true")
    parser.add_argument("--check-a1-workers", action="store_true")
    parser.add_argument("--emit-plane-jet-record", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--emit-genus-three-multiple-audit", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--dimension", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--parameter", type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument("--cusp-vector", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.emit_plane_jet_record:
        if args.dimension is None or args.cusp_vector is None:
            parser.error("the internal jet replay needs --dimension and --cusp-vector")
        vector = tuple(int(value) for value in args.cusp_vector.split(","))
        print(json.dumps(
            _plane_jet_linear_system_record(
                args.dimension, vector, args.parameter
            ),
            sort_keys=True,
        ))
        return
    if args.emit_genus_three_multiple_audit:
        print(json.dumps(
            _genus_three_multiple_audit(args.parameter),
            sort_keys=True,
        ))
        return
    if args.rebuild_a1_workers:
        document = _build_a1_worker_certificate()
        A1_WORKER_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        A1_WORKER_OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({
            "a1_worker_certificate_payload_sha256": document[
                "certificate_payload_sha256"
            ],
            "path": str(A1_WORKER_OUTPUT.relative_to(ROOT)),
        }, indent=2, sort_keys=True))
        return
    if args.check_a1_workers:
        document = _load_a1_worker_certificate()
        print(json.dumps({
            "a1_worker_certificate_payload_sha256": document[
                "certificate_payload_sha256"
            ],
            "path": str(A1_WORKER_OUTPUT.relative_to(ROOT)),
        }, indent=2, sort_keys=True))
        return
    if args.check_archive:
        document = json.loads(args.output.read_text(encoding="utf-8"))
        validate(document)
    else:
        document = build_certificate()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({
        "certificate_payload_sha256": document["certificate_payload_sha256"],
        "totals": document["totals"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
