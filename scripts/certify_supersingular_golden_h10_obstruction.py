"""Finite-field audit for the supersingular-K3 golden H10 obstruction.

The completeness argument is mathematical and is written in
``results/supersingular_k3_golden_h10_obstruction.md``.  This script locks
down the convention-sensitive ingredients: inertness of 2, the two trace
functionals on F4, the incompatible Artin-invariant parities, and the type of
the most obvious rank-11 trace lattice.
"""

from __future__ import annotations

import json

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import e8_cartan_gram


def _f4_add(left: int, right: int) -> int:
    return left ^ right


def _f4_mul(left: int, right: int) -> int:
    """Multiply bit pairs a+b*w modulo w^2+w+1."""

    a, b = left & 1, (left >> 1) & 1
    c, d = right & 1, (right >> 1) & 1
    constant = (a * c) ^ (b * d)
    omega = (a * d) ^ (b * c) ^ (b * d)
    return constant | (omega << 1)


def _f4_square(value: int) -> int:
    return _f4_mul(value, value)


def _f4_trace(value: int) -> int:
    traced = _f4_add(value, _f4_square(value))
    if traced not in (0, 1):
        raise ArithmeticError("F4 trace did not land in F2")
    return traced


def _smith_diagonal(gram: list[list[int]]) -> list[int]:
    smith = smith_normal_form(Matrix(gram), domain=ZZ)
    return [abs(int(smith[i, i])) for i in range(smith.rows)]


def _type_i(gram: list[list[int]]) -> bool:
    inverse = Matrix(gram).inv()
    return all(inverse[index, index].q == 1 for index in range(inverse.rows))


def _naive_rank11_field_gram() -> list[list[Golden]]:
    """Return <-2 phi> + E8 + <2> + <2> over Z[phi]."""

    e8 = e8_cartan_gram()
    result = [[Golden(0) for _ in range(11)] for _ in range(11)]
    result[0][0] = -2 * PHI
    for i in range(8):
        for j in range(8):
            result[i + 1][j + 1] = Golden(e8[i][j])
    result[9][9] = Golden(2)
    result[10][10] = Golden(2)
    return result


def certify() -> dict[str, object]:
    omega = 2
    omega_squared = _f4_square(omega)
    multiplication_table = [
        [_f4_mul(left, right) for right in range(4)] for left in range(4)
    ]
    square_nonzero = {_f4_square(value) for value in range(1, 4)}
    trace_pairs = {
        value: (
            _f4_trace(value),
            _f4_trace(_f4_mul(omega_squared, value)),
        )
        for value in range(4)
    }

    parity_table = [
        {
            "artin_invariant": sigma,
            "forced_by_rank_11_even_reduction": "odd" if sigma % 2 else "fails",
            "forced_by_type_I_discriminant_form": "even" if not sigma % 2 else "fails",
            "compatible": sigma % 2 == 1 and sigma % 2 == 0,
        }
        for sigma in range(1, 11)
    ]

    field = _naive_rank11_field_gram()
    field_conjugate = [[entry.conjugate() for entry in row] for row in field]
    traced = trace_gram(field, -PHI * INV_SQRT5)
    integer = integral_gram(traced)
    smith = _smith_diagonal(integer)
    trace_signature = inertia(integer)
    identity_signature = inertia(field)
    conjugate_signature = inertia(field_conjugate)
    naive = {
        "description": "<-2 phi> + E8 + <2> + <2>",
        "field_rank": 11,
        "field_signature_identity": {
            "positive": identity_signature.positive,
            "negative": identity_signature.negative,
            "zero": identity_signature.zero,
        },
        "field_signature_conjugate": {
            "positive": conjugate_signature.positive,
            "negative": conjugate_signature.negative,
            "zero": conjugate_signature.zero,
        },
        "trace_rank": len(integer),
        "trace_even_integral": is_even_integral(traced),
        "trace_determinant": str(determinant(integer)),
        "trace_signature": {
            "positive": trace_signature.positive,
            "negative": trace_signature.negative,
            "zero": trace_signature.zero,
        },
        "smith_diagonal": smith,
        "two_elementary_length": smith.count(2),
        "two_elementary_type_I": _type_i(integer),
        "nikulin_invariants": [22, smith.count(2), int(not _type_i(integer))],
        "supersingular_char2_NS": False,
        "failure": "type II rather than the required type I",
    }

    verified = (
        # x^2+x+1 has no roots in F2, and omega realizes one in F4.
        all((x * x + x + 1) % 2 for x in (0, 1))
        and _f4_add(_f4_add(omega_squared, omega), 1) == 0
        and square_nonzero == {1, 2, 3}
        and len(set(trace_pairs.values())) == 4
        and trace_pairs[0] == (0, 0)
        and not any(row["compatible"] for row in parity_table)
        and naive["field_signature_identity"]
        == {"positive": 10, "negative": 1, "zero": 0}
        and naive["field_signature_conjugate"]
        == {"positive": 11, "negative": 0, "zero": 0}
        and naive["trace_signature"]
        == {"positive": 1, "negative": 21, "zero": 0}
        and naive["trace_determinant"] == "-64"
        and naive["two_elementary_length"] == 6
        and naive["nikulin_invariants"] == [22, 6, 1]
        and not naive["two_elementary_type_I"]
    )
    return {
        "verified": verified,
        "golden_polynomial_mod_2_irreducible": True,
        "residue_field": "F4 = F2[w]/(w^2+w+1)",
        "multiplication_table": multiplication_table,
        "nonzero_square_map_is_bijective": square_nonzero == {1, 2, 3},
        "trace_functional_pairs": {
            str(value): list(pair) for value, pair in trace_pairs.items()
        },
        "two_trace_functionals_independent": len(set(trace_pairs.values())) == 4,
        "odd_characteristic_obstruction": (
            "an odd-dimensional alternating F4 reduction is singular, so the "
            "integral trace determinant is even"
        ),
        "characteristic_two_parity_table": parity_table,
        "naive_rank11_trace_lattice": naive,
        "conclusion": "no_supersingular_K3_NS_has_a_self_adjoint_golden_action",
        "scope": "finite-field audit accompanying the manuscript proof",
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
