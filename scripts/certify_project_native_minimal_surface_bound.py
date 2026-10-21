#!/usr/bin/env python3
"""Exact checks behind the minimal-surface canonical-carrier theorem.

The surface-classification input is theoretical.  This certifier freezes the
restriction-of-scalars rank arithmetic and every lattice calculation in the
Wehler ideal-triangle calibration.
"""

from __future__ import annotations

import sympy as sp


WEHLER_GRAM = sp.Matrix([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
WEHLER_ROOTS = (
    sp.Matrix([-1, 1, 1]),
    sp.Matrix([1, -1, 1]),
    sp.Matrix([1, 1, -1]),
)


def pairing(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return (left.T * WEHLER_GRAM * right)[0]


def reflection_matrix(root: sp.Matrix) -> sp.Matrix:
    """Return x -> x-2(x,r)/(r,r) r in the H_i basis."""

    square = pairing(root, root)
    return sp.eye(3) - sp.Rational(2, 1) * root * (root.T * WEHLER_GRAM) / square


def trace_signature(field_degree: int, hyperbolic_dimension: int) -> tuple[int, int]:
    """Signature of Tr(q) for one (1,n) and d-1 definite real places."""

    assert field_degree >= 1 and hyperbolic_dimension >= 2
    return (1, field_degree * (hyperbolic_dimension + 1) - 1)


def complex_dimension_ceiling(field_degree: int, picard_ceiling: int = 20) -> int:
    assert field_degree >= 1 and picard_ceiling >= 1
    return picard_ceiling // field_degree - 1


def certify() -> None:
    assert WEHLER_GRAM.det() == 16
    assert WEHLER_GRAM.eigenvals() == {sp.Integer(4): 1, sp.Integer(-2): 2}

    root_gram = sp.Matrix(
        [[pairing(left, right) for right in WEHLER_ROOTS] for left in WEHLER_ROOTS]
    )
    assert root_gram == sp.Matrix([[-4, 4, 4], [4, -4, 4], [4, 4, -4]])

    basis = [sp.eye(3).col(i) for i in range(3)]
    for index, root in enumerate(WEHLER_ROOTS):
        reflection = reflection_matrix(root)
        assert all(entry.q == 1 for entry in reflection)
        assert reflection.T * WEHLER_GRAM * reflection == WEHLER_GRAM
        assert reflection * reflection == sp.eye(3)
        assert reflection * root == -root
        for fixed_index in range(3):
            if fixed_index != index:
                assert reflection * basis[fixed_index] == basis[fixed_index]
        # The deck involution satisfies H_i -> -H_i+2H_j+2H_k.
        assert reflection * basis[index] == basis[index] + 2 * root
        # These are level-two integral orthogonal transformations.
        assert all(int(entry) % 2 == (1 if row == col else 0)
                   for row in range(3) for col in range(3)
                   for entry in [reflection[row, col]])

    # The intersection of the other two mirrors is exactly the H_i ideal ray.
    for index, vertex in enumerate(basis):
        assert pairing(vertex, vertex) == 0
        for root_index, root in enumerate(WEHLER_ROOTS):
            assert (pairing(vertex, root) == 0) == (root_index != index)

    # x^2=4(ab+ac+bc), so no integral vector can have square -2.
    a, b, c = sp.symbols("a b c", integer=True)
    vector = sp.Matrix([a, b, c])
    assert sp.expand(pairing(vector, vector)) == 4 * (a * b + a * c + b * c)

    assert trace_signature(1, 2) == (1, 2)
    assert trace_signature(3, 4) == (1, 14)
    assert complex_dimension_ceiling(1) == 19
    assert complex_dimension_ceiling(2) == 9
    assert complex_dimension_ceiling(6) == 2


def main() -> int:
    certify()
    print("PASS: trace rank and Wehler ideal-triangle calibration are exact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
