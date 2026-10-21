#!/usr/bin/env python3
"""Small exact Q(sqrt(2)) arithmetic and principal-minor PSD necessities.

An element is the pair (a,b) of Fractions representing a+b*sqrt(2).
No numerical tolerances or external arithmetic packages are used.
"""
from fractions import Fraction as F
from itertools import combinations


def q(a=0, b=0):
    return F(a), F(b)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    return x[0]*y[0]+2*x[1]*y[1], x[0]*y[1]+x[1]*y[0]


def inv(x):
    a, b = x
    norm = a*a-2*b*b
    assert norm != 0, 'division by zero in Q(sqrt(2))'
    return a/norm, -b/norm


def div(x, y):
    return mul(x, inv(y))


def sign(x):
    a, b = x
    if not b:
        return (a > 0)-(a < 0)
    if not a:
        return (b > 0)-(b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    difference = a*a-2*b*b
    assert difference != 0  # sqrt(2) is irrational and a,b are nonzero.
    return ((difference > 0)-(difference < 0)) * (1 if a > 0 else -1)


def encode(x):
    return [str(x[0]), str(x[1])]


def determinant(matrix):
    """Exact row elimination; row swaps are tracked, with no sign pivoting."""
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    a = [row[:] for row in matrix]
    out = q(1)
    for k in range(n):
        row = next((i for i in range(k, n) if a[i][k] != q()), None)
        if row is None:
            return q()
        if row != k:
            a[k], a[row] = a[row], a[k]
            out = neg(out)
        pivot = a[k][k]
        out = mul(out, pivot)
        for i in range(k+1, n):
            if a[i][k] == q():
                continue
            factor = div(a[i][k], pivot)
            for j in range(k+1, n):
                a[i][j] = sub(a[i][j], mul(factor, a[k][j]))
            a[i][k] = q()
    return out


def negative_principal_minor(matrix):
    """Return the first strict negative principal minor and checked count.

    A None witness means every principal minor was nonnegative, exactly
    the positive-semidefinite criterion for the real symmetric matrix.
    """
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    assert all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(n))
    checked = 0
    for size in range(1, n+1):
        for ids in combinations(range(n), size):
            value = determinant([[matrix[i][j] for j in ids] for i in ids])
            checked += 1
            if sign(value) < 0:
                return dict(indices=list(ids), determinant=encode(value)), checked
    return None, checked


def recursive_determinant(matrix):
    """Independent Laplace expansion, used only by the small self-check."""
    if not matrix:
        return q(1)
    if len(matrix) == 1:
        return matrix[0][0]
    out = q()
    for j, x in enumerate(matrix[0]):
        if x == q():
            continue
        minor = [[row[k] for k in range(len(matrix)) if k != j] for row in matrix[1:]]
        term = mul(x, recursive_determinant(minor))
        out = add(out, term if j % 2 == 0 else neg(term))
    return out


def self_check():
    assert mul(q(0, 1), q(0, 1)) == q(2)
    assert inv(q(3, 2)) == q(3, -2)
    # Near-zero opposite-sign elements, including the two Pell signs.
    assert sign(q(99, -70)) == 1 and sign(q(239, -169)) == -1
    assert sign(q(-99, 70)) == -1 and sign(q(-239, 169)) == 1
    checked = 0
    for n in range(1, 7):
        for seed in range(12):
            matrix = [[q(((i+1)*(j+3)+seed*(i+2)) % 7-3,
                         ((i+2)*(j+1)+seed*(j+3)) % 5-2)
                       for j in range(n)] for i in range(n)]
            assert determinant(matrix) == recursive_determinant(matrix)
            checked += 1
    assert negative_principal_minor([[q(1), q(0, 1)], [q(0, 1), q(1)]])[0]['determinant'] == ['-1', '0']
    assert negative_principal_minor([[q(2), q(0, 1)], [q(0, 1), q(1)]])[0] is None
    return checked


if __name__ == '__main__':
    print('PASS_EXACT_Q_SQRT2_HELPER', self_check(), 'determinant comparisons')
