"""Exact arithmetic in ``Q(sqrt(2))`` and its K3 trace pairing.

The ring of integers is ``Z[sqrt(2)]``.  For a norm-two field lattice the
totally negative codifferent element

``alpha = -(2 + sqrt(2))/4``

has trace ``-1``.  Hence ``Tr(alpha*h)`` sends a field root of square two to
an integral divisor class of square minus two.  The implementation here is
deliberately small: it supplies the exact operations needed by the silver
finite-core certificate without changing the older golden matrix API.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping, Sequence, TypeAlias

from .matrix import Inertia


RationalLike: TypeAlias = int | Fraction


def _fraction(value: RationalLike | str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, (int, str)):
        return Fraction(value)
    raise TypeError(f"expected an integer, Fraction, or rational string; got {value!r}")


def _rational_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _sign_rational(value: Fraction) -> int:
    return (value > 0) - (value < 0)


@dataclass(frozen=True, slots=True)
class Silver:
    """An exact element ``a + b*sqrt(2)`` of ``Q(sqrt(2))``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", _fraction(self.a))
        object.__setattr__(self, "b", _fraction(self.b))

    @classmethod
    def coerce(cls, value: Silver | RationalLike) -> Silver:
        if isinstance(value, cls):
            return value
        return cls(_fraction(value), Fraction(0))

    @classmethod
    def from_json(cls, value: Mapping[str, Any]) -> Silver:
        if set(value) != {"a", "b"}:
            raise ValueError("a serialized silver element must have keys 'a' and 'b'")
        return cls(_fraction(value["a"]), _fraction(value["b"]))

    def to_json(self) -> dict[str, str]:
        return {"a": _rational_string(self.a), "b": _rational_string(self.b)}

    @property
    def is_rational(self) -> bool:
        return self.b == 0

    @property
    def is_integral(self) -> bool:
        return self.a.denominator == 1 and self.b.denominator == 1

    def conjugate(self) -> Silver:
        return Silver(self.a, -self.b)

    def trace(self) -> Fraction:
        return 2 * self.a

    def norm(self) -> Fraction:
        return self.a * self.a - 2 * self.b * self.b

    def sign(self, embedding: int = 0) -> int:
        """Return the exact sign at ``sqrt(2)`` or ``-sqrt(2)``."""

        if embedding not in (0, 1):
            raise ValueError("embedding must be 0 (identity) or 1 (conjugate)")
        if embedding == 1:
            return self.conjugate().sign()
        rational_sign = _sign_rational(self.a)
        radical_sign = _sign_rational(self.b)
        if rational_sign == 0:
            return radical_sign
        if radical_sign == 0 or rational_sign == radical_sign:
            return rational_sign
        comparison = self.a * self.a - 2 * self.b * self.b
        if comparison > 0:
            return rational_sign
        return radical_sign

    def inverse(self) -> Silver:
        field_norm = self.norm()
        if field_norm == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt(2))")
        return Silver(self.a / field_norm, -self.b / field_norm)

    def __add__(self, other: Silver | RationalLike) -> Silver:
        other = Silver.coerce(other)
        return Silver(self.a + other.a, self.b + other.b)

    def __radd__(self, other: Silver | RationalLike) -> Silver:
        return self + other

    def __sub__(self, other: Silver | RationalLike) -> Silver:
        other = Silver.coerce(other)
        return Silver(self.a - other.a, self.b - other.b)

    def __rsub__(self, other: Silver | RationalLike) -> Silver:
        return Silver.coerce(other) - self

    def __neg__(self) -> Silver:
        return Silver(-self.a, -self.b)

    def __mul__(self, other: Silver | RationalLike) -> Silver:
        other = Silver.coerce(other)
        return Silver(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __rmul__(self, other: Silver | RationalLike) -> Silver:
        return self * other

    def __truediv__(self, other: Silver | RationalLike) -> Silver:
        return self * Silver.coerce(other).inverse()

    def __rtruediv__(self, other: Silver | RationalLike) -> Silver:
        return Silver.coerce(other) / self

    def __pow__(self, exponent: int) -> Silver:
        if not isinstance(exponent, int):
            return NotImplemented
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = SILVER_ONE
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    def __bool__(self) -> bool:
        return self.a != 0 or self.b != 0

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Silver):
            return self.a == other.a and self.b == other.b
        if isinstance(other, (int, Fraction)):
            return self.b == 0 and self.a == other
        return NotImplemented

    def __hash__(self) -> int:
        if self.b == 0:
            return hash(self.a)
        return hash((Silver, self.a, self.b))

    def __lt__(self, other: Silver | RationalLike) -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: Silver | RationalLike) -> bool:
        return (self - other).sign() <= 0

    def __gt__(self, other: Silver | RationalLike) -> bool:
        return (self - other).sign() > 0

    def __ge__(self, other: Silver | RationalLike) -> bool:
        return (self - other).sign() >= 0

    def __str__(self) -> str:
        if self.b == 0:
            return _rational_string(self.a)
        if self.a == 0:
            return f"{_rational_string(self.b)}*sqrt(2)"
        sign = "+" if self.b > 0 else "-"
        return (
            f"{_rational_string(self.a)} {sign} "
            f"{_rational_string(abs(self.b))}*sqrt(2)"
        )


SILVER_ZERO = Silver(0)
SILVER_ONE = Silver(1)
SQRT2 = Silver(0, 1)
SILVER_K3_TRACE_SCALAR = Silver(Fraction(-1, 2), Fraction(-1, 4))
SILVER_RING_BASIS = (SILVER_ONE, SQRT2)


SilverMatrix = Sequence[Sequence[Silver | RationalLike]]


def silver_determinant(matrix: SilverMatrix) -> Silver:
    """Exact determinant over ``Q(sqrt(2))``."""

    work = [[Silver.coerce(entry) for entry in row] for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    result = SILVER_ONE
    sign = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return SILVER_ZERO
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column + 1, size):
                work[row][index] -= factor * work[column][index]
    return sign * result


def silver_inertia(matrix: SilverMatrix, *, embedding: int = 0) -> Inertia:
    """Exact symmetric inertia over either real embedding."""

    work = [[Silver.coerce(entry) for entry in row] for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    if any(work[i][j] != work[j][i] for i in range(size) for j in range(i)):
        raise ValueError("inertia is defined only for symmetric matrices")

    def recurse(current: list[list[Silver]]) -> Inertia:
        order = len(current)
        if order == 0:
            return Inertia(0, 0, 0)
        pivot_index = next(
            (index for index in range(order) if current[index][index]), None
        )
        if pivot_index is not None:
            if pivot_index != 0:
                permutation = [
                    pivot_index,
                    *[index for index in range(order) if index != pivot_index],
                ]
                current = [
                    [current[i][j] for j in permutation] for i in permutation
                ]
            pivot = current[0][0]
            remainder = [
                [
                    current[i][j] - current[i][0] * current[0][j] / pivot
                    for j in range(1, order)
                ]
                for i in range(1, order)
            ]
            tail = recurse(remainder)
            if pivot.sign(embedding) > 0:
                return Inertia(tail.positive + 1, tail.negative, tail.zero)
            return Inertia(tail.positive, tail.negative + 1, tail.zero)
        pair = next(
            (
                (i, j)
                for i in range(order)
                for j in range(i + 1, order)
                if current[i][j]
            ),
            None,
        )
        if pair is None:
            return Inertia(0, 0, order)
        first, second = pair
        permutation = [
            first,
            second,
            *[index for index in range(order) if index not in pair],
        ]
        current = [[current[i][j] for j in permutation] for i in permutation]
        off_diagonal = current[0][1]
        remainder = [
            [
                current[i][j]
                - (
                    current[i][0] * current[1][j]
                    + current[i][1] * current[0][j]
                )
                / off_diagonal
                for j in range(2, order)
            ]
            for i in range(2, order)
        ]
        tail = recurse(remainder)
        return Inertia(tail.positive + 1, tail.negative + 1, tail.zero)

    return recurse(work)


def silver_trace_gram(
    field_gram: SilverMatrix,
    alpha: Silver | RationalLike = SILVER_K3_TRACE_SCALAR,
) -> list[list[Fraction]]:
    """Return ``Tr(alpha*h)`` in the basis ``e_i, sqrt(2)e_i``."""

    gram = [[Silver.coerce(entry) for entry in row] for row in field_gram]
    rank = len(gram)
    if any(len(row) != rank for row in gram):
        raise ValueError("field Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(rank) for j in range(i)):
        raise ValueError("field Gram matrix must be symmetric")
    scalar = Silver.coerce(alpha)
    result = [[Fraction(0) for _ in range(2 * rank)] for _ in range(2 * rank)]
    for i in range(rank):
        for j in range(rank):
            for basis_i, multiplier_i in enumerate(SILVER_RING_BASIS):
                for basis_j, multiplier_j in enumerate(SILVER_RING_BASIS):
                    result[2 * i + basis_i][2 * j + basis_j] = (
                        scalar * multiplier_i * gram[i][j] * multiplier_j
                    ).trace()
    return result

