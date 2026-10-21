"""Exact arithmetic in the real cubic field of conductor seven.

We use

``theta = 2*cos(2*pi/7)``,  ``theta^3 + theta^2 - 2*theta - 1 = 0``.

The polynomial has discriminant ``49`` and ``Z[theta]`` is the full ring of
integers.  Its three roots lie respectively in ``(1, 2)``, ``(-1, 0)`` and
``(-2, -1)``.  The Galois generator used below sends

``theta -> theta^2 - 2 -> 1 - theta - theta^2 -> theta``.

This module deliberately has no floating-point dependency.  In particular,
signs at all three real places are decided by rational interval refinement.
That is the missing primitive needed by a cubic-field Cartan/Vinberg backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
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


Interval = tuple[Fraction, Fraction]


def _interval_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def _interval_multiply(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def _minimal_polynomial(value: Fraction) -> Fraction:
    return value**3 + value**2 - 2 * value - 1


@dataclass(frozen=True, slots=True)
class Cubic7:
    """An exact element ``a + b*theta + c*theta^2``.

    Embedding zero is the distinguished one, with ``1 < theta < 2``.
    Embeddings one and two apply the displayed cyclic Galois generator once
    and twice.  Ordering operators use embedding zero.
    """

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    c: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", _fraction(self.a))
        object.__setattr__(self, "b", _fraction(self.b))
        object.__setattr__(self, "c", _fraction(self.c))

    @classmethod
    def coerce(cls, value: Cubic7 | RationalLike) -> Cubic7:
        if isinstance(value, cls):
            return value
        return cls(_fraction(value))

    @classmethod
    def from_json(cls, value: Mapping[str, Any]) -> Cubic7:
        if set(value) != {"a", "b", "c"}:
            raise ValueError(
                "a serialized conductor-seven element must have keys 'a', 'b', 'c'"
            )
        return cls(_fraction(value["a"]), _fraction(value["b"]), _fraction(value["c"]))

    def to_json(self) -> dict[str, str]:
        return {
            "a": _rational_string(self.a),
            "b": _rational_string(self.b),
            "c": _rational_string(self.c),
        }

    @property
    def is_rational(self) -> bool:
        return self.b == 0 and self.c == 0

    @property
    def is_integral(self) -> bool:
        """Whether the element belongs to the full ring ``Z[theta]``."""

        return all(coefficient.denominator == 1 for coefficient in (self.a, self.b, self.c))

    def conjugate(self, power: int = 1) -> Cubic7:
        """Apply the cyclic Galois generator ``power`` times."""

        power %= 3
        result = self
        for _ in range(power):
            # sigma(theta)=theta^2-2 and sigma(theta)^2=3-theta-theta^2.
            result = Cubic7(
                result.a - 2 * result.b + 3 * result.c,
                -result.c,
                result.b - result.c,
            )
        return result

    def trace(self) -> Fraction:
        # Tr(1)=3, Tr(theta)=-1, Tr(theta^2)=5.
        return 3 * self.a - self.b + 5 * self.c

    def norm(self) -> Fraction:
        product = self * self.conjugate(1) * self.conjugate(2)
        if not product.is_rational:  # pragma: no cover - guards implementation drift
            raise ArithmeticError("the Galois product did not reduce to a rational number")
        return product.a

    def _identity_sign(self) -> int:
        if not self:
            return 0

        # The distinguished root is the unique root of the minimal polynomial
        # in (1,2).  A nonzero quadratic cannot vanish there, since the cubic
        # is irreducible.  Rational interval evaluation therefore terminates.
        lower = Fraction(1)
        upper = Fraction(2)
        while True:
            theta_interval = (lower, upper)
            value_interval = _interval_add(
                (self.a, self.a),
                _interval_multiply(
                    theta_interval,
                    _interval_add(
                        (self.b, self.b),
                        _interval_multiply((self.c, self.c), theta_interval),
                    ),
                ),
            )
            if value_interval[0] > 0:
                return 1
            if value_interval[1] < 0:
                return -1
            midpoint = (lower + upper) / 2
            midpoint_value = _minimal_polynomial(midpoint)
            if midpoint_value == 0:  # impossible for a rational midpoint
                raise ArithmeticError("the irreducible cubic acquired a rational root")
            # p(1)<0<p(2), and the interval contains only the largest root.
            if midpoint_value < 0:
                lower = midpoint
            else:
                upper = midpoint

    @lru_cache(maxsize=None)
    def sign(self, embedding: int = 0) -> int:
        """Return ``-1``, ``0`` or ``1`` at any real embedding, exactly."""

        if embedding not in (0, 1, 2):
            raise ValueError("embedding must be 0, 1, or 2")
        return self.conjugate(embedding)._identity_sign()

    def inverse(self) -> Cubic7:
        field_norm = self.norm()
        if field_norm == 0:
            raise ZeroDivisionError("division by zero in Q(cos(2*pi/7))")
        return self.conjugate(1) * self.conjugate(2) / field_norm

    def __add__(self, other: Cubic7 | RationalLike) -> Cubic7:
        other = Cubic7.coerce(other)
        return Cubic7(self.a + other.a, self.b + other.b, self.c + other.c)

    def __radd__(self, other: Cubic7 | RationalLike) -> Cubic7:
        return self + other

    def __sub__(self, other: Cubic7 | RationalLike) -> Cubic7:
        other = Cubic7.coerce(other)
        return Cubic7(self.a - other.a, self.b - other.b, self.c - other.c)

    def __rsub__(self, other: Cubic7 | RationalLike) -> Cubic7:
        return Cubic7.coerce(other) - self

    def __neg__(self) -> Cubic7:
        return Cubic7(-self.a, -self.b, -self.c)

    def __mul__(self, other: Cubic7 | RationalLike) -> Cubic7:
        other = Cubic7.coerce(other)
        a, b, c = self.a, self.b, self.c
        d, e, f = other.a, other.b, other.c
        # theta^3=1+2*theta-theta^2, theta^4=-1-theta+3*theta^2.
        return Cubic7(
            a * d + b * f + c * e - c * f,
            a * e + b * d + 2 * b * f + 2 * c * e - c * f,
            a * f + b * e + c * d - b * f - c * e + 3 * c * f,
        )

    def __rmul__(self, other: Cubic7 | RationalLike) -> Cubic7:
        return self * other

    def __truediv__(self, other: Cubic7 | RationalLike) -> Cubic7:
        other = Cubic7.coerce(other)
        # Avoid recursively invoking ``inverse`` when dividing by a rational.
        if other.is_rational:
            if other.a == 0:
                raise ZeroDivisionError("division by zero in Q(cos(2*pi/7))")
            return Cubic7(self.a / other.a, self.b / other.a, self.c / other.a)
        return self * other.inverse()

    def __rtruediv__(self, other: Cubic7 | RationalLike) -> Cubic7:
        return Cubic7.coerce(other) / self

    def __pow__(self, exponent: int) -> Cubic7:
        if not isinstance(exponent, int):
            return NotImplemented
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = CUBIC7_ONE
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    def __bool__(self) -> bool:
        return self.a != 0 or self.b != 0 or self.c != 0

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Cubic7):
            return (self.a, self.b, self.c) == (other.a, other.b, other.c)
        if isinstance(other, (int, Fraction)):
            return self.b == 0 and self.c == 0 and self.a == other
        return NotImplemented

    def __hash__(self) -> int:
        if self.is_rational:
            return hash(self.a)
        return hash((Cubic7, self.a, self.b, self.c))

    def __lt__(self, other: Cubic7 | RationalLike) -> bool:
        return (self - other).sign(0) < 0

    def __le__(self, other: Cubic7 | RationalLike) -> bool:
        return (self - other).sign(0) <= 0

    def __gt__(self, other: Cubic7 | RationalLike) -> bool:
        return (self - other).sign(0) > 0

    def __ge__(self, other: Cubic7 | RationalLike) -> bool:
        return (self - other).sign(0) >= 0

    def __str__(self) -> str:
        terms: list[str] = []
        for coefficient, name in ((self.a, ""), (self.b, "theta"), (self.c, "theta^2")):
            if not coefficient:
                continue
            absolute = abs(coefficient)
            body = _rational_string(absolute) if not name else f"{_rational_string(absolute)}*{name}"
            if not terms:
                terms.append(body if coefficient > 0 else f"-{body}")
            else:
                terms.append(("+ " if coefficient > 0 else "- ") + body)
        return " ".join(terms) if terms else "0"


CUBIC7_ZERO = Cubic7(0)
CUBIC7_ONE = Cubic7(1)
THETA7 = Cubic7(0, 1)
"""``2*cos(2*pi/7)``, with sign pattern ``(+,-,-)``."""

RHO7 = Cubic7(-1, 1, 1)
"""``2*cos(pi/7)`` at the distinguished place."""

LABEL7_PRODUCT = Cubic7(2, 1)
"""``4*cos(pi/7)^2 = RHO7^2 = 2 + THETA7``."""

LABEL14_PRODUCT = Cubic7(1, 1, 1)
"""``4*cos(pi/14)^2 = 2 + RHO7 = 1 + THETA7 + THETA7^2``."""

CUBIC7_RING_BASIS = (CUBIC7_ONE, THETA7, THETA7**2)
CUBIC7_MINIMAL_POLYNOMIAL = (Fraction(-1), Fraction(-2), Fraction(1), Fraction(1))
CUBIC7_FIELD_DISCRIMINANT = 49


Cubic7Matrix = Sequence[Sequence[Cubic7 | RationalLike]]


def cubic7_determinant(matrix: Cubic7Matrix) -> Cubic7:
    """Exact determinant over the conductor-seven cubic field."""

    work = [[Cubic7.coerce(entry) for entry in row] for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    result = CUBIC7_ONE
    sign = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return CUBIC7_ZERO
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


def cubic7_inertia(matrix: Cubic7Matrix, *, embedding: int = 0) -> Inertia:
    """Exact symmetric inertia at any of the three real places."""

    work = [[Cubic7.coerce(entry) for entry in row] for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    if any(work[i][j] != work[j][i] for i in range(size) for j in range(i)):
        raise ValueError("inertia is defined only for symmetric matrices")

    def recurse(current: list[list[Cubic7]]) -> Inertia:
        order = len(current)
        if order == 0:
            return Inertia(0, 0, 0)
        pivot_index = next((index for index in range(order) if current[index][index]), None)
        if pivot_index is not None:
            if pivot_index != 0:
                permutation = [
                    pivot_index,
                    *[index for index in range(order) if index != pivot_index],
                ]
                current = [[current[i][j] for j in permutation] for i in permutation]
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
        permutation = [first, second, *[i for i in range(order) if i not in pair]]
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


__all__ = [
    "CUBIC7_FIELD_DISCRIMINANT",
    "CUBIC7_MINIMAL_POLYNOMIAL",
    "CUBIC7_ONE",
    "CUBIC7_RING_BASIS",
    "CUBIC7_ZERO",
    "LABEL7_PRODUCT",
    "LABEL14_PRODUCT",
    "RHO7",
    "THETA7",
    "Cubic7",
    "cubic7_determinant",
    "cubic7_inertia",
]
