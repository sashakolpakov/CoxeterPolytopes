"""Exact arithmetic in Q(phi), where phi^2 = phi + 1."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping, TypeAlias


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
class Golden:
    """An element ``a + b*phi`` of Q(sqrt(5)).

    Coefficients are canonical :class:`fractions.Fraction` objects.  Ordering
    uses the identity real embedding.  Use :meth:`sign` with ``embedding=1``
    when checking the conjugate embedding of an arithmetic lattice.
    """

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", _fraction(self.a))
        object.__setattr__(self, "b", _fraction(self.b))

    @classmethod
    def coerce(cls, value: Golden | RationalLike) -> Golden:
        if isinstance(value, cls):
            return value
        return cls(_fraction(value), Fraction(0))

    @classmethod
    def from_json(cls, value: Mapping[str, Any]) -> Golden:
        if set(value) != {"a", "b"}:
            raise ValueError("a serialized golden element must have exactly keys 'a' and 'b'")
        return cls(_fraction(value["a"]), _fraction(value["b"]))

    def to_json(self) -> dict[str, str]:
        return {"a": _rational_string(self.a), "b": _rational_string(self.b)}

    @property
    def is_rational(self) -> bool:
        return self.b == 0

    @property
    def is_integral(self) -> bool:
        """Whether this element belongs to Z[phi], the full ring of integers."""

        return self.a.denominator == 1 and self.b.denominator == 1

    def conjugate(self) -> Golden:
        """Apply phi -> 1 - phi."""

        return Golden(self.a + self.b, -self.b)

    def trace(self) -> Fraction:
        return 2 * self.a + self.b

    def norm(self) -> Fraction:
        return self.a * self.a + self.a * self.b - self.b * self.b

    def sign(self, embedding: int = 0) -> int:
        """Return -1, 0, or 1 at a real embedding, without approximation."""

        if embedding not in (0, 1):
            raise ValueError("embedding must be 0 (identity) or 1 (conjugate)")
        if embedding == 1:
            return self.conjugate().sign()

        # 2(a + b*phi) = (2a+b) + b*sqrt(5).  If the two terms
        # have opposite signs, compare their squares; equality can occur only
        # when both are zero because sqrt(5) is irrational.
        rational_part = 2 * self.a + self.b
        radical_coefficient = self.b
        rational_sign = _sign_rational(rational_part)
        radical_sign = _sign_rational(radical_coefficient)

        if rational_sign == 0:
            return radical_sign
        if radical_sign == 0 or rational_sign == radical_sign:
            return rational_sign

        comparison = rational_part * rational_part - 5 * radical_coefficient * radical_coefficient
        if comparison > 0:
            return rational_sign
        return radical_sign

    def inverse(self) -> Golden:
        field_norm = self.norm()
        if field_norm == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt(5))")
        conjugate = self.conjugate()
        return Golden(conjugate.a / field_norm, conjugate.b / field_norm)

    def __add__(self, other: Golden | RationalLike) -> Golden:
        other = Golden.coerce(other)
        return Golden(self.a + other.a, self.b + other.b)

    def __radd__(self, other: Golden | RationalLike) -> Golden:
        return self + other

    def __sub__(self, other: Golden | RationalLike) -> Golden:
        other = Golden.coerce(other)
        return Golden(self.a - other.a, self.b - other.b)

    def __rsub__(self, other: Golden | RationalLike) -> Golden:
        return Golden.coerce(other) - self

    def __neg__(self) -> Golden:
        return Golden(-self.a, -self.b)

    def __mul__(self, other: Golden | RationalLike) -> Golden:
        other = Golden.coerce(other)
        # (a+b*phi)(c+d*phi), using phi^2 = phi+1.
        return Golden(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def __rmul__(self, other: Golden | RationalLike) -> Golden:
        return self * other

    def __truediv__(self, other: Golden | RationalLike) -> Golden:
        return self * Golden.coerce(other).inverse()

    def __rtruediv__(self, other: Golden | RationalLike) -> Golden:
        return Golden.coerce(other) / self

    def __pow__(self, exponent: int) -> Golden:
        if not isinstance(exponent, int):
            return NotImplemented
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = ONE
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
        if isinstance(other, Golden):
            return self.a == other.a and self.b == other.b
        if isinstance(other, (int, Fraction)):
            return self.b == 0 and self.a == other
        return NotImplemented

    def __hash__(self) -> int:
        # Cross-type rational equality requires the same hash as Fraction/int.
        if self.b == 0:
            return hash(self.a)
        return hash((Golden, self.a, self.b))

    def __lt__(self, other: Golden | RationalLike) -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: Golden | RationalLike) -> bool:
        return (self - other).sign() <= 0

    def __gt__(self, other: Golden | RationalLike) -> bool:
        return (self - other).sign() > 0

    def __ge__(self, other: Golden | RationalLike) -> bool:
        return (self - other).sign() >= 0

    def __str__(self) -> str:
        if self.b == 0:
            return _rational_string(self.a)
        if self.a == 0:
            return f"{_rational_string(self.b)}*phi"
        sign = "+" if self.b > 0 else "-"
        return f"{_rational_string(self.a)} {sign} {_rational_string(abs(self.b))}*phi"


ZERO = Golden(0)
ONE = Golden(1)
PHI = Golden(0, 1)
