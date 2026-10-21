"""Finite seed manifest for non-unimodular diagonal conductor-seven forms.

Put ``theta = 2*cos(2*pi/7)`` and ``O = Z[theta]``.  The discovery forms are

``diag(-theta, 1, ..., 1, t)``

in ranks five and six, where ``t`` is totally positive.  This module records
one deliberately small *ideal-radical* sweep: every product of distinct atoms
from the inert prime ``(2)``, the ramified prime over 7, and one of the three
primes over the first split rational prime 13, provided its radical norm is at
most 128.  Other prime directions (including the inert primes 3 and 5 and the
split primes from 29 onward) are a second tier.  They make the root-norm
divisor set or the number of jobs larger before the first three genuinely
different split-prime directions have been tested.

The manifest is canonical only for this stated atom/support box; the norm
bound by itself does not define the box, and this is not a classification of
diagonal forms.  The three entries above a split rational prime are
Galois-conjugate, not generally isometric after fixing the
distinguished Lorentzian embedding, so all three are retained.

AlVin's RC7 basis is ``(theta, sigma(theta), sigma^2(theta))`` rather than
``(1, theta, theta^2)``.  Since

``1 = -(theta + sigma(theta) + sigma^2(theta))``, the exact conversion is

``a + b*theta + c*theta^2 -> [b-a-2c, -a-c, -a-2c]``.
"""

from __future__ import annotations

from dataclasses import dataclass

from coxeter_polytopes.exact.cubic7 import Cubic7


@dataclass(frozen=True, slots=True)
class Cubic7DiagonalSeed:
    """One coefficient square class and its intended search priority."""

    seed_id: str
    coefficient: Cubic7
    radical_norm: int
    support: str
    priority: int
    estimated_cost: str
    rationale: str
    calibration: bool = False

    @property
    def field_norm(self) -> int:
        norm = self.coefficient.norm()
        if norm.denominator != 1:
            raise ArithmeticError("an integral seed acquired nonintegral norm")
        return int(norm)

    @property
    def rc7_coefficient(self) -> tuple[int, int, int]:
        return cubic7_to_alvin_rc7(self.coefficient)

    def alvin_qf(self, rank: int) -> str:
        """Return AlVin's exact ``-qf`` argument in the requested rank."""

        if rank not in (5, 6):
            raise ValueError("the initial manifest only targets ranks five and six")
        coefficient = _format_rc7(self.rc7_coefficient)
        return ",".join(("[-1,0,0]", *("1" for _ in range(rank - 2)), coefficient))

    def alvin_command(self, rank: int, binary: str = "./alvin") -> str:
        return f"{binary} -k rc7 -qf '{self.alvin_qf(rank)}' -ip"


def cubic7_to_alvin_rc7(value: Cubic7) -> tuple[int, int, int]:
    """Convert the power basis to AlVin's three-cosine integral basis."""

    if not value.is_integral:
        raise ValueError("AlVin RC7 coefficients must be algebraic integers")
    a, b, c = (int(value.a), int(value.b), int(value.c))
    return b - a - 2 * c, -a - c, -a - 2 * c


def alvin_rc7_to_cubic7(value: tuple[int, int, int]) -> Cubic7:
    """Convert an AlVin three-cosine tuple to ``a+b*theta+c*theta^2``."""

    x, y, z = value
    # theta_1 = theta^2-2 and theta_2 = 1-theta-theta^2.
    return Cubic7(-2 * y + z, x - z, y - z)


def _format_rc7(value: tuple[int, int, int]) -> str:
    return f"[{value[0]},{value[1]},{value[2]}]"


# Reduced representatives modulo multiplication by unit squares.  The
# representatives over 13 are distinguished by the roots 7, 10 and 8 of the
# minimal polynomial modulo 13, respectively.  The mixed representatives
# retain the same order.
_SEEDS = (
    Cubic7DiagonalSeed(
        "unit-control",
        Cubic7(1),
        1,
        "unit",
        0,
        "control",
        "unimodular baseline; detects backend or convention regressions",
        True,
    ),
    Cubic7DiagonalSeed(
        "ramified-7",
        Cubic7(2, -1),
        7,
        "p7",
        1,
        "low",
        "smallest non-unimodular ideal direction; 7 is totally ramified",
    ),
    Cubic7DiagonalSeed(
        "inert-2",
        Cubic7(2),
        8,
        "(2)",
        1,
        "low",
        "smallest inert-prime direction and a distinct determinant square class",
    ),
    Cubic7DiagonalSeed(
        "split-13-r7",
        Cubic7(0, -1, 2),
        13,
        "p13(theta=7)",
        2,
        "medium",
        "first split-prime direction; one of three distinguished real-place classes",
    ),
    Cubic7DiagonalSeed(
        "split-13-r10",
        Cubic7(1, -1, 1),
        13,
        "p13(theta=10)",
        2,
        "medium",
        "first split-prime direction; one of three distinguished real-place classes",
    ),
    Cubic7DiagonalSeed(
        "split-13-r8",
        Cubic7(1, 0, 1),
        13,
        "p13(theta=8)",
        2,
        "medium",
        "first split-prime direction; one of three distinguished real-place classes",
    ),
    Cubic7DiagonalSeed(
        "mixed-2x7",
        Cubic7(2, 2, 2),
        56,
        "(2)*p7",
        3,
        "medium",
        "first mixed-prime direction; tests interaction of inert and ramified support",
    ),
    Cubic7DiagonalSeed(
        "mixed-7x13-r7",
        Cubic7(3, 0, 1),
        91,
        "p7*p13(theta=7)",
        4,
        "high",
        "mixed ramified/split support below the radical-norm cutoff",
    ),
    Cubic7DiagonalSeed(
        "mixed-7x13-r10",
        Cubic7(0, -1, 4),
        91,
        "p7*p13(theta=10)",
        4,
        "high",
        "mixed ramified/split support below the radical-norm cutoff",
    ),
    Cubic7DiagonalSeed(
        "mixed-7x13-r8",
        Cubic7(1, -3, 3),
        91,
        "p7*p13(theta=8)",
        4,
        "high",
        "mixed ramified/split support below the radical-norm cutoff",
    ),
    Cubic7DiagonalSeed(
        "mixed-2x13-r7",
        Cubic7(0, -2, 4),
        104,
        "(2)*p13(theta=7)",
        5,
        "high",
        "mixed inert/split support below the radical-norm cutoff",
    ),
    Cubic7DiagonalSeed(
        "mixed-2x13-r10",
        Cubic7(2, -2, 2),
        104,
        "(2)*p13(theta=10)",
        5,
        "high",
        "mixed inert/split support below the radical-norm cutoff",
    ),
    Cubic7DiagonalSeed(
        "mixed-2x13-r8",
        Cubic7(2, 0, 2),
        104,
        "(2)*p13(theta=8)",
        5,
        "high",
        "mixed inert/split support below the radical-norm cutoff",
    ),
)


def initial_cubic7_diagonal_seeds(
    *, include_control: bool = True
) -> tuple[Cubic7DiagonalSeed, ...]:
    """Return the ordered finite initial sweep."""

    if include_control:
        return _SEEDS
    return tuple(seed for seed in _SEEDS if not seed.calibration)


__all__ = [
    "Cubic7DiagonalSeed",
    "alvin_rc7_to_cubic7",
    "cubic7_to_alvin_rc7",
    "initial_cubic7_diagonal_seeds",
]
