"""Exact screen for the rank-24 modular Gram sequence.

The two endpoint matrices are literal Nebe--Sloane catalogue Grams in their
shared basis.  Put ``D=G1-G0`` and ``G_n=G0+nD``.  The module freezes the
family identity and the complete bounded ``n=4`` silver first-packet screen;
it does not claim Vinberg termination.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

P24_9_SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/P24.9.html"
L24_11_SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/L_24.11.html"
P24_9_GRAM_SHA256 = "1a283c09dbf9dbdb9eb34aec7057ecd594351b9dea51042111efa6d2295c5f50"
L24_11_GRAM_SHA256 = "25f70337de4bcb2bac84e8811e3c658a9a2b27f625f8dea5da7f171679dc75e7"

# The complete 24 by 24 endpoint matrices, row-major.  Each continuation line
# is deliberately just Python whitespace; parsing is invariant under wrapping.
_P24_9_GRAM_TEXT = """
8 -1 2 0 4 -1 -1 -1 1 -1 3 2 2 -4 3 3 0 0 -2 3 0 4 1 1
-1 8 -1 4 -1 -1 0 4 3 -2 0 -1 0 -2 -2 0 -4 0 -4 0 -2 -2 2 2
2 -1 8 -1 2 2 0 -2 -2 -1 -2 0 4 -3 -1 0 -1 -2 3 4 -2 2 -4 3
0 4 -1 8 -1 0 0 2 2 -3 3 1 1 -1 -2 0 -4 1 -3 0 -2 -2 2 0
4 -1 2 -1 8 -2 -1 -1 0 -2 -1 0 3 -4 0 3 1 -1 -2 2 -1 4 1 2
-1 -1 2 0 -2 8 2 0 -2 2 -1 -1 -2 1 -1 -3 0 0 0 3 4 -4 -2 -1
-1 0 0 0 -1 2 8 -2 2 -2 -2 -4 -1 1 -4 -1 -1 -4 -1 2 3 -4 -4 4
-1 4 -2 2 -1 0 -2 8 1 2 1 -1 -2 -1 2 -3 -2 0 -3 -1 -1 -2 4 -2
1 3 -2 2 0 -2 2 1 8 -3 3 -1 -1 0 -1 0 -3 -2 -3 0 0 -2 2 3
-1 -2 -1 -3 -2 2 -2 2 -3 8 0 1 -2 3 4 -3 4 2 2 0 2 -1 0 -4
3 0 -2 3 -1 -1 -2 1 3 0 8 4 0 1 2 1 -1 2 -2 -2 -1 1 4 -2
2 -1 0 1 0 -1 -4 -1 -1 1 4 8 1 0 2 2 0 4 2 -2 -3 4 3 -4
2 0 4 1 3 -2 -1 -2 -1 -2 0 1 8 -2 -2 1 -2 1 2 2 -4 4 0 3
-4 -2 -3 -1 -4 1 1 -1 0 3 1 0 -2 8 1 -3 2 2 2 -1 3 -4 -1 -3
3 -2 -1 -2 0 -1 -4 2 -1 4 2 2 -2 1 8 -1 2 2 1 0 1 2 2 -4
3 0 0 0 3 -3 -1 -3 0 -3 1 2 1 -3 -1 8 1 1 -2 -2 0 4 1 1
0 -4 -1 -4 1 0 -1 -2 -3 4 -1 0 -2 2 2 1 8 1 1 -1 3 1 -2 -3
0 0 -2 1 -1 0 -4 0 -2 2 2 4 1 2 2 1 1 8 1 0 1 2 4 -4
-2 -4 3 -3 -2 0 -1 -3 -3 2 -2 2 2 2 1 -2 1 1 8 0 -2 2 -3 -1
3 0 4 0 2 3 2 -1 0 0 -2 -2 2 -1 0 -2 -1 0 0 8 2 -1 -3 3
0 -2 -2 -2 -1 4 3 -1 0 2 -1 -3 -4 3 1 0 3 1 -2 2 8 -4 -1 -1
4 -2 2 -2 4 -4 -4 -2 -2 -1 1 4 4 -4 2 4 1 2 2 -1 -4 8 2 0
1 2 -4 2 1 -2 -4 4 2 0 4 3 0 -1 2 1 -2 4 -3 -3 -1 2 8 -3
1 2 3 0 2 -1 4 -2 3 -4 -2 -4 3 -3 -4 1 -3 -4 -1 3 -1 0 -3 8
"""

_L24_11_GRAM_TEXT = """
12 -1 3 0 6 -2 -1 -2 2 -2 4 3 3 -6 4 5 0 0 -3 5 0 6 1 2
-1 12 -1 6 -2 -1 0 6 5 -3 1 -1 0 -3 -3 0 -6 0 -6 0 -3 -3 3 3
3 -1 12 -1 3 3 0 -3 -3 -2 -3 0 6 -5 -2 1 -2 -3 4 6 -3 3 -6 5
0 6 -1 12 -2 0 0 3 3 -4 5 2 2 -1 -3 0 -6 2 -4 0 -3 -3 3 0
6 -2 3 -2 12 -3 -1 -2 0 -3 -2 0 4 -6 0 5 2 -2 -3 3 -1 6 1 3
-2 -1 3 0 -3 12 3 0 -3 3 -2 -2 -3 1 -2 -4 0 0 0 4 6 -6 -3 -1
-1 0 0 0 -1 3 12 -3 3 -3 -3 -6 -2 1 -6 -1 -1 -6 -2 3 5 -6 -6 6
-2 6 -3 3 -2 0 -3 12 2 3 2 -1 -3 -1 3 -5 -3 0 -4 -2 -2 -3 6 -3
2 5 -3 3 0 -3 3 2 12 -4 5 -1 -2 0 -1 0 -4 -3 -5 0 0 -3 3 4
-2 -3 -2 -4 -3 3 -3 3 -4 12 0 1 -3 5 6 -5 6 3 3 0 3 -2 0 -6
4 1 -3 5 -2 -2 -3 2 5 0 12 6 0 2 3 1 -2 3 -3 -3 -2 1 6 -3
3 -1 0 2 0 -2 -6 -1 -1 1 6 12 2 0 3 3 -1 6 3 -3 -5 6 5 -6
3 0 6 2 4 -3 -2 -3 -2 -3 0 2 12 -3 -3 2 -3 2 3 3 -6 6 0 4
-6 -3 -5 -1 -6 1 1 -1 0 5 2 0 -3 12 2 -5 3 3 3 -2 4 -6 -1 -5
4 -3 -2 -3 0 -2 -6 3 -1 6 3 3 -3 2 12 -2 3 3 2 0 1 3 3 -6
5 0 1 0 5 -4 -1 -5 0 -5 1 3 2 -5 -2 12 1 1 -3 -2 0 6 1 2
0 -6 -2 -6 2 0 -1 -3 -4 6 -2 -1 -3 3 3 1 12 1 1 -1 5 1 -3 -4
0 0 -3 2 -2 0 -6 0 -3 3 3 6 2 3 3 1 1 12 2 0 1 3 6 -6
-3 -6 4 -4 -3 0 -2 -4 -5 3 -3 3 3 3 2 -3 1 2 12 0 -3 3 -4 -2
5 0 6 0 3 4 3 -2 0 0 -3 -3 3 -2 0 -2 -1 0 0 12 3 -1 -5 5
0 -3 -3 -3 -1 6 5 -2 0 3 -2 -5 -6 4 1 0 5 1 -3 3 12 -6 -2 -1
6 -3 3 -3 6 -6 -6 -3 -3 -2 1 6 6 -6 3 6 1 3 3 -1 -6 12 3 0
1 3 -6 3 1 -3 -6 6 3 0 6 5 0 -1 3 1 -3 6 -4 -5 -2 3 12 -5
2 3 5 0 3 -1 6 -3 4 -6 -3 -6 4 -5 -6 2 -4 -6 -2 5 -1 0 -5 12
"""


def _parse_gram(text: str) -> IntegerMatrix:
    entries = tuple(map(int, text.split()))
    if len(entries) != 24 * 24:
        raise AssertionError("rank-24 Gram text has the wrong length")
    matrix = tuple(tuple(entries[24 * i + j] for j in range(24)) for i in range(24))
    if any(matrix[i][j] != matrix[j][i] for i in range(24) for j in range(24)):
        raise AssertionError("rank-24 Gram text is not symmetric")
    return matrix


@lru_cache(maxsize=1)
def endpoint_grams() -> tuple[IntegerMatrix, IntegerMatrix]:
    return _parse_gram(_P24_9_GRAM_TEXT), _parse_gram(_L24_11_GRAM_TEXT)


def _add(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(left[i][j] + right[i][j] for j in range(24)) for i in range(24))


def _scale(value: int, matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(value * matrix[i][j] for j in range(24)) for i in range(24))


def _subtract(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> IntegerMatrix:
    return _add(left, _scale(-1, right))


def _transpose(matrix: Sequence[Sequence[Fraction | int]]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(Fraction(matrix[j][i]) for j in range(24)) for i in range(24))


def _matmul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(24)) for j in range(24))
        for i in range(24)
    )


@lru_cache(maxsize=None)
def sequence_gram(index: int) -> IntegerMatrix:
    if index < 0:
        raise ValueError("sequence index must be nonnegative")
    g0, g1 = endpoint_grams()
    return _add(g0, _scale(index, _subtract(g1, g0)))


def level(index: int) -> int:
    return index * index + 5 * index + 5


def minimum(index: int) -> int:
    return 8 + 4 * index


def family_matrix_identity() -> bool:
    """Verify ``D G0^-1 D = D-G0/5`` exactly."""

    g0, g1 = endpoint_grams()
    difference = _subtract(g1, g0)
    inverse = _fraction_inverse(tuple(tuple(Fraction(x) for x in row) for row in g0))
    left = _matmul(_matmul(difference, inverse), difference)
    right = tuple(
        tuple(Fraction(difference[i][j]) - Fraction(g0[i][j], 5) for j in range(24))
        for i in range(24)
    )
    return left == right


def _rref_kernel_mod_prime(matrix: IntegerMatrix, prime: int) -> tuple[IntegerMatrix, tuple[int, ...]]:
    work = [[entry % prime for entry in row] for row in matrix]
    pivots: list[int] = []
    row = 0
    for column in range(24):
        pivot = next((i for i in range(row, 24) if work[i][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [(inverse * value) % prime for value in work[row]]
        for i in range(24):
            if i == row or not work[i][column]:
                continue
            factor = work[i][column]
            work[i] = [
                (work[i][j] - factor * work[row][j]) % prime for j in range(24)
            ]
        pivots.append(column)
        row += 1
    free = tuple(i for i in range(24) if i not in pivots)
    kernel = []
    for column in free:
        vector = [0] * 24
        vector[column] = 1
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = (-work[pivot_row][column]) % prime
        kernel.append(tuple(vector))
    return tuple(kernel), tuple(pivots)


@lru_cache(maxsize=1)
def n4_radical_basis() -> IntegerMatrix:
    """Column basis of ``{v:G4*v=0 mod 41}`` in its full integral preimage."""

    kernel, pivots = _rref_kernel_mod_prime(sequence_gram(4), 41)
    columns = list(kernel) + [tuple(41 if i == pivot else 0 for i in range(24)) for pivot in pivots]
    return tuple(tuple(columns[j][i] for j in range(24)) for i in range(24))


@lru_cache(maxsize=1)
def n4_scaled_radical_gram() -> IntegerMatrix:
    basis = n4_radical_basis()
    raw = _matmul(_matmul(_transpose(basis), sequence_gram(4)), basis)
    if any(value.denominator != 1 or value.numerator % 41 for row in raw for value in row):
        raise AssertionError("the mod-41 radical pairing is not divisible by 41")
    return tuple(tuple(int(value / 41) for value in row) for row in raw)


N4_COEFFICIENT = Silver(-3, 5)
N4_PI = Silver(7, 2)
N4_PI_BAR = Silver(7, -2)
N4_FIRST_HEIGHT = Silver(806, 570)
N4_FIRST_CONTROL = Silver(14, 10)
N4_FIRST_TARGET = Silver(1614, 1140)
N4_DUAL_INDEX = 1
N4_DUAL_SQUARE = Fraction(24, 41)
N4_WITNESS_SCALE = Silver(3)
N4_HEIGHT_BOUND = Silver(5000)


@dataclass(frozen=True, slots=True)
class N4NormClass:
    norm: Silver
    pi: bool
    pi_bar: bool

    @property
    def spatial_divisor(self) -> Silver:
        return (N4_PI if self.pi else SILVER_ONE) * (N4_PI_BAR if self.pi_bar else SILVER_ONE)


N4_NORM_CLASSES = tuple(
    N4NormClass(Silver(dyadic) * (N4_PI if pi else 1) * (N4_PI_BAR if pi_bar else 1), pi, pi_bar)
    for dyadic, pi, pi_bar in product((1, 2), (False, True), (False, True))
)


@dataclass(frozen=True, slots=True)
class GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, first: int, second: int) -> int:
        return (
            first * first * self.u_square + second * second * self.v_square
            + 2 * first * second * self.u_dot_v
        )


def _has_short_combination(triple: GramTriple) -> bool:
    determinant = triple.u_square * triple.v_square - triple.u_dot_v**2
    if determinant <= 0:
        return determinant < 0
    first_bound = isqrt((23 * triple.v_square) // determinant) + 1
    second_bound = isqrt((23 * triple.u_square) // determinant) + 1
    return any(
        0 < triple.combination_square(first, second) < 24
        for first in range(-first_bound, first_bound + 1)
        for second in range(-second_bound, second_bound + 1)
        if first or second
    )


@lru_cache(maxsize=None)
def candidate_triples(target: Silver) -> tuple[GramTriple, ...]:
    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational, dot = int(target.a), int(target.b / 2)
    output = []
    for v_square in range(0, rational // 2 + 1, 2):
        u_square = rational - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if (u_square and u_square < 24) or (v_square and v_square < 24):
            continue
        if dot * dot > u_square * v_square:
            continue
        triple = GramTriple(u_square, v_square, dot)
        if not _has_short_combination(triple):
            output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class N4Control:
    norm_class: N4NormClass
    height: Silver
    control: Silver
    target: Silver
    reduced_target: Silver


@lru_cache(maxsize=None)
def n4_controls(height_bound: Silver = N4_HEIGHT_BOUND) -> tuple[N4Control, ...]:
    output: set[N4Control] = set()
    for norm_class in N4_NORM_CLASSES:
        norm = norm_class.norm
        for control in _bounded_for_caps((norm * height_bound / N4_COEFFICIENT, norm / (-N4_COEFFICIENT))):
            if control.sign(0) <= 0 or not (2 * N4_COEFFICIENT * control / norm).is_integral:
                continue
            # At the conjugate split prime the time line is unimodular.
            if norm_class.pi_bar and not (control / N4_PI_BAR).is_integral:
                continue
            target = norm + N4_COEFFICIENT * control * control
            height = N4_COEFFICIENT * control * control / norm
            if target.sign(1) < 0 or height > height_bound:
                continue
            divisor = norm_class.spatial_divisor
            output.add(N4Control(norm_class, height, control, target, target / (divisor * divisor)))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control)))


def n4_gated_triples(control: N4Control) -> tuple[GramTriple, ...]:
    return candidate_triples(control.reduced_target)


SHORT_SHELL_COORDINATE_BOUNDS = {24: 2, 26: 2, 28: 4}


@dataclass(frozen=True, slots=True)
class DirectionMaximum:
    maximum: Silver
    coordinates: tuple[int, int]
    augmented_determinant: Fraction
    retained_count: int


@lru_cache(maxsize=None)
def n4_direction_maximum(triple: GramTriple) -> DirectionMaximum:
    s, t, dot = triple.u_square, triple.v_square, triple.u_dot_v
    p_bound = isqrt((24 * s) // 41) + 1
    q_bound = isqrt((24 * t) // 41) + 1
    retained: list[tuple[Silver, int, int, Fraction]] = []
    for p in range(-p_bound, p_bound + 1):
        for q in range(-q_bound, q_bound + 1):
            determinant = (
                N4_DUAL_SQUARE * (s * t - dot * dot) - t * p * p
                + 2 * dot * p * q - s * q * q
            )
            if determinant < 0:
                continue
            if any(
                norm in SHORT_SHELL_COORDINATE_BOUNDS
                and abs(first * p + second * q) > SHORT_SHELL_COORDINATE_BOUNDS[norm]
                for first in range(-5, 6)
                for second in range(-5, 6)
                if (first or second)
                for norm in (triple.combination_square(first, second),)
            ):
                continue
            retained.append((Silver(p, q), p, q, determinant))
    if not retained:
        raise AssertionError("the exact augmented coordinate screen is empty")
    maximum, p, q, determinant = max(retained)
    return DirectionMaximum(maximum, (p, q), determinant, len(retained))


@dataclass(frozen=True, slots=True)
class N4FirstPacketAudit:
    control_count: int
    distinct_height_count: int
    first_control: N4Control
    triples: tuple[GramTriple, ...]
    maxima: tuple[DirectionMaximum, ...]
    witness_norm: Silver
    margins: tuple[Silver, ...]

    @property
    def strict_hyperideal(self) -> bool:
        return self.witness_norm.sign(0) > 0 and all(value.sign(0) > 0 for value in self.margins)


def n4_first_packet_audit() -> N4FirstPacketAudit:
    controls = n4_controls()
    live = tuple((row, n4_gated_triples(row)) for row in controls if n4_gated_triples(row))
    if not live:
        raise AssertionError("no n=4 necessary packet below the frozen bound")
    row, triples = live[0]
    expected = (GramTriple(798, 408, 570), GramTriple(794, 410, 570))
    if (
        row.height != N4_FIRST_HEIGHT or row.norm_class.norm != 2
        or row.control != N4_FIRST_CONTROL or row.target != N4_FIRST_TARGET
        or triples != expected
    ):
        raise AssertionError("the n=4 first packet drifted")
    maxima = tuple(n4_direction_maximum(triple) for triple in triples)
    witness_norm = N4_WITNESS_SCALE**2 * N4_DUAL_SQUARE - N4_COEFFICIENT
    slab = N4_COEFFICIENT * row.control
    margins = tuple(slab - N4_WITNESS_SCALE * maximum.maximum for maximum in maxima)
    return N4FirstPacketAudit(
        len(controls), len({control.height for control in controls}), row, triples,
        maxima, witness_norm, margins,
    )

