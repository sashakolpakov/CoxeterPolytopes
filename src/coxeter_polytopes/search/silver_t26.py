"""Exact arithmetic and Vinberg controls for the silver ``T26`` host.

``T26`` is the even determinant-three rank-26 lattice in the Nebe--Sloane
catalogue.  This module studies

``<-(2*sqrt(2)-1)> + T26``

over ``Z[sqrt(2)]``.  It supplies the literal catalogue Gram matrix, the
complete eight-class primitive mirror-norm alphabet (including the
discriminant-three gate), exact binary Gram controls, and the fixed dual
direction used by the bounded continuation certificate.  It also freezes an
exact 66-contact vertex of the full minimal-shell polar.  A strict shrink of
that vertex is a much stronger first-packet hyperideal witness than any
catalogue coordinate direction.
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
Vector: TypeAlias = tuple[Silver, ...]

T26_SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb4.html"
T26_MINIMUM = 4
T26_KISSING_NUMBER = 117_936
T26_GROUP_ORDER = 1_268_047_872
T26_SMITH_INVARIANTS = (1,) * 25 + (3,)
T26_RADICAL_THREE_INDEX = 3**25
T26_RADICAL_THREE_MINIMUM = 24
T26_RADICAL_THREE_KISSING_NUMBER = 1_638

T26_LOWER_TRIANGLE: tuple[tuple[int, ...], ...] = (
    (4,),
    (-2, 4),
    (-2, 0, 4),
    (1, 1, 0, 4),
    (-2, 1, 0, -2, 4),
    (-2, 1, 0, 0, 0, 4),
    (1, -2, -1, 0, 0, -1, 4),
    (1, 1, -1, 2, -2, 1, 0, 4),
    (-2, 1, 0, 0, 2, 2, 0, -1, 4),
    (-2, 2, 0, -1, 2, 0, 0, 0, 0, 4),
    (1, 1, -1, 2, -2, 0, -1, 1, -1, 0, 4),
    (1, 1, 0, 2, -2, -1, 0, 1, -1, -1, 1, 4),
    (-2, 0, 2, -2, 2, 0, 0, -1, 0, 2, -2, -2, 4),
    (-1, 2, 1, 2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 4),
    (-2, 2, 1, -1, 2, 0, -1, 0, 0, 2, -1, 0, 2, 0, 4),
    (-1, 2, 1, 2, 0, 0, 0, 1, 0, 1, 0, 2, 0, 2, 1, 4),
    (-2, 2, 0, -1, 1, 2, -1, 1, 1, 1, -1, 0, 1, 0, 2, 1, 4),
    (-1, -1, 1, 0, 1, 1, 0, -1, 1, 0, -1, -2, 1, 0, 0, 0, -1, 4),
    (-2, 0, 0, -2, 1, 2, 0, 0, 1, 1, -1, -2, 1, -1, 1, -1, 2, 1, 4),
    (-2, 1, 2, 1, 1, 1, 0, 0, 1, 1, -1, 0, 1, 2, 1, 2, 0, 2, 0, 4),
    (-2, 0, 2, 0, 0, 2, -1, 0, 1, 0, 0, -1, 1, 0, 0, 0, 1, 1, 1, 1, 4),
    (-1, 2, -1, -1, 1, 0, 0, 1, 0, 2, 0, 0, 1, 0, 2, 0, 2, -2, 1, -1, -1, 4),
    (-1, -1, 1, 0, 0, 0, 2, 0, 0, 1, -1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 4),
    (-2, 0, 0, -1, 2, 2, 0, -1, 2, 1, -1, -2, 1, -1, 1, 0, 1, 2, 2, 1, 1, 0, 1, 4),
    (1, 1, -1, 2, -2, 1, 0, 2, 0, -1, 1, 2, -2, 1, -1, 1, 0, -1, -1, 0, -1, 0, -1, -1, 4),
    (1, -2, 1, 0, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, -1, -1, -2, 1, -1, 0, 1, -2, 0, 0, -1, 4),
)


def _bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(map(int, row)) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("a determinant matrix must be square")
    previous = 1
    sign = 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        value = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * value
                    - work[row][column] * work[column][index]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = value
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def _rank_mod_prime(matrix: IntegerMatrix, prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(row, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [(inverse * value) % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (work[index][j] - factor * work[row][j]) % prime
                for j in range(len(work[0]))
            ]
        row += 1
    return row


@lru_cache(maxsize=1)
def t26_gram() -> IntegerMatrix:
    if tuple(map(len, T26_LOWER_TRIANGLE)) != tuple(range(1, 27)):
        raise AssertionError("the T26 lower triangle has the wrong row lengths")
    return tuple(
        tuple(T26_LOWER_TRIANGLE[max(i, j)][min(i, j)] for j in range(26))
        for i in range(26)
    )


T26_RADICAL_THREE_VECTOR: tuple[int, ...] = (
    2, 1, 1, 2, 2, 0, 0, 1, 0, 0, 2, 0, 0,
    0, 1, 2, 0, 0, 2, 2, 1, 1, 1, 0, 1, 1,
)


@lru_cache(maxsize=1)
def t26_radical_three_basis() -> IntegerMatrix:
    """Column basis of ``{x : (x,T26) is contained in 3Z}``."""

    # The last coordinate of the displayed radical vector is one, so it and
    # 3*e_0,...,3*e_24 form a basis of the full preimage of the mod-three
    # radical.
    return tuple(
        tuple(
            T26_RADICAL_THREE_VECTOR[row] if column == 0
            else (3 if row == column - 1 else 0)
            for column in range(26)
        )
        for row in range(26)
    )


def t26_pair(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    if len(left) != 26 or len(right) != 26:
        raise ValueError("T26 coordinate vectors have length 26")
    gram = t26_gram()
    return sum(
        Fraction(left[i]) * gram[i][j] * Fraction(right[j])
        for i in range(26)
        for j in range(26)
        if gram[i][j]
    )


@lru_cache(maxsize=1)
def dual_direction() -> tuple[Fraction, ...]:
    inverse = _fraction_inverse(
        tuple(tuple(Fraction(value) for value in row) for row in t26_gram())
    )
    # Zero-based coordinate one, as in the certificate brief.
    return tuple(inverse[i][1] for i in range(26))


def t26_catalogue_audit() -> bool:
    gram = t26_gram()
    h = dual_direction()
    radical = T26_RADICAL_THREE_VECTOR
    return (
        all(gram[i][j] == gram[j][i] for i in range(26) for j in range(26))
        and all(gram[i][i] % 2 == 0 for i in range(26))
        and _bareiss_determinant(gram) == 3
        and _rank_mod_prime(gram, 3) == 25
        and all(sum(gram[i][j] * radical[j] for j in range(26)) % 3 == 0 for i in range(26))
        and abs(_bareiss_determinant(t26_radical_three_basis())) == 3**25
        and tuple(sum(Fraction(gram[i][j]) * h[j] for j in range(26)) for i in range(26))
        == (Fraction(0), Fraction(1), *(Fraction(0) for _ in range(24)))
        and t26_pair(h, h) == Fraction(212, 3)
    )


RHO = 2 * SQRT2 - SILVER_ONE
EPSILON = SILVER_ONE + SQRT2
P7_GENERATOR = 3 + SQRT2
FIRST_PACKET_HEIGHT = 5 + 4 * SQRT2
FIRST_PACKET_NORM = Silver(2)
FIRST_PACKET_CONTROL = Silver(2, 1)
FIRST_PACKET_TARGET = Silver(12, 8)
FIRST_PACKET_WALL_COUNT = T26_KISSING_NUMBER


@dataclass(frozen=True, slots=True)
class T26NormClass:
    norm: Silver
    dyadic: bool
    p7: bool
    discriminant_three: bool


NORM_CLASSES = tuple(
    T26NormClass(
        Silver(2 if dyadic else 1)
        * (P7_GENERATOR if p7 else SILVER_ONE)
        * (3 if discriminant_three else 1),
        dyadic,
        p7,
        discriminant_three,
    )
    for dyadic, p7, discriminant_three in product((False, True), repeat=3)
)

COMPLETE_NORM_CLASS_PROOF = (
    "Put pi=sqrt(2), rho=2*pi-1 and delta=3+pi=epsilon*rho. For a "
    "primitive root r=(x,y), q(r) divides 2I(r), with "
    "I(r)=(rho*x,(y,T26)). Away from pi, rho and the determinant prime 3, "
    "unimodularity excludes every prime divisor of q. At pi the spatial "
    "lattice is even and unimodular, so the exponent is zero or two. At "
    "rho it is zero or one; exponent one forces y in rho*T26 and is exact "
    "by primitivity. At the inert rational prime 3 the Smith form 1^25,3 "
    "allows exponent zero or one; exponent one forces x in 3O and both "
    "coefficients of y into S3={v:(v,T26) subset 3Z}. The local choices are "
    "independent. Class number one and U^+=U^2 give the eight canonical "
    "representatives (1 or 2)*(1 or 3+pi)*(1 or 3)."
)


@dataclass(frozen=True, slots=True)
class T26GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, a: int, b: int) -> int:
        return (
            a * a * self.u_square
            + b * b * self.v_square
            + 2 * a * b * self.u_dot_v
        )


def candidate_gram_triples(
    target: Silver, *, discriminant_three: bool = False
) -> tuple[T26GramTriple, ...]:
    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    minimum = T26_RADICAL_THREE_MINIMUM if discriminant_three else T26_MINIMUM
    output = []
    for v_square in range(rational // 2 + 1):
        u_square = rational - 2 * v_square
        if u_square < 0 or u_square % 2 or v_square % 2:
            continue
        if u_square not in (0,) and u_square < minimum:
            continue
        if v_square not in (0,) and v_square < minimum:
            continue
        if dot * dot > u_square * v_square:
            continue
        if discriminant_three and any(value % 3 for value in (u_square, v_square, dot)):
            continue
        triple = T26GramTriple(u_square, v_square, dot)
        if any(
            0 < triple.combination_square(a, b) < minimum
            for a in range(-3, 4)
            for b in range(-3, 4)
            if a or b
        ):
            continue
        output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class T26Control:
    norm_class: T26NormClass
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_controls(height_bound: Silver = FIRST_PACKET_HEIGHT) -> tuple[T26Control, ...]:
    output = set()
    for norm_class in NORM_CLASSES:
        norm = norm_class.norm
        for control in _bounded_for_caps((norm * height_bound / RHO, norm / (-RHO))):
            if control.sign(0) <= 0 or not (2 * RHO * control / norm).is_integral:
                continue
            target = norm + RHO * control * control
            height = RHO * control * control / norm
            if target.sign(1) < 0 or height > height_bound:
                continue
            output.add(T26Control(norm_class, height, control, target))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control)))


def reduced_target(control: T26Control) -> Silver:
    return control.target / (RHO * RHO) if control.norm_class.p7 else control.target


def gated_gram_triples(control: T26Control) -> tuple[T26GramTriple, ...]:
    return candidate_gram_triples(
        reduced_target(control),
        discriminant_three=control.norm_class.discriminant_three,
    )


def height_zero_gated_gram_triples(norm_class: T26NormClass) -> tuple[T26GramTriple, ...]:
    target = norm_class.norm / (RHO * RHO) if norm_class.p7 else norm_class.norm
    return candidate_gram_triples(target, discriminant_three=norm_class.discriminant_three)


WITNESS_SCALE = Silver(Fraction(999, 4000), Fraction(-999, 16000))
WITNESS_NORM = Silver(
    Fraction(190_682_159, 32_000_000),
    Fraction(-33_631_351, 8_000_000),
)


def witness() -> Vector:
    return (SILVER_ONE, *(WITNESS_SCALE * value for value in dual_direction()))


def spatial_pair_silver(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    if len(left) != 26 or len(right) != 26:
        raise ValueError("T26 spatial field vectors have length 26")
    gram = t26_gram()
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(26)
            for j in range(26)
            if gram[i][j]
        ),
        start=Silver(0),
    )


def lorentz_pair(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    if len(left) != 27 or len(right) != 27:
        raise ValueError("T26 host vectors have length 27")
    return -RHO * left[0] * right[0] + spatial_pair_silver(left[1:], right[1:])


def witness_audit() -> bool:
    z = witness()
    return (
        t26_catalogue_audit()
        and WITNESS_SCALE == Silver(Fraction(999, 4000), Fraction(-999, 16000))
        and lorentz_pair(z, z) == WITNESS_NORM
        and WITNESS_NORM.sign(0) > 0
    )


# A vertex of the polar of the complete norm-four shell.  For every minimal
# b in T26 one has |(POLAR_VERTEX,b)|<=1, with equality for 66 sign pairs.
# The freezer independently enumerates and hashes all 58,968 PARI sign
# representatives.  The simple Gram image below also makes its exact square
# inexpensive to audit without loading the shell.
POLAR_VERTEX: tuple[Fraction, ...] = tuple(
    Fraction(value, 30)
    for value in (
        43,
        119,
        50,
        -11,
        1,
        9,
        60,
        5,
        -27,
        -3,
        -20,
        30,
        0,
        -18,
        -43,
        -35,
        6,
        30,
        4,
        -11,
        2,
        -1,
        -16,
        30,
        -19,
        14,
    )
)
POLAR_VERTEX_SQUARE = Fraction(97, 150)
POLAR_VERTEX_CONTACT_COUNT = 66
POLAR_VERTEX_SHELL_MAXIMUM = Fraction(1)
POLAR_WALL_BOUND = Silver(4) - SQRT2  # rho*sqrt(2)
POLAR_STRICT_SCALE = Fraction(999, 1000)
POLAR_BOUNDARY_NORM = Silver(Fraction(316, 25), Fraction(-538, 75))
POLAR_STRICT_NORM = Silver(
    Fraction(315_418_291, 25_000_000),
    Fraction(-44_768_699, 6_250_000),
)
POLAR_STRICT_FIRST_PACKET_MARGIN = Silver(
    Fraction(-1, 500), Fraction(-3, 1000)
)


def polar_vertex_audit() -> bool:
    image = tuple(
        sum(Fraction(t26_gram()[i][j]) * POLAR_VERTEX[j] for j in range(26))
        for i in range(26)
    )
    expected_image = tuple(
        Fraction(value, 10)
        for value in (
            -1, -1, 1, -2, 3, 0, 3, -2, 1, 0, -3, -2, 1,
            -1, 0, -1, 0, 1, 1, 2, 0, 1, 0, 2, 0, 1,
        )
    )
    return (
        image == expected_image
        and t26_pair(POLAR_VERTEX, POLAR_VERTEX) == POLAR_VERTEX_SQUARE
        and POLAR_WALL_BOUND == RHO * SQRT2
        and Silver(POLAR_VERTEX_SQUARE) * POLAR_WALL_BOUND**2 - RHO
        == POLAR_BOUNDARY_NORM
    )


def polar_strict_witness() -> Vector:
    scale = POLAR_STRICT_SCALE * POLAR_WALL_BOUND
    return (SILVER_ONE, *(scale * value for value in POLAR_VERTEX))


def polar_strict_witness_audit() -> bool:
    z = polar_strict_witness()
    # First roots are (2+sqrt(2),(1+sqrt(2))*b).  At shell support one,
    # their maximum pairing is -(1+sqrt(2))*rho*sqrt(2)/1000.
    maximum = (
        -RHO * FIRST_PACKET_CONTROL
        + EPSILON * POLAR_STRICT_SCALE * POLAR_WALL_BOUND
    )
    return (
        polar_vertex_audit()
        and lorentz_pair(z, z) == POLAR_STRICT_NORM
        and POLAR_STRICT_NORM.sign(0) > 0
        and maximum == POLAR_STRICT_FIRST_PACKET_MARGIN
        and maximum.sign(0) < 0
    )


@dataclass(frozen=True, slots=True)
class DualPairingMaximum:
    maximum: Silver
    maximizing_u_pairing: int
    maximizing_v_pairing: int
    retained_count: int


@lru_cache(maxsize=None)
def dual_pairing_maximum(triple: T26GramTriple) -> DualPairingMaximum:
    """PSD upper bound for ``h.(u+sqrt(2)v)``, with ``h^2=212/3``."""

    s, t, dot = triple.u_square, triple.v_square, triple.u_dot_v
    h_square = Fraction(212, 3)
    p_bound = isqrt((212 * s) // 3) + 2
    q_bound = isqrt((212 * t) // 3) + 2
    best: Silver | None = None
    best_pair: tuple[int, int] | None = None
    count = 0
    for p in range(-p_bound, p_bound + 1):
        if p * p > h_square * s:
            continue
        for q in range(-q_bound, q_bound + 1):
            if q * q > h_square * t:
                continue
            determinant = (
                h_square * (s * t - dot * dot)
                - t * p * p
                - s * q * q
                + 2 * dot * p * q
            )
            if determinant < 0:
                continue
            count += 1
            value = Silver(p, q)
            if best is None or value > best:
                best = value
                best_pair = (p, q)
    if best is None or best_pair is None:
        raise AssertionError("the augmented dual Gram enumeration was empty")
    return DualPairingMaximum(best, *best_pair, count)


__all__ = [
    "COMPLETE_NORM_CLASS_PROOF",
    "DualPairingMaximum",
    "EPSILON",
    "FIRST_PACKET_HEIGHT",
    "NORM_CLASSES",
    "P7_GENERATOR",
    "POLAR_BOUNDARY_NORM",
    "POLAR_STRICT_FIRST_PACKET_MARGIN",
    "POLAR_STRICT_NORM",
    "POLAR_STRICT_SCALE",
    "POLAR_VERTEX",
    "POLAR_VERTEX_CONTACT_COUNT",
    "POLAR_VERTEX_SHELL_MAXIMUM",
    "POLAR_VERTEX_SQUARE",
    "POLAR_WALL_BOUND",
    "RHO",
    "T26Control",
    "T26GramTriple",
    "T26NormClass",
    "T26_GROUP_ORDER",
    "T26_KISSING_NUMBER",
    "T26_LOWER_TRIANGLE",
    "T26_MINIMUM",
    "T26_RADICAL_THREE_INDEX",
    "T26_RADICAL_THREE_KISSING_NUMBER",
    "T26_RADICAL_THREE_MINIMUM",
    "T26_RADICAL_THREE_VECTOR",
    "T26_SMITH_INVARIANTS",
    "T26_SOURCE_URL",
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_TARGET",
    "FIRST_PACKET_WALL_COUNT",
    "WITNESS_NORM",
    "WITNESS_SCALE",
    "candidate_gram_triples",
    "dual_direction",
    "dual_pairing_maximum",
    "enumerate_controls",
    "gated_gram_triples",
    "height_zero_gated_gram_triples",
    "lorentz_pair",
    "polar_strict_witness",
    "polar_strict_witness_audit",
    "polar_vertex_audit",
    "reduced_target",
    "spatial_pair_silver",
    "t26_catalogue_audit",
    "t26_gram",
    "t26_pair",
    "t26_radical_three_basis",
    "witness",
    "witness_audit",
]
