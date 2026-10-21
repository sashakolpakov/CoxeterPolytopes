"""Bounded silver screen for the rank-28 lattice ``DualExtremal(28,3)``.

The spatial lattice is the determinant-``3^26``, minimum-eight lattice from
the Nebe--Sloane catalogue.  This module studies

``<-(2*sqrt(2)-1)> + DualExtremal(28,3)``

over ``Z[sqrt(2)]``.  It deliberately certifies a *necessary* first Vinberg
packet.  A strict hyperideal point survives that whole superset, so no shell
realisation assumption is needed for the bounded infinite-volume result.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]

SOURCE_URL = (
    "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/"
    "DualExtremal_28_3.html"
)
DIMENSION = 28
DETERMINANT = 3**26
MINIMUM = 8
KISSING_NUMBER = 1_260
SMITH_INVARIANTS = (1, 1) + (3,) * 26
RADICAL_THREE_INDEX = 9
RADICAL_THREE_MINIMUM = 12
RADICAL_THREE_KISSING_NUMBER = 65_520

_GRAM_TEXT = """
8 4 0 0 2 -2 2 -2 -3 1 -1 3 2 3 4 3 3 0 2 0 4 -2 2 1 -3 -2 3 3
4 8 0 -3 4 -4 -2 2 -3 2 -2 3 -2 0 -1 0 0 -3 4 -3 -1 -4 -2 2 0 -4 0 0
0 0 12 -6 3 -3 3 3 6 -3 -6 0 3 3 -3 3 0 3 -3 -3 3 0 3 3 -3 3 -3 -6
0 -3 -6 12 -6 6 3 -6 -6 -3 3 3 3 0 6 0 0 -3 -3 3 0 3 3 0 3 0 0 6
2 4 3 -6 10 -4 0 5 3 -1 -3 0 -2 2 -3 -1 -2 -1 0 -2 2 0 -3 -1 0 0 3 0
-2 -4 -3 6 -4 10 0 -5 -3 1 0 3 -1 -2 3 -2 2 1 0 2 1 0 3 -2 0 3 -3 0
2 -2 3 3 0 0 10 -3 0 -3 1 3 3 4 2 1 -1 1 -5 -1 6 5 1 3 -3 -1 0 0
-2 2 3 -6 5 -5 -3 10 6 -2 -3 -3 -1 1 -3 -2 -1 -2 0 -1 -2 0 -3 -2 0 0 3 0
-3 -3 6 -6 3 -3 0 6 12 -3 -3 -3 0 0 -3 0 0 3 -3 0 0 0 0 0 -3 3 0 -3
1 2 -3 -3 -1 1 -3 -2 -3 10 0 0 -4 -2 0 -2 2 1 6 -1 1 -3 0 1 0 -3 0 -3
-1 -2 -6 3 -3 0 1 -3 -3 0 10 -3 0 -2 -1 -2 -1 1 1 2 -3 2 -5 -3 0 -1 0 3
3 3 0 3 0 3 3 -3 -3 0 -3 12 -3 0 3 0 3 -3 0 0 3 0 3 3 0 0 0 0
2 -2 3 3 -2 -1 3 -1 0 -4 0 -3 10 2 3 5 1 2 -3 1 2 3 3 2 0 0 0 3
3 0 3 0 2 -2 4 1 0 -2 -2 0 2 8 2 2 1 -1 -2 1 4 2 1 1 -3 -1 3 0
4 -1 -3 6 -3 3 2 -3 -3 0 -1 3 3 2 10 2 4 -1 -1 4 3 1 5 0 0 -2 3 6
3 0 3 0 -1 -2 1 -2 0 -2 -2 0 5 2 2 8 1 2 -2 1 1 -1 4 4 0 -1 0 0
3 0 0 0 -2 2 -1 -1 0 2 -1 3 1 1 4 1 8 1 2 2 2 -2 2 -1 -3 1 0 0
0 -3 3 -3 -1 1 1 -2 3 1 1 -3 2 -1 -1 2 1 8 1 1 1 -1 1 1 -3 2 -3 -3
2 4 -3 -3 0 0 -5 0 -3 6 1 0 -3 -2 -1 -2 2 1 10 -1 -3 -4 -2 0 0 -1 0 0
0 -3 -3 3 -2 2 -1 -1 0 -1 2 0 1 1 4 1 2 1 -1 8 -1 1 2 -1 0 1 3 3
4 -1 3 0 2 1 6 -2 0 1 -3 3 2 4 3 1 2 1 -3 -1 10 3 3 1 -3 0 3 0
-2 -4 0 3 0 0 5 0 0 -3 2 0 3 2 1 -1 -2 -1 -4 1 3 10 -1 0 3 1 3 3
2 -2 3 3 -3 3 1 -3 0 0 -5 3 3 1 5 4 2 1 -2 2 3 -1 10 3 0 2 0 0
1 2 3 0 -1 -2 3 -2 0 1 -3 3 2 1 0 4 -1 1 0 -1 1 0 3 10 0 -3 -3 -3
-3 0 -3 3 0 0 -3 0 -3 0 0 0 0 -3 0 0 -3 -3 0 0 -3 3 0 0 12 0 3 3
-2 -4 3 0 0 3 -1 0 3 -3 -1 0 0 -1 -2 -1 1 2 -1 1 0 1 2 -3 0 10 0 0
3 0 -3 0 3 -3 0 3 0 0 0 0 0 3 3 0 0 -3 0 3 3 3 0 -3 3 0 12 6
3 0 -6 6 0 0 0 0 -3 -3 3 0 3 0 6 0 0 -3 0 3 0 3 0 -3 3 0 6 12
"""


@lru_cache(maxsize=1)
def gram() -> IntegerMatrix:
    rows = tuple(tuple(map(int, row.split())) for row in _GRAM_TEXT.splitlines() if row)
    if len(rows) != DIMENSION or any(len(row) != DIMENSION for row in rows):
        raise AssertionError("the literal catalogue Gram matrix is malformed")
    if any(rows[i][j] != rows[j][i] for i in range(DIMENSION) for j in range(DIMENSION)):
        raise AssertionError("the literal catalogue Gram matrix is not symmetric")
    return rows


def _bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(map(int, row)) for row in matrix]
    previous = 1
    sign = 1
    for column in range(len(work) - 1):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        value = work[column][column]
        for row in range(column + 1, len(work)):
            for index in range(column + 1, len(work)):
                numerator = work[row][index] * value - work[row][column] * work[column][index]
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = value
        for row in range(column + 1, len(work)):
            work[row][column] = 0
    return sign * work[-1][-1]


def _rank_mod_three(matrix: IntegerMatrix) -> int:
    work = [[entry % 3 for entry in row] for row in matrix]
    rank = 0
    for column in range(DIMENSION):
        pivot = next((row for row in range(rank, DIMENSION) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, 3)
        work[rank] = [(inverse * entry) % 3 for entry in work[rank]]
        for row in range(DIMENSION):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (work[row][index] - factor * work[rank][index]) % 3
                for index in range(DIMENSION)
            ]
        rank += 1
    return rank


def catalogue_audit() -> bool:
    matrix = gram()
    return (
        _bareiss_determinant(matrix) == DETERMINANT
        and _rank_mod_three(matrix) == 2
        and all(matrix[i][i] % 2 == 0 for i in range(DIMENSION))
    )


RHO = 2 * SQRT2 - SILVER_ONE
P7 = 3 + SQRT2
FIRST_NECESSARY_HEIGHT = 93 + 66 * SQRT2
FIRST_CUTTING_HEIGHT = 362 + 256 * SQRT2
FIRST_CUTTING_CONTROL = 14 + 10 * SQRT2
FIRST_CUTTING_NORM = Silver(2)


@dataclass(frozen=True, slots=True)
class NormClass:
    norm: Silver
    dyadic: bool
    p7: bool
    discriminant_three: bool


NORM_CLASSES = tuple(
    NormClass(
        Silver(2 if dyadic else 1)
        * (P7 if p7 else SILVER_ONE)
        * (3 if discriminant_three else 1),
        dyadic,
        p7,
        discriminant_three,
    )
    for dyadic, p7, discriminant_three in product((False, True), repeat=3)
)

COMPLETE_NORM_CLASS_PROOF = (
    "For a primitive root q divides 2I, where I=(rho*x,(y,L)). Away from "
    "the ramified dyadic prime, rho above seven, and the determinant prime "
    "three, primitivity excludes all valuations. Evenness gives dyadic "
    "exponent zero or two. The squarefree rho and 3-primary Smith factors "
    "give exponent zero or one; the latter force y in rho*L and respectively "
    "y in S3={v:(v,L) subset 3Z}. Class number one and positive units modulo "
    "squares leave (1 or 2)*(1 or 3+sqrt(2))*(1 or 3)."
)


@dataclass(frozen=True, slots=True)
class Control:
    norm_class: NormClass
    height: Silver
    control: Silver
    target: Silver


@dataclass(frozen=True, slots=True)
class GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, a: int, b: int) -> int:
        return (
            a * a * self.u_square
            + b * b * self.v_square
            + 2 * a * b * self.u_dot_v
        )


FIRST_CUTTING_TRIPLE = GramTriple(362, 182, 256)

# The first catalogue basis vector is a minimum vector.  An affine exact CVP
# slice produces the displayed v with v^2=182 and d.v=37; u=2d+v realizes the
# unique augmented determinant-zero threat at the first cutting height.
FIRST_CUTTING_D = (1,) + (0,) * 27
FIRST_CUTTING_V = (
    4,
    0,
    1,
    3,
    1,
    -1,
    -2,
    1,
    0,
    1,
    1,
    1,
    0,
    0,
    2,
    1,
    -1,
    1,
    0,
    -1,
    0,
    0,
    -1,
    0,
    -1,
    1,
    1,
    -2,
)
FIRST_CUTTING_U = tuple(
    FIRST_CUTTING_V[index] + 2 * FIRST_CUTTING_D[index]
    for index in range(DIMENSION)
)
REPLACEMENT_WITNESS_SCALE = Fraction(12, 25)

# The next accepted cutter.  An exact affine slice gives the displayed
# vector with v^2=828 and d.v=81; u=v+4d realizes the determinant-zero q=2
# threat at H=1629+1152*sqrt(2).
SECOND_CUTTING_HEIGHT = Silver(1629, 1152)
SECOND_CUTTING_CONTROL = Silver(30, 21)
SECOND_CUTTING_TRIPLE = GramTriple(1604, 828, 1152)
SECOND_CUTTING_V = (
    10, 3, 0, 0, -1, 1, 1, 0, 0, 0, 0, -1, -1, 0,
    -1, 0, 1, 1, -1, 0, -1, 1, 1, 0, 0, 0, 1, 0,
)
SECOND_CUTTING_U = tuple(
    SECOND_CUTTING_V[index] + 4 * FIRST_CUTTING_D[index]
    for index in range(DIMENSION)
)
SECOND_CUTTING_W = tuple(
    SECOND_CUTTING_V[index] - 10 * FIRST_CUTTING_D[index]
    for index in range(DIMENSION)
)
SECOND_REPLACEMENT_WITNESS_SCALE = Fraction(479, 1000)

# This is an exact realization of the next raw necessary threat.  Its
# chronological Vinberg acceptance remains open: do not list it among the
# accepted cutters.  Keeping it here prevents a prior nearest-shell vector
# of square 3088 from being confused with the required square-3092 vector.
THIRD_CANDIDATE_HEIGHT = Silver(6149, 4348)
THIRD_CANDIDATE_CONTROL = Silver(58, 41)
THIRD_CANDIDATE_TRIPLE = GramTriple(6116, 3092, 4348)
THIRD_CANDIDATE_W = (
    0, -1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0,
)
THIRD_CANDIDATE_V = tuple(
    THIRD_CANDIDATE_W[index] + 20 * FIRST_CUTTING_D[index]
    for index in range(DIMENSION)
)
THIRD_CANDIDATE_U = tuple(
    THIRD_CANDIDATE_V[index] + 8 * FIRST_CUTTING_D[index]
    for index in range(DIMENSION)
)
THIRD_DIAGNOSTIC_WITNESS_SCALE = Fraction(957, 2000)


@lru_cache(maxsize=None)
def controls(bound: Silver = FIRST_NECESSARY_HEIGHT) -> tuple[Control, ...]:
    output: set[Control] = set()
    for norm_class in NORM_CLASSES:
        q = norm_class.norm
        for x in _bounded_for_caps((q * bound / RHO, q / (-RHO))):
            if x.sign(0) <= 0 or not (2 * RHO * x / q).is_integral:
                continue
            target = q + RHO * x * x
            height = RHO * x * x / q
            if target.sign(1) < 0 or height > bound:
                continue
            output.add(Control(norm_class, height, x, target))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control)))


def reduced_target(control: Control) -> Silver:
    return control.target / (RHO * RHO) if control.norm_class.p7 else control.target


def _has_combination_below(triple: GramTriple, minimum: int) -> bool:
    """Test every binary integral combination below ``minimum``.

    Completing the square gives the displayed finite coefficient bounds in
    the positive-definite case.  In rank one the primitive direction has
    square ``gcd(u^2,v^2,u.v)``.  This replaces the former fixed ``[-3,3]``
    box, which is not a complete minimum gate at later heights.
    """

    u_square = triple.u_square
    v_square = triple.v_square
    dot = triple.u_dot_v
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(u_square, v_square), abs(dot)) < minimum
    u_bound = isqrt(((minimum - 1) * v_square) // determinant) + 1
    v_bound = isqrt(((minimum - 1) * u_square) // determinant) + 1
    return any(
        0 < triple.combination_square(a, b) < minimum
        for a in range(-u_bound, u_bound + 1)
        for b in range(-v_bound, v_bound + 1)
        if a or b
    )


@lru_cache(maxsize=None)
def necessary_triples(control: Control) -> tuple[GramTriple, ...]:
    target = reduced_target(control)
    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    minimum = RADICAL_THREE_MINIMUM if control.norm_class.discriminant_three else MINIMUM
    output: list[GramTriple] = []
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
        if control.norm_class.discriminant_three and any(
            value % 3 for value in (u_square, v_square, dot)
        ):
            continue
        triple = GramTriple(u_square, v_square, dot)
        if _has_combination_below(triple, minimum):
            continue
        output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class PairingMaximum:
    maximum: Silver
    maximizing_pair: tuple[int, int]
    augmented_determinant: int
    retained_extensions: int


@lru_cache(maxsize=None)
def norm_eight_pairing_maximum(triple: GramTriple) -> PairingMaximum:
    """Maximise ``d.(u+sqrt(2)v)`` over a rigorous local Gram superset."""

    u2, v2, uv = triple.u_square, triple.v_square, triple.u_dot_v
    best: Silver | None = None
    best_pair = (0, 0)
    best_determinant = 0
    retained = 0
    for p in range(-isqrt(8 * u2) - 2, isqrt(8 * u2) + 3):
        for q in range(-isqrt(8 * v2) - 2, isqrt(8 * v2) + 3):
            determinant = 8 * (u2 * v2 - uv * uv) - v2 * p * p - u2 * q * q + 2 * uv * p * q
            if determinant < 0:
                continue
            augmented = ((8, p, q), (p, u2, uv), (q, uv, v2))
            if any(
                0
                < sum(
                    coeff[i] * augmented[i][j] * coeff[j]
                    for i in range(3)
                    for j in range(3)
                )
                < MINIMUM
                for coeff in product(range(-3, 4), repeat=3)
                if coeff != (0, 0, 0)
            ):
                continue
            retained += 1
            value = Silver(p, q)
            if best is None or value > best:
                best = value
                best_pair = (p, q)
                best_determinant = determinant
    if best is None:
        raise AssertionError("the augmented pairing enumeration was empty")
    return PairingMaximum(best, best_pair, best_determinant, retained)


def lattice_pair(left: Sequence[int], right: Sequence[int]) -> int:
    matrix = gram()
    return sum(
        int(left[i]) * matrix[i][j] * int(right[j])
        for i in range(DIMENSION)
        for j in range(DIMENSION)
        if left[i] and right[j] and matrix[i][j]
    )


def explicit_root_audit(
    u: Sequence[int],
    v: Sequence[int],
    control: Silver,
    triple: GramTriple,
) -> bool:
    """Check the exact Gram, norm-two identity, integrality and primitivity."""

    observed = GramTriple(
        lattice_pair(u, u),
        lattice_pair(v, v),
        lattice_pair(u, v),
    )
    spatial_norm = Silver(
        observed.u_square + 2 * observed.v_square,
        2 * observed.u_dot_v,
    )
    spatial_pairings = tuple(
        Silver(
            sum(int(u[i]) * gram()[i][j] for i in range(DIMENSION)),
            sum(int(v[i]) * gram()[i][j] for i in range(DIMENSION)),
        )
        for j in range(DIMENSION)
    )
    return (
        observed == triple
        and spatial_norm == Silver(2) + RHO * control**2
        and all(value.is_integral for value in spatial_pairings)
        and any(int(value) % 2 for value in u)
    )


def explicit_root_pairing(
    left_u: Sequence[int],
    left_v: Sequence[int],
    left_control: Silver,
    right_u: Sequence[int],
    right_v: Sequence[int],
    right_control: Silver,
) -> Silver:
    """Lorentz pairing of two roots represented by ``(x,u+sqrt(2)v)``."""

    spatial = Silver(
        lattice_pair(left_u, right_u) + 2 * lattice_pair(left_v, right_v),
        lattice_pair(left_u, right_v) + lattice_pair(left_v, right_u),
    )
    return spatial - RHO * left_control * right_control


def directional_root_pairing(
    scale: Fraction,
    u: Sequence[int],
    v: Sequence[int],
    control: Silver,
) -> Silver:
    """Pair ``(1,scale*d)`` with ``(control,u+sqrt(2)v)``."""

    return Silver(scale) * Silver(
        lattice_pair(FIRST_CUTTING_D, u),
        lattice_pair(FIRST_CUTTING_D, v),
    ) - RHO * control


def second_cutting_root_audit() -> bool:
    return (
        explicit_root_audit(
            SECOND_CUTTING_U,
            SECOND_CUTTING_V,
            SECOND_CUTTING_CONTROL,
            SECOND_CUTTING_TRIPLE,
        )
        and tuple(
            SECOND_CUTTING_V[index] + 4 * FIRST_CUTTING_D[index]
            for index in range(DIMENSION)
        ) == SECOND_CUTTING_U
        and (
            lattice_pair(FIRST_CUTTING_D, SECOND_CUTTING_W),
            lattice_pair(SECOND_CUTTING_W, SECOND_CUTTING_W),
        ) == (1, 8)
    )


def third_candidate_root_audit() -> bool:
    return (
        explicit_root_audit(
            THIRD_CANDIDATE_U,
            THIRD_CANDIDATE_V,
            THIRD_CANDIDATE_CONTROL,
            THIRD_CANDIDATE_TRIPLE,
        )
        and tuple(
            THIRD_CANDIDATE_W[index] + 20 * FIRST_CUTTING_D[index]
            for index in range(DIMENSION)
        ) == THIRD_CANDIDATE_V
        and (
            lattice_pair(FIRST_CUTTING_D, THIRD_CANDIDATE_W),
            lattice_pair(THIRD_CANDIDATE_W, THIRD_CANDIDATE_W),
        ) == (-3, 12)
    )


def second_replacement_witness_norm() -> Silver:
    return Silver(8 * SECOND_REPLACEMENT_WITNESS_SCALE**2) - RHO


def third_diagnostic_witness_norm() -> Silver:
    return Silver(8 * THIRD_DIAGNOSTIC_WITNESS_SCALE**2) - RHO


def first_cutting_root_audit() -> bool:
    """Audit the explicit primitive q=2 wall realizing the first threat."""

    u2 = lattice_pair(FIRST_CUTTING_U, FIRST_CUTTING_U)
    v2 = lattice_pair(FIRST_CUTTING_V, FIRST_CUTTING_V)
    uv = lattice_pair(FIRST_CUTTING_U, FIRST_CUTTING_V)
    target = Silver(u2 + 2 * v2, 2 * uv)
    spatial_pairings = tuple(
        Silver(
            sum(FIRST_CUTTING_U[i] * gram()[i][j] for i in range(DIMENSION)),
            sum(FIRST_CUTTING_V[i] * gram()[i][j] for i in range(DIMENSION)),
        )
        for j in range(DIMENSION)
    )
    return (
        lattice_pair(FIRST_CUTTING_D, FIRST_CUTTING_D) == 8
        and lattice_pair(FIRST_CUTTING_D, FIRST_CUTTING_U) == 53
        and lattice_pair(FIRST_CUTTING_D, FIRST_CUTTING_V) == 37
        and (u2, v2, uv) == (
            FIRST_CUTTING_TRIPLE.u_square,
            FIRST_CUTTING_TRIPLE.v_square,
            FIRST_CUTTING_TRIPLE.u_dot_v,
        )
        and target == FIRST_CUTTING_NORM + RHO * FIRST_CUTTING_CONTROL**2
        and all(value.is_integral for value in spatial_pairings)
        and any(value % 2 for value in FIRST_CUTTING_U)
    )


def first_cutting_old_witness_pairing() -> Silver:
    """Pair ``(1,d/2)`` with the explicit first cutter."""

    return -RHO * FIRST_CUTTING_CONTROL + Silver(53, 37) / 2


def replacement_witness_norm() -> Silver:
    return Silver(8 * REPLACEMENT_WITNESS_SCALE**2) - RHO


def replacement_margin(control: Control, triple: GramTriple) -> Silver:
    """Positive means ``(1,(12/25)d)`` is strict against the type."""

    maximum = norm_eight_pairing_maximum(triple).maximum
    return RHO * control.control - Silver(REPLACEMENT_WITNESS_SCALE) * maximum


@dataclass(frozen=True, slots=True)
class EarlierAcceptanceSpec:
    height: Silver
    norm: Silver
    control: Silver
    triple: GramTriple
    reduced_gram: tuple[tuple[int, int], tuple[int, int]]
    first_coefficient: Silver
    second_coefficient: Silver
    pairing_modulus: int


EARLIER_ACCEPTANCE_SPECS = (
    EarlierAcceptanceSpec(Silver(93, 66), Silver(3), Silver(9, 6), GramTriple(138, 72, 99), ((12, 3), (3, 12)), Silver(-3, -2), Silver(-1, -1), 3),
    EarlierAcceptanceSpec(Silver(93, 66), Silver(6), Silver(12, 9), GramTriple(300, 132, 198), ((12, 6), (6, 36)), Silver(-1, -1), Silver(3, 2), 3),
    EarlierAcceptanceSpec(Silver(93, 66), Silver(6), Silver(12, 9), GramTriple(288, 138, 198), ((18, 0), (0, 30)), Silver(-1, -1), Silver(3, 2), 3),
    EarlierAcceptanceSpec(Silver(93, 66), Silver(6), Silver(12, 9), GramTriple(276, 144, 198), ((24, 6), (6, 24)), Silver(-3, -2), Silver(-1, -1), 3),
    EarlierAcceptanceSpec(Silver(93, 66), Silver(6), Silver(12, 9), GramTriple(264, 150, 198), ((18, -6), (-6, 24)), Silver(-4, -3), Silver(-1, -1), 3),
    EarlierAcceptanceSpec(Silver(93, 66), Silver(6), Silver(12, 9), GramTriple(252, 156, 198), ((12, 6), (6, 12)), Silver(-4, -3), Silver(-1, -1), 3),
    EarlierAcceptanceSpec(Silver(124, 88), Silver(2), Silver(8, 6), GramTriple(122, 64, 88), ((8, 4), (4, 10)), Silver(-4, -3), Silver(3, 2), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(188, 88, 128), ((8, 0), (0, 20)), Silver(-1, -1), Silver(3, 2), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(184, 90, 128), ((10, -2), (-2, 18)), Silver(-1, -1), Silver(3, 2), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(180, 92, 128), ((12, -4), (-4, 16)), Silver(-1, -1), Silver(3, 2), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(176, 94, 128), ((14, 6), (6, 14)), Silver(-3, -2), Silver(-1, -1), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(172, 96, 128), ((12, 4), (4, 12)), Silver(-1, -1), Silver(4, 3), 1),
    EarlierAcceptanceSpec(Silver(181, 128), Silver(2), Silver(10, 7), GramTriple(168, 98, 128), ((8, 0), (0, 10)), Silver(-1, -1), Silver(4, 3), 1),
    EarlierAcceptanceSpec(Silver(186, 132), Silver(6), Silver(18, 12), GramTriple(594, 264, 396), ((66, 0), (0, 0)), Silver(3, 2), Silver(0), 3),
    EarlierAcceptanceSpec(Silver(186, 132), Silver(6), Silver(18, 12), GramTriple(570, 276, 396), ((12, 0), (0, 42)), Silver(-4, -3), Silver(3, 2), 3),
    EarlierAcceptanceSpec(Silver(186, 132), Silver(6), Silver(18, 12), GramTriple(558, 282, 396), ((18, 0), (0, 30)), Silver(-4, -3), Silver(3, 2), 3),
    EarlierAcceptanceSpec(Silver(186, 132), Silver(6), Silver(18, 12), GramTriple(546, 288, 396), ((18, 0), (0, 24)), Silver(3, 2), Silver(-4, -3), 3),
)


@dataclass(frozen=True, slots=True)
class EarlierAcceptanceBound:
    feasible_cross_gram_count: int
    relaxed_maximum: Silver
    relaxed_maximizing_entries: tuple[int, int, int, int]
    positive_relaxed_count: int
    forbidden_positive_extensions: tuple[
        tuple[tuple[int, int, int, int], Silver, tuple[int, int, int, int], int], ...
    ]
    admitted_maximum: Silver
    admitted_maximizing_entries: tuple[int, int, int, int]


@lru_cache(maxsize=None)
def earlier_acceptance_bound(spec: EarlierAcceptanceSpec) -> EarlierAcceptanceBound:
    """Bound the cutter pairing against one earlier necessary binary type.

    The cutter plane has Gram ``[[14,4],[4,26]]``.  Every integral cross
    block satisfying the exact Schur-complement inequalities is enumerated;
    ternary types additionally force every cross entry to be divisible by 3.
    Harmless nonpositive blocks are retained as a superset.  Every positive
    block is rejected only after an explicit norm-below-eight combination is
    exhibited in its four-vector Gram matrix.
    """

    fixed_determinant = 348
    variable = spec.reduced_gram

    def quadratic_numerator(left: int, right: int) -> int:
        return 26 * left * left - 8 * left * right + 14 * right * right

    def columns(norm: int) -> tuple[tuple[int, int, int], ...]:
        left_bound = isqrt(14 * norm)
        right_bound = isqrt(26 * norm)
        modulus = spec.pairing_modulus
        left_start = -((-(-left_bound)) // modulus) * modulus
        right_start = -((-(-right_bound)) // modulus) * modulus
        return tuple(
            (left, right, numerator)
            for left in range(left_start, left_bound + 1, modulus)
            for right in range(right_start, right_bound + 1, modulus)
            if (numerator := quadratic_numerator(left, right))
            <= fixed_determinant * norm
        )

    first = columns(variable[0][0])
    second = columns(variable[1][1])
    fixed_first = Silver(-4, -3)
    fixed_second = Silver(3, 2)
    time_term = RHO * FIRST_CUTTING_CONTROL * spec.control
    weights = (
        fixed_first * spec.first_coefficient,
        fixed_first * spec.second_coefficient,
        fixed_second * spec.first_coefficient,
        fixed_second * spec.second_coefficient,
    )

    def pair_sign(pair: tuple[int, int]) -> int:
        rational, radical = pair
        if not rational:
            return (radical > 0) - (radical < 0)
        if not radical or (rational > 0) == (radical > 0):
            return (rational > 0) - (rational < 0)
        return (
            (rational > 0) - (rational < 0)
            if rational * rational > 2 * radical * radical
            else (radical > 0) - (radical < 0)
        )

    integer_weights = tuple((int(value.a), int(value.b)) for value in weights)
    time_pair = (int(time_term.a), int(time_term.b))
    feasible = 0
    relaxed_best_pair: tuple[int, int] | None = None
    relaxed_entries = (0, 0, 0, 0)
    admitted_best_pair: tuple[int, int] | None = None
    admitted_entries = (0, 0, 0, 0)
    positive = 0
    forbidden = []
    for p, r, first_numerator in first:
        residual_11 = fixed_determinant * variable[0][0] - first_numerator
        for q, s, second_numerator in second:
            residual_22 = fixed_determinant * variable[1][1] - second_numerator
            mixed_numerator = 26 * p * q - 4 * (p * s + r * q) + 14 * r * s
            residual_12 = fixed_determinant * variable[0][1] - mixed_numerator
            if residual_11 * residual_22 < residual_12 * residual_12:
                continue
            feasible += 1
            pairing_pair = (
                integer_weights[0][0] * p
                + integer_weights[1][0] * q
                + integer_weights[2][0] * r
                + integer_weights[3][0] * s
                - time_pair[0],
                integer_weights[0][1] * p
                + integer_weights[1][1] * q
                + integer_weights[2][1] * r
                + integer_weights[3][1] * s
                - time_pair[1],
            )
            entries = (p, q, r, s)
            if relaxed_best_pair is None or pair_sign(
                (pairing_pair[0] - relaxed_best_pair[0], pairing_pair[1] - relaxed_best_pair[1])
            ) > 0:
                relaxed_best_pair = pairing_pair
                relaxed_entries = entries
            if pair_sign(pairing_pair) > 0:
                positive += 1
                four_gram = (
                    (14, 4, p, q),
                    (4, 26, r, s),
                    (p, r, variable[0][0], variable[0][1]),
                    (q, s, variable[1][0], variable[1][1]),
                )
                certificate = next(
                    (
                        (coefficients, square)
                        for coefficients in product(range(-4, 5), repeat=4)
                        if coefficients != (0, 0, 0, 0)
                        and 0
                        < (square := sum(
                            coefficients[i] * four_gram[i][j] * coefficients[j]
                            for i in range(4)
                            for j in range(4)
                        ))
                        < MINIMUM
                    ),
                    None,
                )
                if certificate is None:
                    raise AssertionError("an unresolved positive earlier-wall extension appeared")
                forbidden.append((entries, Silver(*pairing_pair), certificate[0], certificate[1]))
                continue
            if admitted_best_pair is None or pair_sign(
                (pairing_pair[0] - admitted_best_pair[0], pairing_pair[1] - admitted_best_pair[1])
            ) > 0:
                admitted_best_pair = pairing_pair
                admitted_entries = entries
    if relaxed_best_pair is None or admitted_best_pair is None:
        raise AssertionError("the earlier-wall cross-Gram enumeration was empty")
    return EarlierAcceptanceBound(
        feasible,
        Silver(*relaxed_best_pair),
        relaxed_entries,
        positive,
        tuple(forbidden),
        Silver(*admitted_best_pair),
        admitted_entries,
    )


@dataclass(frozen=True, slots=True)
class FirstPrefixAudit:
    control_count: int
    height_count: int
    earlier_nonempty_controls: int
    endpoint_components: tuple[tuple[Silver, Silver, tuple[GramTriple, ...]], ...]
    maxima: tuple[PairingMaximum, ...]
    witness_scale: Silver
    witness_excess: Silver
    margins: tuple[Silver, ...]


@lru_cache(maxsize=1)
def first_prefix_audit() -> FirstPrefixAudit:
    rows = controls()
    nonempty = tuple((row, necessary_triples(row)) for row in rows if necessary_triples(row))
    earlier = tuple(row for row, _ in nonempty if row.height < FIRST_NECESSARY_HEIGHT)
    endpoint = tuple((row, triples) for row, triples in nonempty if row.height == FIRST_NECESSARY_HEIGHT)
    triples = tuple(triple for _, packet in endpoint for triple in packet)
    maxima = tuple(norm_eight_pairing_maximum(triple) for triple in triples)
    scale = Silver(1) / 2
    margins = tuple(
        RHO * row.control - scale * maximum.maximum
        for row, packet in endpoint
        for maximum in (norm_eight_pairing_maximum(triple) for triple in packet)
    )
    audit = FirstPrefixAudit(
        len(rows),
        len({row.height for row in rows}),
        len(earlier),
        tuple((row.norm_class.norm, row.control, packet) for row, packet in endpoint),
        maxima,
        scale,
        8 * scale * scale - RHO,
        margins,
    )
    if len(endpoint) != 2 or len(triples) != 6:
        raise AssertionError("the first necessary packet changed")
    if audit.earlier_nonempty_controls or audit.witness_excess.sign(0) <= 0:
        raise AssertionError("the first-prefix obstruction failed")
    if any(margin.sign(0) <= 0 for margin in margins):
        raise AssertionError("the fixed witness ceased to be strict")
    return audit


__all__ = [
    "COMPLETE_NORM_CLASS_PROOF",
    "DETERMINANT",
    "DIMENSION",
    "EARLIER_ACCEPTANCE_SPECS",
    "FIRST_NECESSARY_HEIGHT",
    "FIRST_CUTTING_CONTROL",
    "FIRST_CUTTING_D",
    "FIRST_CUTTING_HEIGHT",
    "FIRST_CUTTING_NORM",
    "FIRST_CUTTING_TRIPLE",
    "FIRST_CUTTING_U",
    "FIRST_CUTTING_V",
    "KISSING_NUMBER",
    "MINIMUM",
    "NORM_CLASSES",
    "RADICAL_THREE_INDEX",
    "RADICAL_THREE_KISSING_NUMBER",
    "RADICAL_THREE_MINIMUM",
    "RHO",
    "REPLACEMENT_WITNESS_SCALE",
    "SECOND_CUTTING_CONTROL",
    "SECOND_CUTTING_HEIGHT",
    "SECOND_CUTTING_TRIPLE",
    "SECOND_CUTTING_U",
    "SECOND_CUTTING_V",
    "SECOND_CUTTING_W",
    "SECOND_REPLACEMENT_WITNESS_SCALE",
    "SMITH_INVARIANTS",
    "SOURCE_URL",
    "THIRD_CANDIDATE_CONTROL",
    "THIRD_CANDIDATE_HEIGHT",
    "THIRD_CANDIDATE_TRIPLE",
    "THIRD_CANDIDATE_U",
    "THIRD_CANDIDATE_V",
    "THIRD_CANDIDATE_W",
    "THIRD_DIAGNOSTIC_WITNESS_SCALE",
    "catalogue_audit",
    "controls",
    "directional_root_pairing",
    "earlier_acceptance_bound",
    "explicit_root_audit",
    "explicit_root_pairing",
    "first_prefix_audit",
    "first_cutting_old_witness_pairing",
    "first_cutting_root_audit",
    "gram",
    "necessary_triples",
    "norm_eight_pairing_maximum",
    "replacement_margin",
    "replacement_witness_norm",
    "second_cutting_root_audit",
    "second_replacement_witness_norm",
    "third_candidate_root_audit",
    "third_diagnostic_witness_norm",
]
