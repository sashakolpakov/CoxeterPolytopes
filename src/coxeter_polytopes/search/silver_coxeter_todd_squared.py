"""Exact first packet for a new rootless silver host in ``H^24``.

The host is

``<-(12*sqrt(2)-13)> + K12 + K12``

over ``Z[sqrt(2)]``, where ``K12`` is the Coxeter--Todd lattice.  Its
spatial block is even, has minimum four and determinant ``3^12``.  The
coefficient has norm ``-119=-7*17`` and is therefore in a coefficient
squareclass not used by the earlier silver searches in this repository.

The first part proves the complete Vinberg prefix through the first nonempty
height ``67+48*sqrt(2)``.  That packet consists of the norm-four and norm-eight
shells of ``K12^2`` and contains 613,872 walls.  The continuation helpers then
identify the first later height that cuts the initial norm-``4/3`` dual
direction and construct a smaller strict hyperideal replacement through
``393+278*sqrt(2)``.

Both conclusions are deliberately bounded.  The companion freezer certifies
infinite volume through the latter height, not Vinberg termination.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
)
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import (
    LeechGramTriple,
    candidate_leech_gram_triples,
)


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]
Vector: TypeAlias = tuple[Silver, ...]

K12_SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/K12.html"
K12_GRAM_STREAM_SHA256 = (
    "4e25ecfaf7ab6b1e0bb22d0390b1909a0792da7c3db05a9afedc1918a62e7aaa"
)

# The exact catalogue Gram matrix in its displayed integral basis.
K12_GRAM: IntegerMatrix = (
    (4, 0, 0, -2, 0, 0, 2, -1, -1, -1, 2, -1),
    (0, 4, 0, 0, -2, 0, 2, -1, -1, -1, -1, 2),
    (0, 0, 4, 0, 0, -2, 2, 2, 2, -1, -1, -1),
    (-2, 0, 0, 4, 0, 0, -1, -1, 2, 2, -1, -1),
    (0, -2, 0, 0, 4, 0, -1, 2, -1, 2, -1, -1),
    (0, 0, -2, 0, 0, 4, -1, -1, -1, 2, 2, 2),
    (2, 2, 2, -1, -1, -1, 4, 0, 0, -2, 0, 0),
    (-1, -1, 2, -1, 2, -1, 0, 4, 0, 0, -2, 0),
    (-1, -1, 2, 2, -1, -1, 0, 0, 4, 0, 0, -2),
    (-1, -1, -1, 2, 2, 2, -2, 0, 0, 4, 0, 0),
    (2, -1, -1, -1, -1, 2, 0, -2, 0, 0, 4, 0),
    (-1, 2, -1, -1, -1, 2, 0, 0, -2, 0, 0, 4),
)

K12_LEADING_PRINCIPAL_MINORS = (
    4,
    16,
    64,
    192,
    576,
    1728,
    1728,
    1728,
    1728,
    1296,
    972,
    729,
)
K12_SMITH_INVARIANTS = (1,) * 6 + (3,) * 6
K12_MINIMUM = 4
K12_RADICAL_THREE_MINIMUM = 12
K12_SHELL_COUNTS = {4: 756, 6: 4_032, 8: 20_412}

# A column basis of {v in K12 : (v,K12) is contained in 3Z}.  It has index
# 3^6.  Its Gram minimum is independently replayed with PARI by the freezer.
K12_RADICAL_THREE_BASIS: IntegerMatrix = (
    (3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
)

COEFFICIENT = Silver(-13, 12)
EPSILON = SILVER_ONE + SQRT2
EPSILON_SQUARED = Silver(3, 2)
P7_GENERATOR = Silver(3, 1)
P17_GENERATOR = Silver(5, -2)
FIRST_PACKET_HEIGHT = Silver(67, 48)
FIRST_PACKET_Q1_CONTROL = EPSILON_SQUARED
FIRST_PACKET_Q2_CONTROL = SQRT2 * EPSILON_SQUARED
FIRST_PACKET_WALL_COUNT = 613_872

# The first later height at which a wall can meet the strict first-packet
# witness.  The certificate for this statement is deliberately kept in this
# module: it uses only the complete control list, binary Gram data, and the
# four residue classes of the dual coordinate on K12.
FIRST_CUTTING_HEIGHT = Silver(393, 278)
FIRST_CUTTING_NORM = Silver(2)
FIRST_CUTTING_CONTROL = Silver(10, 7)
REPLACEMENT_WITNESS_SCALE = Fraction(7, 4)

COMPLETE_NORM_CLASS_PROOF = (
    "Put pi=sqrt(2), a=12*pi-13, delta7=3+pi and delta17=5-2*pi. "
    "Then a=(pi-1)*delta7*delta17. For a primitive root r=(x,y), "
    "q(r) divides 2I(r), I(r)=(a*x,(y,K12^2)). Away from pi, 3, "
    "delta7 and delta17 the spatial Gram is unimodular, so no prime can "
    "divide q. At each coefficient prime delta7,delta17 the exponent is "
    "zero or one, and exponent one forces y into that prime times K12^2. "
    "At 3 the elementary Smith form 1^12,3^12 gives exponent zero or one; "
    "exponent one forces every spatial pairing to be divisible by 3. At "
    "pi the spatial lattice is even and unimodular, so the only exponents "
    "are zero and two. The four local choices are independent. Since "
    "Z[sqrt(2)] has class number one and every totally-positive unit is a "
    "square, the sixteen products of 1 or 2, optionally 3, delta7 and "
    "delta17, are the complete canonical primitive norm classes."
)


def _bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(map(int, row)) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("a determinant matrix must be square")
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
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
        pivot = next(
            (index for index in range(row, len(work)) if work[index][column]),
            None,
        )
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


def _matmul_integer(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> IntegerMatrix:
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def _transpose(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))


@lru_cache(maxsize=1)
def k12_radical_three_gram() -> IntegerMatrix:
    basis = K12_RADICAL_THREE_BASIS
    return _matmul_integer(_matmul_integer(_transpose(basis), K12_GRAM), basis)


@lru_cache(maxsize=1)
def k12_basis_audit() -> bool:
    inverse = _fraction_inverse(
        tuple(tuple(Fraction(value) for value in row) for row in K12_GRAM)
    )
    leading = tuple(
        _bareiss_determinant(tuple(row[:size] for row in K12_GRAM[:size]))
        for size in range(1, 13)
    )
    radical = k12_radical_three_gram()
    return (
        all(K12_GRAM[i][j] == K12_GRAM[j][i] for i in range(12) for j in range(12))
        and all(K12_GRAM[i][i] % 2 == 0 for i in range(12))
        and leading == K12_LEADING_PRINCIPAL_MINORS
        and _bareiss_determinant(K12_GRAM) == 729
        and _rank_mod_prime(K12_GRAM, 3) == 6
        and all((3 * value).denominator == 1 for row in inverse for value in row)
        and _bareiss_determinant(K12_RADICAL_THREE_BASIS) == 729
        and all(value % 3 == 0 for row in radical for value in row)
    )


def coefficient_factorization_audit() -> bool:
    return (
        COEFFICIENT.sign(0) > 0
        and COEFFICIENT.sign(1) < 0
        and COEFFICIENT.norm() == -119
        and COEFFICIENT == (SQRT2 - 1) * P7_GENERATOR * P17_GENERATOR
        and P7_GENERATOR.norm() == 7
        and P17_GENERATOR.norm() == 17
        and all(P7_GENERATOR.sign(place) > 0 for place in (0, 1))
        and all(P17_GENERATOR.sign(place) > 0 for place in (0, 1))
    )


@dataclass(frozen=True, slots=True)
class CoxeterToddNormClass:
    norm: Silver
    dyadic: bool
    ternary: bool
    p7: bool
    p17: bool

    @property
    def odd_spatial_divisor(self) -> Silver:
        return (
            (P7_GENERATOR if self.p7 else SILVER_ONE)
            * (P17_GENERATOR if self.p17 else SILVER_ONE)
        )


def _build_norm_classes() -> tuple[CoxeterToddNormClass, ...]:
    output = []
    for dyadic, ternary, p7, p17 in product((False, True), repeat=4):
        norm = Silver(2 if dyadic else 1)
        if ternary:
            norm *= 3
        if p7:
            norm *= P7_GENERATOR
        if p17:
            norm *= P17_GENERATOR
        output.append(CoxeterToddNormClass(norm, dyadic, ternary, p7, p17))
    if len({row.norm for row in output}) != 16:
        raise AssertionError("the Coxeter--Todd norm alphabet has collisions")
    return tuple(output)


NORM_CLASSES = _build_norm_classes()


@dataclass(frozen=True, slots=True)
class CoxeterToddControl:
    norm_class: CoxeterToddNormClass
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_controls(
    height_bound: Silver = FIRST_PACKET_HEIGHT,
) -> tuple[CoxeterToddControl, ...]:
    if height_bound.sign(0) < 0:
        raise ValueError("height_bound must be nonnegative")
    output: set[CoxeterToddControl] = set()
    for norm_class in NORM_CLASSES:
        norm = norm_class.norm
        caps = (norm * height_bound / COEFFICIENT, norm / (-COEFFICIENT))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * COEFFICIENT * control / norm).is_integral:
                continue
            target = norm + COEFFICIENT * control * control
            if target.sign(1) < 0:
                continue
            height = COEFFICIENT * control * control / norm
            if height > height_bound:
                continue
            output.add(CoxeterToddControl(norm_class, height, control, target))
    return tuple(
        sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control))
    )


def reduced_target(control: CoxeterToddControl) -> Silver:
    divisor = control.norm_class.odd_spatial_divisor
    return control.target / (divisor * divisor)


def _has_combination_below(triple: LeechGramTriple, minimum: int) -> bool:
    s, t, dot = triple.u_square, triple.v_square, triple.u_dot_v
    determinant = s * t - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(s, t), abs(dot)) < minimum
    # Completing the square gives q(i,j)>=D*i^2/t and >=D*j^2/s.
    i_bound = isqrt(((minimum - 1) * t) // determinant) + 1
    j_bound = isqrt(((minimum - 1) * s) // determinant) + 1
    return any(
        0 < i * i * s + j * j * t + 2 * i * j * dot < minimum
        for i in range(-i_bound, i_bound + 1)
        for j in range(-j_bound, j_bound + 1)
        if i or j
    )


def gated_gram_triples(
    control: CoxeterToddControl,
) -> tuple[tuple[int, int, int], ...]:
    target = reduced_target(control)
    minimum = K12_RADICAL_THREE_MINIMUM if control.norm_class.ternary else K12_MINIMUM
    output = []
    for triple in candidate_leech_gram_triples(target):
        if control.norm_class.ternary and any(
            value % 3
            for value in (triple.u_square, triple.v_square, triple.u_dot_v)
        ):
            continue
        if _has_combination_below(triple, minimum):
            continue
        output.append((triple.u_square, triple.v_square, triple.u_dot_v))
    return tuple(output)


def height_zero_gated_gram_triples(
    norm_class: CoxeterToddNormClass,
) -> tuple[tuple[int, int, int], ...]:
    control = CoxeterToddControl(norm_class, SILVER_ZERO, SILVER_ZERO, norm_class.norm)
    return gated_gram_triples(control)


@dataclass(frozen=True, slots=True)
class FirstPacketComponent:
    norm: Silver
    control: Silver
    target: Silver
    raw_gram_triples: tuple[tuple[int, int, int], ...]
    imprimitive_gram_triples: tuple[tuple[int, int, int], ...]
    impossible_gram_triples: tuple[tuple[int, int, int], ...]
    primitive_gram_triples: tuple[tuple[int, int, int], ...]
    spatial_shell_norm: int
    wall_count: int


def first_packet_components() -> tuple[FirstPacketComponent, ...]:
    rows = tuple(row for row in enumerate_controls() if row.height == FIRST_PACKET_HEIGHT)
    if len(rows) != 2:
        raise AssertionError("the first simultaneous control packet changed")
    by_norm = {row.norm_class.norm: row for row in rows}
    q1 = by_norm[Silver(1)]
    q2 = by_norm[Silver(2)]
    if gated_gram_triples(q1) != ((36, 16, 24),):
        raise AssertionError("the norm-one packet component changed")
    if gated_gram_triples(q2) != ((72, 32, 48), (64, 36, 48)):
        raise AssertionError("the norm-two packet component changed")
    return (
        FirstPacketComponent(
            Silver(1),
            q1.control,
            q1.target,
            ((36, 16, 24),),
            (),
            (),
            ((36, 16, 24),),
            4,
            2 * K12_SHELL_COUNTS[4],
        ),
        FirstPacketComponent(
            Silver(2),
            q2.control,
            q2.target,
            ((64, 36, 48), (68, 34, 48), (72, 32, 48)),
            ((64, 36, 48),),
            ((68, 34, 48),),
            ((72, 32, 48),),
            8,
            2 * K12_SHELL_COUNTS[8] + K12_SHELL_COUNTS[4] ** 2,
        ),
    )


def first_packet_chronology_audit() -> bool:
    controls = enumerate_controls()
    components = first_packet_components()
    return (
        len(controls) == 50
        and len({row.height for row in controls}) == 35
        and sum(row.height < FIRST_PACKET_HEIGHT for row in controls) == 48
        and all(not gated_gram_triples(row) for row in controls if row.height < FIRST_PACKET_HEIGHT)
        and sum(component.wall_count for component in components) == FIRST_PACKET_WALL_COUNT
        and all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES)
    )


def k12_pair(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    if len(left) != 12 or len(right) != 12:
        raise ValueError("K12 coordinate vectors have length twelve")
    return sum(
        Fraction(left[i]) * K12_GRAM[i][j] * Fraction(right[j])
        for i in range(12)
        for j in range(12)
        if K12_GRAM[i][j]
    )


def k12_squared_pair_silver(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    if len(left) != 24 or len(right) != 24:
        raise ValueError("K12^2 coordinate vectors have length twenty-four")
    return sum(
        (
            left[offset + i] * K12_GRAM[i][j] * right[offset + j]
            for offset in (0, 12)
            for i in range(12)
            for j in range(12)
            if K12_GRAM[i][j]
        ),
        start=SILVER_ZERO,
    )


def lorentz_pair(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("host vectors have length twenty-five")
    return -COEFFICIENT * left[0] * right[0] + k12_squared_pair_silver(
        left[1:], right[1:]
    )


def first_packet_root(shell_norm: int, spatial: Sequence[int]) -> Vector:
    if len(spatial) != 24:
        raise ValueError("a K12^2 shell vector has length twenty-four")
    square = k12_squared_pair_silver(
        tuple(Silver(value) for value in spatial),
        tuple(Silver(value) for value in spatial),
    )
    if square != shell_norm or shell_norm not in (4, 8):
        raise ValueError("the spatial vector is not in a first-packet shell")
    if shell_norm == 4:
        time = EPSILON_SQUARED
    else:
        time = SQRT2 * EPSILON_SQUARED
    return (time, *(EPSILON_SQUARED * value for value in spatial))


def first_packet_root_is_reflective(root: Sequence[Silver]) -> bool:
    norm = lorentz_pair(root, root)
    if norm not in (Silver(1), Silver(2)):
        return False
    time_integral = (2 * COEFFICIENT * root[0] / norm).is_integral
    spatial = root[1:]
    pairings = tuple(
        sum(
            (
                spatial[offset + i] * K12_GRAM[i][j]
                for i in range(12)
                if K12_GRAM[i][j]
            ),
            start=SILVER_ZERO,
        )
        for offset in (0, 12)
        for j in range(12)
    )
    return time_integral and all((2 * value / norm).is_integral for value in pairings)


def same_height_acuteness_bounds() -> tuple[tuple[str, Silver], ...]:
    # Minimum four gives d.e<=2 for distinct norm-four vectors,
    # f.g<=6 for distinct norm-eight vectors, and Cauchy gives d.f<=5.
    return (
        ("q1-q1", EPSILON_SQUARED**2 * (Silver(2) - COEFFICIENT)),
        ("q1-q2", EPSILON_SQUARED**2 * (Silver(5) - SQRT2 * COEFFICIENT)),
        ("q2-q2", EPSILON_SQUARED**2 * (Silver(6) - 2 * COEFFICIENT)),
    )


# h=G^-1 e_1 in the first K12 factor.  It belongs to K12^*, has h^2=4/3,
# and h.d is the first integral coordinate of every d in K12.
DUAL_DIRECTION: tuple[Fraction, ...] = (
    Fraction(4, 3),
    Fraction(0),
    Fraction(0),
    Fraction(2, 3),
    Fraction(0),
    Fraction(0),
    Fraction(-2, 3),
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(-1, 3),
    Fraction(-1, 3),
    Fraction(2, 3),
)
WITNESS_SCALE = Fraction(9, 5)


def dual_direction_audit() -> bool:
    products = tuple(
        sum(Fraction(K12_GRAM[i][j]) * DUAL_DIRECTION[j] for j in range(12))
        for i in range(12)
    )
    return products == (Fraction(1),) + (Fraction(0),) * 11 and k12_pair(
        DUAL_DIRECTION, DUAL_DIRECTION
    ) == Fraction(4, 3)


def strict_hyperideal_witness() -> Vector:
    return (
        SILVER_ONE,
        *(Silver(WITNESS_SCALE * value) for value in DUAL_DIRECTION),
        *(SILVER_ZERO for _ in range(12)),
    )


def strict_witness_norm() -> Silver:
    return Silver(WITNESS_SCALE * WITNESS_SCALE * Fraction(4, 3)) - COEFFICIENT


def strict_witness_pairing_bounds() -> tuple[tuple[str, int, Silver], ...]:
    # h.d is integral.  Cauchy gives |h.d|<=2 on d^2=4 and <=3 on
    # f^2<=8.  These bounds cover both factors (the second pairs to zero).
    q1 = EPSILON_SQUARED * (-COEFFICIENT + 2 * WITNESS_SCALE)
    q2 = EPSILON_SQUARED * (-SQRT2 * COEFFICIENT + 3 * WITNESS_SCALE)
    return (("q1 norm-four", 2, q1), ("q2 norm-eight", 3, q2))


def strict_hyperideal_audit() -> bool:
    witness = strict_hyperideal_witness()
    return (
        dual_direction_audit()
        and lorentz_pair(witness, witness) == strict_witness_norm()
        and strict_witness_norm() == Silver(Fraction(433, 25), -12)
        and strict_witness_norm().sign(0) > 0
        and all(bound.sign(0) < 0 for _name, _dot, bound in strict_witness_pairing_bounds())
    )


# ---------------------------------------------------------------------------
# Continuation through the first packet which cuts the 9/5 dual witness.
# ---------------------------------------------------------------------------


_DUAL_COSET_RESIDUAL = {
    0: Fraction(0),
    1: Fraction(13, 4),
    2: Fraction(1),
    3: Fraction(13, 4),
}


def dual_coordinate_coset_minimum(coordinate: int) -> Fraction:
    """Exact minimum in K12 with fixed pairing against ``h``.

    Since ``h=G^-1 e1``, this pairing is the first displayed coordinate.
    The vector ``3h`` is integral, has first coordinate four and norm twelve,
    so the orthogonal residual depends only on the coordinate modulo four.
    The four base residues have minima ``0,4,4,10``.  Consequently

    ``min{x^2:x_1=k} = 3*k^2/4 + rho(k mod 4)``.

    The freezer independently replays the four base residues with PARI.
    """

    coordinate = int(coordinate)
    return Fraction(3 * coordinate * coordinate, 4) + _DUAL_COSET_RESIDUAL[
        coordinate % 4
    ]


@lru_cache(maxsize=None)
def dual_direction_pairing_candidates(
    triple: tuple[int, int, int],
) -> tuple[tuple[int, int], ...]:
    """A finite necessary superset for ``(h.u,h.v)``.

    The 3-by-3 Gram matrix of ``h,u,v`` must be positive semidefinite.  Every
    realization must also obey the exact fixed-coordinate minimum for each of
    the eight combinations ``i*u+j*v`` with ``i,j`` in ``{-1,0,1}``.
    Checking only these eight combinations is intentional: it remains a
    rigorous necessary gate, and is sufficient for all maxima through
    :data:`FIRST_CUTTING_HEIGHT`.
    """

    s, t, dot = map(int, triple)
    a_bound = isqrt((4 * s) // 3) + 2
    b_bound = isqrt((4 * t) // 3) + 2
    output = []
    for a_pair in range(-a_bound, a_bound + 1):
        s4 = 4 * s - 3 * a_pair * a_pair
        if s4 < 0:
            continue
        for b_pair in range(-b_bound, b_bound + 1):
            t4 = 4 * t - 3 * b_pair * b_pair
            cross4 = 4 * dot - 3 * a_pair * b_pair
            if t4 < 0 or s4 * t4 < cross4 * cross4:
                continue
            if any(
                Fraction(i * i * s + j * j * t + 2 * i * j * dot)
                < dual_coordinate_coset_minimum(i * a_pair + j * b_pair)
                for i in (-1, 0, 1)
                for j in (-1, 0, 1)
                if i or j
            ):
                continue
            output.append((a_pair, b_pair))
    return tuple(output)


@dataclass(frozen=True, slots=True)
class DualDirectionLevelBound:
    height: Silver
    norm: Silver
    control: Silver
    gated_gram_type_count: int
    directional_pair_count: int
    maximum_spatial_pairing: Silver
    maximizing_gram_type: tuple[int, int, int]
    maximizing_pair: tuple[int, int]


def dual_direction_level_bound(
    control: CoxeterToddControl,
) -> DualDirectionLevelBound:
    """Maximize pairing with ``h`` over a rigorous necessary Gram superset."""

    triples = gated_gram_triples(control)
    rows = tuple(
        (
            control.norm_class.odd_spatial_divisor * Silver(a_pair, b_pair),
            triple,
            (a_pair, b_pair),
        )
        for triple in triples
        for a_pair, b_pair in dual_direction_pairing_candidates(triple)
    )
    if not rows:
        raise ValueError("the control has no dual-direction Gram candidate")
    maximum, triple, pair = max(rows, key=lambda row: row[0])
    return DualDirectionLevelBound(
        control.height,
        control.norm_class.norm,
        control.control,
        len(triples),
        len(rows),
        maximum,
        triple,
        pair,
    )


def transformed_difference_gram(
    triple: tuple[int, int, int],
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Gram of ``d=u-v, e=2v-u`` for a binary coefficient Gram."""

    s, t, dot = map(int, triple)
    d_square = s + t - 2 * dot
    e_square = s + 4 * t - 4 * dot
    d_dot_e = 3 * dot - s - 2 * t
    return ((d_square, d_dot_e), (d_dot_e, e_square))


# A particularly small realization of the H=393+278*sqrt(2), q=2 type
# (380,204,278).  Here D=U-V and E=2V-U; both lie in the first K12 factor.
FIRST_CUTTING_D: tuple[int, ...] = (
    6,
    0,
    0,
    2,
    0,
    0,
    -3,
    1,
    2,
    -1,
    -2,
    3,
)
FIRST_CUTTING_E: tuple[int, ...] = (
    10,
    1,
    0,
    4,
    -1,
    0,
    -6,
    3,
    2,
    -2,
    -3,
    3,
)
FIRST_CUTTING_U: tuple[int, ...] = tuple(
    2 * d + e for d, e in zip(FIRST_CUTTING_D, FIRST_CUTTING_E, strict=True)
)
FIRST_CUTTING_V: tuple[int, ...] = tuple(
    d + e for d, e in zip(FIRST_CUTTING_D, FIRST_CUTTING_E, strict=True)
)


def first_cutting_root() -> Vector:
    """An explicit primitive accepted root at the first cutting height."""

    spatial = tuple(
        Silver(u, v)
        for u, v in zip(FIRST_CUTTING_U, FIRST_CUTTING_V, strict=True)
    ) + (SILVER_ZERO,) * 12
    return (FIRST_CUTTING_CONTROL, *spatial)


def first_cutting_root_audit() -> bool:
    root = first_cutting_root()
    gram = (
        int(k12_pair(FIRST_CUTTING_U, FIRST_CUTTING_U)),
        int(k12_pair(FIRST_CUTTING_V, FIRST_CUTTING_V)),
        int(k12_pair(FIRST_CUTTING_U, FIRST_CUTTING_V)),
    )
    spatial_pairings = tuple(
        sum(
            (
                root[1 + i] * K12_GRAM[i][j]
                for i in range(12)
                if K12_GRAM[i][j]
            ),
            start=SILVER_ZERO,
        )
        for j in range(12)
    )
    return (
        k12_pair(FIRST_CUTTING_D, FIRST_CUTTING_D) == 28
        and k12_pair(FIRST_CUTTING_E, FIRST_CUTTING_E) == 84
        and k12_pair(FIRST_CUTTING_D, FIRST_CUTTING_E) == 46
        and gram == (380, 204, 278)
        and FIRST_CUTTING_U[0] == 22
        and FIRST_CUTTING_V[0] == 16
        and lorentz_pair(root, root) == FIRST_CUTTING_NORM
        and (COEFFICIENT * FIRST_CUTTING_CONTROL).is_integral
        and all(value.is_integral for value in spatial_pairings)
        # Divisibility by sqrt(2) would force every rational coefficient even.
        and any(value % 2 for value in FIRST_CUTTING_U)
    )


@dataclass(frozen=True, slots=True)
class CrossGramMaximum:
    variable_gram: tuple[tuple[int, int], tuple[int, int]]
    first_column_count: int
    second_column_count: int
    feasible_cross_gram_count: int
    maximum: Silver
    maximizing_cross_entries: tuple[int, int, int, int]


@lru_cache(maxsize=None)
def cross_gram_maximum_against_first_cutter(
    variable_gram: tuple[tuple[int, int], tuple[int, int]],
) -> CrossGramMaximum:
    """Exact PSD maximum against the cutter's ``(D,E)`` coefficient plane.

    If the cross entries are ``A=D.d, B=D.e, C=E.d, D=E.e``, the full
    four-vector Gram matrix must be positive semidefinite.  The fixed block is
    ``[[28,46],[46,84]]`` with determinant 236.  Integer Schur complements
    enumerate every necessary cross block without floating point arithmetic.
    """

    if (
        len(variable_gram) != 2
        or any(len(row) != 2 for row in variable_gram)
        or variable_gram[0][1] != variable_gram[1][0]
    ):
        raise ValueError("variable_gram must be a symmetric 2-by-2 matrix")
    n1 = int(variable_gram[0][0])
    dot = int(variable_gram[0][1])
    n2 = int(variable_gram[1][1])
    determinant = 236

    def quadratic_numerator(left: int, right: int) -> int:
        return 84 * left * left - 92 * left * right + 28 * right * right

    def mixed_numerator(a: int, c: int, b: int, d: int) -> int:
        return 84 * a * b - 46 * (a * d + c * b) + 28 * c * d

    def columns(norm: int) -> tuple[tuple[int, int, int], ...]:
        a_bound = isqrt(28 * norm) + 2
        c_bound = isqrt(84 * norm) + 2
        return tuple(
            (left, right, numerator)
            for left in range(-a_bound, a_bound + 1)
            for right in range(-c_bound, c_bound + 1)
            if (numerator := quadratic_numerator(left, right))
            <= norm * determinant
        )

    first = columns(n1)
    second = columns(n2)
    feasible_count = 0
    best: Silver | None = None
    best_entries: tuple[int, int, int, int] | None = None
    for a_pair, c_pair, q1 in first:
        residual_11 = n1 * determinant - q1
        for b_pair, d_pair, q2 in second:
            residual_22 = n2 * determinant - q2
            residual_12 = dot * determinant - mixed_numerator(
                a_pair, c_pair, b_pair, d_pair
            )
            if residual_11 * residual_22 < residual_12 * residual_12:
                continue
            feasible_count += 1
            value = Silver(
                6 * a_pair + 4 * (b_pair + c_pair) + 3 * d_pair,
                4 * a_pair + 3 * (b_pair + c_pair) + 2 * d_pair,
            )
            if best is None or value > best:
                best = value
                best_entries = (a_pair, b_pair, c_pair, d_pair)
    if best is None or best_entries is None:
        raise AssertionError("the PSD cross-Gram enumeration was empty")
    return CrossGramMaximum(
        variable_gram,
        len(first),
        len(second),
        feasible_count,
        best,
        best_entries,
    )


def replacement_hyperideal_witness() -> Vector:
    return (
        SILVER_ONE,
        *(Silver(REPLACEMENT_WITNESS_SCALE * value) for value in DUAL_DIRECTION),
        *(SILVER_ZERO for _ in range(12)),
    )


def replacement_witness_norm() -> Silver:
    return Silver(
        REPLACEMENT_WITNESS_SCALE
        * REPLACEMENT_WITNESS_SCALE
        * Fraction(4, 3)
    ) - COEFFICIENT


__all__ = [
    "COEFFICIENT",
    "COMPLETE_NORM_CLASS_PROOF",
    "CrossGramMaximum",
    "DUAL_DIRECTION",
    "DualDirectionLevelBound",
    "FIRST_CUTTING_CONTROL",
    "FIRST_CUTTING_D",
    "FIRST_CUTTING_E",
    "FIRST_CUTTING_HEIGHT",
    "FIRST_CUTTING_NORM",
    "FIRST_CUTTING_U",
    "FIRST_CUTTING_V",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_WALL_COUNT",
    "K12_GRAM",
    "K12_GRAM_STREAM_SHA256",
    "K12_RADICAL_THREE_MINIMUM",
    "K12_SHELL_COUNTS",
    "K12_SMITH_INVARIANTS",
    "K12_SOURCE_URL",
    "NORM_CLASSES",
    "P17_GENERATOR",
    "P7_GENERATOR",
    "REPLACEMENT_WITNESS_SCALE",
    "WITNESS_SCALE",
    "coefficient_factorization_audit",
    "cross_gram_maximum_against_first_cutter",
    "dual_coordinate_coset_minimum",
    "dual_direction_audit",
    "dual_direction_level_bound",
    "dual_direction_pairing_candidates",
    "enumerate_controls",
    "first_cutting_root",
    "first_cutting_root_audit",
    "first_packet_chronology_audit",
    "first_packet_components",
    "first_packet_root",
    "first_packet_root_is_reflective",
    "gated_gram_triples",
    "height_zero_gated_gram_triples",
    "k12_basis_audit",
    "k12_radical_three_gram",
    "lorentz_pair",
    "replacement_hyperideal_witness",
    "replacement_witness_norm",
    "same_height_acuteness_bounds",
    "strict_hyperideal_audit",
    "strict_hyperideal_witness",
    "strict_witness_norm",
    "strict_witness_pairing_bounds",
    "transformed_difference_gram",
]
