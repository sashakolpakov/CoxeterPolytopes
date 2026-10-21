"""Exact first packet for the silver prime-seven ``Lambda25`` host.

The Nebe--Sloane catalogue prints the following basis at four times the
quadratic scale used for its integral Gram matrix.  Thus ``G=B*B^t/4``.
Keeping the basis as well as the derived Gram matrix makes the transcription
independently checkable and avoids copying the catalogue's line-wrapped Gram.

The arithmetic lattice is

``<-(2*sqrt(2)-1)> + Lambda25``

over ``Z[sqrt(2)]``.  The positive block has Smith form
``1,1,2^22,16``.  Its dyadic mirror valuations are therefore
``0,2,3,...,10``; the timelike coefficient supplies an independent optional
prime above seven.  This module certifies the complete Vinberg prefix through
the first positive packet at ``H=5+4*sqrt(2)``.

There is one height-zero wall.  The first positive packet is *not* the full
196,656-vector minimal shell: its norm-four reflection gate requires the
minimal vector to lie in the pairing radical modulo two, leaving 98,256
walls.  Every one of them is orthogonal to the height-zero wall, which gives
an exact strict hyperideal witness for this bounded prefix.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd, lcm, prod
from typing import TypeAlias

from sympy import Matrix
from sympy.matrices.normalforms import (
    hermite_normal_form,
    smith_normal_decomp,
    smith_normal_form,
)
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]

LAMBDA25_SOURCE_URL = (
    "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA25.html"
)

P7_LAMBDA25_COEFFICIENT = 2 * SQRT2 - SILVER_ONE
P7_GENERATOR = Silver(3, 1)
EPSILON = SILVER_ONE + SQRT2
FIRST_PACKET_HEIGHT = Silver(5, 4)
FIRST_PACKET_NORM = Silver(4)
FIRST_PACKET_CONTROL = Silver(2, 2)
FIRST_PACKET_TARGET = Silver(24, 16)
FIRST_PACKET_SHELL_NORM = 8
FIRST_PACKET_WALL_COUNT = 98_256
FULL_MINIMAL_SHELL_COUNT = 196_656

# This is the unique vector up to sign in the pairing radical modulo eight
# having square 16.  The chosen sign gives G*w0=-8*e_24.
HEIGHT_ZERO_VECTOR = (
    -2,
    1,
    0,
    1,
    2,
    1,
    1,
    -2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2,
)

# Exact minima of S_m={v in Lambda25 : (v,Lambda25) subset m*Z}.
# The audit script independently reconstructs every S_m and replays qfminim.
PAIRING_RADICAL_MINIMUM = {1: 8, 2: 8, 4: 16, 8: 16, 16: 64}
PAIRING_RADICAL_MINIMAL_COUNT = {
    1: 196_656,
    2: 98_256,
    4: 2,
    8: 2,
    16: 2,
}

LAMBDA25_BASE: IntegerMatrix = (
    (8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (2,2,2,2,0,0,0,0,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0),
    (2,2,0,0,2,2,0,0,2,2,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0),
    (2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0),
    (2,0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,0,0,0,0,0,0,0,0,0),
    (4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0),
    (2,0,2,0,2,0,0,2,2,2,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0),
    (2,0,0,2,2,2,0,0,2,0,2,0,0,0,0,0,2,0,2,0,0,0,0,0,0),
    (2,2,0,0,2,0,2,0,2,0,0,2,0,0,0,0,2,0,0,2,0,0,0,0,0),
    (0,2,2,2,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0,0,0),
    (0,0,0,0,0,0,0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,2,0,0,0),
    (0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0),
    (-3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0),
    (2,0,-2,0,2,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4),
)


def _bareiss_determinant(matrix: IntegerMatrix) -> int:
    work = [list(row) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    previous = 1
    sign = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
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


@lru_cache(maxsize=1)
def lambda25_gram() -> IntegerMatrix:
    gram = []
    for left in LAMBDA25_BASE:
        row = []
        for right in LAMBDA25_BASE:
            numerator = sum(a * b for a, b in zip(left, right))
            if numerator % 4:
                raise AssertionError("catalogue basis pairing is not divisible by four")
            row.append(numerator // 4)
        gram.append(tuple(row))
    return tuple(gram)


def lambda25_catalogue_audit() -> bool:
    gram = lambda25_gram()
    return (
        len(LAMBDA25_BASE) == 25
        and all(len(row) == 25 for row in LAMBDA25_BASE)
        and all(gram[i][j] == gram[j][i] for i in range(25) for j in range(25))
        and all(gram[i][i] % 2 == 0 for i in range(25))
        and min(gram[i][i] for i in range(25)) == 8
        and _bareiss_determinant(gram) == 2**26
    )


@lru_cache(maxsize=1)
def lambda25_smith_invariants() -> tuple[int, ...]:
    smith = smith_normal_form(Matrix(lambda25_gram()), domain=ZZ)
    return tuple(abs(int(smith[index, index])) for index in range(25))


@lru_cache(maxsize=None)
def lambda25_pairing_radical_basis(
    modulus: int,
) -> tuple[tuple[int, ...], ...]:
    """Canonical column basis of ``{v : G*v=0 mod modulus}``."""

    if modulus not in PAIRING_RADICAL_MINIMUM:
        raise ValueError("implemented Lambda25 radical moduli are 1,2,4,8,16")
    if modulus == 1:
        return tuple(tuple(int(i == j) for j in range(25)) for i in range(25))
    gram = Matrix(lambda25_gram())
    diagonal, _left, right = smith_normal_decomp(gram, domain=ZZ)
    factors = tuple(
        modulus // gcd(modulus, abs(int(diagonal[index, index])))
        for index in range(25)
    )
    basis = hermite_normal_form(right * Matrix.diag(*factors))
    if abs(int(basis.det())) != prod(factors):
        raise AssertionError("Lambda25 radical basis has the wrong index")
    if any(int(value) % modulus for value in gram * basis):
        raise AssertionError("Lambda25 radical basis violates its pairing gate")
    return tuple(
        tuple(int(basis[i, j]) for j in range(25)) for i in range(25)
    )


@lru_cache(maxsize=None)
def lambda25_pairing_radical_gram(
    modulus: int,
) -> IntegerMatrix:
    gram = Matrix(lambda25_gram())
    basis = Matrix(lambda25_pairing_radical_basis(modulus))
    radical = basis.T * gram * basis
    return tuple(
        tuple(int(radical[i, j]) for j in range(25)) for i in range(25)
    )


@dataclass(frozen=True, slots=True)
class P7Lambda25NormClass:
    norm: Silver
    dyadic_valuation: int
    p7_gate: bool

    @property
    def coefficient_pairing_moduli(self) -> tuple[int, int]:
        """Pairing moduli for ``y=u+sqrt(2)*v`` after the p7 factor."""

        exponent = max(0, self.dyadic_valuation - 2)
        return 2 ** ((exponent + 1) // 2), 2 ** (exponent // 2)


_DYADIC_CLASSES = tuple(
    (
        (
            Silver(2 ** (valuation // 2))
            if valuation % 2 == 0
            else Silver(2 ** ((valuation - 1) // 2)) * (Silver(2) + SQRT2)
        ),
        valuation,
    )
    for valuation in (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)
)


def _build_norm_classes() -> tuple[P7Lambda25NormClass, ...]:
    output = tuple(
        P7Lambda25NormClass(
            dyadic * (P7_GENERATOR if p7_gate else SILVER_ONE),
            valuation,
            p7_gate,
        )
        for p7_gate in (False, True)
        for dyadic, valuation in _DYADIC_CLASSES
    )
    if len(output) != 20 or len({row.norm for row in output}) != 20:
        raise AssertionError("the p7 Lambda25 norm alphabet has collisions")
    return output


P7_LAMBDA25_NORM_CLASSES = _build_norm_classes()

COMPLETE_P7_LAMBDA25_NORM_PROOF = (
    "Put pi=sqrt(2), rho=2*pi-1 and delta=(1+pi)*rho=3+pi. For a "
    "primitive root r=(x,y), reflection integrality says q(r)|2I(r), "
    "I(r)=(rho*x,(y,Lambda25)). The Smith form 1,1,2^22,16 has largest "
    "pi-primary exponent 8. Evenness excludes v_pi(q)=1. If v_pi(q)=e>2, "
    "then x is divisible by pi^(e-2) and G*y by pi^(e-2); e>10 would make "
    "both x and y divisible by pi. Thus the complete dyadic valuations are "
    "0 and 2 through 10. For y=u+pi*v, writing k=e-2 gives the exact gates "
    "G*u=0 mod 2^ceil(k/2) and G*v=0 mod 2^floor(k/2). Away from pi and "
    "the coefficient prime, local unimodularity excludes every prime. At "
    "delta the exponent is zero or one; in the latter case y lies in "
    "delta*Lambda25. The choices are independent. Class number one and "
    "U^+=U^2 give the twenty canonical generators obtained from "
    "1,2,4+2*pi,4,8+4*pi,8,16+8*pi,16,32+16*pi,32 with an optional "
    "factor delta."
)


@dataclass(frozen=True, slots=True, order=True)
class Lambda25GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int


def _gated_gram_triples_for(
    norm_class: P7Lambda25NormClass,
    target: Silver,
) -> tuple[Lambda25GramTriple, ...]:
    """Necessary and sufficient low-target binary Gram screen in this prefix."""

    if (
        not target.is_integral
        or target.sign(0) < 0
        or target.sign(1) < 0
        or target.b.denominator != 1
    ):
        return ()
    dot = target.b / 2
    if dot.denominator != 1:
        return ()
    rational_target = int(target.a)
    dot_int = int(dot)
    u_modulus, v_modulus = norm_class.coefficient_pairing_moduli
    cross_modulus = lcm(u_modulus, v_modulus)
    output = []
    for v_square in range(0, rational_target // 2 + 1, 2):
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if (
            (u_square and u_square < PAIRING_RADICAL_MINIMUM[u_modulus])
            or (v_square and v_square < PAIRING_RADICAL_MINIMUM[v_modulus])
            or u_square % u_modulus
            or v_square % v_modulus
            or dot_int % cross_modulus
            or dot_int * dot_int > u_square * v_square
        ):
            continue
        output.append(Lambda25GramTriple(u_square, v_square, dot_int))
    return tuple(output)


@dataclass(frozen=True, slots=True)
class P7Lambda25Control:
    norm_class: P7Lambda25NormClass
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_p7_lambda25_controls(
    height_bound: Silver = FIRST_PACKET_HEIGHT,
) -> tuple[P7Lambda25Control, ...]:
    """Enumerate every reflective arithmetic control through ``height_bound``."""

    if height_bound.sign(0) < 0:
        raise ValueError("height_bound must be nonnegative")
    levels: set[P7Lambda25Control] = set()
    for norm_class in P7_LAMBDA25_NORM_CLASSES:
        norm = norm_class.norm
        caps = (
            norm * height_bound / P7_LAMBDA25_COEFFICIENT,
            norm / (-P7_LAMBDA25_COEFFICIENT),
        )
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (
                2 * P7_LAMBDA25_COEFFICIENT * control / norm
            ).is_integral:
                continue
            target = norm + P7_LAMBDA25_COEFFICIENT * control * control
            if target.sign(1) < 0:
                continue
            height = P7_LAMBDA25_COEFFICIENT * control * control / norm
            if height > height_bound:
                continue
            levels.add(P7Lambda25Control(norm_class, height, control, target))
    return tuple(
        sorted(levels, key=lambda row: (row.height, row.norm_class.norm, row.control))
    )


def reduced_spatial_target(level: P7Lambda25Control) -> Silver:
    if level.norm_class.p7_gate:
        return level.target / (P7_GENERATOR * P7_GENERATOR)
    return level.target


def gated_gram_triples(level: P7Lambda25Control) -> tuple[Lambda25GramTriple, ...]:
    return _gated_gram_triples_for(level.norm_class, reduced_spatial_target(level))


def height_zero_gram_triples(
    norm_class: P7Lambda25NormClass,
) -> tuple[Lambda25GramTriple, ...]:
    target = norm_class.norm
    if norm_class.p7_gate:
        target /= P7_GENERATOR * P7_GENERATOR
    return _gated_gram_triples_for(norm_class, target)


def lambda25_pair(
    left: tuple[Silver | int, ...], right: tuple[Silver | int, ...]
) -> Silver:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("Lambda25 vectors have length 25")
    gram = lambda25_gram()
    return sum(
        (
            Silver(gram[i][j]) * Silver.coerce(left[i]) * Silver.coerce(right[j])
            for i in range(25)
            for j in range(25)
            if gram[i][j] and left[i] and right[j]
        ),
        start=SILVER_ZERO,
    )


def p7_lambda25_lorentz_pair(
    left: tuple[Silver | int, ...], right: tuple[Silver | int, ...]
) -> Silver:
    if len(left) != 26 or len(right) != 26:
        raise ValueError("the p7 Lambda25 host has rank 26")
    return (
        -P7_LAMBDA25_COEFFICIENT
        * Silver.coerce(left[0])
        * Silver.coerce(right[0])
        + lambda25_pair(tuple(left[1:]), tuple(right[1:]))
    )


def height_zero_vector_audit() -> bool:
    gram = lambda25_gram()
    image = tuple(
        sum(gram[i][j] * HEIGHT_ZERO_VECTOR[j] for j in range(25))
        for i in range(25)
    )
    square = sum(HEIGHT_ZERO_VECTOR[i] * image[i] for i in range(25))
    return image == (0,) * 24 + (-8,) and square == 16


def first_packet_root(b: tuple[int, ...]) -> tuple[Silver, ...]:
    if len(b) != 25 or lambda25_pair(b, b) != FIRST_PACKET_SHELL_NORM:
        raise ValueError("expected a Lambda25 minimal vector")
    return (FIRST_PACKET_CONTROL, *(EPSILON * coordinate for coordinate in b))


def first_packet_reflection_integral(b: tuple[int, ...]) -> bool:
    root = first_packet_root(b)
    gram = lambda25_gram()
    return all(sum(gram[i][j] * b[j] for j in range(25)) % 2 == 0 for i in range(25)) and all(
        (
            2
            * sum(
                (
                    Silver(gram[i][j]) * root[j + 1]
                    for j in range(25)
                    if gram[i][j]
                ),
                start=SILVER_ZERO,
            )
            / FIRST_PACKET_NORM
        ).is_integral
        for i in range(25)
    )


def first_packet_pairing_from_dot(dot: int) -> Silver:
    return (
        -P7_LAMBDA25_COEFFICIENT * FIRST_PACKET_CONTROL * FIRST_PACKET_CONTROL
        + EPSILON * EPSILON * dot
    )


def strict_hyperideal_witness() -> tuple[Silver, ...]:
    return (SILVER_ONE, *(-coordinate for coordinate in HEIGHT_ZERO_VECTOR))


def strict_witness_norm() -> Silver:
    return Silver(16) - P7_LAMBDA25_COEFFICIENT


def strict_witness_height_zero_pairing() -> Silver:
    return Silver(-16)


def strict_witness_first_packet_pairing() -> Silver:
    # Every reflecting minimal b is orthogonal to HEIGHT_ZERO_VECTOR.
    return -2 * P7_LAMBDA25_COEFFICIENT * EPSILON


__all__ = [
    "COMPLETE_P7_LAMBDA25_NORM_PROOF",
    "EPSILON",
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_SHELL_NORM",
    "FIRST_PACKET_TARGET",
    "FIRST_PACKET_WALL_COUNT",
    "FULL_MINIMAL_SHELL_COUNT",
    "HEIGHT_ZERO_VECTOR",
    "LAMBDA25_BASE",
    "LAMBDA25_SOURCE_URL",
    "Lambda25GramTriple",
    "P7_GENERATOR",
    "P7_LAMBDA25_COEFFICIENT",
    "P7_LAMBDA25_NORM_CLASSES",
    "PAIRING_RADICAL_MINIMAL_COUNT",
    "PAIRING_RADICAL_MINIMUM",
    "enumerate_p7_lambda25_controls",
    "first_packet_pairing_from_dot",
    "first_packet_reflection_integral",
    "first_packet_root",
    "gated_gram_triples",
    "height_zero_gram_triples",
    "height_zero_vector_audit",
    "lambda25_catalogue_audit",
    "lambda25_gram",
    "lambda25_pair",
    "lambda25_pairing_radical_basis",
    "lambda25_pairing_radical_gram",
    "lambda25_smith_invariants",
    "p7_lambda25_lorentz_pair",
    "reduced_spatial_target",
    "strict_hyperideal_witness",
    "strict_witness_first_packet_pairing",
    "strict_witness_height_zero_pairing",
    "strict_witness_norm",
]
