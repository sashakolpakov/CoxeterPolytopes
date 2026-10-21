"""Exact first Vinberg packet for the silver prime-seven ``Lambda_22`` host.

The lattice is

``<-(2*sqrt(2)-1)> + Lambda_22``

over ``Z[sqrt(2)]``.  The positive block is the authoritative integral
Nebe--Sloane Gram matrix.  This module certifies only the complete prefix
through the first nonempty packet at ``H=5+4*sqrt(2)``.

The discriminant of ``Lambda_22`` contributes the dyadic valuations
``0,2,3,4`` and an optional inert-prime-three factor.  The timelike
coefficient contributes an independent optional prime above seven.  Thus
there are sixteen primitive mirror-norm classes.  The prime-seven classes
have spatial component divisible by ``3+sqrt(2)``; the remaining local
pairing gates are applied after that factor is removed.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import lcm
from typing import Sequence

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.cubic7_lambda22 import (
    GRAM_FILE_SHA256,
    GRAM_STREAM_SHA256,
    SOURCE_HTML_SHA256,
    SOURCE_URL,
    lambda22_basis_audit,
    lambda22_gram,
    lambda22_gram_file_sha256,
    lambda22_gram_stream_sha256,
    lambda22_smith_invariants,
)
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import has_short_binary_combination


P7_LAMBDA22_COEFFICIENT = 2 * SQRT2 - SILVER_ONE
P7_GENERATOR = Silver(3, 1)
EPSILON = SILVER_ONE + SQRT2
FIRST_PACKET_HEIGHT = Silver(5, 4)
FIRST_PACKET_NORM = Silver(2)
FIRST_PACKET_CONTROL = Silver(2, 1)
FIRST_PACKET_TARGET = Silver(12, 8)
FIRST_PACKET_SHELL_NORM = 4
FIRST_PACKET_WALL_COUNT = 49_896

# The exact dual direction used by the strict hyperideal witness.  It obeys
# G*k=3*e_22 and k^2=24, so h=k/6 has square 2/3 and h.d=d_22/2.
P7_DUAL_DIRECTION = (
    -4490,
    27287,
    17909,
    13568,
    8412,
    18464,
    12340,
    6324,
    3634,
    1942,
    2575,
    1337,
    544,
    272,
    408,
    214,
    32,
    64,
    38,
    22,
    12,
    8,
)

PAIRING_RADICAL_MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}


@dataclass(frozen=True, slots=True)
class P7Lambda22NormClass:
    norm: Silver
    dyadic_valuation: int
    ternary_gate: bool
    p7_gate: bool

    @property
    def coefficient_pairing_moduli(self) -> tuple[int, int]:
        """Pairing moduli for ``z=u+sqrt(2)*v`` after the p7 factor."""

        u_modulus = v_modulus = 1
        if self.dyadic_valuation == 3:
            u_modulus = 2
        elif self.dyadic_valuation == 4:
            u_modulus = v_modulus = 2
        if self.ternary_gate:
            u_modulus = lcm(u_modulus, 3)
            v_modulus = lcm(v_modulus, 3)
        return u_modulus, v_modulus


_DYADIC_CLASSES = (
    (Silver(1), 0),
    (Silver(2), 2),
    (Silver(4, 2), 3),
    (Silver(4), 4),
)


def _build_norm_classes() -> tuple[P7Lambda22NormClass, ...]:
    output = []
    for p7_gate in (False, True):
        for ternary_gate in (False, True):
            for dyadic, valuation in _DYADIC_CLASSES:
                norm = dyadic
                if ternary_gate:
                    norm *= 3
                if p7_gate:
                    norm *= P7_GENERATOR
                output.append(
                    P7Lambda22NormClass(
                        norm,
                        valuation,
                        ternary_gate,
                        p7_gate,
                    )
                )
    if len({row.norm for row in output}) != 16:
        raise AssertionError("the p7 Lambda22 norm alphabet has collisions")
    return tuple(output)


P7_LAMBDA22_NORM_CLASSES = _build_norm_classes()

COMPLETE_P7_LAMBDA22_NORM_PROOF = (
    "Put pi=sqrt(2), a=2*pi-1, epsilon=1+pi and delta=3+pi=epsilon*a. "
    "For a primitive r=(x,y), reflection integrality is q(r)|2I(r), "
    "where I(r)=(a*x,(y,Lambda22)). Away from pi, 3 and delta the spatial "
    "Gram matrix is unimodular, so a prime divisor of q would divide both "
    "x and y, contradicting primitivity. At delta the exponent is zero or "
    "one; when it is one, y lies in delta*Lambda22. At the inert prime 3, "
    "the Smith factor 6 permits exponent zero or one and the latter requires "
    "G*y divisible by 3. At pi, evenness and the Smith factors 2,6 give the "
    "complete exponents 0,2,3,4. For y=u+pi*v, exponent 3 requires G*u "
    "even and exponent 4 requires both G*u and G*v even. The local choices "
    "are independent. Since Z[sqrt(2)] has class number one and every "
    "totally-positive unit is a square, the sixteen products of the four "
    "dyadic representatives 1,2,4+2*pi,4 with optional factors 3 and delta "
    "are the complete canonical primitive mirror-norm alphabet. For a "
    "delta-class write y=delta*z and apply the pi/3 gates to z."
)


def p7_lambda22_factorization_audit() -> bool:
    return (
        P7_LAMBDA22_COEFFICIENT.sign(0) > 0
        and P7_LAMBDA22_COEFFICIENT.sign(1) < 0
        and P7_LAMBDA22_COEFFICIENT.norm() == -7
        and EPSILON * P7_LAMBDA22_COEFFICIENT == P7_GENERATOR
        and P7_GENERATOR.norm() == 7
        and P7_GENERATOR.sign(0) > 0
        and P7_GENERATOR.sign(1) > 0
    )


@dataclass(frozen=True, slots=True)
class P7Lambda22Control:
    norm_class: P7Lambda22NormClass
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_p7_lambda22_controls(
    height_bound: Silver = FIRST_PACKET_HEIGHT,
) -> tuple[P7Lambda22Control, ...]:
    """Enumerate every reflective control through an algebraic height."""

    if height_bound.sign(0) < 0:
        raise ValueError("height_bound must be nonnegative")
    coefficient = P7_LAMBDA22_COEFFICIENT
    levels: set[P7Lambda22Control] = set()
    for norm_class in P7_LAMBDA22_NORM_CLASSES:
        norm = norm_class.norm
        caps = (norm * height_bound / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if height > height_bound:
                continue
            levels.add(P7Lambda22Control(norm_class, height, control, target))
    return tuple(
        sorted(
            levels,
            key=lambda row: (row.height, row.norm_class.norm, row.control),
        )
    )


def reduced_spatial_target(level: P7Lambda22Control) -> Silver:
    """Return ``z^2`` after writing ``y=delta*z`` when required."""

    if level.norm_class.p7_gate:
        return level.target / (P7_GENERATOR * P7_GENERATOR)
    return level.target


def _gated_gram_triples_for(
    norm_class: P7Lambda22NormClass,
    target: Silver,
) -> tuple[tuple[int, int, int], ...]:
    if not target.is_integral or target.sign(0) < 0 or target.sign(1) < 0:
        return ()
    u_modulus, v_modulus = norm_class.coefficient_pairing_moduli
    cross_modulus = lcm(u_modulus, v_modulus)
    u_minimum = PAIRING_RADICAL_MINIMUM[u_modulus]
    v_minimum = PAIRING_RADICAL_MINIMUM[v_modulus]
    output = []
    for item in candidate_leech_gram_triples(target):
        if has_short_binary_combination(item):
            continue
        u_square, v_square, u_dot_v = (
            item.u_square,
            item.v_square,
            item.u_dot_v,
        )
        if (
            (u_square and u_square < u_minimum)
            or (v_square and v_square < v_minimum)
            or u_square % u_modulus
            or v_square % v_modulus
            or u_dot_v % cross_modulus
        ):
            continue
        output.append((u_square, v_square, u_dot_v))
    return tuple(output)


def ambient_gram_triples(level: P7Lambda22Control) -> tuple[tuple[int, int, int], ...]:
    target = reduced_spatial_target(level)
    if not target.is_integral or target.sign(0) < 0 or target.sign(1) < 0:
        return ()
    return tuple(
        (item.u_square, item.v_square, item.u_dot_v)
        for item in candidate_leech_gram_triples(target)
        if not has_short_binary_combination(item)
    )


def gated_gram_triples(level: P7Lambda22Control) -> tuple[tuple[int, int, int], ...]:
    return _gated_gram_triples_for(level.norm_class, reduced_spatial_target(level))


def height_zero_gram_triples(
    norm_class: P7Lambda22NormClass,
) -> tuple[tuple[int, int, int], ...]:
    target = norm_class.norm
    if norm_class.p7_gate:
        target /= P7_GENERATOR * P7_GENERATOR
    return _gated_gram_triples_for(norm_class, target)


def lambda22_pair_silver(
    left: Sequence[Silver | int], right: Sequence[Silver | int]
) -> Silver:
    if len(left) != 22 or len(right) != 22:
        raise ValueError("Lambda22 vectors have length 22")
    gram = lambda22_gram()
    return sum(
        (
            Silver(gram[i][j]) * Silver.coerce(left[i]) * Silver.coerce(right[j])
            for i in range(22)
            for j in range(22)
            if gram[i][j] and left[i] and right[j]
        ),
        start=SILVER_ZERO,
    )


def p7_lambda22_lorentz_pair(
    left: Sequence[Silver | int], right: Sequence[Silver | int]
) -> Silver:
    if len(left) != 23 or len(right) != 23:
        raise ValueError("the p7 Lambda22 host has rank 23")
    return (
        -P7_LAMBDA22_COEFFICIENT
        * Silver.coerce(left[0])
        * Silver.coerce(right[0])
        + lambda22_pair_silver(left[1:], right[1:])
    )


def first_packet_root(d: Sequence[int]) -> tuple[Silver, ...]:
    if len(d) != 22 or lambda22_pair_silver(d, d) != FIRST_PACKET_SHELL_NORM:
        raise ValueError("expected a Lambda22 minimal vector")
    return (
        FIRST_PACKET_CONTROL,
        *(EPSILON * coordinate for coordinate in d),
    )


def first_packet_reflection_integral(d: Sequence[int]) -> bool:
    root = first_packet_root(d)
    gram = lambda22_gram()
    return (
        (2 * P7_LAMBDA22_COEFFICIENT * root[0] / FIRST_PACKET_NORM).is_integral
        and all(
            (
                2
                * sum(
                    (
                        Silver(gram[i][j]) * root[j + 1]
                        for j in range(22)
                        if gram[i][j]
                    ),
                    start=SILVER_ZERO,
                )
                / FIRST_PACKET_NORM
            ).is_integral
            for i in range(22)
        )
    )


def first_packet_pairing_from_dot(dot: int) -> Silver:
    return (
        -P7_LAMBDA22_COEFFICIENT * FIRST_PACKET_CONTROL * FIRST_PACKET_CONTROL
        + EPSILON * EPSILON * dot
    )


def dual_direction_audit() -> bool:
    gram = lambda22_gram()
    image = tuple(
        sum(gram[i][j] * P7_DUAL_DIRECTION[j] for j in range(22))
        for i in range(22)
    )
    square = sum(
        P7_DUAL_DIRECTION[i] * image[i] for i in range(22)
    )
    return image == (0,) * 21 + (3,) and square == 24


def strict_hyperideal_witness() -> tuple[Silver, ...]:
    """Return ``(1,(4-sqrt(2))*k/8)`` for the pinned dual vector ``k``."""

    scale = (Silver(4) - SQRT2) / 8
    return (SILVER_ONE, *(scale * coordinate for coordinate in P7_DUAL_DIRECTION))


def strict_witness_norm() -> Silver:
    witness = strict_hyperideal_witness()
    return p7_lambda22_lorentz_pair(witness, witness)


def strict_witness_pairing_from_last_coordinate(last_coordinate: int) -> Silver:
    """Pair the witness with the minimal-shell wall indexed by ``d``.

    Since ``G*k=3*e_22``, one has ``(k,d)=3*d_22``.  The exact PARI shell
    census proves ``|d_22|<=2`` over the complete minimal shell.
    """

    ax = P7_LAMBDA22_COEFFICIENT * FIRST_PACKET_CONTROL
    return ax * (Silver(-1) + Silver(3 * last_coordinate) / 8)


def first_packet_control() -> P7Lambda22Control:
    matches = tuple(
        row
        for row in enumerate_p7_lambda22_controls()
        if row.height == FIRST_PACKET_HEIGHT
        and row.norm_class.norm == FIRST_PACKET_NORM
    )
    if len(matches) != 1:
        raise AssertionError("the first p7 Lambda22 packet control drifted")
    return matches[0]


__all__ = [
    "COMPLETE_P7_LAMBDA22_NORM_PROOF",
    "EPSILON",
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_SHELL_NORM",
    "FIRST_PACKET_TARGET",
    "FIRST_PACKET_WALL_COUNT",
    "GRAM_FILE_SHA256",
    "GRAM_STREAM_SHA256",
    "P7_DUAL_DIRECTION",
    "P7_GENERATOR",
    "P7_LAMBDA22_COEFFICIENT",
    "P7_LAMBDA22_NORM_CLASSES",
    "PAIRING_RADICAL_MINIMUM",
    "SOURCE_HTML_SHA256",
    "SOURCE_URL",
    "P7Lambda22Control",
    "P7Lambda22NormClass",
    "ambient_gram_triples",
    "dual_direction_audit",
    "enumerate_p7_lambda22_controls",
    "first_packet_control",
    "first_packet_pairing_from_dot",
    "first_packet_reflection_integral",
    "first_packet_root",
    "gated_gram_triples",
    "height_zero_gram_triples",
    "lambda22_basis_audit",
    "lambda22_gram_file_sha256",
    "lambda22_gram_stream_sha256",
    "lambda22_pair_silver",
    "lambda22_smith_invariants",
    "p7_lambda22_factorization_audit",
    "p7_lambda22_lorentz_pair",
    "reduced_spatial_target",
    "strict_hyperideal_witness",
    "strict_witness_norm",
    "strict_witness_pairing_from_last_coordinate",
]
