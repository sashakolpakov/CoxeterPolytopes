"""Exact first positive Vinberg packet for a silver ``Lambda_22`` host.

The host is ``<-a> + Lambda_22`` over ``Z[sqrt(2)]``, where
``a=sqrt(2)-1``.  The positive block is the exact integral Gram matrix in
the Nebe--Sloane catalogue.  The module certifies only the prefix through
height ``14+10*sqrt(2)``.  It does not claim Vinberg termination.

The nonunimodular point is essential.  The Smith form of ``Lambda_22`` is
``1^20,2,6``.  Consequently the dyadic norm valuation is not merely
squarefree: the complete unit-coefficient alphabet has valuations
``0,2,3,4`` at ``pi=sqrt(2)``.  The inert prime 3 may occur once.  The
corresponding spatial pairing gates are tested before a shell is accepted.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from typing import Sequence

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.cubic7_lambda22 import (
    GRAM_PATH,
    GRAM_FILE_SHA256,
    GRAM_STREAM_SHA256,
    PARI_QFMINIM_CERTIFICATES,
    SOURCE_HTML_SHA256,
    SOURCE_URL,
    lambda22_basis_audit,
    lambda22_gram_file_sha256,
    lambda22_pairing_radical_basis,
    lambda22_pairing_radical_gram,
    lambda22_smith_invariants,
    lambda22_gram,
    lambda22_gram_stream_sha256,
)
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import has_short_binary_combination
from coxeter_polytopes.search.silver_odd_leech import psd_extension_joint_maximum


LAMBDA22_COEFFICIENT = SQRT2 - SILVER_ONE
FIRST_PACKET_HEIGHT = Silver(14, 10)
FIRST_PACKET_NORM = Silver(2)
FIRST_PACKET_CONTROL = Silver(6, 4)
FIRST_PACKET_TARGET = Silver(30, 20)
FIRST_PACKET_SHELL_NORM = 10
FIRST_PACKET_DYADIC_SHELL_NORM = 20
LAMBDA22_NORM10_SHELL_COUNT = 467_596_800
S2_MINIMUM = 12
S3_MINIMUM = 24
S6_MINIMUM = 96
S2_NORM20_SHELL_COUNT = 1_368_576

# pi-adic exponents 0,2,3,4, then optional exponent one at the inert prime 3.
LAMBDA22_MIRROR_NORMS: tuple[Silver, ...] = (
    Silver(1),
    Silver(2),
    Silver(4, 2),
    Silver(4),
    Silver(3),
    Silver(6),
    Silver(12, 6),
    Silver(12),
)

LAMBDA22_NORM_GATE_PROOF = (
    "Put pi=sqrt(2), a=pi-1, and I(r)=(a*x,(y,Lambda22)) for a "
    "primitive vector r=(x,y). Reflection integrality is (q(r))|2I(r). "
    "Away from pi and 3 the spatial Gram is unimodular, hence no prime can "
    "divide q without contradicting primitivity. At the inert prime 3, the "
    "Smith factor 6 gives exponent at most one; if it occurs then "
    "G*y is divisible by 3. At pi, if x is a unit then q=-a*x^2 mod pi and "
    "the exponent is zero. Otherwise the even spatial form excludes "
    "exponent one. The Smith factors 2 and 6 give min(v_pi(I),2)<=2, so "
    "q|2I gives v_pi(q)<=4. Thus the exact dyadic exponents are among "
    "0,2,3,4. Their totally-positive representatives are 1,2,4+2*pi,4; "
    "multiplication by 3 gives the four 3-primary classes. At exponent 3, "
    "G*y divisible by pi means G*u is even for y=u+pi*v. At exponent 4, "
    "G*y divisible by pi^2 means both G*u and G*v are even. For the "
    "3-primary classes both G*u and G*v are divisible by 3. Empty classes remain part of "
    "the complete alphabet, so no realization assumption enters the proof."
)


@dataclass(frozen=True, slots=True)
class SilverLambda22Control:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    dyadic_valuation: int
    ternary_gate: bool


@dataclass(frozen=True, slots=True)
class SilverLambda22NecessaryComponent:
    level: SilverLambda22Control
    triples: tuple[tuple[int, int, int], ...]
    exclusion: str


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


def lambda22_lorentz_pair_silver(
    left: Sequence[Silver | int], right: Sequence[Silver | int]
) -> Silver:
    if len(left) != 23 or len(right) != 23:
        raise ValueError("the silver Lambda22 host has rank 23")
    return (
        -LAMBDA22_COEFFICIENT
        * Silver.coerce(left[0])
        * Silver.coerce(right[0])
        + lambda22_pair_silver(left[1:], right[1:])
    )


def _norm_gates(norm: Silver) -> tuple[int, bool]:
    if norm not in LAMBDA22_MIRROR_NORMS:
        raise ValueError("norm is outside the complete Lambda22 alphabet")
    dyadic = {
        Silver(1): 0,
        Silver(2): 2,
        Silver(4, 2): 3,
        Silver(4): 4,
        Silver(3): 0,
        Silver(6): 2,
        Silver(12, 6): 3,
        Silver(12): 4,
    }[norm]
    ternary = norm in (Silver(3), Silver(6), Silver(12, 6), Silver(12))
    return dyadic, ternary


@lru_cache(maxsize=None)
def enumerate_lambda22_controls(
    height_bound: Silver = FIRST_PACKET_HEIGHT,
) -> tuple[SilverLambda22Control, ...]:
    levels: set[SilverLambda22Control] = set()
    for norm in LAMBDA22_MIRROR_NORMS:
        dyadic, ternary = _norm_gates(norm)
        for control in _bounded_for_caps(
            (norm * height_bound / LAMBDA22_COEFFICIENT, norm / (-LAMBDA22_COEFFICIENT))
        ):
            if control.sign(0) <= 0:
                continue
            if not (2 * LAMBDA22_COEFFICIENT * control / norm).is_integral:
                continue
            target = norm + LAMBDA22_COEFFICIENT * control * control
            if target.sign(1) < 0:
                continue
            height = LAMBDA22_COEFFICIENT * control * control / norm
            if height > height_bound:
                continue
            levels.add(
                SilverLambda22Control(
                    norm, height, control, target, dyadic, ternary
                )
            )
    return tuple(sorted(levels, key=lambda row: (row.height, row.norm, row.control)))


def necessary_components_through_first_packet() -> tuple[SilverLambda22NecessaryComponent, ...]:
    output = []
    for level in enumerate_lambda22_controls():
        triples = tuple(
            (item.u_square, item.v_square, item.u_dot_v)
            for item in candidate_leech_gram_triples(level.target)
            if not has_short_binary_combination(item)
        )
        exclusion = "ordinary minimum-four/Gram exclusion"
        if level.height == Silver(2, Fraction(3, 2)):
            exclusion = "u=v has norm six but Gu even; the dyadic-gate sublattice has minimum twelve"
        elif level.height == Silver(2, 2):
            exclusion = "dyadic gate G*y=0 mod 2 has exact sublattice minimum twelve"
        elif level.height == Silver(6, Fraction(9, 2)):
            exclusion = "u in S6 and v in S3 imply u-v in S3 of norm at most twelve, below min(S3)=24"
        elif level.height == FIRST_PACKET_HEIGHT and level.norm == FIRST_PACKET_NORM:
            exclusion = "realized: u=v=w and w runs over the complete norm-ten shell"
        elif level.height == FIRST_PACKET_HEIGHT and level.norm == Silver(4):
            exclusion = "realized: u=v runs over the complete norm-twenty shell of S2={u:G*u even}"
        elif triples:
            exclusion = "excluded by exact gate congruences/minimum; see frozen audit"
        output.append(SilverLambda22NecessaryComponent(level, triples, exclusion))
    return tuple(output)


def first_packet_root(w: Sequence[int]) -> tuple[Silver, ...]:
    if len(w) != 22 or lambda22_pair_silver(w, w) != FIRST_PACKET_SHELL_NORM:
        raise ValueError("expected a Lambda22 norm-ten shell vector")
    return (FIRST_PACKET_CONTROL, *((SILVER_ONE + SQRT2) * value for value in w))


def first_packet_dyadic_root(u: Sequence[int]) -> tuple[Silver, ...]:
    if len(u) != 22 or lambda22_pair_silver(u, u) != FIRST_PACKET_DYADIC_SHELL_NORM:
        raise ValueError("expected a norm-twenty vector in the dyadic gate")
    gram = lambda22_gram()
    if any(sum(gram[i][j] * u[j] for j in range(22)) % 2 for i in range(22)):
        raise ValueError("the norm-twenty vector is outside S2")
    return (Silver(8, 6), *((SILVER_ONE + SQRT2) * value for value in u))


def first_packet_reflection_integral(w: Sequence[int]) -> bool:
    root = first_packet_root(w)
    gram = lambda22_gram()
    spatial = root[1:]
    return (
        (2 * LAMBDA22_COEFFICIENT * root[0] / FIRST_PACKET_NORM).is_integral
        and all(
            (
                2
                * sum(
                    (Silver(gram[i][j]) * spatial[j] for j in range(22) if gram[i][j]),
                    start=SILVER_ZERO,
                )
                / FIRST_PACKET_NORM
            ).is_integral
            for i in range(22)
        )
    )


def first_packet_dyadic_reflection_integral(u: Sequence[int]) -> bool:
    root = first_packet_dyadic_root(u)
    gram = lambda22_gram()
    return (
        (2 * LAMBDA22_COEFFICIENT * root[0] / Silver(4)).is_integral
        and all(
            (
                sum(
                    (Silver(gram[i][j]) * root[j + 1] for j in range(22) if gram[i][j]),
                    start=SILVER_ZERO,
                )
                / Silver(2)
            ).is_integral
            for i in range(22)
        )
    )


def hyperideal_witness(d: Sequence[int]) -> tuple[Silver, ...]:
    if len(d) != 22 or lambda22_pair_silver(d, d) != 4:
        raise ValueError("expected a Lambda22 minimal vector")
    return (SILVER_ONE, *(Silver(value, 0) / 3 for value in d))


def witness_norm(d: Sequence[int]) -> Silver:
    return lambda22_lorentz_pair_silver(hyperideal_witness(d), hyperideal_witness(d))


def witness_pairing_from_dot(dot: int) -> Silver:
    """Pair the witness with a first-packet root having ``d.w=dot``."""

    return -(SQRT2 - SILVER_ONE) * FIRST_PACKET_CONTROL + (SILVER_ONE + SQRT2) * Silver(dot) / 3


def first_packet_pairing_from_dot(dot: int) -> Silver:
    """Pair two first-packet roots whose underlying norm-ten vectors have dot."""

    return -LAMBDA22_COEFFICIENT * FIRST_PACKET_CONTROL**2 + (Silver(3, 2) * dot)


def dyadic_packet_pairing_from_dot(dot: int) -> Silver:
    return -LAMBDA22_COEFFICIENT * Silver(8, 6) ** 2 + Silver(3, 2) * dot


def cross_packet_pairing_from_dot(dot: int) -> Silver:
    return (
        -LAMBDA22_COEFFICIENT * FIRST_PACKET_CONTROL * Silver(8, 6)
        + Silver(3, 2) * dot
    )


def dyadic_witness_pairing_from_dot(dot: int) -> Silver:
    return -LAMBDA22_COEFFICIENT * Silver(8, 6) + (SILVER_ONE + SQRT2) * Silver(dot) / 3


def minimal_direction_dot_bound() -> int:
    """Return the exact bound ``|d.w|<=5`` for ``d^2=4,w^2=10``.

    Integrality gives an integer dot product.  If ``|d.w|>=6``, one of
    ``w+d`` and ``w-d`` has square at most two.  It cannot be zero because
    the norms differ, contradicting the Lambda22 minimum four.
    """

    return 5


def relaxed_witness_bound(
    u_square: int, v_square: int, u_dot_v: int
) -> tuple[int, int, Silver]:
    """Exact PSD upper bound for ``d.(u+sqrt(2)v)``, with ``d^2=4``."""

    return psd_extension_joint_maximum(
        u_square, v_square, u_dot_v, 4, SQRT2
    )


def gated_relaxed_witness_bound(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    p_modulus: int = 1,
    q_modulus: int = 1,
) -> tuple[int, int, Silver]:
    """PSD/minimum-four bound with exact pairing-gate divisibilities.

    The scan is a necessary superset of realizable triples.  It checks every
    combination ``a*u+b*v+c*d`` with coefficients in ``[-3,3]`` and rejects
    nonzero norms below the Lambda22 minimum four.
    """

    if p_modulus <= 0 or q_modulus <= 0:
        raise ValueError("pairing moduli must be positive")
    best: tuple[int, int, Silver] | None = None
    for p in range(-isqrt(4 * u_square), isqrt(4 * u_square) + 1):
        if p % p_modulus:
            continue
        for q in range(-isqrt(4 * v_square), isqrt(4 * v_square) + 1):
            if q % q_modulus:
                continue
            determinant = (
                4 * (u_square * v_square - u_dot_v * u_dot_v)
                - v_square * p * p
                - u_square * q * q
                + 2 * u_dot_v * p * q
            )
            if determinant < 0:
                continue
            valid = True
            for a in range(-3, 4):
                for b in range(-3, 4):
                    for c in range(-3, 4):
                        if not (a or b or c):
                            continue
                        square = (
                            a * a * u_square
                            + b * b * v_square
                            + 4 * c * c
                            + 2 * a * b * u_dot_v
                            + 2 * a * c * p
                            + 2 * b * c * q
                        )
                        if 0 < square < 4:
                            valid = False
                            break
                    if not valid:
                        break
                if not valid:
                    break
            if not valid:
                continue
            value = Silver(p, q)
            if best is None or value > best[2]:
                best = (p, q, value)
    if best is None:
        raise AssertionError("no gated PSD extension survived")
    return best


__all__ = [
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_SHELL_NORM",
    "FIRST_PACKET_DYADIC_SHELL_NORM",
    "LAMBDA22_NORM10_SHELL_COUNT",
    "S2_MINIMUM",
    "S3_MINIMUM",
    "S6_MINIMUM",
    "S2_NORM20_SHELL_COUNT",
    "FIRST_PACKET_TARGET",
    "GRAM_PATH",
    "GRAM_FILE_SHA256",
    "GRAM_STREAM_SHA256",
    "LAMBDA22_COEFFICIENT",
    "LAMBDA22_MIRROR_NORMS",
    "LAMBDA22_NORM_GATE_PROOF",
    "SOURCE_HTML_SHA256",
    "SOURCE_URL",
    "SilverLambda22Control",
    "enumerate_lambda22_controls",
    "first_packet_pairing_from_dot",
    "first_packet_dyadic_reflection_integral",
    "first_packet_dyadic_root",
    "first_packet_reflection_integral",
    "first_packet_root",
    "hyperideal_witness",
    "lambda22_gram_stream_sha256",
    "lambda22_gram_file_sha256",
    "lambda22_basis_audit",
    "lambda22_smith_invariants",
    "lambda22_pairing_radical_basis",
    "lambda22_pairing_radical_gram",
    "lambda22_lorentz_pair_silver",
    "lambda22_pair_silver",
    "minimal_direction_dot_bound",
    "necessary_components_through_first_packet",
    "relaxed_witness_bound",
    "gated_relaxed_witness_bound",
    "witness_norm",
    "witness_pairing_from_dot",
    "dyadic_packet_pairing_from_dot",
    "cross_packet_pairing_from_dot",
    "dyadic_witness_pairing_from_dot",
]
