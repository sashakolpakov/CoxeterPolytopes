"""Exact first Vinberg packet for a cubic-seven odd-Leech section in H23.

Let ``O24`` be the odd Leech lattice in the catalogue model used by
``silver_odd_leech`` and let ``d=(1,...,1)/sqrt(8)``.  The positive lattice

``S23 = d^perp cap O24``

has rank 23, determinant 3, minimum 3, and a complete shell of 2576 minimal
vectors.  This module studies

``<-theta> + S23`` over ``Z[theta]``, where
``theta=2*cos(2*pi/7)`` is positive at the distinguished real place.

The conclusion is deliberately packet-local.  The first positive packet is
complete and has a strict hyperideal point, so that prefix is not a compact
or finite-volume chamber.  Later mirrors may cut the hyperideal patch.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from typing import Iterator, Sequence, TypeAlias

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp, smith_normal_form

from coxeter_polytopes.exact.cubic7 import CUBIC7_ZERO, THETA7, Cubic7
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps
from coxeter_polytopes.search.silver_odd_leech import (
    ODD_LEECH_BASIS_SCALED,
    bareiss_determinant,
    catalog_golay_codewords,
    odd_leech_pair_scaled,
    odd_leech_model_audit,
    odd_leech_raw_member,
)


RawVector: TypeAlias = tuple[int, ...]

SECTION_NORMAL_SCALED: RawVector = (1,) * 24
SECTION_WITNESS_SCALED: RawVector = (4, -4, *(0 for _ in range(22)))


def _section_pairing_row() -> tuple[int, ...]:
    return tuple(
        odd_leech_pair_scaled(row, SECTION_NORMAL_SCALED)
        for row in ODD_LEECH_BASIS_SCALED
    )


SECTION_PAIRING_ROW = _section_pairing_row()

# The first pairing is one.  Therefore e_i-p_i*e_0, i=1,...,23, is an
# integral basis of the full kernel of pairing with d, not merely a
# finite-index sublattice of it.
SECTION_BASIS_SCALED: tuple[RawVector, ...] = tuple(
    tuple(
        ODD_LEECH_BASIS_SCALED[index][coordinate]
        - SECTION_PAIRING_ROW[index] * ODD_LEECH_BASIS_SCALED[0][coordinate]
        for coordinate in range(24)
    )
    for index in range(1, 24)
)

SECTION_GRAM: tuple[tuple[int, ...], ...] = tuple(
    tuple(odd_leech_pair_scaled(left, right) for right in SECTION_BASIS_SCALED)
    for left in SECTION_BASIS_SCALED
)


def section_smith_invariants() -> tuple[int, ...]:
    smith = smith_normal_form(Matrix(SECTION_GRAM), domain=ZZ)
    return tuple(abs(int(smith[index, index])) for index in range(23))


def iter_minimal_section_vectors_scaled() -> Iterator[RawVector]:
    """Yield the complete norm-three shell of the rank-23 section.

    Odd-Leech norm-three vectors are the Golay sign vectors.  Orthogonality
    to the all-ones normal says exactly that the codeword has weight 12.
    """

    for codeword in catalog_golay_codewords():
        if codeword.bit_count() != 12:
            continue
        vector = tuple(
            -1 if codeword & (1 << coordinate) else 1
            for coordinate in range(24)
        )
        if sum(vector) or not odd_leech_raw_member(vector):
            raise AssertionError("a weight-12 Golay sign vector left the section")
        if odd_leech_pair_scaled(vector, vector) != 3:
            raise AssertionError("a section sign vector lost norm three")
        yield vector


def _rank_mod_prime(rows: Sequence[Sequence[int]], prime: int = 1_000_003) -> int:
    basis: list[list[int]] = []
    pivots: list[int] = []
    for source in rows:
        row = [int(value) % prime for value in source]
        for pivot, old in zip(pivots, basis):
            if row[pivot]:
                factor = row[pivot]
                row = [
                    (value - factor * old[index]) % prime
                    for index, value in enumerate(row)
                ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], prime - 2, prime)
        row = [(value * inverse) % prime for value in row]
        pivots.append(pivot)
        basis.append(row)
    return len(basis)


def _rank_mod_three(matrix: Sequence[Sequence[int]]) -> int:
    work = [[int(value) % 3 for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, 3)
        work[rank] = [(inverse * value) % 3 for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (value - factor * other) % 3
                for value, other in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


@dataclass(frozen=True, slots=True)
class OddLeechSectionModelAudit:
    determinant: int
    smith_invariants: tuple[int, ...]
    leading_principal_minors: tuple[int, ...]
    minimum: int
    minimal_shell_count: int
    minimal_shell_digest: str
    minimal_shell_span_rank: int
    radical_mod_three_rank: int
    radical_generator_coordinates: tuple[int, ...]
    radical_generator_raw_residue: tuple[int, ...]
    radical_minimum_lower_bound: int
    witness_square: int
    witness_shell_pairing_histogram: tuple[tuple[int, int], ...]


RADICAL_NORM_EXCLUSION = (
    "In the displayed section basis, ker(G mod 3) is generated by the last "
    "coordinate, whose raw vector is (-23,1^23), congruent to (1^24) mod 3. "
    "Thus every vector in the pairing radical has all raw coordinates "
    "congruent mod 3. Raw norm 3 means square sum 24: residue 0 is impossible, "
    "while residue +/-1 forces all coordinates +/-1 with one common sign, "
    "contradicting coordinate sum zero. Raw norms 4 and 5 have square sums "
    "32 and 40. In a nonzero residue class, common odd parity makes every "
    "coordinate equal to the least representative +/-1 unless the square "
    "sum is already at least 48, contradicting sum zero; common even parity "
    "has square sum at least 96. Residue zero is impossible modulo 9. "
    "For raw norm 6, square sum 48. "
    "Residue 0 is again impossible. Writing a residue-1 coordinate as 1+3k, "
    "the zero-sum and square-sum equations give sum(k)=-8 and sum(k^2)=8, "
    "hence exactly eight coordinates are -2 and sixteen are 1; this violates "
    "the common coordinate parity of every odd-Leech vector. Residue -1 is "
    "the negative case. Therefore the radical has no vector of norm at most "
    "6 (in particular no norm 3 or 6 vector)."
)


@lru_cache(maxsize=1)
def odd_leech_section_model_audit() -> OddLeechSectionModelAudit:
    parent = odd_leech_model_audit()
    if parent.determinant != 1 or parent.minimal_shell_count != 4096:
        raise AssertionError("the exact odd-Leech parent audit changed")
    if SECTION_PAIRING_ROW[0] != 1:
        raise AssertionError("the section kernel basis lost its primitive pivot")
    if not all(
        sum(row) == 0 and odd_leech_raw_member(row)
        for row in SECTION_BASIS_SCALED
    ):
        raise AssertionError("a section basis row left the exact odd-Leech kernel")
    determinant = bareiss_determinant(SECTION_GRAM)
    minors = tuple(
        bareiss_determinant(tuple(row[:size] for row in SECTION_GRAM[:size]))
        for size in range(1, 24)
    )
    smith = section_smith_invariants()
    if determinant != 3 or smith != (1,) * 22 + (3,) or any(value <= 0 for value in minors):
        raise AssertionError("the section is not the expected positive determinant-three lattice")

    shell = tuple(iter_minimal_section_vectors_scaled())
    digest = sha256()
    histogram: Counter[int] = Counter()
    for vector in shell:
        digest.update(bytes(value + 1 for value in vector))
        histogram[odd_leech_pair_scaled(SECTION_WITNESS_SCALED, vector)] += 1
    if len(shell) != 2576 or _rank_mod_prime(shell) != 23:
        raise AssertionError("the complete weight-12 shell census changed")

    radical_generator = (0,) * 22 + (1,)
    radical_raw = SECTION_BASIS_SCALED[-1]
    if (
        _rank_mod_three(SECTION_GRAM) != 22
        or any(value % 3 for value in Matrix(SECTION_GRAM) * Matrix(radical_generator))
        or tuple(value % 3 for value in radical_raw) != (1,) * 24
    ):
        raise AssertionError("the mod-three pairing radical changed")
    if not odd_leech_raw_member(SECTION_WITNESS_SCALED) or sum(SECTION_WITNESS_SCALED):
        raise AssertionError("the norm-four witness left the section")
    witness_square = odd_leech_pair_scaled(
        SECTION_WITNESS_SCALED, SECTION_WITNESS_SCALED
    )
    if witness_square != 4 or histogram != Counter({0: 1232, -1: 672, 1: 672}):
        raise AssertionError("the section witness audit changed")
    return OddLeechSectionModelAudit(
        determinant=determinant,
        smith_invariants=smith,
        leading_principal_minors=minors,
        minimum=3,
        minimal_shell_count=len(shell),
        minimal_shell_digest=digest.hexdigest(),
        minimal_shell_span_rank=23,
        radical_mod_three_rank=22,
        radical_generator_coordinates=radical_generator,
        radical_generator_raw_residue=tuple(value % 3 for value in radical_raw),
        radical_minimum_lower_bound=7,
        witness_square=witness_square,
        witness_shell_pairing_histogram=tuple(sorted(histogram.items())),
    )


COMPLETE_MIRROR_NORM_PROOF = (
    "For primitive r=(x,y), reflection integrality gives (q(r))|2I(r), "
    "I(r)=(theta*x,(y,S23)). The coefficient theta is a unit and S23 has "
    "Smith invariants 1^22,3. Away from the inert primes (2),(3), local "
    "unimodularity and primitivity force v(q)=0. At (2), I(r) is a unit and "
    "v_2(q)<=v_2(2)=1. At (3), q-divisibility forces x divisible by 3 and y "
    "into the one-dimensional pairing radical; its primitive pairing ideal "
    "has valuation exactly one, so v_3(q)<=1. The Minkowski bound 14/9 gives "
    "class number one, and -1, theta, theta^2-2 have full unit-signature rank, "
    "so totally positive units are squares. The complete canonical norms are "
    "therefore 1,2,3,6. For norms 3 and 6 every coefficient of y lies in the "
    "mod-three pairing radical."
)

CANONICAL_MIRROR_NORMS: tuple[Cubic7, ...] = tuple(
    Cubic7(value) for value in (1, 2, 3, 6)
)


@dataclass(frozen=True, slots=True, order=True)
class OddLeechSectionGramType:
    u_square: int
    v_square: int
    w_square: int
    u_dot_v: int
    u_dot_w: int
    v_dot_w: int
    determinant: int

    @property
    def matrix(self) -> tuple[tuple[int, int, int], ...]:
        return (
            (self.u_square, self.u_dot_v, self.u_dot_w),
            (self.u_dot_v, self.v_square, self.v_dot_w),
            (self.u_dot_w, self.v_dot_w, self.w_square),
        )

    @property
    def rank(self) -> int:
        return int(Matrix(self.matrix).rank())


def _rank_three_determinant(values: tuple[int, int, int, int, int, int]) -> int:
    u, v, w, uv, uw, vw = values
    return (
        u * v * w
        + 2 * uv * uw * vw
        - u * vw * vw
        - v * uw * uw
        - w * uv * uv
    )


def _positive_quotient(
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    raw = Matrix(matrix)
    if raw.det():
        return matrix
    diagonal, left, right = smith_normal_decomp(raw, domain=ZZ)
    if left * raw * right != diagonal:
        raise AssertionError("the Smith transport convention changed")
    nonzero = tuple(
        index for index in range(raw.rows) if diagonal[index, index] != 0
    )
    zero = tuple(index for index in range(raw.rows) if index not in nonzero)
    if zero and raw * right[:, zero] != Matrix.zeros(raw.rows, len(zero)):
        raise AssertionError("the Smith zero columns do not span the integral kernel")
    complement = right[:, nonzero]
    quotient = complement.T * raw * complement
    if quotient.det() <= 0:
        raise AssertionError("the coefficient Gram quotient is not positive definite")
    return tuple(
        tuple(int(quotient[i, j]) for j in range(len(nonzero)))
        for i in range(len(nonzero))
    )


def _integer_gram_minimum_at_least(
    matrix: tuple[tuple[int, ...], ...], minimum: int
) -> bool:
    quotient = _positive_quotient(matrix)
    if not quotient:
        return False
    inverse = Matrix(quotient).inv()
    bounds = []
    for index in range(len(quotient)):
        diagonal = inverse[index, index]
        cap = Fraction(int(diagonal.p), int(diagonal.q)) * (minimum - 1)
        bounds.append(isqrt(cap.numerator // cap.denominator))

    def visit(index: int, vector: tuple[int, ...]) -> bool:
        if index == len(bounds):
            if not any(vector):
                return True
            square = sum(
                quotient[i][j] * vector[i] * vector[j]
                for i in range(len(vector))
                for j in range(len(vector))
            )
            return not (0 < square < minimum)
        return all(
            visit(index + 1, (*vector, value))
            for value in range(-bounds[index], bounds[index] + 1)
        )

    return visit(0, ())


@lru_cache(maxsize=None)
def section_rank_three_gram_types(
    target: Cubic7,
    coefficient_minimum: int = 3,
    pairing_modulus: int = 1,
) -> tuple[OddLeechSectionGramType, ...]:
    """Enumerate every necessary Gram type for ``y=u+theta*v+theta^2*w``."""

    target = Cubic7.coerce(target)
    if not target.is_integral or coefficient_minimum <= 0 or pairing_modulus <= 0:
        return ()
    a, b, c = (int(target.a), int(target.b), int(target.c))
    trace = int(target.trace())
    if trace < 0:
        return ()
    output = []
    for v_square in range(2 * trace // 7 + 1):
        if 0 < v_square < coefficient_minimum:
            continue
        for w_square in range(2 * trace // 7 + 1):
            if 0 < w_square < coefficient_minimum:
                continue
            cauchy = isqrt(v_square * w_square)
            for v_dot_w in range(-cauchy, cauchy + 1):
                u_square = a - 2 * v_dot_w + w_square
                if (
                    u_square < 0
                    or u_square > trace
                    or 0 < u_square < coefficient_minimum
                ):
                    continue
                twice_u_dot_v = b + w_square - 4 * v_dot_w
                twice_u_dot_w = c - v_square - 3 * w_square + 2 * v_dot_w
                if twice_u_dot_v % 2 or twice_u_dot_w % 2:
                    continue
                u_dot_v = twice_u_dot_v // 2
                u_dot_w = twice_u_dot_w // 2
                if (
                    u_dot_v * u_dot_v > u_square * v_square
                    or u_dot_w * u_dot_w > u_square * w_square
                ):
                    continue
                values = (
                    u_square,
                    v_square,
                    w_square,
                    u_dot_v,
                    u_dot_w,
                    v_dot_w,
                )
                determinant = _rank_three_determinant(values)
                if determinant < 0 or any(value % pairing_modulus for value in values):
                    continue
                gram = (
                    (u_square, u_dot_v, u_dot_w),
                    (u_dot_v, v_square, v_dot_w),
                    (u_dot_w, v_dot_w, w_square),
                )
                if not _integer_gram_minimum_at_least(gram, coefficient_minimum):
                    continue
                output.append(OddLeechSectionGramType(*values, determinant))
    return tuple(sorted(output))


FIRST_PACKET_NORM = Cubic7(2)
FIRST_PACKET_CONTROL = Cubic7(2) + THETA7
FIRST_PACKET_SCALE = Cubic7(1) + THETA7
FIRST_PACKET_TARGET = Cubic7(3, 6, 3)
FIRST_PACKET_HEIGHT = Cubic7(Fraction(1, 2), 3, Fraction(3, 2))


@dataclass(frozen=True, slots=True, order=True)
class OddLeechSectionControlLevel:
    height: Cubic7
    norm: Cubic7
    control: Cubic7
    target: Cubic7
    pairing_modulus: int
    coefficient_minimum: int
    necessary_gram_types: tuple[OddLeechSectionGramType, ...]


@lru_cache(maxsize=1)
def section_controls_through_first_packet() -> tuple[OddLeechSectionControlLevel, ...]:
    levels = set()
    # The endpoint is below seven at the distinguished embedding.  The
    # integer cap seven is a complete Minkowski superset; exact comparison
    # removes its overrun.
    for norm in CANONICAL_MIRROR_NORMS:
        caps = (norm * 7 / THETA7, norm / (-THETA7), norm / (-THETA7))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * THETA7 * control / norm).is_integral:
                continue
            target = norm + THETA7 * control * control
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = THETA7 * control * control / norm
            if (FIRST_PACKET_HEIGHT - height).sign(0) < 0:
                continue
            modulus = 3 if norm in (Cubic7(3), Cubic7(6)) else 1
            minimum = 7 if modulus == 3 else 3
            levels.add(
                OddLeechSectionControlLevel(
                    height,
                    norm,
                    control,
                    target,
                    modulus,
                    minimum,
                    section_rank_three_gram_types(target, minimum, modulus),
                )
            )
    return tuple(sorted(levels))


@dataclass(frozen=True, slots=True)
class OddLeechSectionFirstPacketAudit:
    model: OddLeechSectionModelAudit
    canonical_mirror_norms: tuple[Cubic7, ...]
    height_zero_ambient_types: tuple[tuple[Cubic7, tuple[OddLeechSectionGramType, ...]], ...]
    height_zero_gated_types: tuple[tuple[Cubic7, tuple[OddLeechSectionGramType, ...]], ...]
    height_zero_wall_count: int
    controls: tuple[OddLeechSectionControlLevel, ...]
    first_height: Cubic7
    first_norm: Cubic7
    first_control: Cubic7
    first_scale: Cubic7
    first_target: Cubic7
    first_gram_type: OddLeechSectionGramType
    first_packet_root_count: int
    distinct_spatial_pairings: tuple[int, ...]
    distinct_root_pairings: tuple[Cubic7, ...]
    witness_square: int
    witness_maximum_pairing: int
    witness_minimum_wall_margin: Cubic7
    witness_hyperideal_excess: Cubic7
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=1)
def odd_leech_section_first_packet_audit() -> OddLeechSectionFirstPacketAudit:
    model = odd_leech_section_model_audit()
    ambient = tuple(
        (norm, section_rank_three_gram_types(norm, 3, 1))
        for norm in CANONICAL_MIRROR_NORMS
    )
    gated = tuple(
        (
            norm,
            section_rank_three_gram_types(
                norm,
                7 if norm in (Cubic7(3), Cubic7(6)) else 3,
                3 if norm in (Cubic7(3), Cubic7(6)) else 1,
            ),
        )
        for norm in CANONICAL_MIRROR_NORMS
    )
    if tuple(len(types) for _norm, types in ambient) != (0, 0, 1, 1):
        raise AssertionError("the height-zero ambient Gram census changed")
    if any(types for _norm, types in gated):
        raise AssertionError("a height-zero reflecting Gram type survived")

    controls = section_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_types)
    if len(controls) != 5 or len(nonempty) != 1:
        raise AssertionError("the complete first-prefix control census changed")
    first = nonempty[0]
    expected_gram = OddLeechSectionGramType(3, 3, 0, 3, 0, 0, 0)
    if (
        first.height != FIRST_PACKET_HEIGHT
        or first.norm != FIRST_PACKET_NORM
        or first.control != FIRST_PACKET_CONTROL
        or first.target != FIRST_PACKET_TARGET
        or first.necessary_gram_types != (expected_gram,)
    ):
        raise AssertionError("the first section packet changed")
    if (
        -THETA7 * FIRST_PACKET_CONTROL * FIRST_PACKET_CONTROL
        + 3 * FIRST_PACKET_SCALE * FIRST_PACKET_SCALE
        != FIRST_PACKET_NORM
        or abs(FIRST_PACKET_CONTROL.norm()) != 1
        or abs(FIRST_PACKET_SCALE.norm()) != 1
    ):
        raise AssertionError("the displayed packet roots lost norm or primitivity")

    # For distinct norm-three vectors b,c, minimum three applied to b-c and
    # b+c gives b.c in {-1,0,1}, apart from the opposite case b.c=-3.
    spatial_pairings = (-3, -1, 0, 1)
    root_pairings = tuple(
        -THETA7 * FIRST_PACKET_CONTROL * FIRST_PACKET_CONTROL
        + FIRST_PACKET_SCALE * FIRST_PACKET_SCALE * Cubic7(pairing)
        for pairing in spatial_pairings
    )
    if any(pairing.sign(0) >= 0 for pairing in root_pairings):
        raise AssertionError("the simultaneous first packet ceased to be acute")

    maximum = max(pairing for pairing, _count in model.witness_shell_pairing_histogram)
    margin = THETA7 * FIRST_PACKET_CONTROL - FIRST_PACKET_SCALE * maximum
    excess = Cubic7(model.witness_square) - THETA7
    if margin.sign(0) <= 0 or excess.sign(0) <= 0:
        raise AssertionError("the strict hyperideal witness ceased to work")
    return OddLeechSectionFirstPacketAudit(
        model=model,
        canonical_mirror_norms=CANONICAL_MIRROR_NORMS,
        height_zero_ambient_types=ambient,
        height_zero_gated_types=gated,
        height_zero_wall_count=0,
        controls=controls,
        first_height=first.height,
        first_norm=first.norm,
        first_control=first.control,
        first_scale=FIRST_PACKET_SCALE,
        first_target=first.target,
        first_gram_type=expected_gram,
        first_packet_root_count=model.minimal_shell_count,
        distinct_spatial_pairings=spatial_pairings,
        distinct_root_pairings=root_pairings,
        witness_square=model.witness_square,
        witness_maximum_pairing=maximum,
        witness_minimum_wall_margin=margin,
        witness_hyperideal_excess=excess,
        finite_volume=False,
        compact=False,
    )


__all__ = [
    "CANONICAL_MIRROR_NORMS",
    "COMPLETE_MIRROR_NORM_PROOF",
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_SCALE",
    "FIRST_PACKET_TARGET",
    "RADICAL_NORM_EXCLUSION",
    "SECTION_BASIS_SCALED",
    "SECTION_GRAM",
    "SECTION_NORMAL_SCALED",
    "SECTION_PAIRING_ROW",
    "SECTION_WITNESS_SCALED",
    "OddLeechSectionControlLevel",
    "OddLeechSectionFirstPacketAudit",
    "OddLeechSectionGramType",
    "OddLeechSectionModelAudit",
    "iter_minimal_section_vectors_scaled",
    "odd_leech_section_first_packet_audit",
    "odd_leech_section_model_audit",
    "section_controls_through_first_packet",
    "section_rank_three_gram_types",
    "section_smith_invariants",
]
