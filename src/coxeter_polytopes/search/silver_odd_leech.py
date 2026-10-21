"""Exact first Vinberg packet for a silver odd-Leech host in H24.

The host is

``<-a> + O_24`` over ``Z[sqrt(2)]``, with ``a=-3+4*sqrt(2)`` and
``O_24`` the odd Leech lattice.  This module deliberately certifies only the
complete prefix through the first nonempty positive-height packet.  The
prefix has a strict hyperideal patch, so it is not a finite-volume chamber;
later mirrors can still cut that patch.

The odd Leech lattice is represented in coordinates scaled by ``sqrt(8)``.
For the catalogue Golay code ``C`` and
``Z_+^24={z in Z^24 : sum(z) is even}``, its two raw-coordinate cosets are

``2*C + 4*Z_+^24`` and ``1 + 2*C + 4*Z_+^24``.

Thus the actual Euclidean pairing is the raw dot product divided by eight.
The shifted coset has norm at least three, with equality exactly at the 4096
Golay sign vectors.  The unshifted coset has norm at least four.  This proves
the full minimal-shell census; it does not assume that an arbitrary set of
4096 sign patterns is the odd Leech shell.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from typing import Iterator, Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


RawVector: TypeAlias = tuple[int, ...]

EPSILON = SILVER_ONE + SQRT2
ODD_LEECH_COEFFICIENT = Silver(-3, 4)
P23_TOTALLY_POSITIVE_GENERATOR = Silver(5, 1)

# The dyadic prime is (sqrt(2)); canonical totally-positive generators for
# its exponents 0,1,2 are 1, 2+sqrt(2), 2.  The independent p23 exponent is
# 0 or 1.  Multiplying gives the six normalized mirror-norm candidates.
ODD_LEECH_MIRROR_NORMS: tuple[Silver, ...] = (
    Silver(1),
    Silver(2, 1),
    Silver(2),
    Silver(5, 1),
    Silver(12, 7),
    Silver(10, 2),
)
P23_MIRROR_NORMS = frozenset(ODD_LEECH_MIRROR_NORMS[3:])

ODD_LEECH_NORM_GATE_PROOF = (
    "Put pi=sqrt(2), a=-3+4*pi, and p=(a), of ideal norm 23. For a "
    "primitive r=(x,y), reflection integrality is (q(r))|2I(r), with "
    "I(r)=(a*x,(y,O24)). Away from (pi) and p, any prime dividing q "
    "would divide x and, by unimodularity, y, contradicting primitivity. "
    "At (pi), a is a unit and primitivity makes I(r) a unit, so "
    "v_pi(q)<=v_pi(2)=2; oddness supplies no even-lattice lower bound. "
    "Thus the dyadic exponents are among 0,1,2. At p, a positive "
    "q-valuation and q|2(y,O24) force y in p*O24. Then x is a p-unit, "
    "I(r) has valuation one, and q=-a*x^2 mod p^2 has valuation exactly "
    "one. Hence the p-exponent is 0 or 1, and the exponent-one classes "
    "have the spatial gate y in a*O24. Modulo totally-positive unit "
    "squares the six generators are 1,2+pi,2,5+pi,12+7*pi,10+2*pi."
)


# Integral basis from the Nebe--Sloane Catalogue of Lattices, on the raw
# sqrt(8) scale.  Its raw Gram divided by eight is integral and unimodular.
ODD_LEECH_BASIS_SCALED: tuple[RawVector, ...] = (
    (8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0),
    (4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0),
    (2, 0, 2, 0, 2, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0),
    (2, 0, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0),
    (2, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0),
    (0, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0),
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
)


def _raw_word(row: Sequence[int]) -> int:
    """Return the mod-four codeword encoded by an even raw row."""

    if len(row) != 24 or any(value % 2 for value in row):
        raise ValueError("expected 24 even raw coordinates")
    return sum(((value // 2) & 1) << index for index, value in enumerate(row))


# The eleven octad rows in the catalogue basis, together with the all-ones
# word obtained by doubling the shifted-coset basis vector, generate C.
CATALOG_GOLAY_GENERATORS: tuple[int, ...] = tuple(
    _raw_word(row)
    for row in ODD_LEECH_BASIS_SCALED
    if all(value % 2 == 0 for value in row)
    and any(value % 4 == 2 for value in row)
) + ((1 << 24) - 1,)


@lru_cache(maxsize=1)
def catalog_golay_codewords() -> tuple[int, ...]:
    """Enumerate the exact [24,12,8] code used by the catalogue basis."""

    if len(CATALOG_GOLAY_GENERATORS) != 12:
        raise AssertionError("the catalogue Golay generator count drifted")
    words = []
    for message in range(1 << 12):
        word = 0
        for bit, generator in enumerate(CATALOG_GOLAY_GENERATORS):
            if message & (1 << bit):
                word ^= generator
        words.append(word)
    words = sorted(set(words))
    if len(words) != 4096:
        raise AssertionError("the catalogue code lost dimension twelve")
    distribution = Counter(word.bit_count() for word in words)
    if distribution != Counter({0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}):
        raise AssertionError("the catalogue code lost the Golay weight enumerator")
    if any((left & right).bit_count() % 2 for left in words for right in CATALOG_GOLAY_GENERATORS):
        raise AssertionError("the catalogue code ceased to be self-orthogonal")
    return tuple(words)


@lru_cache(maxsize=1)
def _catalog_golay_codeword_set() -> frozenset[int]:
    return frozenset(catalog_golay_codewords())


def odd_leech_raw_member(vector: Sequence[int]) -> bool:
    """Test membership in the two-coset odd-Leech Construction-B model."""

    if len(vector) != 24 or any(not isinstance(value, int) for value in vector):
        return False
    parity = vector[0] & 1
    if any((value & 1) != parity for value in vector):
        return False
    if parity == 0:
        codeword = sum(((value // 2) & 1) << i for i, value in enumerate(vector))
        residual = tuple(
            (value - 2 * ((codeword >> i) & 1)) // 4
            for i, value in enumerate(vector)
        )
        if any(
            value != 2 * ((codeword >> i) & 1) + 4 * residual[i]
            for i, value in enumerate(vector)
        ):
            return False
    else:
        codeword = sum(
            (((value - 1) // 2) & 1) << i for i, value in enumerate(vector)
        )
        residual = tuple(
            (value - 1 - 2 * ((codeword >> i) & 1)) // 4
            for i, value in enumerate(vector)
        )
        if any(
            value != 1 + 2 * ((codeword >> i) & 1) + 4 * residual[i]
            for i, value in enumerate(vector)
        ):
            return False
    return codeword in _catalog_golay_codeword_set() and sum(residual) % 2 == 0


def odd_leech_pair_scaled(left: Sequence[int], right: Sequence[int]) -> int:
    """Pair raw sqrt(8)-scaled odd-Leech vectors exactly."""

    if len(left) != 24 or len(right) != 24:
        raise ValueError("odd-Leech raw vectors have length 24")
    numerator = sum(x * y for x, y in zip(left, right))
    if numerator % 8:
        raise ValueError("the raw pairing is not integral after division by eight")
    return numerator // 8


ODD_LEECH_GRAM: tuple[tuple[int, ...], ...] = tuple(
    tuple(odd_leech_pair_scaled(left, right) for right in ODD_LEECH_BASIS_SCALED)
    for left in ODD_LEECH_BASIS_SCALED
)


def bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Exact fraction-free determinant for a small integral matrix."""

    work = [list(row) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    if not size:
        return 1
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * pivot_value
                    - work[row][column] * work[column][index]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = pivot_value
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


@lru_cache(maxsize=1)
def odd_leech_leading_principal_minors() -> tuple[int, ...]:
    return tuple(
        bareiss_determinant(tuple(row[:size] for row in ODD_LEECH_GRAM[:size]))
        for size in range(1, 25)
    )


def iter_minimal_odd_leech_vectors_scaled() -> Iterator[RawVector]:
    """Yield the complete norm-three shell as proven Golay sign vectors."""

    for codeword in catalog_golay_codewords():
        vector = tuple(
            -1 if codeword & (1 << index) else 1 for index in range(24)
        )
        if not odd_leech_raw_member(vector):
            raise AssertionError("a Golay sign vector left the shifted coset")
        if odd_leech_pair_scaled(vector, vector) != 3:
            raise AssertionError("a Golay sign vector lost norm three")
        yield vector


ODD_LEECH_NORM4_WITNESS_SCALED = ODD_LEECH_BASIS_SCALED[1]


@dataclass(frozen=True, slots=True)
class OddLeechModelAudit:
    determinant: int
    leading_principal_minors: tuple[int, ...]
    golay_weight_distribution: tuple[tuple[int, int], ...]
    minimal_shell_count: int
    minimal_shell_digest: str
    witness_square: int
    witness_pairing_histogram: tuple[tuple[int, int], ...]


@lru_cache(maxsize=1)
def odd_leech_model_audit() -> OddLeechModelAudit:
    """Audit the basis, construction, minimum, and complete minimal shell."""

    if not all(odd_leech_raw_member(row) for row in ODD_LEECH_BASIS_SCALED):
        raise AssertionError("the catalogue basis left the Construction-B model")
    determinant = bareiss_determinant(ODD_LEECH_GRAM)
    minors = odd_leech_leading_principal_minors()
    if determinant != 1 or any(value <= 0 for value in minors):
        raise AssertionError("the odd-Leech Gram matrix is not positive unimodular")

    digest = sha256()
    histogram: Counter[int] = Counter()
    total = 0
    for vector in iter_minimal_odd_leech_vectors_scaled():
        total += 1
        digest.update(bytes(value + 1 for value in vector))
        histogram[odd_leech_pair_scaled(vector, ODD_LEECH_NORM4_WITNESS_SCALED)] += 1
    if total != 4096:
        raise AssertionError("the odd-Leech minimal-shell count drifted")
    if set(histogram) - {-1, 0, 1}:
        raise AssertionError("the norm-four witness acquired a large pairing")
    return OddLeechModelAudit(
        determinant=determinant,
        leading_principal_minors=minors,
        golay_weight_distribution=tuple(
            sorted(Counter(word.bit_count() for word in catalog_golay_codewords()).items())
        ),
        minimal_shell_count=total,
        minimal_shell_digest=digest.hexdigest(),
        witness_square=odd_leech_pair_scaled(
            ODD_LEECH_NORM4_WITNESS_SCALED,
            ODD_LEECH_NORM4_WITNESS_SCALED,
        ),
        witness_pairing_histogram=tuple(sorted(histogram.items())),
    )


@dataclass(frozen=True, slots=True)
class OddLeechGramTriple:
    u_square: int
    v_square: int
    u_dot_v: int


def necessary_odd_leech_gram_triples(target: Silver) -> tuple[OddLeechGramTriple, ...]:
    """Necessary binary Gram data using integrality, minimum three, and Cauchy."""

    if not target.is_integral:
        return ()
    rational = int(target.a)
    radical = int(target.b)
    if radical % 2 or rational < 0:
        return ()
    dot = radical // 2
    output = []
    for v_square in range(rational // 2 + 1):
        u_square = rational - 2 * v_square
        if u_square not in (0,) and u_square < 3:
            continue
        if v_square not in (0,) and v_square < 3:
            continue
        if dot * dot > u_square * v_square:
            continue
        output.append(OddLeechGramTriple(u_square, v_square, dot))
    return tuple(output)


@dataclass(frozen=True, slots=True)
class OddLeechControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    spatial_square_target: Silver
    p23_spatial_gate: bool
    reduced_target: Silver
    necessary_gram_triples: tuple[OddLeechGramTriple, ...]


FIRST_PACKET_HEIGHT = Silver(7, 6) / 2
FIRST_PACKET_CONTROL = EPSILON
FIRST_PACKET_NORM = Silver(2)
FIRST_PACKET_TARGET = Silver(9, 6)
HEIGHT_ZERO_CANONICAL_NORMS: tuple[Silver, ...] = ODD_LEECH_MIRROR_NORMS[:3]
HEIGHT_ZERO_TRACE_MINIMUM = 6
H45_ENDPOINT = Silver(45, 32) / 2


@lru_cache(maxsize=1)
def odd_leech_controls_through_first_packet() -> tuple[OddLeechControlLevel, ...]:
    """Exhaust normalized controls through ``7/2+3*sqrt(2)``."""

    output = set()
    # The exact endpoint is below eight at the distinguished embedding.
    # The Minkowski box at integer height eight is therefore a safe complete
    # superset; the final exact comparison removes its overrun.
    for norm in ODD_LEECH_MIRROR_NORMS:
        caps = (norm * 8 / ODD_LEECH_COEFFICIENT, norm / (-ODD_LEECH_COEFFICIENT))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * ODD_LEECH_COEFFICIENT * control / norm).is_integral:
                continue
            target = norm + ODD_LEECH_COEFFICIENT * control * control
            if target.sign(1) < 0:
                continue
            height = ODD_LEECH_COEFFICIENT * control * control / norm
            if (FIRST_PACKET_HEIGHT - height).sign(0) < 0:
                continue
            gated = norm in P23_MIRROR_NORMS
            reduced = target / (ODD_LEECH_COEFFICIENT ** 2) if gated else target
            output.add(
                OddLeechControlLevel(
                    norm=norm,
                    height=height,
                    control=control,
                    spatial_square_target=target,
                    p23_spatial_gate=gated,
                    reduced_target=reduced,
                    necessary_gram_triples=necessary_odd_leech_gram_triples(reduced),
                )
            )
    return tuple(sorted(output, key=lambda level: (level.height, level.norm, level.control)))


def _binary_minimum(triple: OddLeechGramTriple) -> int:
    """Exact minimum of a positive-semidefinite binary integral form."""

    aa = triple.u_square
    cc = triple.v_square
    bb = triple.u_dot_v
    if aa < 0 or cc < 0 or bb * bb > aa * cc:
        raise ValueError("expected a positive-semidefinite Gram triple")
    if aa == 0:
        return cc
    if cc == 0:
        return aa
    # Gauss-reduce the Gram [[aa,bb],[bb,cc]].  Integral basis operations
    # preserve the represented norms.  In a reduced binary positive form the
    # smaller diagonal is the least nonzero represented value.  The same
    # loop handles rank one unless it discovers an integral null vector.
    while True:
        if aa > cc:
            aa, cc = cc, aa
        if aa == 0:
            return cc
        quotient = (bb + aa // 2) // aa
        if quotient:
            cc, bb = cc - 2 * quotient * bb + quotient * quotient * aa, bb - quotient * aa
            continue
        if abs(2 * bb) > aa:
            quotient = 1 if bb > 0 else -1
            cc, bb = cc - 2 * quotient * bb + quotient * quotient * aa, bb - quotient * aa
            continue
        if cc == 0:
            return 0
        return aa


def _realizable_minimum_three_triples(
    target: Silver,
) -> tuple[OddLeechGramTriple, ...]:
    """Apply every universal binary consequence of odd-Leech minimum three."""

    output = []
    for triple in necessary_odd_leech_gram_triples(target):
        determinant = (
            triple.u_square * triple.v_square - triple.u_dot_v * triple.u_dot_v
        )
        minimum = _binary_minimum(triple)
        if minimum < 3:
            continue
        if determinant == 0:
            # A rank-one pair is u=m*d, v=n*d with primitive integers m,n.
            # Dividing all three Gram entries by gcd gives d^2.  The gcd is
            # the square of gcd(m,n) times d^2, hence it is precisely d^2
            # after the primitive reduction and must be at least three.
            from math import gcd

            base_square = gcd(
                triple.u_square,
                gcd(triple.v_square, abs(triple.u_dot_v)),
            )
            if base_square < 3:
                continue
        output.append(triple)
    return tuple(output)


@lru_cache(maxsize=1)
def odd_leech_controls_through_h45() -> tuple[OddLeechControlLevel, ...]:
    """Complete normalized control/necessary-Gram census through H45."""

    output = set()
    # H45 is below 46 at the distinguished embedding.  Enumerate the safe
    # integer cap and then enforce the exact endpoint.
    for norm in ODD_LEECH_MIRROR_NORMS:
        caps = (norm * 46 / ODD_LEECH_COEFFICIENT, norm / (-ODD_LEECH_COEFFICIENT))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * ODD_LEECH_COEFFICIENT * control / norm).is_integral:
                continue
            target = norm + ODD_LEECH_COEFFICIENT * control * control
            if target.sign(1) < 0:
                continue
            height = ODD_LEECH_COEFFICIENT * control * control / norm
            if (H45_ENDPOINT - height).sign(0) < 0:
                continue
            gated = norm in P23_MIRROR_NORMS
            reduced = target / (ODD_LEECH_COEFFICIENT ** 2) if gated else target
            output.add(
                OddLeechControlLevel(
                    norm=norm,
                    height=height,
                    control=control,
                    spatial_square_target=target,
                    p23_spatial_gate=gated,
                    reduced_target=reduced,
                    necessary_gram_triples=_realizable_minimum_three_triples(reduced),
                )
            )
    return tuple(sorted(output, key=lambda level: (level.height, level.norm, level.control)))


@dataclass(frozen=True, slots=True)
class OddLeechH45GramBound:
    original_triple: OddLeechGramTriple
    transformed_description: str
    fixed_direction_maximum_pairing: Silver
    wall_bound: Silver
    strict_margin: Silver


@dataclass(frozen=True, slots=True)
class OddLeechH45PrefixAudit:
    controls: tuple[OddLeechControlLevel, ...]
    nonempty_necessary_levels: tuple[OddLeechControlLevel, ...]
    terminal_height: Silver
    terminal_norm: Silver
    terminal_control: Silver
    terminal_target: Silver
    terminal_triples: tuple[OddLeechGramTriple, ...]
    fixed_direction_bounds: tuple[OddLeechH45GramBound, ...]
    witness_scale: Fraction
    witness_square: int
    minimum_wall_margin: Silver
    hyperideal_excess: Silver


def psd_extension_joint_maximum(
    first_square: int,
    second_square: int,
    first_dot_second: int,
    fixed_square: int,
    second_coefficient: Silver,
) -> tuple[int, int, Silver]:
    """Maximize ``m+c*n`` over all integral positive-semidefinite extensions.

    Here ``m`` and ``n`` are the pairings of the fixed direction with the
    first and second binary vectors.  The two 2x2 principal-minor bounds make
    the enumeration finite, and the displayed determinant formula is exact.
    """

    if (
        min(first_square, second_square, fixed_square) < 0
        or first_dot_second * first_dot_second > first_square * second_square
    ):
        raise ValueError("expected a positive-semidefinite binary Gram matrix")
    best: tuple[int, int, Silver] | None = None
    for first_pairing in range(
        -isqrt(first_square * fixed_square),
        isqrt(first_square * fixed_square) + 1,
    ):
        for second_pairing in range(
            -isqrt(second_square * fixed_square),
            isqrt(second_square * fixed_square) + 1,
        ):
            determinant = (
                first_square * second_square * fixed_square
                + 2 * first_dot_second * first_pairing * second_pairing
                - first_square * second_pairing * second_pairing
                - second_square * first_pairing * first_pairing
                - fixed_square * first_dot_second * first_dot_second
            )
            if determinant < 0:
                continue
            value = Silver(first_pairing) + second_coefficient * second_pairing
            if best is None or (value - best[2]).sign(0) > 0:
                best = first_pairing, second_pairing, value
    if best is None:
        raise AssertionError("a PSD binary Gram matrix lost all extensions")
    return best


@lru_cache(maxsize=1)
def odd_leech_h45_prefix_audit() -> OddLeechH45PrefixAudit:
    """Certify a scaled norm-four direction through every possible H45 wall.

    The terminal packet may contain billions of ordinary shell pairs.  No
    realization expansion is needed: positive-semidefiniteness of the 3x3
    Gram matrix with the fixed norm-four direction gives rigorous pairing
    bounds for a superset of every actual odd-Leech pair.
    """

    controls = odd_leech_controls_through_h45()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    if len(controls) != 22 or len(nonempty) != 2:
        raise AssertionError("the H45 control census drifted")
    first, terminal = nonempty
    if (
        first.height != FIRST_PACKET_HEIGHT
        or terminal.height != H45_ENDPOINT
        or terminal.norm != Silver(2)
        or terminal.control != Silver(3, 2)
        or terminal.reduced_target != Silver(47, 32)
    ):
        raise AssertionError("the H45 endpoints changed")
    triples = (
        OddLeechGramTriple(27, 10, 16),
        OddLeechGramTriple(25, 11, 16),
        OddLeechGramTriple(23, 12, 16),
    )
    if terminal.necessary_gram_triples != triples:
        raise AssertionError("the H45 necessary Gram types changed")

    wall_bound = ODD_LEECH_COEFFICIENT * terminal.control
    # T1: put b=u-2v (b^2=3), w=v (w^2=10), b.w=-4, so
    # y=b+(2+sqrt2)w.  A complete integer scan of PSD Gram(b,w,d)
    # has joint maximum at (b.d,w.d)=(-2,6).
    # T2: put d'=u-v (d'^2=4), v^2=11, d'.v=5.  The analogous complete
    # scan has joint maximum (d'.d,v.d)=(3,6).
    # T3: put b=u-v (b^2=3), w=v (w^2=12), b.w=4.  Its joint maximum is
    # again (b.d,w.d)=(3,6).
    #
    # These are bounds for d itself.  The strict hyperideal witness through
    # H45 is z=(5/6)d; scale the maxima before comparing with the walls.
    witness_scale = Fraction(5, 6)
    transformed_grams = (
        (3, 10, -4, Silver(2, 1)),
        (4, 11, 5, EPSILON),
        (3, 12, 4, EPSILON),
    )
    joint_maxima = tuple(
        psd_extension_joint_maximum(*gram[:3], 4, gram[3])
        for gram in transformed_grams
    )
    if joint_maxima != (
        (-2, 6, Silver(10, 6)),
        (3, 6, Silver(9, 6)),
        (3, 6, Silver(9, 6)),
    ):
        raise AssertionError("the H45 joint PSD maxima drifted")
    maxima = tuple(witness_scale * maximum for _m, _n, maximum in joint_maxima)
    descriptions = (
        "b=u-2v: (b^2,w^2,b.w)=(3,10,-4), y=b+(2+sqrt(2))*w; complete integer 3x3 PSD scan gives maximum y.d=10+6sqrt(2), then z=(5/6)d",
        "d'=u-v: (d'^2,v^2,d'.v)=(4,11,5), y=d'+(1+sqrt(2))*v; complete integer 3x3 PSD scan gives maximum y.d=9+6sqrt(2), then z=(5/6)d",
        "b=u-v: (b^2,w^2,b.w)=(3,12,4), y=b+(1+sqrt(2))*w; complete integer 3x3 PSD scan gives maximum y.d=9+6sqrt(2), then z=(5/6)d",
    )
    bounds = tuple(
        OddLeechH45GramBound(
            original_triple=triple,
            transformed_description=description,
            fixed_direction_maximum_pairing=maximum,
            wall_bound=wall_bound,
            strict_margin=wall_bound - maximum,
        )
        for triple, description, maximum in zip(triples, descriptions, maxima)
    )
    if any(bound.strict_margin.sign(0) <= 0 for bound in bounds):
        raise AssertionError("the H45 fixed direction lost strict clearance")
    first_audit = odd_leech_first_packet_audit()
    first_scaled_margin = (
        ODD_LEECH_COEFFICIENT
        - witness_scale * first_audit.witness_maximum_pairing
    )
    minimum_margin = min(
        (first_scaled_margin, *(bound.strict_margin for bound in bounds))
    )
    scaled_excess = (
        Silver(first_audit.witness_square) * witness_scale * witness_scale
        - ODD_LEECH_COEFFICIENT
    )
    if scaled_excess.sign(0) <= 0:
        raise AssertionError("the scaled H45 witness ceased to be hyperideal")
    return OddLeechH45PrefixAudit(
        controls=controls,
        nonempty_necessary_levels=nonempty,
        terminal_height=terminal.height,
        terminal_norm=terminal.norm,
        terminal_control=terminal.control,
        terminal_target=terminal.reduced_target,
        terminal_triples=triples,
        fixed_direction_bounds=bounds,
        witness_scale=witness_scale,
        witness_square=first_audit.witness_square,
        minimum_wall_margin=minimum_margin,
        hyperideal_excess=scaled_excess,
    )


@dataclass(frozen=True, slots=True)
class OddLeechFirstPacketAudit:
    model: OddLeechModelAudit
    mirror_norms: tuple[Silver, ...]
    height_zero_wall_count: int
    controls: tuple[OddLeechControlLevel, ...]
    first_height: Silver
    first_norm: Silver
    first_control: Silver
    first_target: Silver
    first_gram_triple: OddLeechGramTriple
    first_packet_root_count: int
    distinct_spatial_pairings: tuple[int, ...]
    distinct_root_pairings: tuple[Silver, ...]
    witness_square: int
    witness_maximum_pairing: int
    witness_minimum_wall_margin: Silver
    witness_hyperideal_excess: Silver


@lru_cache(maxsize=1)
def odd_leech_first_packet_audit() -> OddLeechFirstPacketAudit:
    """Certify the complete first packet and a strict hyperideal witness."""

    model = odd_leech_model_audit()
    # At height zero, a p23-supported norm would force y in p*O24 and hence
    # contradict primitivity.  The remaining normalized norms have traces
    # 2,4,4.  But for y=u+sqrt(2)v, Tr(y^2)=2(u^2+2v^2)>=6 by the odd-Leech
    # minimum, so no height-zero mirror exists.
    if max(int(norm.trace()) for norm in HEIGHT_ZERO_CANONICAL_NORMS) >= HEIGHT_ZERO_TRACE_MINIMUM:
        raise AssertionError("the height-zero trace exclusion stopped being strict")
    controls = odd_leech_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    if len(controls) != 9 or len(nonempty) != 1:
        raise AssertionError("the exact first-prefix control census drifted")
    first = nonempty[0]
    if (
        first.height != FIRST_PACKET_HEIGHT
        or first.norm != FIRST_PACKET_NORM
        or first.control != FIRST_PACKET_CONTROL
        or first.spatial_square_target != FIRST_PACKET_TARGET
        or first.necessary_gram_triples != (OddLeechGramTriple(3, 3, 3),)
    ):
        raise AssertionError("the first odd-Leech packet changed")

    # Distinct codewords differ by a nonzero Golay word of weight
    # 8,12,16,24.  Their sign-vector pairings are 3-weight/4.
    spatial_pairings = (1, 0, -1, -3)
    root_pairings = tuple(
        EPSILON ** 2 * (-ODD_LEECH_COEFFICIENT + pairing)
        for pairing in spatial_pairings
    )
    if any(value.sign(0) >= 0 for value in root_pairings):
        raise AssertionError("the first packet ceased to be Vinberg-acute")

    maximum = max(pairing for pairing, _count in model.witness_pairing_histogram)
    margin = ODD_LEECH_COEFFICIENT - maximum
    excess = Silver(model.witness_square) - ODD_LEECH_COEFFICIENT
    if margin.sign(0) <= 0 or excess.sign(0) <= 0:
        raise AssertionError("the odd-Leech prefix lost its hyperideal witness")
    return OddLeechFirstPacketAudit(
        model=model,
        mirror_norms=ODD_LEECH_MIRROR_NORMS,
        height_zero_wall_count=0,
        controls=controls,
        first_height=first.height,
        first_norm=first.norm,
        first_control=first.control,
        first_target=first.spatial_square_target,
        first_gram_triple=first.necessary_gram_triples[0],
        first_packet_root_count=model.minimal_shell_count,
        distinct_spatial_pairings=spatial_pairings,
        distinct_root_pairings=root_pairings,
        witness_square=model.witness_square,
        witness_maximum_pairing=maximum,
        witness_minimum_wall_margin=margin,
        witness_hyperideal_excess=excess,
    )


__all__ = [
    "CATALOG_GOLAY_GENERATORS",
    "EPSILON",
    "FIRST_PACKET_CONTROL",
    "FIRST_PACKET_HEIGHT",
    "FIRST_PACKET_NORM",
    "FIRST_PACKET_TARGET",
    "H45_ENDPOINT",
    "HEIGHT_ZERO_CANONICAL_NORMS",
    "HEIGHT_ZERO_TRACE_MINIMUM",
    "ODD_LEECH_BASIS_SCALED",
    "ODD_LEECH_COEFFICIENT",
    "ODD_LEECH_GRAM",
    "ODD_LEECH_MIRROR_NORMS",
    "ODD_LEECH_NORM4_WITNESS_SCALED",
    "ODD_LEECH_NORM_GATE_PROOF",
    "OddLeechControlLevel",
    "OddLeechFirstPacketAudit",
    "OddLeechGramTriple",
    "OddLeechH45GramBound",
    "OddLeechH45PrefixAudit",
    "OddLeechModelAudit",
    "bareiss_determinant",
    "catalog_golay_codewords",
    "iter_minimal_odd_leech_vectors_scaled",
    "necessary_odd_leech_gram_triples",
    "odd_leech_controls_through_first_packet",
    "odd_leech_first_packet_audit",
    "odd_leech_h45_prefix_audit",
    "odd_leech_leading_principal_minors",
    "odd_leech_model_audit",
    "odd_leech_controls_through_h45",
    "odd_leech_pair_scaled",
    "odd_leech_raw_member",
    "psd_extension_joint_maximum",
]
