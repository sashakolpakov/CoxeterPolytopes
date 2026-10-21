"""Exact full three-node AGCat cycle packets over ``Q(sqrt(5))``.

The universal golden frame alphabet has 6,875 length-two products, 6,869
of them dotted.  A direct cubic enumeration of its three-edge multisets is
unnecessary.  This module groups the 6,874 nonzero values by principal-ideal
squareclass, joins three buckets by XOR, and streams precisely the products
whose square is a golden integer square.

For an acute symmetrizable Cartan triangle with edge products ``c1,c2,c3``,
the directed three-cycle is the negative identity-positive root ``b`` of

``b^2=c1*c2*c3``.

The normalized Gram determinant is

``8 - 2*(c1+c2+c3) + 2*b``.

At the conjugate place every independent frame triple is positive definite,
so its determinant must be positive.  At the selected place the independent
triple can be positive definite, degenerate, or Lorentzian.  Degenerate
packets stay in the universal census; compactness removes affine/parabolic
subdiagrams because a cocompact reflection group has no parabolic elements.
The result here is a complete local three-node packet tranche, not a
compatible high-rank packet enumeration.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
from math import comb
from typing import Iterator

from coxeter_polytopes.classification.golden_agcat_manifest import (
    golden_coxeter_frame_pair_alphabet,
)
from coxeter_polytopes.exact.golden import ONE, PHI, Golden


PrimeIdealLabel = tuple[int, int]
SquareclassKey = tuple[PrimeIdealLabel, ...]


def _prime_factors(value: int) -> tuple[tuple[int, int], ...]:
    """Factor a positive rational integer by exact trial division."""

    if value <= 0:
        raise ValueError("expected a positive integer")
    output = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor:
            divisor = 3 if divisor == 2 else divisor + 2
            continue
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        output.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        output.append((value, 1))
    return tuple(output)


def _roots_mod_prime(prime: int) -> tuple[int, ...]:
    return tuple(
        residue
        for residue in range(prime)
        if (residue * residue - residue - 1) % prime == 0
    )


def _hensel_root(residue: int, prime: int, exponent: int) -> int:
    """Lift a simple root of ``X^2-X-1`` modulo ``p`` to ``p^exponent``."""

    modulus = prime
    for _ in range(1, exponent):
        quotient = (residue * residue - residue - 1) // modulus
        correction = (
            -quotient * pow((2 * residue - 1) % prime, -1, prime)
        ) % prime
        residue += correction * modulus
        modulus *= prime
    if (residue * residue - residue - 1) % modulus:
        raise AssertionError("Hensel lift ceased to be a polynomial root")
    return residue


def _valuation(value: int, prime: int, cap: int) -> int:
    if not value:
        return cap
    exponent = 0
    while exponent < cap and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


@dataclass(frozen=True, slots=True, order=True)
class GoldenElementSquareclass:
    ideal_parity: SquareclassKey
    conjugate_sign: int
    unit_parity: int


@lru_cache(maxsize=None)
def golden_element_squareclass(value: Golden) -> GoldenElementSquareclass:
    """Return the element squareclass of a positive golden integer.

    The label ``(p,0)`` denotes the unique inert or ramified prime above
    ``p``.  For a split prime, ``(p,r)`` denotes the prime at which
    ``phi`` reduces to the smaller or larger polynomial root ``r`` modulo
    ``p``.  Hensel lifting extracts exact valuations even when the rational
    norm has repeated factors.
    """

    if not value.is_integral or not value:
        raise ValueError("squareclass factorization needs a nonzero golden integer")
    field_norm = value.norm()
    if field_norm.denominator != 1:
        raise AssertionError("an integral golden element acquired fractional norm")
    absolute_norm = abs(int(field_norm))
    rational, phi_coefficient = int(value.a), int(value.b)
    parity: list[PrimeIdealLabel] = []
    for prime, norm_exponent in _prime_factors(absolute_norm):
        if prime == 5:
            ramified_generator = 2 + PHI
            quotient = value
            valuation = 0
            while (quotient / ramified_generator).is_integral:
                quotient /= ramified_generator
                valuation += 1
            if valuation != norm_exponent:
                raise AssertionError("the ramified valuation disagrees with the norm")
            if valuation % 2:
                parity.append((prime, 0))
            continue

        roots = _roots_mod_prime(prime)
        if not roots:
            valuation = min(
                _valuation(rational, prime, norm_exponent),
                _valuation(phi_coefficient, prime, norm_exponent),
            )
            if 2 * valuation != norm_exponent:
                raise AssertionError("an inert valuation disagrees with the norm")
            if valuation % 2:
                parity.append((prime, 0))
            continue

        valuations = []
        for residue in roots:
            lifted = _hensel_root(residue, prime, norm_exponent)
            valuation = _valuation(
                rational + phi_coefficient * lifted,
                prime,
                norm_exponent,
            )
            valuations.append(valuation)
            if valuation % 2:
                parity.append((prime, residue))
        if sum(valuations) != norm_exponent:
            raise AssertionError("split-prime valuations disagree with the norm")
    ideal_parity = tuple(sorted(parity))

    # The unit group is +/-phi^Z.  Once ideal parity and the two real signs
    # are recorded, no further squareclass bit remains: -1 is detected by
    # the selected sign (always positive in this alphabet), and phi parity is
    # detected by the norm sign, equivalently the conjugate sign here.
    return GoldenElementSquareclass(
        ideal_parity=ideal_parity,
        conjugate_sign=value.sign(1),
        unit_parity=0,
    )


def _golden_square_root_coefficients(value: Golden) -> Golden | None:
    """Return the identity-positive golden square root, if it exists."""

    from math import isqrt

    if not value.is_integral or value.sign(0) <= 0 or value.sign(1) <= 0:
        return None
    a, b = int(value.a), int(value.b)
    coefficient = 4 * a + 2 * b
    discriminant = coefficient * coefficient - 20 * b * b
    if discriminant < 0:
        return None
    discriminant_root = isqrt(discriminant)
    if discriminant_root * discriminant_root != discriminant:
        return None
    candidates = []
    # From A=m^2+n^2 and B=2mn+n^2, the integer t=n^2 solves
    # 5t^2-(4A+2B)t+B^2=0.
    for t_numerator in (
        coefficient - discriminant_root,
        coefficient + discriminant_root,
    ):
        if t_numerator < 0 or t_numerator % 10:
            continue
        n_square = t_numerator // 10
        n_abs = isqrt(n_square)
        if n_abs * n_abs != n_square:
            continue
        for n in ({n_abs, -n_abs} if n_abs else {0}):
            if n:
                m_numerator = b - n * n
                if m_numerator % (2 * n):
                    continue
                candidates.append(Golden(m_numerator // (2 * n), n))
            elif b == 0:
                m_abs = isqrt(a)
                if m_abs * m_abs == a:
                    candidates.extend((Golden(m_abs), Golden(-m_abs)))
    roots = tuple(sorted(set(root for root in candidates if root * root == value)))
    positive = tuple(root for root in roots if root.sign(0) > 0)
    return positive[0] if len(roots) == 2 and roots[0] == -roots[1] and len(positive) == 1 else None


def _squareclass_xor(
    left: GoldenElementSquareclass, right: GoldenElementSquareclass
) -> GoldenElementSquareclass:
    return GoldenElementSquareclass(
        tuple(sorted(set(left.ideal_parity) ^ set(right.ideal_parity))),
        left.conjugate_sign * right.conjugate_sign,
        left.unit_parity ^ right.unit_parity,
    )


def golden_square_root_from_squareclass_product(
    left: Golden, middle: Golden, right: Golden
) -> Golden:
    """Return the positive golden square root after the squareclass gate."""

    product = left * middle * right
    if product.sign(0) <= 0 or product.sign(1) <= 0:
        raise AssertionError("an edge product left the totally-positive cone")
    root = _golden_square_root_coefficients(product)
    if root is None:
        raise AssertionError("squareclass-compatible product lost its golden square root")
    # Total positivity of the product forces the two embeddings of a root
    # either to have equal or opposite signs.  A totally-positive unit
    # multiple cannot change this.  Cartan acuteness fixes only sigma0(root)
    # positive before taking the negative directed cycle.
    if root.sign(1) == 0:
        raise AssertionError("a nonzero field root vanished at one embedding")
    return root


@dataclass(frozen=True, slots=True)
class GoldenTriangleStratum:
    edge_kinds: tuple[str, str, str]
    root_conjugate_sign: int
    selected_determinant_sign: int
    conjugate_determinant_sign: int
    count: int


@dataclass(frozen=True, slots=True)
class GoldenTrianglePacketCensus:
    alphabet_count: int
    zero_edge_value_count: int
    nonzero_edge_value_count: int
    squareclass_bucket_count: int
    maximum_squareclass_bucket_size: int
    compatible_squareclass_bucket_triples: int
    raw_unordered_nonzero_edge_multisets: int
    square_product_edge_multisets: int
    conjugate_positive_packet_count: int
    selected_positive_packet_count: int
    selected_lorentzian_packet_count: int
    raw_selected_zero_determinant_count: int
    retained_selected_zero_determinant_count: int
    conjugate_nonpositive_packet_count: int
    strata: tuple[GoldenTriangleStratum, ...]
    canonical_stream_sha256: str
    retained_stream_sha256: str


def _edge_kind(value: Golden) -> str:
    return "finite" if (value - 4).sign(0) < 0 else "dotted"


def _encoded_golden(value: Golden) -> bytes:
    if not value.is_integral:
        raise AssertionError("triangle stream value left Z[phi]")
    return f"{int(value.a)},{int(value.b)}".encode()


def _pair_mul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def _pair_sign(value: tuple[int, int], embedding: int) -> int:
    a, b = value
    if embedding:
        a, b = a + b, -b
    rational_part = 2 * a + b
    if not b:
        return (rational_part > 0) - (rational_part < 0)
    if not rational_part:
        return (b > 0) - (b < 0)
    rational_sign = (rational_part > 0) - (rational_part < 0)
    radical_sign = (b > 0) - (b < 0)
    if rational_sign == radical_sign:
        return rational_sign
    comparison = rational_part * rational_part - 5 * b * b
    return rational_sign if comparison > 0 else radical_sign


def _pair_square_root(value: tuple[int, int]) -> tuple[int, int]:
    """Fast raw-coefficient version of the audited golden root solver."""

    from math import isqrt

    a, b = value
    coefficient = 4 * a + 2 * b
    discriminant = coefficient * coefficient - 20 * b * b
    if discriminant < 0:
        raise AssertionError("squareclass join emitted negative root discriminant")
    discriminant_root = isqrt(discriminant)
    if discriminant_root * discriminant_root != discriminant:
        raise AssertionError("squareclass join emitted nonsquare root discriminant")
    roots = []
    for numerator in (coefficient - discriminant_root, coefficient + discriminant_root):
        if numerator < 0 or numerator % 10:
            continue
        n_square = numerator // 10
        n_abs = isqrt(n_square)
        if n_abs * n_abs != n_square:
            continue
        for n in ({n_abs, -n_abs} if n_abs else {0}):
            if n:
                m_numerator = b - n * n
                if m_numerator % (2 * n):
                    continue
                root = (m_numerator // (2 * n), n)
            else:
                m_abs = isqrt(a)
                if b or m_abs * m_abs != a:
                    continue
                root = (m_abs, 0)
            if _pair_mul(root, root) == value and _pair_sign(root, 0) > 0:
                roots.append(root)
    roots = list(set(roots))
    if len(roots) != 1:
        raise AssertionError("squareclass join lost its unique identity-positive root")
    return roots[0]


@lru_cache(maxsize=1)
def golden_triangle_packet_census() -> GoldenTrianglePacketCensus:
    """Stream the complete nonorthogonal three-edge packet census."""

    alphabet = golden_coxeter_frame_pair_alphabet()
    nonzero = tuple(value for value in alphabet if value)
    buckets: dict[GoldenElementSquareclass, list[Golden]] = defaultdict(list)
    for value in nonzero:
        buckets[golden_element_squareclass(value)].append(value)
    for values in buckets.values():
        values.sort()
    raw_buckets = {
        key: tuple((int(value.a), int(value.b)) for value in values)
        for key, values in buckets.items()
    }
    keys = tuple(sorted(buckets))
    key_index = {key: index for index, key in enumerate(keys)}

    bucket_triples = []
    square_product_count = 0
    for first_index, first_key in enumerate(keys):
        first_count = len(buckets[first_key])
        for second_index in range(first_index, len(keys)):
            second_key = keys[second_index]
            third_key = _squareclass_xor(first_key, second_key)
            third_index = key_index.get(third_key, -1)
            if third_index < second_index:
                continue
            second_count = len(buckets[second_key])
            third_count = len(buckets[third_key])
            bucket_triples.append((first_index, second_index, third_index))
            if first_index == second_index == third_index:
                square_product_count += comb(first_count + 2, 3)
            elif first_index == second_index:
                square_product_count += comb(first_count + 1, 2) * third_count
            elif second_index == third_index:
                square_product_count += first_count * comb(second_count + 1, 2)
            else:
                square_product_count += first_count * second_count * third_count

    strata: Counter[tuple[tuple[str, str, str], int, int, int]] = Counter()
    canonical_digest = sha256()
    retained_digest = sha256()
    raw_selected_zero = 0
    retained_selected_zero = 0
    conjugate_nonpositive = 0
    emitted = 0
    retained = 0

    for first_index, second_index, third_index in bucket_triples:
        first_values = raw_buckets[keys[first_index]]
        second_values = raw_buckets[keys[second_index]]
        third_values = raw_buckets[keys[third_index]]
        if first_index == second_index == third_index:
            iterator = combinations_with_replacement(
                first_values, 3
            )
        elif first_index == second_index:
            iterator = (
                (left, middle, right)
                for left, middle in combinations_with_replacement(first_values, 2)
                for right in third_values
            )
        elif second_index == third_index:
            iterator = (
                (left, middle, right)
                for left in first_values
                for middle, right in combinations_with_replacement(second_values, 2)
            )
        else:
            iterator = (
                (left, middle, right)
                for left in first_values
                for middle in second_values
                for right in third_values
            )

        for raw_edges in iterator:
            edges = tuple(sorted(raw_edges))
            product = _pair_mul(_pair_mul(edges[0], edges[1]), edges[2])
            root = _pair_square_root(product)
            cycle = (-root[0], -root[1])
            determinant = (
                8 - 2 * sum(edge[0] for edge in edges) + 2 * cycle[0],
                -2 * sum(edge[1] for edge in edges) + 2 * cycle[1],
            )
            selected_sign = _pair_sign(determinant, 0)
            conjugate_sign = _pair_sign(determinant, 1)
            root_conjugate_sign = _pair_sign(root, 1)
            edge_kinds = tuple(
                sorted(
                    "finite" if _pair_sign((edge[0] - 4, edge[1]), 0) < 0 else "dotted"
                    for edge in edges
                )
            )
            record = (
                ";".join(
                    f"{value[0]},{value[1]}"
                    for value in (*edges, cycle, determinant)
                )
                + "\n"
            ).encode()
            canonical_digest.update(record)
            emitted += 1
            if selected_sign == 0:
                raw_selected_zero += 1
                if conjugate_sign > 0:
                    retained_selected_zero += 1
            if conjugate_sign <= 0:
                conjugate_nonpositive += 1
            else:
                retained_digest.update(record)
                retained += 1
            strata[(edge_kinds, root_conjugate_sign, selected_sign, conjugate_sign)] += 1

    if emitted != square_product_count:
        raise AssertionError("the streamed square-product count changed")
    if retained + conjugate_nonpositive != emitted:
        raise AssertionError("the conjugate determinant partition changed")

    strata_rows = tuple(
        GoldenTriangleStratum(*key, count)
        for key, count in sorted(strata.items())
    )
    selected_positive = sum(
        row.count
        for row in strata_rows
        if row.conjugate_determinant_sign > 0 and row.selected_determinant_sign > 0
    )
    selected_lorentzian = sum(
        row.count
        for row in strata_rows
        if row.conjugate_determinant_sign > 0 and row.selected_determinant_sign < 0
    )
    return GoldenTrianglePacketCensus(
        alphabet_count=len(alphabet),
        zero_edge_value_count=len(alphabet) - len(nonzero),
        nonzero_edge_value_count=len(nonzero),
        squareclass_bucket_count=len(buckets),
        maximum_squareclass_bucket_size=max(map(len, buckets.values())),
        compatible_squareclass_bucket_triples=len(bucket_triples),
        raw_unordered_nonzero_edge_multisets=comb(len(nonzero) + 2, 3),
        square_product_edge_multisets=square_product_count,
        conjugate_positive_packet_count=retained,
        selected_positive_packet_count=selected_positive,
        selected_lorentzian_packet_count=selected_lorentzian,
        raw_selected_zero_determinant_count=raw_selected_zero,
        retained_selected_zero_determinant_count=retained_selected_zero,
        conjugate_nonpositive_packet_count=conjugate_nonpositive,
        strata=strata_rows,
        canonical_stream_sha256=canonical_digest.hexdigest(),
        retained_stream_sha256=retained_digest.hexdigest(),
    )


__all__ = [name for name in globals() if not name.startswith("_")]
