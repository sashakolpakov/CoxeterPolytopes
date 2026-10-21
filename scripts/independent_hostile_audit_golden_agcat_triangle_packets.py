#!/usr/bin/env python3
"""Independent hostile audit of the full golden AGCat triangle tranche.

The length-two alphabet is treated as the declared input.  This replay does
not import ``golden_triangle_packets``: it independently implements golden
integer arithmetic, prime-ideal parity, square-root reconstruction, bucket
joining, determinant signs, multiplicities, strata, and both stream hashes.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from functools import cmp_to_key
from hashlib import sha256
from itertools import combinations_with_replacement
from math import comb, isqrt
from pathlib import Path

from coxeter_polytopes.classification.golden_agcat_manifest import (
    golden_coxeter_frame_pair_alphabet,
)


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = ROOT / "data/enumerations/golden_agcat_triangle_packets.json"
OUTPUT = (
    ROOT
    / "data/enumerations/golden_agcat_triangle_packets_independent_hostile_audit.json"
)

Pair = tuple[int, int]
PrimeLabel = tuple[int, int]
Squareclass = tuple[tuple[PrimeLabel, ...], int]


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def multiply(left: Pair, right: Pair) -> Pair:
    """Multiply ``a+b*phi`` pairs using ``phi^2=phi+1``."""

    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def norm(value: Pair) -> int:
    a, b = value
    return a * a + a * b - b * b


def sign(value: Pair, embedding: int = 0) -> int:
    a, b = value
    if embedding == 1:
        a, b = a + b, -b
    elif embedding != 0:
        raise ValueError("embedding must be zero or one")
    rational = 2 * a + b
    sr = (rational > 0) - (rational < 0)
    sb = (b > 0) - (b < 0)
    if not sr:
        return sb
    if not sb or sr == sb:
        return sr
    return sr if rational * rational - 5 * b * b > 0 else sb


def compare(left: Pair, right: Pair) -> int:
    return sign((left[0] - right[0], left[1] - right[1]))


def exact_divide(left: Pair, right: Pair) -> Pair | None:
    """Divide in ``Z[phi]`` and return ``None`` if the quotient is nonintegral."""

    a, b = left
    c, d = right
    denominator = c * c + c * d - d * d
    first = a * (c + d) - b * d
    second = b * c - a * d
    if not denominator or first % denominator or second % denominator:
        return None
    return first // denominator, second // denominator


def prime_factors(value: int) -> tuple[tuple[int, int], ...]:
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


def polynomial_roots_mod_prime(prime: int) -> tuple[int, ...]:
    return tuple(
        residue
        for residue in range(prime)
        if (residue * residue - residue - 1) % prime == 0
    )


def hensel_root(residue: int, prime: int, exponent: int) -> int:
    modulus = prime
    for _ in range(1, exponent):
        quotient = (residue * residue - residue - 1) // modulus
        derivative_inverse = pow((2 * residue - 1) % prime, -1, prime)
        residue += (-quotient * derivative_inverse % prime) * modulus
        modulus *= prime
    if (residue * residue - residue - 1) % modulus:
        raise AssertionError("independent Hensel lift failed")
    return residue


def integer_valuation(value: int, prime: int, cap: int) -> int:
    if not value:
        return cap
    exponent = 0
    while exponent < cap and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def squareclass(value: Pair) -> Squareclass:
    """Principal-ideal parity plus the non-square unit/sign bit."""

    if value == (0, 0) or sign(value) <= 0:
        raise ValueError("expected an identity-positive nonzero golden integer")
    absolute_norm = abs(norm(value))
    a, b = value
    parity: list[PrimeLabel] = []
    for prime, norm_exponent in prime_factors(absolute_norm):
        if prime == 5:
            quotient = value
            valuation = 0
            while (next_quotient := exact_divide(quotient, (2, 1))) is not None:
                quotient = next_quotient
                valuation += 1
            if valuation != norm_exponent:
                raise AssertionError("ramified valuation/norm disagreement")
            if valuation % 2:
                parity.append((5, 0))
            continue

        roots = polynomial_roots_mod_prime(prime)
        if not roots:
            valuation = min(
                integer_valuation(a, prime, norm_exponent),
                integer_valuation(b, prime, norm_exponent),
            )
            if 2 * valuation != norm_exponent:
                raise AssertionError("inert valuation/norm disagreement")
            if valuation % 2:
                parity.append((prime, 0))
            continue

        valuations = []
        for residue in roots:
            lifted = hensel_root(residue, prime, norm_exponent)
            valuation = integer_valuation(a + b * lifted, prime, norm_exponent)
            valuations.append(valuation)
            if valuation % 2:
                parity.append((prime, residue))
        if sum(valuations) != norm_exponent:
            raise AssertionError("split valuation/norm disagreement")

    # All input edge values are positive at the selected place.  The
    # remaining unit squareclass is phi^Z mod squares, detected by the
    # conjugate sign because phi has signs (+,-).
    return tuple(sorted(parity)), sign(value, 1)


def xor(left: Squareclass, right: Squareclass) -> Squareclass:
    return (
        tuple(sorted(set(left[0]) ^ set(right[0]))),
        left[1] * right[1],
    )


def identity_positive_square_root(value: Pair) -> Pair:
    """Reconstruct the unique identity-positive integral root and verify it."""

    a, b = value
    coefficient = 4 * a + 2 * b
    discriminant = coefficient * coefficient - 20 * b * b
    if discriminant < 0:
        raise AssertionError("squareclass-compatible product has negative discriminant")
    delta = isqrt(discriminant)
    if delta * delta != discriminant:
        raise AssertionError("squareclass-compatible product has nonsquare discriminant")
    candidates = set()
    for numerator in (coefficient - delta, coefficient + delta):
        if numerator < 0 or numerator % 10:
            continue
        n_square = numerator // 10
        n_abs = isqrt(n_square)
        if n_abs * n_abs != n_square:
            continue
        for n in ((n_abs, -n_abs) if n_abs else (0,)):
            if n:
                m_numerator = b - n * n
                if m_numerator % (2 * n):
                    continue
                root = m_numerator // (2 * n), n
            else:
                m_abs = isqrt(a)
                if b or m_abs * m_abs != a:
                    continue
                root = m_abs, 0
            if multiply(root, root) == value and sign(root) > 0:
                candidates.add(root)
    if len(candidates) != 1:
        raise AssertionError("compatible product lost its unique positive root")
    return candidates.pop()


def expected_multiset_count(
    first_index: int,
    second_index: int,
    third_index: int,
    first_size: int,
    second_size: int,
    third_size: int,
) -> int:
    if first_index == second_index == third_index:
        return comb(first_size + 2, 3)
    if first_index == second_index:
        return comb(first_size + 1, 2) * third_size
    if second_index == third_index:
        return first_size * comb(second_size + 1, 2)
    return first_size * second_size * third_size


def replay() -> dict[str, object]:
    alphabet = tuple((int(value.a), int(value.b)) for value in golden_coxeter_frame_pair_alphabet())
    nonzero = tuple(value for value in alphabet if value != (0, 0))
    if len(alphabet) != 6_875 or len(nonzero) != 6_874:
        raise AssertionError("declared edge alphabet drifted")
    if any(sign(value) <= 0 or sign(value, 1) <= 0 for value in nonzero):
        raise AssertionError("a nonzero frame edge left the totally-positive cone")

    buckets: dict[Squareclass, list[Pair]] = defaultdict(list)
    for value in nonzero:
        buckets[squareclass(value)].append(value)
    for values in buckets.values():
        values.sort(key=cmp_to_key(compare))
    keys = tuple(sorted(buckets))
    key_index = {key: index for index, key in enumerate(keys)}

    # A representative product inside each bucket must be an actual square;
    # this is a direct computational check of the ideal/sign classification.
    for values in buckets.values():
        representative = values[0]
        for value in values:
            identity_positive_square_root(multiply(representative, value))

    bucket_triples = []
    formula_count = 0
    for first_index, first_key in enumerate(keys):
        for second_index in range(first_index, len(keys)):
            third_index = key_index.get(xor(first_key, keys[second_index]), -1)
            if third_index < second_index:
                continue
            bucket_triples.append((first_index, second_index, third_index))
            formula_count += expected_multiset_count(
                first_index,
                second_index,
                third_index,
                len(buckets[first_key]),
                len(buckets[keys[second_index]]),
                len(buckets[keys[third_index]]),
            )

    canonical = sha256()
    retained_digest = sha256()
    strata: Counter[tuple[tuple[str, str, str], int, int, int]] = Counter()
    emitted = retained = rejected = 0
    raw_zero = retained_zero = 0
    zero_records = []

    for first_index, second_index, third_index in bucket_triples:
        first_values = buckets[keys[first_index]]
        second_values = buckets[keys[second_index]]
        third_values = buckets[keys[third_index]]
        expected = expected_multiset_count(
            first_index,
            second_index,
            third_index,
            len(first_values), len(second_values), len(third_values)
        )
        local_count = 0
        if first_index == second_index == third_index:
            iterator = combinations_with_replacement(first_values, 3)
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
            local_count += 1
            # The production stream declares coefficient-pair lexical order;
            # this is a canonical total order even though it differs from the
            # identity-real order used within squareclass buckets.
            edges = tuple(sorted(raw_edges))
            product = multiply(multiply(edges[0], edges[1]), edges[2])
            root = identity_positive_square_root(product)
            cycle = -root[0], -root[1]
            determinant = (
                8 - 2 * sum(edge[0] for edge in edges) + 2 * cycle[0],
                -2 * sum(edge[1] for edge in edges) + 2 * cycle[1],
            )
            selected_sign = sign(determinant)
            conjugate_sign = sign(determinant, 1)
            edge_kinds = tuple(
                sorted(
                    "finite" if sign((edge[0] - 4, edge[1])) < 0 else "dotted"
                    for edge in edges
                )
            )
            record = (
                ";".join(f"{a},{b}" for a, b in (*edges, cycle, determinant)) + "\n"
            ).encode()
            canonical.update(record)
            emitted += 1
            if selected_sign == 0:
                raw_zero += 1
                zero_records.append(
                    {
                        "edges": [list(edge) for edge in edges],
                        "cycle": list(cycle),
                        "determinant": list(determinant),
                        "conjugate_sign": conjugate_sign,
                    }
                )
                if conjugate_sign > 0:
                    retained_zero += 1
            if conjugate_sign > 0:
                retained += 1
                retained_digest.update(record)
            else:
                rejected += 1
            strata[(edge_kinds, sign(root, 1), selected_sign, conjugate_sign)] += 1
        if local_count != expected:
            raise AssertionError("bucket iterator multiplicity disagrees with closed formula")

    if emitted != formula_count or emitted != retained + rejected:
        raise AssertionError("independent stream partition failed")
    selected_positive = sum(
        count for (_kinds, _root, selected, conjugate), count in strata.items()
        if selected > 0 and conjugate > 0
    )
    selected_lorentzian = sum(
        count for (_kinds, _root, selected, conjugate), count in strata.items()
        if selected < 0 and conjugate > 0
    )
    return {
        "alphabet_count": len(alphabet),
        "nonzero_count": len(nonzero),
        "bucket_count": len(buckets),
        "maximum_bucket_size": max(map(len, buckets.values())),
        "compatible_bucket_triples": len(bucket_triples),
        "raw_unordered_nonzero_multisets": comb(len(nonzero) + 2, 3),
        "square_product_packets": emitted,
        "formula_square_product_packets": formula_count,
        "conjugate_positive_packets": retained,
        "conjugate_nonpositive_packets": rejected,
        "retained_selected_positive": selected_positive,
        "retained_selected_lorentzian": selected_lorentzian,
        "raw_selected_degenerate": raw_zero,
        "retained_selected_degenerate": retained_zero,
        "zero_records": zero_records,
        "strata": [
            {
                "edge_kinds": list(key[0]),
                "root_conjugate_sign": key[1],
                "selected_determinant_sign": key[2],
                "conjugate_determinant_sign": key[3],
                "count": count,
            }
            for key, count in sorted(strata.items())
        ],
        "canonical_sha256": canonical.hexdigest(),
        "retained_sha256": retained_digest.hexdigest(),
    }


def audit() -> dict[str, object]:
    replayed = replay()
    expected = {
        "alphabet_count": 6_875,
        "nonzero_count": 6_874,
        "bucket_count": 3_661,
        "maximum_bucket_size": 109,
        "compatible_bucket_triples": 27_088,
        "raw_unordered_nonzero_multisets": 54_158_527_500,
        "square_product_packets": 4_768_141,
        "formula_square_product_packets": 4_768_141,
        "conjugate_positive_packets": 2_633_885,
        "conjugate_nonpositive_packets": 2_134_256,
        "retained_selected_positive": 0,
        "retained_selected_lorentzian": 2_633_885,
        "raw_selected_degenerate": 1,
        "retained_selected_degenerate": 0,
        "canonical_sha256": "41386c1bfb64c497acc2ddb967c9a0d646f82e67dfc8907b8bc7d3b50e85bbe3",
        "retained_sha256": "8be12b4371f49e1a0dd35f2c290c4e5735e64da56dc209655d1e0af5b49476b6",
    }
    for key, value in expected.items():
        if replayed[key] != value:
            raise AssertionError(f"independent {key} drifted: {replayed[key]!r}")
    if replayed["zero_records"] != [
        {
            "edges": [[1, 0], [1, 0], [1, 0]],
            "cycle": [-1, 0],
            "determinant": [0, 0],
            "conjugate_sign": 0,
        }
    ]:
        raise AssertionError("the sole degenerate record is not finite (3,3,3)")

    production = json.loads(PRODUCTION.read_text())
    counts = production["counts"]
    production_comparison = {
        "edge_alphabet": replayed["alphabet_count"],
        "nonzero_edge_values": replayed["nonzero_count"],
        "element_squareclass_buckets": replayed["bucket_count"],
        "maximum_bucket_size": replayed["maximum_bucket_size"],
        "compatible_bucket_triples": replayed["compatible_bucket_triples"],
        "raw_unordered_nonzero_edge_multisets": replayed["raw_unordered_nonzero_multisets"],
        "integral_golden_square_product_packets": replayed["square_product_packets"],
        "conjugate_positive_definite_packets": replayed["conjugate_positive_packets"],
        "conjugate_nonpositive_rejections": replayed["conjugate_nonpositive_packets"],
        "retained_selected_positive_definite": replayed["retained_selected_positive"],
        "retained_selected_lorentzian": replayed["retained_selected_lorentzian"],
        "raw_selected_degenerate": replayed["raw_selected_degenerate"],
        "retained_selected_degenerate": replayed["retained_selected_degenerate"],
    }
    for key, value in production_comparison.items():
        if counts.get(key) != value:
            raise AssertionError(f"production count {key} is stale or incorrect")
    digests = production["stream_digests"]
    if (
        digests["all_integral_square_product_packets_sha256"]
        != replayed["canonical_sha256"]
        or digests["conjugate_positive_packets_sha256"]
        != replayed["retained_sha256"]
    ):
        raise AssertionError("production stream digest disagrees")

    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": (
            "Independent replay of the complete nonorthogonal three-node local "
            "packet tranche; no high-rank compatibility or emptiness claim."
        ),
        "independence": (
            "The declared length-two alphabet is imported, but the owner triangle "
            "module is not. Golden arithmetic, prime-ideal squareclasses, roots, "
            "bucket joins, multiplicities, determinants, strata, and digests are "
            "reimplemented."
        ),
        "production_artifact_sha256": sha256(PRODUCTION.read_bytes()).hexdigest(),
        "replay": replayed,
        "logic": {
            "squareclass_completeness": (
                "Z[phi] has class number one and units +/-phi^n. Ideal valuation "
                "parities plus conjugate sign therefore characterize squares of "
                "identity-positive elements; every bucket value was also checked "
                "against its representative by explicit square reconstruction."
            ),
            "multiplicity": (
                "Every compatible squareclass multiset i<=j<=k with XOR zero is "
                "generated exactly once. Closed multiset formulas are checked "
                "against the actual iterator count separately for all 27088 joins."
            ),
            "determinant": (
                "For diagonal 2 and directed cycle b=-sqrt(c12*c23*c31), exact "
                "expansion gives 8-2*sum(cij)+2*b. At the conjugate place all "
                "2x2 minors are positive, so determinant >0 is exactly the missing "
                "three-vector positive-definiteness condition."
            ),
            "degenerate_scope": (
                "Raw selected degeneracy and retained selected degeneracy are "
                "separate. The sole raw zero is the finite affine (3,3,3) record; "
                "it is also conjugate-degenerate and hence never retained."
            ),
            "compact_scope": (
                "This local census neither forces a complete triangle in a rank-23--30 "
                "frame nor checks overlapping cycles, gauges, orders, conductors, or "
                "Vinberg termination. Q(sqrt(5)) dimensions 22--29 remain open."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = audit()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "verdict": payload["verdict"],
                "square_product_packets": payload["replay"]["square_product_packets"],
                "retained_packets": payload["replay"]["conjugate_positive_packets"],
                "canonical_sha256": payload["replay"]["canonical_sha256"],
                "retained_sha256": payload["replay"]["retained_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
