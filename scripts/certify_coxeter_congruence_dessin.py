#!/usr/bin/env python3
"""Exact finite congruence-dessin calibration on the E10 H9 chamber.

The general number-field theorem is proved in the accompanying report.  This
certificate deliberately uses only integer arithmetic.  It verifies the
finite-order obstruction witnesses, a prescribed finite set of nonrelations,
and the height-to-one-level Gram reconstruction at the same (coincidentally
sufficient) rational prime 5.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import resource
import sys
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_rank10_two_full_nef_simplices.json"
SOURCE_SHA256 = "11b937f0a78091452d1a50cb2b4cb554949485b1469b4ad2616ee629e24320f8"
DEFAULT_OUTPUT = ROOT / "data/enumerations/coxeter_congruence_dessin_e10_h9.json"
RSS_CEILING_BYTES = 1 << 29

Matrix = tuple[int, ...]


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("congruence-dessin replay exceeded 0.5 GiB RSS")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
            _rss_guard()
    return digest.hexdigest()


def _identity(rank: int) -> Matrix:
    return tuple(1 if row == column else 0 for row in range(rank) for column in range(rank))


def _multiply(left: Matrix, right: Matrix, rank: int) -> Matrix:
    return tuple(
        sum(left[rank * row + k] * right[rank * k + column] for k in range(rank))
        for row in range(rank)
        for column in range(rank)
    )


def _power(matrix: Matrix, exponent: int, rank: int) -> Matrix:
    result = _identity(rank)
    base = matrix
    while exponent:
        if exponent & 1:
            result = _multiply(result, base, rank)
        exponent >>= 1
        if exponent:
            base = _multiply(base, base, rank)
    return result


def _transpose(matrix: Matrix, rank: int) -> Matrix:
    return tuple(matrix[rank * column + row] for row in range(rank) for column in range(rank))


def _difference(left: Matrix, right: Matrix) -> Matrix:
    return tuple(a - b for a, b in zip(left, right, strict=True))


def _content(matrix: Iterable[int]) -> int:
    result = 0
    for entry in matrix:
        result = math.gcd(result, abs(entry))
    if result == 0:
        raise AssertionError("a claimed nonrelation was the zero matrix")
    return result


def _prime_factors(value: int) -> set[int]:
    value = abs(value)
    factors: set[int] = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        factors.add(value)
    return factors


def _is_prime(value: int) -> bool:
    return value >= 2 and all(value % divisor for divisor in range(2, math.isqrt(value) + 1))


def _valuation(value: int, prime: int) -> int:
    if value == 0:
        raise AssertionError("valuation of zero is not used")
    answer = 0
    while value % prime == 0:
        answer += 1
        value //= prime
    return answer


def _determinant_bareiss(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows]
    rank = len(matrix)
    sign = 1
    denominator = 1
    for column in range(rank - 1):
        pivot = next((row for row in range(column, rank) if matrix[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign *= -1
        pivot_value = matrix[column][column]
        for row in range(column + 1, rank):
            for other in range(column + 1, rank):
                numerator = matrix[row][other] * pivot_value - matrix[row][column] * matrix[column][other]
                if numerator % denominator:
                    raise AssertionError("Bareiss division was not exact")
                matrix[row][other] = numerator // denominator
        denominator = pivot_value
        for row in range(column + 1, rank):
            matrix[row][column] = 0
    return sign * matrix[-1][-1]


def _reflection(gram: list[list[int]], index: int) -> Matrix:
    """Reflection in basis root e_index for a square-minus-two Gram matrix."""

    rank = len(gram)
    matrix = list(_identity(rank))
    # Column vectors: s_i(x)=x+(x,e_i)e_i, so only row i changes.
    for column in range(rank):
        matrix[rank * index + column] += gram[column][index]
    return tuple(matrix)


def _equal_mod(left: Matrix, right: Matrix, modulus: int) -> bool:
    return all((a - b) % modulus == 0 for a, b in zip(left, right, strict=True))


def _order_mod(matrix: Matrix, modulus: int, rank: int, limit: int) -> int | None:
    identity = _identity(rank)
    current = identity
    for exponent in range(1, limit + 1):
        current = _multiply(current, matrix, rank)
        if _equal_mod(current, identity, modulus):
            return exponent
    return None


def _matrix_digest(matrix: Matrix) -> str:
    payload = ",".join(str(entry) for entry in matrix).encode()
    return hashlib.sha256(payload).hexdigest()


def build() -> dict[str, object]:
    actual_sha = _sha256(SOURCE)
    if actual_sha != SOURCE_SHA256:
        raise AssertionError(f"source SHA changed: {actual_sha} != {SOURCE_SHA256}")
    source = json.loads(SOURCE.read_text())
    case = next(item for item in source["cases"] if item["id"] == "R10-UE8-NEF-H9-10")
    gram = case["matrix_certificate"]["wall_gram"]
    rank = len(gram)
    assert rank == 10 and all(len(row) == rank for row in gram)
    assert gram == [list(row) for row in zip(*gram, strict=True)]
    assert [gram[index][index] for index in range(rank)] == [-2] * rank
    determinant = _determinant_bareiss(gram)
    assert determinant == -1

    identity = _identity(rank)
    gram_matrix = tuple(entry for row in gram for entry in row)
    reflections = tuple(_reflection(gram, index) for index in range(rank))
    assert all(_multiply(matrix, matrix, rank) == identity for matrix in reflections)
    assert all(matrix != identity for matrix in reflections)
    assert len(set(reflections)) == rank
    assert all(
        _multiply(_multiply(_transpose(matrix, rank), gram_matrix, rank), matrix, rank)
        == gram_matrix
        for matrix in reflections
    )

    pair_rows: list[dict[str, int]] = []
    obstruction_contents: list[dict[str, object]] = []
    for index, matrix in enumerate(reflections):
        obstruction_contents.append(
            {"kind": "generator_nonidentity", "indices": [index], "content": _content(_difference(matrix, identity))}
        )
    for first in range(rank):
        for second in range(first + 1, rank):
            obstruction_contents.append(
                {
                    "kind": "marked_generators_distinct",
                    "indices": [first, second],
                    "content": _content(_difference(reflections[first], reflections[second])),
                }
            )
            product = _multiply(reflections[first], reflections[second], rank)
            expected_order = 3 if gram[first][second] == 1 else 2
            assert gram[first][second] in (0, 1)
            assert _power(product, expected_order, rank) == identity
            proper_power_contents: dict[str, int] = {}
            for prime_divisor in sorted(_prime_factors(expected_order)):
                exponent = expected_order // prime_divisor
                witness = _difference(_power(product, exponent, rank), identity)
                witness_content = _content(witness)
                proper_power_contents[str(exponent)] = witness_content
                obstruction_contents.append(
                    {
                        "kind": "pair_proper_power_nonidentity",
                        "indices": [first, second],
                        "expected_order": expected_order,
                        "exponent": exponent,
                        "content": witness_content,
                    }
                )
            pair_rows.append(
                {
                    "first": first,
                    "second": second,
                    "gram_entry": gram[first][second],
                    "exact_order": expected_order,
                    "proper_power_witness_content": min(proper_power_contents.values()),
                }
            )

    # Exact-order preservation itself does not require tameness.  We impose
    # p not dividing 2*lcm(m_ij) for the clean algebraic deck interpretation.
    bad_primes = _prime_factors(2 * math.lcm(*(row["exact_order"] for row in pair_rows)))
    for obstruction in obstruction_contents:
        bad_primes.update(_prime_factors(int(obstruction["content"])))

    degree = 1
    height_bound = max(abs(entry) for row in gram for entry in row)
    reconstruction_threshold = (2 * height_bound) ** degree
    prime = next(
        candidate
        for candidate in range(2, 100)
        if _is_prime(candidate) and candidate not in bad_primes and candidate > reconstruction_threshold
    )
    assert prime == 5

    reduced_reflections = {tuple(entry % prime for entry in matrix) for matrix in reflections}
    assert len(reduced_reflections) == rank
    assert all(not _equal_mod(matrix, identity, prime) for matrix in reflections)
    for row in pair_rows:
        product = _multiply(reflections[row["first"]], reflections[row["second"]], rank)
        assert _order_mod(product, prime, rank, row["exact_order"]) == row["exact_order"]

    # A finite list of genuine matrix nonrelations in the fixed-p tower.
    coxeter_element = identity
    for reflection in reflections:
        coxeter_element = _multiply(coxeter_element, reflection, rank)
    nonrelation_rows: list[dict[str, object]] = []
    for exponent in range(1, 25):
        power = _power(coxeter_element, exponent, rank)
        difference = _difference(power, identity)
        content = _content(difference)
        nonrelation_rows.append(
            {
                "word": f"(s0*s1*...*s9)^{exponent}",
                "power": exponent,
                "difference_content": content,
                "p_adic_valuation": _valuation(content, prime),
                "matrix_sha256": _matrix_digest(power),
            }
        )
    tower_exponent = 1 + max(int(row["p_adic_valuation"]) for row in nonrelation_rows)
    modulus = prime**tower_exponent
    assert all(
        not _equal_mod(_power(coxeter_element, int(row["power"]), rank), identity, modulus)
        for row in nonrelation_rows
    )

    # Height-to-one-level lifting is logically separate from order
    # preservation.  Here k=Q, B=2, and 5>(2B)^1=4.
    alphabet = list(range(-height_bound, height_bound + 1))
    residue_to_integer = {entry % prime: entry for entry in alphabet}
    assert len(residue_to_integer) == len(alphabet)
    reconstructed_gram = [
        [residue_to_integer[entry % prime] for entry in row]
        for row in gram
    ]
    assert reconstructed_gram == gram

    reduced_payload = json.dumps(
        [[entry % prime for entry in matrix] for matrix in reflections],
        separators=(",", ":"),
    ).encode()
    certificate: dict[str, object] = {
        "all_checks_pass": True,
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": SOURCE_SHA256,
            "case_id": case["id"],
        },
        "calibration": {
            "id": "E10-H9-10",
            "hyperbolic_dimension": case["hyperbolic_dimension"],
            "facet_count": case["facet_count"],
            "lattice": case["lattice"],
            "gram_rank": rank,
            "gram_determinant": determinant,
            "finite_pair_order_histogram": {
                "2": sum(row["exact_order"] == 2 for row in pair_rows),
                "3": sum(row["exact_order"] == 3 for row in pair_rows),
            },
        },
        "finite_order_level": {
            "principal_congruence_ideal": f"({prime})",
            "rational_prime": prime,
            "tame_excluded_primes": sorted(_prime_factors(2 * math.lcm(*(row["exact_order"] for row in pair_rows)))),
            "all_obstruction_contents": sorted({int(row["content"]) for row in obstruction_contents}),
            "marked_involution_count": len(reduced_reflections),
            "exact_pair_orders_mod_level": True,
            "order_2_pair_count": sum(row["exact_order"] == 2 for row in pair_rows),
            "order_3_pairs": [
                [row["first"], row["second"]]
                for row in pair_rows
                if row["exact_order"] == 3
            ],
            "pair_data_sha256": hashlib.sha256(
                json.dumps(pair_rows, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
            "reduced_marked_reflections_sha256": hashlib.sha256(reduced_payload).hexdigest(),
        },
        "finite_nonrelation_tower": {
            "fixed_prime_ideal": f"({prime})",
            "tested_level": f"({prime}^{tower_exponent})",
            "tower_exponent": tower_exponent,
            "prescribed_nonrelation_count": len(nonrelation_rows),
            "all_prescribed_nonrelations_separated": True,
            "prescribed_powers": [int(row["power"]) for row in nonrelation_rows],
            "difference_content_histogram": {
                str(content): sum(row["difference_content"] == content for row in nonrelation_rows)
                for content in sorted({int(row["difference_content"]) for row in nonrelation_rows})
            },
            "nonrelation_data_sha256": hashlib.sha256(
                json.dumps(nonrelation_rows, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
        },
        "height_one_level_reconstruction": {
            "logical_role": "separate from finite-order preservation",
            "field": "Q",
            "field_degree": degree,
            "conjugate_height_bound_B": height_bound,
            "strict_norm_threshold": reconstruction_threshold,
            "prime_ideal_norm": prime,
            "strict_inequality": prime > reconstruction_threshold,
            "bounded_integral_alphabet": alphabet,
            "residues_are_injective": True,
            "exact_gram_reconstructed": True,
        },
        "scope": {
            "determinant_line_extension": "block-diagonal extension by the identity on the auxiliary line preserves this certificate",
            "single_finite_level_infinite_relations": "not_detected",
            "stack_vs_neat_cover": "finite deck monodromy retains marked involutions; a neat cover has no inertia",
        },
        "resource_contract": {
            "in_process_rss_ceiling_bytes": RSS_CEILING_BYTES,
            "outer_watchdog_command": "python3 scripts/run_with_rss_limit.py --memory-gib 0.5 -- python3 scripts/certify_coxeter_congruence_dessin.py --check-output data/enumerations/coxeter_congruence_dessin_e10_h9.json",
        },
    }
    _rss_guard()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    certificate = build()
    if args.check_output is not None:
        archived = json.loads(args.check_output.read_text())
        if archived != certificate:
            raise AssertionError(f"archive mismatch: {args.check_output}")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
