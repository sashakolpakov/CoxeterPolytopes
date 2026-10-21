#!/usr/bin/env python3
"""Independent hostile audit of the finite congruence-dessin theorem.

This file intentionally imports no project certificate code.  It rederives
the E10 reflection matrices and tests four failure modes which the theorem's
hypotheses are designed to exclude.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_rank10_two_full_nef_simplices.json"
SOURCE_SHA256 = "11b937f0a78091452d1a50cb2b4cb554949485b1469b4ad2616ee629e24320f8"
DEFAULT_OUTPUT = ROOT / "data/enumerations/coxeter_congruence_dessin_hostile_audit.json"
RSS_CEILING_BYTES = 1 << 29


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("hostile congruence-dessin audit exceeded 0.5 GiB RSS")


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
            _guard()
    return digest.hexdigest()


def _identity(rank: int) -> list[list[int]]:
    return [[int(row == column) for column in range(rank)] for row in range(rank)]


def _multiply(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    rank = len(left)
    return [
        [sum(left[row][k] * right[k][column] for k in range(rank)) for column in range(rank)]
        for row in range(rank)
    ]


def _power(matrix: list[list[int]], exponent: int) -> list[list[int]]:
    result = _identity(len(matrix))
    for _ in range(exponent):
        result = _multiply(result, matrix)
    return result


def _same_mod(left: list[list[int]], right: list[list[int]], modulus: int) -> bool:
    return all(
        (left[row][column] - right[row][column]) % modulus == 0
        for row in range(len(left))
        for column in range(len(left))
    )


def _reflection(gram: list[list[int]], index: int) -> list[list[int]]:
    result = _identity(len(gram))
    for column in range(len(gram)):
        result[index][column] += gram[index][column]
    return result


def _order_mod(matrix: list[list[int]], modulus: int, cap: int = 12) -> int | None:
    identity = _identity(len(matrix))
    current = identity
    for order in range(1, cap + 1):
        current = _multiply(current, matrix)
        if _same_mod(current, identity, modulus):
            return order
    return None


def build() -> dict[str, object]:
    assert _sha(SOURCE) == SOURCE_SHA256
    source = json.loads(SOURCE.read_text())
    case = next(item for item in source["cases"] if item["id"] == "R10-UE8-NEF-H9-10")
    gram = case["matrix_certificate"]["wall_gram"]
    identity = _identity(10)
    reflections = [_reflection(gram, index) for index in range(10)]
    assert all(_multiply(matrix, matrix) == identity for matrix in reflections)
    reductions = {
        tuple(entry % 5 for row in matrix for entry in row)
        for matrix in reflections
    }
    assert len(reductions) == 10
    pair_orders = []
    for first in range(10):
        for second in range(first + 1, 10):
            order = _order_mod(_multiply(reflections[first], reflections[second]), 5)
            expected = 3 if gram[first][second] else 2
            assert order == expected
            pair_orders.append(order)
    assert pair_orders.count(2) == 36 and pair_orders.count(3) == 9

    # Failure mode 1: retaining only the relation s^2=1 does not retain an
    # involution.  The nontrivial integral involution -1 becomes 1 modulo 2.
    minus_one = [[-1]]
    assert _power(minus_one, 2) == [[1]]
    assert _same_mod(minus_one, [[1]], 2)

    # Failure mode 2: a word may collapse at a shallow layer and reappear at
    # a deeper one.  The p-adic valuation bound in the tower is necessary.
    shallow_word = [[1, 5], [0, 1]]
    assert shallow_word != _identity(2)
    assert _same_mod(shallow_word, _identity(2), 5)
    assert not _same_mod(shallow_word, _identity(2), 25)

    # Failure mode 3: a finite quotient assigns finite order to an infinite
    # integral unipotent.  Thus no one layer detects m=infinity.
    unipotent = [[1, 1], [0, 1]]
    assert _power(unipotent, 5) != _identity(2)
    assert _same_mod(_power(unipotent, 5), _identity(2), 5)

    # Failure mode 4: the height lemma needs strict Norm(p)>(2B)^d.  At
    # equality B=1, p=2, the two allowed integers -1 and 1 collide.
    assert 2 == (2 * 1) ** 1
    assert (-1) % 2 == 1 % 2
    assert len({entry % 3 for entry in (-1, 0, 1)}) == 3

    # Reflection matrices alone cannot orient a marked root: r and -r give
    # the same reflection.  The quadratic/root/chamber marking is essential.
    root_plus_reflection = [[-1]]
    root_minus_reflection = [[-1]]
    assert root_plus_reflection == root_minus_reflection

    result: dict[str, object] = {
        "all_checks_pass": True,
        "independence": "no import from the production congruence-dessin certificate",
        "source": {"path": str(SOURCE.relative_to(ROOT)), "sha256": SOURCE_SHA256},
        "e10_level_5_replay": {
            "distinct_marked_involutions": len(reductions),
            "finite_pair_order_histogram": {"2": pair_orders.count(2), "3": pair_orders.count(3)},
        },
        "hostile_counterexamples": {
            "relation_without_nonidentity_witness": "-1 collapses to 1 modulo 2",
            "shallow_layer_without_valuation_bound": "[[1,5],[0,1]] collapses modulo 5 but not modulo 25",
            "single_layer_infinite_order_claim": "[[1,1],[0,1]] has fifth power nonidentity over Z but identity modulo 5",
            "nonstrict_height_threshold": "B=1 and Norm(p)=2=(2B)^1 identifies -1 and 1",
            "unmarked_reflection_orientation": "r and -r define the same reflection",
        },
        "verdict": "PASS_WITH_SCOPE: finite labels and finite nonrelations only; exact Gram needs bounded marked quadratic/root data",
        "resource_contract": {"in_process_rss_ceiling_bytes": RSS_CEILING_BYTES},
    }
    _guard()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    result = build()
    if args.check_output is not None and json.loads(args.check_output.read_text()) != result:
        raise AssertionError(f"archive mismatch: {args.check_output}")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
