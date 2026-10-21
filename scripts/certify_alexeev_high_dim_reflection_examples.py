#!/usr/bin/env python3
"""Exact audit of three high-dimensional matrices attached to Alexeev 2023.

Download the v4 arXiv ancillary file ``all_cases_in_2el_paper.txt`` and run

    python scripts/certify_alexeev_high_dim_reflection_examples.py PATH

The script parses only the three named Sage matrix assignments.  Rank and
inertia are then checked over the rationals.  No Sage installation is needed.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from sympy import Matrix


SOURCE_SHA256 = "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"

EXPECTED = {
    (15, 9, 1): {
        "matrix_size": 98,
        "rank": 15,
        "diagonal_counts": {-2: 84, -4: 14},
        "inertia": (1, 14, 83),
        "leading_block_determinant": 512,
    },
    (20, 2, 1): {
        "matrix_size": 50,
        "rank": 20,
        "diagonal_counts": {-2: 45, -4: 5},
        "inertia": (1, 19, 30),
        "leading_block_determinant": -4,
    },
    (22, 2, 0): {
        "matrix_size": 210,
        "rank": 22,
        "diagonal_counts": {-2: 42, -4: 168},
        "inertia": (1, 21, 188),
        "leading_block_determinant": -4,
    },
}


def _key_text(key: tuple[int, int, int]) -> str:
    return "(" + ",".join(str(value) for value in key) + ")"


def _extract_matrix(source: str, key: tuple[int, int, int]) -> Matrix:
    pattern = (
        r"^coxeter_matrix\["
        + re.escape(_key_text(key))
        + r"\] = matrix\((\[.*\])\)$"
    )
    match = re.search(pattern, source, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"matrix assignment {_key_text(key)} was not found")
    rows = ast.literal_eval(match.group(1))
    return Matrix(rows)


def _audit_matrix(matrix: Matrix, expected: dict[str, Any]) -> dict[str, Any]:
    size = expected["matrix_size"]
    rank_expected = expected["rank"]
    if matrix.shape != (size, size):
        raise ValueError(f"expected a {size} by {size} matrix, got {matrix.shape}")

    rank = int(matrix.rank())
    leading_determinants = [
        int(matrix[:index, :index].det())
        for index in range(1, rank_expected + 1)
    ]
    if any(value == 0 for value in leading_determinants):
        raise ValueError("the advertised leading principal block has a zero pivot")

    previous = 1
    pivot_signs: list[int] = []
    for determinant in leading_determinants:
        pivot_signs.append(1 if determinant * previous > 0 else -1)
        previous = determinant

    positive = pivot_signs.count(1)
    negative = pivot_signs.count(-1)
    inertia = (positive, negative, size - rank)
    diagonal_counts = dict(sorted(Counter(int(matrix[i, i]) for i in range(size)).items()))

    checks = {
        "symmetric": matrix == matrix.T,
        "rank": rank == rank_expected,
        "diagonal_counts": diagonal_counts == expected["diagonal_counts"],
        "inertia": inertia == expected["inertia"],
        "leading_block_determinant": (
            leading_determinants[-1] == expected["leading_block_determinant"]
        ),
    }
    return {
        "matrix_size": size,
        "rank": rank,
        "diagonal_counts": {str(key): value for key, value in diagonal_counts.items()},
        "leading_principal_determinants": leading_determinants,
        "ldl_pivot_signs": pivot_signs,
        "inertia_positive_negative_zero": list(inertia),
        "checks": checks,
        "verified": all(checks.values()),
    }


def certify(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    source = raw.decode("utf-8")
    examples = {
        _key_text(key): _audit_matrix(_extract_matrix(source, key), expected)
        for key, expected in EXPECTED.items()
    }
    return {
        "source": str(path),
        "source_sha256": digest,
        "expected_source_sha256": SOURCE_SHA256,
        "source_digest_matches": digest == SOURCE_SHA256,
        "examples": examples,
        "verified": digest == SOURCE_SHA256
        and all(example["verified"] for example in examples.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Alexeev's attached Sage matrix file")
    args = parser.parse_args()
    print(json.dumps(certify(args.source), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
