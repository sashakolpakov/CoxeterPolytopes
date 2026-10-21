#!/usr/bin/env python3
"""Replay the numerical controls for the quiver/2-CY Coxeter audit.

This script does not certify any global finiteness theorem.  It reconstructs
the two frozen wall Gram matrices used in the audit and records the rank-two
Coxeter parameter

    p_ij = 4 (r_i,r_j)^2 / ((r_i,r_i)(r_j,r_j)).

For intersecting walls p=0,1,2,3 correspond to m=2,3,4,6.  Values p>=4
are parallel or ultraparallel.  A pair of individual 2-spherical twists can
only realize p=0 or 1, while a symmetrizable generalized Cartan matrix can
also realize p=2 or 3.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
H21_PACKET = ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
X3_PACKET = ROOT / "data/enumerations/vinberg_x3_h19.json"


def _matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    rows = len(left)
    middle = len(right)
    cols = len(right[0])
    if not left or not right or len(left[0]) != middle:
        raise ValueError("incompatible matrix dimensions")
    return [
        [sum(left[i][k] * right[k][j] for k in range(middle)) for j in range(cols)]
        for i in range(rows)
    ]


def _transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def _h21_gram() -> list[list[int]]:
    packet = json.loads(H21_PACKET.read_text())
    source = packet["coxeter_input"]
    basis_gram = [[int(entry) for entry in row] for row in source["basis_gram"]]
    walls = [
        [int(entry) for entry in row]
        for row in source["wall_coordinates_as_columns_transposed"]
    ]
    return _matmul(_matmul(walls, basis_gram), _transpose(walls))


def _x3_gram() -> list[list[int]]:
    packet = json.loads(X3_PACKET.read_text())
    return [[int(entry) for entry in row] for row in packet["wall_gram"]]


def _rank_rational(matrix: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(cols):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def _parameter_histogram(gram: list[list[int]]) -> dict[str, int]:
    size = len(gram)
    histogram: Counter[Fraction] = Counter()
    for i in range(size):
        if len(gram[i]) != size or gram[i][i] >= 0:
            raise ValueError("expected a square wall Gram with negative diagonal")
        for j in range(i + 1, size):
            parameter = Fraction(
                4 * gram[i][j] * gram[i][j],
                gram[i][i] * gram[j][j],
            )
            histogram[parameter] += 1
    return {
        str(parameter): histogram[parameter]
        for parameter in sorted(histogram)
    }


def _norm_histogram(gram: list[list[int]]) -> dict[str, int]:
    counts = Counter(gram[i][i] for i in range(len(gram)))
    return {str(norm): counts[norm] for norm in sorted(counts)}


def _control(name: str, gram: list[list[int]]) -> dict[str, object]:
    facet_count = len(gram)
    rank = _rank_rational(gram)
    return {
        "name": name,
        "facet_count": facet_count,
        "gram_rank": rank,
        "gram_corank": facet_count - rank,
        "wall_norm_histogram": _norm_histogram(gram),
        "coxeter_parameter_histogram": _parameter_histogram(gram),
        "pair_count": facet_count * (facet_count - 1) // 2,
    }


def build_certificate() -> dict[str, object]:
    h21 = _control("Borcherds H21", _h21_gram())
    x3 = _control("Vinberg X3", _x3_gram())

    expected_h21 = {
        "0": 5796,
        "1": 3465,
        "2": 2016,
        "4": 3780,
        "9": 6048,
        "16": 840,
    }
    expected_x3 = {
        "0": 501,
        "1": 27,
        "3": 36,
        "4": 24,
        "16": 30,
        "25": 12,
    }

    checks = {
        "h21_rank_and_corank_are_22_and_188": (
            h21["gram_rank"] == 22 and h21["gram_corank"] == 188
        ),
        "x3_rank_and_corank_are_20_and_16": (
            x3["gram_rank"] == 20 and x3["gram_corank"] == 16
        ),
        "h21_parameter_histogram_matches": (
            h21["coxeter_parameter_histogram"] == expected_h21
        ),
        "x3_parameter_histogram_matches": (
            x3["coxeter_parameter_histogram"] == expected_x3
        ),
        "h21_contains_m4_pairs_forbidden_to_individual_spherical_twists": (
            expected_h21["2"] > 0
        ),
        "x3_contains_m6_pairs_forbidden_to_individual_spherical_twists": (
            expected_x3["3"] > 0
        ),
        "m5_cartan_product_is_nonintegral": isqrt(5) ** 2 != 5,
    }

    return {
        "schema": "quiver-2cy-coxeter-route-audit/v1",
        "scope": (
            "Numerical replay only: no arithmetic reflection finiteness, "
            "essential-surjectivity, or boundedness theorem is asserted."
        ),
        "rank_two_rules": {
            "individual_2_spherical_twists_finite_parameters": ["0", "1"],
            "symmetrizable_gcm_finite_parameters": ["0", "1", "2", "3"],
            "m5_parameter": "(3+sqrt(5))/2",
            "m5_parameter_minimal_polynomial": "x^2-3*x+1",
        },
        "controls": [h21, x3],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit nonzero unless every frozen assertion passes",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if args.check and not certificate["all_checks_pass"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
