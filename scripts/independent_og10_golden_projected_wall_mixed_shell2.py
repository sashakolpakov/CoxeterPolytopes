#!/usr/bin/env python3
"""Independent exact replay of the mixed-E8 ``b^2 <= 4`` wall box.

This file deliberately imports no project code.  In particular, it does not
use the translated-sphere Schur recursion in
``search_og10_golden_projected_wall.py``.  Instead it uses the necessary bound

    a^2 <= 2 (w')^2 + 2 (phi')^2 b^2
        < 4 + 8 (phi')^2 < 8

for ``w' = a + phi' b``, ``(w')^2 < 2``, and ``b^2 <= 4``.  Since E8 is even,
every possible ``a`` is therefore in the exact shells of square 0, 2, 4, or
6.  Those shells and the square-2/square-4 ``b`` shells are generated as
Weyl orbits and pinned against the E8 theta-series coefficients.  Exhausting
their Cartesian product is an independent completeness proof for this box.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE_PATH = (
    ROOT
    / "data/candidates/og10_golden_projected_wall_mixed_shell2_no_go.json"
)

# Bourbaki-style simple-root Gram matrix in the exact coordinate convention
# used by the archived OG10 fixed slice.
E8_GRAM: tuple[tuple[int, ...], ...] = (
    (2, -1, 0, 0, 0, 0, 0, 0),
    (-1, 2, -1, 0, 0, 0, 0, 0),
    (0, -1, 2, -1, 0, 0, 0, -1),
    (0, 0, -1, 2, -1, 0, 0, 0),
    (0, 0, 0, -1, 2, -1, 0, 0),
    (0, 0, 0, 0, -1, 2, -1, 0),
    (0, 0, 0, 0, 0, -1, 2, 0),
    (0, 0, -1, 0, 0, 0, 0, 2),
)

# One seed in each nonzero shell needed below.  Simple-root reflections close
# each seed into a full Weyl orbit; the resulting cardinalities equal the
# exact theta-series coefficients 240 sigma_3(m).
SHELL_SEEDS: dict[int, tuple[int, ...]] = {
    2: (-1, -1, -1, -1, -1, -1, -1, -1),
    4: (-1, -1, -1, -1, -1, 0, -1, -1),
    6: (-1, -1, -1, -1, -1, -1, -1, 1),
}
EXPECTED_SHELL_COUNTS = {2: 240, 4: 2160, 6: 6720}

# These are respectively 400 times the archived invariant interior and 13
# times its v61 escape ray, restricted to the E8 coordinates.
INTERIOR_E8_SCALED = (184, 364, 540, 440, 336, 228, 116, 272)
RAY_E8_PHI_COEFFICIENTS_SCALED = (26, 52, 78, 64, 48, 32, 16, 40)
INTERIOR_XU_SCALED = ((601, 199), (2, 0))
RAY_XU_SCALED = ((13, 0), (0, 1))

Pair = tuple[int, int]
Vector = tuple[int, ...]
MixedRecord = tuple[
    Pair,
    tuple[Pair, ...],
    Pair,
    Pair,
    Vector,
    Vector,
]
CoefficientRecord = tuple[Pair, Pair, Pair, Pair, Pair]

ZERO_PAIR: Pair = (0, 0)
INV_PHI_PAIR: Pair = (-1, 1)


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def neg(value: Pair) -> Pair:
    return -value[0], -value[1]


def mul(left: Pair, right: Pair) -> Pair:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def scale(value: Pair, scalar: int) -> Pair:
    return scalar * value[0], scalar * value[1]


def conjugate(value: Pair) -> Pair:
    return value[0] + value[1], -value[1]


def sign(value: Pair) -> int:
    """Exact sign of ``a+b*phi`` using integer comparisons only."""

    a, b = value
    rational_part = 2 * a + b
    if rational_part == 0:
        return (b > 0) - (b < 0)
    if b == 0 or (rational_part > 0) == (b > 0):
        return (rational_part > 0) - (rational_part < 0)
    comparison = rational_part * rational_part - 5 * b * b
    if comparison > 0:
        return (rational_part > 0) - (rational_part < 0)
    return (b > 0) - (b < 0)


def vector_square(vector: Vector) -> int:
    return sum(
        vector[i] * E8_GRAM[i][j] * vector[j]
        for i in range(8)
        for j in range(8)
    )


def gram_vector(vector: Vector) -> Vector:
    return tuple(
        sum(E8_GRAM[i][j] * vector[j] for j in range(8)) for i in range(8)
    )


def inner(left: Vector, gram_right: Vector) -> int:
    return sum(left[i] * gram_right[i] for i in range(8))


@lru_cache(maxsize=None)
def weyl_shell(square: int) -> tuple[Vector, ...]:
    """Generate one E8 shell by simple-root Weyl reflections."""

    seed = SHELL_SEEDS[square]
    if vector_square(seed) != square:
        raise AssertionError(f"bad E8 shell seed for square {square}")
    orbit = {seed}
    frontier = [seed]
    while frontier:
        vector = frontier.pop()
        covector = gram_vector(vector)
        for index in range(8):
            reflected = list(vector)
            reflected[index] -= covector[index]
            candidate = tuple(reflected)
            if candidate not in orbit:
                if vector_square(candidate) != square:
                    raise AssertionError("simple-root reflection changed the norm")
                orbit.add(candidate)
                frontier.append(candidate)
    if len(orbit) != EXPECTED_SHELL_COUNTS[square]:
        raise AssertionError(
            f"E8 shell {square} has {len(orbit)} vectors, expected "
            f"{EXPECTED_SHELL_COUNTS[square]}"
        )
    return tuple(sorted(orbit))


def digest_record(hasher: object, record: object) -> None:
    hasher.update(repr(record).encode("ascii"))
    hasher.update(b"\n")


def is_signature_admissible(projected_square_twice: Pair) -> bool:
    if sign(projected_square_twice) <= 0:
        return False
    schur_twice = add((4, 0), neg(projected_square_twice))
    return sign(schur_twice) > 0 and sign(conjugate(schur_twice)) > 0


def coefficient_records(bound: int) -> tuple[CoefficientRecord, ...]:
    box = tuple(
        (a, b)
        for a in range(-bound, bound + 1)
        for b in range(-bound, bound + 1)
    )
    x_data: list[tuple[Pair, Pair]] = []
    for x_covector in box:
        term = neg(mul(mul(x_covector, x_covector), INV_PHI_PAIR))
        if sign(conjugate(add((4, 0), neg(term)))) > 0:
            x_data.append((x_covector, term))
    u_data: list[tuple[Pair, Pair]] = []
    for u_covector in box:
        term = mul(u_covector, u_covector)
        if sign(conjugate(add((4, 0), neg(term)))) > 0:
            u_data.append((u_covector, term))
    records = []
    for x_covector, x_term in x_data:
        for u_covector, u_term in u_data:
            records.append(
                (
                    x_covector,
                    u_covector,
                    add(x_term, u_term),
                    add(
                        mul(x_covector, INTERIOR_XU_SCALED[0]),
                        mul(u_covector, INTERIOR_XU_SCALED[1]),
                    ),
                    add(
                        mul(x_covector, RAY_XU_SCALED[0]),
                        mul(u_covector, RAY_XU_SCALED[1]),
                    ),
                )
            )
    return tuple(records)


@lru_cache(maxsize=1)
def independent_mixed_records() -> tuple[MixedRecord, ...]:
    """Enumerate every eligible ``(a,b)`` without translated Schur pruning."""

    # Exact check of 4 + 8(phi')^2 < 8.  Since (phi')^2=2-phi,
    # the positive difference is -12+8phi.
    if sign((-12, 8)) <= 0:
        raise AssertionError("the independent a^2 < 8 bound failed")

    zero = (0,) * 8
    a_vectors = (zero,) + tuple(
        vector
        for square in (2, 4, 6)
        for vector in weyl_shell(square)
    )
    b_vectors = tuple(
        vector for square in (2, 4) for vector in weyl_shell(square)
    )
    if len(a_vectors) != 9121 or len(b_vectors) != 2400:
        raise AssertionError("unexpected independent E8 box cardinality")

    a_data = tuple(
        (vector, vector_square(vector), gram_vector(vector)) for vector in a_vectors
    )
    records: list[MixedRecord] = []
    for b_vector in b_vectors:
        b_square = vector_square(b_vector)
        gram_b = gram_vector(b_vector)
        b_on_interior = sum(
            gram_b[i] * INTERIOR_E8_SCALED[i] for i in range(8)
        )
        for a_vector, a_square, gram_a in a_data:
            product = inner(a_vector, gram_b)
            # w^2=(a^2+b^2)+(b^2+2(a,b))*phi.
            square = (a_square + b_square, b_square + 2 * product)
            conjugate_square = conjugate(square)
            if sign(conjugate_square) <= 0:
                continue
            if sign(add((2, 0), neg(conjugate_square))) <= 0:
                continue
            if sign(square) <= 0:
                raise AssertionError("positive-definite identity-place square failed")

            covector = tuple((gram_a[i], gram_b[i]) for i in range(8))
            a_on_interior = sum(
                gram_a[i] * INTERIOR_E8_SCALED[i] for i in range(8)
            )
            on_interior = (a_on_interior, b_on_interior)
            # (x+y*phi)*(r*phi) = yr + (x+y)r*phi.
            on_ray = (
                sum(
                    gram_b[i] * RAY_E8_PHI_COEFFICIENTS_SCALED[i]
                    for i in range(8)
                ),
                sum(
                    (gram_a[i] + gram_b[i])
                    * RAY_E8_PHI_COEFFICIENTS_SCALED[i]
                    for i in range(8)
                ),
            )
            records.append(
                (
                    scale(square, 2),
                    covector,
                    on_interior,
                    on_ray,
                    a_vector,
                    b_vector,
                )
            )
    return tuple(sorted(records, key=lambda record: (record[5], record[4])))


@lru_cache(maxsize=None)
def run_independent_replay(coefficient_bound: int = 30) -> dict[str, object]:
    records = independent_mixed_records()
    eligible_hasher = hashlib.sha256()
    by_square: dict[Pair, list[MixedRecord]] = {}
    a_square_histogram: Counter[int] = Counter()
    b_square_histogram: Counter[int] = Counter()
    for record in records:
        by_square.setdefault(record[0], []).append(record)
        a_square_histogram[vector_square(record[4])] += 1
        b_square_histogram[vector_square(record[5])] += 1
        digest_record(
            eligible_hasher,
            (record[4], record[5], record[0], record[1], record[2], record[3]),
        )

    coefficients = coefficient_records(coefficient_bound)
    admissible_by_square: dict[Pair, tuple[CoefficientRecord, ...]] = {}
    signature_parameter_count = 0
    for e8_square_twice in sorted(by_square):
        admissible = tuple(
            data
            for data in coefficients
            if is_signature_admissible(add(data[2], e8_square_twice))
        )
        admissible_by_square[e8_square_twice] = admissible
        signature_parameter_count += len(admissible)

    wall_hasher = hashlib.sha256()
    wall_test_count = 0
    cutting_records = []
    separation_sign_histogram: Counter[str] = Counter()
    for (
        e8_square_twice,
        _,
        e8_on_interior,
        e8_on_ray,
        a_vector,
        b_vector,
    ) in records:
        for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
            admissible_by_square[e8_square_twice]
        ):
            projected_square_twice = add(non_e8_square, e8_square_twice)
            if not is_signature_admissible(projected_square_twice):
                raise AssertionError("an enumerated wall failed the exact signature gate")
            on_interior = add(non_e8_interior, e8_on_interior)
            on_ray = add(non_e8_ray, e8_on_ray)
            signs = (sign(on_interior), sign(on_ray))
            separation_sign_histogram[f"{signs[0]},{signs[1]}"] += 1
            wall_record = (
                a_vector,
                b_vector,
                A,
                C,
                projected_square_twice,
                on_interior,
                on_ray,
            )
            digest_record(wall_hasher, wall_record)
            wall_test_count += 1
            if signs[0] * signs[1] == -1:
                cutting_records.append(wall_record)

    return {
        "replay": "standalone_E8_Weyl_shell_Cartesian_product",
        "imports_project_code": False,
        "coefficient_bound": coefficient_bound,
        "necessity_bound": "a^2<4+8*(phi')^2<8; E8 even, hence a^2 in {0,2,4,6}",
        "a_vector_count": 9121,
        "b_vector_count": 2400,
        "raw_a_b_pair_count": 9121 * 2400,
        "theta_pinned_shell_counts": {
            str(square): EXPECTED_SHELL_COUNTS[square] for square in (2, 4, 6)
        },
        "eligible_w_count": len(records),
        "eligible_w_square_class_count": len(by_square),
        "eligible_w_a_square_histogram": {
            str(key): a_square_histogram[key] for key in sorted(a_square_histogram)
        },
        "eligible_w_b_square_histogram": {
            str(key): b_square_histogram[key] for key in sorted(b_square_histogram)
        },
        "eligible_w_records_sha256": eligible_hasher.hexdigest(),
        "signature_admissible_parameter_pair_count_summed_over_w_square_classes": signature_parameter_count,
        "exact_wall_test_count": wall_test_count,
        "separation_sign_histogram": {
            key: separation_sign_histogram[key]
            for key in sorted(separation_sign_histogram)
        },
        "canonical_wall_records_sha256": wall_hasher.hexdigest(),
        "cutting_covector_count": len(cutting_records),
        "cutting_records": cutting_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coefficient-bound", type=int, default=30)
    parser.add_argument(
        "--compare-certificate",
        type=Path,
        default=None,
        help="compare the independent replay block to a frozen shell2 certificate",
    )
    args = parser.parse_args()
    replay = run_independent_replay(args.coefficient_bound)
    if args.compare_certificate is not None:
        payload = json.loads(args.compare_certificate.read_text())
        expected = payload["independent_replay"]
        if replay != expected:
            raise AssertionError("independent shell2 replay differs from certificate")
        print(json.dumps({"independent_replay_matches": True}, sort_keys=True))
        return
    print(json.dumps(replay, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
