#!/usr/bin/env python3
"""Standalone exact replay of the mixed-E8 ``b^2 <= 6`` wall box.

This file imports no project code and does not use production's translated
Schur recursion.  If ``w'=a+phi'*b``, ``(w')^2<2``, and ``b^2<=6``, then

    a^2 <= 2(w')^2 + 2(phi')^2 b^2
        < 4 + 12(phi')^2 = 28 - 12phi < 10.

Evenness of E8 therefore restricts ``a^2`` to 0,2,4,6,8.  Exact Weyl orbits,
pinned to the E8 theta series, generate these shells.  Each of the b-shells
2,4,6 is one Weyl orbit.  We test the complete a-shell list at one b seed and
transport the eligible list to every b by an explicit integral Weyl matrix.
This is a bijective symmetry quotient of all 242,965,920 raw (a,b) pairs.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterator


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE_PATH = (
    ROOT
    / "data/candidates/og10_golden_projected_wall_mixed_shell3_no_go.json"
)

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

# Square 8 has two Weyl orbits: a primitive orbit of size 17,280 and the
# 240 doubled roots.  The lower three nonzero shells are single orbits.
SHELL_ORBIT_SEEDS: dict[int, tuple[tuple[int, ...], ...]] = {
    2: ((-1, -1, -1, -1, -1, -1, -1, -1),),
    4: ((-1, -1, -1, -1, -1, 0, -1, -1),),
    6: ((-1, -1, -1, -1, -1, -1, -1, 1),),
    8: (
        (-1, 0, -1, 0, 0, 0, 0, 1),
        (-2, 0, 0, 0, 0, 0, 0, 0),
    ),
}
EXPECTED_ORBIT_COUNTS = {
    2: (240,),
    4: (2160,),
    6: (6720,),
    8: (17280, 240),
}
EXPECTED_SHELL_COUNTS = {2: 240, 4: 2160, 6: 6720, 8: 17520}
EXPECTED_SEED_ELIGIBLE_COUNTS = {2: 58, 4: 80, 6: 57}

INTERIOR_E8_SCALED = (184, 364, 540, 440, 336, 228, 116, 272)
RAY_E8_PHI_COEFFICIENTS_SCALED = (26, 52, 78, 64, 48, 32, 16, 40)
INTERIOR_XU_SCALED = ((601, 199), (2, 0))
RAY_XU_SCALED = ((13, 0), (0, 1))

Pair = tuple[int, int]
Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
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


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(8)) for i in range(8)
    )


def reflected_matrix(matrix: Matrix, simple_index: int) -> Matrix:
    rows = [list(row) for row in matrix]
    rows[simple_index] = [
        matrix[simple_index][j]
        - sum(E8_GRAM[simple_index][k] * matrix[k][j] for k in range(8))
        for j in range(8)
    ]
    return tuple(tuple(row) for row in rows)


def reflect_vector(vector: Vector, simple_index: int) -> Vector:
    covector = gram_vector(vector)
    reflected = list(vector)
    reflected[simple_index] -= covector[simple_index]
    return tuple(reflected)


@lru_cache(maxsize=None)
def weyl_orbit(seed: Vector) -> tuple[Vector, ...]:
    square = vector_square(seed)
    orbit = {seed}
    frontier = [seed]
    while frontier:
        vector = frontier.pop()
        for index in range(8):
            candidate = reflect_vector(vector, index)
            if candidate not in orbit:
                if vector_square(candidate) != square:
                    raise AssertionError("a simple-root reflection changed the norm")
                orbit.add(candidate)
                frontier.append(candidate)
    return tuple(sorted(orbit))


@lru_cache(maxsize=None)
def weyl_shell(square: int) -> tuple[Vector, ...]:
    shell: set[Vector] = set()
    observed_orbit_counts = []
    for seed in SHELL_ORBIT_SEEDS[square]:
        if vector_square(seed) != square:
            raise AssertionError(f"bad E8 shell seed for square {square}")
        orbit = set(weyl_orbit(seed))
        if shell & orbit:
            raise AssertionError(f"declared E8 square-{square} orbits overlap")
        shell.update(orbit)
        observed_orbit_counts.append(len(orbit))
    if tuple(observed_orbit_counts) != EXPECTED_ORBIT_COUNTS[square]:
        raise AssertionError(
            f"E8 square-{square} orbit counts changed: {observed_orbit_counts}"
        )
    if len(shell) != EXPECTED_SHELL_COUNTS[square]:
        raise AssertionError(
            f"E8 shell {square} has {len(shell)} vectors, expected "
            f"{EXPECTED_SHELL_COUNTS[square]} from 240*sigma_3({square // 2})"
        )
    return tuple(sorted(shell))


@lru_cache(maxsize=None)
def transport_orbit(square: int) -> dict[Vector, Matrix]:
    """Map the shell seed to every vector by an explicit Weyl matrix."""

    seeds = SHELL_ORBIT_SEEDS[square]
    if len(seeds) != 1:
        raise AssertionError("transport_orbit is used only for one-orbit b-shells")
    seed = seeds[0]
    identity: Matrix = tuple(
        tuple(int(i == j) for j in range(8)) for i in range(8)
    )
    transports = {seed: identity}
    frontier = [seed]
    while frontier:
        vector = frontier.pop()
        matrix = transports[vector]
        for index in range(8):
            candidate = reflect_vector(vector, index)
            if candidate in transports:
                continue
            candidate_matrix = reflected_matrix(matrix, index)
            if matrix_vector(candidate_matrix, seed) != candidate:
                raise AssertionError("Weyl transport matrix has the wrong image")
            transports[candidate] = candidate_matrix
            frontier.append(candidate)
    if set(transports) != set(weyl_shell(square)):
        raise AssertionError(f"transport orbit does not exhaust square {square}")
    return transports


@lru_cache(maxsize=1)
def a_vectors() -> tuple[Vector, ...]:
    # 10-(4+12(phi')^2)=12phi-18 is strictly positive.
    if sign((-18, 12)) <= 0:
        raise AssertionError("the independent a^2<10 bound failed")
    vectors = ((0,) * 8,) + tuple(
        vector
        for square in (2, 4, 6, 8)
        for vector in weyl_shell(square)
    )
    if len(vectors) != 26641 or len(set(vectors)) != len(vectors):
        raise AssertionError("unexpected independent a-shell cardinality")
    return vectors


@lru_cache(maxsize=None)
def seed_eligible_a(b_square: int) -> tuple[Vector, ...]:
    seed = SHELL_ORBIT_SEEDS[b_square][0]
    gram_b = gram_vector(seed)
    eligible = []
    for a_vector in a_vectors():
        a_square = vector_square(a_vector)
        product = inner(a_vector, gram_b)
        square = (a_square + b_square, b_square + 2 * product)
        conjugate_square = conjugate(square)
        if sign(conjugate_square) <= 0:
            continue
        if sign(add((2, 0), neg(conjugate_square))) <= 0:
            continue
        eligible.append(a_vector)
    if len(eligible) != EXPECTED_SEED_ELIGIBLE_COUNTS[b_square]:
        raise AssertionError(
            f"seed b^2={b_square} has {len(eligible)} eligible a-vectors"
        )
    return tuple(sorted(eligible))


def digest_record(hasher: object, record: object) -> None:
    hasher.update(repr(record).encode("ascii"))
    hasher.update(b"\n")


def is_signature_admissible(projected_square_twice: Pair) -> bool:
    if sign(projected_square_twice) <= 0:
        return False
    schur_twice = add((4, 0), neg(projected_square_twice))
    return sign(schur_twice) > 0 and sign(conjugate(schur_twice)) > 0


@lru_cache(maxsize=None)
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


def iter_independent_mixed_records() -> Iterator[MixedRecord]:
    """Yield all eligible records in production's canonical ``(b,a)`` order."""

    transports: dict[Vector, tuple[int, Matrix]] = {}
    for b_square in (2, 4, 6):
        for b_vector, matrix in transport_orbit(b_square).items():
            if b_vector in transports:
                raise AssertionError("distinct b-shells overlap")
            transports[b_vector] = (b_square, matrix)
    if len(transports) != 9120:
        raise AssertionError("unexpected independent b-shell cardinality")

    a_shell_set = set(a_vectors())
    for b_vector in sorted(transports):
        b_square, matrix = transports[b_vector]
        gram_b = gram_vector(b_vector)
        b_on_interior = sum(
            gram_b[i] * INTERIOR_E8_SCALED[i] for i in range(8)
        )
        b_on_ray = sum(
            gram_b[i] * RAY_E8_PHI_COEFFICIENTS_SCALED[i]
            for i in range(8)
        )
        transported_a = tuple(
            sorted(matrix_vector(matrix, vector) for vector in seed_eligible_a(b_square))
        )
        if len(set(transported_a)) != len(transported_a):
            raise AssertionError("Weyl transport is not injective on eligible a")
        for a_vector in transported_a:
            if a_vector not in a_shell_set:
                raise AssertionError("Weyl transport left the complete a-shell box")
            a_square = vector_square(a_vector)
            gram_a = gram_vector(a_vector)
            product = inner(a_vector, gram_b)
            square = (a_square + b_square, b_square + 2 * product)
            conjugate_square = conjugate(square)
            if sign(conjugate_square) <= 0 or sign(
                add((2, 0), neg(conjugate_square))
            ) <= 0:
                raise AssertionError("transported record failed its exact leaf gate")
            if sign(square) <= 0:
                raise AssertionError("identity-place E8 square is not positive")
            covector = tuple((gram_a[i], gram_b[i]) for i in range(8))
            on_interior = (
                sum(
                    gram_a[i] * INTERIOR_E8_SCALED[i] for i in range(8)
                ),
                b_on_interior,
            )
            on_ray = (
                b_on_ray,
                sum(
                    (gram_a[i] + gram_b[i])
                    * RAY_E8_PHI_COEFFICIENTS_SCALED[i]
                    for i in range(8)
                ),
            )
            yield (
                scale(square, 2),
                covector,
                on_interior,
                on_ray,
                a_vector,
                b_vector,
            )


@lru_cache(maxsize=None)
def run_independent_replay(coefficient_bound: int = 30) -> dict[str, object]:
    eligible_hasher = hashlib.sha256()
    eligible_count = 0
    by_square: Counter[Pair] = Counter()
    a_square_histogram: Counter[int] = Counter()
    b_square_histogram: Counter[int] = Counter()
    for record in iter_independent_mixed_records():
        eligible_count += 1
        by_square[record[0]] += 1
        a_square_histogram[vector_square(record[4])] += 1
        b_square_histogram[vector_square(record[5])] += 1
        digest_record(
            eligible_hasher,
            (record[4], record[5], record[0], record[1], record[2], record[3]),
        )

    coefficients = coefficient_records(coefficient_bound)
    admissible_by_square: dict[Pair, tuple[CoefficientRecord, ...]] = {}
    signature_parameter_count = 0
    signature_count_by_square: dict[str, int] = {}
    for e8_square_twice in sorted(by_square):
        admissible = tuple(
            data
            for data in coefficients
            if is_signature_admissible(add(data[2], e8_square_twice))
        )
        admissible_by_square[e8_square_twice] = admissible
        signature_parameter_count += len(admissible)
        signature_count_by_square[repr(e8_square_twice)] = len(admissible)

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
    ) in iter_independent_mixed_records():
        for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
            admissible_by_square[e8_square_twice]
        ):
            projected_square_twice = add(non_e8_square, e8_square_twice)
            if not is_signature_admissible(projected_square_twice):
                raise AssertionError("an enumerated wall failed the signature gate")
            on_interior = add(non_e8_interior, e8_on_interior)
            on_ray = add(non_e8_ray, e8_on_ray)
            signs = sign(on_interior), sign(on_ray)
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
        "replay": "standalone_E8_Weyl_orbit_transport",
        "imports_project_code": False,
        "coefficient_bound": coefficient_bound,
        "necessity_bound": (
            "a^2<4+12*(phi')^2=28-12*phi<10; E8 even, hence "
            "a^2 in {0,2,4,6,8}"
        ),
        "a_vector_count": 26641,
        "b_vector_count": 9120,
        "raw_a_b_pair_count_before_Weyl_quotient": 26641 * 9120,
        "seed_pair_test_count_after_Weyl_quotient": 26641 * 3,
        "theta_pinned_shell_counts": {
            str(square): EXPECTED_SHELL_COUNTS[square]
            for square in (2, 4, 6, 8)
        },
        "shell_Weyl_orbit_counts": {
            str(square): list(EXPECTED_ORBIT_COUNTS[square])
            for square in (2, 4, 6, 8)
        },
        "eligible_a_count_at_b_seed_by_b_square": {
            str(square): EXPECTED_SEED_ELIGIBLE_COUNTS[square]
            for square in (2, 4, 6)
        },
        "eligible_w_count": eligible_count,
        "eligible_w_square_class_count": len(by_square),
        "eligible_w_square_class_histogram": {
            repr(key): by_square[key] for key in sorted(by_square)
        },
        "eligible_w_a_square_histogram": {
            str(key): a_square_histogram[key] for key in sorted(a_square_histogram)
        },
        "eligible_w_b_square_histogram": {
            str(key): b_square_histogram[key] for key in sorted(b_square_histogram)
        },
        "eligible_w_records_sha256": eligible_hasher.hexdigest(),
        "signature_admissible_parameter_pair_count_summed_over_w_square_classes": signature_parameter_count,
        "signature_parameter_count_by_w_square": signature_count_by_square,
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
    parser.add_argument("--compare-certificate", type=Path, default=None)
    args = parser.parse_args()
    replay = run_independent_replay(args.coefficient_bound)
    if args.compare_certificate is not None:
        payload = json.loads(args.compare_certificate.read_text())
        if replay != payload["independent_replay"]:
            raise AssertionError("independent shell3 replay differs from certificate")
        print(json.dumps({"independent_replay_matches": True}, sort_keys=True))
        return
    print(json.dumps(replay, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
