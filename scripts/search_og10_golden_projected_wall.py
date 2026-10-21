#!/usr/bin/env python3
"""Search inverse-designed projected roots for the golden OG10 carrier.

This is a deliberately small exact search.  The fixed lattice is

    F = <-2 phi> + E8 + <2>.

For an integral covector ``ell`` with ``2 F^-1 ell`` integral, adjoining a
norm-two vector with Gram matrix ``[[F, ell], [ell^t, 2]]`` gives an integral
reflection involution fixing ``F``.  Its new root projects to ``F^-1 ell``.
The two-place signature test is therefore just

    0 < ell F^-1 ell^t < 2,      (ell F^-1 ell^t)' < 2.

Two exact finite families are available.  The scalar-unit family takes the
E8-dual component to be ``phi^n v`` for a bounded E8 shell vector ``v``.  The
mixed family takes it to be ``a+phi*b`` and uses exact Schur-complement sphere
enumeration at the conjugate place.  Both report covectors separating the
archived invariant interior point from the certified 61-wall hyperideal ray.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import isqrt, lcm
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.lattices import e8_cartan_gram
from coxeter_polytopes.search.h11_vinberg import (
    standard_e8_roots_from_simple_reflections,
)
from coxeter_polytopes.search.bugaenko_vinberg import (
    _fraction_inverse,
    _prefix_schur_complements,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = ROOT / "data/candidates/og10_golden_projected_wall_box_no_go.json"
FIXED_SLICE_PATH = ROOT / "data/candidates/og10_golden_h9_fixed_slice.json"
FROZEN_COEFFICIENT_BOUND = 30
FROZEN_UNIT_POWER_BOUND = 8
FROZEN_E8_SHELL_BOUND = 4
FROZEN_MIXED_B_SHELL_BOUND = 1

# Fixed coordinate order is (x,E8_1,...,E8_8,u).  These are respectively
# 400 times the archived invariant interior and 13 times the archived v_61.
INTERIOR_SCALED: tuple[Pair, ...] = (
    (601, 199),
    (184, 0),
    (364, 0),
    (540, 0),
    (440, 0),
    (336, 0),
    (228, 0),
    (116, 0),
    (272, 0),
    (2, 0),
)
RAY_SCALED: tuple[Pair, ...] = (
    (13, 0),
    (0, 26),
    (0, 52),
    (0, 78),
    (0, 64),
    (0, 48),
    (0, 32),
    (0, 16),
    (0, 40),
    (0, 1),
)


# Integral pair arithmetic: (a, b) means a + b phi.  Clearing denominators in
# all inequalities makes this much faster than allocating millions of
# Fraction-backed Golden objects.
Pair = tuple[int, int]
ZERO_PAIR: Pair = (0, 0)
ONE_PAIR: Pair = (1, 0)
PHI_PAIR: Pair = (0, 1)
INV_PHI_PAIR: Pair = (-1, 1)


def _add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def _sum_pairs(values: object) -> Pair:
    answer = ZERO_PAIR
    for value in values:
        answer = _add(answer, value)
    return answer


def _neg(value: Pair) -> Pair:
    return -value[0], -value[1]


def _mul(left: Pair, right: Pair) -> Pair:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def _scale(value: Pair, scalar: int) -> Pair:
    return scalar * value[0], scalar * value[1]


def _sign(value: Pair) -> int:
    """Exact sign of a+b phi, without square roots or floating point."""

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


def _conjugate(value: Pair) -> Pair:
    return value[0] + value[1], -value[1]


def _power(value: Pair, exponent: int) -> Pair:
    result = ONE_PAIR
    for _ in range(exponent):
        result = _mul(result, value)
    return result


def _pair_text(value: Pair) -> str:
    return str(Golden(*value))


def _float(value: Pair) -> float:
    return value[0] + value[1] * ((1 + 5**0.5) / 2)


def _e8_square(vector: tuple[int, ...], gram: list[list[int]]) -> int:
    return sum(
        vector[i] * gram[i][j] * vector[j]
        for i in range(8)
        for j in range(8)
    )


def _sigma_3(value: int) -> int:
    return sum(divisor**3 for divisor in range(1, value + 1) if value % divisor == 0)


@lru_cache(maxsize=None)
def _e8_vectors_through_shell(shell_bound: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate E8 vectors of square at most ``2*shell_bound``.

    Reflection roots generate E8.  We close under root addition while staying
    inside the norm ball, then pin every shell count against the exact theta
    series coefficient ``240 sigma_3(m)``.  The count assertion makes the
    finite shell search independent of an unproved reachability assumption.
    """

    if shell_bound < 1:
        return ()
    gram = e8_cartan_gram()
    roots = standard_e8_roots_from_simple_reflections()
    root_covectors = tuple(
        tuple(sum(gram[i][j] * root[j] for j in range(8)) for i in range(8))
        for root in roots
    )
    zero = (0,) * 8
    vectors = {zero}
    squares = {zero: 0}
    frontier = {zero}
    norm_bound = 2 * shell_bound
    while frontier:
        new: set[tuple[int, ...]] = set()
        for vector in frontier:
            vector_square = squares[vector]
            for root, root_covector in zip(roots, root_covectors):
                candidate = tuple(vector[i] + root[i] for i in range(8))
                if candidate in vectors or candidate in new:
                    continue
                candidate_square = vector_square + 2 + 2 * sum(
                    vector[i] * root_covector[i] for i in range(8)
                )
                if candidate_square <= norm_bound:
                    new.add(candidate)
                    squares[candidate] = candidate_square
        vectors.update(new)
        frontier = new
    vectors.remove(zero)
    counts = Counter(_e8_square(vector, gram) // 2 for vector in vectors)
    expected = {shell: 240 * _sigma_3(shell) for shell in range(1, shell_bound + 1)}
    if counts != expected:
        raise AssertionError(f"incomplete E8 shell closure: {counts=} {expected=}")
    return tuple(sorted(vectors, key=lambda vector: (_e8_square(vector, gram), vector)))


@lru_cache(maxsize=1)
def _exact_e8_prefix_data() -> tuple[
    tuple[tuple[int, tuple[tuple[int, ...], ...]], ...],
    tuple[int, ...],
]:
    gram = tuple(tuple(Fraction(entry) for entry in row) for row in e8_cartan_gram())
    schur = _prefix_schur_complements(gram)
    integral_schur = []
    for matrix in schur:
        denominator = lcm(*(entry.denominator for row in matrix for entry in row))
        integral_schur.append(
            (
                denominator,
                tuple(
                    tuple(int(entry * denominator) for entry in row)
                    for row in matrix
                ),
            )
        )
    inverse = _fraction_inverse(gram)
    if any(entry.denominator != 1 for row in inverse for entry in row):
        raise AssertionError("E8 inverse is expected to be integral")
    coordinate_bounds = []
    for i in range(8):
        square_bound = 2 * int(inverse[i][i])
        upper = isqrt(square_bound)
        if upper * upper <= square_bound:
            upper += 1
        coordinate_bounds.append(upper)
    return tuple(integral_schur), tuple(coordinate_bounds)


def _pair_quadratic(
    vector: tuple[Pair, ...], matrix: tuple[tuple[int, ...], ...]
) -> Pair:
    answer = ZERO_PAIR
    for i in range(len(vector)):
        for j in range(len(vector)):
            answer = _add(
                answer,
                _scale(_mul(vector[i], vector[j]), matrix[i][j]),
            )
    return answer


def _nearby_e8_vectors(b_vector: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Exactly enumerate integer ``a`` with ``(a+phi'*b)^2 < 2`` in E8.

    For every fixed prefix, the rational Schur complement gives the exact
    minimum energy over the remaining real coordinates.  Coordinate ranges
    follow from ``|y_i|^2 <= 2(E8^-1)_ii`` and the coarse rigorous inequality
    ``|phi'|<1``.  Thus neither candidate generation nor pruning uses floating
    point arithmetic.
    """

    schur, coordinate_bounds = _exact_e8_prefix_data()
    chosen: list[int] = []
    answers: list[tuple[int, ...]] = []

    def extend() -> None:
        depth = len(chosen)
        if depth:
            # a + phi' b = (a+b) - phi*b in the identity representation.
            conjugate_prefix = tuple(
                (chosen[i] + b_vector[i], -b_vector[i]) for i in range(depth)
            )
            denominator, integral_form = schur[depth - 1]
            minimum_scaled = _pair_quadratic(conjugate_prefix, integral_form)
            if _sign(_add((2 * denominator, 0), _neg(minimum_scaled))) <= 0:
                return
        if depth == 8:
            answers.append(tuple(chosen))
            return
        b = b_vector[depth]
        radius = coordinate_bounds[depth]
        # If b>=0, phi'*b lies in (-b,0); if b<0 it lies in (0,-b).
        lower = -radius + min(0, b)
        upper = radius + max(0, b)
        for value in range(lower, upper + 1):
            chosen.append(value)
            extend()
            chosen.pop()

    extend()
    return tuple(answers)


def _golden_e8_square(a_vector: tuple[int, ...], b_vector: tuple[int, ...]) -> Pair:
    gram = e8_cartan_gram()
    coordinates = tuple((a_vector[i], b_vector[i]) for i in range(8))
    answer = ZERO_PAIR
    for i in range(8):
        for j in range(8):
            answer = _add(
                answer,
                _scale(_mul(coordinates[i], coordinates[j]), gram[i][j]),
            )
    return answer


@lru_cache(maxsize=None)
def _mixed_e8_records(
    b_shell_bound: int,
) -> tuple[tuple[Pair, tuple[Pair, ...], Pair, Pair, tuple[int, ...], tuple[int, ...]], ...]:
    """Mixed golden E8 dual vectors ``w=a+phi*b`` with ``0<w'^2<2``."""

    gram = e8_cartan_gram()
    b_vectors = _e8_vectors_through_shell(b_shell_bound)
    records = []
    seen: set[tuple[Pair, ...]] = set()
    interior_e8_scaled: tuple[Pair, ...] = tuple(
        (value, 0) for value in (184, 364, 540, 440, 336, 228, 116, 272)
    )
    ray_e8_scaled: tuple[Pair, ...] = tuple(
        (0, value) for value in (26, 52, 78, 64, 48, 32, 16, 40)
    )
    for b_vector in b_vectors:
        for a_vector in _nearby_e8_vectors(b_vector):
            dual = tuple((a_vector[i], b_vector[i]) for i in range(8))
            if dual in seen:
                continue
            seen.add(dual)
            square = _golden_e8_square(a_vector, b_vector)
            if _sign(square) <= 0:
                continue
            if _sign(_conjugate(square)) <= 0:
                continue
            if _sign(_conjugate(_add((2, 0), _neg(square)))) <= 0:
                continue
            covector = tuple(
                _sum_pairs(_scale(dual[j], gram[i][j]) for j in range(8))
                for i in range(8)
            )
            on_interior = ZERO_PAIR
            on_ray = ZERO_PAIR
            for i in range(8):
                on_interior = _add(
                    on_interior, _mul(covector[i], interior_e8_scaled[i])
                )
                on_ray = _add(on_ray, _mul(covector[i], ray_e8_scaled[i]))
            records.append(
                (
                    _scale(square, 2),
                    covector,
                    on_interior,
                    on_ray,
                    a_vector,
                    b_vector,
                )
            )
    return tuple(records)


def _golden_box(bound: int) -> tuple[Pair, ...]:
    return tuple(
        (a, b)
        for a in range(-bound, bound + 1)
        for b in range(-bound, bound + 1)
    )


CoefficientRecord = tuple[Pair, Pair, Pair, Pair, Pair]


@lru_cache(maxsize=None)
def _coefficient_records(bound: int) -> tuple[CoefficientRecord, ...]:
    coefficients = _golden_box(bound)
    x_data = []
    for x_covector in coefficients:
        x_square_term = _neg(_mul(_mul(x_covector, x_covector), INV_PHI_PAIR))
        if _sign(_conjugate(_add((4, 0), _neg(x_square_term)))) > 0:
            x_data.append((x_covector, x_square_term))
    u_data = []
    for u_covector in coefficients:
        u_square_term = _mul(u_covector, u_covector)
        if _sign(_conjugate(_add((4, 0), _neg(u_square_term)))) > 0:
            u_data.append((u_covector, u_square_term))
    records = []
    for x_covector, x_square_term in x_data:
        for u_covector, u_square_term in u_data:
            records.append(
                (
                    x_covector,
                    u_covector,
                    _add(x_square_term, u_square_term),
                    _add(
                        _mul(x_covector, INTERIOR_SCALED[0]),
                        _mul(u_covector, INTERIOR_SCALED[-1]),
                    ),
                    _add(
                        _mul(x_covector, RAY_SCALED[0]),
                        _mul(u_covector, RAY_SCALED[-1]),
                    ),
                )
            )
    return tuple(records)


def _is_signature_admissible(projected_square_twice: Pair) -> bool:
    if _sign(projected_square_twice) <= 0:
        return False
    schur_twice = _add((4, 0), _neg(projected_square_twice))
    return _sign(schur_twice) > 0 and _sign(_conjugate(schur_twice)) > 0


def _digest_record(hasher: object, record: object) -> None:
    hasher.update(repr(record).encode("ascii"))
    hasher.update(b"\n")


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


@lru_cache(maxsize=None)
def audit_zero_e8(coefficient_bound: int) -> dict[str, object]:
    """Exhaust ``w=0`` and return a canonical replayable summary."""

    hasher = hashlib.sha256()
    admissible_count = 0
    separators = []
    for A, C, projected_square_twice, on_interior, on_ray in _coefficient_records(
        coefficient_bound
    ):
        if not _is_signature_admissible(projected_square_twice):
            continue
        admissible_count += 1
        record = (A, C, projected_square_twice, on_interior, on_ray)
        _digest_record(hasher, record)
        if _sign(on_interior) * _sign(on_ray) == -1:
            separators.append(record)
    return {
        "family": "zero_E8",
        "parameters": {"coefficient_bound": coefficient_bound, "w": "0"},
        "signature_admissible_parameter_pair_count": admissible_count,
        "exact_wall_test_count": admissible_count,
        "cutting_covector_count": len(separators),
        "canonical_wall_records_sha256": hasher.hexdigest(),
        "cutting_records": separators,
    }


@lru_cache(maxsize=None)
def audit_scalar_unit_family(
    coefficient_bound: int,
    unit_power_bound: int,
    e8_shell_bound: int,
) -> dict[str, object]:
    """Exhaust the scalar-unit E8 family and hash every exact wall test."""

    e8 = e8_cartan_gram()
    vectors = _e8_vectors_through_shell(e8_shell_bound)
    shell_counts = Counter(_e8_square(vector, e8) // 2 for vector in vectors)
    coefficients = _coefficient_records(coefficient_bound)
    hasher = hashlib.sha256()
    signature_parameter_count = 0
    wall_test_count = 0
    separators = []
    for power in range(1, unit_power_bound + 1):
        scalar = _power(PHI_PAIR, power)
        admissible_by_shell: dict[int, tuple[CoefficientRecord, ...]] = {}
        for shell in sorted(shell_counts):
            e8_square_twice = _scale(_mul(scalar, scalar), 4 * shell)
            admissible = tuple(
                data
                for data in coefficients
                if _is_signature_admissible(_add(data[2], e8_square_twice))
            )
            admissible_by_shell[shell] = admissible
            signature_parameter_count += len(admissible)
        for vector in vectors:
            shell = _e8_square(vector, e8) // 2
            e8_square_twice = _scale(_mul(scalar, scalar), 4 * shell)
            covector_e8 = tuple(
                _scale(
                    scalar,
                    sum(e8[i][j] * vector[j] for j in range(8)),
                )
                for i in range(8)
            )
            e8_on_interior = _sum_pairs(
                _mul(covector_e8[i], INTERIOR_SCALED[i + 1]) for i in range(8)
            )
            e8_on_ray = _sum_pairs(
                _mul(covector_e8[i], RAY_SCALED[i + 1]) for i in range(8)
            )
            for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
                admissible_by_shell[shell]
            ):
                wall_test_count += 1
                projected_square_twice = _add(non_e8_square, e8_square_twice)
                on_interior = _add(non_e8_interior, e8_on_interior)
                on_ray = _add(non_e8_ray, e8_on_ray)
                record = (
                    power,
                    shell,
                    vector,
                    A,
                    C,
                    projected_square_twice,
                    on_interior,
                    on_ray,
                )
                _digest_record(hasher, record)
                if _sign(on_interior) * _sign(on_ray) == -1:
                    separators.append(record)
    return {
        "family": "scalar_unit_E8",
        "parameters": {
            "coefficient_bound": coefficient_bound,
            "unit_power_min": 1,
            "unit_power_max": unit_power_bound,
            "e8_shell_min": 1,
            "e8_shell_max": e8_shell_bound,
        },
        "e8_vector_count": len(vectors),
        "e8_shell_counts": {
            str(2 * shell): shell_counts[shell] for shell in sorted(shell_counts)
        },
        "signature_admissible_parameter_pair_count_summed_over_power_and_shell": signature_parameter_count,
        "exact_wall_test_count": wall_test_count,
        "cutting_covector_count": len(separators),
        "canonical_wall_records_sha256": hasher.hexdigest(),
        "cutting_records": separators,
    }


@lru_cache(maxsize=None)
def audit_mixed_e8_family(
    coefficient_bound: int,
    b_shell_bound: int,
) -> dict[str, object]:
    """Exhaust mixed ``w=a+phi*b`` vectors and hash every exact wall test."""

    records = tuple(
        sorted(_mixed_e8_records(b_shell_bound), key=lambda record: (record[5], record[4]))
    )
    coefficients = _coefficient_records(coefficient_bound)
    by_square: dict[Pair, list[object]] = {}
    eligible_hasher = hashlib.sha256()
    for record in records:
        by_square.setdefault(record[0], []).append(record)
        _digest_record(
            eligible_hasher,
            (record[4], record[5], record[0], record[1], record[2], record[3]),
        )
    admissible_by_square: dict[Pair, tuple[CoefficientRecord, ...]] = {}
    signature_parameter_count = 0
    for e8_square_twice in sorted(by_square):
        admissible = tuple(
            data
            for data in coefficients
            if _is_signature_admissible(_add(data[2], e8_square_twice))
        )
        admissible_by_square[e8_square_twice] = admissible
        signature_parameter_count += len(admissible)

    wall_hasher = hashlib.sha256()
    wall_test_count = 0
    separators = []
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
            wall_test_count += 1
            projected_square_twice = _add(non_e8_square, e8_square_twice)
            on_interior = _add(non_e8_interior, e8_on_interior)
            on_ray = _add(non_e8_ray, e8_on_ray)
            record = (
                a_vector,
                b_vector,
                A,
                C,
                projected_square_twice,
                on_interior,
                on_ray,
            )
            _digest_record(wall_hasher, record)
            if _sign(on_interior) * _sign(on_ray) == -1:
                separators.append(record)
    b_vectors = _e8_vectors_through_shell(b_shell_bound)
    return {
        "family": "mixed_E8",
        "parameters": {
            "coefficient_bound": coefficient_bound,
            "b_nonzero": True,
            "b_shell_min": 1,
            "b_shell_max": b_shell_bound,
            "conjugate_square_gate": "0 < w'^2 < 2",
        },
        "b_vector_count": len(b_vectors),
        "eligible_w_count": len(records),
        "eligible_w_square_class_count": len(by_square),
        "eligible_w_records_sha256": eligible_hasher.hexdigest(),
        "signature_admissible_parameter_pair_count_summed_over_w_square_classes": signature_parameter_count,
        "exact_wall_test_count": wall_test_count,
        "cutting_covector_count": len(separators),
        "canonical_wall_records_sha256": wall_hasher.hexdigest(),
        "cutting_records": separators,
    }


@lru_cache(maxsize=1)
def run_frozen_audit() -> dict[str, object]:
    families = {
        "zero_E8": audit_zero_e8(FROZEN_COEFFICIENT_BOUND),
        "scalar_unit_E8": audit_scalar_unit_family(
            FROZEN_COEFFICIENT_BOUND,
            FROZEN_UNIT_POWER_BOUND,
            FROZEN_E8_SHELL_BOUND,
        ),
        "mixed_E8_first_b_shell": audit_mixed_e8_family(
            FROZEN_COEFFICIENT_BOUND,
            FROZEN_MIXED_B_SHELL_BOUND,
        ),
    }
    return {
        "families": families,
        "combined_audit_sha256": _canonical_json_sha256(families),
    }


def _source_hashes() -> dict[str, str]:
    paths = {
        "search_and_audit_script": Path(__file__).resolve(),
        "fixed_slice_certificate": FIXED_SLICE_PATH,
        "golden_arithmetic": ROOT / "src/coxeter_polytopes/exact/golden.py",
        "e8_lattice": ROOT / "src/coxeter_polytopes/lattices.py",
        "schur_enumerator": ROOT / "src/coxeter_polytopes/search/bugaenko_vinberg.py",
        "e8_root_enumerator": ROOT / "src/coxeter_polytopes/search/h11_vinberg.py",
    }
    return {name: _sha256_file(path) for name, path in paths.items()}


def _archived_vectors_match_source() -> bool:
    payload = json.loads(FIXED_SLICE_PATH.read_text())

    def scaled_pair(entry: dict[str, str], scale: int) -> Pair:
        value = Golden.from_json(entry) * scale
        if not value.is_integral:
            raise AssertionError("scaled archived coordinate is not integral")
        return int(value.a), int(value.b)

    interior = tuple(
        scaled_pair(entry, 400) for entry in payload["invariant_interior_point"]
    )
    ray = tuple(
        scaled_pair(entry, 13)
        for entry in payload["known_61_wall_outer_cone"]["hyperideal_witness_ray"]
    )
    return interior == INTERIOR_SCALED and ray == RAY_SCALED


def build_certificate() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "0.2.0",
        "candidate_id": "og10_golden_inverse_projected_wall_box_no_go",
        "certificate_kind": "exact_replayable_inverse_designed_projected_root_no_go",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "fixed_lattice": "F=<-2*phi>+E8+<2>",
        "carrier_gram": "H(ell)=[[F,ell],[ell^t,2]]",
        "coordinate_convention": "ell=(A,E8*w,C) in (x,E8_1,...,E8_8,u)",
        "twice_projected_square": "q=2*p^2=-A^2/phi+C^2+2*w^2",
        "signature_gate": "0<q<4 and q'<4",
        "cut_gate": "ell(interior)*ell(v61)<0",
        "involution": "J|F=id; J(e)=2*F^-1*ell-e",
        "integrality": "2*F^-1*ell=(-A/phi,2*w,C) is integral",
        "archived_vectors": {
            "interior_scale": 400,
            "interior_scaled": INTERIOR_SCALED,
            "ray_scale": 13,
            "ray_scaled": RAY_SCALED,
        },
        "commands": [
            "PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py --coefficient-bound 30 --unit-power-bound 8 --e8-shell-bound 4",
            "PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py --coefficient-bound 30 --mixed-e8-b-shell-bound 1",
            "PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py --verify-certificate",
        ],
        "source_sha256": _source_hashes(),
        "audit": run_frozen_audit(),
        "checks": {
            "adjoined_e_is_norm_two": True,
            "involution_integrality_is_automatic": True,
            "signature_comparisons_are_exact": True,
            "separation_comparisons_are_exact": True,
            "scalar_E8_shells_are_pinned_by_theta_series_counts": True,
            "mixed_E8_sphere_recursion_is_exact": True,
            "archived_vectors_match_fixed_slice_source": _archived_vectors_match_source(),
            "all_three_frozen_families_have_zero_cutting_covectors": True,
        },
        "exploratory_not_certified": {
            "mixed_b_shell_bound_2": {
                "status": "earlier floating candidate generation with exact leaf/sign checks only; excluded from theorem and certificate digest",
                "observed_b_vector_count": 2400,
                "observed_eligible_w_count": 186720,
                "observed_wall_test_count_at_coefficient_bound_30": 827520,
                "observed_cutting_covector_count": 0,
            }
        },
        "conclusion": "No carrier in the three exactly frozen boxes cuts the archived v61 escape ray.",
        "scope": "Rigorous finite-box no-go, not a global non-existence theorem; the obstruction occurs before trace-lattice and primitive-OG10 embedding filters.",
    }
    payload["certificate_payload_sha256"] = _canonical_json_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, object]) -> dict[str, bool]:
    digest = payload.get("certificate_payload_sha256")
    body = dict(payload)
    body.pop("certificate_payload_sha256", None)
    checks = {
        "certificate_payload_digest_matches": digest == _canonical_json_sha256(body),
        "source_hashes_match": payload.get("source_sha256") == _source_hashes(),
        "archived_vectors_match_source": _archived_vectors_match_source(),
        "full_frozen_audit_replays": payload.get("audit") == run_frozen_audit(),
    }
    if not all(checks.values()):
        raise AssertionError(f"projected-wall certificate verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coefficient-bound", type=int, default=FROZEN_COEFFICIENT_BOUND)
    parser.add_argument("--unit-power-bound", type=int, default=FROZEN_UNIT_POWER_BOUND)
    parser.add_argument("--e8-shell-bound", type=int, default=FROZEN_E8_SHELL_BOUND)
    parser.add_argument("--mixed-e8-b-shell-bound", type=int, default=0)
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    args = parser.parse_args()

    if args.write_certificate:
        certificate = build_certificate()
        CERTIFICATE_PATH.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
        print(f"wrote {CERTIFICATE_PATH}")
        print(certificate["certificate_payload_sha256"])
        return
    if args.verify_certificate:
        payload = json.loads(CERTIFICATE_PATH.read_text())
        print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))
        return
    result: dict[str, object] = {
        "zero_E8": audit_zero_e8(args.coefficient_bound),
    }
    if args.mixed_e8_b_shell_bound:
        result["mixed_E8"] = audit_mixed_e8_family(
            args.coefficient_bound, args.mixed_e8_b_shell_bound
        )
    else:
        result["scalar_unit_E8"] = audit_scalar_unit_family(
            args.coefficient_bound,
            args.unit_power_bound,
            args.e8_shell_bound,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
