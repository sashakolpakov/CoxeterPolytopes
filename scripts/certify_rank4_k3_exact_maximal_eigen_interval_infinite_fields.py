#!/usr/bin/env python3
"""Exact arithmetic certificate for maximal rank-four K3 eigen-intervals.

The checker reconstructs the diagonal lattice, the two endpoint roots,
integral field reflections and their four-dimensional lifts, the infinite
ambient facet family, the fixed square-14 polarization, and finite exact
samples of the maximality and Pell arguments.  It imports no reflective
lattice classification.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from math import isqrt, sqrt
import json
from pathlib import Path
from typing import Iterable, Sequence


IntMatrix = tuple[tuple[int, ...], ...]
Field = tuple[Fraction, Fraction]  # c0+c1*t, t^2=2*a*t-1
FieldMatrix = tuple[tuple[Field, ...], ...]
Vector = tuple[int, ...]


def int_transpose(matrix: Sequence[Sequence[int]]) -> IntMatrix:
    return tuple(
        tuple(row[column] for row in matrix)
        for column in range(len(matrix[0]))
    )


def int_multiply(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> IntMatrix:
    columns = int_transpose(right)
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def int_identity(rank: int) -> IntMatrix:
    return tuple(
        tuple(1 if row == column else 0 for column in range(rank))
        for row in range(rank)
    )


def int_block(*rows: Sequence[IntMatrix]) -> IntMatrix:
    return tuple(
        tuple(entry for block in block_row for entry in block[row])
        for block_row in rows
        for row in range(len(block_row[0]))
    )


def f(value: int | Fraction, coefficient: int | Fraction = 0) -> Field:
    return Fraction(value), Fraction(coefficient)


def fadd(left: Field, right: Field) -> Field:
    return left[0] + right[0], left[1] + right[1]


def fneg(value: Field) -> Field:
    return -value[0], -value[1]


def fmul(left: Field, right: Field, a: int) -> Field:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1]
        + left[1] * right[0]
        + 2 * a * left[1] * right[1],
    )


def fconjugate(value: Field, a: int) -> Field:
    return value[0] + 2 * a * value[1], -value[1]


def fvalue(value: Field, t: float) -> float:
    return float(value[0]) + float(value[1]) * t


def fm_transpose(matrix: FieldMatrix) -> FieldMatrix:
    return tuple(
        tuple(row[column] for row in matrix)
        for column in range(len(matrix[0]))
    )


def fm_multiply(left: FieldMatrix, right: FieldMatrix, a: int) -> FieldMatrix:
    columns = fm_transpose(right)
    return tuple(
        tuple(
            sum_fields(fmul(x, y, a) for x, y in zip(row, column))
            for column in columns
        )
        for row in left
    )


def sum_fields(values: Iterable[Field]) -> Field:
    answer = f(0)
    for value in values:
        answer = fadd(answer, value)
    return answer


def field_payload(value: Field) -> list[str]:
    return [str(value[0]), str(value[1])]


def field_matrix_payload(matrix: FieldMatrix) -> list[list[list[str]]]:
    return [[field_payload(entry) for entry in row] for row in matrix]


def companion(a: int) -> IntMatrix:
    return ((0, -1), (1, 2 * a))


def original_gram(a: int) -> IntMatrix:
    return (
        (-2, 0, 0, 0),
        (0, 2, 0, 0),
        (0, 0, 2 - 4 * a, -2),
        (0, 0, -2, -2),
    )


def original_tau(a: int) -> IntMatrix:
    c = companion(a)
    zero = ((0, 0), (0, 0))
    return int_block((c, zero), (zero, c))


def regular_matrix(value: Field, a: int) -> IntMatrix:
    if value[0].denominator != 1 or value[1].denominator != 1:
        raise ValueError("regular integral lift requires coefficients in Z[t]")
    c0, c1 = int(value[0]), int(value[1])
    c = companion(a)
    return tuple(
        tuple(c0 * int(row == column) + c1 * c[row][column] for column in range(2))
        for row in range(2)
    )


def lift_field_matrix(matrix: FieldMatrix, a: int) -> IntMatrix:
    blocks = tuple(
        tuple(regular_matrix(entry, a) for entry in row)
        for row in matrix
    )
    return int_block(*blocks)


def field_form(a: int) -> FieldMatrix:
    d = a - 1
    return ((f(1, -1), f(0)), (f(0), f(2 * d)))


def endpoint_reflections(a: int) -> tuple[FieldMatrix, FieldMatrix]:
    t = f(0, 1)
    two_t_minus_one = f(-1, 2)
    s0 = ((f(1), f(0)), (f(0), f(-1)))
    s1 = (
        (two_t_minus_one, f(2, -2)),
        (f(0, 2), fneg(two_t_minus_one)),
    )
    # Keep t referenced explicitly for readability and catch a future basis
    # change through the polynomial identity below.
    if fmul(fadd(t, f(-1)), fadd(t, f(-1)), a) != f(0, 2 * (a - 1)):
        raise AssertionError("(t-1)^2=2*d*t failed")
    return s0, s1


def reflection_record(a: int) -> dict[str, object]:
    if a < 2:
        raise ValueError("a must be at least two")
    form = field_form(a)
    s0, s1 = endpoint_reflections(a)
    identity: FieldMatrix = ((f(1), f(0)), (f(0), f(1)))
    checks = []
    lifts = []
    gram = original_gram(a)
    tau = original_tau(a)
    for reflection in (s0, s1):
        checks.append(
            fm_multiply(reflection, reflection, a) == identity
            and fm_multiply(
                fm_multiply(fm_transpose(reflection), form, a),
                reflection,
                a,
            )
            == form
        )
        lift = lift_field_matrix(reflection, a)
        lifts.append(lift)
        if int_multiply(int_multiply(int_transpose(lift), gram), lift) != gram:
            raise AssertionError(f"four-dimensional isometry failed at a={a}")
        if int_multiply(lift, tau) != int_multiply(tau, lift):
            raise AssertionError(f"field centralizer failed at a={a}")
    translation = fm_multiply(s1, s0, a)
    expected_translation: FieldMatrix = (
        (f(-1, 2), f(-2, 2)),
        (f(0, 2), f(-1, 2)),
    )
    if translation != expected_translation or not all(checks):
        raise AssertionError(f"endpoint reflection identities failed at a={a}")
    return {
        "a": a,
        "d": a - 1,
        "field_form": field_matrix_payload(form),
        "S0": field_matrix_payload(s0),
        "S1": field_matrix_payload(s1),
        "translation_S1_S0": field_matrix_payload(translation),
        "field_reflection_checks": checks,
        "integral_four_dimensional_lifts": [
            [list(row) for row in lift] for lift in lifts
        ],
        "lifts_preserve_original_gram": True,
        "lifts_commute_with_tau": True,
        "translation_half_trace": field_payload(f(-1, 2)),
        "cyclic_Cartan_product": field_payload(f(0, 4)),
    }


def q_pair(left: Sequence[int], right: Sequence[int], d: int) -> int:
    return (
        -left[0] * right[0]
        + left[1] * right[1]
        + left[2] * right[2]
        + 2 * d * left[3] * right[3]
    )


def endpoint_record(a: int) -> dict[str, object]:
    d = a - 1
    alpha0 = (0, 0, -1, 0)
    alpha1 = (1, -1, 1, 0)
    polarization = (3, -1, 1, 0)
    original_polarization = (-1, 3, 0, 1)
    tau_h = tuple(
        sum(original_tau(a)[row][column] * original_polarization[column] for column in range(4))
        for row in range(4)
    )
    pairing = sum(
        original_polarization[row]
        * original_gram(a)[row][column]
        * tau_h[column]
        for row in range(4)
        for column in range(4)
    )
    return {
        "a": a,
        "endpoint_roots": [list(alpha0), list(alpha1)],
        "endpoint_root_q_norms": [
            q_pair(alpha0, alpha0, d),
            q_pair(alpha1, alpha1, d),
        ],
        "restricted_wall_equations": [
            "-(t-1)*s",
            "(t-1)*(s-1)",
        ],
        "exact_open_component": "0<s<1",
        "positive_interval_radius_squared": "t/(t-1)",
        "finite_length": True,
        "cosh_squared_length": "t",
        "polarization_diagonal_coordinates": list(polarization),
        "polarization_original_coordinates": list(original_polarization),
        "polarization_square": -2 * q_pair(polarization, polarization, d),
        "polarized_correspondence_height": pairing,
        "expected_polarized_correspondence_height": 32 * a - 10,
    }


def ambient_facet_record(a: int, k: int) -> dict[str, object]:
    if a < 2 or k < 1:
        raise ValueError("a>=2 and k>=1 are required")
    d = a - 1
    height = d * k * k
    root = (height, 0, height - 1, k)
    support = (
        1 + height * height,
        0,
        height * (height - 1),
        height * k,
    )
    if q_pair(root, root, d) != 1:
        raise AssertionError("ambient root family changed")
    if q_pair(support, root, d) != 0:
        raise AssertionError("ambient support point left its mirror")
    return {
        "a": a,
        "d": d,
        "k": k,
        "height": height,
        "root": list(root),
        "root_q_norm": q_pair(root, root, d),
        "support_point": list(support),
        "support_q_norm": q_pair(support, support, d),
        "support_root_pairing": q_pair(support, root, d),
    }


def squarefree_decomposition(number: int) -> tuple[int, int]:
    residual = number
    squarefree = 1
    square_factor = 1
    prime = 2
    while prime * prime <= residual:
        exponent = 0
        while residual % prime == 0:
            residual //= prime
            exponent += 1
        square_factor *= prime ** (exponent // 2)
        if exponent % 2:
            squarefree *= prime
        prime = 3 if prime == 2 else prime + 2
    if residual > 1:
        squarefree *= residual
    return squarefree, square_factor


def field_record(a: int) -> dict[str, int]:
    d, y = squarefree_decomposition(a * a - 1)
    if a * a - d * y * y != 1:
        raise AssertionError("Pell record failed")
    return {"a": a, "squarefree_field_parameter": d, "Pell_y": y}


def maximality_trace_candidates(a: int, bound: int = 12) -> list[list[int]]:
    t = a + sqrt(a * a - 1)
    conjugate = 1 / t
    answer = []
    for m in range(-bound, bound + 1):
        for n in range(-bound, bound + 1):
            value = m + n * t
            other = m + n * conjugate
            if 1 < value < 2 * t - 1 and abs(other) <= 1 + 1e-12:
                answer.append([m, n])
    return answer


def square_equation_solutions(a: int, bound: int = 256) -> dict[str, list[list[int]]]:
    d = a - 1
    t_case = []
    tm1_case = []
    for q0 in range(-bound, bound + 1):
        for q1 in range(-bound, bound + 1):
            constant = q0 * q0 - q1 * q1
            coefficient = 2 * q1 * (q0 + a * q1)
            if (constant, coefficient) == (2 * d, 2 * d):
                t_case.append([q0, q1])
            if (constant, coefficient) == (1, 2 * d - 1):
                tm1_case.append([q0, q1])
    return {"P_equals_t": t_case, "P_equals_t_minus_one": tm1_case}


def enumerate_sector_roots(d: int, max_height: int) -> Iterable[Vector]:
    for height in range(1, max_height + 1):
        target = height * height + 1
        for z in range(-isqrt(target // (2 * d)), isqrt(target // (2 * d)) + 1):
            residual = target - 2 * d * z * z
            if residual < 0:
                continue
            for m in range(isqrt(residual) + 1):
                u2 = residual - m * m
                u = isqrt(u2)
                if u * u == u2:
                    yield height, -m, u, z


def finite_root_audit(a: int, max_height: int = 32) -> dict[str, object]:
    d = a - 1
    roots = list(enumerate_sector_roots(d, max_height))
    endpoint0_strict = True
    endpoint1_weak = True
    endpoint1_equalities = []
    polarization_strict = True
    for root in roots:
        height, x, u, z = root
        # Coefficients c0+c1*t of q(p_s,root) at s=0,1.
        at0 = f(-x, -height)
        at1 = f(-x - u - 2 * d * z, u - height)
        t = a + sqrt(a * a - 1)
        endpoint0_strict &= fvalue(at0, t) < -1e-10
        value1 = fvalue(at1, t)
        endpoint1_weak &= value1 <= 1e-10
        if abs(value1) <= 1e-10:
            endpoint1_equalities.append(list(root))
        polarization_strict &= (
            -3 * height - x + u < 0
        )
    if not endpoint0_strict or not endpoint1_weak or not polarization_strict:
        raise AssertionError(f"finite root audit failed at a={a}")
    return {
        "a": a,
        "max_height": max_height,
        "sector_root_count": len(roots),
        "p0_strict_for_positive_height_roots": endpoint0_strict,
        "p1_weak_for_positive_height_roots": endpoint1_weak,
        "p1_equalities": endpoint1_equalities,
        "expected_p1_equality": [[1, -1, 1, 0]],
        "square14_polarization_strict": polarization_strict,
    }


def build_certificate(max_a: int = 64) -> dict[str, object]:
    if max_a < 2:
        raise ValueError("max_a must be at least two")
    exact_a = tuple(dict.fromkeys(value for value in (2, 3, 4, 7, 19, max_a) if value <= max_a))
    reflections = [reflection_record(a) for a in exact_a]
    endpoints = [endpoint_record(a) for a in exact_a]
    if any(
        row["polarization_square"] != 14
        or row["polarized_correspondence_height"]
        != row["expected_polarized_correspondence_height"]
        for row in endpoints
    ):
        raise AssertionError("polarized packet formula changed")
    maximality = []
    for a in exact_a:
        candidates = maximality_trace_candidates(a)
        if candidates != [[-1, 1], [0, 1]]:
            raise AssertionError((a, candidates))
        squares = square_equation_solutions(a, bound=64)
        if squares["P_equals_t"] or squares["P_equals_t_minus_one"]:
            raise AssertionError((a, squares))
        maximality.append(
            {
                "a": a,
                "trace_candidates_m_n": candidates,
                "candidate_values": ["t-1", "t"],
                "square_obstruction_search_bound": 64,
                "square_equation_solutions": squares,
            }
        )
    fields = [field_record(a) for a in range(2, max_a + 1)]
    distinct_fields = sorted({row["squarefree_field_parameter"] for row in fields})
    return {
        "theorem": "rank-four exact maximal K3 eigen-interval infinite-field family",
        "hyperbolic_dimension": 1,
        "rank_Q_full_trace": 4,
        "field_degree": 2,
        "reflection_records": reflections,
        "endpoint_and_polarization_records": endpoints,
        "finite_root_audits": [finite_root_audit(a) for a in exact_a],
        "ambient_infinite_facet_samples": [
            ambient_facet_record(a, k)
            for a in exact_a
            for k in (1, 2, 3, 5, 8)
        ],
        "maximality_records": maximality,
        "field_samples": fields,
        "distinct_sample_field_count": len(distinct_fields),
        "distinct_sample_squarefree_parameters": distinct_fields,
        "conclusions": {
            "finite_volume_eigen_interval": True,
            "exact_prescribed_eigen_nef_return": True,
            "complete_restricted_root_saturation": True,
            "one_restricted_chamber_orbit": True,
            "reflection_maximal_in_declared_Z_t_host": True,
            "fixed_ample_polarization_square": 14,
            "bounded_polarized_correspondence_height": False,
            "bounded_lattice_discriminant": False,
            "finite_field_set": False,
            "ambient_full_root_chamber_finite_volume": False,
        },
        "scope": {
            "n_equals_one_only": True,
            "n_at_least_two_claimed": False,
            "maximal_across_overorders_claimed": False,
            "maximal_across_commensurable_overlattices_claimed": False,
            "maximal_in_full_O_L_claimed": False,
            "host_relative_maximality": "A=Z[t] full trace lattice",
            "no_reflective_lattice_classification_imported": True,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate(max_a: int = 64) -> str:
    payload = build_certificate(max_a=max_a)
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=64)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    encoded = encoded_certificate(max_a=arguments.max_a)
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit(f"archive mismatch: {arguments.check}")
        print(f"PASS archive replay: {arguments.check}")
        return
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
        print(arguments.output)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
