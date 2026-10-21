#!/usr/bin/env python3
"""Exact replay for the power-coupled rank-six K3 pentagon obstruction."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
from typing import Sequence


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def pentagon_closure_record(
    first: int | Fraction,
    second: int | Fraction,
) -> dict[str, object]:
    """Replay the five rank-three principal-minor identities exactly."""
    A = Fraction(first)
    C = Fraction(second)
    delta = A * C - A - C
    if A <= 1 or C <= 1 or delta <= 0:
        raise ValueError("compact right-pentagon parameters need A,C>1 and Delta>0")
    B = (A - 1) * (C - 1)
    D = A * (C - 1) / delta
    E = C * (A - 1) / delta
    minors = (
        C * E - C - D - E + 1,
        -A + B * E - B - E + 1,
        B * D - B - C - D + 1,
        A * D - A - D - E + 1,
        A * C - A - B - C + 1,
    )
    if any(minors):
        raise AssertionError("a rank-three pentagon principal minor did not vanish")
    if not (B > 1 and D > 1 and E > 1):
        raise AssertionError("sample is not in the compact ultraparallel region")
    return {
        "A": fraction_text(A),
        "C": fraction_text(C),
        "Delta": fraction_text(delta),
        "B": fraction_text(B),
        "D": fraction_text(D),
        "E": fraction_text(E),
        "principal_four_by_four_minors": [fraction_text(value) for value in minors],
        "all_five_principal_minors_vanish": True,
    }


def power_coupling_record(
    first: int | Fraction,
    exponent: int,
) -> dict[str, object]:
    A = Fraction(first)
    if A <= 1:
        raise ValueError("the distinguished A-coordinate must exceed one")
    if exponent < 2:
        raise ValueError("the power exponent must be at least two")
    C = A**exponent
    u = A**exponent - A ** (exponent - 1) - 1
    if u <= 0:
        raise ValueError("the chosen sample does not have positive pentagon denominator")
    closure = pentagon_closure_record(A, C)
    D = Fraction(closure["D"])
    E = Fraction(closure["E"])
    v = 4 * (E - 1)
    if E != 1 + 1 / u:
        raise AssertionError("the E closure formula failed")
    if D != 1 + A ** (exponent - 1) / u:
        raise AssertionError("the D closure formula failed")
    if u * v != 4:
        raise AssertionError("the fixed divisor identity failed")
    return {
        "A": fraction_text(A),
        "q": exponent,
        "C": fraction_text(C),
        "u": fraction_text(u),
        "D": fraction_text(D),
        "E": fraction_text(E),
        "v": fraction_text(v),
        "Delta_equals_A_times_u": Fraction(closure["Delta"]) == A * u,
        "u_times_v": fraction_text(u * v),
    }


def definite_place_record(
    conjugate: int | Fraction,
    exponent: int,
) -> dict[str, object]:
    A1 = Fraction(conjugate)
    if not 0 < A1 < 1:
        raise ValueError("the definite-place squared correlation must lie in (0,1)")
    if exponent < 2:
        raise ValueError("the power exponent must be at least two")
    u1 = A1**exponent - A1 ** (exponent - 1) - 1
    v1 = 4 / u1
    if not -2 < u1 < -1:
        raise AssertionError("the definite-place u interval failed")
    if not -4 < v1 < -2:
        raise AssertionError("the definite-place v interval failed")
    return {
        "A_conjugate": fraction_text(A1),
        "q": exponent,
        "u_conjugate": fraction_text(u1),
        "v_conjugate": fraction_text(v1),
        "u_interval": "(-2,-1)",
        "v_interval": "(-4,-2)",
    }


def proper_divisors(value: int) -> tuple[int, ...]:
    if value <= 0:
        raise ValueError("value must be positive")
    return tuple(candidate for candidate in range(1, value + 1) if value % candidate == 0)


def norm_divisor_records() -> list[dict[str, object]]:
    records = []
    for norm_magnitude in proper_divisors(16):
        records.append(
            {
                "negative_norm_magnitude": norm_magnitude,
                "divides_norm_of_four": True,
                "distinguished_u_strict_upper_bound": fraction_text(
                    Fraction(16, norm_magnitude)
                ),
                "uniform_u_strict_upper_bound": 16,
            }
        )
    return records


def squarefree_part(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    result = 1
    prime = 2
    remaining = value
    while prime * prime <= remaining:
        parity = 0
        while remaining % prime == 0:
            remaining //= prime
            parity ^= 1
        if parity:
            result *= prime
        prime += 1
    if remaining > 1:
        result *= remaining
    return result


def exact_height_inequality(trace: int, norm: int) -> bool:
    """Check (T+sqrt(D))/2 < (1+sqrt(69))/2 without floats."""
    discriminant = trace * trace - 4 * norm
    if trace < 2 or discriminant <= 0:
        return False
    r = trace - 1
    leftover = 69 - r * r - discriminant
    return leftover > 0 and leftover * leftover > 4 * r * r * discriminant


def quadratic_candidate_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for trace in range(2, 6):
        for norm in range(1, 5):
            discriminant = trace * trace - 4 * norm
            if discriminant <= 0 or isqrt(discriminant) ** 2 == discriminant:
                continue
            # A_1=(T-sqrt(D))/2 lies in (0,1).  Positivity follows from N>0.
            if discriminant <= (trace - 2) ** 2:
                continue
            if not exact_height_inequality(trace, norm):
                continue
            radicand = squarefree_part(discriminant)
            records.append(
                {
                    "trace": trace,
                    "norm": norm,
                    "polynomial": f"x^2-{trace}x+{norm}",
                    "discriminant": discriminant,
                    "field_squarefree_radicand": radicand,
                    "field": f"Q(sqrt({radicand}))",
                    "distinguished_root_above_one": True,
                    "conjugate_root_in_zero_one": True,
                    "distinguished_root_below_height_envelope": True,
                }
            )
    return records


def build_certificate() -> dict[str, object]:
    candidates = quadratic_candidate_records()
    expected = (
        (3, 1, 5, 5),
        (4, 1, 12, 3),
        (4, 2, 8, 2),
        (5, 2, 17, 17),
        (5, 3, 13, 13),
    )
    observed = tuple(
        (
            record["trace"],
            record["norm"],
            record["discriminant"],
            record["field_squarefree_radicand"],
        )
        for record in candidates
    )
    if observed != expected:
        raise AssertionError("the exact quadratic field list changed")
    return {
        "theorem": "power-coupled exact rank-six K3 right pentagons have a finite field set",
        "pentagon_closure_records": [
            pentagon_closure_record(2, 3),
            pentagon_closure_record(2, 4),
            pentagon_closure_record(Fraction(3, 2), 4),
            pentagon_closure_record(Fraction(5, 2), Fraction(11, 3)),
        ],
        "power_coupling_records": [
            power_coupling_record(2, exponent) for exponent in range(2, 9)
        ],
        "definite_place_records": [
            definite_place_record(conjugate, exponent)
            for conjugate in (Fraction(1, 10), Fraction(1, 3), Fraction(1, 2), Fraction(9, 10))
            for exponent in (2, 3, 7, 19)
        ],
        "negative_norm_divisor_records": norm_divisor_records(),
        "quadratic_field_candidates": candidates,
        "conclusions": {
            "pentagon_closure_denominator": "AC-A-C",
            "power_coupling_fixed_divisor_identity": "u*v=4",
            "negative_norm_magnitude_divides_16": True,
            "distinguished_A_upper_bound": "(1+sqrt(69))/2",
            "field_generating_pair_cycle": "xi=4A",
            "field_generating_pair_cycle_upper_bound": "2(1+sqrt(69))",
            "possible_field_radicands": [2, 3, 5, 13, 17],
            "finite_field_set_inside_power_coupled_ansatz": True,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_rank_six_required": True,
            "finite_projective_rank_three_order_module_required": True,
            "exact_same_compact_right_pentagon_eigen_nef_chamber_required": True,
            "complete_saturated_integral_reflection_packet_required": True,
            "reflection_maximality_in_declared_order_host_retained": True,
            "additional_power_coupling_C_equals_A_to_q_required": True,
            "additional_A_algebraic_integer_field_generator_required": True,
            "reflection_maximality_used_to_derive_height_bound": False,
            "arbitrary_independent_A_C_coupling_closed": False,
            "existence_of_a_carrier_in_the_ansatz_asserted": False,
            "uniform_mixed_correspondence_degree_derived": False,
            "uses_reflective_lattice_classification": False,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args(argv)
    encoded = encoded_certificate()
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
