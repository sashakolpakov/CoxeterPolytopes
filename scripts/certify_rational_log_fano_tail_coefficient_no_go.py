#!/usr/bin/env python3
"""Exact replay for the hostile log-Fano tail coefficient no-go."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_log_fano_tail_coefficient_no_go.json"
)


def as_fraction(value: int | str | Fraction) -> Fraction:
    """Convert exact JSON/test input and reject floats."""
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError("coefficients must be exact rational numbers")
    try:
        return Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("coefficients must be exact rational numbers") from exc


def tail_energy(
    length: int,
    chain_coefficients: Iterable[int | str | Fraction] = (),
    endpoint_coefficient: int | str | Fraction = 0,
) -> dict[str, Any]:
    """Return the exact exceptional coordinates and square loss."""
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError("tail length must be a positive integer")
    coefficients = [as_fraction(value) for value in chain_coefficients]
    if not coefficients and length > 1:
        coefficients = [Fraction(0)] * (length - 1)
    if len(coefficients) != length - 1:
        raise ValueError("a length-m tail has exactly m-1 chain coefficients")
    endpoint = as_fraction(endpoint_coefficient)
    if any(value < 0 or value > 1 for value in coefficients + [endpoint]):
        raise ValueError("boundary coefficients must lie in [0,1]")

    if length == 1:
        coordinates = [1 + endpoint]
    else:
        coordinates = [1 + coefficients[0]]
        coordinates.extend(
            1 - coefficients[index - 1] + coefficients[index]
            for index in range(1, length - 1)
        )
        coordinates.append(1 - coefficients[-1] + endpoint)

    coordinate_sum = sum(coordinates, Fraction(0))
    energy = sum((value * value for value in coordinates), Fraction(0))
    cauchy_floor = coordinate_sum * coordinate_sum / length
    return {
        "length": length,
        "chain_coefficients": [str(value) for value in coefficients],
        "endpoint_coefficient": str(endpoint),
        "exceptional_coordinates": [str(value) for value in coordinates],
        "coordinate_sum": str(coordinate_sum),
        "expected_coordinate_sum": str(length + endpoint),
        "square_loss": str(energy),
        "cauchy_floor": str(cauchy_floor),
        "loss_at_least_length": energy >= length,
        "local_tail_intersections": [
            str(coordinates[index] - coordinates[index + 1])
            for index in range(length - 1)
        ]
        + [str(coordinates[-1])],
    }


def h1_nef_intersections(
    first: int | str | Fraction,
    second: int | str | Fraction,
    gamma_first: int | str | Fraction = 0,
    gamma_second: int | str | Fraction = 0,
) -> tuple[Fraction, Fraction]:
    """Return (K+Delta).B_i for the minimal H1 incidence."""
    a_1, a_2 = as_fraction(first), as_fraction(second)
    g_1, g_2 = as_fraction(gamma_first), as_fraction(gamma_second)
    return 2 - 4 * a_1 + 2 * a_2 + g_1, 2 + 2 * a_1 - 4 * a_2 + g_2


def h2_nef_intersections(
    boundary: int | str | Fraction,
    minus_two: int | str | Fraction,
    gamma_boundary: int | str | Fraction = 0,
    gamma_minus_two: int | str | Fraction = 0,
) -> tuple[Fraction, Fraction]:
    """Return (K+Delta).(B,R) for the minimal H2 incidence."""
    a, b = as_fraction(boundary), as_fraction(minus_two)
    g_b = as_fraction(gamma_boundary)
    g_r = as_fraction(gamma_minus_two)
    return 2 - 4 * a + 2 * b + g_b, 2 * a - 2 * b + g_r


def physical_relative_coefficients(
    length: int,
    chain_coefficients: Iterable[int | str | Fraction] = (),
    endpoint_coefficient: int | str | Fraction = 0,
) -> list[Fraction]:
    """Coefficients in K+Delta-pi^*(K0+Delta0), in the prime chain."""
    row = tail_energy(length, chain_coefficients, endpoint_coefficient)
    chain = [Fraction(value) for value in row["chain_coefficients"]]
    endpoint = Fraction(row["endpoint_coefficient"])
    return [Fraction(index + 1) + value for index, value in enumerate(chain)] + [
        Fraction(length) + endpoint
    ]


def build_certificate() -> dict[str, Any]:
    samples: list[dict[str, Any]] = []
    for length in (1, 2, 3, 8):
        samples.append(tail_energy(length))
        samples.append(
            tail_energy(length, [Fraction(1, 2)] * (length - 1), Fraction(1, 2))
        )
        samples.append(tail_energy(length, [1] * (length - 1), 1))

    grid = [Fraction(index, 4) for index in range(5)]
    h1_feasible = []
    h2_feasible = []
    for first in grid:
        for second in grid:
            h1_values = h1_nef_intersections(first, second)
            if max(h1_values) <= 0:
                h1_feasible.append([str(first), str(second)])
            h2_values = h2_nef_intersections(first, second)
            if max(h2_values) <= 0:
                h2_feasible.append([str(first), str(second)])

    payload = {
        "schema_version": 2,
        "scope": {
            "method": "surface blowup, adjunction, canonical section, exact intersections",
            "absolute_tail_bound_claimed": False,
            "conditional_bound": "m<[-(K0+Delta0)]^2 if -(K+Delta) is nef and big",
            "minimal_hostile_H1_H2_weak_log_Fano_exists": False,
            "higher_parallel_dotted_H1_H2_claimed": False,
            "classification_used": False,
        },
        "tail_energy_identity": {
            "coordinate_sum": "sum(x_j)=m+b",
            "square_identity": "L_m^2=L_0^2-sum(x_j^2)",
            "sharp_inequality": "L_m^2<=L_0^2-(m+b)^2/m<=L_0^2-m",
            "equality_choice": "all tail coefficients zero",
            "sample_rows": samples,
        },
        "hostile_coefficient_rigidity": {
            "H1_feasible_quarter_grid": h1_feasible,
            "H2_feasible_quarter_grid": h2_feasible,
            "exact_H1_conclusion": "a1=a2=1 and Gamma.(B1+B2)=0",
            "exact_H2_conclusion": "a=b=1 and Gamma.(B+R)=0",
            "positive_uniform_klt_gap_possible": False,
        },
        "relative_discrepancy_control": {
            "zero_boundary_rows": {
                str(length): [str(value) for value in physical_relative_coefficients(length)]
                for length in (1, 2, 3, 8)
            },
            "crepant_effective_tail_boundary_exists": False,
            "formal_crepant_coefficients": "a_i=-i, b=-m",
        },
        "log_Chern_changes": {
            "tail_omitted_from_reduced_boundary": {
                "log_c1_square_change": "-m",
                "log_c2_change": "+m",
            },
            "full_reduced_tail_in_boundary": {
                "log_c1_square_change": "-(m+3)",
                "log_c2_relative_to_core": "-1, independent of m",
                "scope": "m>=1; centers outside the reduced SNC core boundary",
            },
            "BMY_direction_bounds_m_without_nef_bigness": False,
        },
    }
    certify(payload)
    payload["all_checks_pass"] = True
    return payload


def certify(payload: dict[str, Any]) -> None:
    for row in payload["tail_energy_identity"]["sample_rows"]:
        if row["coordinate_sum"] != row["expected_coordinate_sum"]:
            raise AssertionError("the telescoping tail-coordinate sum failed")
        if not row["loss_at_least_length"]:
            raise AssertionError("the tail energy fell below its sharp floor")
        if Fraction(row["square_loss"]) < Fraction(row["cauchy_floor"]):
            raise AssertionError("Cauchy--Schwarz failed")

    rigidity = payload["hostile_coefficient_rigidity"]
    if rigidity["H1_feasible_quarter_grid"] != [["1", "1"]]:
        raise AssertionError("the H1 coefficient-one gate changed")
    if rigidity["H2_feasible_quarter_grid"] != [["1", "1"]]:
        raise AssertionError("the H2 coefficient-one gate changed")

    for length, coefficients in payload["relative_discrepancy_control"][
        "zero_boundary_rows"
    ].items():
        expected = [str(index) for index in range(1, int(length) + 1)]
        if coefficients != expected:
            raise AssertionError("the relative discrepancy ladder changed")


def render(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write and args.check:
        parser.error("choose at most one of --write and --check")

    payload = build_certificate()
    rendered = render(payload)
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(rendered, encoding="utf-8")
    elif args.check:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen archive is stale; rerun with --write")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
