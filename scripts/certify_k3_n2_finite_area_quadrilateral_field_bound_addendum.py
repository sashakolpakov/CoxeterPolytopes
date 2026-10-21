#!/usr/bin/env python3
"""Certificate for the ideal-vertex extension of the quadrilateral theorem."""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/k3_n2_finite_area_quadrilateral_field_bound_addendum.json"
PINS = {
    "compact_report": (
        "results/k3_n2_compact_quadrilateral_field_bound.md",
        "7f300ef5b54e9dbf0ca47140b5ca7134119046fc6f1ee39bd46bce58a6028528",
    ),
    "compact_certifier": (
        "scripts/certify_k3_n2_compact_quadrilateral_field_bound.py",
        "25ccdd5e21ae2ad5ff5d43d079fb9308b3dbacb4780fa0ba8638f7c27f9d884d",
    ),
    "compact_archive": (
        "data/enumerations/k3_n2_compact_quadrilateral_field_bound.json",
        "a61e4a28960e082380af1fc7e0dd3f8b311ab2871b6e361268f1f4e50b6213c4",
    ),
    "compact_tests": (
        "tests/test_k3_n2_compact_quadrilateral_field_bound.py",
        "0517df99c24402c881e32c6666e91abf67c0d4e3f84fb5be2cc857aaf25e072c",
    ),
    "independent_report": (
        "results/k3_n2_compact_quadrilateral_field_bound_independent_hostile_audit.md",
        "484388eeba8d7f6f73b0309b711abfc6e9fc68e4886ae003b54cc7b4477fbbbf",
    ),
    "independent_script": (
        "scripts/audit_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py",
        "79d428bbc5ef6baf09a6bb6e85c5cb0ed0946cbfad60649b59839c91b78ba32f",
    ),
    "independent_tests": (
        "tests/test_k3_n2_compact_quadrilateral_field_bound_independent_hostile.py",
        "0abeb7bb7183b92a87e8d46970d28849bdecc4732163804189ee79bf47e06eea",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        if actual != expected:
            raise ArithmeticError((label, expected, actual))
        observed[label] = actual
    return observed


def endpoint_coefficient_record() -> dict[str, object]:
    maxima = [0, 0, 0, 0]
    endpoint_count = 0
    for a, b, c, d in product((0, 1), repeat=4):
        P = a * b + c * d
        Q = a * d + b * c
        R = a * c + b * d
        K = (
            -a * a
            - b * b
            - c * c
            - d * d
            + a * a * c * c
            + b * b * d * d
            - 2 * a * b * c * d
        )
        assert 0 <= P <= 2 and 0 <= Q <= 2 and 0 <= R <= 2
        assert abs(K) <= 8
        maxima = [
            max(maxima[0], P),
            max(maxima[1], Q),
            max(maxima[2], R),
            max(maxima[3], abs(K)),
        ]
        endpoint_count += 1
    return {
        "hypercube_endpoints_checked": endpoint_count,
        "endpoint_maxima_P_Q_R_absK": maxima,
        "adjacent_magnitudes": "0<=a,b,c,d<=1",
        "weak_coefficient_bounds": "0<=P,Q,R<=2 and |K|<=8",
        "weak_determinant_growth": "(U^2-1)(V^2-1)<=4U+4V+4UV+8",
    }


def diagonal_record() -> dict[str, object]:
    trial = 32_772
    rational = 537_034_767
    surd = 8_389_888
    square_difference = rational**2 - 4097 * surd**2
    assert square_difference == 17_606_953_271_521 > 0
    return {
        "strict_norm_floor": "A-1,C-1>1/4096",
        "opposite_magnitude_strict_lower_endpoint": "sqrt(4097)/64",
        "coefficient_functions_strictly_decrease_for_V_gt_1": True,
        "weak_growth_plus_strict_forbidden_endpoint_recovers_strict_quadratic": True,
        "strict_quadratic": "U^2<(16384+256sqrt(4097))U+(32769+256sqrt(4097))",
        "trial": trial,
        "trial_remainder": f"{rational}-{surd}sqrt(4097)",
        "remainder_square_difference": square_difference,
        "both_opposite_magnitudes_strict_upper_bound": trial,
        "both_opposite_squared_ratios_strict_upper_bound": trial**2,
    }


def field_record() -> dict[str, object]:
    opposite_bound = 32_772
    H = 16 * opposite_bound**4
    difference = 32 * H
    assert H == 18_455_752_922_365_956_096
    assert difference == 590_584_093_515_710_595_072
    return {
        "field_degree_upper_bound": 6,
        "simple_cycle_weak_conjugate_height": H,
        "complete_cyclic_packet_generates_field_required": True,
        "selected_generator_count_upper_bound": 2,
        "primitive_grid": "theta1+c*theta2 for c=0,...,15",
        "Vandermonde_difference_bound": difference,
        "Vandermonde_exponent": 30,
        "field_discriminant_bound": f"{difference}^30",
        "quadratic_field_discriminant_upper_bound": 4_857_950_597,
    }


def build_certificate() -> dict[str, object]:
    return {
        "schema": "k3_n2_finite_area_quadrilateral_field_bound_addendum/v1",
        "theorem": "finite-area full-trace K3 Coxeter quadrilateral cyclic ground fields form a finite set",
        "dependency_pins": verify_pins(),
        "ideal_endpoint_coefficients": endpoint_coefficient_record(),
        "diagonal_bound": diagonal_record(),
        "field_bound": field_record(),
        "conclusions": {
            "compact_quadrilateral_theorem_extended_to_finite_or_ideal_vertices": True,
            "same_all_degree_discriminant_bound": True,
            "same_quadratic_discriminant_bound": True,
            "noncompact_infinite_area_or_ultraideal_cases_closed": False,
            "arbitrary_n2_polygon_closed": False,
        },
    }


def digest(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args(argv)
    encoded = encoded_certificate()
    if arguments.check:
        if ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("finite-area quadrilateral addendum archive mismatch")
        print("PASS finite-area quadrilateral addendum archive")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
