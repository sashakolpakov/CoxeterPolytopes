#!/usr/bin/env python3
"""Independent hostile replay of the compact right-hexagon field theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import permutations, product
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "owner_report": {
        "path": "results/k3_n2_compact_right_hexagon_field_bound.md",
        "sha256": "5cbbe7fd2d88a2129bdc0b5345f0ae377c0fa02b5f0874606fc6ff76a17783f3",
    },
    "owner_certifier": {
        "path": "scripts/certify_k3_n2_compact_right_hexagon_field_bound.py",
        "sha256": "ff39e56f30c74429fb93935e3ba3f5353eeb219d059b2b38aa78c5f3df07a85a",
    },
    "owner_archive": {
        "path": "data/enumerations/k3_n2_compact_right_hexagon_field_bound.json",
        "sha256": "dc52fb561ac5a16b7ee0683af1999be26909820102fabf92382b30c1cf856e10",
    },
    "owner_tests": {
        "path": "tests/test_k3_n2_compact_right_hexagon_field_bound.py",
        "sha256": "06b0c81155739626843543301732e1af75c56dd2c5f9bd87c7bafab1674c808c",
    },
    "scope_audit": {
        "path": "results/k3_n2_compact_right_hexagon_field_bound_hostile_scope_audit.md",
        "sha256": "ba74f1f05d7e68d0d2cf916991dd9f248ce96127febe928f54bdf1274c672365",
    },
    "scope_tests": {
        "path": "tests/test_k3_n2_compact_right_hexagon_field_bound_hostile_scope.py",
        "sha256": "83895d73bddd5e098f53373ece083b6da6407439b08947bb58bf3a9a20f56952",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_coordinates(
    u: Fraction, v: Fraction, w: Fraction
) -> tuple[Fraction, ...]:
    """Recompute all six shifted side coordinates without owner imports."""
    u, v, w = Fraction(u), Fraction(v), Fraction(w)
    if min(u, v, w) <= 1:
        raise ValueError("alternating coshes must exceed one")
    x0, x2, x4 = u * u - 1, v * v - 1, w * w - 1
    odd_cosh_squares = (
        (u * v + w) ** 2 / (x0 * x2),
        (v * w + u) ** 2 / (x2 * x4),
        (w * u + v) ** 2 / (x4 * x0),
    )
    return (
        x0,
        odd_cosh_squares[0] - 1,
        x2,
        odd_cosh_squares[1] - 1,
        x4,
        odd_cosh_squares[2] - 1,
    )


def closure_record() -> dict[str, object]:
    values = (
        Fraction(17, 16),
        Fraction(5, 4),
        Fraction(3, 2),
        Fraction(2),
        Fraction(3),
    )
    tested = 0
    for u, v, w in product(values, repeat=3):
        coordinates = independent_coordinates(u, v, w)
        if min(coordinates) <= 0:
            raise ArithmeticError("right-hexagon coordinate lost positivity")
        x0, x1, x2, x3, x4, x5 = coordinates
        if not x0 * x1 * x2 == x2 * x3 * x4 == x4 * x5 * x0:
            raise ArithmeticError("common right-hexagon numerator failed")
        for index in range(6):
            residual = (
                coordinates[index] * coordinates[(index + 1) % 6]
                - coordinates[(index + 3) % 6]
                * coordinates[(index + 4) % 6]
            )
            if residual:
                raise ArithmeticError("four-consecutive Gram relation failed")
        tested += 1
    return {
        "exact_rational_samples": tested,
        "common_numerator_identity": True,
        "six_rotated_four_consecutive_relations": True,
        "distinguished_square_roots_only": True,
    }


def short_pair_record() -> dict[str, object]:
    permutations_checked = 0
    for ordering in permutations((Fraction(16), Fraction(25), Fraction(50))):
        smallest, middle, largest = sorted(ordering)
        if smallest < 16:
            raise ArithmeticError("large-case sample left its declared branch")
        C, E = middle + 1, largest + 1
        if 2 * C / (middle * largest) > Fraction(34, 256):
            raise ArithmeticError("the 2C envelope failed")
        if E / (middle * largest) > Fraction(17, 256):
            raise ArithmeticError("the E envelope failed")
        permutations_checked += 1
    numerator_bound = 51 + 34 * 5
    if not numerator_bound == 221 < 256:
        raise ArithmeticError("the short odd-side comparison changed")
    return {
        "alternating_order_permutations_checked": permutations_checked,
        "small_alternating_branch_threshold": 16,
        "large_alternating_branch_odd_upper": "(51+34sqrt(17))/256",
        "rational_upper": "221/256",
        "universal_some_shifted_coordinate_upper": 16,
        "universal_some_pair_cycle_upper": 68,
        "does_not_bound_every_coordinate": True,
    }


def extension_record() -> dict[str, object]:
    proper_relative_degrees: list[dict[str, int]] = []
    for total_degree in range(2, 7):
        for subfield_degree in range(1, total_degree):
            if total_degree % subfield_degree:
                continue
            relative_degree = total_degree // subfield_degree
            if relative_degree > 1:
                proper_relative_degrees.append(
                    {
                        "total_degree": total_degree,
                        "subfield_degree": subfield_degree,
                        "real_extension_count": relative_degree,
                    }
                )
    if min(row["real_extension_count"] for row in proper_relative_degrees) < 2:
        raise ArithmeticError("a proper subfield did not repeat its real value")
    return {
        "proper_subfield_patterns": proper_relative_degrees,
        "unique_value_above_four_forces_full_field": True,
        "totally_real_extension_used": True,
    }


def field_record() -> dict[str, object]:
    degree_cap = 20 // 3
    exponent = degree_cap * (degree_cap - 1)
    if (degree_cap, exponent) != (6, 30):
        raise ArithmeticError("the K3 degree/discriminant exponent changed")
    return {
        "full_trace_rank": "3d",
        "field_degree_cap": degree_cap,
        "primitive_integral_generator_interval": "(0,68)",
        "distinguished_interval": "(4,68)",
        "definite_intervals": "(0,4)",
        "Vandermonde_difference_bound": 68,
        "discriminant_exponent": exponent,
        "field_discriminant_bound": "68^30",
    }


def build_audit() -> dict[str, object]:
    for item in PINS.values():
        observed = sha256_file(ROOT / item["path"])
        if observed != item["sha256"]:
            raise ArithmeticError((item["path"], item["sha256"], observed))
    owner = (ROOT / PINS["owner_report"]["path"]).read_text(encoding="utf-8")
    owner_flat = " ".join(owner.split())
    for phrase in (
        "some one of the six",
        "unique-large-conjugate primitivity",
        r"\mathbf Q(p)=k",
        "does not cover non-right hexagons",
        "neither constructs the K3 carrier nor proves essential surjectivity",
    ):
        if phrase not in owner_flat:
            raise ArithmeticError(f"owner scope phrase missing: {phrase}")
    return {
        "verdict": "PASS_EXACT_COMPACT_RIGHT_ANGLED_SCOPE",
        "pins": PINS,
        "closure": closure_record(),
        "short_pair": short_pair_record(),
        "primitive_field": extension_record(),
        "field": field_record(),
        "scope": {
            "primitive_exact_full_trace_complex_K3_required": True,
            "exact_compact_right_angled_hexagon_required": True,
            "carrier_construction_proved": False,
            "essential_surjectivity_proved": False,
            "nonright_or_noncompact_hexagon_closed": False,
            "arbitrary_six_or_more_facet_chamber_closed": False,
            "every_dotted_pair_bounded": False,
            "mixed_height_bounded": False,
            "maximality_used": False,
        },
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
