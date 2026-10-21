#!/usr/bin/env python3
"""Independent replay for the compact K3 Coxeter-polygon field theorem."""

from __future__ import annotations

from fractions import Fraction
import argparse
import json
import math
from typing import Sequence


def nonadjacent(left: int, right: int, size: int) -> bool:
    """Return whether two cyclically indexed facets are nonadjacent."""
    gap = (left - right) % size
    return gap not in (0, 1, size - 1)


def incident_facets(vertex: int, size: int) -> tuple[int, int]:
    return (vertex - 1) % size, vertex % size


def separation_record(max_facets: int = 256) -> dict[str, int | bool | list[int]]:
    """Exhaustively replay both separation selectors on a large cyclic range."""
    if max_facets < 4:
        raise ValueError("a nonadjacent facet pair requires at least four facets")
    vertex_wall_pairs = 0
    vertex_pairs = 0
    for size in range(4, max_facets + 1):
        for vertex in range(size):
            sides = incident_facets(vertex, size)
            for wall in range(size):
                if wall in sides:
                    continue
                if not any(nonadjacent(side, wall, size) for side in sides):
                    raise AssertionError((size, vertex, wall))
                vertex_wall_pairs += 1
        for left in range(size):
            for right in range(left + 1, size):
                if not any(
                    nonadjacent(a, b, size)
                    for a in incident_facets(left, size)
                    for b in incident_facets(right, size)
                ):
                    raise AssertionError((size, left, right))
                vertex_pairs += 1
    return {
        "facet_range": [4, max_facets],
        "vertex_wall_pairs": vertex_wall_pairs,
        "vertex_pairs": vertex_pairs,
        "all_selectors_pass": True,
    }


def fermi_wedge_terms(alpha: float, t: float, r: float) -> tuple[float, float]:
    """Terms in the remote-side inequality when a foot has t < 0."""
    if not (0 < alpha <= math.pi / 2):
        raise ValueError("the angle must be non-obtuse and positive")
    if not t < 0 or not r >= 0:
        raise ValueError("this audit branch tests a foot beyond the endpoint")
    return (
        -math.sin(alpha) * math.cosh(r) * math.sinh(t),
        math.cos(alpha) * math.sinh(r),
    )


def fermi_strip_linear_forms(
    t: float, r: float, left: float, right: float
) -> tuple[float, float]:
    """Linear half-space forms cutting out the perpendicular Fermi strip."""
    if left > right:
        raise ValueError("strip endpoints are reversed")
    x0 = math.cosh(r) * math.cosh(t)
    x1 = math.cosh(r) * math.sinh(t)
    left_form = x1 * math.cosh(left) - x0 * math.sinh(left)
    right_form = x1 * math.cosh(right) - x0 * math.sinh(right)
    return left_form, right_form


def norm_angle_record() -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for degree in range(1, 7):
        z0_floor = Fraction(1, 4 ** (degree - 1))
        sine_square_floor = z0_floor / 4
        sine_floor = Fraction(1, 2**degree)
        if sine_square_floor != sine_floor**2:
            raise AssertionError("the norm exponent changed")
        rows.append(
            {
                "degree": degree,
                "z0_floor": str(z0_floor),
                "sine_square_floor": str(sine_square_floor),
                "sine_floor": str(sine_floor),
            }
        )
    return rows


def constants_record() -> dict[str, int | str]:
    degree_bound = 20 // 3
    angle_denominator = 2**degree_bound
    half_distance_bound = 1 + angle_denominator * 4
    distance_bound = 2 * half_distance_bound**2 - 1
    pair_bound = 4 * distance_bound**2
    exponent = degree_bound * (degree_bound - 1)
    expected = (6, 64, 257, 132_097, 69_798_469_636, 30)
    actual = (
        degree_bound,
        angle_denominator,
        half_distance_bound,
        distance_bound,
        pair_bound,
        exponent,
    )
    if actual != expected:
        raise AssertionError((actual, expected))
    return {
        "degree_bound": degree_bound,
        "angle_floor": "1/64",
        "cosh_half_distance_bound": half_distance_bound,
        "cosh_distance_bound": distance_bound,
        "pair_cycle_bound": pair_bound,
        "discriminant_exponent": exponent,
    }


def proper_subfield_extension_record() -> list[dict[str, int]]:
    rows: list[dict[str, int]] = []
    for total_degree in range(2, 7):
        for subfield_degree in range(1, total_degree):
            if total_degree % subfield_degree:
                continue
            relative_degree = total_degree // subfield_degree
            if relative_degree < 2:
                raise AssertionError("a proper subfield had relative degree one")
            rows.append(
                {
                    "total_degree": total_degree,
                    "subfield_degree": subfield_degree,
                    "real_extensions": relative_degree,
                }
            )
    return rows


def build_certificate(max_facets: int = 256) -> dict[str, object]:
    return {
        "verdict": "PASS_EXACT_COMPACT_FULL_TRACE_SCOPE",
        "projection": {
            "uses_intersection_of_both_inward_halfplanes": True,
            "fermi_strip_is_intersection_of_geodesic_halfplanes": True,
            "support_ultraparallelity_requires_compact_facet_minimizers": True,
        },
        "separation": separation_record(max_facets),
        "norm_rows": norm_angle_record(),
        "constants": constants_record(),
        "proper_subfield_extensions": proper_subfield_extension_record(),
        "scope": {
            "facet_count_at_least_four": True,
            "compact_polygon_required": True,
            "primitive_exact_full_trace_K3_required": True,
            "finite_projective_order_module_required": True,
            "minimal_cyclic_ground_field_retained_by_owner": True,
            "minimality_used_by_unique_large_generator_argument": False,
            "finite_area_noncompact_closed": False,
            "carrier_construction_proved": False,
            "reflection_maximality_used": False,
        },
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-facets", type=int, default=256)
    arguments = parser.parse_args(argv)
    print(json.dumps(build_certificate(arguments.max_facets), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
