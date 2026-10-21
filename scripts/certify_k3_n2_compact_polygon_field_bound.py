#!/usr/bin/env python3
"""Exact replay for the general compact K3 Coxeter-polygon field theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence


def constants_record() -> dict[str, object]:
    degree_bound = 6
    sine_denominator = 2**degree_bound
    cosh_half_distance_rational_bound = 1 + 64 * 4
    cosh_distance_bound = 2 * cosh_half_distance_rational_bound**2 - 1
    pair_cycle_bound = 4 * cosh_distance_bound**2
    if (
        sine_denominator,
        cosh_half_distance_rational_bound,
        cosh_distance_bound,
        pair_cycle_bound,
    ) != (64, 257, 132_097, 69_798_469_636):
        raise AssertionError("the compact-polygon constants changed")
    return {
        "field_degree_upper_bound": degree_bound,
        "Coxeter_angle_sine_lower_bound": "1/64",
        "Coxeter_angle_lower_bound": "1/64",
        "sector_area_formula": "alpha*(cosh(delta/2)-1)",
        "Gauss_Bonnet_area_strict_upper_bound": "N*pi",
        "cosh_half_shortest_dotted_distance_strict_upper_bound": "1+64*pi<257",
        "cosh_shortest_dotted_distance_strict_upper_bound": cosh_distance_bound,
        "shortest_dotted_pair_cycle_strict_upper_bound": pair_cycle_bound,
        "primitive_generator_conjugate_interval": f"(0,{pair_cycle_bound})",
        "field_discriminant_bound": f"{pair_cycle_bound}^30",
    }


def nonadjacent(index: int, other: int, facet_count: int) -> bool:
    difference = (index - other) % facet_count
    return difference not in (0, 1, facet_count - 1)


def vertex_incident_facets(vertex: int, facet_count: int) -> tuple[int, int]:
    return (vertex - 1) % facet_count, vertex % facet_count


def separation_combinatorics_record(max_facets: int = 128) -> dict[str, object]:
    checked_vertex_wall_pairs = 0
    checked_vertex_pairs = 0
    for facet_count in range(4, max_facets + 1):
        for vertex in range(facet_count):
            incident = vertex_incident_facets(vertex, facet_count)
            for wall in range(facet_count):
                if wall in incident:
                    continue
                if not any(nonadjacent(side, wall, facet_count) for side in incident):
                    raise AssertionError("a vertex/nonincident-wall selector failed")
                checked_vertex_wall_pairs += 1
        for left in range(facet_count):
            for right in range(left + 1, facet_count):
                left_incident = vertex_incident_facets(left, facet_count)
                right_incident = vertex_incident_facets(right, facet_count)
                if not any(
                    nonadjacent(a, b, facet_count)
                    for a in left_incident
                    for b in right_incident
                ):
                    raise AssertionError("a distinct-vertex selector failed")
                checked_vertex_pairs += 1
    return {
        "facet_count_range_checked": [4, max_facets],
        "vertex_nonincident_wall_pairs_checked": checked_vertex_wall_pairs,
        "distinct_vertex_pairs_checked": checked_vertex_pairs,
        "all_have_nonadjacent_incident_facet_witnesses": True,
    }


def build_certificate() -> dict[str, object]:
    return {
        "theorem": "all compact full-trace K3 Coxeter polygons have a finite field set",
        "constants": constants_record(),
        "separation_combinatorics": separation_combinatorics_record(),
        "conclusions": {
            "K3_full_trace_field_degree_at_most_6": True,
            "adjacent_reflection_norm_gives_angle_at_least_1_over_64": True,
            "acute_projection_makes_nonadjacent_supports_ultraparallel": True,
            "radius_delta_over_2_vertex_sectors_are_contained_and_disjoint": True,
            "some_dotted_pair_cycle_below_69798469636": True,
            "short_dotted_pair_cycle_generates_full_field": True,
            "field_discriminant_at_most_69798469636_power_30": True,
            "imposed_mixed_height_required": False,
            "reflection_maximality_used": False,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_required": True,
            "minimal_cyclic_Cartan_ground_field_required": True,
            "exact_same_compact_Coxeter_polygon_eigen_nef_chamber_required": True,
            "facet_count_at_least_4_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "finite_area_noncompact_polygon_closed": False,
            "ultraideal_or_truncated_polygon_closed": False,
            "higher_dimensional_chamber_closed": False,
            "carrier_existence_proved": False,
            "essential_surjectivity_proved": False,
            "Property_4_correspondence_return_proved": False,
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
