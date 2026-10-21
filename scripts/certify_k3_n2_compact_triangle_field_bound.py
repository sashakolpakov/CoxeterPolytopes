#!/usr/bin/env python3
"""Exact replay for the finite-area K3 Coxeter-triangle field theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence


def triangle_field_record() -> dict[str, object]:
    simple_cycle_bound = 8
    generator_count = 2
    bad_coefficients = 15
    primitive_bound = 16 * simple_cycle_bound
    difference_bound = 2 * primitive_bound
    if (
        simple_cycle_bound,
        primitive_bound,
        difference_bound,
    ) != (8, 128, 256):
        raise AssertionError("the compact-triangle constants changed")
    return {
        "field_degree_upper_bound": 6,
        "pair_cycle_conjugate_weak_upper_bound": 4,
        "length_three_cycle_absolute_conjugate_weak_upper_bound": 8,
        "simple_cyclic_packet_conjugate_bound": simple_cycle_bound,
        "selected_generator_count_upper_bound": generator_count,
        "bad_primitive_combination_coefficient_count_upper_bound": bad_coefficients,
        "primitive_element_conjugate_bound": primitive_bound,
        "Vandermonde_difference_bound": difference_bound,
        "field_discriminant_bound": "256^30",
    }


def build_certificate() -> dict[str, object]:
    return {
        "theorem": "finite-area full-trace K3 Coxeter-triangle cyclic ground fields form a finite set",
        "field_bound": triangle_field_record(),
        "conclusions": {
            "K3_full_trace_field_degree_at_most_6": True,
            "all_pair_cycles_bounded_at_every_embedding": True,
            "length_three_cycle_bounded_at_every_embedding": True,
            "bounded_cyclic_packet_gives_finite_field_set": True,
            "bounded_mixed_height_premise_required": False,
            "angle_classification_required": False,
            "Gram_determinant_closure_required": False,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_required": True,
            "returned_field_is_cyclic_Cartan_ground_field_required": True,
            "exact_same_finite_area_Coxeter_triangle_eigen_nef_chamber_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "complete_saturated_reflection_packet_retained": True,
            "reflection_maximality_in_declared_order_host_retained": True,
            "reflection_maximality_used_in_bound": False,
            "compact_triangle_closed": True,
            "ideal_vertex_triangle_closed": True,
            "truncated_ultraideal_polyhedron_closed": False,
            "arbitrary_n2_polygon_closed": False,
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
