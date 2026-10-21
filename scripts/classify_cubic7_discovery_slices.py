#!/usr/bin/env python3
"""Print the exact conductor-seven calibration/discovery census as JSON."""

from __future__ import annotations

import argparse
import json

from coxeter_polytopes.classification.cubic7_packets import (
    enumerate_equal_norm_compact_simplices,
    enumerate_symmetric_compact_prisms,
)
from coxeter_polytopes.search.cubic7_diagonal import enumerate_diagonal_root_pool


def _element(value):
    return value.to_json()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=10)
    parser.add_argument("--indent", type=int, default=2)
    arguments = parser.parse_args()

    triangles = enumerate_equal_norm_compact_simplices(3)
    tetrahedra = enumerate_equal_norm_compact_simplices(4)
    prisms = enumerate_symmetric_compact_prisms()
    root_pool = enumerate_diagonal_root_pool(3, arguments.height)
    record = {
        "schema": "coxeter.cubic7.discovery_slices.v1",
        "field": {
            "name": "Q(2*cos(2*pi/7))",
            "generator_polynomial": "x^3+x^2-2*x-1",
            "ring_of_integers": "Z[theta]",
            "discriminant": 49,
        },
        "equal_norm_simplex_box": {
            "labels": [2, 3, 7],
            "compact_cubic_triangles": [list(packet.labels) for packet in triangles],
            "compact_cubic_tetrahedra": [list(packet.labels) for packet in tetrahedra],
        },
        "symmetric_base_prism_box": [
            {
                "base_labels": list(packet.base_labels),
                "lateral_labels": list(packet.lateral_labels),
                "dotted_gram_entry": _element(packet.dotted_entry),
                "cosh_base_distance": _element(packet.cosh_base_distance),
                "primitive_radical_vector": [
                    _element(value) for value in packet.radical_vector
                ],
                "base_swap_root": [
                    _element(value) for value in packet.base_swap_root
                ],
                "novelty_status": "arithmetic compact double; not a new maximal group",
            }
            for packet in prisms
        ],
        "diagonal_host": {
            "lattice": "<-theta> + <1>^3",
            "height_bound": arguments.height,
            "root_norms_exhausted": [1, 2],
            "control_counts_by_norm": dict(root_pool.control_counts_by_norm),
            "spatial_alphabet_counts_by_norm": dict(
                root_pool.spatial_alphabet_counts_by_norm
            ),
            "oriented_positive_height_root_count": len(root_pool.roots),
            "root_counts_by_norm": {
                str(norm): sum(root.norm == norm for root in root_pool.roots)
                for norm in (1, 2)
            },
            "status": (
                "complete bounded mirror pool; irredundant Vinberg chamber "
                "selection and finite-volume stopping remain"
            ),
        },
    }
    print(json.dumps(record, indent=arguments.indent, sort_keys=True))


if __name__ == "__main__":
    main()
