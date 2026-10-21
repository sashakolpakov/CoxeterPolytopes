#!/usr/bin/env python3
"""Archive the exact H=100 ``<-(theta-1)>+N(A1^24)`` prefix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_a1_golay import (
    COMPLETE_MIRROR_NORM_PROOF,
    UNIT_A1_GOLAY_SEED,
    unit_prefix_through_h100,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _level(level, pattern) -> dict[str, object]:
    counts: dict[str, int] = {}
    for entry in pattern.entries:
        counts[str(entry)] = counts.get(str(entry), 0) + 1
    return {
        "norm": _value(level.norm),
        "height": _value(level.height),
        "control": _value(level.control),
        "coordinate_square_target": _value(level.coordinate_square_target),
        "spatial_pattern": counts,
        "residue_profile_000_to_111": list(pattern.residue_profile),
        "labeled_wall_count": pattern.labeled_count,
    }


def _geometry(geometry) -> dict[str, object]:
    return {
        "packet_count": geometry.packet_count,
        "total_wall_count": geometry.total_wall_count,
        "extreme_ray_count": geometry.extreme_ray_count,
        "finite_vertex_count": geometry.finite_vertex_count,
        "cusp_count": geometry.cusp_count,
        "hyperideal_ray_count": geometry.hyperideal_ray_count,
        "uniform_coordinate": _value(geometry.uniform_coordinate),
        "uniform_norm": _value(geometry.uniform_norm),
        "active_wall_count": geometry.active_wall_count,
        "active_spatial_rank": geometry.active_spatial_rank,
        "finite_volume": geometry.finite_volume,
        "compact": geometry.compact,
    }


def build_audit() -> dict[str, object]:
    audit = unit_prefix_through_h100()
    return {
        "schema_version": "0.1.0",
        "lattice": "<-(theta-1)> + N(A1^24, extended Golay glue)",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "height_bound": audit.complete_height_bound,
        "canonical_mirror_norms": [
            norm.to_json() for norm in UNIT_A1_GOLAY_SEED.canonical_mirror_norms
        ],
        "complete_mirror_norm_proof": COMPLETE_MIRROR_NORM_PROOF,
        "control_level_count": audit.control_level_count,
        "raw_square_partition_count": audit.raw_square_partition_count,
        "nonempty_level_count": audit.nonempty_level_count,
        "first_packet": _level(audit.first_level, audit.first_pattern),
        "first_packet_geometry": _geometry(audit.first_geometry),
        "second_packet": _level(audit.second_level, audit.second_pattern),
        "second_packet_pairing_audits": {
            "with_first_packet": [_value(value) for value in audit.first_to_second_pairings],
            "same_height_distinct_roots": [
                _value(value) for value in audit.second_same_height_pairings
            ],
        },
        "second_packet_geometry": _geometry(audit.second_geometry),
        "scope": (
            "Every primitive mirror norm and every dominant norm equation through H=100 "
            "is included. The seven XOR Golay-weight tests reject all nonzero-residue "
            "patterns; both survivors have zero residue and explicit root-lattice "
            "labelings. Nonclosure is asserted only for this exhaustive prefix."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_a1_golay_theta_minus_one_h100.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["second_packet_geometry"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
