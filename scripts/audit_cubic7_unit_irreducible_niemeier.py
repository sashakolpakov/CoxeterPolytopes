#!/usr/bin/env python3
"""Archive the unit-coefficient D24+ and A24 first alcoves."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_a1_golay import COMPLETE_MIRROR_NORM_PROOF
from coxeter_polytopes.search.cubic7_a24_plus_unit import unit_a24_plus_first_packet
from coxeter_polytopes.search.cubic7_d24_plus_unit import unit_d24_plus_first_packet


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def build_audit() -> dict[str, object]:
    d24 = unit_d24_plus_first_packet()
    a24 = unit_a24_plus_first_packet()
    return {
        "schema_version": "0.1.0",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "timelike_coefficient": "theta-1",
        "canonical_mirror_norms": ["1", "2"],
        "complete_mirror_norm_proof": COMPLETE_MIRROR_NORM_PROOF,
        "shared_first_control": _value(d24.first_control),
        "shared_first_height": _value(d24.first_height),
        "D24_plus": {
            "root_norm": _value(d24.root_norm),
            "spatial_coset": d24.spatial_coset,
            "glue_basis_integral": d24.glue_basis_integral,
            "reflection_integral": d24.reflection_integral,
            "affine_marks_verified": d24.affine_marks_verified,
            "cap": _value(d24.cap),
            "extreme_ray_count": len(d24.rays),
            "finite_vertex_count": d24.finite_vertex_count,
            "cusp_count": d24.cusp_count,
            "hyperideal_ray_count": d24.hyperideal_ray_count,
            "finite_volume": d24.finite_volume,
            "rays": [
                {
                    "node": ray.node,
                    "mark": ray.mark,
                    "inverse_cartan_diagonal": ray.inverse_cartan_diagonal,
                    "norm": _value(ray.norm),
                    "kind": ray.kind,
                }
                for ray in d24.rays
            ],
        },
        "A24_Niemeier": {
            "root_norm": _value(a24.root_norm),
            "spatial_coset": a24.spatial_coset,
            "reflection_integral": a24.reflection_integral,
            "cap": _value(a24.cap),
            "extreme_ray_count": len(a24.rays),
            "finite_vertex_count": a24.finite_vertex_count,
            "cusp_count": a24.cusp_count,
            "hyperideal_ray_count": a24.hyperideal_ray_count,
            "finite_volume": a24.finite_volume,
            "rays": [
                {
                    "node": ray.node,
                    "inverse_cartan_diagonal": str(ray.inverse_cartan_diagonal),
                    "norm": _value(ray.norm),
                    "kind": ray.kind,
                }
                for ray in a24.rays
            ],
        },
        "scope": (
            "These are exact affine-alcove first-packet classifications. "
            "They prove neither global nonreflectivity nor persistence of an individual ray."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_unit_irreducible_niemeier_first_packet.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: payload[key] for key in ("D24_plus", "A24_Niemeier")}, indent=2)[:2000])


if __name__ == "__main__":
    main()
