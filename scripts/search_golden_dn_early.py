#!/usr/bin/env python3
"""Archive golden ``<-2phi> + D_n`` walls/rays through ``x=phi^2``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import ONE, PHI, Golden
from coxeter_polytopes.search.golden_dn import (
    analyze_golden_dn_cone,
    golden_dn_first_cap_geometry,
    run_golden_dn_search,
)


def _golden(value: Golden) -> dict[str, str]:
    return value.to_json()


def _vector(vector: tuple[Golden, ...]) -> list[dict[str, str]]:
    return [_golden(value) for value in vector]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/golden_dn_h22_h29_through_phi2.json"),
    )
    args = parser.parse_args()
    controls = (ONE, PHI, PHI * PHI)
    records = []
    for dimension in range(22, 30):
        search = run_golden_dn_search(dimension, controls)
        geometry = analyze_golden_dn_cone(search)
        first = golden_dn_first_cap_geometry(dimension)
        records.append(
            {
                "lattice": f"<-2phi> + D{dimension}",
                "hyperbolic_dimension": dimension,
                "root_norm": 2,
                "controls": [_golden(control) for control in controls],
                "levels": [
                    {
                        "control": _golden(level.control),
                        "height": _golden(level.height),
                        "target": _golden(level.target),
                        "candidate_count": level.candidate_count,
                        "accepted_root_count": len(level.roots),
                        "accepted_roots": [
                            {
                                "coset": root.spatial_coset,
                                "vector": _vector(root.vector),
                            }
                            for root in level.roots
                        ],
                    }
                    for level in search.levels
                ],
                "wall_count": len(search.walls),
                "first_cap_geometry": {
                    "extreme_ray_count": first.extreme_ray_count,
                    "finite_vertex_count": first.finite_vertex_count,
                    "cusp_count": first.cusp_count,
                    "hyperideal_ray_count": first.hyperideal_ray_count,
                },
                "through_phi_squared_geometry": {
                    "extreme_ray_count": len(geometry.extreme_rays),
                    "finite_vertex_count": geometry.finite_vertex_count,
                    "cusp_count": geometry.cusp_count,
                    "hyperideal_ray_count": geometry.hyperideal_ray_count,
                    "irredundant": geometry.irredundant,
                    "finite_volume": geometry.finite_volume,
                    "compact": geometry.compact,
                    "rays": [
                        {
                            "vector": _vector(ray.vector),
                            "active_facets": list(ray.active_facets),
                            "norm": _golden(ray.norm),
                            "kind": ray.kind,
                        }
                        for ray in geometry.extreme_rays
                    ],
                },
            }
        )
    d24_plus = run_golden_dn_search(24, controls, d24_plus=True)
    d24_plus_geometry = analyze_golden_dn_cone(d24_plus)
    payload = {
        "schema": "golden-dn-early-vinberg-0.1.0",
        "field": "Q(phi), phi^2=phi+1",
        "scope": (
            "Exact norm-two D_n/D24+ root search for controls 1,phi,phi^2 "
            "and exact cone rays. The D_n nonunimodular mirror-norm audit is "
            "not complete; later norm classes or controls may add walls."
        ),
        "records": records,
        "d24_plus_glue_check": {
            "lattice": "<-2phi> + D24+",
            "glue_cosets": "all lambda in O/2, including 1,phi,1+phi spinor residues",
            "level_candidate_counts": [level.candidate_count for level in d24_plus.levels],
            "level_accepted_counts": [len(level.roots) for level in d24_plus.levels],
            "wall_count": len(d24_plus.walls),
            "extreme_ray_count": len(d24_plus_geometry.extreme_rays),
            "finite_vertex_count": d24_plus_geometry.finite_vertex_count,
            "cusp_count": d24_plus_geometry.cusp_count,
            "hyperideal_ray_count": d24_plus_geometry.hyperideal_ray_count,
            "irredundant": d24_plus_geometry.irredundant,
            "finite_volume": d24_plus_geometry.finite_volume,
        },
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
