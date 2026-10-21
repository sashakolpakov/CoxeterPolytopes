"""Independently classify the exact wall cone of the bounded H11 search."""

from __future__ import annotations

import argparse
import json

from coxeter_polytopes.exact.golden import Golden, ZERO
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse, _mat_vec
from coxeter_polytopes.search.h11_vinberg import (
    _standard_e8_gram,
    h11_basis_gram,
    run_h11_vinberg_search,
)
from coxeter_polytopes.verify.chamber import (
    diagnose_chamber_cone,
    verify_hyperideal_extreme_ray,
)


IDENTITY_BOUND = 20
HYPERIDEAL_WITNESS = (0, 3, 4, 5, 6, 7, 8, 11, 85)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full-census",
        action="store_true",
        help="enumerate every extreme ray, rather than verify only the obstruction",
    )
    args = parser.parse_args()
    search = run_h11_vinberg_search(IDENTITY_BOUND)
    inverse_e8 = _fraction_inverse(_standard_e8_gram())
    rho = _mat_vec(inverse_e8, [Golden(1)] * 8)
    inward = (ZERO, Golden(1), *rho)
    interior = tuple(
        search.control_vector[i] + inward[i] / 10 for i in range(10)
    )
    witness = verify_hyperideal_extreme_ray(
        search.roots, h11_basis_gram(), HYPERIDEAL_WITNESS, dimension=9
    )
    result = {
        "identity_bound": IDENTITY_BOUND,
        "facet_count": len(search.roots),
        "level_solution_counts": [len(level.roots) for level in search.levels],
        "visited_node_count": sum(
            sum(level.visited_by_depth) for level in search.levels
        ),
        "status": "bounded_prefix_not_finite_volume",
        "hyperideal_witness_active_facets": list(witness.active_facets),
        "hyperideal_witness_ray": [entry.to_json() for entry in witness.ray],
        "hyperideal_witness_norm": witness.norm.to_json(),
    }
    if args.full_census:
        geometry = diagnose_chamber_cone(
            search.roots, h11_basis_gram(), interior, dimension=9
        )
        result["full_geometry"] = {
            "finite_ray_count": len(geometry.finite_vertices),
            "cusp_count": len(geometry.cusps),
            "hyperideal_ray_count": len(geometry.hyperideal_rays),
            "opposite_sheet_ray_count": len(geometry.opposite_sheet_rays),
        }
    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
