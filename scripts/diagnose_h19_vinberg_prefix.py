"""Verify a concise hyperideal obstruction for the bounded H19 prefix."""

from __future__ import annotations

import json

from coxeter_polytopes.search.h19_vinberg import (
    h19_basis_gram,
    run_h19_vinberg_search,
)
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray


HYPERIDEAL_WITNESS = (0, 1, 2, 3, 4, 5, 8, 9, 13)


def main() -> None:
    search = run_h19_vinberg_search(5)
    witness = verify_hyperideal_extreme_ray(
        search.roots, h19_basis_gram(), HYPERIDEAL_WITNESS, dimension=9
    )
    print(
        json.dumps(
            {
                "identity_bound": 5,
                "facet_count": len(search.roots),
                "hyperideal_witness_active_facets": list(witness.active_facets),
                "hyperideal_witness_ray": [entry.to_json() for entry in witness.ray],
                "hyperideal_witness_norm": witness.norm.to_json(),
                "conclusion": (
                    "the bounded prefix is neither compact nor finite-volume; "
                    "later roots may cut this ray"
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
