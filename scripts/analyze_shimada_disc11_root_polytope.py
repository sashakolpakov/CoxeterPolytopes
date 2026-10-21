"""Analyze the distinct geometric root walls found in the K3 traversal.

This is intentionally separate from the global completeness proof.  It asks
what polyhedron the currently accumulated effective ``(-2)`` walls bound and
certifies its Gram matrix and projective extreme rays exactly.  Compactness of
this *outer wall polyhedron* does not by itself prove that no further divisor
wall cuts it; that is the job of the global induced-chamber traversal.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import inertia, rank
from coxeter_polytopes.geometric.golden_halfspaces import inverse_golden_matrix
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.explore_shimada_disc11_global_eigen_nef import explore_global_cell
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
)


def _pair(left, gram, right):
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        Golden(0),
    )


def analyze_root_polytope(path: Path, *, traversal_states: int = 10) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    traversal = explore_global_cell(path, max_states=traversal_states)
    functionals = tuple(
        tuple(
            Golden.from_json(entry)
            for entry in record["normalized_restricted_functional"]
        )
        for record in traversal["visible_root_rows"]
    )
    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    inverse = inverse_golden_matrix(gram)

    # The archived functional is positive in the chamber.  The independent
    # chamber verifier uses outward roots satisfying (ample, root) < 0.
    roots = tuple(
        tuple(
            -sum(
                (inverse[i][j] * functional[j] for j in range(10)),
                Golden(0),
            )
            for i in range(10)
        )
        for functional in functionals
    )
    wall_gram = tuple(
        tuple(_pair(left, gram, right) for right in roots) for left in roots
    )

    relation_records = []
    nonacute_pairs = []
    all_relations_coxeter = True
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            product = wall_gram[i][j]
            if product.sign() > 0:
                nonacute_pairs.append((i, j, product))
                all_relations_coxeter = False
                continue
            relation = classify_root_pair(
                wall_gram[i][i], product, wall_gram[j][j]
            )
            if relation.kind == "finite_unrecognized":
                all_relations_coxeter = False
            relation_records.append((i, j, relation))

    geometry = None
    if not nonacute_pairs:
        geometry = diagnose_chamber_cone(
            roots,
            gram,
            tuple(AMPLE_FIELD_COEFFICIENTS),
            9,
        )

    relation_counts = Counter(
        (
            relation.kind,
            str(relation.label) if relation.label is not None else "none",
            str(relation.q),
        )
        for _, _, relation in relation_records
    )
    exact_outer_compact = bool(
        geometry is not None
        and not geometry.cusps
        and not geometry.hyperideal_rays
        and not geometry.opposite_sheet_rays
    )
    return {
        "verified": (
            traversal["verified_exact_local_reductions"]
            and len(functionals) == traversal["unique_visible_root_rays"]
            and rank(wall_gram) == 10
            and inertia(gram).positive == 9
            and inertia(gram).negative == 1
            and all_relations_coxeter
            and geometry is not None
        ),
        "source_sha256": digest,
        "traversal_states": traversal["processed_physical_states"],
        "traversal_globally_complete": traversal["global_traversal_complete"],
        "distinct_visible_root_walls": len(functionals),
        "field_gram_inertia": {
            "positive": inertia(gram).positive,
            "negative": inertia(gram).negative,
            "zero": inertia(gram).zero,
        },
        "wall_gram_rank": rank(wall_gram),
        "wall_gram_inertia": {
            "positive": inertia(wall_gram).positive,
            "negative": inertia(wall_gram).negative,
            "zero": inertia(wall_gram).zero,
        },
        "wall_norms": [wall_gram[i][i].to_json() for i in range(len(roots))],
        "wall_gram": [
            [entry.to_json() for entry in row] for row in wall_gram
        ],
        "coxeter_relation_counts": [
            {
                "kind": key[0],
                "label": key[1],
                "q": key[2],
                "count": count,
            }
            for key, count in sorted(relation_counts.items())
        ],
        "nonacute_pairs": [
            {"facets": [i, j], "product": product.to_json()}
            for i, j, product in nonacute_pairs
        ],
        "all_pairs_are_coxeter": all_relations_coxeter,
        "outer_wall_polytope": (
            {
                "finite_vertices": [list(value) for value in geometry.finite_vertices],
                "cusps": [list(value) for value in geometry.cusps],
                "hyperideal_rays": [list(value) for value in geometry.hyperideal_rays],
                "opposite_sheet_rays": [
                    list(value) for value in geometry.opposite_sheet_rays
                ],
                "compact": exact_outer_compact,
            }
            if geometry is not None
            else None
        ),
        "scope": (
            "exact_geometry_of_accumulated_root_walls_not_global_wall_completeness"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--traversal-states", type=int, default=10)
    args = parser.parse_args()
    print(
        json.dumps(
            analyze_root_polytope(
                args.disc11,
                traversal_states=args.traversal_states,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
