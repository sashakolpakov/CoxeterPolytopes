"""Certify the cyclic order-three fold of the rank-15 K3 wall chamber.

This is deliberately separate from the involution search.  In addition to
the optimized Vinberg implementation it runs the literal all-principal-
subsets formulation of CoxIter Proposition 1.10 and compares both endpoint
sets with an exact extreme-ray enumeration.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from coxeter_polytopes.verify.vinberg_volume import (
    VinbergVolumeVerification,
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import (
    _coxeter_label,
    _enumerate_vertices,
    _literal_coxiter_prop_1_10_audit,
    _orbit_sum_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = (
    ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"
)
ARCHIVE_PATH = (
    ROOT / "data" / "enumerations" / "k3_rank15_order3_h12_fold.json"
)

# One generator; its inverse has the same two wall orbits and hence the same
# folded chamber.
ORDER_THREE_CYCLES = (
    (14, 15, 16),
    (17, 18, 19),
)


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _permutation(size: int, cycles: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    permutation = list(range(size))
    for cycle in cycles:
        for index, label in enumerate(cycle):
            permutation[label - 1] = cycle[(index + 1) % len(cycle)] - 1
    return tuple(permutation)


def _orbits(permutation: tuple[int, ...]) -> list[list[int]]:
    result: list[list[int]] = []
    seen: set[int] = set()
    for start in range(len(permutation)):
        if start in seen:
            continue
        orbit: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            orbit.append(current + 1)
            current = permutation[current]
        result.append(orbit)
    return result


def _is_order_exactly_three(permutation: tuple[int, ...]) -> bool:
    square = tuple(permutation[permutation[index]] for index in range(len(permutation)))
    cube = tuple(permutation[square[index]] for index in range(len(permutation)))
    return permutation != tuple(range(len(permutation))) and cube == tuple(
        range(len(permutation))
    )


def _ns_action(
    case: dict[str, Any], wall_gram: sp.Matrix, permutation: tuple[int, ...]
) -> dict[str, Any]:
    basis = [int(label) - 1 for label in case["basis_indices"]]
    basis_gram = wall_gram.extract(basis, basis)
    inverse_gram = basis_gram.inv()
    coordinates = [
        inverse_gram * wall_gram.extract(basis, [wall])
        for wall in range(wall_gram.rows)
    ]
    action = sp.Matrix.hstack(*(coordinates[permutation[index]] for index in basis))
    identity = sp.eye(len(basis))
    square = action * action
    cube = square * action

    def integral(matrix: sp.Matrix) -> bool:
        return all(sp.Rational(entry).q == 1 for entry in matrix)

    discriminant_identity = integral((action - identity) * inverse_gram)
    discriminant_minus_identity = integral((action + identity) * inverse_gram)
    checks = {
        "all_wall_coordinates_integral": all(integral(vector) for vector in coordinates),
        "action_integral": integral(action),
        "action_unimodular": abs(int(action.det())) == 1,
        "action_has_exact_order_three": cube == identity and action != identity,
        "action_isometric": action.T * basis_gram * action == basis_gram,
        "action_permutes_every_wall": all(
            action * coordinates[index] == coordinates[permutation[index]]
            for index in range(wall_gram.rows)
        ),
    }
    return {
        "basis_wall_labels": [index + 1 for index in basis],
        "basis_gram_determinant": int(basis_gram.det()),
        "action_matrix_columns": _matrix_payload(action),
        "determinant": int(action.det()),
        "characteristic_polynomial": str(sp.factor(action.charpoly().as_expr())),
        "fixed_rank": len(basis) - (action - identity).rank(),
        "moving_rank": (action - identity).rank(),
        "discriminant_action_is_identity": discriminant_identity,
        "discriminant_action_is_minus_identity": discriminant_minus_identity,
        "discriminant_action_has_exact_order_three": (
            not discriminant_identity and cube == identity
        ),
        "generic_K3_realization": {
            "status": "not_realized_on_general_K3_special_period_unassessed",
            "reason": (
                "Kondo Section 4.5 excludes the S3 factor from the symplectic "
                "action on the general surface.  Exact discriminant action is "
                "nontrivial of order three, so it cannot glue to the identity "
                "transcendental action of a symplectic order-three automorphism "
                "at a general period.  A special period with a matching "
                "non-scalar transcendental Hodge isometry is not assessed."
            ),
            "source": {
                "author": "Shigeyuki Kondo",
                "doi": "10.1017/S0027763000001653",
                "locator": "Section 4.5, Figure 5",
            },
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _coxeter_audit(orbit_gram: sp.Matrix, orbits: list[list[int]]) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    edges: list[dict[str, Any]] = []
    failures: list[str] = []
    for left in range(orbit_gram.rows):
        for right in range(left + 1, orbit_gram.rows):
            left_norm = sp.Rational(orbit_gram[left, left])
            right_norm = sp.Rational(orbit_gram[right, right])
            product = sp.Rational(orbit_gram[left, right])
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                failures.append(f"{left + 1},{right + 1}: {error}")
                continue
            counts[label] += 1
            if label != "2":
                ratio = sp.factor(product**2 / (left_norm * right_norm))
                edges.append(
                    {
                        "orbit_indices": [left + 1, right + 1],
                        "wall_orbits": [orbits[left], orbits[right]],
                        "norms": [int(left_norm), int(right_norm)],
                        "product": int(product),
                        "normalized_square_ratio": str(ratio),
                        "label": label,
                    }
                )
    for label in ("2", "3", "4", "6", "infinity", "ultraparallel"):
        counts.setdefault(label, 0)
    return {
        "all_pairs_are_exact_coxeter_relations": not failures,
        "edge_counts": dict(sorted(counts.items())),
        "nonright_edges": edges,
        "failures": failures,
    }


def _reflection_audit(
    wall_gram: sp.Matrix, orbit_gram: sp.Matrix, orbits: list[list[int]]
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for index, orbit in enumerate(orbits):
        internal = [
            int(wall_gram[left - 1, right - 1])
            for position, left in enumerate(orbit)
            for right in orbit[position + 1 :]
        ]
        norm = int(orbit_gram[index, index])
        pairwise_disjoint = all(value == 0 for value in internal)
        passes = pairwise_disjoint and norm == -2 * len(orbit)
        records.append(
            {
                "orbit_index": index + 1,
                "wall_orbit": orbit,
                "orbit_size": len(orbit),
                "internal_intersections": internal,
                "norm": norm,
                "commuting_root_product_restricts_to_integral_orbit_reflection": passes,
            }
        )
    return {
        "principle": (
            "for k pairwise-disjoint (-2)-roots in one orbit, their commuting "
            "reflection product restricts to reflection in the orbit sum of norm -2k"
        ),
        "orbits": records,
        "all_fold_reflections_integral": all(
            record[
                "commuting_root_product_restricts_to_integral_orbit_reflection"
            ]
            for record in records
        ),
    }


def _optimized_digest(result: VinbergVolumeVerification) -> str:
    return _digest(
        {
            "finite_vertices": result.finite_vertices,
            "ideal_vertices": [
                {
                    "facets": endpoint.facets,
                    "affine_components": endpoint.affine_components,
                }
                for endpoint in result.ideal_vertices
            ],
            "edge_audits": [
                {
                    "facets": audit.facets,
                    "finite_endpoints": audit.finite_endpoints,
                    "ideal_endpoints": audit.ideal_endpoints,
                }
                for audit in result.edge_audits
            ],
            "parabolic_neighborhood_audits": [
                {
                    "affine_component": audit.affine_component,
                    "zero_neighborhood": audit.zero_neighborhood,
                    "neighborhood_components": audit.neighborhood_components,
                    "passes": audit.is_rank_dimension_minus_one_parabolic,
                }
                for audit in result.parabolic_neighborhood_audits
            ],
        }
    )


def _literal_digest(literal: dict[str, Any]) -> str:
    return _digest(
        {
            "spherical_edges": literal["spherical_edges"],
            "spherical_vertices": literal["spherical_vertices"],
            "parabolic_vertices": literal["parabolic_vertices"],
            "edge_extensions": literal["edge_extensions"],
            "bad_edges": literal["bad_edges"],
        }
    )


def build_certificate() -> dict[str, Any]:
    source = _load(SOURCE_PATH)
    case = next(case for case in source["cases"] if int(case["rho"]) == 15)
    wall_gram = build_gram(case)
    permutation = _permutation(wall_gram.rows, ORDER_THREE_CYCLES)
    inverse = tuple(permutation[permutation[index]] for index in range(len(permutation)))
    orbits = _orbits(permutation)
    inverse_orbits = _orbits(inverse)
    orbit_sums = _orbit_sum_matrix(wall_gram.rows, orbits)
    orbit_gram = orbit_sums.T * wall_gram * orbit_sums
    dimension = orbit_gram.rank() - 1
    coxeter = _coxeter_audit(orbit_gram, orbits)
    reflections = _reflection_audit(wall_gram, orbit_gram, orbits)
    ns_action = _ns_action(case, wall_gram, permutation)

    vertices, basis, basis_gram = _enumerate_vertices(
        orbit_gram, list(range(orbit_gram.rows)), []
    )
    finite = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal = [vertex for vertex in vertices if vertex["type"] == "hyperideal"]
    vectors = [sp.Matrix(vertex["basis_coordinates"]) for vertex in vertices]
    pairing_matrix = orbit_gram.extract(range(orbit_gram.rows), basis)
    interior = sum(vectors, sp.zeros(orbit_gram.rank(), 1))
    interior_square = int((interior.T * basis_gram * interior)[0])
    strict_interior = interior_square > 0 and all(
        value > 0 for value in pairing_matrix * interior
    )
    same_time_component = all(
        (vectors[left].T * basis_gram * vectors[right])[0] > 0
        for left in range(len(vectors))
        for right in range(left + 1, len(vectors))
    )
    support_checks: dict[str, dict[str, Any]] = {}
    for facet in range(orbit_gram.rows):
        equality = [
            vector
            for vector, vertex in zip(vectors, vertices, strict=True)
            if facet + 1 in vertex["active_facets"]
        ]
        equality_rank = sp.Matrix.hstack(*equality).rank() if equality else 0
        strict_ray = any(
            (orbit_gram.extract([facet], basis) * vector)[0] > 0
            for vector in vectors
        )
        support_checks[str(facet + 1)] = {
            "equality_ray_span_rank": equality_rank,
            "strict_ray_exists": strict_ray,
            "is_genuine_facet": equality_rank == dimension and strict_ray,
        }

    ideal_link_checks = []
    for vertex in ideal:
        active = [label - 1 for label in vertex["active_facets"]]
        gram = orbit_gram.extract(active, active)
        inertia = exact_inertia(gram)
        ideal_link_checks.append(
            {
                "active_facets": vertex["active_facets"],
                "rank": gram.rank(),
                "inertia": list(inertia),
                "parabolic_rank_dimension_minus_one": (
                    gram.rank() == dimension - 1
                    and inertia == (0, dimension - 1, len(active) - dimension + 1)
                ),
            }
        )

    positive_normal_gram = [
        [-int(orbit_gram[left, right]) for right in range(orbit_gram.cols)]
        for left in range(orbit_gram.rows)
    ]
    optimized = verify_vinberg_volume_criterion(
        positive_normal_gram, dimension, simple_roots_complete=True
    )
    literal = _literal_coxiter_prop_1_10_audit(
        positive_normal_gram, dimension
    )
    optimized_finite = set(optimized.finite_vertices)
    optimized_ideal = {endpoint.facets for endpoint in optimized.ideal_vertices}
    literal_finite = set(literal["spherical_vertices"])
    literal_ideal = set(literal["parabolic_vertices"])
    extreme_finite = {
        tuple(label - 1 for label in vertex["active_facets"]) for vertex in finite
    }
    extreme_ideal = {
        tuple(label - 1 for label in vertex["active_facets"]) for vertex in ideal
    }
    literal_extension_histogram = {
        str(count): multiplicity
        for count, multiplicity in sorted(literal["extension_histogram"].items())
    }
    literal_endpoint_histogram = {
        f"{finite_count},{ideal_count}": multiplicity
        for (finite_count, ideal_count), multiplicity in sorted(
            literal["endpoint_type_histogram"].items()
        )
    }
    checks = {
        "permutation_has_exact_order_three": _is_order_exactly_three(permutation),
        "inverse_has_same_wall_orbits": {
            frozenset(orbit) for orbit in inverse_orbits
        }
        == {frozenset(orbit) for orbit in orbits},
        "permutation_preserves_weighted_wall_gram": all(
            wall_gram[left, right]
            == wall_gram[permutation[left], permutation[right]]
            for left in range(wall_gram.rows)
            for right in range(wall_gram.rows)
        ),
        "integral_NS_action": ns_action["all_checks_pass"],
        "fixed_space_is_H12": exact_inertia(orbit_gram)
        == (1, 12, orbit_gram.rows - 13),
        "seventeen_spacelike_orbit_normals": orbit_gram.rows == 17
        and all(orbit_gram[index, index] < 0 for index in range(orbit_gram.rows)),
        "all_exact_coxeter_relations": coxeter[
            "all_pairs_are_exact_coxeter_relations"
        ],
        "all_fold_reflections_integral": reflections[
            "all_fold_reflections_integral"
        ],
        "all_vertices_finite_or_ideal": not hyperideal,
        "strict_timelike_interior": strict_interior,
        "one_time_component": same_time_component,
        "all_walls_irredundant": all(
            record["is_genuine_facet"] for record in support_checks.values()
        ),
        "all_extreme_ideal_links_parabolic": all(
            record["parabolic_rank_dimension_minus_one"]
            for record in ideal_link_checks
        ),
        "optimized_vinberg_finite_volume": optimized.finite_volume is True,
        "optimized_vinberg_noncompact": optimized.compact is False,
        "optimized_no_bad_edges": not optimized.finite_volume_bad_edges,
        "optimized_full_parabolic_neighborhoods": all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in optimized.parabolic_neighborhood_audits
        ),
        "literal_finite_covolume": literal["finite_covolume"] is True,
        "literal_every_edge_has_two_endpoints": literal_extension_histogram
        == {"2": 623},
        "literal_no_bad_edges": not literal["bad_edges"],
        "literal_matches_optimized_endpoints": (
            literal_finite == optimized_finite and literal_ideal == optimized_ideal
        ),
        "both_diagram_audits_match_extreme_rays": (
            optimized_finite == extreme_finite and optimized_ideal == extreme_ideal
        ),
    }

    return {
        "schema_version": 1,
        "status_date": "2026-08-06",
        "artifact_id": "k3_rank15_order3_h12_fold",
        "artifact_status": "exact_literal_and_optimized_finite_volume_certified",
        "source_artifact": str(SOURCE_PATH.relative_to(ROOT)),
        "source_sha256": _sha256(SOURCE_PATH),
        "ambient_rho": 15,
        "ambient_lattice": case["lattice"],
        "permutation_images_one_based": [value + 1 for value in permutation],
        "inverse_permutation_images_one_based": [value + 1 for value in inverse],
        "wall_orbits": orbits,
        "orbit_gram": _matrix_payload(orbit_gram),
        "orbit_gram_sha256": _digest(_matrix_payload(orbit_gram)),
        "orbit_gram_rank": orbit_gram.rank(),
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "hyperbolic_dimension": dimension,
        "facet_count": orbit_gram.rows,
        "coxeter": coxeter,
        "fold_reflections": reflections,
        "neron_severi_action": ns_action,
        "vertex_basis_orbits": [orbits[index] for index in basis],
        "vertices": vertices,
        "finite_vertex_count": len(finite),
        "ideal_vertex_count": len(ideal),
        "hyperideal_vertex_count": len(hyperideal),
        "interior_square": interior_square,
        "facet_support_checks": support_checks,
        "ideal_link_checks": ideal_link_checks,
        "optimized_vinberg": {
            "decision": optimized.decision,
            "simple_roots_complete": optimized.simple_roots_complete,
            "elliptic_edge_count": optimized.elliptic_edge_count,
            "finite_vertex_count": optimized.finite_vertex_count,
            "ideal_vertex_count": optimized.ideal_vertex_count,
            "parabolic_neighborhood_audit_count": len(
                optimized.parabolic_neighborhood_audits
            ),
            "finite_volume_bad_edges": [
                list(edge) for edge in optimized.finite_volume_bad_edges
            ],
            "ideal_vertices": [
                {
                    "facets": [index + 1 for index in endpoint.facets],
                    "affine_components": [
                        [index + 1 for index in component]
                        for component in endpoint.affine_components
                    ],
                }
                for endpoint in optimized.ideal_vertices
            ],
            "canonical_full_audit_sha256": _optimized_digest(optimized),
        },
        "literal_coxiter_prop_1_10": {
            "implementation": (
                "all principal rank-11 and rank-12 spherical subsets plus "
                "all principal subsets of orders 12 through 17 tested for "
                "componentwise parabolic rank 11"
            ),
            "spherical_rank_11_edge_count": len(literal["spherical_edges"]),
            "spherical_rank_12_vertex_count": len(literal["spherical_vertices"]),
            "parabolic_rank_11_vertex_count": len(literal["parabolic_vertices"]),
            "edge_extension_histogram": literal_extension_histogram,
            "endpoint_type_histogram": literal_endpoint_histogram,
            "parabolic_vertices_one_based": [
                [index + 1 for index in vertex]
                for vertex in literal["parabolic_vertices"]
            ],
            "bad_edges": [[index + 1 for index in edge] for edge in literal["bad_edges"]],
            "finite_covolume": literal["finite_covolume"],
            "matches_optimized_endpoint_sets": (
                literal_finite == optimized_finite and literal_ideal == optimized_ideal
            ),
            "canonical_full_literal_audit_sha256": _literal_digest(literal),
        },
        "novelty_status": "not_assessed; no claim that this Coxeter polytope is new",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "finite_volume_certified": all(checks.values()),
        "compact": optimized.compact,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="recompute the full literal certificate and compare to the archive",
    )
    arguments = parser.parse_args()
    certificate = build_certificate()
    if arguments.check_archive and certificate != _load(ARCHIVE_PATH):
        raise SystemExit("recomputed certificate differs from archive")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
