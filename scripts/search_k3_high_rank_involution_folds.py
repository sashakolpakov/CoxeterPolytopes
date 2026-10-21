"""Exhaust and certify weighted-wall-graph involution folds in ranks 15--17.

The source diagrams are the published complete (-2)-curve configurations of
the high-rank finite-automorphism K3 surfaces.  This script does two logically
separate things:

1. it enumerates *all* automorphisms of each exact weighted wall Gram matrix
   by a self-contained backtracking search, and selects every nonidentity
   involution;
2. for every conjugacy class whose fixed hyperbolic dimension is greater
   than nine, it performs exact Coxeter-angle, folded-reflection,
   irredundancy, extreme-ray, and Vinberg finite-volume checks.

A weighted graph symmetry need not be induced by an automorphism of the K3
surface.  The output records Kondō's realization result separately; in
particular, the two high-dimensional finite-volume folds found here are
graph/lattice symmetries not induced by the published automorphism actions on
the general K3 surfaces.  Possible realization at a special period is not
asserted either way.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp

from coxeter_polytopes.verify.vinberg_volume import (
    VinbergVolumeVerification,
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import (
    _coxeter_label,
    _enumerate_vertices,
    _orbit_sum_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = (
    ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"
)
ARCHIVE_PATH = (
    ROOT / "data" / "enumerations" / "k3_high_rank_involution_folds_rho15_17.json"
)
TARGET_RHOS = (15, 16, 17)

Permutation = tuple[int, ...]


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_text(value: sp.Rational) -> str:
    value = sp.Rational(value)
    return str(int(value)) if value.q == 1 else f"{value.p}/{value.q}"


def _identity(size: int) -> Permutation:
    return tuple(range(size))


def _compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left after right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def _inverse(permutation: Permutation) -> Permutation:
    inverse = [0] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(inverse)


def _conjugate(action: Permutation, element: Permutation) -> Permutation:
    return _compose(_compose(action, element), _inverse(action))


def _cycles(permutation: Permutation) -> list[list[int]]:
    cycles: list[list[int]] = []
    seen: set[int] = set()
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current + 1)
            current = permutation[current]
        cycles.append(cycle)
    return cycles


def _is_involution(permutation: Permutation) -> bool:
    return all(permutation[permutation[index]] == index for index in range(len(permutation)))


def _weighted_graph_automorphisms(
    gram: sp.Matrix,
) -> tuple[tuple[Permutation, ...], int]:
    """Enumerate the exact weighted complete-graph automorphism group.

    Initial vertex colours are full row-weight histograms.  At each search
    node the next source vertex is chosen by the minimum remaining domain,
    and a target is allowed exactly when every weight to an already mapped
    vertex agrees.  A completed bijection therefore preserves every entry of
    the Gram matrix, including omitted/zero intersections.  Conversely every
    Gram automorphism follows one branch, so the enumeration is exhaustive.
    """

    if gram.rows != gram.cols or gram != gram.T:
        raise ValueError("weighted wall Gram must be square and symmetric")
    size = gram.rows
    weights = _matrix_payload(gram)
    signatures = [tuple(sorted(Counter(row).items())) for row in weights]
    colour_classes = {
        signature: tuple(
            index for index, candidate in enumerate(signatures) if candidate == signature
        )
        for signature in set(signatures)
    }

    mapping: dict[int, int] = {}
    used_targets: set[int] = set()
    automorphisms: list[Permutation] = []
    search_node_count = 0

    def visit() -> None:
        nonlocal search_node_count
        search_node_count += 1
        if len(mapping) == size:
            permutation = tuple(mapping[index] for index in range(size))
            if any(
                weights[source][other]
                != weights[permutation[source]][permutation[other]]
                for source in range(size)
                for other in range(size)
            ):
                raise AssertionError("completed backtracking branch is not an automorphism")
            automorphisms.append(permutation)
            return

        best_source: int | None = None
        best_domain: list[int] | None = None
        for source in range(size):
            if source in mapping:
                continue
            domain = [
                target
                for target in colour_classes[signatures[source]]
                if target not in used_targets
                and all(
                    weights[source][known_source]
                    == weights[target][known_target]
                    for known_source, known_target in mapping.items()
                )
            ]
            if not domain:
                return
            if best_domain is None or len(domain) < len(best_domain):
                best_source = source
                best_domain = domain

        if best_source is None or best_domain is None:
            raise AssertionError("incomplete branch has no source vertex")
        for target in best_domain:
            mapping[best_source] = target
            used_targets.add(target)
            visit()
            used_targets.remove(target)
            del mapping[best_source]

    visit()
    result = tuple(sorted(set(automorphisms)))
    return result, search_node_count


def _group_checks(group: Sequence[Permutation]) -> dict[str, bool]:
    if not group:
        return {
            "contains_identity": False,
            "closed_under_composition": False,
            "closed_under_inverse": False,
        }
    group_set = set(group)
    return {
        "contains_identity": _identity(len(group[0])) in group_set,
        "closed_under_composition": all(
            _compose(left, right) in group_set for left in group for right in group
        ),
        "closed_under_inverse": all(_inverse(element) in group_set for element in group),
    }


def _conjugacy_classes(
    group: Sequence[Permutation], elements: Iterable[Permutation]
) -> tuple[tuple[Permutation, ...], ...]:
    remaining = set(elements)
    classes: list[tuple[Permutation, ...]] = []
    while remaining:
        representative = min(remaining)
        conjugates = {
            _conjugate(action, representative)
            for action in group
            if _conjugate(action, representative) in remaining
        }
        if not conjugates:
            raise AssertionError("conjugacy orbit unexpectedly empty")
        classes.append(tuple(sorted(conjugates)))
        remaining.difference_update(conjugates)
    return tuple(sorted(classes, key=lambda conjugacy_class: conjugacy_class[0]))


def _find_conjugator(
    group: Sequence[Permutation], source: Permutation, target: Permutation
) -> Permutation:
    for action in group:
        if _conjugate(action, source) == target:
            return action
    raise ValueError("target does not lie in the source conjugacy class")


def _orbit_relabeling(
    source_cycles: Sequence[Sequence[int]],
    target_cycles: Sequence[Sequence[int]],
    conjugator: Permutation,
) -> tuple[int, ...]:
    target_lookup = {
        frozenset(cycle): index for index, cycle in enumerate(target_cycles)
    }
    relabeling = []
    for cycle in source_cycles:
        image = frozenset(conjugator[label - 1] + 1 for label in cycle)
        relabeling.append(target_lookup[image])
    return tuple(relabeling)


def _surface_action_status(
    rho: int,
    permutation: Permutation,
    gram: sp.Matrix,
    group: Sequence[Permutation],
) -> dict[str, Any]:
    """Separate Kondō-realized actions from abstract graph symmetries."""

    cycles = _cycles(permutation)
    fixed_wall_count = sum(len(cycle) == 1 for cycle in cycles)
    internal_intersections = sorted(
        int(gram[cycle[0] - 1, cycle[1] - 1])
        for cycle in cycles
        if len(cycle) == 2
    )
    central = all(_compose(permutation, action) == _compose(action, permutation) for action in group)

    if rho == 15:
        if central:
            status = "realized_by_symplectic_surface_involution"
            reason = (
                "Kondo Section 4.5 realizes the central graph reflection by a "
                "symplectic involution."
            )
        else:
            status = "not_realized_on_general_K3_special_period_unassessed"
            reason = (
                "Kondo Section 4.5 proves that the S3 factor of the graph "
                "group is not represented symplectically; the other Z/2 "
                "factor in the automorphism group of the general surface "
                "acts trivially on NS.  A special-period realization is not "
                "decided here."
            )
    elif rho == 16:
        if 2 in internal_intersections:
            status = "realized_by_symplectic_surface_involution"
            reason = (
                "This is Kondo's tau_1: it swaps the two components E1,E2 "
                "of the square-zero pencil |E1+E2| and is realized "
                "symplectically (Section 4.4)."
            )
        else:
            status = "not_realized_on_general_K3_special_period_unassessed"
            reason = (
                "Kondo Section 4.4 excludes the second graph reflection "
                "tau_2 because it fixes fifteen rational curves; multiplying "
                "by the realized tau_1 excludes tau_1*tau_2 from the "
                "symplectic graph image.  The general nonsymplectic factor "
                "acts trivially on NS; a special-period realization is not "
                "decided here."
            )
    elif rho == 17:
        status = "realized_by_symplectic_surface_involution"
        reason = (
            "Kondo Section 4.3 realizes the unique nonidentity graph "
            "symmetry by a symplectic involution."
        )
    else:
        raise ValueError(f"unsupported rho {rho}")

    return {
        "status": status,
        "reason": reason,
        "graph_action_is_central": central,
        "fixed_wall_count": fixed_wall_count,
        "two_cycle_internal_intersections": internal_intersections,
        "source": {
            "author": "Shigeyuki Kondo",
            "title": "Algebraic K3 surfaces with finite automorphism groups",
            "doi": "10.1017/S0027763000001653",
            "locator": f"Section 4.{20-rho}",
        },
    }


def _neron_severi_action_audit(
    case: dict[str, Any], permutation: Permutation, wall_gram: sp.Matrix
) -> dict[str, Any]:
    """Recover and certify the induced integral action in a curve basis.

    The archived ``basis_indices`` are an integral Néron--Severi basis.  The
    pairings of any wall with that basis determine its coordinates exactly.
    This also lets us test whether the discriminant action is ``+1`` or
    ``-1``: in basis coordinates this is equivalent to integrality of
    ``(T-I)B^{-1}`` or ``(T+I)B^{-1}``, respectively.
    """

    basis_indices = [int(label) - 1 for label in case["basis_indices"]]
    basis_gram = wall_gram.extract(basis_indices, basis_indices)
    basis_gram_inverse = basis_gram.inv()
    wall_coordinates = [
        basis_gram_inverse * wall_gram.extract(basis_indices, [wall])
        for wall in range(wall_gram.rows)
    ]
    all_wall_coordinates_integral = all(
        sp.Rational(entry).q == 1
        for coordinates in wall_coordinates
        for entry in coordinates
    )
    action = sp.Matrix.hstack(
        *(wall_coordinates[permutation[basis_wall]] for basis_wall in basis_indices)
    )
    identity = sp.eye(basis_gram.rows)
    action_on_discriminant_is_identity = all(
        sp.Rational(entry).q == 1
        for entry in (action - identity) * basis_gram_inverse
    )
    action_on_discriminant_is_minus_identity = all(
        sp.Rational(entry).q == 1
        for entry in (action + identity) * basis_gram_inverse
    )
    permutes_every_wall = all(
        action * wall_coordinates[wall] == wall_coordinates[permutation[wall]]
        for wall in range(wall_gram.rows)
    )
    checks = {
        "archived_wall_coordinates_are_integral": all_wall_coordinates_integral,
        "action_matrix_is_integral": all(
            sp.Rational(entry).q == 1 for entry in action
        ),
        "action_is_unimodular": abs(int(action.det())) == 1,
        "action_is_an_involution": action * action == identity,
        "action_preserves_basis_gram": action.T * basis_gram * action == basis_gram,
        "action_permutes_every_archived_wall": permutes_every_wall,
    }
    return {
        "basis_wall_labels": [index + 1 for index in basis_indices],
        "basis_gram_determinant": int(basis_gram.det()),
        "action_matrix_columns": _matrix_payload(action),
        "determinant": int(action.det()),
        "fixed_rank": basis_gram.rows - (action - identity).rank(),
        "anti_fixed_rank": basis_gram.rows - (action + identity).rank(),
        "discriminant_action_is_identity": action_on_discriminant_is_identity,
        "discriminant_action_is_minus_identity": (
            action_on_discriminant_is_minus_identity
        ),
        "discriminant_action_is_scalar_plus_or_minus_one": (
            action_on_discriminant_is_identity
            or action_on_discriminant_is_minus_identity
        ),
        "checks": checks,
        "all_integral_lattice_checks_pass": all(checks.values()),
    }


def _coxeter_audit(
    orbit_gram: sp.Matrix, facet_indices: Sequence[int], cycles: Sequence[Sequence[int]]
) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    edges: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for left_position, left in enumerate(facet_indices):
        for right in facet_indices[left_position + 1 :]:
            left_norm = sp.Rational(orbit_gram[left, left])
            right_norm = sp.Rational(orbit_gram[right, right])
            product = sp.Rational(orbit_gram[left, right])
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                failures.append(
                    {
                        "orbit_indices": [left + 1, right + 1],
                        "wall_orbits": [list(cycles[left]), list(cycles[right])],
                        "error": str(error),
                    }
                )
                continue
            counts[label] += 1
            if label == "2":
                continue
            ratio = sp.factor(product**2 / (left_norm * right_norm))
            edges.append(
                {
                    "orbit_indices": [left + 1, right + 1],
                    "wall_orbits": [list(cycles[left]), list(cycles[right])],
                    "norms": [int(left_norm), int(right_norm)],
                    "product": int(product),
                    "normalized_square_ratio": _rational_text(sp.Rational(ratio)),
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


def _fold_reflection_audit(
    wall_gram: sp.Matrix,
    orbit_gram: sp.Matrix,
    facet_indices: Sequence[int],
    cycles: Sequence[Sequence[int]],
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for orbit_index in facet_indices:
        orbit = list(cycles[orbit_index])
        norm = int(orbit_gram[orbit_index, orbit_index])
        if len(orbit) == 1:
            passes = norm == -2
            kind = "fixed_minus_two_root"
            intersection = None
        elif len(orbit) == 2:
            intersection = int(wall_gram[orbit[0] - 1, orbit[1] - 1])
            passes = intersection == 0 and norm == -4
            kind = "disjoint_two_cycle"
        else:
            intersection = None
            passes = False
            kind = "unsupported_orbit"
        checks.append(
            {
                "orbit_index": orbit_index + 1,
                "wall_orbit": orbit,
                "kind": kind,
                "norm": norm,
                "internal_intersection": intersection,
                "integral_reflection_on_fixed_lattice": passes,
            }
        )
    return {
        "all_folded_reflections_are_integral": all(
            check["integral_reflection_on_fixed_lattice"] for check in checks
        ),
        "checks": checks,
    }


def _canonical_vinberg_audit_digest(result: VinbergVolumeVerification) -> str:
    payload = {
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
                "is_rank_dimension_minus_one_parabolic": (
                    audit.is_rank_dimension_minus_one_parabolic
                ),
            }
            for audit in result.parabolic_neighborhood_audits
        ],
    }
    return _json_digest(payload)


def _geometry_certificate(
    rho: int,
    cycles: list[list[int]],
    orbit_gram: sp.Matrix,
    coxeter_audit: dict[str, Any],
    reflection_audit: dict[str, Any],
    neron_severi_action: dict[str, Any],
) -> dict[str, Any]:
    rank = orbit_gram.rank()
    dimension = rank - 1
    facet_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] < 0
    ]
    isotropic_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] == 0
    ]
    timelike_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] > 0
    ]
    facet_gram = orbit_gram.extract(facet_indices, facet_indices)

    vertices, basis, basis_gram = _enumerate_vertices(
        orbit_gram, facet_indices, timelike_indices
    )
    finite_vertices = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal_vertices = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal_vertices = [
        vertex for vertex in vertices if vertex["type"] == "hyperideal"
    ]
    pairing_matrix = orbit_gram.extract(range(orbit_gram.rows), basis)
    vertex_vectors = [sp.Matrix(vertex["basis_coordinates"]) for vertex in vertices]
    same_time_component = all(
        (vertex_vectors[left].T * basis_gram * vertex_vectors[right])[0] > 0
        for left in range(len(vertex_vectors))
        for right in range(left + 1, len(vertex_vectors))
    )
    interior = sum(vertex_vectors, sp.zeros(rank, 1))
    interior_pairings = pairing_matrix * interior
    interior_square = int((interior.T * basis_gram * interior)[0])
    strict_timelike_interior = (
        interior_square > 0
        and all(interior_pairings[index] > 0 for index in facet_indices)
        and all(interior_pairings[index] > 0 for index in timelike_indices)
    )

    facet_support_checks: dict[str, dict[str, Any]] = {}
    for facet in facet_indices:
        equality_vectors = [
            sp.Matrix(vertex["basis_coordinates"])
            for vertex in vertices
            if facet + 1 in vertex["active_facets"]
        ]
        equality_rank = (
            sp.Matrix.hstack(*equality_vectors).rank() if equality_vectors else 0
        )
        strict_exists = any(
            (orbit_gram.extract([facet], basis) * vector)[0] > 0
            for vector in vertex_vectors
        )
        facet_support_checks[str(facet + 1)] = {
            "equality_ray_span_rank": equality_rank,
            "strict_ray_exists": strict_exists,
            "is_genuine_facet": equality_rank == dimension and strict_exists,
        }

    ideal_link_checks: list[dict[str, Any]] = []
    for vertex in ideal_vertices:
        active = [index - 1 for index in vertex["active_facets"]]
        active_gram = orbit_gram.extract(active, active)
        inertia = exact_inertia(active_gram)
        ideal_link_checks.append(
            {
                "active_facets": vertex["active_facets"],
                "inertia": list(inertia),
                "rank": active_gram.rank(),
                "parabolic_of_rank_dimension_minus_one": (
                    inertia[0] == 0
                    and inertia[1] == dimension - 1
                    and active_gram.rank() == dimension - 1
                ),
            }
        )

    positive_normal_gram = [
        [-int(facet_gram[left, right]) for right in range(facet_gram.cols)]
        for left in range(facet_gram.rows)
    ]
    vinberg = verify_vinberg_volume_criterion(
        positive_normal_gram,
        dimension,
        # The source curve graph is the complete wall set.  Restriction to a
        # graph-isometry fixed space makes all inequalities in an orbit
        # identical, and the facet audit above proves none of the resulting
        # negative orbit inequalities is redundant.
        simple_roots_complete=True,
    )

    facet_labels = tuple(index + 1 for index in facet_indices)
    vinberg_finite = {
        tuple(facet_labels[position] for position in vertex)
        for vertex in vinberg.finite_vertices
    }
    vinberg_ideal = {
        tuple(facet_labels[position] for position in endpoint.facets)
        for endpoint in vinberg.ideal_vertices
    }
    extreme_finite = {tuple(vertex["active_facets"]) for vertex in finite_vertices}
    extreme_ideal = {tuple(vertex["active_facets"]) for vertex in ideal_vertices}
    extension_histogram = Counter(audit.endpoint_count for audit in vinberg.edge_audits)
    endpoint_type_histogram = Counter(
        (len(audit.finite_endpoints), len(audit.ideal_endpoints))
        for audit in vinberg.edge_audits
    )

    checks = {
        "dimension_is_greater_than_nine": dimension > 9,
        "fixed_space_is_lorentzian": exact_inertia(orbit_gram)
        == (1, dimension, orbit_gram.rows - rank),
        "all_orbit_normals_are_spacelike": not isotropic_indices and not timelike_indices,
        "all_wall_pairs_are_exact_coxeter_relations": coxeter_audit[
            "all_pairs_are_exact_coxeter_relations"
        ],
        "folded_reflections_are_integral": reflection_audit[
            "all_folded_reflections_are_integral"
        ],
        "graph_action_is_integral_NS_isometry": neron_severi_action[
            "all_integral_lattice_checks_pass"
        ],
        "all_vertices_are_finite_or_ideal": not hyperideal_vertices,
        "all_vertices_lie_in_one_time_component": same_time_component,
        "strict_timelike_interior_exists": strict_timelike_interior,
        "every_orbit_wall_is_irredundant": all(
            check["is_genuine_facet"] for check in facet_support_checks.values()
        ),
        "all_ideal_links_are_parabolic": all(
            check["parabolic_of_rank_dimension_minus_one"]
            for check in ideal_link_checks
        ),
        "vinberg_finite_volume": vinberg.finite_volume is True,
        "vinberg_noncompact": vinberg.compact is False,
        "vinberg_has_no_bad_edges": not vinberg.finite_volume_bad_edges,
        "every_elliptic_edge_has_two_endpoints": set(extension_histogram) == {2},
        "vinberg_vertices_match_extreme_rays": (
            vinberg_finite == extreme_finite and vinberg_ideal == extreme_ideal
        ),
        "all_full_parabolic_neighborhoods_pass": all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in vinberg.parabolic_neighborhood_audits
        ),
    }

    return {
        "ambient_rho": rho,
        "fixed_rank": rank,
        "hyperbolic_dimension": dimension,
        "facet_count": len(facet_indices),
        "orbit_cycles": cycles,
        "orbit_gram": _matrix_payload(orbit_gram),
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "facet_gram_rank": facet_gram.rank(),
        "facet_gram_inertia": list(exact_inertia(facet_gram)),
        "facet_orbits": [cycles[index] for index in facet_indices],
        "isotropic_orbits": [cycles[index] for index in isotropic_indices],
        "timelike_orbits": [cycles[index] for index in timelike_indices],
        "coxeter": coxeter_audit,
        "fold_reflections": reflection_audit,
        "neron_severi_action": neron_severi_action,
        "vertex_basis_orbits": [cycles[index] for index in basis],
        "vertices": vertices,
        "finite_vertex_count": len(finite_vertices),
        "ideal_vertex_count": len(ideal_vertices),
        "hyperideal_vertex_count": len(hyperideal_vertices),
        "interior_square": interior_square,
        "facet_support_checks": facet_support_checks,
        "ideal_link_checks": ideal_link_checks,
        "vinberg_volume_criterion": {
            "simple_roots_complete": vinberg.simple_roots_complete,
            "completeness_evidence": (
                "published_complete_K3_wall_graph_plus_exact_fixed_orbit_"
                "restriction_and_irredundancy"
            ),
            "decision": vinberg.decision,
            "elliptic_edge_count": vinberg.elliptic_edge_count,
            "finite_vertex_count": vinberg.finite_vertex_count,
            "ideal_vertex_count": vinberg.ideal_vertex_count,
            "edge_extension_histogram": {
                str(count): multiplicity
                for count, multiplicity in sorted(extension_histogram.items())
            },
            "endpoint_type_histogram": {
                f"{finite},{ideal}": multiplicity
                for (finite, ideal), multiplicity in sorted(
                    endpoint_type_histogram.items()
                )
            },
            "finite_volume_bad_edges": [
                list(edge) for edge in vinberg.finite_volume_bad_edges
            ],
            "ideal_vertices": [
                {
                    "facets": [facet_labels[position] for position in endpoint.facets],
                    "wall_orbits": [
                        cycles[facet_indices[position]] for position in endpoint.facets
                    ],
                    "affine_components": [
                        [facet_labels[position] for position in component]
                        for component in endpoint.affine_components
                    ],
                }
                for endpoint in vinberg.ideal_vertices
            ],
            "parabolic_neighborhood_audit_count": len(
                vinberg.parabolic_neighborhood_audits
            ),
            "canonical_full_audit_sha256": _canonical_vinberg_audit_digest(vinberg),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "finite_volume_certified": all(checks.values()),
        "compact": vinberg.compact,
    }


def _fast_involution_record(
    case: dict[str, Any],
    candidate_index: int,
    permutation: Permutation,
    wall_gram: sp.Matrix,
    group: Sequence[Permutation],
) -> tuple[dict[str, Any], sp.Matrix, dict[str, Any], dict[str, Any]]:
    rho = int(case["rho"])
    cycles = _cycles(permutation)
    if any(len(cycle) > 2 for cycle in cycles):
        raise AssertionError("selected permutation is not an involution")
    orbit_sums = _orbit_sum_matrix(wall_gram.rows, cycles)
    orbit_gram = orbit_sums.T * wall_gram * orbit_sums
    rank = orbit_gram.rank()
    dimension = rank - 1
    facet_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] < 0
    ]
    isotropic_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] == 0
    ]
    timelike_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] > 0
    ]
    coxeter = _coxeter_audit(orbit_gram, facet_indices, cycles)
    reflections = _fold_reflection_audit(
        wall_gram, orbit_gram, facet_indices, cycles
    )
    neron_severi_action = _neron_severi_action_audit(
        case, permutation, wall_gram
    )
    payload = _matrix_payload(orbit_gram)
    record = {
        "candidate_id": f"rho{rho}_involution_{candidate_index}",
        "permutation_images_one_based": [value + 1 for value in permutation],
        "cycles": cycles,
        "orbit_count": len(cycles),
        "orbit_gram": payload,
        "orbit_gram_sha256": _json_digest(payload),
        "fixed_rank": rank,
        "hyperbolic_dimension": dimension,
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "negative_orbit_count": len(facet_indices),
        "isotropic_orbits": [cycles[index] for index in isotropic_indices],
        "timelike_orbits": [cycles[index] for index in timelike_indices],
        "coxeter_edge_counts": coxeter["edge_counts"],
        "all_negative_orbit_pairs_are_exact_coxeter_relations": coxeter[
            "all_pairs_are_exact_coxeter_relations"
        ],
        "all_negative_fold_reflections_are_integral": reflections[
            "all_folded_reflections_are_integral"
        ],
        "neron_severi_action": neron_severi_action,
        "dimension_greater_than_nine": dimension > 9,
        "surface_action": _surface_action_status(
            rho, permutation, wall_gram, group
        ),
    }
    return record, orbit_gram, coxeter, reflections


def _analyze_rho(case: dict[str, Any], inventory_only: bool) -> dict[str, Any]:
    rho = int(case["rho"])
    wall_gram = build_gram(case)
    group, search_node_count = _weighted_graph_automorphisms(wall_gram)
    group_set = set(group)
    identity = _identity(wall_gram.rows)
    involutions = tuple(
        permutation
        for permutation in group
        if permutation != identity and _is_involution(permutation)
    )
    conjugacy_classes = _conjugacy_classes(group, involutions)

    fast_by_permutation: dict[Permutation, dict[str, Any]] = {}
    gram_by_permutation: dict[Permutation, sp.Matrix] = {}
    coxeter_by_permutation: dict[Permutation, dict[str, Any]] = {}
    reflection_by_permutation: dict[Permutation, dict[str, Any]] = {}
    for index, permutation in enumerate(involutions, start=1):
        record, orbit_gram, coxeter, reflections = _fast_involution_record(
            case, index, permutation, wall_gram, group
        )
        fast_by_permutation[permutation] = record
        gram_by_permutation[permutation] = orbit_gram
        coxeter_by_permutation[permutation] = coxeter
        reflection_by_permutation[permutation] = reflections

    class_records: list[dict[str, Any]] = []
    geometry_certificates: list[dict[str, Any]] = []
    for class_index, conjugacy_class in enumerate(conjugacy_classes, start=1):
        representative = conjugacy_class[0]
        representative_record = fast_by_permutation[representative]
        class_id = f"rho{rho}_class_{class_index}"
        member_transports = []
        representative_cycles = _cycles(representative)
        representative_gram = gram_by_permutation[representative]
        for member in conjugacy_class:
            conjugator = _find_conjugator(group, representative, member)
            target_cycles = _cycles(member)
            relabeling = _orbit_relabeling(
                representative_cycles, target_cycles, conjugator
            )
            target_gram = gram_by_permutation[member]
            congruence = all(
                representative_gram[left, right]
                == target_gram[relabeling[left], relabeling[right]]
                for left in range(representative_gram.rows)
                for right in range(representative_gram.cols)
            )
            if not congruence:
                raise AssertionError("conjugacy transport does not preserve orbit Gram")
            fast_by_permutation[member]["conjugacy_class_id"] = class_id
            fast_by_permutation[member]["geometry_representative_candidate_id"] = (
                representative_record["candidate_id"]
            )
            member_transports.append(
                {
                    "candidate_id": fast_by_permutation[member]["candidate_id"],
                    "conjugator_images_one_based": [value + 1 for value in conjugator],
                    "representative_orbit_to_member_orbit_one_based": [
                        value + 1 for value in relabeling
                    ],
                    "orbit_gram_congruence_verified": congruence,
                }
            )

        eligible = (
            representative_record["dimension_greater_than_nine"]
            and representative_record[
                "all_negative_orbit_pairs_are_exact_coxeter_relations"
            ]
            and representative_record["all_negative_fold_reflections_are_integral"]
            and representative_record["neron_severi_action"][
                "all_integral_lattice_checks_pass"
            ]
            and not representative_record["isotropic_orbits"]
            and not representative_record["timelike_orbits"]
        )
        geometry_id: str | None = None
        if eligible and not inventory_only:
            geometry_id = f"rho{rho}_class_{class_index}_geometry"
            geometry = _geometry_certificate(
                rho,
                representative_cycles,
                representative_gram,
                coxeter_by_permutation[representative],
                reflection_by_permutation[representative],
                representative_record["neron_severi_action"],
            )
            geometry["geometry_certificate_id"] = geometry_id
            geometry["representative_candidate_id"] = representative_record[
                "candidate_id"
            ]
            geometry["conjugacy_class_id"] = class_id
            geometry["surface_automorphism_realization"] = representative_record[
                "surface_action"
            ]
            geometry["novelty_status"] = (
                "not_assessed; no claim that the abstract Coxeter polytope is new"
            )
            geometry_certificates.append(geometry)
        for member in conjugacy_class:
            fast_by_permutation[member]["geometry_certificate_id"] = geometry_id

        class_records.append(
            {
                "conjugacy_class_id": class_id,
                "size": len(conjugacy_class),
                "representative_candidate_id": representative_record["candidate_id"],
                "member_candidate_ids": [
                    fast_by_permutation[member]["candidate_id"]
                    for member in conjugacy_class
                ],
                "eligible_for_full_geometry_gate": eligible,
                "member_transports": member_transports,
            }
        )

    involution_records = [fast_by_permutation[element] for element in involutions]
    group_checks = _group_checks(group)
    if not all(group_checks.values()) or set(group) != group_set:
        raise AssertionError("weighted graph automorphisms failed group checks")
    order_histogram = Counter()
    for element in group:
        power = identity
        for order in range(1, len(group) + 1):
            power = _compose(element, power)
            if power == identity:
                order_histogram[order] += 1
                break
        else:
            raise AssertionError("element order does not divide group order")

    return {
        "rho": rho,
        "ambient_lattice": case["lattice"],
        "wall_count": wall_gram.rows,
        "wall_gram_rank": wall_gram.rank(),
        "wall_gram_inertia": list(exact_inertia(wall_gram)),
        "weighted_graph_automorphism_search": {
            "algorithm": (
                "exact_row_weight_colour_backtracking_with_minimum_remaining_domain"
            ),
            "search_node_count": search_node_count,
            "automorphism_count": len(group),
            "nonidentity_involution_count": len(involutions),
            "conjugacy_class_count": len(conjugacy_classes),
            "element_order_histogram": {
                str(order): multiplicity
                for order, multiplicity in sorted(order_histogram.items())
            },
            "group_checks": group_checks,
            "all_permutations_preserve_full_weighted_gram": all(
                all(
                    wall_gram[source, other]
                    == wall_gram[permutation[source], permutation[other]]
                    for source in range(wall_gram.rows)
                    for other in range(wall_gram.rows)
                )
                for permutation in group
            ),
        },
        "conjugacy_classes": class_records,
        "involutions": involution_records,
        "geometry_certificates": geometry_certificates,
    }


def build_certificate(
    rhos: set[int] | None = None, *, inventory_only: bool = False
) -> dict[str, Any]:
    source = _load(SOURCE_PATH)
    selected = [
        case
        for case in source["cases"]
        if int(case["rho"]) in TARGET_RHOS
        and (rhos is None or int(case["rho"]) in rhos)
    ]
    cases = [_analyze_rho(case, inventory_only) for case in selected]
    all_geometry = [
        geometry for case in cases for geometry in case["geometry_certificates"]
    ]
    inventory_payload = [
        {
            "rho": case["rho"],
            "automorphism_count": case["weighted_graph_automorphism_search"][
                "automorphism_count"
            ],
            "involutions": [
                {
                    "candidate_id": involution["candidate_id"],
                    "permutation": involution["permutation_images_one_based"],
                    "cycles": involution["cycles"],
                    "orbit_gram_sha256": involution["orbit_gram_sha256"],
                    "fixed_rank": involution["fixed_rank"],
                    "dimension": involution["hyperbolic_dimension"],
                    "surface_action_status": involution["surface_action"]["status"],
                }
                for involution in case["involutions"]
            ],
        }
        for case in cases
    ]
    full_run = rhos is None or rhos == set(TARGET_RHOS)
    certificate = {
        "schema_version": 1,
        "status_date": "2026-08-06",
        "artifact_id": "k3_high_rank_involution_folds_rho15_17",
        "artifact_status": (
            "exact_inventory_only" if inventory_only else "exact_search_and_geometry_certified"
        ),
        "source_artifact": str(SOURCE_PATH.relative_to(ROOT)),
        "source_sha256": _sha256(SOURCE_PATH),
        "scope": (
            "all nonidentity weighted-wall-graph involutions in rho=15,16,17; "
            "full exact geometry for conjugacy classes with fixed H-dimension >9"
        ),
        "method": {
            "graph_exhaustion": (
                "self-contained exact weighted-Gram automorphism backtracking"
            ),
            "fold": "orbit-sum Gram on the graph-isometry fixed subspace",
            "divisor_lattice_gate": (
                "recover the graph action in the archived integral NS curve "
                "basis; verify unimodularity, order two, isometry, and every "
                "wall image"
            ),
            "coxeter_gate": (
                "exact normalized products: 0, 1/4, 1/2, 3/4, 1, or >1"
            ),
            "integral_reflection_gate": (
                "fixed (-2)-root, or disjoint two-cycle with orbit norm -4"
            ),
            "irredundancy_gate": (
                "all projective extreme rays enumerated exactly; equality rays "
                "span each facet hyperplane and a strict ray exists"
            ),
            "finite_volume_gate": (
                "Vinberg elliptic-edge endpoint criterion with full parabolic "
                "zero-neighborhood audit and complete simple-root hard gate"
            ),
        },
        "surface_realization_warning": (
                "A symmetry of the complete weighted curve graph induces an exact "
                "rational isometry of the divisor span and hence a Coxeter fixed "
                "slice, but need not be induced by an automorphism of the selected "
                "K3 surface.  Special-period realization requires a separate "
                "Torelli/discriminant-form argument."
        ),
        "source_action_evidence": {
            "kondo_doi": "10.1017/S0027763000001653",
            "locator": "Sections 4.3--4.5, Figures 3--5",
            "rank15": (
                "Aut(weighted graph)=S3 x Z/2; only the central Z/2 graph "
                "reflection is symplectically realized."
            ),
            "rank16": (
                "Aut(weighted graph)=(Z/2)^2; tau_1 is realized, while tau_2 "
                "is excluded because it fixes fifteen rational curves."
            ),
            "rank17": "the unique graph involution is symplectically realized",
        },
        "selected_rhos": [case["rho"] for case in cases],
        "inventory_sha256": _json_digest(inventory_payload),
        "cases": cases,
        "summary": {
            "weighted_graph_automorphism_counts": {
                str(case["rho"]): case["weighted_graph_automorphism_search"][
                    "automorphism_count"
                ]
                for case in cases
            },
            "nonidentity_involution_counts": {
                str(case["rho"]): case["weighted_graph_automorphism_search"][
                    "nonidentity_involution_count"
                ]
                for case in cases
            },
            "dimension_greater_than_nine_involution_counts": {
                str(case["rho"]): sum(
                    involution["dimension_greater_than_nine"]
                    for involution in case["involutions"]
                )
                for case in cases
            },
            "fully_analyzed_geometry_conjugacy_class_count": len(all_geometry),
            "certified_finite_volume_fold_count_up_to_graph_conjugacy": sum(
                geometry["finite_volume_certified"] for geometry in all_geometry
            ),
            "certified_finite_volume_involution_count": sum(
                len(conjugacy_class["member_candidate_ids"])
                for case in cases
                for conjugacy_class in case["conjugacy_classes"]
                if conjugacy_class["eligible_for_full_geometry_gate"]
                and any(
                    geometry["conjugacy_class_id"]
                    == conjugacy_class["conjugacy_class_id"]
                    and geometry["finite_volume_certified"]
                    for geometry in case["geometry_certificates"]
                )
            ),
            "certified_dimensions_up_to_graph_conjugacy": sorted(
                geometry["hyperbolic_dimension"]
                for geometry in all_geometry
                if geometry["finite_volume_certified"]
            ),
            "all_certified_positive_folds_are_surface_automorphism_fixed": all(
                geometry["surface_automorphism_realization"]["status"]
                == "realized_by_symplectic_surface_involution"
                for geometry in all_geometry
                if geometry["finite_volume_certified"]
            ),
            "novelty_status": (
                "unresolved; these exact folded polytopes are not claimed new"
            ),
        },
        "all_checks_pass": (
            (not full_run or len(cases) == len(TARGET_RHOS))
            and all(
                all(case["weighted_graph_automorphism_search"]["group_checks"].values())
                and case["weighted_graph_automorphism_search"][
                    "all_permutations_preserve_full_weighted_gram"
                ]
                for case in cases
            )
            and (inventory_only or all(geometry["all_checks_pass"] for geometry in all_geometry))
        ),
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rho",
        type=int,
        action="append",
        choices=TARGET_RHOS,
        help="analyze only this Picard rank (repeatable)",
    )
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="skip expensive extreme-ray and Vinberg geometry certificates",
    )
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="perform a full recomputation and compare it to the frozen archive",
    )
    arguments = parser.parse_args()
    selected = set(arguments.rho) if arguments.rho else None
    if arguments.check_archive and (selected is not None or arguments.inventory_only):
        raise SystemExit("--check-archive requires the full non-inventory run")
    certificate = build_certificate(selected, inventory_only=arguments.inventory_only)
    if arguments.check_archive:
        archived = _load(ARCHIVE_PATH)
        if certificate != archived:
            raise SystemExit("recomputed certificate differs from archive")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
