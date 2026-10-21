"""Enumerate exact subgroup folds of the five complete high-rank K3 nef chambers.

This is an example-finding program, not a novelty certificate.  It first
computes the full automorphism group of each *weighted* complete curve graph.
For every subgroup it forms the orbit-sum Gram matrix of the fixed cone and
applies cheap exact gates before invoking Vinberg's finite-volume criterion.

The geometric reflection gate used here is deliberately strict: every
negative-square curve orbit must consist of pairwise disjoint (-2)-curves.
Then the product of their commuting root reflections restricts to reflection
in the orbit sum on the invariant space.  Graph symmetries that pass this gate
still require a separate Torelli/discriminant argument before being called
automorphisms of an algebraic K3 surface.
"""

from __future__ import annotations

from collections import Counter, deque
from itertools import combinations
import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import sympy as sp

from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import (
    _coxeter_label,
    _enumerate_vertices,
    _literal_coxiter_prop_1_10_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"

Permutation = tuple[int, ...]
Subgroup = frozenset[Permutation]


def _compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left after right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def _identity(size: int) -> Permutation:
    return tuple(range(size))


def _order(permutation: Permutation) -> int:
    identity = _identity(len(permutation))
    power = identity
    for exponent in range(1, 10_000):
        power = _compose(permutation, power)
        if power == identity:
            return exponent
    raise AssertionError("permutation order guard exceeded")


def _weighted_graph(gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    for left in range(gram.rows):
        for right in range(left):
            value = int(gram[left, right])
            if value:
                graph.add_edge(left, right, weight=value)
    return graph


def _automorphisms(gram: sp.Matrix) -> tuple[Permutation, ...]:
    graph = _weighted_graph(gram)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("weight", None),
    )
    permutations = {
        tuple(int(mapping[index]) for index in range(gram.rows))
        for mapping in matcher.isomorphisms_iter()
    }
    for permutation in permutations:
        permutation_matrix = sp.zeros(gram.rows)
        for source, target in enumerate(permutation):
            permutation_matrix[target, source] = 1
        if permutation_matrix.T * gram * permutation_matrix != gram:
            raise AssertionError("weighted graph automorphism does not preserve Gram")
    return tuple(sorted(permutations))


def _generated_subgroup(generators: Iterable[Permutation]) -> Subgroup:
    generators = tuple(generators)
    if not generators:
        raise ValueError("need a permutation size")
    subgroup = {_identity(len(generators[0])), *generators}
    queue = deque(subgroup)
    while queue:
        left = queue.popleft()
        for right in tuple(subgroup):
            for product in (_compose(left, right), _compose(right, left)):
                if product not in subgroup:
                    subgroup.add(product)
                    queue.append(product)
    return frozenset(subgroup)


def _all_subgroups(group: tuple[Permutation, ...]) -> tuple[Subgroup, ...]:
    identity_subgroup = frozenset({_identity(len(group[0]))})
    subgroups = {identity_subgroup}
    queue = deque([identity_subgroup])
    while queue:
        subgroup = queue.popleft()
        for element in group:
            if element in subgroup:
                continue
            generated = _generated_subgroup((*subgroup, element))
            if generated not in subgroups:
                subgroups.add(generated)
                queue.append(generated)
    return tuple(
        sorted(
            subgroups,
            key=lambda subgroup: (len(subgroup), tuple(sorted(subgroup))),
        )
    )


def _orbits(subgroup: Subgroup, size: int) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(size))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(sorted({element[start] for element in subgroup}))
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def _orbit_gram(
    gram: sp.Matrix, orbits: tuple[tuple[int, ...], ...]
) -> sp.Matrix:
    sums = sp.zeros(gram.rows, len(orbits))
    for column, orbit in enumerate(orbits):
        for row in orbit:
            sums[row, column] = 1
    return sums.T * gram * sums


def _pairwise_disjoint(gram: sp.Matrix, orbit: tuple[int, ...]) -> bool:
    return all(gram[left, right] == 0 for left, right in combinations(orbit, 2))


def _coxeter_edges(
    orbit_gram: sp.Matrix, facets: list[int]
) -> tuple[bool, Counter[str], list[dict[str, Any]], str | None]:
    counts: Counter[str] = Counter()
    edges: list[dict[str, Any]] = []
    try:
        for left_position, left in enumerate(facets):
            for right in facets[left_position + 1 :]:
                label = _coxeter_label(
                    sp.Rational(orbit_gram[left, left]),
                    sp.Rational(orbit_gram[right, right]),
                    sp.Rational(orbit_gram[left, right]),
                )
                counts[label] += 1
                if label != "2":
                    edges.append(
                        {
                            "facets": [left + 1, right + 1],
                            "label": label,
                        }
                    )
    except ValueError as error:
        return False, counts, edges, str(error)
    return True, counts, edges, None


def _certify_fold(
    gram: sp.Matrix,
    subgroup: Subgroup,
    orbits: tuple[tuple[int, ...], ...],
    *,
    run_vinberg: bool = True,
    full_postcheck: bool = True,
) -> dict[str, Any]:
    folded = _orbit_gram(gram, orbits)
    rank = folded.rank()
    inertia = exact_inertia(folded)
    dimension = rank - 1
    facets = [index for index in range(folded.rows) if folded[index, index] < 0]
    isotropic = [index for index in range(folded.rows) if folded[index, index] == 0]
    timelike = [index for index in range(folded.rows) if folded[index, index] > 0]
    record: dict[str, Any] = {
        "subgroup_order": len(subgroup),
        "element_orders": sorted(Counter(_order(element) for element in subgroup).elements()),
        "orbits": [[value + 1 for value in orbit] for orbit in orbits],
        "fixed_rank": rank,
        "hyperbolic_dimension": dimension,
        "folded_inertia": list(inertia),
        "facet_orbit_indices": [index + 1 for index in facets],
        "isotropic_orbit_indices": [index + 1 for index in isotropic],
        "timelike_orbit_indices": [index + 1 for index in timelike],
    }
    gates = {
        "dimension_above_nine": dimension > 9,
        "lorentzian_fixed_span": inertia == (1, rank - 1, folded.rows - rank),
        "no_isotropic_orbit_normal": not isotropic,
        "enough_facets": len(facets) >= dimension + 1,
        "facet_orbits_pairwise_disjoint": all(
            _pairwise_disjoint(gram, orbits[index]) for index in facets
        ),
    }
    record["gates"] = gates
    if not all(gates.values()):
        record["status"] = "rejected_pre_coxeter"
        return record

    is_coxeter, edge_counts, edges, error = _coxeter_edges(folded, facets)
    record["coxeter_edge_counts"] = dict(sorted(edge_counts.items()))
    record["coxeter_edges"] = edges
    record["coxeter_error"] = error
    record["gates"]["exact_coxeter_angles"] = is_coxeter
    if not is_coxeter:
        record["status"] = "rejected_non_coxeter"
        return record
    if not run_vinberg:
        record["status"] = "exact_coxeter_candidate"
        return record

    wall_gram = [
        [-int(folded[left, right]) for right in facets] for left in facets
    ]
    volume = verify_vinberg_volume_criterion(
        wall_gram, dimension, simple_roots_complete=True
    )
    record["vinberg"] = {
        "decision": volume.decision,
        "finite_volume": volume.finite_volume,
        "compact": volume.compact,
        "elliptic_edge_count": volume.elliptic_edge_count,
        "finite_vertex_count": volume.finite_vertex_count,
        "ideal_vertex_count": volume.ideal_vertex_count,
        "bad_edge_count": len(volume.finite_volume_bad_edges),
    }
    if volume.finite_volume is not True:
        record["status"] = "rejected_infinite_volume"
        return record
    if not full_postcheck:
        record["status"] = "vinberg_finite_volume_candidate"
        return record

    vertices, basis, basis_gram = _enumerate_vertices(folded, facets, timelike)
    finite = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal = [vertex for vertex in vertices if vertex["type"] == "hyperideal"]
    vertex_vectors = [sp.Matrix(vertex["basis_coordinates"]) for vertex in vertices]
    interior = sum(vertex_vectors, sp.zeros(rank, 1))
    pairing_matrix = folded.extract(range(folded.rows), basis)
    interior_pairings = pairing_matrix * interior
    strict_interior = bool(
        vertex_vectors
        and (interior.T * basis_gram * interior)[0] > 0
        and all(interior_pairings[index] > 0 for index in facets)
        and all(interior_pairings[index] > 0 for index in timelike)
    )
    facet_support: dict[str, bool] = {}
    for facet in facets:
        equality_vectors = [
            sp.Matrix(vertex["basis_coordinates"])
            for vertex in vertices
            if facet + 1 in vertex["active_facets"]
        ]
        equality_rank = (
            sp.Matrix.hstack(*equality_vectors).rank() if equality_vectors else 0
        )
        strict_exists = any(
            (pairing_matrix.extract([facet], range(rank)) * vector)[0] > 0
            for vector in vertex_vectors
        )
        facet_support[str(facet + 1)] = equality_rank == dimension and strict_exists

    literal = _literal_coxiter_prop_1_10_audit(wall_gram, dimension)
    record.update(
        {
            "folded_gram": [
                [int(folded[left, right]) for right in range(folded.cols)]
                for left in range(folded.rows)
            ],
            "vertices": vertices,
            "finite_vertex_count": len(finite),
            "ideal_vertex_count": len(ideal),
            "hyperideal_vertex_count": len(hyperideal),
            "strict_timelike_interior": strict_interior,
            "facet_support": facet_support,
            "literal_prop_1_10": {
                "spherical_edge_count": len(literal["spherical_edges"]),
                "spherical_vertex_count": len(literal["spherical_vertices"]),
                "parabolic_vertex_count": len(literal["parabolic_vertices"]),
                "bad_edges": [list(edge) for edge in literal["bad_edges"]],
                "finite_covolume": literal["finite_covolume"],
            },
        }
    )
    final_checks = {
        "no_hyperideal_vertices": not hyperideal,
        "strict_timelike_interior": strict_interior,
        "every_wall_is_a_facet": all(facet_support.values()),
        "literal_vinberg_finite_covolume": literal["finite_covolume"] is True,
    }
    record["final_checks"] = final_checks
    record["status"] = (
        "certified_finite_volume_fold"
        if all(final_checks.values())
        else "rejected_geometric_postcheck"
    )
    return record


def analyze(
    *,
    selected_rhos: set[int] | None = None,
    run_vinberg: bool = True,
    full_postcheck: bool = True,
) -> dict[str, Any]:
    source = json.loads(SOURCE.read_text())
    cases: list[dict[str, Any]] = []
    for case in source["cases"]:
        if selected_rhos is not None and int(case["rho"]) not in selected_rhos:
            continue
        gram = build_gram(case)
        group = _automorphisms(gram)
        subgroups = _all_subgroups(group)
        folds = [
            _certify_fold(
                gram,
                subgroup,
                _orbits(subgroup, gram.rows),
                run_vinberg=run_vinberg,
                full_postcheck=full_postcheck,
            )
            for subgroup in subgroups
            if len(subgroup) > 1
        ]
        cases.append(
            {
                "rho": int(case["rho"]),
                "wall_count": gram.rows,
                "sourced_general_automorphism_group": case[
                    "general_automorphism_group"
                ],
                "weighted_graph_automorphism_group_order": len(group),
                "element_order_histogram": dict(
                    sorted(Counter(_order(element) for element in group).items())
                ),
                "subgroup_count_including_identity": len(subgroups),
                "fold_status_histogram": dict(
                    sorted(Counter(fold["status"] for fold in folds).items())
                ),
                "folds": folds,
            }
        )
    return {
        "artifact_id": "k3_high_rank_nef_subgroup_fold_search",
        "artifact_status": "exact_exhaustive_weighted_graph_subgroup_search",
        "scope_warning": (
            "A weighted wall-graph subgroup is not called an algebraic K3 "
            "automorphism until a separate Torelli/discriminant realization gate passes."
        ),
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--screen-only",
        action="store_true",
        help="stop finite-volume survivors before extreme-ray/literal audits",
    )
    parser.add_argument(
        "--pre-vinberg",
        action="store_true",
        help="stop after exact Lorentzian/Coxeter/reflection gates",
    )
    parser.add_argument(
        "--rho",
        action="append",
        type=int,
        help="restrict to a Picard rank (repeatable)",
    )
    arguments = parser.parse_args()
    print(
        json.dumps(
            analyze(
                selected_rhos=set(arguments.rho) if arguments.rho else None,
                run_vinberg=not arguments.pre_vinberg,
                full_postcheck=not arguments.screen_only,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
