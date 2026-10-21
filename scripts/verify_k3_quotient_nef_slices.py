"""Exactly verify the automorphism-fixed quotient-nef Coxeter chambers.

The input consists of published complete K3 curve graphs and published
surface-automorphism actions.  The verifier folds each wall system by the
chosen involution, recognizes every dihedral angle exactly, and enumerates
all projective extreme rays of the folded nef cone.  No floating-point
geometry is used.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache, reduce
from itertools import combinations
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import (
    build_gram,
    exact_inertia,
)
from coxeter_polytopes.verify.vinberg_volume import verify_vinberg_volume_criterion
from coxeter_polytopes.exact.diagram import is_elliptic, is_parabolic, parabolic_rank


ROOT = Path(__file__).resolve().parents[1]
HIGH_RANK_PATH = (
    ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"
)
SLICE_PATH = ROOT / "data" / "known_examples" / "k3_quotient_nef_slices.json"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _permutation_from_cycles(
    size: int, cycles: Iterable[Iterable[int]]
) -> tuple[int, ...]:
    permutation = list(range(size))
    seen: set[int] = set()
    for raw_cycle in cycles:
        cycle = [int(value) - 1 for value in raw_cycle]
        if not cycle:
            raise ValueError("empty permutation cycle")
        if any(value < 0 or value >= size for value in cycle):
            raise ValueError(f"cycle outside 1,...,{size}: {raw_cycle}")
        if seen.intersection(cycle):
            raise ValueError(f"cycles are not disjoint: {raw_cycle}")
        seen.update(cycle)
        for index, value in enumerate(cycle):
            permutation[value] = cycle[(index + 1) % len(cycle)]
    if seen != set(range(size)):
        raise ValueError("cycles do not cover every wall label")
    return tuple(permutation)


def _orbit_sum_matrix(size: int, cycles: list[list[int]]) -> sp.Matrix:
    matrix = sp.zeros(size, len(cycles))
    for column, cycle in enumerate(cycles):
        for label in cycle:
            matrix[int(label) - 1, column] = 1
    return matrix


def _nondegenerate_principal_basis(matrix: sp.Matrix) -> tuple[int, ...]:
    rank = matrix.rank()
    for indices in combinations(range(matrix.rows), rank):
        if matrix.extract(indices, indices).det() != 0:
            return indices
    raise ValueError("nondegenerate symmetric span has no principal basis")


def _positive_primitive(vector: sp.Matrix) -> tuple[int, ...]:
    rationals = [sp.Rational(value) for value in vector]
    denominator_lcm = int(sp.ilcm(*[value.q for value in rationals]))
    integers = [int(value * denominator_lcm) for value in rationals]
    common = reduce(gcd, (abs(value) for value in integers if value), 0)
    if common == 0:
        raise ValueError("zero projective vector")
    return tuple(value // common for value in integers)


def _coxeter_label(
    left_norm: sp.Rational, right_norm: sp.Rational, product: sp.Rational
) -> str:
    if left_norm >= 0 or right_norm >= 0:
        raise ValueError("Coxeter facets require spacelike normals")
    if product < 0:
        raise ValueError("folded nef walls are not acute-angled")
    if product == 0:
        return "2"
    ratio = sp.factor(product**2 / (left_norm * right_norm))
    labels = {
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "infinity",
    }
    if ratio in labels:
        return labels[ratio]
    if ratio > 1:
        return "ultraparallel"
    raise ValueError(f"non-Coxeter normalized product {ratio}")


def _literal_coxiter_prop_1_10_audit(
    wall_gram: list[list[int]], dimension: int
) -> dict[str, Any]:
    """Enumerate the extensions in CoxIter Proposition 1.10 literally.

    This deliberately does not use the optimized full-zero-neighborhood
    construction in ``verify_vinberg_volume_criterion``.  It scans every
    spherical rank-``n-1`` and rank-``n`` principal subdiagram and every
    componentwise-parabolic rank-``n-1`` principal subdiagram, then counts
    containment extensions exactly as the proposition states.
    """

    facet_count = len(wall_gram)

    def subgram(facets: tuple[int, ...]) -> list[list[int]]:
        return [[wall_gram[left][right] for right in facets] for left in facets]

    spherical_edges = tuple(
        facets
        for facets in combinations(range(facet_count), dimension - 1)
        if is_elliptic(subgram(facets))
    )
    spherical_vertices = tuple(
        facets
        for facets in combinations(range(facet_count), dimension)
        if is_elliptic(subgram(facets))
    )

    # A nonempty parabolic diagram of rank n-1 and c affine components has
    # order n-1+c, hence order at least n.  Scanning from n therefore remains
    # a literal exhaustive subset test while avoiding irrelevant small sets.
    parabolic_vertices: list[tuple[int, ...]] = []
    for order in range(dimension, facet_count + 1):
        for facets in combinations(range(facet_count), order):
            candidate = subgram(facets)
            if is_parabolic(candidate) and parabolic_rank(candidate) == dimension - 1:
                parabolic_vertices.append(facets)

    spherical_vertex_sets = tuple(
        (facets, set(facets)) for facets in spherical_vertices
    )
    parabolic_vertex_sets = tuple(
        (facets, set(facets)) for facets in parabolic_vertices
    )
    edge_extensions: list[dict[str, Any]] = []
    for edge in spherical_edges:
        edge_set = set(edge)
        finite_extensions = tuple(
            facets for facets, facets_set in spherical_vertex_sets if edge_set <= facets_set
        )
        ideal_extensions = tuple(
            facets for facets, facets_set in parabolic_vertex_sets if edge_set <= facets_set
        )
        edge_extensions.append(
            {
                "edge": edge,
                "finite_extensions": finite_extensions,
                "ideal_extensions": ideal_extensions,
                "extension_count": len(finite_extensions) + len(ideal_extensions),
            }
        )

    histogram = Counter(audit["extension_count"] for audit in edge_extensions)
    endpoint_type_histogram = Counter(
        (
            len(audit["finite_extensions"]),
            len(audit["ideal_extensions"]),
        )
        for audit in edge_extensions
    )
    has_required_vertex = bool(spherical_vertices or parabolic_vertices)
    bad_edges = tuple(
        audit["edge"] for audit in edge_extensions if audit["extension_count"] != 2
    )
    return {
        "spherical_edges": spherical_edges,
        "spherical_vertices": spherical_vertices,
        "parabolic_vertices": tuple(parabolic_vertices),
        "edge_extensions": tuple(edge_extensions),
        "extension_histogram": dict(sorted(histogram.items())),
        "endpoint_type_histogram": dict(sorted(endpoint_type_histogram.items())),
        "has_required_vertex": has_required_vertex,
        "bad_edges": bad_edges,
        "finite_covolume": has_required_vertex and not bad_edges,
    }


def _enumerate_vertices(
    orbit_gram: sp.Matrix,
    facet_indices: list[int],
    timelike_indices: list[int],
) -> tuple[list[dict[str, Any]], tuple[int, ...], sp.Matrix]:
    """Enumerate every projective extreme ray of the folded nef cone."""

    rank = orbit_gram.rank()
    dimension = rank - 1
    basis = _nondegenerate_principal_basis(orbit_gram)
    basis_gram = orbit_gram.extract(basis, basis)
    pairing_matrix = orbit_gram.extract(range(orbit_gram.rows), basis)
    vertices: dict[tuple[int, ...], dict[str, Any]] = {}

    for active_candidate in combinations(facet_indices, dimension):
        equations = orbit_gram.extract(active_candidate, basis)
        if equations.rank() != dimension:
            continue
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            raise AssertionError("codimension-d active set has wrong nullity")

        oriented: tuple[sp.Matrix, sp.Matrix] | None = None
        for vector in (nullspace[0], -nullspace[0]):
            pairings = pairing_matrix * vector
            if all(pairings[index] >= 0 for index in facet_indices) and all(
                pairings[index] > 0 for index in timelike_indices
            ):
                oriented = vector, pairings
                break
        if oriented is None:
            continue

        vector, pairings = oriented
        primitive = _positive_primitive(vector)
        primitive_vector = sp.Matrix(primitive)
        primitive_pairings = pairing_matrix * primitive_vector
        square = int((primitive_vector.T * basis_gram * primitive_vector)[0])
        active = tuple(
            index for index in facet_indices if primitive_pairings[index] == 0
        )
        if orbit_gram.extract(active, basis).rank() != dimension:
            raise AssertionError("recorded active set does not define an extreme ray")

        record = {
            "basis_coordinates": list(primitive),
            "active_facets": [index + 1 for index in active],
            "square": square,
            "type": "finite" if square > 0 else "ideal" if square == 0 else "hyperideal",
        }
        previous = vertices.get(primitive)
        if previous is not None and previous != record:
            raise AssertionError("inconsistent duplicate projective ray")
        vertices[primitive] = record

    ordered = sorted(vertices.values(), key=lambda item: item["basis_coordinates"])
    return ordered, basis, basis_gram


def analyze_case(
    slice_case: dict[str, Any], high_rank_case: dict[str, Any]
) -> dict[str, Any]:
    wall_gram = build_gram(high_rank_case)
    cycles = [[int(value) for value in cycle] for cycle in slice_case["involution_cycles"]]
    permutation = _permutation_from_cycles(wall_gram.rows, cycles)
    permutation_matrix = sp.zeros(wall_gram.rows)
    for source, target in enumerate(permutation):
        permutation_matrix[target, source] = 1

    orbit_sums = _orbit_sum_matrix(wall_gram.rows, cycles)
    orbit_gram = orbit_sums.T * wall_gram * orbit_sums
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
    # Curve classes have negative square in the surface convention.  Negating
    # gives the positive-normal Coxeter Gram convention used by Vinberg's
    # diagram criterion.  The published complete curve graph and the full
    # orbit fold provide a complete (not prefix) list of fixed-slice facets.
    vinberg_wall_gram = [
        [-int(facet_gram[left, right]) for right in range(facet_gram.cols)]
        for left in range(facet_gram.rows)
    ]
    vinberg_volume = verify_vinberg_volume_criterion(
        vinberg_wall_gram,
        orbit_gram.rank() - 1,
        simple_roots_complete=True,
    )
    literal_prop_1_10 = _literal_coxiter_prop_1_10_audit(
        vinberg_wall_gram,
        orbit_gram.rank() - 1,
    )
    fold_reflection_checks: dict[str, bool] = {}
    for orbit_index in facet_indices:
        orbit = cycles[orbit_index]
        if len(orbit) == 1:
            # A fixed (-2)-root reflection commutes with the involution and
            # restricts integrally to the invariant lattice.
            fold_reflection_checks[str(orbit_index + 1)] = True
        elif len(orbit) == 2:
            left, right = orbit[0] - 1, orbit[1] - 1
            # For two disjoint roots the commuting product s_left*s_right
            # restricts to reflection in their orbit sum (of norm -4).
            fold_reflection_checks[str(orbit_index + 1)] = (
                wall_gram[left, right] == 0
                and orbit_gram[orbit_index, orbit_index] == -4
            )
        else:
            fold_reflection_checks[str(orbit_index + 1)] = False

    edge_counts = {"3": 0, "4": 0, "6": 0, "infinity": 0, "ultraparallel": 0}
    coxeter_edges: list[dict[str, Any]] = []
    for left_position, left in enumerate(facet_indices):
        for right in facet_indices[left_position + 1 :]:
            label = _coxeter_label(
                sp.Rational(orbit_gram[left, left]),
                sp.Rational(orbit_gram[right, right]),
                sp.Rational(orbit_gram[left, right]),
            )
            if label != "2":
                edge_counts[label] += 1
                coxeter_edges.append(
                    {"facets": [left + 1, right + 1], "label": label}
                )

    vertices, basis, basis_gram = _enumerate_vertices(
        orbit_gram, facet_indices, timelike_indices
    )
    finite = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal = [vertex for vertex in vertices if vertex["type"] == "hyperideal"]
    vertex_vectors = [sp.Matrix(vertex["basis_coordinates"]) for vertex in vertices]
    same_time_component = all(
        (vertex_vectors[left].T * basis_gram * vertex_vectors[right])[0] > 0
        for left in range(len(vertex_vectors))
        for right in range(left + 1, len(vertex_vectors))
    )
    interior = sum(vertex_vectors, sp.zeros(orbit_gram.rank(), 1))
    interior_pairings = orbit_gram.extract(range(orbit_gram.rows), basis) * interior
    interior_square = (interior.T * basis_gram * interior)[0]
    strict_timelike_interior = (
        interior_square > 0
        and all(interior_pairings[index] > 0 for index in facet_indices)
        and all(interior_pairings[index] > 0 for index in timelike_indices)
    )

    # Every inequality is a genuine facet: its equality rays span a
    # hyperplane, and some other extreme ray makes the inequality strict.
    facet_support_checks: dict[str, bool] = {}
    for facet in facet_indices:
        equality_vectors = [
            sp.Matrix(vertex["basis_coordinates"])
            for vertex in vertices
            if facet + 1 in vertex["active_facets"]
        ]
        strict_exists = any(
            (orbit_gram.extract([facet], basis) * sp.Matrix(vertex["basis_coordinates"]))[0]
            > 0
            for vertex in vertices
        )
        equality_rank = (
            sp.Matrix.hstack(*equality_vectors).rank() if equality_vectors else 0
        )
        facet_support_checks[str(facet + 1)] = (
            equality_rank == orbit_gram.rank() - 1 and strict_exists
        )

    ideal_link_checks = []
    for vertex in ideal:
        active = [index - 1 for index in vertex["active_facets"]]
        active_gram = orbit_gram.extract(active, active)
        inertia = exact_inertia(active_gram)
        ideal_link_checks.append(
            {
                "active_facets": vertex["active_facets"],
                "inertia": list(inertia),
                "rank": active_gram.rank(),
                "negative_semidefinite_of_rank_dimension_minus_one": (
                    inertia[0] == 0
                    and inertia[1] == orbit_gram.rank() - 2
                    and active_gram.rank() == orbit_gram.rank() - 2
                ),
            }
        )

    expected = slice_case["expected"]
    facet_labels = tuple(index + 1 for index in facet_indices)
    vinberg_finite_vertices = {
        tuple(facet_labels[position] for position in vertex)
        for vertex in vinberg_volume.finite_vertices
    }
    vinberg_ideal_vertices = {
        tuple(facet_labels[position] for position in endpoint.facets)
        for endpoint in vinberg_volume.ideal_vertices
    }
    literal_finite_vertices = {
        tuple(facet_labels[position] for position in vertex)
        for vertex in literal_prop_1_10["spherical_vertices"]
    }
    literal_ideal_vertices = {
        tuple(facet_labels[position] for position in vertex)
        for vertex in literal_prop_1_10["parabolic_vertices"]
    }
    literal_summary = {
        "spherical_rank_n_minus_one_count": len(
            literal_prop_1_10["spherical_edges"]
        ),
        "spherical_rank_n_count": len(literal_prop_1_10["spherical_vertices"]),
        "parabolic_rank_n_minus_one_count": len(
            literal_prop_1_10["parabolic_vertices"]
        ),
        "extension_histogram": {
            str(count): multiplicity
            for count, multiplicity in literal_prop_1_10[
                "extension_histogram"
            ].items()
        },
        "endpoint_type_histogram": {
            f"{finite_count},{ideal_count}": multiplicity
            for (finite_count, ideal_count), multiplicity in literal_prop_1_10[
                "endpoint_type_histogram"
            ].items()
        },
    }
    expected_literal_summary = expected["literal_coxiter_prop_1_10"]
    extreme_finite_vertices = {
        tuple(vertex["active_facets"]) for vertex in finite
    }
    extreme_ideal_vertices = {
        tuple(vertex["active_facets"]) for vertex in ideal
    }
    recorded_edge_counts = {
        key: edge_counts[key] for key in ("3", "4", "infinity")
    }
    checks = {
        "published_wall_gram_preserved": (
            permutation_matrix.T * wall_gram * permutation_matrix == wall_gram
        ),
        "action_is_an_involution": all(
            permutation[permutation[index]] == index
            for index in range(len(permutation))
        ),
        "fixed_rank": orbit_gram.rank() == int(expected["fixed_rank"]),
        "lorentzian_fixed_space": exact_inertia(orbit_gram)
        == (
            1,
            int(expected["fixed_rank"]) - 1,
            orbit_gram.rows - int(expected["fixed_rank"]),
        ),
        "hyperbolic_dimension": orbit_gram.rank() - 1
        == int(expected["hyperbolic_dimension"]),
        "no_isotropic_orbit_normals": not isotropic_indices,
        "facet_count": len(facet_indices) == int(expected["facet_count"]),
        "timelike_orbit_count": len(timelike_indices)
        == int(expected["redundant_timelike_orbit_count"]),
        "facet_gram_inertia": list(exact_inertia(facet_gram))
        == expected["facet_gram_inertia"],
        "coxeter_edge_counts": recorded_edge_counts
        == expected["coxeter_edge_counts"],
        "no_ultraparallel_edges": edge_counts["ultraparallel"] == 0,
        "folded_reflections_are_integral": all(fold_reflection_checks.values()),
        "finite_vertex_count": len(finite) == int(expected["finite_vertex_count"]),
        "ideal_vertex_count": len(ideal) == int(expected["ideal_vertex_count"]),
        "total_vertex_count": len(vertices) == int(expected["total_vertex_count"]),
        "no_hyperideal_vertices": not hyperideal,
        "all_vertices_in_one_time_component": same_time_component,
        "strict_timelike_interior": strict_timelike_interior,
        "every_wall_is_a_facet": all(facet_support_checks.values()),
        "all_ideal_links_are_parabolic": all(
            check["negative_semidefinite_of_rank_dimension_minus_one"]
            for check in ideal_link_checks
        ),
        "vinberg_finite_volume_criterion": vinberg_volume.finite_volume is True,
        "vinberg_noncompact_criterion": vinberg_volume.compact is False,
        "vinberg_finite_vertices_match_extreme_rays": (
            vinberg_finite_vertices == extreme_finite_vertices
        ),
        "vinberg_ideal_vertices_match_extreme_rays": (
            vinberg_ideal_vertices == extreme_ideal_vertices
        ),
        "vinberg_full_parabolic_neighborhoods": all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in vinberg_volume.parabolic_neighborhood_audits
        ),
        "literal_coxiter_prop_1_10_counts": (
            literal_summary == expected_literal_summary
        ),
        "literal_coxiter_prop_1_10_finite_covolume": (
            literal_prop_1_10["finite_covolume"] is True
            and literal_prop_1_10["has_required_vertex"] is True
            and not literal_prop_1_10["bad_edges"]
        ),
        "literal_prop_1_10_matches_optimized_endpoint_sets": (
            literal_finite_vertices == vinberg_finite_vertices
            and literal_ideal_vertices == vinberg_ideal_vertices
            and len(literal_prop_1_10["spherical_edges"])
            == vinberg_volume.elliptic_edge_count
        ),
    }
    if "facet_gram_determinant" in expected:
        checks["facet_gram_determinant"] = int(facet_gram.det()) == int(
            expected["facet_gram_determinant"]
        )
    if "facet_gram_rank" in expected:
        checks["facet_gram_rank"] = facet_gram.rank() == int(
            expected["facet_gram_rank"]
        )

    return {
        "ambient_rho": int(slice_case["ambient_rho"]),
        "ambient_lattice": slice_case["ambient_lattice"],
        "fixed_rank": orbit_gram.rank(),
        "hyperbolic_dimension": orbit_gram.rank() - 1,
        "orbit_cycles": cycles,
        "orbit_gram": [list(map(int, orbit_gram.row(index))) for index in range(orbit_gram.rows)],
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "facet_orbits": [[cycles[index] for index in facet_indices][position] for position in range(len(facet_indices))],
        "redundant_timelike_orbits": [cycles[index] for index in timelike_indices],
        "facet_gram_rank": facet_gram.rank(),
        "facet_gram_inertia": list(exact_inertia(facet_gram)),
        "coxeter_edges": coxeter_edges,
        "coxeter_edge_counts": recorded_edge_counts,
        "vertex_basis_orbits": [cycles[index] for index in basis],
        "vertices": vertices,
        "finite_vertex_count": len(finite),
        "ideal_vertex_count": len(ideal),
        "ideal_link_checks": ideal_link_checks,
        "vinberg_volume_criterion": {
            "simple_roots_complete": vinberg_volume.simple_roots_complete,
            "decision": vinberg_volume.decision,
            "elliptic_edge_count": vinberg_volume.elliptic_edge_count,
            "finite_vertex_count": vinberg_volume.finite_vertex_count,
            "ideal_vertex_count": vinberg_volume.ideal_vertex_count,
            "finite_volume_bad_edges": [
                list(edge) for edge in vinberg_volume.finite_volume_bad_edges
            ],
            "ideal_vertices": [
                {
                    "facets": [facet_labels[position] for position in endpoint.facets],
                    "affine_components": [
                        [facet_labels[position] for position in component]
                        for component in endpoint.affine_components
                    ],
                }
                for endpoint in vinberg_volume.ideal_vertices
            ],
            "parabolic_neighborhood_audit_count": len(
                vinberg_volume.parabolic_neighborhood_audits
            ),
        },
        "literal_coxiter_prop_1_10": {
            **literal_summary,
            "has_required_vertex": literal_prop_1_10["has_required_vertex"],
            "parabolic_vertices": [
                [facet_labels[position] for position in vertex]
                for vertex in literal_prop_1_10["parabolic_vertices"]
            ],
            "bad_edges": [
                [facet_labels[position] for position in edge]
                for edge in literal_prop_1_10["bad_edges"]
            ],
            "finite_covolume": literal_prop_1_10["finite_covolume"],
            "matches_optimized_endpoint_sets": (
                literal_finite_vertices == vinberg_finite_vertices
                and literal_ideal_vertices == vinberg_ideal_vertices
            ),
        },
        "facet_support_checks": facet_support_checks,
        "fold_reflection_checks": fold_reflection_checks,
        "finite_volume": not hyperideal and all(
            check["negative_semidefinite_of_rank_dimension_minus_one"]
            for check in ideal_link_checks
        ),
        "compact": not ideal and not hyperideal,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


@lru_cache(maxsize=1)
def analyze() -> dict[str, Any]:
    high_rank = _load(HIGH_RANK_PATH)
    slices = _load(SLICE_PATH)
    high_rank_by_rho = {int(case["rho"]): case for case in high_rank["cases"]}
    cases = [
        analyze_case(case, high_rank_by_rho[int(case["ambient_rho"])])
        for case in slices["cases"]
    ]
    return {
        "artifact_id": slices["artifact_id"],
        "artifact_status": slices["artifact_status"],
        "cases": cases,
        "all_checks_pass": all(case["all_checks_pass"] for case in cases),
    }


def main() -> None:
    print(json.dumps(analyze(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
