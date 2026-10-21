#!/usr/bin/env python3
"""Hostile, polytope-first recognition audit for one-boundary Coble carriers.

This replay deliberately does not import the surface-first Coble certifier.
It starts with archived odd-unimodular H10 wall presentations, searches for a
characteristic norm-minus-one boundary facet, reconstructs the orthogonal
Coble--Mukai root chamber and its cusps, and only then applies Kondo's finite
classification as a labelled final gate.

Failure is only an obstruction to this narrowly specified Coble carrier.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


ROOT = Path(__file__).resolve().parents[1]
D5 = ROOT / "data/enumerations/shimada_disc16_d5_h10_30_lineage.json"
R12_GRAPH = (
    ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
)
R12_59 = (
    ROOT
    / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/polytope_first_coble_h10_hostile_audit.json"
)

SCHEMA = "coxeter.polytope_first_coble_h10_hostile.v1"

# Vertex order is the order of the twelve bullets in Kondo's Figure Gamma_II.
KONDO_TYPE_II_EDGES = (
    (5, 1),
    (7, 1),
    (8, 9),
    (8, 11),
    (9, 10),
    (11, 10),
    (0, 4),
    (0, 6),
    (4, 2),
    (6, 2),
    (0, 8),
    (1, 10),
    (3, 7),
    (2, 3),
    (3, 5),
)

GRAPH_IDS = {
    ("R12_UE8A1x2", 2): "R12-E8A1x2-C2-H10-12",
    ("R12_UD8A1x2", 2): "R12-D8A1x2-C2a-H10-15",
    ("R12_UD8A1x2", 3): "R12-D8A1x2-C2b-H10-16",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix was expected")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _vector_payload(vector: sp.Matrix) -> list[int]:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("an integral vector was expected")
    return [int(value) for value in vector]


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("wall coordinates left the integral lattice")
    entries = [int(value) for value in vector]
    content = math.gcd(*(abs(value) for value in entries))
    if content == 0:
        raise ArithmeticError("zero wall")
    # Divide by positive content: the archived cone orientation is retained.
    return sp.Matrix([value // content for value in entries])


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i]), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i)
                    if work[i, j]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            i, j = pair
            transform = sp.eye(work.rows)
            transform[:, i] = sp.eye(work.rows)[:, i] + sp.eye(work.rows)[:, j]
            work = transform.T * work * transform
            pivot = i
        if pivot:
            order = [pivot, *range(pivot), *range(pivot + 1, work.rows)]
            work = work.extract(order, order)
        value = work[0, 0]
        if value > 0:
            positive += 1
        elif value < 0:
            negative += 1
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return positive, negative, 0


def _characteristic(gram: sp.Matrix, vector: sp.Matrix) -> bool:
    diagonal = sp.Matrix([gram[i, i] for i in range(gram.rows)])
    return all(int(value) % 2 == 0 for value in gram * vector - diagonal)


def _kernel_basis(row: sp.Matrix) -> sp.Matrix:
    diagonal, _, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(row).convert_to(ZZ)
    )
    smith = diagonal.to_Matrix()
    transform = right.to_Matrix()
    rank = row.rank()
    if rank != 1 or row * transform != smith:
        raise ArithmeticError("Smith kernel replay failed")
    return transform[:, rank:]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i]
    ]


def _oriented_gram(gram: sp.Matrix, walls: Sequence[sp.Matrix]) -> sp.Matrix:
    wall_matrix = sp.Matrix.hstack(*walls)
    result = wall_matrix.T * gram * wall_matrix
    if any(
        result[i, j] < 0
        for i in range(result.rows)
        for j in range(i)
    ):
        raise ArithmeticError("the archived wall orientation is not acute")
    return result


def _root_graph(root_gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(root_gram.rows))
    for left in range(root_gram.rows):
        for right in range(left):
            pairing = int(root_gram[left, right])
            if pairing:
                graph.add_edge(right, left, pairing=pairing)
    return graph


def _type_ii_match(root_gram: sp.Matrix) -> tuple[list[int] | None, int]:
    if root_gram.shape != (12, 12):
        return None, 0
    source = _root_graph(root_gram)
    target = nx.Graph()
    target.add_nodes_from(range(12))
    target.add_edges_from(
        (left, right, {"pairing": 1})
        for left, right in KONDO_TYPE_II_EDGES
    )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        source,
        target,
        edge_match=lambda left, right: left["pairing"] == right["pairing"],
    )
    mappings = [tuple(mapping[i] for i in range(12)) for mapping in matcher.isomorphisms_iter()]
    return (list(min(mappings)) if mappings else None), len(mappings)


def _cone_rays(
    gram: sp.Matrix, boundary: sp.Matrix, roots: Sequence[sp.Matrix]
) -> list[dict[str, Any]]:
    rays: dict[tuple[int, ...], tuple[int, ...]] = {}
    for facets in itertools.combinations(range(len(roots)), 9):
        equations = sp.Matrix.vstack(
            boundary.T * gram,
            *(roots[index].T * gram for index in facets),
        )
        if equations.rank() != 10:
            continue
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            raise ArithmeticError("an extreme-ray system has wrong nullity")
        vector = nullspace[0]
        denominator = int(sp.ilcm(*(value.q for value in vector)))
        entries = [int(value * denominator) for value in vector]
        content = math.gcd(*(abs(value) for value in entries))
        vector = sp.Matrix([value // content for value in entries])
        pairings = [int((vector.T * gram * root)[0]) for root in roots]
        if all(value <= 0 for value in pairings):
            vector = -vector
            pairings = [-value for value in pairings]
        if not all(value >= 0 for value in pairings):
            continue
        key = tuple(_vector_payload(vector))
        rays[key] = tuple(
            index for index, value in enumerate(pairings) if value == 0
        )
    result = []
    for coordinates, active in sorted(rays.items()):
        vector = sp.Matrix(coordinates)
        square = int((vector.T * gram * vector)[0])
        result.append(
            {
                "vector": list(coordinates),
                "square": square,
                "kind": (
                    "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
                ),
                "active_root_positions_one_based": [index + 1 for index in active],
            }
        )
    return result


def _type_ii_saturation_uniqueness(
    complement_gram: sp.Matrix,
    root_coordinates: sp.Matrix,
    root_gram: sp.Matrix,
) -> dict[str, Any]:
    """Audit the integral U-gate omitted by a mere real Gram comparison."""

    root_basis = hermite_normal_form(root_coordinates)
    if root_basis.shape != (10, 10):
        raise ArithmeticError("the Type-II roots do not span rank ten")
    generated_gram = root_basis.T * complement_gram * root_basis
    generator_coordinates = root_basis.inv() * root_coordinates
    if any(value.q != 1 for value in generator_coordinates):
        raise ArithmeticError("HNF did not generate the root lattice")

    discriminant_classes: list[tuple[tuple[int, ...], int, bool]] = []
    for bits in itertools.product((0, 1), repeat=10):
        vector = sp.Matrix([sp.Rational(value, 2) for value in bits])
        if all(value.q == 1 for value in generated_gram * vector):
            square = int((vector.T * generated_gram * vector)[0])
            discriminant_classes.append((bits, square, square % 2 == 0))
    nonzero_isotropic = [
        bits
        for bits, _, isotropic in discriminant_classes
        if any(bits) and isotropic
    ]

    graph = _root_graph(root_gram)
    pivots = list(generator_coordinates.rref()[1])
    independent = generator_coordinates[:, pivots]
    actions = set()
    for mapping in nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        edge_match=lambda left, right: left["pairing"] == right["pairing"],
    ).isomorphisms_iter():
        images = sp.Matrix.hstack(
            *(generator_coordinates[:, mapping[index]] for index in pivots)
        )
        action = images * independent.inv()
        if not (
            all(value.q == 1 for value in action)
            and abs(int(action.det())) == 1
            and action.T * generated_gram * action == generated_gram
        ):
            raise ArithmeticError("a diagram symmetry did not lift integrally")
        image_keys = []
        for bits in nonzero_isotropic:
            image = action * sp.Matrix(
                [sp.Rational(value, 2) for value in bits]
            )
            image_keys.append(
                tuple(int(2 * (value - sp.floor(value))) % 2 for value in image)
            )
        actions.add(tuple(image_keys))

    transitive = bool(nonzero_isotropic) and {
        action[0] for action in actions
    } == set(nonzero_isotropic)
    return {
        "root_span_index_in_E10": abs(int(root_basis.det())),
        "root_lattice_determinant": int(generated_gram.det()),
        "root_lattice_smith_diagonal": _smith_diagonal(generated_gram),
        "discriminant_group": "(Z/2)^2",
        "discriminant_classes": [
            {
                "twice_coordinates_mod_two": list(bits),
                "square_representative": square,
                "isotropic_for_even_overlattice": isotropic,
            }
            for bits, square, isotropic in discriminant_classes
        ],
        "nonzero_isotropic_index_two_extensions": len(nonzero_isotropic),
        "weighted_diagram_automorphism_order": len(
            list(
                nx.algorithms.isomorphism.GraphMatcher(
                    graph,
                    graph,
                    edge_match=lambda left, right: left["pairing"]
                    == right["pairing"],
                ).isomorphisms_iter()
            )
        ),
        "diagram_group_transitive_on_even_unimodular_extensions": transitive,
        "unique_E10_saturation_up_to_diagram_automorphism": (
            len(nonzero_isotropic) == 2 and transitive
        ),
    }


def _candidate_audit(
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    wall_gram: sp.Matrix,
    boundary_index: int,
) -> dict[str, Any]:
    boundary = walls[boundary_index]
    root_indices: list[int] = []
    exceptional_indices: list[int] = []
    profile: Counter[tuple[int, int]] = Counter()
    failures = []
    for index, wall in enumerate(walls):
        if index == boundary_index:
            continue
        square = int(wall_gram[index, index])
        pairing = int(wall_gram[boundary_index, index])
        profile[(pairing, square)] += 1
        if (pairing, square) == (0, -2):
            root_indices.append(index)
        elif (pairing, square) == (1, -1):
            exceptional_indices.append(index)
        else:
            failures.append(index)

    record: dict[str, Any] = {
        "boundary_facet_one_based": boundary_index + 1,
        "boundary_vector": _vector_payload(boundary),
        "incidence_profile_pairing_norm": {
            f"{pairing},{square}": count
            for (pairing, square), count in sorted(profile.items())
        },
        "coble_incidence_pass": not failures,
        "forbidden_facet_indices_one_based": [index + 1 for index in failures],
    }
    if failures:
        return record

    complement_basis = _kernel_basis(boundary.T * gram)
    complement_gram = complement_basis.T * gram * complement_basis
    splitting = sp.Matrix.hstack(boundary, complement_basis)
    roots = [walls[index] for index in root_indices]
    exceptional = [walls[index] for index in exceptional_indices]
    root_matrix = sp.Matrix.hstack(*roots)
    root_coordinates = complement_gram.inv() * (
        complement_basis.T * gram * root_matrix
    )
    if any(value.q != 1 for value in root_coordinates):
        raise ArithmeticError("root left b-perp")
    root_gram = root_matrix.T * gram * root_matrix

    rays = _cone_rays(gram, boundary, roots)
    ray_histogram = Counter(ray["kind"] for ray in rays)
    cusp_vectors = {
        tuple(ray["vector"]) for ray in rays if ray["kind"] == "ideal"
    }
    extracted = {
        tuple(_vector_payload(wall + boundary)): index
        for index, wall in zip(exceptional_indices, exceptional)
    }
    cusp_match = cusp_vectors == set(extracted)
    type_ii_permutation, type_ii_isomorphism_count = _type_ii_match(root_gram)

    if type_ii_permutation is not None:
        saturation = _type_ii_saturation_uniqueness(
            complement_gram, root_coordinates, root_gram
        )
    else:
        saturation = None

    if type_ii_permutation is not None and cusp_match:
        classification = {
            "status": "nonempty_complex_Coble_carrier",
            "kondo_type": "II",
            "first_surface_input": (
                "Kondo's existence/classification of the complex one-boundary "
                "Type-II Coble surface"
            ),
        }
    elif cusp_match and len(roots) != 12:
        classification = {
            "status": "Coble_carrier_obstructed",
            "reason": (
                "the finite-volume CM root chamber would force finite Aut, "
                "whereas Kondo's complex one-boundary finite-Aut catalogue "
                "has exactly 12 roots in both Types I and II"
            ),
            "observed_root_count": len(roots),
        }
    else:
        classification = {
            "status": "Coble_carrier_obstructed",
            "reason": "the exact cusp regeneration or Kondo diagram gate failed",
        }

    record.update(
        {
            "root_facet_indices_one_based": [index + 1 for index in root_indices],
            "exceptional_facet_indices_one_based": [
                index + 1 for index in exceptional_indices
            ],
            "b_perp": {
                "basis_columns": _matrix_payload(complement_basis),
                "gram": _matrix_payload(complement_gram),
                "determinant": int(complement_gram.det()),
                "inertia": list(_exact_inertia(complement_gram)),
                "even": all(
                    int(complement_gram[i, i]) % 2 == 0 for i in range(10)
                ),
                "direct_splitting_determinant": int(splitting.det()),
                "is_E10": (
                    complement_gram.det() == -1
                    and _exact_inertia(complement_gram) == (1, 9, 0)
                    and all(
                        int(complement_gram[i, i]) % 2 == 0
                        for i in range(10)
                    )
                ),
            },
            "root_gram": _matrix_payload(root_gram),
            "root_span_smith_diagonal_in_b_perp": _smith_diagonal(
                root_coordinates
            ),
            "cm_cone": {
                "ray_count": len(rays),
                "kind_histogram": dict(sorted(ray_histogram.items())),
                "rays": rays,
            },
            "cusp_regeneration": {
                "extracted_E_plus_b_count": len(extracted),
                "complete_primitive_isotropic_cusp_count": len(cusp_vectors),
                "sets_equal": cusp_match,
                "cusp_to_exceptional_facet_one_based": {
                    ",".join(map(str, cusp)): extracted[cusp] + 1
                    for cusp in sorted(cusp_vectors)
                    if cusp in extracted
                },
            },
            "kondo_type_II_root_gram_match": type_ii_permutation is not None,
            "root_to_kondo_type_II_zero_based_permutation": type_ii_permutation,
            "type_II_weighted_graph_isomorphism_count": type_ii_isomorphism_count,
            "integral_saturation_uniqueness": saturation,
            "complete_package_regenerated_from_P": cusp_match,
            "surface_classification_gate": classification,
        }
    )
    return record


def _presentation_audit(
    chamber_id: str,
    presentation_id: str,
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    source: str,
) -> dict[str, Any]:
    if not (
        abs(int(gram.det())) == 1
        and _exact_inertia(gram) == (1, 10, 0)
        and any(int(gram[i, i]) % 2 for i in range(11))
    ):
        raise ArithmeticError("presentation is not odd unimodular I_(1,10)")
    wall_gram = _oriented_gram(gram, walls)
    if any(wall_gram[i, i] not in (-1, -2) for i in range(len(walls))):
        raise ArithmeticError("primitive wall alphabet is not {-1,-2}")
    characteristic = [
        index
        for index, wall in enumerate(walls)
        if wall_gram[index, index] == -1 and _characteristic(gram, wall)
    ]
    candidates = [
        _candidate_audit(gram, walls, wall_gram, index)
        for index in characteristic
    ]
    return {
        "chamber_id": chamber_id,
        "presentation_id": presentation_id,
        "source": source,
        "facet_count": len(walls),
        "lattice_gram": _matrix_payload(gram),
        "primitive_norm_histogram": dict(
            sorted(Counter(int(wall_gram[i, i]) for i in range(len(walls))).items())
        ),
        "all_off_diagonal_pairings_nonnegative": True,
        "characteristic_norm_minus_one_facet_count": len(characteristic),
        "candidate_audits": candidates,
    }


def _fold_presentations() -> list[dict[str, Any]]:
    presentations = []
    graph_data = json.loads(R12_GRAPH.read_text(encoding="utf-8"))
    for stock in graph_data["stocks"]:
        for survivor in stock["survivors"]:
            chamber_id = GRAPH_IDS[(stock["stock_key"], survivor["class_index"])]
            fixed = survivor["fixed_lattice"]
            fixed_roots = sp.Matrix(fixed["orbit_root_coordinate_columns"])
            for position, overclosure in enumerate(
                fixed["maximal_stable_overlattices"]
            ):
                gram = sp.Matrix(overclosure["gram"])
                basis = sp.Matrix(overclosure["basis_columns_in_fixed_lattice"])
                coordinates = basis.inv() * fixed_roots
                walls = [_primitive(coordinates[:, index]) for index in range(coordinates.cols)]
                presentations.append(
                    _presentation_audit(
                        chamber_id,
                        f"{chamber_id}:overclosure-{position}",
                        gram,
                        walls,
                        str(R12_GRAPH.relative_to(ROOT)),
                    )
                )

    data_59 = json.loads(R12_59.read_text(encoding="utf-8"))
    for case in data_59["cases"]:
        fixed = case["fixed_lattice"]
        fixed_roots = sp.Matrix(fixed["orbit_root_coordinate_columns"])
        for position, overclosure in enumerate(fixed["maximal_stable_overlattices"]):
            gram = sp.Matrix(overclosure["gram"])
            basis = sp.Matrix(overclosure["basis_columns_in_fixed_lattice"])
            coordinates = basis.inv() * fixed_roots
            walls = [_primitive(coordinates[:, index]) for index in range(coordinates.cols)]
            presentations.append(
                _presentation_audit(
                    case["id"],
                    f"{case['id']}:overclosure-{position}",
                    gram,
                    walls,
                    str(R12_59.relative_to(ROOT)),
                )
            )
    return presentations


def _d5_presentation() -> dict[str, Any]:
    data = json.loads(D5.read_text(encoding="utf-8"))
    gram = sp.Matrix(data["saturated_fixed_lattice"]["target_half_scaled_gram"])
    walls = []
    for record in data["primitive_mirrors"]["records"]:
        vector = sp.Matrix(record["primitive_I_1_10_coordinates"])
        sign = 1 if sp.Rational(record["presented_equals_scale_times_primitive"]) > 0 else -1
        walls.append(sign * vector)
    return _presentation_audit(
        "D16-D5-H10-30",
        "D16-D5-H10-30:canonical-half-scale",
        gram,
        walls,
        str(D5.relative_to(ROOT)),
    )


def build_certificate() -> dict[str, Any]:
    presentations = [*_fold_presentations(), _d5_presentation()]
    chamber_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    candidates = []
    for presentation in presentations:
        chamber_groups[presentation["chamber_id"]].append(presentation)
        candidates.extend(presentation["candidate_audits"])

    incidence_pass = [item for item in candidates if item["coble_incidence_pass"]]
    regenerated = [
        item for item in incidence_pass if item["complete_package_regenerated_from_P"]
    ]
    positives = [
        item
        for item in regenerated
        if item["surface_classification_gate"]["status"]
        == "nonempty_complex_Coble_carrier"
    ]
    obstructions = [
        item
        for item in regenerated
        if item["surface_classification_gate"]["status"]
        == "Coble_carrier_obstructed"
    ]

    chamber_summary = {}
    for chamber_id, group in sorted(chamber_groups.items()):
        group_candidates = [
            candidate
            for presentation in group
            for candidate in presentation["candidate_audits"]
        ]
        chamber_summary[chamber_id] = {
            "presentation_count": len(group),
            "characteristic_candidate_count": len(group_candidates),
            "incidence_pass_count": sum(
                candidate["coble_incidence_pass"] for candidate in group_candidates
            ),
            "regenerated_package_count": sum(
                candidate.get("complete_package_regenerated_from_P", False)
                for candidate in group_candidates
            ),
            "classification_statuses": sorted(
                {
                    candidate["surface_classification_gate"]["status"]
                    for candidate in group_candidates
                    if candidate.get("complete_package_regenerated_from_P", False)
                }
            ),
        }

    expected = {
        "six_chamber_types": len(chamber_groups) == 6,
        "twenty_four_integral_presentations": len(presentations) == 24,
        "one_hundred_fifty_four_characteristic_candidates": len(candidates) == 154,
        "two_pass_incidence": len(incidence_pass) == 2,
        "two_regenerate_complete_package": len(regenerated) == 2,
        "one_kondo_positive": len(positives) == 1,
        "one_kondo_obstruction": len(obstructions) == 1,
        "positive_is_D16_D5": any(
            presentation["chamber_id"] == "D16-D5-H10-30"
            and any(
                candidate.get("surface_classification_gate", {}).get("status")
                == "nonempty_complex_Coble_carrier"
                for candidate in presentation["candidate_audits"]
            )
            for presentation in presentations
        ),
        "D16_boundary_is_unique_facet_26": next(
            presentation
            for presentation in presentations
            if presentation["chamber_id"] == "D16-D5-H10-30"
        )["candidate_audits"][0]["boundary_facet_one_based"]
        == 26,
        "D16_unique_E10_saturation": next(
            candidate
            for presentation in presentations
            if presentation["chamber_id"] == "D16-D5-H10-30"
            for candidate in presentation["candidate_audits"]
        )["integral_saturation_uniqueness"][
            "unique_E10_saturation_up_to_diagram_automorphism"
        ],
    }
    if not all(expected.values()):
        raise ArithmeticError(f"hostile audit changed: {expected}")

    return {
        "schema": SCHEMA,
        "scope": {
            "input": (
                "all pre-existing local H10 chamber artifacts with an explicit "
                "odd-unimodular I_(1,10) presentation or reflection-stable "
                "odd-unimodular overclosure"
            ),
            "chamber_types": sorted(chamber_groups),
            "failure_means": (
                "obstruction only to a complex finite-Aut one-boundary Coble "
                "full-nef carrier"
            ),
            "not_claimed": (
                "universality for Coble surfaces, other surfaces, or abstract "
                "arithmetic Coxeter polytopes"
            ),
        },
        "logic_layers": {
            "P_only": [
                "odd-unimodular normalization and primitive acute walls",
                "characteristic norm-minus-one boundary search",
                "b-perp=E10 test",
                "root/exceptional incidence split",
                "complete exact CM ray enumeration",
                "cusp-to-E regeneration and integral saturation uniqueness",
            ],
            "first_surface_classification_input": (
                "Kondo's theorem that the complex finite-Aut one-boundary "
                "Coble types are exactly I and II, both with 12 roots, plus "
                "his explicit nonempty Type-II example"
            ),
            "no_surface_first_wall_list_imported": True,
            "no_project_script_imports": True,
            "no_import_cycle": True,
        },
        "summary": {
            "chamber_type_count": len(chamber_groups),
            "integral_presentation_count": len(presentations),
            "characteristic_norm_minus_one_candidate_count": len(candidates),
            "Coble_incidence_pass_count": len(incidence_pass),
            "complete_P_package_regeneration_count": len(regenerated),
            "Kondo_nonempty_complex_carrier_count": len(positives),
            "Kondo_catalogue_obstruction_count": len(obstructions),
            "decision": (
                "D16-D5-H10-30 is the unique positive in this frozen local "
                "denominator; P10/12 is the only other arithmetic survivor "
                "and is obstructed by its ten-root CM chamber"
            ),
        },
        "chamber_summary": chamber_summary,
        "presentations": presentations,
        "checks": expected,
        "all_checks_pass": all(expected.values()),
        "sources": {
            str(D5.relative_to(ROOT)): _sha256(D5),
            str(R12_GRAPH.relative_to(ROOT)): _sha256(R12_GRAPH),
            str(R12_59.relative_to(ROOT)): _sha256(R12_59),
            "kondo_primary": {
                "arxiv": "2102.00591",
                "machine_role": "classification label only; not parsed by replay",
                "claims_used": [
                    "Main Theorem/Table 2",
                    "Example TypeII",
                    "CME10, NefCone, Vinberg2, Halphen, fibration-isotropic",
                ],
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
