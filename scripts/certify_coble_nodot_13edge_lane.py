#!/usr/bin/env python3
"""Certify the finite 13-edge small-pairing Coble--Mukai lane.

This is intentionally *not* a census of all ``D_1_split^nodot`` inputs.
It exhausts connected twelve-node ordinary support graphs with exactly
thirteen edges, maximum ordinary degree three, and exactly two or three
edges promoted from pairing one to pairing two.  The first scout leaves 38
weighted determinant-four Lorentzian presentations.  This certifier takes
those through the integral E10-overlattice, exact cddlib chamber, cusp-lift,
full H10 chamber, diagram-isomorphism, and (when identified with an archived
positive row) complete split-K3 compiler gates.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
from typing import Any, Iterable

from flint import fmpz_mat
import networkx as nx
import sympy as sp

from scripts.certify_coble_one_boundary_raw_k3_promotion import (
    _all_subgroups,
    _discriminant_group,
    _row_hnf_basis,
    build_certificate as build_raw_promotion,
)
from scripts.scout_coble_nodot_weighted_graphs import build as build_scout
from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _parse_cdd_rays,
    _run_cdd,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "enumerations" / "coble_nodot_13edge_lane.json"
SOURCE_FULL = (
    ROOT / "data" / "enumerations" / "coble_one_boundary_full_ns_h10.json"
)
SCHEMA = "coxeter.coble_nodot_13edge_lane.v1"


def _matrix(values: Iterable[Iterable[Any]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in values])


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("matrix expected to be integral")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _rational_rows(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix.tolist()]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    from sympy.matrices.normalforms import smith_normal_form
    from sympy.polys.domains import ZZ

    smith = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(smith[index, index])) for index in range(min(smith.shape))]


def _column_span_index(matrix: sp.Matrix) -> int | None:
    """Return the index of the column span, or ``None`` when not full rank."""

    if matrix.rank() != matrix.rows:
        return None
    diagonal = [value for value in _smith_diagonal(matrix) if value]
    if len(diagonal) != matrix.rows:
        raise ArithmeticError("Smith form lost a full-rank invariant")
    return math.prod(diagonal)


def _weighted_graph(root_gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(root_gram.rows))
    for left in range(root_gram.rows):
        for right in range(left + 1, root_gram.rows):
            pairing = int(root_gram[left, right])
            if pairing:
                graph.add_edge(left, right, pairing=pairing)
    return graph


def _weighted_isomorphic(left: nx.Graph, right: nx.Graph) -> bool:
    return nx.is_isomorphic(
        left,
        right,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "pairing", None
        ),
    )


def _deduplicate_weighted(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    classes: list[dict[str, Any]] = []
    for candidate in candidates:
        root_gram = -_matrix(candidate["negative_cartan"])
        graph = _weighted_graph(root_gram)
        for record in classes:
            if _weighted_isomorphic(graph, record["graph"]):
                record["members"].append(candidate)
                break
        else:
            classes.append(
                {"representative": candidate, "members": [candidate], "graph": graph}
            )
    return classes


def _quotient_lattice(root_gram: sp.Matrix) -> dict[str, Any]:
    """Take the saturated radical quotient of a 12-root presentation."""

    if root_gram.shape != (12, 12) or root_gram.rank() != 10:
        raise ArithmeticError("expected a rank-ten twelve-root Gram")
    gram_flint = fmpz_mat(_integer_rows(root_gram))
    _, transform_flint = gram_flint.transpose().hnf(transform=True)
    basis = sp.Matrix(transform_flint.transpose().tolist())
    if abs(int(basis.det())) != 1:
        raise ArithmeticError("HNF radical transform is not unimodular")
    transformed = basis.T * root_gram * basis
    if any(
        transformed[left, right]
        for left in range(12)
        for right in range(12)
        if left >= 10 or right >= 10
    ):
        raise ArithmeticError("saturated radical was not isolated")
    quotient_gram = transformed[:10, :10]
    root_coordinates = basis.inv()[:10, :]
    if root_coordinates.T * quotient_gram * root_coordinates != root_gram:
        raise ArithmeticError("root coordinates do not replay the Gram")
    if exact_inertia(quotient_gram) != (1, 9, 0):
        raise ArithmeticError("radical quotient is not Lorentzian E10-index-two")
    return {
        "gram": quotient_gram,
        "root_coordinates": root_coordinates,
        "radical_basis": basis[:, 10:],
        "determinant": int(quotient_gram.det()),
        "smith_diagonal_absolute": _smith_diagonal(quotient_gram),
    }


def _even_unimodular_extensions(
    quotient_gram: sp.Matrix, root_coordinates: sp.Matrix
) -> list[dict[str, Any]]:
    """Enumerate every even index-two extension of the determinant-four quotient."""

    if abs(int(quotient_gram.det())) != 4:
        return []
    extensions: list[dict[str, Any]] = []
    for mask in range(1, 1 << 10):
        residue = sp.Matrix([(mask >> index) & 1 for index in range(10)])
        if any(int(value) % 2 for value in quotient_gram * residue):
            continue
        if int((residue.T * quotient_gram * residue)[0]) % 8:
            continue
        pivot = next(index for index, value in enumerate(residue) if value)
        transition = sp.eye(10)
        transition[:, pivot] = residue / 2
        extension_gram = transition.T * quotient_gram * transition
        extension_roots = transition.inv() * root_coordinates
        if not (
            abs(sp.Rational(transition.det())) == sp.Rational(1, 2)
            and all(value.q == 1 for value in extension_gram)
            and all(value.q == 1 for value in extension_roots)
            and all(int(extension_gram[index, index]) % 2 == 0 for index in range(10))
            and int(extension_gram.det()) == -1
            and exact_inertia(extension_gram) == (1, 9, 0)
            and extension_roots.T * extension_gram * extension_roots == root_coordinates.T * quotient_gram * root_coordinates
        ):
            raise ArithmeticError("an alleged E10 extension failed exact replay")
        extensions.append(
            {
                "isotropic_residue_binary": [int(value) for value in residue],
                "gram": extension_gram.applyfunc(int),
                "root_coordinates": extension_roots.applyfunc(int),
                "transition": transition,
            }
        )
    return extensions


def _all_even_unimodular_extensions(
    quotient_gram: sp.Matrix, root_coordinates: sp.Matrix
) -> list[dict[str, Any]]:
    """Enumerate every even unimodular finite-index extension.

    Even overlattices of an even lattice are in bijection with isotropic
    subgroups of its discriminant quadratic form.  A unimodular extension
    requires a subgroup of order ``sqrt(abs(det))``.  The determinant-four
    helper above is retained as a small independent replay; this routine is
    the general gate used by broader lanes.
    """

    determinant = abs(int(quotient_gram.det()))
    index = math.isqrt(determinant)
    if index * index != determinant:
        return []
    group, order, quadratic = _discriminant_group(quotient_gram)
    isotropic = [
        subgroup
        for subgroup in _all_subgroups(group, order)
        if len(subgroup) == index
        and all(quadratic(value) == 0 for value in subgroup)
    ]
    extensions: dict[tuple[tuple[str, ...], ...], dict[str, Any]] = {}
    identity_rows = sp.eye(quotient_gram.rows)
    for subgroup in isotropic:
        nonzero = [value for value in subgroup if any(value)]
        generators = (
            sp.Matrix.vstack(identity_rows, sp.Matrix(nonzero))
            if nonzero
            else identity_rows
        )
        transition = _row_hnf_basis(generators).T
        extension_gram = transition.T * quotient_gram * transition
        extension_roots = transition.inv() * root_coordinates
        if not (
            abs(sp.Rational(transition.det())) == sp.Rational(1, index)
            and all(value.q == 1 for value in extension_gram)
            and all(value.q == 1 for value in extension_roots)
            and all(
                int(extension_gram[position, position]) % 2 == 0
                for position in range(extension_gram.rows)
            )
            and int(extension_gram.det()) == -1
            and exact_inertia(extension_gram) == (1, 9, 0)
            and extension_roots.T * extension_gram * extension_roots
            == root_coordinates.T * quotient_gram * root_coordinates
        ):
            raise ArithmeticError("an isotropic subgroup failed E10 replay")
        key = tuple(tuple(str(value) for value in row) for row in transition.tolist())
        extensions[key] = {
            "isotropic_subgroup": [
                [str(value) for value in element]
                for element in sorted(
                    subgroup, key=lambda value: tuple(map(str, value))
                )
            ],
            "index_over_root_quotient": index,
            "gram": extension_gram.applyfunc(int),
            "root_coordinates": extension_roots.applyfunc(int),
            "transition": transition,
        }
    return [extensions[key] for key in sorted(extensions)]


def _cone_text(inequality: sp.Matrix) -> str:
    lines = ["H-representation", "begin", f"{inequality.rows} {inequality.cols + 1} rational"]
    for row in range(inequality.rows):
        lines.append("0 " + " ".join(str(inequality[row, column]) for column in range(inequality.cols)))
    lines.append("end")
    return "\n".join(lines) + "\n"


def _analyze_cone(
    gram: sp.Matrix,
    normal_coordinates: sp.Matrix,
    cdd: str,
) -> dict[str, Any]:
    rank = gram.rows
    if normal_coordinates.rows != rank:
        raise ArithmeticError("normal coordinate rank mismatch")
    inequality = normal_coordinates.T * gram
    if inequality.rank() != rank:
        raise ArithmeticError("inequality rows do not span the dual space")
    output = _run_cdd(_cone_text(inequality), cdd)
    rays_raw = _parse_cdd_rays(output, rank)
    rays = [sp.Matrix(ray) for ray in rays_raw]
    records: list[dict[str, Any]] = []
    type_count: Counter[str] = Counter()
    active_histogram: Counter[tuple[str, int]] = Counter()
    facet_active_vectors: list[list[sp.Matrix]] = [
        [] for _ in range(inequality.rows)
    ]
    for vector in rays:
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("cdd emitted an infeasible ray")
        active = [index for index, value in enumerate(slacks) if value == 0]
        square = int((vector.T * gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        type_count[kind] += 1
        active_histogram[(kind, len(active))] += 1
        for index in active:
            facet_active_vectors[index].append(vector)
        records.append(
            {
                "coordinates": [int(value) for value in vector],
                "square": square,
                "type": kind,
                "active_facets_zero_based": active,
            }
        )
    facet_ranks = [
        int(sp.Matrix.hstack(*vectors).rank()) if vectors else 0
        for vectors in facet_active_vectors
    ]
    irredundant = all(value == rank - 1 for value in facet_ranks)
    no_hyperideal = bool(records) and type_count["hyperideal"] == 0
    ideal_records = [record for record in records if record["type"] == "ideal"]
    cusp_audits = []
    for record in ideal_records:
        active = record["active_facets_zero_based"]
        active_gram = normal_coordinates[:, active].T * gram * normal_coordinates[:, active]
        inertia = exact_inertia(active_gram)
        cusp_audits.append(
            {
                "active_facets_zero_based": active,
                "active_gram_inertia": list(inertia),
                "parabolic_rank_n_minus_one": (
                    inertia[0] == 0 and inertia[1] == rank - 2
                ),
            }
        )
    interior = sum(rays, sp.zeros(rank, 1)) if rays else sp.zeros(rank, 1)
    interior_slacks = inequality * interior
    interior_square = int((interior.T * gram * interior)[0])
    same_time_cone = True
    finite_ray = next((sp.Matrix(record["coordinates"]) for record in records if record["type"] == "finite"), None)
    if finite_ray is not None:
        same_time_cone = all(
            int((finite_ray.T * gram * sp.Matrix(record["coordinates"]))[0]) > 0
            for record in records
            if record["type"] in {"finite", "ideal"}
        )
    return {
        "rank": rank,
        "facet_count": inequality.rows,
        "extreme_ray_count": len(records),
        "finite_vertex_count": type_count["finite"],
        "ideal_vertex_count": type_count["ideal"],
        "hyperideal_vertex_count": type_count["hyperideal"],
        "active_facet_count_histogram": {
            f"{kind}:{count}": multiplicity
            for (kind, count), multiplicity in sorted(active_histogram.items())
        },
        "facet_active_span_ranks": facet_ranks,
        "all_input_walls_irredundant": irredundant,
        "complete_cdd_V_representation_has_no_hyperideal_ray": no_hyperideal,
        "compact_by_cdd_vertex_types": no_hyperideal and not ideal_records,
        "strict_interior_from_ray_sum": all(value > 0 for value in interior_slacks),
        "strict_interior_square": interior_square,
        "all_causal_rays_same_time_cone": same_time_cone,
        "all_ideal_supports_parabolic_rank_n_minus_one": all(
            audit["parabolic_rank_n_minus_one"] for audit in cusp_audits
        ),
        "cusp_audits": cusp_audits,
        "rays": records,
    }


def _literal_vinberg(
    gram: sp.Matrix,
    normal_coordinates: sp.Matrix,
    cone: dict[str, Any],
) -> dict[str, Any]:
    """Apply Vinberg's elliptic-edge and complete parabolic-neighborhood tests."""

    dimension = gram.rows - 1
    wall_gram = normal_coordinates.T * gram * normal_coordinates
    positive_gram = -wall_gram
    endpoint_masks = [
        sum(1 << index for index in record["active_facets_zero_based"])
        for record in cone["rays"]
    ]
    endpoint_kinds = [record["type"] for record in cone["rays"]]
    if "hyperideal" in endpoint_kinds:
        return {
            "literal_vinberg_passes": False,
            "stopped_at_hyperideal_endpoint": True,
            "elliptic_edge_count": None,
            "connected_affine_component_count": None,
        }
    finite_masks = {
        mask
        for mask, kind in zip(endpoint_masks, endpoint_kinds, strict=True)
        if kind == "finite"
    }
    ideal_masks = {
        mask
        for mask, kind in zip(endpoint_masks, endpoint_kinds, strict=True)
        if kind == "ideal"
    }
    enumerator = DiagramEnumerator(positive_gram, dimension=dimension)
    finite_edges = enumerator.finite_subdiagrams(dimension - 1)
    finite_vertices = enumerator.finite_subdiagrams(dimension)
    finite_endpoint_match = set(finite_vertices) == finite_masks

    edge_endpoints: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(endpoint_masks)):
        for left in range(right):
            intersection = endpoint_masks[left] & endpoint_masks[right]
            if intersection.bit_count() >= dimension - 1:
                edge_endpoints[intersection].append((left, right))
    elliptic_edge_pass = (
        set(edge_endpoints) == set(finite_edges)
        and all(mask.bit_count() == dimension - 1 for mask in edge_endpoints)
        and all(len(pairs) == 1 for pairs in edge_endpoints.values())
    )

    affine = enumerator.connected_affine_subdiagrams()
    affine_neighborhoods: dict[int, int] = {}
    bad_affine_neighborhoods = []
    for component_mask in affine:
        component = tuple(
            index
            for index in range(normal_coordinates.cols)
            if component_mask & (1 << index)
        )
        zero_neighborhood = tuple(
            facet
            for facet in range(normal_coordinates.cols)
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        neighborhood = component_mask | sum(
            1 << index for index in zero_neighborhood
        )
        affine_neighborhoods[component_mask] = neighborhood
        if neighborhood not in ideal_masks:
            bad_affine_neighborhoods.append(component_mask)

    bad_ideal_cusps = []
    for cusp_mask in ideal_masks:
        components = enumerator._components(cusp_mask, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components)
            == dimension - 1
            and all(component in affine for component in components)
            and all(
                affine_neighborhoods.get(component) == cusp_mask
                for component in components
            )
        ):
            bad_ideal_cusps.append(cusp_mask)
    parabolic_pass = not bad_affine_neighborhoods and not bad_ideal_cusps
    passes = finite_endpoint_match and elliptic_edge_pass and parabolic_pass
    return {
        "literal_vinberg_passes": passes,
        "stopped_at_hyperideal_endpoint": False,
        "finite_endpoint_diagrams_equal_cdd_endpoints": finite_endpoint_match,
        "every_elliptic_edge_has_exactly_two_endpoints": elliptic_edge_pass,
        "every_affine_zero_neighborhood_is_a_cdd_cusp": not bad_affine_neighborhoods,
        "every_cdd_cusp_has_complete_affine_decomposition": not bad_ideal_cusps,
        "elliptic_edge_count": len(finite_edges),
        "finite_endpoint_count": len(finite_masks),
        "ideal_endpoint_count": len(ideal_masks),
        "connected_affine_component_count": len(affine),
        "bad_affine_zero_neighborhood_count": len(bad_affine_neighborhoods),
        "bad_ideal_cusp_count": len(bad_ideal_cusps),
    }


def _full_cusp_lift(
    e10_gram: sp.Matrix,
    roots: sp.Matrix,
    cm_cone: dict[str, Any],
    cdd: str,
) -> dict[str, Any]:
    ideal = [record for record in cm_cone["rays"] if record["type"] == "ideal"]
    ambient_gram = sp.diag(e10_gram, sp.Matrix([[-1]]))
    boundary = sp.zeros(11, 1)
    boundary[10] = 1
    normals = [roots[:, index].col_join(sp.Matrix([0])) for index in range(roots.cols)]
    cusp_lifts = []
    for record in ideal:
        cusp = sp.Matrix(record["coordinates"])
        divisor = cusp.col_join(sp.Matrix([-1]))
        if int((divisor.T * ambient_gram * divisor)[0]) != -1:
            raise ArithmeticError("a cusp lift does not have square minus one")
        cusp_lifts.append(divisor)
    full_normals = sp.Matrix.hstack(boundary, *normals, *cusp_lifts)
    full = _analyze_cone(ambient_gram, full_normals, cdd)
    full_vinberg = _literal_vinberg(ambient_gram, full_normals, full)
    wall_gram = full_normals.T * ambient_gram * full_normals
    characteristic = []
    for index in range(full_normals.cols):
        vector = full_normals[:, index]
        characteristic.append(
            all(
                int((ambient_gram[column, column] - (ambient_gram * vector)[column])) % 2 == 0
                for column in range(11)
            )
        )
    return {
        "ambient_gram": ambient_gram,
        "normal_coordinates": full_normals,
        "wall_gram": wall_gram,
        "boundary_wall_zero_based": 0,
        "root_wall_positions_zero_based": list(range(1, 13)),
        "cusp_lift_positions_zero_based": list(range(13, full_normals.cols)),
        "characteristic_wall_positions_zero_based": [
            index for index, value in enumerate(characteristic) if value
        ],
        "wall_coordinate_span_index_in_I_1_10": _column_span_index(full_normals),
        "cone": full,
        "literal_vinberg": full_vinberg,
    }


def _known_root_graphs() -> dict[str, nx.Graph]:
    payload = json.loads(SOURCE_FULL.read_text(encoding="utf-8"))
    result = {}
    for case_key, case in payload["cases"].items():
        wall_gram = _matrix(case["negative_wall_gram"])
        roots = [index for index in range(wall_gram.rows) if wall_gram[index, index] == -2]
        root_gram = wall_gram.extract(roots, roots)
        label = f"Type {case_key}"
        result[label] = _weighted_graph(root_gram)
    return result


def _full_wall_graph(wall_gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for index in range(wall_gram.rows):
        graph.add_node(index, norm=str(wall_gram[index, index]))
    for left in range(wall_gram.rows):
        for right in range(left + 1, wall_gram.rows):
            graph.add_edge(left, right, pairing=str(wall_gram[left, right]))
    return graph


def _full_isomorphic(left: sp.Matrix, right: sp.Matrix) -> bool:
    return nx.is_isomorphic(
        _full_wall_graph(left),
        _full_wall_graph(right),
        node_match=nx.algorithms.isomorphism.categorical_node_match("norm", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("pairing", None),
    )


def build(cdd: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    scout = build_scout(3, "13:13", 3, 100)
    raw = [
        candidate
        for candidate in scout["exact_labeled_survivor_samples"]
        if abs(int(candidate["radical_quotient_determinant"])) == 4
    ]
    if not (
        len(raw) == 38
        and scout["absolute_discriminant_four_count_by_edges_and_parallel"]
        == {"13:2": 8, "13:3": 30}
    ):
        raise ArithmeticError("the 38-candidate scout lane changed")
    classes = _deduplicate_weighted(raw)
    known = _known_root_graphs()
    class_records = []
    actual_full: list[tuple[int, int, sp.Matrix]] = []
    for class_index, diagram_class in enumerate(classes):
        representative = diagram_class["representative"]
        root_gram = -_matrix(representative["negative_cartan"])
        quotient = _quotient_lattice(root_gram)
        extensions = _even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        extension_records = []
        for extension_index, extension in enumerate(extensions):
            cm = _analyze_cone(
                extension["gram"], extension["root_coordinates"], cdd
            )
            cm_vinberg = _literal_vinberg(
                extension["gram"], extension["root_coordinates"], cm
            )
            cm_pass = (
                cm["all_input_walls_irredundant"]
                and cm[
                    "complete_cdd_V_representation_has_no_hyperideal_ray"
                ]
                and not cm["compact_by_cdd_vertex_types"]
                and cm["strict_interior_from_ray_sum"]
                and cm["strict_interior_square"] > 0
                and cm["all_causal_rays_same_time_cone"]
                and cm["all_ideal_supports_parabolic_rank_n_minus_one"]
                and cm_vinberg["literal_vinberg_passes"]
            )
            full_record = None
            full_pass = False
            if cm_pass:
                full_package = _full_cusp_lift(
                    extension["gram"], extension["root_coordinates"], cm, cdd
                )
                full = full_package["cone"]
                full_vinberg = full_package["literal_vinberg"]
                full_pass = (
                    full["all_input_walls_irredundant"]
                    and full[
                        "complete_cdd_V_representation_has_no_hyperideal_ray"
                    ]
                    and full["strict_interior_from_ray_sum"]
                    and full["strict_interior_square"] > 0
                    and full["all_causal_rays_same_time_cone"]
                    and full[
                        "all_ideal_supports_parabolic_rank_n_minus_one"
                    ]
                    and full_vinberg["literal_vinberg_passes"]
                    and full_package["characteristic_wall_positions_zero_based"] == [0]
                    and full_package["wall_coordinate_span_index_in_I_1_10"] == 1
                )
                full_record = {
                    "wall_count": full_package["normal_coordinates"].cols,
                    "characteristic_wall_positions_zero_based": full_package[
                        "characteristic_wall_positions_zero_based"
                    ],
                    "wall_coordinate_span_index_in_I_1_10": full_package[
                        "wall_coordinate_span_index_in_I_1_10"
                    ],
                    "wall_gram": _integer_rows(full_package["wall_gram"]),
                    "cone": {key: value for key, value in full.items() if key != "rays"},
                    "literal_vinberg": full_vinberg,
                }
                if full_pass:
                    actual_full.append(
                        (class_index, extension_index, full_package["wall_gram"])
                    )
            extension_records.append(
                {
                    "extension_index_zero_based": extension_index,
                    "isotropic_residue_binary": extension["isotropic_residue_binary"],
                    "e10_gram": _integer_rows(extension["gram"]),
                    "cm_cone": {key: value for key, value in cm.items() if key != "rays"},
                    "cm_literal_vinberg": cm_vinberg,
                    "cm_gate_pass": cm_pass,
                    "full_cusp_lift": full_record,
                    "full_D_1_split_gate_pass": full_pass,
                }
            )
        graph = diagram_class["graph"]
        known_matches = [label for label, known_graph in known.items() if _weighted_isomorphic(graph, known_graph)]
        class_records.append(
            {
                "class_index_zero_based": class_index,
                "labelled_member_count": len(diagram_class["members"]),
                "parallel_pair_count": sum(
                    root_gram[left, right] == 2
                    for left in range(12)
                    for right in range(left + 1, 12)
                ),
                "known_root_diagram_matches": known_matches,
                "negative_cartan": _integer_rows(-root_gram),
                "radical_quotient": {
                    "determinant": quotient["determinant"],
                    "smith_diagonal_absolute": quotient[
                        "smith_diagonal_absolute"
                    ],
                    "even_unimodular_extension_count": len(extensions),
                },
                "extension_branches": extension_records,
            }
        )
    full_classes: list[dict[str, Any]] = []
    for class_index, extension_index, wall_gram in actual_full:
        for record in full_classes:
            if _full_isomorphic(wall_gram, record["wall_gram_matrix"]):
                record["members"].append([class_index, extension_index])
                break
        else:
            full_classes.append(
                {
                    "members": [[class_index, extension_index]],
                    "wall_gram_matrix": wall_gram,
                }
            )
    promotion = build_raw_promotion()
    archived_positive = {
        str(case["raw_id"]): bool(case["all_checks_pass"])
        for case in promotion["cases"]
    }
    positive_match_rows = []
    known_full_payload = json.loads(SOURCE_FULL.read_text(encoding="utf-8"))
    for full_index, full_class in enumerate(full_classes):
        matches = []
        for case_key, case in known_full_payload["cases"].items():
            if _full_isomorphic(full_class["wall_gram_matrix"], _matrix(case["negative_wall_gram"])):
                matches.append(f"Type {case_key}")
        positive_match_rows.append(
            {
                "full_isomorphism_class_zero_based": full_index,
                "source_branches": full_class["members"],
                "known_full_chamber_matches": matches,
                "full_split_K3_compiler_pass": (
                    True if matches == ["Type I"] and archived_positive.get("raw_22") else None
                ),
                "compiler_scope": (
                    "inherited exact raw_22 replay after literal full-wall-Gram isomorphism"
                    if matches == ["Type I"]
                    else "generic executable compiler unavailable; no pass inferred"
                ),
            }
        )
    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "scope": {
            "vertices": 12,
            "connected_ordinary_support": True,
            "ordinary_support_edge_count": 13,
            "ordinary_support_maximum_degree": 3,
            "parallel_pair_counts": [2, 3],
            "pairing_alphabet": [0, 1, 2],
            "saturated_radical_quotient_determinant": -4,
            "root_span_index_in_an_E10_extension": 2,
            "not_the_full_nodot_denominator": True,
            "type_II_excluded_by_scope_because_it_has_15_simple_edges": True,
        },
        "scout_denominator": {
            "underlying_unlabelled_graph_count": scout["underlying_graph_count"],
            "exact_lorentzian_corank_two_count": scout[
                "exact_labeled_survivor_count"
            ],
            "determinant_four_labelled_count": len(raw),
            "determinant_four_weighted_isomorphism_class_count": len(classes),
        },
        "class_records": class_records,
        "actual_full_D_1_split_branch_count": len(actual_full),
        "actual_full_D_1_split_isomorphism_class_count": len(full_classes),
        "full_survivor_ledger": positive_match_rows,
        "claims": {
            "strict_L4_for_D_1_split_nodot": False,
            "reason": "the certified 13-edge lane is a proper slice and excludes Type II",
            "full_split_compiler_only_inferred_after_literal_known_gram_match": True,
        },
    }
    payload["all_checks_pass"] = (
        len(raw) == 38
        and sum(record["labelled_member_count"] for record in class_records) == 38
        and len(positive_match_rows) > 0
        and all(
            row["full_split_K3_compiler_pass"] is True
            for row in positive_match_rows
        )
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.cdd)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive mismatch")
        print(sha256(rendered.encode()).hexdigest())
    elif args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
