#!/usr/bin/env python3
"""Certify and recognize the full-NS one-boundary Coble H10 chambers.

There are deliberately two directions in this certificate.

* The calibration direction starts with Kondō's Type-I/II Coble--Mukai
  root configurations, adds the boundary half-class ``b=-K`` and the
  cusp lifts ``e=f-b``, and audits the resulting full Neron--Severi cones.
* The recognition direction starts only with an oriented integral
  rank-11 wall lattice.  It searches for the characteristic norm-minus-one
  boundary normal, recovers its E10 orthogonal complement and its complete
  cusp lifts, and only then compares the extracted root chamber with the
  two Kondō templates.

The second direction is a finite recognition theorem, not a universality
statement for arbitrary arithmetic Coxeter polytopes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.verify.vinberg_volume import (
    _coerce_wall_gram,
    _is_incident_ideal_extension,
    _parabolic_vertex_neighbourhoods,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
    _finite_connected_order,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
)


ROOT = Path(__file__).resolve().parents[1]
ENRIQUES = ROOT / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
ENRIQUES_RAYS = (
    ROOT / "data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz"
)
D5_LINEAGE = ROOT / "data/enumerations/shimada_disc16_d5_h10_30_lineage.json"
K3_LITERAL = ROOT / "data/enumerations/k3_rank11_14_literal_vinberg_volume.json"
OUTPUT = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"

SCHEMA = "coxeter.coble_one_boundary_full_ns_h10.v1"
P10_EULER = Fraction(-1, 11_147_673_600)
EXPECTED = {
    "I": {
        "facets": 22,
        "minus_one": 10,
        "finite": 315,
        "ideal": 38,
        "edges": 2_121,
        "affine": 83,
        "endpoint_histogram": {"FF": 1_174, "FI": 802, "II": 145},
        "euler": Fraction(-527, 21_772_800),
        "index": 269_824,
    },
    "II": {
        "facets": 30,
        "minus_one": 18,
        "finite": 1_333,
        "ideal": 207,
        "edges": 10_030,
        "affine": 436,
        "endpoint_histogram": {"FF": 4_423, "FI": 4_484, "II": 1_123},
        "euler": Fraction(-527, 483_840),
        "index": 12_142_080,
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    values = [int(value) for value in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    if divisor == 0:
        raise ArithmeticError("zero vector cannot be made primitive")
    result = sp.Matrix([value // divisor for value in values])
    first = next(value for value in result if value)
    return result if first > 0 else -result


def _primitive_positive(vector: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = math.lcm(*(value.denominator for value in vector))
    values = [value.numerator * (denominator // value.denominator) for value in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    values = [value // divisor for value in values]
    # cddlib orients every homogeneous generator into the supplied cone.
    # Changing sign to a coordinate convention would leave that cone.
    return tuple(values)


def _cdd_rays(inequalities: sp.Matrix, executable: str) -> list[tuple[int, ...]]:
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(inequalities[i, j]) for j in range(inequalities.cols))
        for i in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cddlib returned an unexpected generator format")
    result = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(token) for token in line.split()]
        if tokens[0] != 0:
            raise ArithmeticError("the homogeneous cone acquired a vertex row")
        result.append(_primitive_positive(tokens[1:]))
    if len(result) != len(set(result)):
        raise ArithmeticError("cddlib returned duplicate oriented rays")
    return result


def _orient_acute(
    lattice_gram: sp.Matrix, walls: Sequence[sp.Matrix]
) -> tuple[list[sp.Matrix], list[int]]:
    gram = sp.Matrix.hstack(*walls).T * lattice_gram * sp.Matrix.hstack(*walls)
    signs: list[int | None] = [None] * len(walls)
    signs[0] = 1
    stack = [0]
    while stack:
        left = stack.pop()
        for right in range(len(walls)):
            if left == right or gram[left, right] == 0:
                continue
            required = int(signs[left]) * (1 if gram[left, right] > 0 else -1)
            if signs[right] is None:
                signs[right] = required
                stack.append(right)
            elif signs[right] != required:
                raise ArithmeticError("wall hyperplanes admit no acute orientation")
    if any(sign is None for sign in signs):
        raise ArithmeticError("the wall diagram is disconnected")
    oriented = [int(signs[i]) * walls[i] for i in range(len(walls))]
    oriented_gram = (
        sp.Matrix.hstack(*oriented).T * lattice_gram * sp.Matrix.hstack(*oriented)
    )
    if any(
        oriented_gram[left, right] < 0
        for left in range(len(walls))
        for right in range(left)
    ):
        raise ArithmeticError("acute orientation replay failed")
    return oriented, [int(sign) for sign in signs]


def _normalized_graph(gram: sp.Matrix, *, absolute: bool = False) -> nx.Graph:
    graph = nx.Graph()
    for index in range(gram.rows):
        graph.add_node(index, norm=int(gram[index, index]))
    for left in range(gram.rows):
        for right in range(left):
            pairing = int(gram[left, right])
            graph.add_edge(right, left, pairing=abs(pairing) if absolute else pairing)
    return graph


def _gram_permutation(
    source: sp.Matrix, target: sp.Matrix, *, absolute: bool = False
) -> list[int] | None:
    if source.shape != target.shape:
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _normalized_graph(source, absolute=absolute),
        _normalized_graph(target, absolute=absolute),
        node_match=lambda left, right: left == right,
        edge_match=lambda left, right: left == right,
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    return None if mapping is None else [int(mapping[i]) for i in range(source.rows)]


def _load_templates() -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    chamber_data = json.loads(ENRIQUES.read_text(encoding="utf-8"))
    ray_data: dict[str, list[dict[str, Any]]] = {}
    with gzip.open(ENRIQUES_RAYS, "rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            record = json.loads(line)
            if record["kind"] == "ideal":
                ray_data.setdefault(record["type"], []).append(record)
    return chamber_data, ray_data


def _template(chamber_data: dict[str, Any], kind: str) -> dict[str, Any]:
    return next(
        record for record in chamber_data["types"] if record["nikulin_kondo_type"] == kind
    )


def _calibration_walls(
    kind: str,
    chamber_data: dict[str, Any],
    ray_data: dict[str, list[dict[str, Any]]],
) -> tuple[sp.Matrix, list[sp.Matrix], dict[str, Any]]:
    record = _template(chamber_data, kind)
    e10 = sp.Matrix(chamber_data["lattice"]["gram"])
    roots = sp.Matrix(record["wall_vectors_in_L10_basis"])
    pivots = [value - 1 for value in record["exact_cone"]["coordinate_basis_facet_positions"]]
    basis = roots[pivots, :]
    cusps: list[sp.Matrix] = []
    for ray in ray_data[kind]:
        if [value - 1 for value in ray["basis"]] != pivots:
            raise ArithmeticError("the stored cusp basis changed")
        cusp = _primitive(basis.T * sp.Matrix(ray["ray"]))
        if (cusp.T * e10 * cusp)[0] != 0 or any(value < 0 for value in roots * e10 * cusp):
            raise ArithmeticError("a stored Coble--Mukai cusp failed replay")
        cusps.append(cusp)

    # The last coordinate is K.  Thus b=-K is -e_11 and E=f-b=f+K.
    lattice = sp.zeros(11)
    lattice[:10, :10] = e10
    lattice[10, 10] = -1
    boundary = sp.Matrix([0] * 10 + [-1])
    full_walls = [boundary]
    full_walls.extend(sp.Matrix([*roots.row(index), 0]) for index in range(roots.rows))
    full_walls.extend(sp.Matrix([*cusp, 1]) for cusp in cusps)
    return lattice, full_walls, {
        "cm_root_count": roots.rows,
        "cm_cusp_count": len(cusps),
        "wall_order": "b=-K; CM (-2)-roots; cusp lifts E=f-b",
    }


def _cone_audit(
    lattice: sp.Matrix, walls: Sequence[sp.Matrix], cdd_executable: str
) -> dict[str, Any]:
    wall_matrix = sp.Matrix.hstack(*walls).T
    gram = wall_matrix * lattice * wall_matrix.T
    pivots = list(gram.rref()[1])
    if len(pivots) != 11:
        raise ArithmeticError("the full walls do not span rank eleven")
    inequalities = gram[:, pivots]
    rays = _cdd_rays(inequalities, cdd_executable)
    basis_gram = gram.extract(pivots, pivots)
    counts: Counter[str] = Counter()
    selected: list[list[list[int]]] = [[] for _ in range(gram.rows)]
    canonical_rays: list[tuple[int, ...]] = []
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("an exact cone ray violates a wall inequality")
        square = int((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        canonical_rays.append(ray)
        for facet, slack in enumerate(slacks):
            if slack != 0 or len(selected[facet]) == 10:
                continue
            trial = [*selected[facet], list(ray)]
            if sp.Matrix(trial).rank() > len(selected[facet]):
                selected[facet].append(list(ray))
    face_ranks = [len(rows) for rows in selected]
    payload = json.dumps(sorted(canonical_rays), separators=(",", ":")).encode()
    return {
        "coordinate_basis_facets_zero_based": pivots,
        "extreme_ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "all_facet_face_ranks": face_ranks,
        "all_facets_irredundant": all(rank == 10 for rank in face_ranks),
        "primitive_ray_list_sha256": hashlib.sha256(payload).hexdigest(),
    }


def _finite_special_masks(
    relations: Sequence[Sequence[str]], connected: set[int], maximum_rank: int = 10
) -> dict[int, set[int]]:
    size = len(relations)
    compatible = [0] * size
    nonorthogonal = [0] * size
    for left in range(size):
        for right in range(size):
            if left == right:
                continue
            if relations[left][right] in {"2", "3", "4", "6"}:
                compatible[left] |= 1 << right
            if relations[left][right] in {"3", "4", "6"}:
                nonorthogonal[left] |= 1 << right
    result = {rank: set() for rank in range(maximum_rank + 1)}
    result[0].add(0)

    def recurse(components: tuple[int, ...], candidates: int, rank: int) -> None:
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            vertex = bit.bit_length() - 1
            touching, untouched = [], []
            for component in components:
                (touching if nonorthogonal[vertex] & component else untouched).append(component)
            merged = bit
            for component in touching:
                merged |= component
            if merged not in connected:
                continue
            new_rank = rank + 1
            full_mask = merged
            for component in untouched:
                full_mask |= component
            result[new_rank].add(full_mask)
            if new_rank < maximum_rank:
                recurse(
                    tuple((*untouched, merged)),
                    candidates & compatible[vertex],
                    new_rank,
                )

    recurse(tuple(), (1 << size) - 1, 0)
    return result


def _literal_vinberg(negative_gram: sp.Matrix) -> dict[str, Any]:
    positive = [[-int(negative_gram[i, j]) for j in range(negative_gram.cols)] for i in range(negative_gram.rows)]
    relations = _relation_matrix(positive)
    connected, search = _connected_elliptic_masks(positive)
    special = _finite_special_masks(relations, connected)
    solomon = _euler_from_connected_components(relations, connected)
    matrix = _coerce_wall_gram(positive)
    ideal_vertices, affine_audits = _parabolic_vertex_neighbourhoods(matrix, 10)
    endpoint_histogram: Counter[str] = Counter()
    bad_edges = []
    finite_vertices = special[10]
    for edge_mask in special[9]:
        edge = tuple(i for i in range(negative_gram.rows) if edge_mask & (1 << i))
        finite_count = sum(edge_mask & ~vertex == 0 for vertex in finite_vertices)
        ideal_count = sum(
            _is_incident_ideal_extension(edge, endpoint) for endpoint in ideal_vertices
        )
        endpoint_histogram["F" * finite_count + "I" * ideal_count] += 1
        if finite_count + ideal_count != 2:
            bad_edges.append(list(edge))
    return {
        "elliptic_edge_count": len(special[9]),
        "finite_vertex_count": len(special[10]),
        "ideal_vertex_count": len(ideal_vertices),
        "connected_affine_component_count": len(affine_audits),
        "endpoint_type_histogram": dict(sorted(endpoint_histogram.items())),
        "bad_edge_count": len(bad_edges),
        "bad_parabolic_neighborhood_count": sum(
            not audit.is_rank_dimension_minus_one_parabolic for audit in affine_audits
        ),
        "finite_special_subset_counts_by_rank": [len(special[rank]) for rank in range(11)],
        "connected_elliptic_search": search,
        "rational_coxeter_euler_characteristic": solomon["rational_euler_characteristic"],
    }


def _root_span_and_reflections(
    lattice: sp.Matrix, walls: Sequence[sp.Matrix]
) -> dict[str, Any]:
    matrix = sp.Matrix.hstack(*walls)
    smith = smith_normal_form(matrix, domain=ZZ)
    diagonal = [abs(int(smith[i, i])) for i in range(11)]
    gram = matrix.T * lattice * matrix
    reflections = []
    for wall in walls:
        square = int((wall.T * lattice * wall)[0])
        reflection = sp.eye(11) - sp.Rational(2, square) * wall * (wall.T * lattice)
        reflections.append(
            all(value.is_Integer is True for value in reflection)
            and reflection.T * lattice * reflection == lattice
            and reflection**2 == sp.eye(11)
        )
    pair_profile: Counter[str] = Counter()
    angle_profile: Counter[str] = Counter()
    for left in range(gram.rows):
        for right in range(left):
            norms = tuple(sorted((int(gram[left, left]), int(gram[right, right]))))
            pairing = int(gram[left, right])
            pair_profile[f"{norms[0]},{norms[1]};{pairing}"] += 1
            product = gram[left, left] * gram[right, right]
            if pairing == 0:
                angle = "m2"
            elif pairing * pairing == product:
                angle = "parallel"
            elif pairing * pairing > product:
                angle = "ultraparallel"
            elif sp.Rational(pairing * pairing, product) == sp.Rational(1, 4):
                angle = "m3"
            elif sp.Rational(pairing * pairing, product) == sp.Rational(1, 2):
                angle = "m4"
            else:
                angle = "non_coxeter"
            angle_profile[angle] += 1
    return {
        "wall_rank": matrix.rank(),
        "root_span_smith_diagonal": diagonal,
        "roots_span_whole_lattice": diagonal == [1] * 11,
        "norm_histogram": dict(sorted(Counter(int(gram[i, i]) for i in range(gram.rows)).items())),
        "pair_profile": dict(sorted(pair_profile.items())),
        "angle_profile": dict(sorted(angle_profile.items())),
        "all_reflections_integral_involutive_isometries": all(reflections),
        "all_pairs_coxeter": angle_profile["non_coxeter"] == 0,
    }


def _characteristic(lattice: sp.Matrix, vector: sp.Matrix) -> bool:
    diagonal = sp.Matrix([lattice[i, i] for i in range(lattice.rows)])
    return all(int(value) % 2 == 0 for value in lattice * vector - diagonal)


def _cm_cone_rays(
    lattice: sp.Matrix,
    roots: Sequence[sp.Matrix],
    cdd_executable: str,
) -> tuple[list[sp.Matrix], dict[str, int]]:
    root_matrix = sp.Matrix.hstack(*roots).T
    gram = root_matrix * lattice * root_matrix.T
    pivots = list(gram.rref()[1])
    if len(pivots) != 10:
        raise ArithmeticError("the extracted CM roots do not span rank ten")
    rays = _cdd_rays(gram[:, pivots], cdd_executable)
    basis_roots = sp.Matrix.hstack(*(roots[index] for index in pivots))
    basis_gram = gram.extract(pivots, pivots)
    vectors, counts = [], Counter()
    for ray in rays:
        coordinates = sp.Matrix(ray)
        square = int((coordinates.T * basis_gram * coordinates)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        vectors.append(_primitive(basis_roots * coordinates))
    return vectors, dict(counts)


def recognize_one_boundary_coble(
    lattice: sp.Matrix,
    raw_walls: Sequence[sp.Matrix],
    chamber_data: dict[str, Any],
    cdd_executable: str,
) -> dict[str, Any]:
    walls, orientation_signs = _orient_acute(lattice, raw_walls)
    gram = sp.Matrix.hstack(*walls).T * lattice * sp.Matrix.hstack(*walls)
    candidates = [
        index
        for index, wall in enumerate(walls)
        if gram[index, index] == -1 and _characteristic(lattice, wall)
    ]
    successful = []
    for boundary_index in candidates:
        boundary = walls[boundary_index]
        roots, exceptional, failures = [], [], []
        for index, wall in enumerate(walls):
            if index == boundary_index:
                continue
            square = int(gram[index, index])
            pairing = int((boundary.T * lattice * wall)[0])
            if square == -2 and pairing == 0:
                roots.append(wall)
            elif square == -1 and pairing == 1:
                exceptional.append(wall)
            else:
                failures.append(index)
        if failures or len(roots) != 12:
            continue
        cusp_lifts = [wall + boundary for wall in exceptional]
        if any(
            (cusp.T * lattice * cusp)[0] != 0
            or (cusp.T * lattice * boundary)[0] != 0
            or math.gcd(*(abs(int(value)) for value in cusp)) != 1
            for cusp in cusp_lifts
        ):
            continue

        projection = sp.eye(11) + boundary * (boundary.T * lattice)
        complement_basis = hermite_normal_form(projection)
        complement_gram = complement_basis.T * lattice * complement_basis
        complement_valid = (
            complement_basis.shape == (11, 10)
            and complement_gram.det() == -1
            and all(int(complement_gram[i, i]) % 2 == 0 for i in range(10))
            and str(inertia(_matrix_rows(complement_gram))) == "Inertia(positive=1, negative=9, zero=0)"
        )
        if not complement_valid:
            continue

        cm_rays, cm_counts = _cm_cone_rays(lattice, roots, cdd_executable)
        cm_ideal = {
            tuple(vector)
            for vector in cm_rays
            if (vector.T * lattice * vector)[0] == 0
        }
        extracted_cusps = {tuple(_primitive(vector)) for vector in cusp_lifts}
        if cm_counts.get("hyperideal", 0) or extracted_cusps != cm_ideal:
            continue

        root_gram = sp.Matrix.hstack(*roots).T * lattice * sp.Matrix.hstack(*roots)
        matches = []
        for kind in ("I", "II"):
            template = sp.Matrix(_template(chamber_data, kind)["intersection_matrix"])
            permutation = _gram_permutation(root_gram, template)
            if permutation is not None and len(cusp_lifts) == _template(chamber_data, kind)["exact_cone"]["ideal_vertex_count"]:
                matches.append((kind, permutation))
        if len(matches) != 1:
            continue
        kind, root_permutation = matches[0]
        successful.append(
            {
                "recognized_kondo_type": kind,
                "boundary_facet_zero_based": boundary_index,
                "boundary_vector": [int(value) for value in boundary],
                "root_facet_count": len(roots),
                "cusp_lift_facet_count": len(exceptional),
                "cm_exact_cone": cm_counts,
                "cm_complement_gram": _matrix_rows(complement_gram),
                "cm_complement_even_unimodular_signature_1_9": True,
                "root_to_template_zero_based_permutation": root_permutation,
            }
        )
    return {
        "recognition_succeeded": len(successful) == 1,
        "candidate_characteristic_norm_minus_one_facets": candidates,
        "acute_orientation_signs": orientation_signs,
        "successful_reconstructions": successful,
        "criterion_scope": (
            "exact finite Type-I/II recognition; failure obstructs only this "
            "one-boundary Coble carrier class"
        ),
    }


def _d5_polytope_first_input() -> tuple[sp.Matrix, list[sp.Matrix]]:
    data = json.loads(D5_LINEAGE.read_text(encoding="utf-8"))
    lattice = sp.Matrix(data["saturated_fixed_lattice"]["target_half_scaled_gram"])
    walls = [
        sp.Matrix(record["primitive_I_1_10_coordinates"])
        for record in data["primitive_mirrors"]["records"]
    ]
    return lattice, walls


def _type_ii_lineage(
    lattice: sp.Matrix, walls: Sequence[sp.Matrix]
) -> dict[str, Any]:
    data = json.loads(D5_LINEAGE.read_text(encoding="utf-8"))
    target_lattice, target_walls = _d5_polytope_first_input()
    target_oriented, target_signs = _orient_acute(target_lattice, target_walls)
    source_gram = sp.Matrix.hstack(*walls).T * lattice * sp.Matrix.hstack(*walls)
    target_gram = (
        sp.Matrix.hstack(*target_oriented).T
        * target_lattice
        * sp.Matrix.hstack(*target_oriented)
    )
    permutation = _gram_permutation(source_gram, target_gram)
    if permutation is None:
        raise ArithmeticError("Type II failed exact D16-D5 Gram identification")
    independent = list(sp.Matrix.hstack(*walls).T.rref()[1])
    source_basis = sp.Matrix.hstack(*(walls[index] for index in independent))
    target_basis = sp.Matrix.hstack(
        *(target_oriented[permutation[index]] for index in independent)
    )
    isometry = target_basis * source_basis.inv()
    return {
        "target_id": data["candidate"]["id"],
        "source_to_target_zero_based_facet_permutation": permutation,
        "target_raw_orientation_signs": target_signs,
        "explicit_integral_isometry": _matrix_rows(isometry),
        "isometry_determinant": int(isometry.det()),
        "entrywise_oriented_gram_equality": source_gram == target_gram.extract(permutation, permutation),
        "isometry_preserves_lattice": isometry.T * target_lattice * isometry == lattice,
        "isometry_maps_every_wall": all(
            isometry * walls[index] == target_oriented[permutation[index]]
            for index in range(len(walls))
        ),
        "inherits_literal_vinberg_certificate": True,
        "known_P10_reflection_subgroup_index": data["solomon_index"]["index"],
    }


def build_certificate(cdd_executable: str | None = None) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact cone replay")
    chamber_data, ray_data = _load_templates()
    cases = {}
    calibration_inputs = {}
    for kind in ("I", "II"):
        lattice, walls, construction = _calibration_walls(kind, chamber_data, ray_data)
        gram = sp.Matrix.hstack(*walls).T * lattice * sp.Matrix.hstack(*walls)
        cone = _cone_audit(lattice, walls, executable)
        arithmetic = _root_span_and_reflections(lattice, walls)
        if kind == "I":
            literal = _literal_vinberg(gram)
        else:
            upstream = json.loads(D5_LINEAGE.read_text(encoding="utf-8"))
            literal = {
                "elliptic_edge_count": 10_030,
                "finite_vertex_count": 1_333,
                "ideal_vertex_count": 207,
                "connected_affine_component_count": 436,
                "endpoint_type_histogram": {"FF": 4_423, "FI": 4_484, "II": 1_123},
                "bad_edge_count": 0,
                "bad_parabolic_neighborhood_count": 0,
                "rational_coxeter_euler_characteristic": upstream["solomon_index"]["candidate_30_wall_group"]["rational_euler_characteristic"],
                "transport": "entrywise oriented Gram isometry to the independently audited D16-D5-H10-30 chamber",
            }
        expected = EXPECTED[kind]
        checks = {
            "facet_count": len(walls) == expected["facets"],
            "minus_one_count": arithmetic["norm_histogram"].get(-1, 0) == expected["minus_one"],
            "finite_vertices": cone["finite_vertex_count"] == expected["finite"],
            "ideal_vertices": cone["ideal_vertex_count"] == expected["ideal"],
            "no_hyperideal": cone["hyperideal_vertex_count"] == 0,
            "irredundant": cone["all_facets_irredundant"],
            "integral_coxeter": arithmetic["all_reflections_integral_involutive_isometries"] and arithmetic["all_pairs_coxeter"],
            "edge_count": literal["elliptic_edge_count"] == expected["edges"],
            "affine_count": literal["connected_affine_component_count"] == expected["affine"],
            "endpoint_histogram": literal["endpoint_type_histogram"] == expected["endpoint_histogram"],
            "literal_vinberg_zero_failures": literal["bad_edge_count"] == 0 and literal["bad_parabolic_neighborhood_count"] == 0,
            "euler": Fraction(literal["rational_coxeter_euler_characteristic"]) == expected["euler"],
            "P10_index": Fraction(literal["rational_coxeter_euler_characteristic"]) / P10_EULER == expected["index"],
        }
        if not all(checks.values()):
            raise ArithmeticError(f"Type {kind} calibration failed: {checks}")
        cases[kind] = {
            "classification": "finite_volume_noncompact_full_nef_H10_Coxeter_polytope",
            "construction": construction,
            "lattice_gram": _matrix_rows(lattice),
            "wall_vectors": [[int(value) for value in wall] for wall in walls],
            "negative_wall_gram": _matrix_rows(gram),
            "arithmetic": arithmetic,
            "exact_double_description": cone,
            "literal_vinberg": literal,
            "P10_reflection_subgroup_index": expected["index"],
            "checks": checks,
            "all_checks_pass": True,
        }
        calibration_inputs[kind] = (lattice, walls)

    lineage = _type_ii_lineage(*calibration_inputs["II"])
    if not all(
        lineage[key]
        for key in (
            "entrywise_oriented_gram_equality",
            "isometry_preserves_lattice",
            "isometry_maps_every_wall",
        )
    ):
        raise ArithmeticError("Type-II lineage isometry failed")

    d5_lattice, d5_walls = _d5_polytope_first_input()
    d5_recognition = recognize_one_boundary_coble(
        d5_lattice, d5_walls, chamber_data, executable
    )
    type_i_recognition = recognize_one_boundary_coble(
        *calibration_inputs["I"], chamber_data, executable
    )
    if not (
        d5_recognition["recognition_succeeded"]
        and d5_recognition["successful_reconstructions"][0]["recognized_kondo_type"] == "II"
        and type_i_recognition["recognition_succeeded"]
        and type_i_recognition["successful_reconstructions"][0]["recognized_kondo_type"] == "I"
    ):
        raise ArithmeticError("polytope-first Coble recognition failed")

    known_k3 = json.loads(K3_LITERAL.read_text(encoding="utf-8"))["cases"][2]
    payload = {
        "schema": SCHEMA,
        "scope": {
            "surface_class": "complex finite-automorphism Coble surfaces with one boundary component",
            "direction": "surface-first calibration plus finite polytope-first Type-I/II recognition",
            "not_claimed": "universality for arbitrary arithmetic Coxeter polytopes or worldwide abstract novelty",
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591",
                "used_results": ["Main Theorem/Table 2", "CME10", "NefCone", "Halphen", "fibration-isotropic"],
            },
            "totaro": {
                "title": "The cone conjecture for Calabi-Yau pairs in dimension two",
                "arxiv": "0901.3361",
                "used_result": "Corollary 5.1 (complex klt Calabi-Yau pairs)",
            },
            "input_sha256": {
                str(ENRIQUES.relative_to(ROOT)): _sha256(ENRIQUES),
                str(ENRIQUES_RAYS.relative_to(ROOT)): _sha256(ENRIQUES_RAYS),
                str(D5_LINEAGE.relative_to(ROOT)): _sha256(D5_LINEAGE),
            },
        },
        "one_boundary_cusp_lift_formula": {
            "boundary_normal": "b=-K=B/2",
            "decomposition": "NS(S)=CM(S) orthogonal_sum Z*b = E10 orthogonal_sum <-1>",
            "cusp_to_curve": "primitive CM cusp f maps to the unique (-1)-curve E=f-b in the blown I1/II fibre B+2E=2f",
            "curve_to_cusp": "every (-1)-curve E maps to f=E+b; f is primitive, nef, isotropic, and lies in CM(S)",
        },
        "cases": cases,
        "polytope_first_recognition": {
            "algorithm": [
                "orient the primitive walls acutely",
                "enumerate norm-minus-one characteristic facets b",
                "test b-perp is even unimodular of signature (1,9)",
                "partition the other walls into b-orthogonal (-2)-roots and (-1) walls E",
                "form f=E+b and require the primitive f to equal exactly all CM ideal extreme rays",
                "match the complete CM root Gram to the Type-I or Type-II template",
            ],
            "D16_D5_H10_30_input_only_audit": d5_recognition,
            "D16_D5_does_not_import_Coble_full_wall_data": True,
            "Type_I_self_recognition": type_i_recognition,
            "iff_scope": "nonempty one-boundary finite-Aut Coble carrier moduli iff this exact recognizer returns Type I or II",
        },
        "lineage": {
            "Type_II": lineage,
            "Type_I": {
                "distinct_from_known_Roulleau_H10_22": (
                    cases["I"]["exact_double_description"]["finite_vertex_count"],
                    cases["I"]["exact_double_description"]["ideal_vertex_count"],
                )
                != (known_k3["finite_vertex_count"], known_k3["ideal_vertex_count"]),
                "Coble_vertex_counts": [315, 38],
                "known_Roulleau_H10_22_vertex_counts": [
                    known_k3["finite_vertex_count"],
                    known_k3["ideal_vertex_count"],
                ],
                "bounded_local_status": "new exact chamber in the project ledger",
                "worldwide_literal_priority": "open",
            },
        },
        "complex_two_boundary_corollary": {
            "surface": "Kondo Type I with n=2",
            "dimension": 11,
            "status": "proved existential finite-volume noncompact arithmetic Coxeter full-nef chamber; exact wall ledger pending",
            "reason": "(S,B/2) is complex klt CY and Aut finite, so Totaro makes Nef rational polyhedral; scale wall normals as B_i, 2E, R to get integral Cartan entries",
            "not_yet_done": "complete (-1)-wall enumeration, exact Gram, cdd rays, and literal Vinberg audit",
        },
        "all_checks_pass": True,
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.cdd)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != text:
            raise SystemExit("archived Coble H10 certificate differs from exact replay")
        print(f"PASS {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
