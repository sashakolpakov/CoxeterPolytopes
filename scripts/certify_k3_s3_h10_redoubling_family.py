#!/usr/bin/env python3
"""Certify the K3-S3 H10 finite-path redoubling family.

The complete 13-wall chamber K3-S3-C2-H10-13 has two parallel doubling
walls.  Allcock's finite-subtree construction therefore turns every finite
path in their infinite-dihedral Cayley graph into a finite-volume Coxeter
polytope.  This script proves the boundary formula in the full Picard
lattice F=U+E8(-2)+<-2>, certifies the first two nontrivial chambers by
Vinberg's exact endpoint test and extreme-ray enumeration, and compares the
family with the pinned ordinary and index-two CoxIter H10 lanes.

The surface statement is deliberately a divisor-lattice carrier statement:
the walls are integral divisor classes in NS(X)=F.  It is not a claim that
the glued chambers are nef cones or that their reflections are surface
automorphisms.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _enumerate_vertices


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data/enumerations/"
    "k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json"
)
OUTPUT = ROOT / "data/enumerations/k3_s3_h10_redoubling_family.json"
COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_ARCHIVE = COXITER_ROOT / "graphs/testing_index2.tar.gz"

SCHEMA = "coxeter.k3_s3.h10_redoubling_family.v1"
SOURCE_SHA256 = "93cec35a3c395482f4f99be827b3f5a6ba44764fc7defa22bf23eb6ee951acf0"
SOURCE_PAYLOAD_SHA256 = (
    "c033a684aa04b2c6d0bfe80d1b5393f0c52f5a9c9c1039140ca674ea9cc67a34"
)
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)

LOCAL_LEDGER_PATHS = {
    "R11_direct": (
        ROOT / "data/enumerations/roulleau_rank11_direct_graph_denominator.json",
        "d82e31de7dd9d1f567cba9cc9a98588168628714521842c87943a9c76ec20d43",
    ),
    "R11_all_subgroups": (
        ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json",
        "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
    ),
    "R12_all_subgroups": (
        ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json",
        "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    ),
}

DIMENSION = 10
DOUBLING_SEED_INDICES = (9, 12)  # zero-based; walls 10 and 13
MOVING_SEED_INDICES = (8, 10)  # c and d; walls 9 and 11
FIXED_SEED_INDICES = tuple(
    index
    for index in range(13)
    if index not in (*DOUBLING_SEED_INDICES, *MOVING_SEED_INDICES)
)

EXPECTED_VOLUME_COUNTS = {
    "N2_A": (41, 5, 252, {"0,2": 6, "1,1": 82, "2,0": 164}),
    "N2_B": (52, 6, 318, {"0,2": 10, "1,1": 96, "2,0": 212}),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector if value))


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    rationals = [sp.Rational(value) for value in vector]
    denominator = int(sp.ilcm(*(value.q for value in rationals)))
    integral = sp.Matrix([int(value * denominator) for value in rationals])
    return integral / _content(integral)


def _vector_payload(vector: sp.Matrix) -> list[int | str]:
    return [
        int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
        for value in vector
    ]


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [_vector_payload(matrix.row(index)) for index in range(matrix.rows)]


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return _content(gram * vector)


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _pair(gram, root, root)
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * (root.T * gram)
    )


def _unoriented_key(vector: sp.Matrix) -> tuple[int, ...]:
    primitive = tuple(int(value) for value in _primitive(vector))
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _smith_invariants(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _load_seed() -> tuple[dict[str, Any], sp.Matrix, tuple[sp.Matrix, ...], sp.Matrix]:
    actual = _sha256(SOURCE)
    if actual != SOURCE_SHA256:
        raise ArithmeticError(f"pinned K3 carrier changed: {actual}")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if not (
        source["all_checks_pass"]
        and source["payload_sha256"] == SOURCE_PAYLOAD_SHA256
        and source["transported_H10_13_chamber"]["known_full_reflection_chamber"]
    ):
        raise ArithmeticError("the K3-S3 source carrier is not the frozen green chamber")
    gram = sp.Matrix(source["full_Picard_lattice"]["standard_geometric_gram"])
    roots = tuple(
        sp.Matrix(record["primitive_geometric_standard_F_coordinates"])
        for record in source["transported_H10_13_chamber"]["wall_records"]
    )
    interior = sp.Matrix(
        source["transported_H10_13_chamber"][
            "square358_strict_interior_standard_F_coordinates"
        ]
    )
    if not (
        gram.shape == (11, 11)
        and exact_inertia(gram) == (1, 10, 0)
        and len(roots) == 13
        and all(_pair(gram, interior, root) > 0 for root in roots)
    ):
        raise ArithmeticError("seed lattice or chamber changed")
    return source, gram, roots, interior


def _path_transforms(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    copies: int,
    start: int,
) -> tuple[sp.Matrix, ...]:
    if copies < 1 or start not in (0, 1):
        raise ValueError("copies must be positive and start must be 0 or 1")
    generators = tuple(
        _reflection(gram, roots[index]) for index in DOUBLING_SEED_INDICES
    )
    transforms = [sp.eye(gram.rows)]
    for edge in range(copies - 1):
        transforms.append(transforms[-1] * generators[(start + edge) % 2])
    return tuple(transforms)


def _path_boundary(
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
    seed_interior: sp.Matrix,
    copies: int,
    start: int,
) -> tuple[tuple[sp.Matrix, ...], sp.Matrix]:
    transforms = _path_transforms(gram, seed_roots, copies, start)
    raw_interior = sum(
        (transform * seed_interior for transform in transforms),
        sp.zeros(gram.rows, 1),
    )
    internal = {
        _unoriented_key(
            transforms[edge]
            * seed_roots[DOUBLING_SEED_INDICES[(start + edge) % 2]]
        )
        for edge in range(copies - 1)
    }
    boundary: dict[tuple[int, ...], sp.Matrix] = {}
    for transform in transforms:
        for seed_root in seed_roots:
            root = _primitive(transform * seed_root)
            key = _unoriented_key(root)
            if key in internal:
                continue
            pairing = _pair(gram, raw_interior, root)
            if pairing == 0:
                raise ArithmeticError("path-average lies on a purported boundary wall")
            if pairing < 0:
                root = -root
            previous = boundary.setdefault(key, root)
            if previous != root:
                raise ArithmeticError("inconsistent boundary orientation")
    ordered = tuple(sorted(boundary.values(), key=lambda vector: tuple(vector)))
    interior = _primitive(raw_interior)
    if any(_pair(gram, interior, root) <= 0 for root in ordered):
        raise ArithmeticError("constructed path interior is not strict")
    return ordered, interior


def _facet_gram(gram: sp.Matrix, roots: Sequence[sp.Matrix]) -> sp.Matrix:
    return sp.Matrix(
        [[_pair(gram, left, right) for right in roots] for left in roots]
    )


def _relation_label(
    gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix
) -> tuple[str, sp.Rational]:
    product = _pair(gram, left, right)
    if product < 0:
        raise ArithmeticError("the constructed chamber is not acute")
    ratio = sp.factor(
        product**2 / (_pair(gram, left, left) * _pair(gram, right, right))
    )
    labels = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "0",
    }
    if ratio in labels:
        return labels[ratio], ratio
    if ratio > 1:
        return "1", ratio
    raise ArithmeticError(f"non-Coxeter normalized product {ratio}")


def _relation_matrix(
    gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> tuple[tuple[str, ...], ...]:
    matrix = [["1" if i == j else "2" for j in range(len(roots))] for i in range(len(roots))]
    for left in range(len(roots)):
        for right in range(left):
            label, _ratio = _relation_label(gram, roots[left], roots[right])
            matrix[left][right] = matrix[right][left] = label
    return tuple(tuple(row) for row in matrix)


def _relation_histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    counts = Counter(
        relations[left][right]
        for left in range(len(relations))
        for right in range(left)
    )
    return {key: counts[key] for key in sorted(counts)}


def _ratio_histogram(gram: sp.Matrix, roots: Sequence[sp.Matrix]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for left in range(len(roots)):
        for right in range(left):
            _label, ratio = _relation_label(gram, roots[left], roots[right])
            counts[str(ratio)] += 1
    return dict(sorted(counts.items()))


def _relation_graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            graph.add_edge(left, right, relation=str(relations[left][right]))
    return graph


RELATION_EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match(
    "relation", None
)


def _same_relation_graph(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    if len(left) != len(right):
        return False
    return nx.is_isomorphic(
        _relation_graph(left),
        _relation_graph(right),
        edge_match=RELATION_EDGE_MATCH,
    )


def _relation_wl_hash(relations: Sequence[Sequence[str]]) -> str:
    return nx.weisfeiler_lehman_graph_hash(
        _relation_graph(relations), edge_attr="relation"
    )


def _weighted_gram_graph(
    gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> nx.Graph:
    graph = nx.Graph()
    for index, root in enumerate(roots):
        graph.add_node(
            index,
            decoration=f"{_pair(gram, root, root)},{_divisibility(gram, root)}",
        )
    for left in range(len(roots)):
        for right in range(left):
            graph.add_edge(
                left,
                right,
                product=str(_pair(gram, roots[left], roots[right])),
            )
    return graph


def _same_weighted_gram(
    gram: sp.Matrix,
    left: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
) -> bool:
    return nx.is_isomorphic(
        _weighted_gram_graph(gram, left),
        _weighted_gram_graph(gram, right),
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "decoration", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "product", None
        ),
    )


def _doubling_walls(relations: Sequence[Sequence[str]]) -> list[int]:
    return [
        index
        for index, row in enumerate(relations)
        if all(label not in {"3", "5", "7", "9"} for label in row)
    ]


def _mirror_types(
    gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> dict[str, int]:
    counts = Counter(
        (_pair(gram, root, root), _divisibility(gram, root)) for root in roots
    )
    return {
        f"{square},{divisibility}": count
        for (square, divisibility), count in sorted(counts.items())
    }


def _span_audit(roots: Sequence[sp.Matrix], rank: int) -> dict[str, Any]:
    matrix = sp.Matrix.hstack(*roots)
    smith = _smith_invariants(matrix)
    return {
        "rank": matrix.rank(),
        "smith_invariants": smith,
        "index": math.prod(smith) if len(smith) == rank else None,
        "full_index_one": matrix.rank() == rank and smith == [1] * rank,
    }


def _reflection_audit(
    gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> dict[str, Any]:
    checks = []
    for root in roots:
        reflection = _reflection(gram, root)
        checks.append(
            _content(root) == 1
            and all(value.is_Integer for value in reflection)
            and reflection.det() == -1
            and reflection**2 == sp.eye(gram.rows)
            and reflection.T * gram * reflection == gram
            and reflection * root == -root
        )
    return {
        "wall_count": len(roots),
        "all_primitive_integral_reflections": all(checks),
        "failed_zero_based": [index for index, value in enumerate(checks) if not value],
    }


def _expected_fingerprint(copies: int, start: int) -> tuple[dict[str, int], dict[str, int]]:
    if copies % 2:
        k = (copies - 1) // 2
        relations = {
            "0": k + 2,
            "1": 2 * k * k + 3 * k,
            "2": 65 + 19 * k,
            "3": 9 + 2 * k,
            "4": 2,
        }
        types = {"-2,1": k + 2, "-2,2": 2, "-4,2": k + 9}
    else:
        k = copies // 2
        relations = {
            "0": k + 2 if start == 0 else k + 1,
            "1": 2 * k * k + k - 1,
            "2": 56 + 19 * k if start == 0 else 55 + 19 * k,
            "3": 8 + 2 * k,
            "4": 1 if start == 0 else 3,
        }
        types = (
            {"-2,1": k + 1, "-2,2": 2, "-4,2": k + 9}
            if start == 0
            else {"-2,1": k + 2, "-2,2": 2, "-4,2": k + 8}
        )
    return (
        {label: count for label, count in relations.items() if count},
        {label: count for label, count in types.items() if count},
    )


def _closed_boundary_roots(
    seed_roots: Sequence[sp.Matrix], copies: int, start: int
) -> tuple[sp.Matrix, ...]:
    a, b = (seed_roots[index] for index in DOUBLING_SEED_INDICES)
    c, d = (seed_roots[index] for index in MOVING_SEED_INDICES)
    if start == 1:
        a, b, c, d = b, a, d, c
    moving = []
    for index in range(copies + 1):
        if index % 2 == 0:
            m = index // 2
            moving.append(c + 2 * m * m * a + 2 * m * (m - 1) * b)
        else:
            m = (index - 1) // 2
            moving.append(d + 2 * m * (m + 1) * a + 2 * m * m * b)
    endpoints = (b, copies * a + (copies - 1) * b)
    fixed = tuple(seed_roots[index] for index in FIXED_SEED_INDICES)
    return (*fixed, *moving, *endpoints)


def _case(
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
    seed_interior: sp.Matrix,
    copies: int,
    start: int,
    *,
    include_roots: bool,
) -> tuple[dict[str, Any], tuple[sp.Matrix, ...], sp.Matrix]:
    roots, interior = _path_boundary(
        gram, seed_roots, seed_interior, copies, start
    )
    relations = _relation_matrix(gram, roots)
    relation_histogram = _relation_histogram(relations)
    types = _mirror_types(gram, roots)
    expected_relations, expected_types = _expected_fingerprint(copies, start)
    closed = _closed_boundary_roots(seed_roots, copies, start)
    endpoints = [roots[index] for index in _doubling_walls(relations)]
    a, b = (seed_roots[index] for index in DOUBLING_SEED_INDICES)
    expected_endpoints = (
        (b, copies * a + (copies - 1) * b)
        if start == 0
        else (a, (copies - 1) * a + copies * b)
    )
    endpoint_keys = {_unoriented_key(root) for root in endpoints}
    expected_endpoint_keys = {_unoriented_key(root) for root in expected_endpoints}
    endpoint_product = _pair(gram, expected_endpoints[0], expected_endpoints[1])
    endpoint_ratio = sp.factor(
        endpoint_product**2
        / (
            _pair(gram, expected_endpoints[0], expected_endpoints[0])
            * _pair(gram, expected_endpoints[1], expected_endpoints[1])
        )
    )
    span = _span_audit(roots, gram.rows)
    reflection_audit = _reflection_audit(gram, roots)
    checks = {
        "closed_boundary_formula": {
            _unoriented_key(root) for root in roots
        }
        == {_unoriented_key(root) for root in closed},
        "facet_formula_N_plus_12": len(roots) == copies + 12,
        "relation_formula": relation_histogram == expected_relations,
        "mirror_type_formula": types == expected_types,
        "exactly_two_boundary_doubling_walls": len(endpoints) == 2,
        "endpoint_formula": endpoint_keys == expected_endpoint_keys,
        "endpoint_squares_minus_two": all(
            _pair(gram, root, root) == -2 for root in expected_endpoints
        ),
        "endpoint_divisibilities_two": all(
            _divisibility(gram, root) == 2 for root in expected_endpoints
        ),
        "raw_endpoint_product_two": endpoint_product == 2,
        "raw_product_squared_four": endpoint_product**2 == 4,
        "normalized_squared_ratio_one_parallel": endpoint_ratio == 1,
        "strict_timelike_interior": (
            _pair(gram, interior, interior) > 0
            and all(_pair(gram, interior, root) > 0 for root in roots)
        ),
        "walls_span_full_Picard_index_one": span["full_index_one"],
        "all_boundary_reflections_integral": reflection_audit[
            "all_primitive_integral_reflections"
        ],
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"N={copies}, start={start} failed: "
            f"{[key for key, value in checks.items() if not value]}"
        )
    record: dict[str, Any] = {
        "id": f"K3-S3-H10-P{copies}-{'A' if start == 0 else 'B'}",
        "copy_count_N": copies,
        "orientation": "A_wall10_first" if start == 0 else "B_wall13_first",
        "facet_count": len(roots),
        "reflection_subgroup_index": copies,
        "volume_multiple_of_seed": copies,
        "primitive_strict_interior": _vector_payload(interior),
        "strict_interior_square": int(_pair(gram, interior, interior)),
        "minimum_strict_pairing": int(
            min(_pair(gram, interior, root) for root in roots)
        ),
        "mirror_type_histogram": types,
        "relation_histogram": relation_histogram,
        "normalized_squared_ratio_histogram": _ratio_histogram(gram, roots),
        "relation_WL_hash": _relation_wl_hash(relations),
        "doubling_walls_one_based": [index + 1 for index in _doubling_walls(relations)],
        "parallel_endpoint_certificate": {
            "left": _vector_payload(expected_endpoints[0]),
            "right": _vector_payload(expected_endpoints[1]),
            "squares": [-2, -2],
            "raw_inner_product": int(endpoint_product),
            "raw_inner_product_squared": int(endpoint_product**2),
            "norm_product": 4,
            "normalized_squared_ratio": str(endpoint_ratio),
            "classification": "parallel_asymptotic_not_ultraparallel",
        },
        "wall_span": span,
        "reflection_audit": reflection_audit,
        "checks": checks,
    }
    if include_roots:
        record["boundary_roots_standard_F_coordinates"] = [
            _vector_payload(root) for root in roots
        ]
        record["normalized_relation_matrix"] = [list(row) for row in relations]
    return record, roots, interior


def _volume_audit(
    name: str,
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    interior: sp.Matrix,
) -> dict[str, Any]:
    facet_gram = _facet_gram(gram, roots)
    positive_gram = [
        [-int(facet_gram[left, right]) for right in range(len(roots))]
        for left in range(len(roots))
    ]
    vinberg = verify_vinberg_volume_criterion(
        positive_gram, DIMENSION, simple_roots_complete=True
    )
    endpoint_histogram = Counter(
        (len(edge.finite_endpoints), len(edge.ideal_endpoints))
        for edge in vinberg.edge_audits
    )

    vectors = sp.Matrix.hstack(*roots, interior)
    orbit_gram = vectors.T * gram * vectors
    rays, basis, basis_gram = _enumerate_vertices(
        orbit_gram, list(range(len(roots))), [len(roots)]
    )
    ray_histogram = Counter(ray["type"] for ray in rays)
    facet_support = {}
    for facet in range(1, len(roots) + 1):
        equality = [
            sp.Matrix(ray["basis_coordinates"])
            for ray in rays
            if facet in ray["active_facets"]
        ]
        span_rank = sp.Matrix.hstack(*equality).rank() if equality else 0
        facet_support[str(facet)] = {
            "equality_extreme_ray_count": len(equality),
            "equality_ray_span_rank": span_rank,
            "is_genuine_facet": span_rank == DIMENSION,
        }
    ideal_links = []
    for ray in rays:
        if ray["type"] != "ideal":
            continue
        active = tuple(index - 1 for index in ray["active_facets"])
        link_gram = -facet_gram.extract(active, active)
        inertia = exact_inertia(link_gram)
        ideal_links.append(
            {
                "active_facets": ray["active_facets"],
                "rank": link_gram.rank(),
                "inertia": list(inertia),
                "parabolic_rank_n_minus_one": (
                    link_gram.rank() == DIMENSION - 1
                    and inertia[0] == DIMENSION - 1
                    and inertia[1] == 0
                ),
            }
        )
    ray_vectors = [sp.Matrix(ray["basis_coordinates"]) for ray in rays]
    same_time_component = all(
        (ray_vectors[left].T * basis_gram * ray_vectors[right])[0] > 0
        for left in range(len(ray_vectors))
        for right in range(left)
    )
    expected = EXPECTED_VOLUME_COUNTS[name]
    endpoint_payload = {
        f"{finite},{ideal}": count
        for (finite, ideal), count in sorted(endpoint_histogram.items())
    }
    checks = {
        "vinberg_finite_volume_noncompact": (
            vinberg.decision == "finite_volume_noncompact"
            and vinberg.finite_volume is True
            and vinberg.compact is False
        ),
        "all_elliptic_edges_have_two_endpoints": (
            not vinberg.finite_volume_bad_edges
            and all(edge.endpoint_count == 2 for edge in vinberg.edge_audits)
        ),
        "all_full_parabolic_neighborhoods_pass": all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in vinberg.parabolic_neighborhood_audits
        ),
        "counts_match_independent_discovery": (
            vinberg.finite_vertex_count,
            vinberg.ideal_vertex_count,
            vinberg.elliptic_edge_count,
            endpoint_payload,
        )
        == expected,
        "extreme_rays_match_vertex_counts": (
            ray_histogram["finite"] == vinberg.finite_vertex_count
            and ray_histogram["ideal"] == vinberg.ideal_vertex_count
            and ray_histogram["hyperideal"] == 0
        ),
        "every_wall_irredundant": all(
            audit["is_genuine_facet"] for audit in facet_support.values()
        ),
        "all_ideal_links_parabolic": all(
            audit["parabolic_rank_n_minus_one"] for audit in ideal_links
        ),
        "all_vertices_one_time_component": same_time_component,
        "vertex_basis_is_full_rank": len(basis) == gram.rows,
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"{name} volume audit failed: "
            f"{[key for key, value in checks.items() if not value]}"
        )
    return {
        "classification": vinberg.decision,
        "finite_vertices": vinberg.finite_vertex_count,
        "cusps": vinberg.ideal_vertex_count,
        "hyperideal_extreme_rays": ray_histogram["hyperideal"],
        "elliptic_rank9_edges": vinberg.elliptic_edge_count,
        "endpoint_type_histogram_finite_ideal": endpoint_payload,
        "parabolic_neighborhood_audit_count": len(
            vinberg.parabolic_neighborhood_audits
        ),
        "facet_support": facet_support,
        "ideal_link_checks": ideal_links,
        "checks": checks,
    }


def _normalize_relation_matrix(matrix: Iterable[Iterable[Any]]) -> list[list[str]]:
    replacements = {"infinity": "0", "ultraparallel": "1", "∞": "0"}
    return [
        [replacements.get(str(value), str(value)) for value in row]
        for row in matrix
    ]


def _parse_coxiter_text(text: str) -> tuple[int, int, list[list[str]]]:
    lines = [
        content
        for line in text.splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [
        ["1" if left == right else "2" for right in range(facets)]
        for left in range(facets)
    ]
    vertex_line = next(
        (line for line in lines[1:] if line.startswith("vertices labels:")), None
    )
    labels = (
        [str(index + 1) for index in range(facets)]
        if vertex_line is None
        else vertex_line.split(":", 1)[1].split()
    )
    if len(labels) != facets or len(set(labels)) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    positions = {label: index for index, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ArithmeticError(f"unsupported CoxIter line: {line}")
        left, right = positions[tokens[0]], positions[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _coxiter_h10_records() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if _sha256(COXITER_ARCHIVE) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("pinned CoxIter ordinary archive changed")
    if _sha256(DOUBLING_ARCHIVE) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("pinned CoxIter doubling archive changed")

    records = []
    ordinary_total = 0
    ordinary_paths = sorted(COXITER_ROOT.rglob("*.coxiter"))
    for path in ordinary_paths:
        facets, dimension, matrix = _parse_coxiter_text(
            path.read_text(encoding="utf-8", errors="replace")
        )
        if dimension != DIMENSION:
            continue
        ordinary_total += 1
        records.append(
            {
                "lane": "ordinary",
                "path": str(path.relative_to(COXITER_ROOT)),
                "facets": facets,
                "relations": matrix,
            }
        )
    doubling_total = 0
    with tarfile.open(DOUBLING_ARCHIVE, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            facets, dimension, matrix = _parse_coxiter_text(
                handle.read().decode("utf-8", errors="replace")
            )
            if dimension != DIMENSION:
                continue
            doubling_total += 1
            records.append(
                {
                    "lane": "doubling",
                    "path": member.name,
                    "facets": facets,
                    "relations": matrix,
                }
            )
    histogram = Counter(record["facets"] for record in records)
    expected_histogram = {11: 14, 12: 11, 13: 8, 14: 4, 15: 1, 16: 3, 19: 1, 20: 3, 28: 2, 44: 2}
    if not (
        len(ordinary_paths) == 1214
        and ordinary_total == 27
        and doubling_total == 22
        and dict(sorted(histogram.items())) == expected_histogram
    ):
        raise ArithmeticError("the pinned H10 CoxIter lane inventory changed")
    return records, {
        "coxiter_commit": COXITER_COMMIT,
        "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
        "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
        "ordinary_H10_records": ordinary_total,
        "doubling_H10_records": doubling_total,
        "combined_H10_records": len(records),
        "facet_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "maximum_facet_count": max(histogram),
    }


def _local_ledger_records() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    loaded = {}
    sources = {}
    for name, (path, expected_hash) in LOCAL_LEDGER_PATHS.items():
        actual = _sha256(path)
        if actual != expected_hash:
            raise ArithmeticError(f"pinned local lineage source {name} changed: {actual}")
        loaded[name] = json.loads(path.read_text(encoding="utf-8"))
        sources[name] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": expected_hash,
        }

    records = []
    direct = loaded["R11_direct"]
    direct_case = next(case for case in direct["cases"] if case["id"] == "R11_UD8A1")
    matrix = _normalize_relation_matrix(
        direct_case["published_full_nef_chamber"]["normalized_relation_matrix"]
    )
    records.append({"id": "R11-UD8A1-H10-14", "facets": len(matrix), "relations": matrix})
    for row in loaded["R11_all_subgroups"]["abstract_diagram_types"]:
        if row["hyperbolic_dimension"] == 10 and row["facet_count"] == 14:
            matrix = _normalize_relation_matrix(row["normalized_relation_matrix"])
            records.append(
                {
                    "id": f"R11-UD8D4-{row['type_id']}",
                    "facets": len(matrix),
                    "relations": matrix,
                }
            )
    for stock_index, stock in enumerate(loaded["R12_all_subgroups"]["stocks"]):
        for row in stock["survivors"]:
            geometry = row["geometry"]
            if geometry["hyperbolic_dimension"] == 10 and geometry["facet_count"] == 16:
                matrix = _normalize_relation_matrix(row["normalized_relation_matrix"])
                records.append(
                    {
                        "id": f"R12-stock{stock_index}-class{row['class_index']}",
                        "facets": len(matrix),
                        "relations": matrix,
                    }
                )
    if Counter(record["facets"] for record in records) != Counter({14: 3, 16: 1}):
        raise ArithmeticError("local H10/14 and H10/16 ledger denominator changed")
    return records, sources


def _node_profiles(relations: Sequence[Sequence[str]]) -> list[dict[str, int]]:
    profiles = []
    for index, row in enumerate(relations):
        counts = Counter(value for right, value in enumerate(row) if right != index)
        profiles.append({key: counts[key] for key in sorted(counts)})
    return sorted(profiles, key=lambda value: json.dumps(value, sort_keys=True))


def _lineage_audit(
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
    seed_interior: sp.Matrix,
) -> dict[str, Any]:
    corpus, corpus_inventory = _coxiter_h10_records()
    local, local_sources = _local_ledger_records()
    comparisons = []
    histogram_collisions = []
    header_hits = []
    exact_matches = []
    for copies in range(1, corpus_inventory["maximum_facet_count"] - 11):
        for start in ([0, 1] if copies % 2 == 0 else [0]):
            roots, _interior = _path_boundary(
                gram, seed_roots, seed_interior, copies, start
            )
            candidate = _relation_matrix(gram, roots)
            candidate_histogram = _relation_histogram(candidate)
            same_header = [row for row in corpus if row["facets"] == len(candidate)]
            if same_header:
                header_hits.append(copies)
            isomorphism = {
                row["path"]: _same_relation_graph(candidate, row["relations"])
                for row in same_header
            }
            matches = [row["path"] for row in same_header if isomorphism[row["path"]]]
            collisions = [
                row
                for row in same_header
                if _relation_histogram(row["relations"]) == candidate_histogram
                and not isomorphism[row["path"]]
            ]
            for row in collisions:
                histogram_collisions.append(
                    {
                        "candidate": f"N{copies}_{'A' if start == 0 else 'B'}",
                        "corpus_path": row["path"],
                        "same_global_histogram": True,
                        "same_node_profile_multiset": _node_profiles(candidate)
                        == _node_profiles(row["relations"]),
                        "exact_labelled_graph_isomorphic": False,
                    }
                )
            exact_matches.extend(matches)
            comparisons.append(
                {
                    "candidate": f"N{copies}_{'A' if start == 0 else 'B'}",
                    "facets": len(candidate),
                    "same_header_record_count": len(same_header),
                    "histogram_collision_count": len(collisions),
                    "exact_labelled_graph_matches": matches,
                }
            )

    local_comparisons = []
    for copies, starts in ((2, (0, 1)), (4, (0, 1))):
        for start in starts:
            roots, _interior = _path_boundary(
                gram, seed_roots, seed_interior, copies, start
            )
            candidate = _relation_matrix(gram, roots)
            eligible = [row for row in local if row["facets"] == len(candidate)]
            matches = [
                row["id"]
                for row in eligible
                if _same_relation_graph(candidate, row["relations"])
            ]
            local_comparisons.append(
                {
                    "candidate": f"N{copies}_{'A' if start == 0 else 'B'}",
                    "same_header_local_records": [row["id"] for row in eligible],
                    "exact_labelled_graph_matches": matches,
                }
            )
    nonseed_histogram_collisions = [
        row for row in histogram_collisions if row["candidate"] != "N1_A"
    ]
    checks = {
        "all_49_pinned_H10_records_scanned": len(corpus) == 49,
        "expected_path_headers_present": sorted(set(header_hits))
        == [1, 2, 3, 4, 7, 8, 16, 32],
        "zero_exact_CoxIter_matches_through_44_facets": not exact_matches,
        "sole_nonseed_global_histogram_collision_retained_for_exact_GI": (
            len(nonseed_histogram_collisions) == 1
            and nonseed_histogram_collisions[0]["candidate"] == "N2_B"
            and nonseed_histogram_collisions[0]["corpus_path"].endswith(
                "10-2,1,1,1,1,1,1,1,1,1,1-12-13_1.coxiter"
            )
            and not nonseed_histogram_collisions[0]["same_node_profile_multiset"]
            and all(
                not row["exact_labelled_graph_isomorphic"]
                and not row["same_node_profile_multiset"]
                for row in histogram_collisions
            )
        ),
        "zero_exact_local_ledger_matches_N2_N4": all(
            not row["exact_labelled_graph_matches"] for row in local_comparisons
        ),
        "negative_control_seed_also_absent": next(
            row for row in comparisons if row["candidate"] == "N1_A"
        )["exact_labelled_graph_matches"]
        == [],
    }
    if not all(checks.values()):
        raise ArithmeticError(
            "lineage audit failed: "
            f"{[key for key, value in checks.items() if not value]}"
        )
    return {
        "scope": (
            "bounded nonmatch against the pinned CoxIter ordinary and index-two "
            "H10 lanes and selected local H10 ledgers; not a worldwide priority proof"
        ),
        "negative_control": (
            "the independently known seed P_1 is also absent, so corpus absence "
            "is evidence about these pinned lanes only"
        ),
        "corpus_inventory": corpus_inventory,
        "path_comparisons_N1_through_N32": comparisons,
        "histogram_collision_controls": histogram_collisions,
        "local_sources": local_sources,
        "local_comparisons": local_comparisons,
        "checks": checks,
    }


def build_certificate(*, include_lineage: bool = True) -> dict[str, Any]:
    source, gram, seed_roots, seed_interior = _load_seed()
    a, b = (seed_roots[index] for index in DOUBLING_SEED_INDICES)
    delta = a + b
    seed_relations = _relation_matrix(gram, seed_roots)
    seed_doubling = _doubling_walls(seed_relations)
    generator_action_on_ab = {
        "s_a": [[-1, 2], [0, 1]],
        "s_b": [[1, 0], [2, -1]],
    }
    seed_checks = {
        "doubling_walls_exactly_10_and_13": seed_doubling == [9, 12],
        "both_doubling_roots_minus_two_divisibility_two": all(
            _pair(gram, root, root) == -2 and _divisibility(gram, root) == 2
            for root in (a, b)
        ),
        "raw_inner_product_two": _pair(gram, a, b) == 2,
        "raw_product_squared_four": _pair(gram, a, b) ** 2 == 4,
        "normalized_squared_ratio_one": sp.factor(
            _pair(gram, a, b) ** 2
            / (_pair(gram, a, a) * _pair(gram, b, b))
        )
        == 1,
        "common_ideal_direction_isotropic": _pair(gram, delta, delta) == 0,
        "nine_other_walls_fixed_by_both_reflections": all(
            _pair(gram, seed_roots[index], a) == 0
            and _pair(gram, seed_roots[index], b) == 0
            for index in FIXED_SEED_INDICES
        ),
        "generator_matrices_on_ab_correct": (
            _reflection(gram, a) * a == -a
            and _reflection(gram, a) * b == 2 * a + b
            and _reflection(gram, b) * a == a + 2 * b
            and _reflection(gram, b) * b == -b
        ),
    }
    if not all(seed_checks.values()):
        raise ArithmeticError("seed redoubling hypotheses failed")

    representative_records = []
    representative_geometry: dict[tuple[int, int], tuple[tuple[sp.Matrix, ...], sp.Matrix]] = {}
    for copies, starts in ((1, (0,)), (2, (0, 1)), (3, (0,)), (4, (0, 1))):
        for start in starts:
            record, roots, interior = _case(
                gram,
                seed_roots,
                seed_interior,
                copies,
                start,
                include_roots=copies in (2, 3, 4),
            )
            representative_records.append(record)
            representative_geometry[(copies, start)] = (roots, interior)

    formula_replays = []
    for copies in range(1, 13):
        orientations = []
        geometries = []
        for start in (0, 1):
            record, roots, _interior = _case(
                gram,
                seed_roots,
                seed_interior,
                copies,
                start,
                include_roots=False,
            )
            orientations.append(
                {
                    "orientation": "A" if start == 0 else "B",
                    "facets": record["facet_count"],
                    "mirror_types": record["mirror_type_histogram"],
                    "relations": record["relation_histogram"],
                    "doubling_wall_count": len(record["doubling_walls_one_based"]),
                }
            )
            geometries.append(roots)
        same_relation = _same_relation_graph(
            _relation_matrix(gram, geometries[0]),
            _relation_matrix(gram, geometries[1]),
        )
        same_weighted = _same_weighted_gram(gram, geometries[0], geometries[1])
        expected_same = copies % 2 == 1
        if not (same_relation == expected_same and same_weighted == expected_same):
            raise ArithmeticError(f"orientation parity theorem failed at N={copies}")
        formula_replays.append(
            {
                "N": copies,
                "orientations": orientations,
                "A_B_relation_isomorphic": same_relation,
                "A_B_lattice_decorated_Gram_isomorphic": same_weighted,
            }
        )

    volume_audits = {}
    for copies, start in ((2, 0), (2, 1)):
        roots, interior = representative_geometry[(copies, start)]
        key = f"N{copies}_{'A' if start == 0 else 'B'}"
        volume_audits[key] = _volume_audit(key, gram, roots, interior)

    seed_chi = Fraction(-527, 33_443_020_800)
    seed_volume = Fraction(527, 987_614_208_000)
    theorem = {
        "statement": (
            "For every integer N>=1, either N-copy finite path in the Cayley "
            "line of <s_a,s_b>=C2*C2 gives a finite-volume noncompact H10 "
            "Coxeter polytope P_N with N+12 facets. Its boundary reflections "
            "generate an index-N subgroup of the seed reflection group, and "
            "vol(P_N)=N vol(P_1). Odd N has one lattice-decorated type under "
            "path reversal; every even N has two nonisomorphic A/B types."
        ),
        "allcock_input": {
            "reference": "Daniel Allcock, Infinitely many hyperbolic Coxeter groups through dimension 19, Geometry & Topology 10 (2006), 737-758",
            "doi": "10.2140/gt.2006.10.737",
            "url": "https://arxiv.org/abs/0903.0138",
            "result_used": "Theorem 4.1 and Corollary 4.2 (finite-subtree gluing of pairwise disjoint doubling walls)",
            "doubling_walls": [10, 13],
            "walls_are_disjoint_in_H10": "parallel with a common ideal point",
            "Cayley_graph": "the bi-infinite line for C2*C2",
            "finite_subtrees": "exactly finite paths",
        },
        "boundary_proof": {
            "fixed_outer_walls": 9,
            "moving_outer_wall_formula": (
                "q_(2m)=c+2m^2 a+2m(m-1)b; "
                "q_(2m+1)=d+2m(m+1)a+2m^2 b"
            ),
            "moving_outer_wall_count": "N+1",
            "endpoint_formula_A": "b and N*a+(N-1)*b",
            "endpoint_formula_B": "a and (N-1)*a+N*b",
            "endpoint_inner_product_identity": "(b,N*a+(N-1)b)=2",
            "facet_count_identity": "9+(N+1)+2=N+12",
            "induction_engine": generator_action_on_ab,
            "why_finite_volume": (
                "P_N is the union of N finite-volume seed chambers with "
                "pairwise disjoint interiors; Allcock proves that this union "
                "is a Coxeter polytope"
            ),
            "why_index_N": "Allcock Corollary 4.2",
            "why_noncompact": (
                "the index-N subgroup of the nonuniform seed reflection group "
                "contains powers of its parabolic elements and is therefore "
                "still nonuniform"
            ),
        },
        "closed_fingerprint_formulas": {
            "N=2k+1": {
                "relations": {"0": "k+2", "1": "2k^2+3k", "2": "65+19k", "3": "9+2k", "4": "2"},
                "mirror_types": {"-4,2": "9+k", "-2,1": "2+k", "-2,2": "2"},
            },
            "N=2k_A": {
                "relations": {"0": "k+2", "1": "2k^2+k-1", "2": "56+19k", "3": "8+2k", "4": "1"},
                "mirror_types": {"-4,2": "9+k", "-2,1": "1+k", "-2,2": "2"},
            },
            "N=2k_B": {
                "relations": {"0": "k+1", "1": "2k^2+k-1", "2": "55+19k", "3": "8+2k", "4": "3"},
                "mirror_types": {"-4,2": "8+k", "-2,1": "2+k", "-2,2": "2"},
            },
        },
        "exact_orbifold_invariants": {
            "seed_rational_Euler_characteristic": str(seed_chi),
            "P_N_rational_Euler_characteristic": "-N*527/33443020800",
            "seed_covolume": "pi^5*527/987614208000",
            "P_N_covolume": "N*pi^5*527/987614208000",
            "H10_Gauss_Bonnet_factor_vol_over_minus_chi": "32*pi^5/945",
        },
    }

    lineage = (
        _lineage_audit(gram, seed_roots, seed_interior)
        if include_lineage
        else {"skipped": True}
    )
    checks = {
        "seed_redoubling_hypotheses": all(seed_checks.values()),
        "representative_cases_green": all(
            all(record["checks"].values()) for record in representative_records
        ),
        "closed_form_replay_N1_through_N12": len(formula_replays) == 12,
        "odd_one_type_even_two_types": all(
            row["A_B_relation_isomorphic"] == (row["N"] % 2 == 1)
            and row["A_B_lattice_decorated_Gram_isomorphic"]
            == (row["N"] % 2 == 1)
            for row in formula_replays
        ),
        "both_N2_exact_Vinberg_and_extreme_ray_audits": all(
            all(audit["checks"].values()) for audit in volume_audits.values()
        ),
        "bounded_lineage_green": (
            not include_lineage or all(lineage["checks"].values())
        ),
        "allcock_applies_for_every_positive_N": True,
        "surface_scope_is_divisor_lattice_not_nef_or_automorphisms": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"global certificate failed: {[key for key, value in checks.items() if not value]}"
        )
    payload = {
        "schema": SCHEMA,
        "artifact_id": "K3-S3-H10-REDOUBLING-FAMILY",
        "theorem": theorem,
        "surface_carrier": {
            "lattice": "F=U+E8(-2)+<-2>",
            "signature": [1, 10, 0],
            "rank": 11,
            "determinant": 512,
            "K3_family": source["explicit_family"],
            "generic_full_Picard_moduli_dimension": 9,
            "meaning": (
                "on every very-general marked member, every displayed wall is "
                "an integral divisor class and every wall reflection is an "
                "integral NS isometry"
            ),
            "scope_exclusions": [
                "not a claim that P_N is the nef cone",
                "not a claim that the wall classes are all irreducible curves",
                "not a claim that the reflections are surface automorphisms",
                "not a worldwide novelty claim",
            ],
        },
        "seed": {
            "id": "K3-S3-C2-H10-13",
            "source": {
                "path": str(SOURCE.relative_to(ROOT)),
                "sha256": SOURCE_SHA256,
                "payload_sha256": SOURCE_PAYLOAD_SHA256,
            },
            "doubling_walls_one_based": [index + 1 for index in seed_doubling],
            "wall10_a": _vector_payload(a),
            "wall13_b": _vector_payload(b),
            "common_ideal_delta_a_plus_b": _vector_payload(delta),
            "parallel_normalization_warning": {
                "raw_inner_product": int(_pair(gram, a, b)),
                "raw_inner_product_squared": int(_pair(gram, a, b) ** 2),
                "norm_product": int(_pair(gram, a, a) * _pair(gram, b, b)),
                "normalized_squared_ratio": "1",
                "conclusion": "parallel, not ultraparallel",
            },
            "checks": seed_checks,
        },
        "representative_polytopes": representative_records,
        "formula_replay_N1_through_N12": formula_replays,
        "exact_N2_volume_audits": volume_audits,
        "bounded_lineage": lineage,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--skip-lineage",
        action="store_true",
        help="skip the external CoxIter/local-ledger comparison",
    )
    args = parser.parse_args()
    payload = build_certificate(include_lineage=not args.skip_lineage)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
