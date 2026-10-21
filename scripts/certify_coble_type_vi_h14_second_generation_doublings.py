#!/usr/bin/env python3
"""Certify every second-generation Type-VI Coble H14 doubling.

The two frozen single doubles H14/100 and H14/101 have respectively seven
and ten intrinsic Coxeter-diagram orbits on their 68 doubling walls.  This
certifier doubles one representative of every orbit and checks the resulting
17 index-four chambers from the exact divisor lattices and transported seed
rays.

The output deliberately distinguishes an integral divisor-lattice chamber
from a nef/effective chamber.  All new normals are primitive integral
reflection images in the physical Neron--Severi lattice, but after crossing
the original nef wall they are not asserted effective on the Coble surface.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.audit_jacobian_kummer_h12_parallel_cusps import _affine_component
from scripts.certify_coble_type_vi_h14_single_doublings import (
    _content,
    _divisibility,
    _pair,
    _reflection,
    _unoriented_key,
)
from scripts.enumerate_jacobian_kummer_h10_class0007_components import (
    _components,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SINGLE = ROOT / "data/enumerations/coble_type_vi_h14_single_doublings.json"
SOURCE = ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
SEED_RAYS = (
    ROOT
    / "data/enumerations/shimada_disc16_d5_h14_independent_exact_cone_rays.jsonl.gz"
)
OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_h14_second_generation_doublings.json"
)

PINNED = {
    "single": "d32c9dbd96d1bb97c49736bd0d3d4fb6b967eb3754a771fc61d557d9b3a7d88d",
    "source": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "seed_rays": "5dd8c3b98a2e1afa2a8e6b887bfe8ad6d614de573ee277e3eefb8bab4ba01aab",
}

PARENT_KEYS = ("H14/100", "H14/101")
PARENT_PREFIX = {"H14/100": "A", "H14/101": "B"}
EXPECTED_PARENT_AUTOMORPHISM_ORDERS = {"H14/100": 16, "H14/101": 12}
EXPECTED_ORBITS = {
    "H14/100": [(8, 12), (16, 20), (8, 21), (16, 22), (4, 26), (8, 31), (8, 44)],
    "H14/101": [
        (12, 12),
        (6, 16),
        (6, 17),
        (12, 20),
        (6, 21),
        (6, 25),
        (6, 31),
        (6, 32),
        (6, 39),
        (2, 45),
    ],
}


@dataclass(frozen=True, slots=True)
class RayRecord:
    kind: str
    active: frozenset[int]


@dataclass(slots=True)
class ParentData:
    key: str
    roots: list[sp.Matrix]
    physical_roots: list[sp.Matrix]
    interior: sp.Matrix
    wall_gram: sp.Matrix
    rays: list[RayRecord]
    automorphism_order: int
    doubling_orbits: list[list[int]]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rational(value: Any) -> sp.Rational:
    return sp.Rational(str(value))


def _vector(payload: Sequence[Any]) -> sp.Matrix:
    return sp.Matrix([_rational(value) for value in payload])


def _vector_integers(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1:
        raise ValueError("expected a column vector")
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("an integral vector acquired a denominator")
    return [int(value) for value in vector]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    return [
        abs(int(value))
        for value in smith_normal_form(matrix, domain=ZZ).diagonal()
    ]


def _wall_ratio(wall_gram: sp.Matrix, left: int, right: int) -> Fraction:
    numerator = int(wall_gram[left, right]) ** 2
    denominator = int(wall_gram[left, left] * wall_gram[right, right])
    return Fraction(numerator, denominator)


def _doubling_walls_from_gram(wall_gram: sp.Matrix) -> list[int]:
    return [
        left
        for left in range(wall_gram.rows)
        if all(
            right == left
            or _wall_ratio(wall_gram, left, right) != Fraction(1, 4)
            for right in range(wall_gram.rows)
        )
    ]


def _intrinsic_graph(wall_gram: sp.Matrix) -> nx.Graph:
    """Normalized Coxeter graph, with no divisor-length decoration."""

    graph = nx.Graph()
    graph.add_nodes_from(range(wall_gram.rows))
    for left in range(wall_gram.rows):
        for right in range(left):
            ratio = _wall_ratio(wall_gram, left, right)
            if ratio:
                graph.add_edge(
                    left,
                    right,
                    ratio=f"{ratio.numerator}/{ratio.denominator}",
                )
    return graph


def _intrinsic_automorphism_orbits(
    wall_gram: sp.Matrix, candidates: Sequence[int]
) -> tuple[int, list[list[int]]]:
    graph = _intrinsic_graph(wall_gram)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "ratio", None
        ),
    )
    automorphisms = list(matcher.isomorphisms_iter())
    unseen = set(candidates)
    orbits: list[list[int]] = []
    while unseen:
        representative = min(unseen)
        orbit = sorted(mapping[representative] for mapping in automorphisms)
        orbit = sorted(set(orbit))
        if not set(orbit) <= set(candidates):
            raise ArithmeticError("an intrinsic symmetry left the doubling set")
        orbits.append(orbit)
        unseen -= set(orbit)
    return len(automorphisms), orbits


def _intrinsic_refinement_hash(
    wall_gram: sp.Matrix, rounds: int = 10
) -> dict[str, Any]:
    """An abstract normalized-Gram isometry invariant.

    This is color refinement on the complete edge-colored Coxeter graph.  A
    different final multiset hash proves nonisometry.  Equality would not by
    itself prove isometry, but no equality occurs in the 17-row census.
    """

    size = wall_gram.rows
    ratios = [
        [_wall_ratio(wall_gram, left, right) for right in range(size)]
        for left in range(size)
    ]
    colors = ["facet"] * size
    for _ in range(rounds):
        new_colors: list[str] = []
        for left in range(size):
            signature = Counter(
                (
                    ratios[left][right].numerator,
                    ratios[left][right].denominator,
                    colors[right],
                )
                for right in range(size)
                if right != left
            )
            encoded = tuple(
                sorted(
                    (numerator, denominator, color, multiplicity)
                    for (
                        numerator,
                        denominator,
                        color,
                    ), multiplicity in signature.items()
                )
            )
            new_colors.append(hashlib.sha256(repr(encoded).encode()).hexdigest())
        colors = new_colors
    histogram = Counter(colors)
    return {
        "rounds": rounds,
        "multiset_sha256": hashlib.sha256(
            repr(sorted(histogram.items())).encode()
        ).hexdigest(),
        "stable_color_class_count": len(histogram),
        "largest_color_class": max(histogram.values()),
        "uses_only_normalized_gram": True,
    }


def _seed_ray_records(
    source: dict[str, Any],
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
) -> list[RayRecord]:
    permutation = source["lineage"]["type_VI_to_Shimada_zero_based"]
    inverse = {shimada: type_vi for type_vi, shimada in enumerate(permutation)}
    with gzip.open(SEED_RAYS, "rt", encoding="utf-8") as handle:
        header = json.loads(next(handle))
        raw = [json.loads(line) for line in handle]
    if not (
        header["ray_count"] == 28_872
        and len(raw) == 28_872
        and header["facet_indexing"] == "one_based"
    ):
        raise ArithmeticError("the seed ray archive changed")
    basis = [
        inverse[index - 1] for index in header["coordinate_basis_facet_indices"]
    ]
    basis_matrix = sp.Matrix.hstack(*(seed_roots[index] for index in basis))
    records: list[RayRecord] = []
    for vector, active in raw:
        ray = basis_matrix * sp.Matrix(vector)
        square = _pair(gram, ray, ray)
        if square < 0:
            kind = "hyperideal"
        elif square == 0:
            kind = "ideal"
        else:
            kind = "finite"
        records.append(
            RayRecord(
                kind,
                frozenset(inverse[index - 1] for index in active),
            )
        )
    if Counter(record.kind for record in records) != Counter(
        {"finite": 27_845, "ideal": 1_027}
    ):
        raise ArithmeticError("the seed ray types changed")
    return records


def _transport_parent_rays(
    seed_records: Sequence[RayRecord],
    seed_wall_gram: sp.Matrix,
    glued_seed_wall: int,
    boundary_origins: Sequence[Sequence[Sequence[int]]],
) -> list[RayRecord]:
    origin_to_parent = {
        (copy, seed): parent
        for parent, origins in enumerate(boundary_origins)
        for copy, seed in origins
    }
    output: list[RayRecord] = []
    for record in seed_records:
        if glued_seed_wall not in record.active:
            for copy in range(2):
                output.append(
                    RayRecord(
                        record.kind,
                        frozenset(
                            origin_to_parent[(copy, seed)]
                            for seed in record.active
                        ),
                    )
                )
            continue
        survives = any(
            seed_wall_gram[glued_seed_wall, other] != 0
            for other in record.active
            if other != glued_seed_wall
        )
        if not survives:
            continue
        output.append(
            RayRecord(
                record.kind,
                frozenset(
                    origin_to_parent[(copy, seed)]
                    for seed in record.active - {glued_seed_wall}
                    for copy in range(2)
                ),
            )
        )
    return output


def _prepare_parents() -> tuple[
    dict[str, Any], dict[str, Any], sp.Matrix, sp.Matrix, list[ParentData]
]:
    hashes = {
        "single": _sha256(SINGLE),
        "source": _sha256(SOURCE),
        "seed_rays": _sha256(SEED_RAYS),
    }
    if hashes != PINNED:
        raise ArithmeticError(f"a pinned source changed: {hashes}")
    single = json.loads(SINGLE.read_text())
    source = json.loads(SOURCE.read_text())
    gram = sp.Matrix(source["arithmetic"]["gram"])
    ns_gram = sp.Matrix(source["full_NS_lattice"]["gram"])
    seed_roots = [sp.Matrix(vector) for vector in source["arithmetic"]["wall_vectors"]]
    seed_matrix = sp.Matrix.hstack(*seed_roots)
    seed_wall_gram = seed_matrix.T * gram * seed_matrix
    seed_rays = _seed_ray_records(source, gram, seed_roots)
    parents: list[ParentData] = []
    for key in PARENT_KEYS:
        payload = single["positive_doubles"][key]
        construction = payload["construction"]
        roots = [
            sp.Matrix(vector)
            for vector in construction["wall_vectors_in_arithmetic_basis"]
        ]
        physical_roots = [
            sp.Matrix(vector)
            for vector in construction["physical_NS_divisor_carrier"][
                "wall_vectors_in_NS_basis"
            ]
        ]
        interior = _vector(construction["strict_interior"])
        wall_matrix = sp.Matrix.hstack(*roots)
        wall_gram = wall_matrix.T * gram * wall_matrix
        rays = _transport_parent_rays(
            seed_rays,
            seed_wall_gram,
            construction["selected_seed_wall_zero_based"],
            construction["boundary_origin_seed_facets"],
        )
        expected_rays = payload["ray_census"]
        if Counter(record.kind for record in rays) != Counter(
            {
                "finite": expected_rays["finite"],
                "ideal": expected_rays["ideal"],
            }
        ):
            raise ArithmeticError(f"the reconstructed {key} rays changed")
        candidates = _doubling_walls_from_gram(wall_gram)
        automorphism_order, orbits = _intrinsic_automorphism_orbits(
            wall_gram, candidates
        )
        expected_orbits = EXPECTED_ORBITS[key]
        if not (
            automorphism_order == EXPECTED_PARENT_AUTOMORPHISM_ORDERS[key]
            and [(len(orbit), orbit[0] + 1) for orbit in orbits]
            == expected_orbits
            and set().union(*(set(orbit) for orbit in orbits))
            == set(candidates)
        ):
            raise ArithmeticError(f"the intrinsic {key} orbit census changed")
        parents.append(
            ParentData(
                key,
                roots,
                physical_roots,
                interior,
                wall_gram,
                rays,
                automorphism_order,
                orbits,
            )
        )
    return single, source, gram, ns_gram, parents


def _construct_boundary(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    interior: sp.Matrix,
    selected: int,
) -> tuple[list[dict[str, Any]], sp.Matrix, sp.Matrix]:
    selected_root = roots[selected]
    reflection = _reflection(gram, selected_root)
    doubled_interior = interior + reflection * interior
    if _pair(gram, doubled_interior, doubled_interior) <= 0:
        raise ArithmeticError("a second-generation interior is not timelike")
    boundary: dict[tuple[int, ...], dict[str, Any]] = {}
    for copy in range(2):
        for parent_index, parent_root in enumerate(roots):
            if parent_index == selected:
                continue
            image = parent_root if copy == 0 else reflection * parent_root
            if _pair(gram, doubled_interior, image) < 0:
                image = -image
            if _pair(gram, doubled_interior, image) <= 0:
                raise ArithmeticError("a child wall is not strictly oriented")
            witness = (
                interior
                - _pair(gram, interior, parent_root)
                / _pair(gram, parent_root, parent_root)
                * parent_root
            )
            if copy:
                witness = reflection * witness
            key = _unoriented_key(image)
            record = boundary.setdefault(
                key,
                {"root": image, "witness": witness, "origins": []},
            )
            record["origins"].append((copy, parent_index))
    return list(boundary.values()), doubled_interior, reflection


def _construct_physical_boundary(
    ns_gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    selected: int,
) -> tuple[list[dict[str, Any]], sp.Matrix]:
    selected_root = roots[selected]
    reflection = _reflection(ns_gram, selected_root)
    boundary: dict[tuple[int, ...], dict[str, Any]] = {}
    for copy in range(2):
        for parent_index, parent_root in enumerate(roots):
            if parent_index == selected:
                continue
            image = parent_root if copy == 0 else reflection * parent_root
            key = _unoriented_key(image)
            record = boundary.setdefault(key, {"root": image, "origins": []})
            record["origins"].append((copy, parent_index))
    return list(boundary.values()), reflection


def _transport_child_rays(
    parent: ParentData,
    selected: int,
    boundary_origins: Sequence[Sequence[tuple[int, int]]],
) -> tuple[list[RayRecord], dict[str, Any]]:
    origin_to_child = {
        origin: child
        for child, origins in enumerate(boundary_origins)
        for origin in origins
    }
    output: list[RayRecord] = []
    on_glue: Counter[str] = Counter()
    retained: Counter[str] = Counter()
    flattened: Counter[str] = Counter()
    for record in parent.rays:
        if selected not in record.active:
            for copy in range(2):
                output.append(
                    RayRecord(
                        record.kind,
                        frozenset(
                            origin_to_child[(copy, facet)]
                            for facet in record.active
                        ),
                    )
                )
            continue
        on_glue[record.kind] += 1
        survives = any(
            parent.wall_gram[selected, other] != 0
            for other in record.active
            if other != selected
        )
        if not survives:
            flattened[record.kind] += 1
            continue
        retained[record.kind] += 1
        output.append(
            RayRecord(
                record.kind,
                frozenset(
                    origin_to_child[(copy, facet)]
                    for facet in record.active - {selected}
                    for copy in range(2)
                ),
            )
        )
    counts = Counter(record.kind for record in output)
    if counts["hyperideal"]:
        raise ArithmeticError("a second-generation chamber gained a hyperideal ray")
    if len(set(output)) != len(output):
        raise ArithmeticError("transport produced duplicate projective ray incidences")
    return output, {
        "finite": counts["finite"],
        "ideal": counts["ideal"],
        "hyperideal": counts["hyperideal"],
        "parent_rays_on_glue": dict(sorted(on_glue.items())),
        "retained_on_glue": dict(sorted(retained.items())),
        "flattened_on_glue": dict(sorted(flattened.items())),
        "complete_transport_rule": (
            "two reflected copies off the glue; one retained copy on the glue "
            "iff another active wall is nonorthogonal; otherwise the ray flattens"
        ),
        "transport_lemma_proof": (
            "Off the glue a neighborhood of the convex double lies in one copy, "
            "so every extreme ray is a reflected parent extreme ray. On the glue, "
            "every possible extreme ray is already extreme in the glue face of the "
            "parent. Removing the glued normal h drops its active-normal span by h. "
            "If every other active normal is h-orthogonal, reflection fixes that "
            "span and the ray is flattened. Otherwise s_h(a)-a is a nonzero "
            "multiple of h for some active a, so the two reflected active sets "
            "recover full rank and the ray remains extreme. Reflection and fixing "
            "the glue preserve Lorentzian causal type."
        ),
    }


def _literal_endpoint_audit(
    child_rays: Sequence[RayRecord], wall_gram: sp.Matrix
) -> dict[str, Any]:
    """Literal edge/cusp audit from the complete transported ray incidences."""

    dimension = 14
    positive_gram = [
        [-int(wall_gram[left, right]) for right in range(wall_gram.cols)]
        for left in range(wall_gram.rows)
    ]
    finite_vertices = {
        tuple(sorted(record.active))
        for record in child_rays
        if record.kind == "finite"
    }
    ideal_vertices = {
        tuple(sorted(record.active))
        for record in child_rays
        if record.kind == "ideal"
    }
    if not (
        len(finite_vertices)
        == sum(record.kind == "finite" for record in child_rays)
        and len(ideal_vertices)
        == sum(record.kind == "ideal" for record in child_rays)
        and {len(vertex) for vertex in finite_vertices} == {dimension}
    ):
        raise ArithmeticError("the child finite endpoint incidences are malformed")

    endpoint_counts: defaultdict[tuple[int, ...], list[int]] = defaultdict(
        lambda: [0, 0]
    )
    for vertex in finite_vertices:
        for deleted in range(dimension):
            endpoint_counts[vertex[:deleted] + vertex[deleted + 1 :]][0] += 1

    affine_cache: dict[tuple[int, ...], bool] = {}
    component_signature_histogram: Counter[tuple[int, ...]] = Counter()
    invalid_neighborhood_count = 0
    full_zero_neighborhood_failures = 0
    ideal_endpoint_incidence_count = 0
    for neighborhood in ideal_vertices:
        components = _components(positive_gram, neighborhood)
        component_signature_histogram[
            tuple(sorted(len(component) for component in components))
        ] += 1
        valid = sum(len(component) - 1 for component in components) == dimension - 1
        for component in components:
            if component not in affine_cache:
                affine_cache[component] = _affine_component(
                    [
                        [positive_gram[left][right] for right in component]
                        for left in component
                    ]
                )
            valid &= affine_cache[component]
            component_set = set(component)
            full_neighborhood = tuple(
                facet
                for facet in range(len(positive_gram))
                if facet in component_set
                or all(positive_gram[facet][node] == 0 for node in component)
            )
            if full_neighborhood != neighborhood:
                full_zero_neighborhood_failures += 1
        if not valid:
            invalid_neighborhood_count += 1
        for deleted_nodes in product(*components):
            deleted = set(deleted_nodes)
            edge = tuple(
                sorted(
                    facet
                    for component in components
                    for facet in component
                    if facet not in deleted
                )
            )
            if len(edge) != dimension - 1:
                raise ArithmeticError("a cusp deletion produced a non-edge order")
            endpoint_counts[edge][1] += 1
            ideal_endpoint_incidence_count += 1

    endpoint_count_histogram: Counter[int] = Counter()
    endpoint_type_histogram: Counter[tuple[int, int]] = Counter()
    for finite, ideal in endpoint_counts.values():
        endpoint_count_histogram[finite + ideal] += 1
        endpoint_type_histogram[(finite, ideal)] += 1
    bad_edge_count = sum(
        multiplicity
        for endpoints, multiplicity in endpoint_count_histogram.items()
        if endpoints != 2
    )
    if not (
        endpoint_count_histogram == Counter({2: len(endpoint_counts)})
        and bad_edge_count == 0
        and invalid_neighborhood_count == 0
        and full_zero_neighborhood_failures == 0
    ):
        raise ArithmeticError("a literal second-generation Vinberg audit failed")
    return {
        "method": (
            "literal rank-13 edge endpoint and full affine zero-neighborhood "
            "audit on the complete transported exact ray incidence list"
        ),
        "extreme_ray_dependency": (
            "uses the independently complete union/ray transport certificate; "
            "it is not advertised as a no-ray-import diagram enumeration"
        ),
        "finite_vertex_count": len(finite_vertices),
        "ideal_vertex_count": len(ideal_vertices),
        "elliptic_rank_13_edge_count": len(endpoint_counts),
        "finite_endpoint_incidence_count": len(finite_vertices) * dimension,
        "ideal_endpoint_incidence_count": ideal_endpoint_incidence_count,
        "endpoint_count_histogram": {
            str(key): value for key, value in sorted(endpoint_count_histogram.items())
        },
        "endpoint_type_histogram": {
            f"{finite},{ideal}": value
            for (finite, ideal), value in sorted(endpoint_type_histogram.items())
        },
        "cusp_component_order_signature_histogram": {
            "+".join(map(str, signature)): value
            for signature, value in sorted(component_signature_histogram.items())
        },
        "distinct_affine_components_checked": len(affine_cache),
        "invalid_full_parabolic_neighborhood_count": invalid_neighborhood_count,
        "full_zero_neighborhood_failure_count": full_zero_neighborhood_failures,
        "bad_edge_count": bad_edge_count,
        "decision": "finite_volume_noncompact",
    }


def _root_and_relation_audit(
    gram: sp.Matrix, roots: Sequence[sp.Matrix], wall_gram: sp.Matrix
) -> tuple[dict[str, Any], list[list[str]]]:
    types: Counter[tuple[int, int]] = Counter()
    relations: Counter[str] = Counter()
    ratio_matrix: list[list[str]] = [["0/1"] * len(roots) for _ in roots]
    reflection_checks = []
    for left, root in enumerate(roots):
        types[(int(_pair(gram, root, root)), _divisibility(gram, root))] += 1
        reflection = _reflection(gram, root)
        reflection_checks.append(
            _content(root) == 1
            and all(value.q == 1 for value in reflection)
            and reflection.T * gram * reflection == gram
            and reflection**2 == sp.eye(gram.rows)
        )
        ratio_matrix[left][left] = "1/1"
        for right, other in enumerate(roots[:left]):
            product = wall_gram[left, right]
            if product < 0:
                raise ArithmeticError("a child wall pairing has the wrong orientation")
            ratio = _wall_ratio(wall_gram, left, right)
            text = f"{ratio.numerator}/{ratio.denominator}"
            ratio_matrix[left][right] = ratio_matrix[right][left] = text
            if ratio == 0:
                relations["m2"] += 1
            elif ratio == Fraction(1, 4):
                relations["m3"] += 1
            elif ratio == Fraction(1, 2):
                relations["m4"] += 1
            elif ratio == 1:
                relations["parallel"] += 1
            elif ratio > 1:
                relations["ultraparallel"] += 1
            else:
                relations["invalid_finite"] += 1
    if not all(reflection_checks) or relations["invalid_finite"]:
        raise ArithmeticError("a child Coxeter/reflection gate failed")
    return {
        "primitive_square_divisibility_histogram": {
            f"{square},{divisibility}": value
            for (square, divisibility), value in sorted(types.items())
        },
        "relation_histogram": dict(sorted(relations.items())),
        "all_roots_primitive_and_reflections_integral": True,
        "all_off_diagonal_pairings_nonnegative": True,
        "invalid_finite_angle_count": 0,
    }, ratio_matrix


def certify_case(
    parent: ParentData,
    orbit_index: int,
    gram: sp.Matrix,
    ns_gram: sp.Matrix,
    *,
    literal: bool = True,
) -> dict[str, Any]:
    orbit = parent.doubling_orbits[orbit_index]
    selected = orbit[0]
    if selected not in _doubling_walls_from_gram(parent.wall_gram):
        raise ArithmeticError("the selected parent wall is not a doubling wall")
    arithmetic_records, interior, reflection = _construct_boundary(
        gram, parent.roots, parent.interior, selected
    )
    physical_records, physical_reflection = _construct_physical_boundary(
        ns_gram, parent.physical_roots, selected
    )
    arithmetic_by_origins = {
        tuple(record["origins"]): record for record in arithmetic_records
    }
    physical_by_origins = {
        tuple(record["origins"]): record for record in physical_records
    }
    if set(arithmetic_by_origins) != set(physical_by_origins):
        raise ArithmeticError("arithmetic and physical child walls misalign")
    origins = [tuple(record["origins"]) for record in arithmetic_records]
    roots = [record["root"] for record in arithmetic_records]
    physical_roots = [physical_by_origins[origin]["root"] for origin in origins]
    wall_matrix = sp.Matrix.hstack(*roots)
    physical_matrix = sp.Matrix.hstack(*physical_roots)
    wall_gram = wall_matrix.T * gram * wall_matrix
    physical_wall_gram = physical_matrix.T * ns_gram * physical_matrix

    selected_physical = parent.physical_roots[selected]
    if not (
        _pair(ns_gram, selected_physical, selected_physical) == -1
        and all(value.q == 1 for value in physical_reflection)
        and physical_reflection.T * ns_gram * physical_reflection == ns_gram
        and physical_reflection**2 == sp.eye(ns_gram.rows)
    ):
        raise ArithmeticError("the selected physical divisor reflection failed")

    wall_smith = _smith_diagonal(wall_matrix)
    physical_smith = _smith_diagonal(physical_matrix)
    if not (
        wall_matrix.rank() == physical_matrix.rank() == gram.rows == 15
        and wall_smith == physical_smith == [1] * 15
        and all(_content(root) == 1 for root in physical_roots)
    ):
        raise ArithmeticError("a child wall stock failed primitive full span")

    root_audit, ratio_matrix = _root_and_relation_audit(
        gram, roots, wall_gram
    )
    for left in range(len(physical_roots)):
        for right in range(left + 1):
            if _wall_ratio(physical_wall_gram, left, right) != _wall_ratio(
                wall_gram, left, right
            ):
                raise ArithmeticError("physical/arithmetic normalized Grams differ")

    minimum_witness_slack: sp.Rational | None = None
    for target, record in enumerate(arithmetic_records):
        witness = record["witness"]
        slacks = list(witness.T * gram * wall_matrix)
        if slacks[target] != 0 or any(
            slack <= 0 for index, slack in enumerate(slacks) if index != target
        ):
            raise ArithmeticError("a transported facet witness failed")
        local_minimum = min(
            slack for index, slack in enumerate(slacks) if index != target
        )
        minimum_witness_slack = (
            local_minimum
            if minimum_witness_slack is None
            else min(minimum_witness_slack, local_minimum)
        )

    child_rays, ray_census = _transport_child_rays(parent, selected, origins)
    literal_audit = (
        _literal_endpoint_audit(child_rays, wall_gram)
        if literal
        else {"status": "not_replayed_in_this_call"}
    )

    redoubling = _doubling_walls_from_gram(wall_gram)
    parallel_pairs = [
        (right, left)
        for position, left in enumerate(redoubling)
        for right in redoubling[:position]
        if _wall_ratio(wall_gram, left, right) == 1
    ]
    if not parallel_pairs:
        raise ArithmeticError("a second-generation chamber is not redoublable")
    pair = parallel_pairs[0]
    case_number = orbit_index + 1
    prefix = PARENT_PREFIX[parent.key]
    identifier = (
        f"COBLE-VI3-D4{prefix}{case_number:02d}-H14-{len(roots)}"
    )
    intrinsic = _intrinsic_refinement_hash(wall_gram)
    checks = {
        "selected_wall_is_intrinsic_doubling_wall": True,
        "physical_glued_divisor_square_minus_one": True,
        "physical_glued_reflection_integral": True,
        "all_boundary_normals_primitive_integral_divisor_classes": True,
        "arithmetic_and_physical_normalized_grams_equal": True,
        "all_reflections_integral": root_audit[
            "all_roots_primitive_and_reflections_integral"
        ],
        "all_facets_irredundant": minimum_witness_slack is not None,
        "strict_timelike_interior": bool(
            _pair(gram, interior, interior) > 0
        ),
        "zero_hyperideal_rays": ray_census["hyperideal"] == 0,
        "has_ideal_rays": ray_census["ideal"] > 0,
        "literal_vinberg_edges_pass": (
            not literal or literal_audit["bad_edge_count"] == 0
        ),
        "literal_vinberg_cusps_pass": (
            not literal
            or literal_audit["invalid_full_parabolic_neighborhood_count"] == 0
            and literal_audit["full_zero_neighborhood_failure_count"] == 0
        ),
        "redoublable": bool(parallel_pairs),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"{identifier} failed: {checks}")
    return {
        "id": identifier,
        "parent": parent.key,
        "parent_intrinsic_automorphism_order": parent.automorphism_order,
        "doubling_orbit_index": orbit_index,
        "doubling_orbit_size": len(orbit),
        "doubling_orbit_one_based": [index + 1 for index in orbit],
        "selected_parent_wall_zero_based": selected,
        "selected_parent_wall_one_based": selected + 1,
        "hyperbolic_dimension": 14,
        "facet_count": len(roots),
        "reflection_subgroup_index_in_parent": 2,
        "reflection_subgroup_index_in_H14_60_seed": 4,
        "volume_multiple_of_H14_60_seed": 4,
        "compact": False,
        "ambient_inertia": list(exact_inertia(gram)),
        "wall_gram_inertia_by_full_span": [1, 14, len(roots) - 15],
        "wall_span_rank": wall_matrix.rank(),
        "wall_span_smith": wall_smith,
        "physical_NS_divisor_carrier": {
            "NS": "I_(1,14)",
            "selected_parent_divisor_square": -1,
            "selected_reflection_integral": True,
            "wall_span_smith": physical_smith,
            "wall_vectors_sha256": _json_hash(
                [_vector_integers(root) for root in physical_roots]
            ),
            "normalized_gram_equals_arithmetic_entrywise": True,
            "effectivity_scope": (
                "primitive integral divisor classes only; reflected classes "
                "are not asserted effective and the chamber is not asserted nef"
            ),
        },
        "strict_interior": [str(value) for value in interior],
        "strict_interior_square": str(_pair(gram, interior, interior)),
        "minimum_transported_facet_witness_slack": str(minimum_witness_slack),
        "all_facets_irredundant": True,
        "root_and_relation_audit": root_audit,
        "arithmetic_wall_vectors_sha256": _json_hash(
            [_vector_integers(root) for root in roots]
        ),
        "exact_normalized_squared_gram_sha256": _json_hash(ratio_matrix),
        "ray_census": ray_census,
        "literal_vinberg_endpoint_cusp_audit": literal_audit,
        "abstract_isometry_invariant": intrinsic,
        "redoubling_frontier": {
            "doubling_wall_count": len(redoubling),
            "parallel_pairs_among_doubling_walls": len(parallel_pairs),
            "explicit_parallel_pair_zero_based": list(pair),
            "explicit_parallel_pair_roots": [
                [int(value) for value in roots[index]] for index in pair
            ],
            "allcock_iteration_available": True,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def _factor_archive_cases(
    cases: Sequence[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Factor repeated facet-count audits without weakening the replay.

    Every full audit is performed before this compact archival projection.
    Equality is checked exactly, so a future case that merely shares a facet
    count but not its root/ray/Vinberg profile will stop the certificate.
    """

    profiles: dict[str, Any] = {}
    archived: list[dict[str, Any]] = []
    for case in cases:
        profile_key = f"F{case['facet_count']}"
        ray_profile = dict(case["ray_census"])
        ray_profile.pop("transport_lemma_proof")
        redoubling = case["redoubling_frontier"]
        literal = case["literal_vinberg_endpoint_cusp_audit"]
        literal_profile = (
            literal
            if "bad_edge_count" not in literal
            else {
                key: literal[key]
                for key in (
                    "method",
                    "extreme_ray_dependency",
                    "finite_vertex_count",
                    "ideal_vertex_count",
                    "elliptic_rank_13_edge_count",
                    "endpoint_count_histogram",
                    "endpoint_type_histogram",
                    "distinct_affine_components_checked",
                    "invalid_full_parabolic_neighborhood_count",
                    "full_zero_neighborhood_failure_count",
                    "bad_edge_count",
                    "decision",
                )
            }
        )
        profile = {
            "facet_count": case["facet_count"],
            "root_and_relation_audit": case["root_and_relation_audit"],
            "ray_census": ray_profile,
            "literal_vinberg_endpoint_cusp_audit": literal_profile,
            "redoubling": {
                "doubling_wall_count": redoubling["doubling_wall_count"],
                "parallel_pairs_among_doubling_walls": redoubling[
                    "parallel_pairs_among_doubling_walls"
                ],
            },
        }
        if profile_key in profiles and profiles[profile_key] != profile:
            raise ArithmeticError(
                f"{profile_key} cases no longer share one exact audit profile"
            )
        profiles[profile_key] = profile
        physical = case["physical_NS_divisor_carrier"]
        archived.append(
            {
                "id": case["id"],
                "parent": case["parent"],
                "doubling_orbit_index": case["doubling_orbit_index"],
                "doubling_orbit_size": case["doubling_orbit_size"],
                "doubling_orbit_one_based": case[
                    "doubling_orbit_one_based"
                ],
                "selected_parent_wall_one_based": case[
                    "selected_parent_wall_one_based"
                ],
                "facet_count": case["facet_count"],
                "audit_profile": profile_key,
                "wall_gram_inertia_by_full_span": case[
                    "wall_gram_inertia_by_full_span"
                ],
                "strict_interior_sha256": _json_hash(case["strict_interior"]),
                "strict_interior_square": case["strict_interior_square"],
                "minimum_facet_witness_slack": case[
                    "minimum_transported_facet_witness_slack"
                ],
                "arithmetic_wall_vectors_sha256": case[
                    "arithmetic_wall_vectors_sha256"
                ],
                "physical_NS_wall_vectors_sha256": physical[
                    "wall_vectors_sha256"
                ],
                "exact_normalized_squared_gram_sha256": case[
                    "exact_normalized_squared_gram_sha256"
                ],
                "abstract_isometry_invariant_sha256": case[
                    "abstract_isometry_invariant"
                ]["multiset_sha256"],
                "explicit_parallel_pair_zero_based": redoubling[
                    "explicit_parallel_pair_zero_based"
                ],
                "explicit_parallel_pair_roots_sha256": _json_hash(
                    redoubling["explicit_parallel_pair_roots"]
                ),
                "full_case_audit_sha256": _json_hash(case),
                "all_checks_pass": True,
            }
        )
    return dict(sorted(profiles.items())), archived


def build_certificate(*, literal: bool = True) -> dict[str, Any]:
    single, source, gram, ns_gram, parents = _prepare_parents()
    cases = [
        certify_case(parent, orbit_index, gram, ns_gram, literal=literal)
        for parent in parents
        for orbit_index in range(len(parent.doubling_orbits))
    ]
    invariant_hashes = [
        case["abstract_isometry_invariant"]["multiset_sha256"] for case in cases
    ]
    if len(set(invariant_hashes)) != len(cases):
        raise ArithmeticError("two second-generation abstract invariants collided")
    facet_histogram = Counter(case["facet_count"] for case in cases)
    if facet_histogram != Counter({180: 3, 181: 4, 182: 4, 183: 6}):
        raise ArithmeticError(f"the child facet histogram changed: {facet_histogram}")
    checks = {
        "source_certificates_pass": (
            single["all_checks_pass"] is True and source["all_checks_pass"] is True
        ),
        "two_parent_chambers": len(parents) == 2,
        "parent_intrinsic_automorphism_orders_16_12": [
            parent.automorphism_order for parent in parents
        ]
        == [16, 12],
        "parent_doubling_orbit_counts_7_10": [
            len(parent.doubling_orbits) for parent in parents
        ]
        == [7, 10],
        "seventeen_cases": len(cases) == 17,
        "all_cases_pass": all(case["all_checks_pass"] for case in cases),
        "all_pairwise_nonisometric_by_abstract_invariant": (
            len(set(invariant_hashes)) == 17
        ),
        "all_zero_hyperideal": all(
            case["ray_census"]["hyperideal"] == 0 for case in cases
        ),
        "all_noncompact": all(case["ray_census"]["ideal"] > 0 for case in cases),
        "all_literal_vinberg_audits_pass": (
            not literal
            or all(
                case["literal_vinberg_endpoint_cusp_audit"]["bad_edge_count"] == 0
                and case["literal_vinberg_endpoint_cusp_audit"][
                    "invalid_full_parabolic_neighborhood_count"
                ]
                == 0
                for case in cases
            )
        ),
        "all_physical_integral_divisor_carriers": all(
            case["physical_NS_divisor_carrier"]["selected_reflection_integral"]
            for case in cases
        ),
        "all_irredundant": all(case["all_facets_irredundant"] for case in cases),
        "all_redoublable": all(
            case["redoubling_frontier"]["allcock_iteration_available"]
            for case in cases
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the global second-generation checks failed: {checks}")
    profiles, archived_cases = _factor_archive_cases(cases)
    return {
        "schema": "coxeter.coble_type_vi_h14_second_generation_doublings.v1",
        "artifact_id": "coble_type_vi_h14_second_generation_doublings",
        "status_date": "2026-08-09",
        "direction": "SURFACE -> POLYTOPE positive divisor-lattice factory",
        "inputs": {
            "single_doubling_certificate": {
                "path": str(SINGLE.relative_to(ROOT)),
                "sha256": PINNED["single"],
            },
            "type_vi_full_ns_certificate": {
                "path": str(SOURCE.relative_to(ROOT)),
                "sha256": PINNED["source"],
            },
            "seed_ray_archive": {
                "path": str(SEED_RAYS.relative_to(ROOT)),
                "sha256": PINNED["seed_rays"],
            },
        },
        "classification": {
            "parent_intrinsic_automorphism_orders": {
                parent.key: parent.automorphism_order for parent in parents
            },
            "doubling_wall_count_per_parent": 68,
            "intrinsic_doubling_orbit_counts": {
                parent.key: len(parent.doubling_orbits) for parent in parents
            },
            "second_generation_type_count": len(cases),
            "facet_count_histogram": {
                str(key): value for key, value in sorted(facet_histogram.items())
            },
            "pairwise_nonisometry_proof": (
                "all 17 exact normalized-Gram color-refinement multiset hashes "
                "are distinct; the invariant uses no lattice root-length decoration"
            ),
            "all_reflection_subgroup_index_in_seed": 4,
            "all_volume_multiple_of_seed": 4,
            "ambient_lattice_inertia": [1, 14, 0],
            "every_arithmetic_and_physical_wall_span": {
                "rank": 15,
                "smith_diagonal": [1] * 15,
            },
        },
        "audit_profiles_by_facet_count": profiles,
        "cases": archived_cases,
        "novelty": {
            "new_to_project_abstract_polytopes": 17,
            "worldwide_priority": "open",
            "reason": (
                "Allcock already gives infinite H14 families; no comprehensive "
                "publication-level isometry census was available"
            ),
        },
        "algebraic_scope": {
            "positive": (
                "every wall is a primitive integral physical NS divisor class "
                "and every mirror reflection is integral"
            ),
            "not_claimed": (
                "effectivity of reflected classes, nef-cone equality, or a "
                "surface-automorphism fundamental domain after wall crossing"
            ),
        },
        "polyhedral_union_ray_transport_lemma": {
            "statement": (
                "For a convex Coxeter double C union s_h(C), the complete extreme-"
                "ray list consists of two reflected copies of each parent ray off "
                "h, plus each parent ray on h whose active set contains a wall not "
                "orthogonal to h; an on-h ray with all other active walls "
                "orthogonal to h is flattened. No other extreme rays occur."
            ),
            "proof": (
                "Away from h the union is locally one parent copy. On h any "
                "candidate is an extreme ray of the parent glue face, hence a "
                "parent extreme ray. Deleting h lowers active-normal rank by at "
                "most one. If all remaining active normals are "
                "h-orthogonal, their reflected span is unchanged and has deficient "
                "rank. If some active a is nonorthogonal, s_h(a)-a is a nonzero "
                "multiple of h, so the doubled active normals restore full rank."
            ),
            "causal_type_consequence": (
                "The retained rays are reflected copies or fixed glue rays, so "
                "their finite/ideal/hyperideal type is unchanged; flattening only "
                "removes rays."
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-literal", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    artifact = build_certificate(literal=not arguments.no_literal)
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
            raise SystemExit(f"archived artifact differs: {OUTPUT}")
    if arguments.write:
        OUTPUT.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
