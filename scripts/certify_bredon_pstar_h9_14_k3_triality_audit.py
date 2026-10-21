#!/usr/bin/env python3
"""Hostile polytope-first audit of Bredon's exceptional ADEG polytope P_*.

The script starts from the four-arm theta Coxeter diagram printed in Figure 5
of arXiv:2507.05153.  It reconstructs both crystallographic root-length
orientations, builds the canonical order-three simply-laced unfolding, and
then compares the unique K3-admissible branch with the already frozen
``A18_4_0_C07`` factory and K3-nef-promotion certificates.

It also records the elementary symmetrizability obstruction which prevents a
single D4-triality/simply-laced unfolding theorem from covering all 24 ADEG
polytopes: Bredon's one-6 and three-6 triangles already fail the cycle test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
PROMOTION = (
    ROOT / "data" / "enumerations" / "a18_c01_c07_k3_nef_promotions.json"
)
CONSTRUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "bredon_pstar_h9_14_k3_triality_audit.json"
)

EXPECTED_SHA256 = {
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "promotion": "2ab7acc379e48ac55944e92b8d614c0992d7c482af1300c670f2bf8e899ba1b1",
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
}

SCHEMA = "coxeter.bredon_pstar_h9_14_k3_triality_audit.v1"

# Canonical Figure-5 ordering:
#
#       a_j --(6)-- b_j -- c_j
#      /                       \
#     L                         R        (j=0,1,2,3),
#
# where L--a_j, b_j--c_j, and c_j--R are ordinary (m=3) edges.
# The order is L, a0,b0,c0, ..., a3,b3,c3, R.
LEFT = 0
RIGHT = 13
ARMS = tuple((1 + 3 * j, 2 + 3 * j, 3 + 3 * j) for j in range(4))
CANONICAL_TO_FACTORY = (2, 1, 12, 7, 3, 4, 5, 8, 9, 13, 10, 0, 11, 6)

# One-based source-wall order giving an exact equivariant identification of
# the canonical 24-node unfolding with the 24 short A18 source roots.
UNFOLD_TO_SOURCE_WALL = (
    4,
    16,
    12,
    3,
    17,
    11,
    30,
    9,
    5,
    32,
    22,
    6,
    7,
    23,
    15,
    13,
    14,
    31,
    25,
    18,
    24,
    1,
    19,
    8,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(i)] for i in range(matrix.rows)]


def smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(smith[i, i])) for i in range(min(matrix.shape))]


def canonical_integral_gram(a_component_long: bool) -> sp.Matrix:
    """Return one of the two primitive crystallographic symmetrizations."""

    a_component = {LEFT, *(arm[0] for arm in ARMS)}
    long = a_component if a_component_long else set(range(14)) - a_component
    gram = sp.zeros(14)
    for vertex in range(14):
        gram[vertex, vertex] = 6 if vertex in long else 2
    for a, b, c in ARMS:
        edges = ((LEFT, a, 3), (a, b, 6), (b, c, 3), (c, RIGHT, 3))
        for left, right, label in edges:
            if label == 6:
                value = -3
            elif left in long:
                if right not in long:
                    raise AssertionError("an m=3 edge changed root length")
                value = -3
            else:
                if right in long:
                    raise AssertionError("an m=3 edge changed root length")
                value = -1
            gram[left, right] = gram[right, left] = value
    return gram


def normalized_gram(integral: sp.Matrix) -> list[list[str]]:
    result: list[list[str]] = []
    for i in range(integral.rows):
        row = []
        for j in range(integral.cols):
            if i == j:
                row.append("1")
            elif integral[i, j] == 0:
                row.append("0")
            else:
                ratio = sp.simplify(
                    integral[i, j]
                    / sp.sqrt(integral[i, i] * integral[j, j])
                )
                row.append(str(ratio))
        result.append(row)
    return result


def relation_edges(integral: sp.Matrix) -> list[dict[str, int]]:
    edges = []
    for i in range(integral.rows):
        for j in range(i + 1, integral.cols):
            if integral[i, j] == 0:
                continue
            squared = sp.Rational(
                integral[i, j] ** 2,
                integral[i, i] * integral[j, j],
            )
            label = {sp.Rational(1, 4): 3, sp.Rational(3, 4): 6}.get(squared)
            if label is None:
                raise AssertionError("the canonical diagram ceased to be ADEG")
            edges.append({"left": i, "right": j, "label": label})
    return edges


def simply_laced_unfold(
    integral: sp.Matrix,
) -> tuple[sp.Matrix, list[tuple[int, int]], sp.Matrix, list[int]]:
    """Unfold norm-six vertices into three orthogonal norm-two vertices."""

    long = {i for i in range(integral.rows) if integral[i, i] == 6}
    vertices: list[tuple[int, int]] = []
    for vertex in range(integral.rows):
        vertices.extend(
            [(vertex, phase) for phase in range(3)]
            if vertex in long
            else [(vertex, 0)]
        )
    position = {vertex: index for index, vertex in enumerate(vertices)}
    gram = 2 * sp.eye(len(vertices))
    for i in range(integral.rows):
        for j in range(i + 1, integral.cols):
            if integral[i, j] == 0:
                continue
            if i in long and j in long:
                if integral[i, j] != -3:
                    raise AssertionError("long-long edge is not simply laced after folding")
                for phase in range(3):
                    left, right = position[i, phase], position[j, phase]
                    gram[left, right] = gram[right, left] = -1
            elif (i in long) != (j in long):
                if integral[i, j] != -3:
                    raise AssertionError("unequal lengths do not form an m=6 edge")
                long_vertex = i if i in long else j
                short_vertex = j if i in long else i
                for phase in range(3):
                    left = position[long_vertex, phase]
                    right = position[short_vertex, 0]
                    gram[left, right] = gram[right, left] = -1
            else:
                if integral[i, j] != -1:
                    raise AssertionError("short-short edge is not simply laced")
                left, right = position[i, 0], position[j, 0]
                gram[left, right] = gram[right, left] = -1

    orbit_sum = sp.zeros(integral.rows, len(vertices))
    c3_permutation = list(range(len(vertices)))
    for vertex in range(integral.rows):
        if vertex in long:
            for phase in range(3):
                orbit_sum[vertex, position[vertex, phase]] = 1
                c3_permutation[position[vertex, phase]] = position[
                    vertex, (phase + 1) % 3
                ]
        else:
            orbit_sum[vertex, position[vertex, 0]] = 1
    if orbit_sum * gram * orbit_sum.T != integral:
        raise AssertionError("the orbit sums do not fold to the input Gram")
    return gram, vertices, orbit_sum, c3_permutation


def wall_span_index(rows: sp.Matrix) -> int:
    hnf = hermite_normal_form(rows.T)
    if hnf.rows != hnf.cols:
        raise AssertionError("the walls do not have full lattice rank")
    return abs(int(hnf.det()))


def symmetrizability_test(
    vertex_count: int,
    m3_edges: Iterable[tuple[int, int]],
    m6_edges: Iterable[tuple[int, int]],
) -> dict[str, Any]:
    """Test the necessary and sufficient 1/3 root-length cocycle condition."""

    parent = list(range(vertex_count))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for left, right in m3_edges:
        union(left, right)
    components = {find(vertex) for vertex in range(vertex_count)}
    color: dict[int, int] = {}
    adjacency: dict[int, list[int]] = {component: [] for component in components}
    internal_m6 = []
    for left, right in m6_edges:
        left_component, right_component = find(left), find(right)
        if left_component == right_component:
            internal_m6.append((left, right))
        else:
            adjacency[left_component].append(right_component)
            adjacency[right_component].append(left_component)
    odd_cycle = False
    if not internal_m6:
        for start in components:
            if start in color:
                continue
            color[start] = 0
            stack = [start]
            while stack:
                current = stack.pop()
                for neighbour in adjacency[current]:
                    expected = 1 - color[current]
                    if neighbour in color and color[neighbour] != expected:
                        odd_cycle = True
                    elif neighbour not in color:
                        color[neighbour] = expected
                        stack.append(neighbour)
    passes = not internal_m6 and not odd_cycle
    return {
        "passes": passes,
        "m3_component_count": len(components),
        "internal_m6_edges": [list(edge) for edge in internal_m6],
        "odd_m6_cycle_after_m3_contraction": odd_cycle,
        "criterion": (
            "contract every m=3 component; the m=6 quotient multigraph must "
            "have no loop and must be bipartite"
        ),
    }


def build() -> dict[str, Any]:
    actual_hashes = {
        "factory": sha256(FACTORY),
        "promotion": sha256(PROMOTION),
        "construction": sha256(CONSTRUCTION),
    }
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError("a pinned P-star source certificate changed")

    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    promotion = json.loads(PROMOTION.read_text(encoding="utf-8"))
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    t3 = next(row for row in factory["constructions"] if row["id"] == "A18_4_0_C07")
    c07 = next(row for row in promotion["cases"] if row["id"] == "A18_4_0_C07")
    a18 = next(
        row for row in construction["inventory"]["examples"] if row["key"] == [18, 4, 0]
    )

    canonical = canonical_integral_gram(a_component_long=True)
    reverse = canonical_integral_gram(a_component_long=False)
    factory_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in t3["section"]["restricted_mirrors"]["positive_wall_gram"]]
    )
    if factory_gram.extract(CANONICAL_TO_FACTORY, CANONICAL_TO_FACTORY) != canonical:
        raise AssertionError("the factory T3 Gram is not Bredon's Figure-5 Gram")

    active = [
        row
        for row in t3["section"]["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    ]
    active_coordinates = sp.Matrix(
        [row["primitive_integral_same_direction_coordinates"] for row in active]
    )
    presentation_signs = sp.diag(
        *[sp.Rational(row["presentation_scale_from_primitive"]) for row in active]
    )
    oriented_active_coordinates = presentation_signs * active_coordinates
    fixed_gram = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in t3["section"]["saturated_fixed_lattice"]["gram"]
        ]
    )
    if not (
        active_coordinates.rank() == 10
        and wall_span_index(active_coordinates) == 1
        and oriented_active_coordinates
        * (-fixed_gram)
        * oriented_active_coordinates.T
        == factory_gram
    ):
        raise AssertionError("the fourteen facets do not span the fixed lattice integrally")

    unfolded, unfolded_vertices, orbit_sum, c3_permutation = simply_laced_unfold(canonical)
    reverse_unfolded, _, _, _ = simply_laced_unfold(reverse)

    ambient_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in a18["ambient"]["gram"]]
    )
    source_walls = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in a18["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    source_indices = [index - 1 for index in UNFOLD_TO_SOURCE_WALL]
    source_coordinates = source_walls.extract(source_indices, range(18))
    source_gram = source_coordinates * ambient_gram * source_coordinates.T
    if not (
        source_gram == -unfolded
        and source_coordinates.rank() == 18
        and wall_span_index(source_coordinates) == 1
    ):
        raise AssertionError("the canonical unfolding is not the complete A18 short-root span")

    generator = next(
        row for row in a18["subgroup_classes"] if row["class_index"] == 7
    )["representative_generators"][0]
    source_position = {source: index for index, source in enumerate(source_indices)}
    induced = [source_position[generator[source]] for source in source_indices]
    if induced != c3_permutation:
        raise AssertionError("the canonical triality does not equal the archived C07 action")

    canonical_edges = relation_edges(canonical)
    m3_edges = [
        (edge["left"], edge["right"])
        for edge in canonical_edges
        if edge["label"] == 3
    ]
    m6_edges = [
        (edge["left"], edge["right"])
        for edge in canonical_edges
        if edge["label"] == 6
    ]
    pstar_symmetrizable = symmetrizability_test(14, m3_edges, m6_edges)
    triangle_one_6 = symmetrizability_test(
        3, ((1, 2), (2, 0)), ((0, 1),)
    )
    triangle_three_6 = symmetrizability_test(
        3, (), ((0, 1), (1, 2), (2, 0))
    )
    if not (
        pstar_symmetrizable["passes"]
        and not triangle_one_6["passes"]
        and not triangle_three_6["passes"]
    ):
        raise AssertionError("the ADEG symmetrizability controls changed")

    literal = t3["literal_vinberg"]
    cone = t3["exact_cone"]["exact_double_description"]
    if not (
        cone["finite_vertex_count"] == 128
        and cone["ideal_vertex_count"] == 6
        and cone["hyperideal_vertex_count"] == 0
        and literal["elliptic_edge_count"] == 671
        and literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighborhood_count"] == 0
        and literal["literal_vinberg_passes"]
    ):
        raise AssertionError("the literal finite-volume certificate is not green")

    nef = c07["nef_and_quotient_cone"]
    period = c07["period_and_torelli"]
    coinvariant = c07["ambient_NS"]["coinvariant_lattice"]
    if not (
        promotion["all_checks_pass"]
        and period["picard_lattice_is_exactly_N"]
        and period["invariant_kahler_chamber_exists"]
        and coinvariant["minus_two_rootless"]
        and nef["invariant_nef_section_equality"]
        == "factory cone = P_r intersect Fix(g) = Nef(X) intersect Fix(g)"
    ):
        raise AssertionError("the K3 exact-return certificate is not green")

    return {
        "schema": SCHEMA,
        "source": {
            "primary_source": "Naomi Bredon, On ADEG-polyhedra in hyperbolic spaces, arXiv:2507.05153v1, Figure 5 and Theorem 4.2",
            "primary_source_url": "https://arxiv.org/abs/2507.05153",
            "primary_source_pdf_sha256_downloaded_2026_08_09": "8155bf44b9e0d9888266d71edb6f5706ca86a1edac95fe060dc63147f9f24989",
            "primary_source_tex_sha256_downloaded_2026_08_09": "873a6be52641754c374a520f1163b6e3babccf97b9851f4add7c3e45e7e7cc8b",
            "upstream_certificate_sha256": actual_hashes,
        },
        "pstar": {
            "dimension": 9,
            "facets": 14,
            "canonical_vertex_order": [
                "L",
                "a0",
                "b0",
                "c0",
                "a1",
                "b1",
                "c1",
                "a2",
                "b2",
                "c2",
                "a3",
                "b3",
                "c3",
                "R",
            ],
            "four_arm_symbols": [[3, 6, 3, 3] for _ in range(4)],
            "relation_edges": canonical_edges,
            "normalized_gram_Q_sqrt3": normalized_gram(canonical),
            "published_f_vector": [134, 671, 1480, 1909, 1606, 917, 356, 91, 14],
            "independently_replayed_f_entries": {
                "vertices": 134,
                "edges": 671,
                "facets": 14,
            },
        },
        "integral_orientations": {
            "k3_admissible": {
                "long_vertices": [0, 1, 4, 7, 10],
                "primitive_norm_histogram": {"2": 9, "6": 5},
                "gram": matrix_payload(canonical),
                "rank_inertia": [canonical.rank(), list(exact_inertia(canonical))],
                "smith_diagonal_absolute": smith_diagonal(canonical),
                "fixed_lattice": "U + A2(-1)^4",
                "fixed_lattice_determinant": -81,
                "walls_span_fixed_lattice_index": wall_span_index(active_coordinates),
                "factory_facet_permutation_zero_based": list(CANONICAL_TO_FACTORY),
            },
            "reversed": {
                "long_vertices": sorted(set(range(14)) - {0, 1, 4, 7, 10}),
                "primitive_norm_histogram": {"2": 5, "6": 9},
                "gram": matrix_payload(reverse),
                "rank_inertia": [reverse.rank(), list(exact_inertia(reverse))],
                "smith_diagonal_absolute": smith_diagonal(reverse),
                "candidate_lattice": "U(3) + A2(-1)^4",
                "candidate_lattice_determinant": -729,
            },
        },
        "finite_volume": {
            "classification": "finite_volume_noncompact",
            "finite_ideal_hyperideal": [128, 6, 0],
            "all_extreme_rays_same_time_cone": t3["time_orientation"]["all_extreme_rays_same_time_cone"],
            "strict_interior_pairing_range": [
                t3["time_orientation"]["minimum_strict_interior_pairing"],
                t3["time_orientation"]["maximum_strict_interior_pairing"],
            ],
            "elliptic_edges": literal["elliptic_edge_count"],
            "connected_affine_components": literal["connected_affine_component_count"],
            "bad_edges": literal["bad_edge_count"],
            "bad_parabolic_neighborhoods": literal["bad_parabolic_neighborhood_count"],
            "ray_archive": t3["artifacts"],
        },
        "triality_unfolding": {
            "canonical_node_count": len(unfolded_vertices),
            "canonical_rank_inertia": [unfolded.rank(), list(exact_inertia(unfolded))],
            "canonical_smith_diagonal_absolute": smith_diagonal(unfolded),
            "canonical_characteristic_polynomial": str(
                sp.factor(unfolded.charpoly().as_expr())
            ),
            "folded_gram_equals_integral_pstar_gram": orbit_sum * unfolded * orbit_sum.T == canonical,
            "ambient_lattice": "Nikulin (18,4,0), determinant -16",
            "unfolded_roots_span_ambient_lattice_index": wall_span_index(source_coordinates),
            "unfold_to_source_wall_one_based": list(UNFOLD_TO_SOURCE_WALL),
            "source_gram_equals_negative_unfolded_gram": source_gram == -unfolded,
            "triality_action_equals_archived_C07_action": induced == c3_permutation,
            "reversed_node_count": reverse_unfolded.rows,
            "reversed_rank_inertia": [
                reverse_unfolded.rank(),
                list(exact_inertia(reverse_unfolded)),
            ],
            "reversed_smith_diagonal_absolute": smith_diagonal(reverse_unfolded),
            "reversed_k3_obstruction": "rank 28 and three negative directions, so it is neither a K3 Picard lattice nor Lorentzian",
        },
        "k3_and_surface_return": {
            "ambient_NS": "Nikulin (18,4,0), rank 18, determinant -16",
            "transcendental_lattice": c07["transcendental_and_glue"]["T"],
            "primitive_k3_graph_glue": True,
            "automorphism": "genuine purely non-symplectic order-three K3 automorphism",
            "coinvariant_lattice": coinvariant["isometry_type"],
            "coinvariant_minimum_positive_norm": coinvariant["minimum_positive_norm"],
            "coinvariant_minus_two_rootless": coinvariant["minus_two_rootless"],
            "exact_picard_lattice": period["picard_lattice_is_exactly_N"],
            "exact_invariant_nef_return": nef["invariant_nef_section_equality"],
            "all_restricted_long_walls_redundant": nef["all_restricted_long_walls_redundant"],
            "normal_quotient_nef_return": nef["quotient_nef_conclusion"],
            "surface_realization_novelty": "project-new exact K3 invariant-nef and normal-quotient-nef realization of Bredon's published-new abstract P_star",
        },
        "adeg24_uniformity": {
            "general_criterion": "contract all m=3 components; an admissible simply-laced C3/D4 fold exists exactly when the m=6 quotient has no loop and is bipartite",
            "pstar": pstar_symmetrizable,
            "bredon_H2_triangle_with_one_m6": triangle_one_6,
            "bredon_H2_triangle_with_three_m6": triangle_three_6,
            "whole_24_member_uniform_simply_laced_unfolding": False,
            "proof_scope": "the two displayed H2 counterexamples already refute class-wide universality; this is not an exhaustive per-row K3-carrier census",
        },
        "inverse_score": {
            "single_Pstar": {
                "C": "FAIL: a singleton calibration is not a non-singleton surface-free class",
                "Car": "PASS",
                "Return": "PASS",
                "AG": "PASS: K3 embedding, Torelli, and order-three classification are genuinely AG-native",
                "Pushback": "FAIL: no new theorem on every member of a non-singleton polytope class",
                "level": "L3 one-datum exact-return calibration, not L4",
            },
            "whole_ADEG24": {
                "C": "PASS: Bredon's ADEG definition is surface-free and non-singleton",
                "Car": "FAIL: the one-m6 and three-m6 triangles obstruct every admissible symmetric unfolding",
                "Return": "FAIL class-wide; PASS only for P_star",
                "AG": "conditional PASS on the P_star branch only",
                "Pushback": "FAIL: no AG-native sharpening beyond Bredon's published 24-row classification and n<=11 bound",
                "level": "not an inverse theorem",
            },
        },
        "lineage": {
            "abstract_polytope": "published-new in Bredon 2025",
            "arithmetic_lattice_and_full_reflection_chamber": "known Scharlau-Kirschmer II_(9,1)(3^4)",
            "order_three_K3_deformation_type": "classical",
            "project_new": "exact marked D4-triality reconstruction plus K3 invariant-nef and normal-quotient-nef realization",
            "not_claimed": [
                "a new abstract Coxeter polytope",
                "a new arithmetic reflection group",
                "an exact new volume formula beyond Bredon's commensurability statement",
                "a strict L4 polytope-to-surface theorem",
            ],
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    rendered = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.check_archive:
        if OUTPUT.read_bytes() != rendered:
            raise AssertionError("the frozen P-star audit changed")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(rendered)
    print(json.dumps({
        "all_checks_pass": payload["all_checks_pass"],
        "finite_ideal_hyperideal": payload["finite_volume"]["finite_ideal_hyperideal"],
        "canonical_unfold_rank": payload["triality_unfolding"]["canonical_rank_inertia"][0],
        "reversed_unfold_inertia": payload["triality_unfolding"]["reversed_rank_inertia"][1],
        "ADEG24_uniform_unfolding": payload["adeg24_uniformity"]["whole_24_member_uniform_simply_laced_unfolding"],
    }, sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
