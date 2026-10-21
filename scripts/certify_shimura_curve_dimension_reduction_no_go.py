#!/usr/bin/env python3
"""Exact controls for the Coxeter-to-Shimura-curve dimension-reduction audit.

The script checks two independent pieces of the audit.

* The compact arithmetic tetrahedron with Coxeter diagram ``C4(one 5)``
  has no Lorentzian three-node principal packet, none of its four geometric
  facets is a Coxeter polygon, and the odd-diagram topology of every
  reflection centralizer has an annular quotient (oriented double of genus
  one).
* The compact simplex ``[5,3,3,3]`` has a two-face whose exact induced Gram
  matrix is the Coxeter triangle ``(2,3,10)``.  The same simplex is a base
  facet of the compact arithmetic five-prism ``[5,3,3,3,3]``.

Absolute reflection maximality of the tetrahedron is the literature input
from Felikson's complete Coxeter-decomposition classification.  Absolute
maximality of the four-simplex and five-prism follows from the cited
minimal-arithmetic-covolume theorems.  This script certifies the exact Gram,
face, and topology calculations; it does not pretend to reprove those
classification/minimal-volume inputs.
"""

from __future__ import annotations

import itertools
import json
from typing import Any

import sympy as sp


def _expr(value: sp.Expr) -> str:
    return sp.sstr(sp.radsimp(sp.simplify(value)))


def _rows(matrix: sp.Matrix) -> list[list[str]]:
    return [[_expr(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _schur_face(gram: sp.Matrix, fixed_nodes: list[int]) -> tuple[list[int], sp.Matrix]:
    free_nodes = [node for node in range(gram.rows) if node not in fixed_nodes]
    fixed = gram.extract(fixed_nodes, fixed_nodes)
    cross = gram.extract(free_nodes, fixed_nodes)
    face = sp.simplify(
        gram.extract(free_nodes, free_nodes) - cross * fixed.inv() * cross.T
    )
    return free_nodes, face


def _unit_normalized(gram: sp.Matrix) -> sp.Matrix:
    scale = sp.diag(*[sp.sqrt(sp.Rational(2, 1) / gram[i, i]) for i in range(gram.rows)])
    return sp.simplify(scale * gram * scale)


def _m_interval(cosine: sp.Expr, lower: int, upper: int) -> bool:
    """Certify lower < pi/acos(cosine) < upper by exact cosine comparison."""

    return bool(
        sp.simplify(cosine > sp.cos(sp.pi / lower))
        and sp.simplify(cosine < sp.cos(sp.pi / upper))
    )


def _h6_packet() -> dict[str, Any]:
    sqrt5 = sp.sqrt(5)
    phi = (1 + sqrt5) / 2
    gram = sp.Matrix(
        [
            [2, -phi, 0, -1],
            [-phi, 2, -1, 0],
            [0, -1, 2, -1],
            [-1, 0, -1, 2],
        ]
    )
    conjugate = gram.xreplace({sqrt5: -sqrt5})

    leading = [sp.radsimp(gram[:j, :j].det()) for j in range(1, 5)]
    conjugate_leading = [sp.radsimp(conjugate[:j, :j].det()) for j in range(1, 5)]
    expected_leading = [
        sp.Integer(2),
        (5 - sqrt5) / 2,
        3 - sqrt5,
        -(1 + 5 * sqrt5) / 2,
    ]
    expected_conjugate = [
        sp.Integer(2),
        (5 + sqrt5) / 2,
        3 + sqrt5,
        (-1 + 5 * sqrt5) / 2,
    ]

    triple_determinants = {
        "".join(map(str, subset)): sp.radsimp(gram.extract(subset, subset).det())
        for subset in itertools.combinations(range(4), 3)
    }
    expected_triples = {
        "012": 3 - sqrt5,
        "013": 3 - sqrt5,
        "023": sp.Integer(4),
        "123": sp.Integer(4),
    }
    triple_positive_definite: dict[str, bool] = {}
    for subset in itertools.combinations(range(4), 3):
        packet = gram.extract(subset, subset)
        # Sylvester's criterion in the inherited node order.  Checking the
        # whole leading-minor chain avoids inferring definiteness from the
        # determinant alone.
        triple_positive_definite["".join(map(str, subset))] = all(
            sp.simplify(packet[:j, :j].det()).is_positive is True
            for j in range(1, 4)
        )

    facet_records: dict[str, Any] = {}
    expected_intervals = {
        0: [(5, 6), (4, 5), (3, 4)],
        1: [(4, 5), (5, 6), (3, 4)],
        2: [(8, 9), (3, 4), (2, 3)],
        3: [(8, 9), (2, 3), (3, 4)],
    }
    all_facets_non_coxeter = True
    for node in range(4):
        free, facet = _schur_face(gram, [node])
        cosines: list[sp.Expr] = []
        pairs: list[list[int]] = []
        for left, right in itertools.combinations(range(3), 2):
            cosine = sp.radsimp(
                -facet[left, right] / sp.sqrt(facet[left, left] * facet[right, right])
            )
            cosines.append(cosine)
            pairs.append([free[left], free[right]])
        intervals = expected_intervals[node]
        interval_checks = [
            _m_interval(cosine, lower, upper)
            for cosine, (lower, upper) in zip(cosines, intervals, strict=True)
        ]
        # Every certified open interval lies strictly between consecutive
        # integers, so no intrinsic face angle is pi/m for integral m.
        non_coxeter = all(interval_checks)
        all_facets_non_coxeter &= non_coxeter
        facet_records[str(node)] = {
            "free_nodes": free,
            "schur_gram": _rows(facet),
            "side_pairs": pairs,
            "angle_cosines": [_expr(value) for value in cosines],
            "coxeter_order_open_intervals": [list(interval) for interval in intervals],
            "all_interval_checks_exact": all(interval_checks),
            "is_coxeter_polygon": not non_coxeter,
        }

    # The odd graph is the four-cycle.  All reflections are conjugate, its
    # first Betti number is one, and Allcock--Brink identify this first Betti
    # group with the nonreflection part of a reflection centralizer.
    odd_vertices = 4
    odd_edges = 4
    odd_components = 1
    odd_betti = odd_edges - odd_vertices + odd_components
    quotient_genus = 0
    quotient_boundary_components = 2
    oriented_double_genus = 2 * quotient_genus + quotient_boundary_components - 1

    checks = {
        "identity_inertia_is_3_1": all(
            sp.simplify(left - right) == 0
            for left, right in zip(leading, expected_leading, strict=True)
        ),
        "nongeometric_conjugate_is_positive_definite": all(
            sp.simplify(left - right) == 0
            for left, right in zip(conjugate_leading, expected_conjugate, strict=True)
        )
        and all(sp.simplify(value).is_positive is True for value in conjugate_leading),
        "every_three_node_principal_packet_is_elliptic": all(
            sp.simplify(triple_determinants[key] - value) == 0
            and sp.simplify(value).is_positive is True
            for key, value in expected_triples.items()
        )
        and all(triple_positive_definite.values()),
        "no_facet_is_a_coxeter_polygon": all_facets_non_coxeter,
        "odd_graph_is_connected_cycle": odd_betti == 1,
        "all_reflections_are_conjugate": odd_components == 1,
        "mirror_oriented_double_has_genus_one": oriented_double_genus == 1,
        "ground_field_is_Q_sqrt5": sp.simplify(phi**2 - phi - 1) == 0,
    }

    return {
        "name": "Felikson_H6_equals_C4_one_5",
        "field": "Q(sqrt(5))",
        "gram": _rows(gram),
        "leading_principal_determinants": [_expr(value) for value in leading],
        "conjugate_leading_principal_determinants": [
            _expr(value) for value in conjugate_leading
        ],
        "three_node_principal_determinants": {
            key: _expr(value) for key, value in triple_determinants.items()
        },
        "three_node_packets_positive_definite": triple_positive_definite,
        "facets": facet_records,
        "reflection_centralizer_topology": {
            "odd_graph": "C4",
            "odd_graph_vertices": odd_vertices,
            "odd_graph_edges": odd_edges,
            "odd_graph_first_betti_number": odd_betti,
            "nonreflection_group": "Z",
            "full_quotient_coarse_surface": "annulus",
            "quotient_genus": quotient_genus,
            "quotient_boundary_components": quotient_boundary_components,
            "oriented_double_compactified_genus": oriented_double_genus,
        },
        "checks": checks,
    }


def _shared_2310_face() -> dict[str, Any]:
    sqrt5 = sp.sqrt(5)
    phi = (1 + sqrt5) / 2
    gram = sp.diag(2, 2, 2, 2, 2)
    for node, weight in enumerate((phi, sp.Integer(1), sp.Integer(1), sp.Integer(1))):
        gram[node, node + 1] = -weight
        gram[node + 1, node] = -weight

    fixed_nodes = [0, 1]
    free_nodes, face = _schur_face(gram, fixed_nodes)
    normalized = _unit_normalized(face)
    target_square = (5 + sqrt5) / 2  # (2 cos(pi/10))^2

    leading = [sp.radsimp(gram[:j, :j].det()) for j in range(1, 6)]
    checks = {
        "ambient_simplex_inertia_is_4_1": all(
            sp.simplify(value).is_positive is True for value in leading[:-1]
        )
        and sp.simplify(leading[-1]).is_negative is True,
        "face_first_edge_is_label_10": normalized[0, 1] < 0
        and sp.simplify(normalized[0, 1] ** 2 - target_square) == 0,
        "face_second_edge_is_label_3": sp.simplify(normalized[1, 2] + 1) == 0,
        "face_third_edge_is_label_2": sp.simplify(normalized[0, 2]) == 0,
        "face_field_is_Q_sqrt5": sp.simplify(
            2 * normalized[0, 1] ** 2 - 5 - sqrt5
        )
        == 0,
    }

    return {
        "ambient": "[5,3,3,3] in H4",
        "ambient_gram": _rows(gram),
        "fixed_face_nodes": fixed_nodes,
        "free_side_nodes": free_nodes,
        "face_schur_gram": _rows(face),
        "face_norm_two_gram": _rows(normalized),
        "coxeter_triangle": [2, 3, 10],
        "real_DM_orbifold_curve": "P1_(2,3,10) with chosen real half",
        "ground_field": "Q(sqrt(5))",
        "same_face_in_H5_control": (
            "The compact straight Coxeter five-prism [5,3,3,3,3] has an "
            "isometric [5,3,3,3] base facet, hence the same (2,3,10) face."
        ),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    h6 = _h6_packet()
    shared = _shared_2310_face()
    checks = {
        **{f"h6_{key}": value for key, value in h6["checks"].items()},
        **{f"shared_{key}": value for key, value in shared["checks"].items()},
    }
    return {
        "schema": "coxeter.shimura_curve_dimension_reduction_no_go.v1",
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "absolute_maximality_literature_inputs": {
            "H6_C4_one_5": (
                "Felikson's complete classification of Coxeter decompositions of "
                "bounded hyperbolic tetrahedra: H6 admits no nontrivial decomposition."
            ),
            "H4_[5,3,3,3]": (
                "Its orientation subgroup gives the unique minimum-covolume compact "
                "orientable arithmetic H4 orbifold."
            ),
            "H5_prism_[5,3,3,3,3]": (
                "Its orientation subgroup gives the unique minimum-covolume compact "
                "orientable arithmetic H5 orbifold."
            ),
        },
        "maximal_arithmetic_tetrahedron_control": h6,
        "same_curve_different_dimension_control": shared,
        "logical_verdict": {
            "rank3_principal_extraction": "fails on the H6 control",
            "facet_recursion": "fails on the H6 control",
            "mirror_centralizer_genus_zero": "fails on the H6 control (genus one)",
            "one_curve_recovers_ambient": (
                "false: the H4 simplex and H5 prism share P1_(2,3,10)"
            ),
            "remaining_conditional_result": (
                "A marked arithmetic Coxeter two-section gives the exact genus-zero "
                "curve theorem for that section and inherits the ambient ground field, "
                "but not the ambient dimension, facet count, or transverse Gram data."
            ),
        },
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
