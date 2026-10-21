#!/usr/bin/env python3
"""Exact Gram obstruction for one surviving record-30 ``a=2`` labeling.

The combinatorial/local Coxeter gates leave the labeling frozen below alive.
This certificate performs the missing global metric test.  Its key witness is
an 11-by-11 principal minor of the twice-normalized Gram matrix which contains
no dotted entry and has determinant ``7*(3-sqrt(5)) != 0``.  Thus every choice
of the three dotted entries has rank at least 11, whereas a Coxeter chamber in
``H^9`` has Gram rank at most 10.

The scope is intentionally one exact labeling and one exact placement of two
extra ``m=3`` edges.  It is not a census of the other locally admissible label
states on this nerve.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    is_finite_coxeter_labeling,
    is_lanner_coxeter_labeling,
    pair,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/h9_14_record30_a2_orbit10_gram_obstruction.json"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "a3b5996c71cc616fb18dd99f6c54befab01196d9b7dfb23b964ceb2d7f7eef17"
)

FACET_COUNT = 14
TARGET_HYPERBOLIC_DIMENSION = 9
TARGET_GRAM_RANK = TARGET_HYPERBOLIC_DIMENSION + 1

# The seven inclusion-minimal nonfaces of the surviving abstract nerve.
DOTTED_PAIRS = ((0, 11), (1, 12), (2, 13))
HIGHER_NONFACES = (
    (0, 3, 4, 5),
    (6, 7, 8, 11),
    (3, 4, 5, 9, 10),
    (6, 7, 8, 9, 10),
)
MINIMAL_NONFACES = DOTTED_PAIRS + HIGHER_NONFACES

# This is sparse terminal state 0, followed by the locally admissible pair of
# cross-component edges (0,1) and (2,11).  Every unlisted meeting pair is m=2.
FINITE_NONRIGHT_LABELS = (
    ((0, 1), 3),
    ((0, 5), 3),
    ((2, 11), 3),
    ((3, 5), 3),
    ((3, 10), 3),
    ((4, 5), 5),
    ((6, 8), 3),
    ((6, 9), 3),
    ((7, 8), 5),
    ((8, 11), 3),
    ((9, 10), 3),
)

PRINCIPAL_ORDER = tuple(range(11))
ANCHOR_ORDER = tuple(range(10))

Quadratic = tuple[Fraction, Fraction]  # a + b*sqrt(5)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return _canonical_sha256(payload)


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _finite_labels() -> dict[tuple[int, int], int]:
    dotted = {pair(*edge) for edge in DOTTED_PAIRS}
    labels = {
        pair(left, right): 2
        for left, right in combinations(range(FACET_COUNT), 2)
        if pair(left, right) not in dotted
    }
    labels.update({pair(*edge): label for edge, label in FINITE_NONRIGHT_LABELS})
    return labels


def _actual_vertices() -> tuple[tuple[int, ...], ...]:
    supports = tuple(set(face) for face in MINIMAL_NONFACES)
    return tuple(
        vertex
        for vertex in combinations(range(FACET_COUNT), TARGET_HYPERBOLIC_DIMENSION)
        if not any(support <= set(vertex) for support in supports)
    )


def _twice_gram(
    order: Sequence[int], labels: Mapping[tuple[int, int], int]
) -> sp.Matrix:
    values = {
        2: sp.Integer(0),
        3: -sp.Integer(1),
        4: -sp.sqrt(2),
        5: -(1 + sp.sqrt(5)) / 2,
    }
    return sp.Matrix(
        [
            [
                sp.Integer(2)
                if left == right
                else values[labels[pair(left, right)]]
                for right in order
            ]
            for left in order
        ]
    )


def _qadd(left: Quadratic, right: Quadratic) -> Quadratic:
    return left[0] + right[0], left[1] + right[1]


def _qmul(left: Quadratic, right: Quadratic) -> Quadratic:
    a, b = left
    c, d = right
    return a * c + 5 * b * d, a * d + b * c


def _qscale(value: Quadratic, scalar: int) -> Quadratic:
    return value[0] * scalar, value[1] * scalar


def _qstr(value: Quadratic) -> str:
    a, b = value
    expression = sp.Rational(a.numerator, a.denominator) + sp.Rational(
        b.numerator, b.denominator
    ) * sp.sqrt(5)
    return sp.sstr(expression)


def _fixed_principal_edges() -> tuple[tuple[int, int, int], ...]:
    support = set(PRINCIPAL_ORDER)
    return tuple(
        (edge[0], edge[1], label)
        for edge, label in FINITE_NONRIGHT_LABELS
        if set(edge) <= support
    )


def _forest_matching_expansion() -> tuple[Quadratic, tuple[Quadratic, ...]]:
    """Compute the fixed determinant independently by the forest formula.

    For a weighted forest with diagonal 2, only fixed points and transpositions
    contribute to the determinant.  Hence

        det(2I-A) = sum_M (-1)^|M| 2^(n-2|M|) prod_{e in M} w_e^2,

    where the sum runs over graph matchings.  The two m=5 edges have
    ``w^2=((1+sqrt(5))/2)^2=(3+sqrt(5))/2``; m=3 edges have square one.
    """

    edges = _fixed_principal_edges()
    vertices = set(PRINCIPAL_ORDER)

    # A tiny union-find check makes the premise of the matching formula
    # executable rather than merely asserted.
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right, _label in edges:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            raise AssertionError("the fixed principal non-right graph is not a forest")
        parent[root_right] = root_left

    square_weight: dict[int, Quadratic] = {
        3: (Fraction(1), Fraction(0)),
        5: (Fraction(3, 2), Fraction(1, 2)),
    }
    sums: list[Quadratic] = []
    determinant: Quadratic = (Fraction(0), Fraction(0))
    for size in range(len(PRINCIPAL_ORDER) // 2 + 1):
        weighted_sum: Quadratic = (Fraction(0), Fraction(0))
        for selection in combinations(edges, size):
            endpoints = [vertex for edge in selection for vertex in edge[:2]]
            if len(set(endpoints)) != 2 * size:
                continue
            product: Quadratic = (Fraction(1), Fraction(0))
            for _left, _right, label in selection:
                product = _qmul(product, square_weight[label])
            weighted_sum = _qadd(weighted_sum, product)
        sums.append(weighted_sum)
        determinant = _qadd(
            determinant,
            _qscale(
                weighted_sum,
                (-1) ** size * 2 ** (len(PRINCIPAL_ORDER) - 2 * size),
            ),
        )
    return determinant, tuple(sums)


def build() -> dict[str, object]:
    labels = _finite_labels()
    vertices = _actual_vertices()
    bad_vertices = tuple(
        vertex
        for vertex in vertices
        if not is_finite_coxeter_labeling(vertex, labels)
    )
    bad_lanners = tuple(
        face
        for face in HIGHER_NONFACES
        if not is_lanner_coxeter_labeling(face, labels)
    )
    if len(vertices) != 128 or bad_vertices or bad_lanners:
        raise AssertionError("the frozen local Coxeter survivor changed")

    principal = _twice_gram(PRINCIPAL_ORDER, labels)
    anchor = _twice_gram(ANCHOR_ORDER, labels)
    determinant = sp.factor(principal.det(), extension=sp.sqrt(5))
    anchor_determinant = sp.factor(anchor.det(), extension=sp.sqrt(5))
    expected_determinant = 7 * (3 - sp.sqrt(5))
    expected_anchor = 52 - 24 * sp.sqrt(5)
    if sp.simplify(determinant - expected_determinant) != 0:
        raise AssertionError("the decisive eleven-minor changed")
    if sp.simplify(anchor_determinant - expected_anchor) != 0:
        raise AssertionError("the ten-facet Schur anchor changed")
    if anchor_determinant == 0:
        raise AssertionError("the Schur anchor unexpectedly became singular")

    coupling = principal[:10, 10]
    schur_residual = sp.factor(
        principal[10, 10] - (coupling.T * anchor.inv() * coupling)[0],
        extension=sp.sqrt(5),
    )
    expected_schur = -7 * (9 + 5 * sp.sqrt(5)) / 44
    if sp.simplify(schur_residual - expected_schur) != 0:
        raise AssertionError("the exact Schur residual changed")
    if sp.simplify(anchor_determinant * schur_residual - determinant) != 0:
        raise AssertionError("the Schur determinant identity failed")

    independent_determinant, matching_sums = _forest_matching_expansion()
    if independent_determinant != (Fraction(21), Fraction(-7)):
        raise AssertionError("the independent forest determinant changed")
    if sp.simplify(
        determinant
        - (
            sp.Rational(independent_determinant[0].numerator, independent_determinant[0].denominator)
            + sp.Rational(independent_determinant[1].numerator, independent_determinant[1].denominator)
            * sp.sqrt(5)
        )
    ) != 0:
        raise AssertionError("the two exact determinant computations disagree")
    if not (sp.sqrt(5) < 3):
        raise AssertionError("the exact positivity proof sqrt(5)<3 failed")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record30_a2_orbit10_exact_gram_obstruction",
        "claim": (
            "The frozen locally admissible record-30 a=2 orbit-10 labeling, "
            "with added m=3 edges (0,1) and (2,11), has no Gram realization "
            "in H^9 for any values of its three dotted entries."
        ),
        "candidate": {
            "record_index": 30,
            "addition_schedule": 2,
            "positive_A4_orbit_index": 10,
            "facet_count": FACET_COUNT,
            "minimal_nonfaces": [list(face) for face in MINIMAL_NONFACES],
            "minimal_nonface_masks": [_mask(face) for face in MINIMAL_NONFACES],
            "higher_lanner_supports": [list(face) for face in HIGHER_NONFACES],
            "dotted_pairs": [
                {
                    "pair": list(edge),
                    "twice_gram_variable": name,
                    "required_inequality": f"{name} < -2",
                }
                for edge, name in zip(DOTTED_PAIRS, ("x", "y", "z"), strict=True)
            ],
            "unlisted_meeting_pair_label": 2,
            "finite_nonright_labels": [
                {"pair": list(edge), "label": label}
                for edge, label in FINITE_NONRIGHT_LABELS
            ],
            "rank_nine_vertex_count": len(vertices),
        },
        "local_replay": {
            "all_rank_nine_vertex_diagrams_finite": not bad_vertices,
            "bad_rank_nine_vertices": [list(vertex) for vertex in bad_vertices],
            "all_order_four_and_five_minimal_nonfaces_lanner": not bad_lanners,
            "bad_higher_lanner_supports": [list(face) for face in bad_lanners],
            "warning": (
                "Passing these local tests is necessary but not sufficient; "
                "the global Gram-rank test below is decisive."
            ),
        },
        "metric_obstruction": {
            "normalization": (
                "Twice the normalized Gram: diagonal 2 and finite entry "
                "-2*cos(pi/m)."
            ),
            "fixed_principal_order": list(PRINCIPAL_ORDER),
            "fixed_principal_nonright_edges": [
                {"pair": [left, right], "label": label}
                for left, right, label in _fixed_principal_edges()
            ],
            "contains_no_dotted_entry": True,
            "principal_determinant": sp.sstr(determinant),
            "principal_determinant_expanded": sp.sstr(sp.expand(determinant)),
            "positivity_proof": (
                "3-sqrt(5)>0 because both sides are positive and 3^2=9>5; "
                "therefore 7*(3-sqrt(5)) is nonzero."
            ),
            "schur_elimination": {
                "anchor_order": list(ANCHOR_ORDER),
                "anchor_determinant": sp.sstr(anchor_determinant),
                "last_pivot_residual": sp.sstr(schur_residual),
                "identity": "det(G[0..10]) = det(G[0..9]) * residual",
            },
            "independent_forest_matching_expansion": {
                "fixed_nonright_graph_is_a_forest": True,
                "weighted_matching_sums_by_size": [
                    _qstr(value) for value in matching_sums
                ],
                "determinant": _qstr(independent_determinant),
            },
            "full_gram_rank_lower_bound": 11,
            "H9_gram_rank_upper_bound": TARGET_GRAM_RANK,
            "rank_condition_passes": False,
            "inertia_9_1_4_possible": False,
        },
        "conclusion": {
            "this_exact_labeling_excluded": True,
            "reason": (
                "A principal 11-by-11 minor is nonsingular independently of "
                "the dotted variables, so the full Gram rank is at least 11. "
                "Fourteen wall normals in R^(9,1) have Gram rank at most 10."
            ),
            "new_hyperbolic_polytope": False,
        },
        "scope": {
            "single_sparse_terminal_state": True,
            "single_two_cross_edge_placement": [[0, 1], [2, 11]],
            "all_values_of_the_three_dotted_entries": True,
            "other_sparse_terminal_states_excluded": False,
            "other_cross_edge_placements_excluded": False,
            "three_or_more_cross_edges_excluded": False,
            "whole_orbit10_nerve_excluded": False,
            "whole_record30_a2_schedule_excluded": False,
            "finite_volume_noncompact_statement": False,
        },
    }
    document["certificate_sha256"] = _certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record30_a2_orbit10_exact_gram_obstruction":
        raise ValueError("unexpected record-30 a2 orbit-10 certificate kind")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("the record-30 a2 orbit-10 Gram certificate did not replay")
    digest = _certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("the record-30 a2 orbit-10 certificate digest is invalid")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("the record-30 a2 orbit-10 certificate is not frozen")
    metric = rebuilt["metric_obstruction"]
    return {
        "rank_nine_vertices": rebuilt["candidate"]["rank_nine_vertex_count"],
        "principal_order": metric["fixed_principal_order"],
        "principal_determinant": metric["principal_determinant"],
        "full_gram_rank_lower_bound": metric["full_gram_rank_lower_bound"],
        "H9_gram_rank_upper_bound": metric["H9_gram_rank_upper_bound"],
        "this_exact_labeling_excluded": rebuilt["conclusion"][
            "this_exact_labeling_excluded"
        ],
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--print", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        print(json.dumps(verify(document), indent=2, sort_keys=True))
        return
    document = build()
    if args.print:
        print(json.dumps(document, indent=2, sort_keys=True))
    else:
        print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
