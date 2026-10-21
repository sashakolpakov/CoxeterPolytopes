#!/usr/bin/env python3
r"""Certify an infinite stabilizer of the complete golden H9 chamber.

The complete chamber

    P9 \subset H(<-2*phi> + E8 + <2>)

has Bugaenko's compact H8 chamber as its complete ``u=0`` facet.  The
nontrivial reversal of the Bugaenko diagram extends by the identity on the
last coordinate.  Since it preserves both the complete facet and its chosen
side, it preserves P9 globally.  Its product with the already certified
Pell-fold involution has trace greater than the dimension, and hence has
infinite order.

This script freezes all finite exact linear-algebra inputs to that argument.
The only non-computational step is the standard facet-adjacency lemma: an
arrangement chamber has exactly one adjacent chamber on either side of the
relative interior of a complete facet.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from scripts import certify_og10_golden_h9_fixed_slice as fixed_slice
from scripts import certify_og10_golden_h9_complete_chamber_involution as pell_fold


ROOT = Path(__file__).resolve().parents[1]
EXTENDED_PREFIX = (
    ROOT / "data/candidates/og10_golden_vinberg_through_4_plus_5phi.json"
)
FIXED_COORDINATES = (0, *range(2, 11))

# The nontrivial automorphism of Bugaenko's eleven-node labeled diagram.
BUGAENKO_REVERSAL = (8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9)

Vector = tuple[Golden, ...]
Matrix = tuple[Vector, ...]


def _sum(values: Iterable[Golden]) -> Golden:
    return sum(values, start=ZERO)


def _identity(size: int) -> Matrix:
    return tuple(
        tuple(ONE if i == j else ZERO for j in range(size))
        for i in range(size)
    )


def _compose(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            _sum(left[i][k] * right[k][j] for k in range(len(left)))
            for j in range(len(right))
        )
        for i in range(len(left))
    )


def _apply(matrix: Matrix, vector: Sequence[Golden]) -> Vector:
    return tuple(
        _sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return _sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def _inverse(matrix: Matrix) -> Matrix:
    """Exact Gauss--Jordan inverse over Q(phi)."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("inverse requires a square matrix")
    work = [
        list(matrix[i])
        + [ONE if i == j else ZERO for j in range(size)]
        for i in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            raise ArithmeticError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][j] - scale * work[column][j]
                for j in range(2 * size)
            ]
    return tuple(
        tuple(work[i][size:])
        for i in range(size)
    )


def _matrix_isometry(matrix: Matrix, gram: Matrix) -> bool:
    columns = tuple(
        tuple(matrix[i][j] for i in range(len(matrix)))
        for j in range(len(matrix))
    )
    return all(
        _pair(columns[i], gram, columns[j]) == gram[i][j]
        for i in range(len(matrix))
        for j in range(len(matrix))
    )


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _encoded_matrix(matrix: Matrix) -> list[list[dict[str, str]]]:
    return [_encoded_vector(row) for row in matrix]


def _bugaenko_reversal(
    fixed_roots: Sequence[Vector],
) -> tuple[Matrix, tuple[Vector, ...]]:
    """Construct the facet reversal from the actual fixed-slice wall basis."""

    u_zero_roots = tuple(root for root in fixed_roots if not root[-1])
    ordered = tuple(
        u_zero_roots[index]
        for index in fixed_slice.BUGAENKO_REORDERING
    )
    basis_roots = tuple(root[:-1] for root in ordered[:9])
    basis_matrix = tuple(
        tuple(basis_roots[j][i] for j in range(9))
        for i in range(9)
    )
    image_matrix = tuple(
        tuple(basis_matrix[i][BUGAENKO_REVERSAL[j]] for j in range(9))
        for i in range(9)
    )
    reversal_nine = _compose(image_matrix, _inverse(basis_matrix))
    reversal = tuple(
        tuple(
            reversal_nine[i][j]
            if i < 9 and j < 9
            else ONE if i == j else ZERO
            for j in range(10)
        )
        for i in range(10)
    )
    return reversal, ordered


def build() -> dict[str, object]:
    ambient_gram = og10_field_gram()
    gram = tuple(
        tuple(ambient_gram[i][j] for j in FIXED_COORDINATES)
        for i in FIXED_COORDINATES
    )
    archived = fixed_slice._load_roots(EXTENDED_PREFIX, 1144)
    roots = tuple(
        tuple(root[i] for i in FIXED_COORDINATES)
        for root in archived
        if fixed_slice._apply_h(root) == root
    )

    reversal, bugaenko_roots = _bugaenko_reversal(roots)
    facet_normal = (ZERO, *(ZERO for _ in range(8)), Golden(-1))

    g = pell_fold._block_matrix((PHI, Golden(-1), PHI, -PHI))
    product = _compose(reversal, g)
    product_trace = _sum(product[i][i] for i in range(10))
    trace_margin = product_trace - Golden(10)

    reversal_wall_images = tuple(
        _apply(reversal, root) for root in bugaenko_roots
    )
    expected_wall_images = tuple(
        bugaenko_roots[BUGAENKO_REVERSAL[i]] for i in range(11)
    )
    reversal_last_row = reversal[-1]
    reversal_last_column = tuple(row[-1] for row in reversal)

    # These two archived theorems supply the global, rather than prefix-level,
    # inputs.  Their own independent tests rebuild their exact certificates.
    fixed_dependency = json.loads(
        (ROOT / "data/candidates/og10_golden_h9_fixed_slice.json").read_text()
    )
    pell_dependency = json.loads(
        (
            ROOT
            / "data/candidates/og10_golden_h9_complete_chamber_involution.json"
        ).read_text()
    )

    checks = {
        "fixed_archive_has_61_roots": len(roots) == 61,
        "complete_bugaenko_facet_dependency_is_verified": (
            fixed_dependency["verified"] is True
            and fixed_dependency["right_angled_bugaenko_facet"][
                "complete_facet_wall_count"
            ]
            == 11
            and fixed_dependency["right_angled_bugaenko_facet"][
                "verified_integral_lattice_identification"
            ]
            is True
        ),
        "pell_fold_complete_chamber_dependency_is_verified": (
            pell_dependency["verified"] is True
            and pell_dependency["checks"][
                "g_is_integral_isometric_involution"
            ]
            is True
            and pell_dependency["checks"][
                "folded_endpoint_is_strict_for_all_61_walls"
            ]
            is True
        ),
        "bugaenko_reversal_is_O_integral": all(
            entry.is_integral for row in reversal for entry in row
        ),
        "bugaenko_reversal_is_an_isometric_involution": (
            _matrix_isometry(reversal, gram)
            and _compose(reversal, reversal) == _identity(10)
        ),
        "bugaenko_reversal_maps_all_11_oriented_facet_walls": (
            len(bugaenko_roots) == 11
            and reversal_wall_images == expected_wall_images
        ),
        "bugaenko_reversal_fixes_facet_normal_and_preserves_u_side": (
            _apply(reversal, facet_normal) == facet_normal
            and reversal_last_row
            == (ZERO, *(ZERO for _ in range(8)), ONE)
            and reversal_last_column
            == (ZERO, *(ZERO for _ in range(8)), ONE)
        ),
        "product_trace_is_8_plus_8phi": product_trace == Golden(8, 8),
        "product_trace_exceeds_real_dimension": trace_margin.sign() > 0,
    }

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_OG10_golden_H9_infinite_chamber_stabilizer",
        "candidate_id": "og10_golden_h9_infinite_stabilizer",
        "verified": all(checks.values()),
        "checks": checks,
        "fixed_lattice": {
            "decomposition": "<-2*phi> + E8 + <2>",
            "rank_over_Z_phi": 10,
            "hyperbolic_dimension": 9,
            "conjugate_form": "positive_definite",
            "arithmetic_group": "uniform",
        },
        "bugaenko_facet_reversal": {
            "diagram_permutation": list(BUGAENKO_REVERSAL),
            "matrix_row_major": _encoded_matrix(reversal),
            "facet_normal": _encoded_vector(facet_normal),
            "complete_chamber_argument": (
                "The reversal is an integral isometry, hence permutes the full "
                "norm-two arrangement.  It maps the complete Bugaenko facet "
                "setwise and preserves its u-side.  At a relative-interior "
                "point of a complete arrangement facet there is exactly one "
                "adjacent chamber on each side.  Therefore the reversal "
                "stabilizes the complete P9 chamber."
            ),
        },
        "infinite_order_product": {
            "pell_fold_g_row_major": _encoded_matrix(g),
            "product_order": "A_times_g",
            "product_matrix_row_major": _encoded_matrix(product),
            "trace": product_trace.to_json(),
            "trace_minus_dimension": trace_margin.to_json(),
            "proof": (
                "A finite-order real 10 by 10 matrix has all eigenvalues on "
                "the unit circle and therefore absolute trace at most 10.  "
                "The identity-embedding trace is 8+8*phi>10, so the product "
                "has infinite order."
            ),
        },
        "group_theoretic_conclusion": {
            "complete_chamber_stabilizer": "infinite",
            "norm_two_reflection_group_index": "infinite",
            "P9_hyperbolic_volume": "infinite",
            "argument": (
                "The norm-two Weyl group is normal in the uniform arithmetic "
                "group and acts simply transitively on arrangement chambers.  "
                "The chamber stabilizer therefore identifies with the "
                "arithmetic quotient by the Weyl group.  Its infinite subgroup "
                "generated by the two involutions forces infinite index and "
                "infinite chamber volume."
            ),
        },
        "ambient_consequence": {
            "golden_H10_chamber_volume": "infinite",
            "enhanced_OG10_birational_centralizer": "infinite",
            "argument": (
                "The established fixed-facet gate says finite H10 volume would "
                "make P9 compact.  Infinite P9 volume rules this out.  The "
                "enhanced-covolume theorem then makes Bir(X,tau) infinite."
            ),
            "scope": (
                "The explicit infinite-order H9 stabilizer need not itself lift "
                "to the archived golden H10 carrier.  Infinitude of the enhanced "
                "birational centralizer is a consequence of the independent "
                "finite-covolume equivalence, not an explicit lift claim."
            ),
        },
        "conclusion": (
            "The complete golden H9 norm-two chamber and the golden OG10 H10 "
            "norm-two chamber both have infinite hyperbolic volume.  This "
            "excludes the present algebraic carrier as a finite-volume or "
            "compact Coxeter-polytope example."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded)


if __name__ == "__main__":
    main()
