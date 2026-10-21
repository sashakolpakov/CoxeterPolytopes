#!/usr/bin/env python3
"""Exact linear-algebra checks for the integral-packet dimension bridge.

This is deliberately not a Coxeter-polytope census.  It verifies the symbolic
star-Gram obstruction showing that bounded four-node packet types do not bound
global Lorentzian rank, and freezes the rank inequalities used by the bridge
theorem.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/integral_packet_dimension_bridge.json"


def star_gram(leaves: int) -> sp.Matrix:
    """Return -2I plus the adjacency matrix of the star K_{1,leaves}."""

    assert leaves >= 1
    matrix = -2 * sp.eye(leaves + 1)
    for leaf in range(1, leaves + 1):
        matrix[0, leaf] = 1
        matrix[leaf, 0] = 1
    return matrix


def star_characteristic_coefficients(leaves: int) -> list[int]:
    lam = sp.symbols("lambda")
    actual = star_gram(leaves).charpoly(lam).as_expr()
    expected = (lam + 2) ** (leaves - 1) * ((lam + 2) ** 2 - leaves)
    assert sp.expand(actual - expected) == 0
    return [int(coefficient) for coefficient in sp.Poly(actual, lam).all_coeffs()]


def inertia_from_formula(leaves: int) -> tuple[int, int, int]:
    """Return (positive, negative, zero) for the star Gram."""

    if leaves > 4:
        return (1, leaves, 0)
    if leaves == 4:
        return (0, leaves, 1)
    return (0, leaves + 1, 0)


def exact_matrix_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    eigenvalues = matrix.eigenvals()
    positive = negative = zero = 0
    for eigenvalue, multiplicity in eigenvalues.items():
        sign = sp.signsimp(eigenvalue)
        if sign.is_positive:
            positive += multiplicity
        elif sign.is_negative:
            negative += multiplicity
        elif sign.is_zero:
            zero += multiplicity
        else:
            raise AssertionError(f"undecided exact eigenvalue sign: {eigenvalue}")
    return positive, negative, zero


def build_archive() -> dict[str, object]:
    samples = []
    for leaves in range(1, 13):
        matrix = star_gram(leaves)
        if leaves == 4:
            assert matrix.det() == 0
        else:
            assert matrix.det() != 0
        expected_inertia = inertia_from_formula(leaves)
        assert exact_matrix_inertia(matrix) == expected_inertia
        samples.append(
            {
                "leaves": leaves,
                "order": leaves + 1,
                "rank": int(matrix.rank()),
                "inertia_positive_negative_zero": list(expected_inertia),
                "characteristic_coefficients_descending": star_characteristic_coefficients(
                    leaves
                ),
            }
        )

    # Every at-most-four-node principal submatrix is either edgeless or a
    # single star K_{1,j}, j<=3.  (A principal subgraph containing the centre
    # cannot also have isolated vertices.)  The star blocks have eigenvalues
    # -2+sqrt(j), -2-sqrt(j), and -2, hence all listed types are negative
    # definite.  Keep the disconnected edgeless types explicit: this field is
    # an exhaustive isomorphism-type list, not merely a connected-type list.
    local_types = []
    for order in range(1, 5):
        matrix = -2 * sp.eye(order)
        assert exact_matrix_inertia(matrix)[0] == 0
        assert exact_matrix_inertia(matrix)[2] == 0
        local_types.append(
            {
                "kind": "edgeless",
                "order": order,
                "inertia_positive_negative_zero": list(exact_matrix_inertia(matrix)),
            }
        )
    for local_leaves in range(1, 4):
        matrix = star_gram(local_leaves)
        assert exact_matrix_inertia(matrix)[0] == 0
        assert exact_matrix_inertia(matrix)[2] == 0
        local_types.append(
            {
                "kind": "star",
                "star_leaves": local_leaves,
                "order": local_leaves + 1,
                "inertia_positive_negative_zero": list(exact_matrix_inertia(matrix)),
            }
        )

    # The order-six row (five leaves) is the finite-volume boundary case:
    # deleting a leaf leaves the parabolic K_{1,4}.  From six leaves onward,
    # the corresponding deletion is already indefinite.
    deletion_gate = []
    for leaves in range(5, 13):
        delete_center = -2 * sp.eye(leaves)
        delete_leaf = star_gram(leaves - 1)
        deletion_gate.append(
            {
                "leaves": leaves,
                "delete_center_inertia": list(exact_matrix_inertia(delete_center)),
                "delete_leaf_inertia": list(exact_matrix_inertia(delete_leaf)),
                "finite_or_ideal_vertex_deletion_gate": leaves == 5,
            }
        )

    return {
        "schema": "integral_packet_dimension_bridge_exact_v1",
        "linear_rank_theorem": {
            "root_trace_rank": "d*(n+1)",
            "faithful_packet_inequality": "d*(n+1) <= rank(H_Z) <= B",
            "dimension_consequence": "n <= floor(B/d)-1 <= B-1",
            "facet_degree_consequence": "if sum H.D_i <= D then f <= D and n <= D-1",
            "interval_generator_consequence": "d <= 23^t; d <= 23 if one coordinate is primitive",
        },
        "star_counterfamily": {
            "formula": "G_m=-2*I+Adjacency(K_{1,m})",
            "characteristic_polynomial": "(lambda+2)^(m-1)*((lambda+2)^2-m)",
            "entry_alphabet": [-2, 0, 1],
            "at_most_four_node_types": local_types,
            "all_four_node_principal_submatrices_negative_definite": True,
            "lorentzian_for_all_m_at_least_5": True,
            "unbounded_rank": "rank(G_m)=m+1 for m>=5",
            "samples": samples,
            "global_vertex_deletion_gate": deletion_gate,
            "scope": (
                "Lorentzian Coxeter Grams; only m=5 passes every vertex-deletion "
                "elliptic/parabolic gate, and no polytope/maximality promotion is used"
            ),
        },
        "surface_counterfamily": {
            "surface": "Bl_{p_1,...,p_r}(P^2)",
            "picard_rank": "r+1",
            "exceptional_curve_intersection_matrix": "-I_r",
            "every_at_most_four_exceptional_curve_packet": "a principal submatrix of -I_4",
            "conclusion": "fixed surface dimension and bounded individual curve degree do not bound global cohomological rank",
        },
        "open_global_gate": (
            "A uniform bounded cohomological descent must inject the full restriction-of-scalars "
            "root module into one global packet with an independently certified rank/cell ceiling; "
            "four-node type finiteness alone cannot supply this."
        ),
    }


def canonical_json(data: dict[str, object]) -> str:
    return json.dumps(data, sort_keys=True, indent=2) + "\n"


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    payload = canonical_json(build_archive())
    if args.write:
        ARCHIVE.write_text(payload, encoding="utf-8")
    if args.check or not args.write:
        assert ARCHIVE.read_text(encoding="utf-8") == payload, "archive is stale"
    print("PASS: faithful-rank inequality frozen; star local-to-global obstruction exact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
