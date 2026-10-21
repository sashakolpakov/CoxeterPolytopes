#!/usr/bin/env python3
"""Exact certificate for the [3,3,6] counterexample to a Lanner reduction.

Nikulin's Proposition 3 in arXiv:0708.3991v2 says in its ``u<2``
branch that a minimal hyperbolic subgraph is Lanner.  The marked Coxeter
orthoscheme [3,3,6] is instead quasi-Lanner.  This script checks that claim
using exact SymPy arithmetic only; it imports no reflection-group
classification or dimension theorem.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/nikulin_336_quasilanner_counterexample.json"
SCHEMA = "coxeter.nikulin_336_quasilanner_counterexample.v1"


def _expr(value: sp.Expr) -> str:
    return sp.sstr(sp.factor(value))


def _matrix_rows(matrix: sp.Matrix) -> list[list[str]]:
    return [[_expr(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _principal_determinants(matrix: sp.Matrix) -> list[sp.Expr]:
    return [sp.factor(matrix[:j, :j].det()) for j in range(1, matrix.rows + 1)]


def _negative_definite(matrix: sp.Matrix) -> bool:
    return all(
        sp.simplify((-1) ** j * determinant).is_positive is True
        for j, determinant in enumerate(_principal_determinants(matrix), start=1)
    )


def _negative_semidefinite(matrix: sp.Matrix) -> bool:
    indices = range(matrix.rows)
    for size in range(1, matrix.rows + 1):
        for subset in itertools.combinations(indices, size):
            minor = matrix.extract(subset, subset).det()
            if sp.simplify((-1) ** size * minor).is_nonnegative is not True:
                return False
    return True


def build_certificate() -> dict[str, Any]:
    sqrt3 = sp.sqrt(3)
    gram = sp.Matrix(
        [
            [-2, 1, 0, 0],
            [1, -2, 1, 0],
            [0, 1, -2, sqrt3],
            [0, 0, sqrt3, -2],
        ]
    )

    leading = _principal_determinants(gram)
    if leading != [-2, 3, -4, -1]:
        raise ArithmeticError(f"unexpected leading determinants: {leading}")
    pivots = [leading[0]] + [sp.factor(leading[j] / leading[j - 1]) for j in range(1, 4)]
    if pivots != [-2, -sp.Rational(3, 2), -sp.Rational(4, 3), sp.Rational(1, 4)]:
        raise ArithmeticError(f"unexpected LDL pivots: {pivots}")

    expected_deletions = {
        1: {
            "determinant": 0,
            "rank": 2,
            "charpoly_coefficients": [1, 6, 8, 0],
            "inertia": {"positive": 0, "negative": 2, "zero": 1},
            "type": "affine_[3,6]",
        },
        2: {
            "determinant": -2,
            "rank": 3,
            "charpoly_coefficients": [1, 6, 9, 2],
            "inertia": {"positive": 0, "negative": 3, "zero": 0},
            "type": "elliptic",
        },
        3: {
            "determinant": -6,
            "rank": 3,
            "charpoly_coefficients": [1, 6, 11, 6],
            "inertia": {"positive": 0, "negative": 3, "zero": 0},
            "type": "elliptic",
        },
        4: {
            "determinant": -4,
            "rank": 3,
            "charpoly_coefficients": [1, 6, 10, 4],
            "inertia": {"positive": 0, "negative": 3, "zero": 0},
            "type": "elliptic",
        },
    }

    deletion_records: dict[str, Any] = {}
    for node in range(1, 5):
        deletion = gram.minor_submatrix(node - 1, node - 1)
        expected = expected_deletions[node]
        coefficients = [int(c) for c in deletion.charpoly().all_coeffs()]
        if deletion.det() != expected["determinant"]:
            raise ArithmeticError((node, deletion.det()))
        if deletion.rank() != expected["rank"]:
            raise ArithmeticError((node, deletion.rank()))
        if coefficients != expected["charpoly_coefficients"]:
            raise ArithmeticError((node, coefficients))
        if node == 1:
            null_vector = sp.Matrix([1, 2, sqrt3])
            if deletion * null_vector != sp.zeros(3, 1):
                raise ArithmeticError("affine null-vector check failed")
            if not _negative_semidefinite(deletion):
                raise ArithmeticError("affine deletion is not negative semidefinite")
        elif not _negative_definite(deletion):
            raise ArithmeticError(f"deletion {node} is not negative definite")

        deletion_records[str(node)] = {
            **expected,
            "matrix": _matrix_rows(deletion),
            "negative_semidefinite": _negative_semidefinite(deletion),
            "negative_definite": _negative_definite(deletion),
        }
        if node == 1:
            deletion_records[str(node)]["primitive_null_vector"] = ["1", "2", "sqrt(3)"]

    # Nodes 2 and 3 are the two endpoint-exclusive walls.  Their deletions
    # are elliptic, so both endpoints of the marked edge are finite.
    marked_u = gram[1, 2]
    if marked_u != 1:
        raise ArithmeticError(marked_u)

    # The diagram is a tree, so all field-generating cyclic products reduce
    # to the three nonzero two-cycles a_ij^2.
    edges = [(1, 2), (2, 3), (3, 4)]
    edge_squares = [sp.factor(gram[i - 1, j - 1] ** 2) for i, j in edges]
    if edge_squares != [1, 1, 3]:
        raise ArithmeticError(edge_squares)
    z = sp.symbols("z")
    coefficient_minpoly = sp.Poly(sp.minimal_polynomial(sqrt3, z), z)
    if coefficient_minpoly.all_coeffs() != [1, 0, -3]:
        raise ArithmeticError(coefficient_minpoly)

    checks = {
        "full_determinant_is_minus_one": gram.det() == -1,
        "full_inertia_is_one_positive_three_negative": pivots[-1] > 0
        and all(pivot < 0 for pivot in pivots[:-1]),
        "all_vertex_deletions_are_negative_semidefinite": all(
            record["negative_semidefinite"] for record in deletion_records.values()
        ),
        "node_1_deletion_is_affine": deletion_records["1"]["type"] == "affine_[3,6]",
        "nodes_2_and_3_are_finite_endpoint_deletions": deletion_records["2"][
            "negative_definite"
        ]
        and deletion_records["3"]["negative_definite"],
        "marked_u_is_one_and_less_than_two": marked_u == 1 and marked_u < 2,
        "minimal_indefinite_but_not_lanner": all(
            record["negative_semidefinite"] for record in deletion_records.values()
        )
        and not deletion_records["1"]["negative_definite"],
        "cyclic_products_generate_Q": edge_squares == [1, 1, 3],
        "affine_null_vector_is_exact": deletion_records["1"]["primitive_null_vector"]
        == ["1", "2", "sqrt(3)"],
    }
    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise ArithmeticError(checks)

    return {
        "schema": SCHEMA,
        "status": "exact_counterexample_to_lanner_sentence_field_conclusion_repairable",
        "source": {
            "paper": "V. V. Nikulin, On ground fields of arithmetic hyperbolic reflection groups",
            "arxiv": "0708.3991v2",
            "location": "Proposition 3, printed pages 8--9, u<2 branch",
            "source_pdf_sha256": "bf5e66703c20245ab229b199aa464b3955f3645716d7ae1b967566cadce082d0",
        },
        "diagram": {
            "coxeter_chain": [3, 3, 6],
            "gram_normalization": "diagonal -2; edge entry 2*cos(pi/m)",
            "gram_matrix": _matrix_rows(gram),
            "determinant": int(gram.det()),
            "leading_principal_determinants": [_expr(value) for value in leading],
            "ldl_pivots": [_expr(value) for value in pivots],
            "inertia": {"positive": 1, "negative": 3, "zero": 0},
            "vertex_deletions": deletion_records,
        },
        "marked_finite_edge": {
            "endpoint_exclusive_nodes": [2, 3],
            "u": int(marked_u),
            "u_less_than_2": True,
            "both_endpoint_deletions_negative_definite": True,
        },
        "arithmetic_packet": {
            "graph_is_tree": True,
            "nonzero_two_cycle_products": [int(value) for value in edge_squares],
            "ground_field": "Q",
            "coefficient_field": "Q(sqrt(3))",
            "coefficient_minimal_polynomial_coefficients": [1, 0, -3],
            "all_cyclic_products_are_algebraic_integers": True,
            "v_arithmetic_nonidentity_ground_field_embedding_condition": "vacuous because ground field is Q",
        },
        "geometric_verdict": {
            "minimal_indefinite": True,
            "lanner": False,
            "quasi_lanner": True,
            "finite_volume": True,
            "ideal_vertices": 1,
            "finite_vertices": 3,
            "counterexample_scope": "the printed implication, not the proposition's repaired field-set conclusion",
        },
        "repair": {
            "replacement": "minimal hyperbolic subgraph is Lanner or quasi-Lanner",
            "quasi_lanner_branch": "finite-volume noncompact arithmetic simplex, hence ground field Q",
            "downstream_field_sets_change": False,
            "reason": "Q already belongs to FL^4 and FT",
        },
        "checks": checks,
        "all_checks_pass": True,
        "methodological_scope": {
            "proved_here": "an exact four-node counterexample and its field-preserving repair",
            "not_proved_here": "the five-shape enumeration, N(14), or any global dimension theorem",
            "role_of_nikulin": "hostile regression target, not a project dependency",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the frozen JSON archive")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
