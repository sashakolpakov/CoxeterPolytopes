#!/usr/bin/env python3
"""Certify the explicit E7(-2) near miss behind the H10/45 obstruction.

This is deliberately independent of the universal primitive-glue theorem.
It verifies one primitive E7(-2) embedding in
``T=U(2)^2+D6(-1)``, its integral orthogonal reversal, the exact mismatch
of its discriminant action with the archived Picard transposition, and an
index-four false positive whose E7(-2) is not primitive.  It also records
the deterministic bounded-search observations without treating them as an
exhaustion theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    same_row_lattice,
    saturated_row_kernel,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = (
    ROOT
    / "data"
    / "enumerations"
    / "roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
)
SEARCH = ROOT / "scripts" / "search_r12_h10_45_special_period_embedding.py"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "r12_h10_45_e7_embedding_search_boundary.json"
)

EXPECTED_SHA256 = {
    "upstream": "448bcc47937a2763bbd4ea97466fbc48924dc862e88ba8719e5ab98b0cc4c916",
    "search": "80dc5a15f7132de3260992a311e32528657bdccd130a34c6638586db935269ed",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def matrix_digest(matrix: sp.Matrix) -> str:
    payload = json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    work = [[int(value) & 1 for value in row] for row in rows]
    rank = 0
    if not work:
        return 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and work[index][column]:
                work[index] = [
                    left ^ right for left, right in zip(work[index], work[rank])
                ]
        rank += 1
    return rank


def gf2_nullspace(matrix: sp.Matrix) -> sp.Matrix:
    work = np.array(matrix.tolist(), dtype=int) & 1
    row_count, column_count = work.shape
    rank = 0
    pivots: list[int] = []
    for column in range(column_count):
        choices = np.flatnonzero(work[rank:, column])
        if not len(choices):
            continue
        pivot = rank + int(choices[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in range(row_count):
            if row != rank and work[row, column]:
                work[row] ^= work[rank]
        pivots.append(column)
        rank += 1
    free = [column for column in range(column_count) if column not in pivots]
    rows = []
    for free_column in free:
        vector = [0] * column_count
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = int(work[row, free_column])
        rows.append(vector)
    return sp.Matrix(rows)


def gf2_coordinates(vector: Sequence[int], basis: sp.Matrix) -> list[int]:
    target = tuple(int(value) & 1 for value in vector)
    for bits in itertools.product((0, 1), repeat=basis.rows):
        image = tuple(
            sum(bits[index] * int(basis[index, column]) for index in range(basis.rows))
            & 1
            for column in range(basis.cols)
        )
        if image == target:
            return list(bits)
    raise ArithmeticError("vector is not in the displayed GF(2) span")


def discriminant_action(action: sp.Matrix, half_basis: sp.Matrix) -> sp.Matrix:
    rows = []
    for row in range(half_basis.rows):
        image = half_basis.row(row) * action
        rows.append(gf2_coordinates([int(value) & 1 for value in image], half_basis))
    return sp.Matrix(rows)


def q2(bits: Sequence[int], half_basis: sp.Matrix, gram: sp.Matrix) -> int:
    representative = sp.Matrix([bits]) * half_basis
    numerator = int((representative * gram * representative.T)[0])
    if numerator % 2:
        raise ArithmeticError("half-vector representative is not dual")
    return (numerator // 2) % 4


def fixed_q2_histogram(
    action: sp.Matrix, half_basis: sp.Matrix, gram: sp.Matrix
) -> dict[str, int]:
    histogram = {str(value): 0 for value in range(4)}
    for bits in itertools.product((0, 1), repeat=half_basis.rows):
        row = sp.Matrix([bits])
        if (row * action - row).applyfunc(lambda value: int(value) & 1) != sp.zeros(
            1, half_basis.rows
        ):
            continue
        histogram[str(q2(bits, half_basis, gram))] += 1
    return histogram


def build_certificate() -> dict[str, Any]:
    actual_hashes = {"upstream": sha256(UPSTREAM), "search": sha256(SEARCH)}
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError("a pinned H10/45 source or search script changed")

    upstream = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    case = next(
        row
        for row in upstream["cases"]
        if row["id"] == "R12-D4x2A1x2-S6transposition-H10-45"
    )
    volume = case["vinberg_finite_volume_audit"]
    if not (
        upstream["all_checks_pass"]
        and case["all_checks_pass"]
        and volume["all_checks_pass"]
        and volume["dimension"] == 10
        and volume["facet_count"] == 45
        and volume["elliptic_edge_count"] == 7332
        and volume["finite_vertex_count"] == 1048
        and volume["ideal_vertex_count"] == 93
        and volume["bad_edge_count"] == 0
        and volume["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the upstream H10/45 Vinberg theorem changed")

    n_gram = sp.Matrix(upstream["coordinate_reconstruction"]["basis_gram"])
    n_action = sp.eye(12)
    n_action[1, 1] = n_action[2, 2] = 0
    n_action[1, 2] = n_action[2, 1] = 1
    anti_n = sp.zeros(1, 12)
    anti_n[0, 1], anti_n[0, 2] = 1, -1
    n_half_basis = gf2_nullspace(n_gram)
    n_discr_action = discriminant_action(n_action, n_half_basis)
    n_fixed_histogram = fixed_q2_histogram(n_discr_action, n_half_basis, n_gram)
    if not (
        n_gram.det() == -64
        and exact_inertia(n_gram) == (1, 11, 0)
        and n_action * n_gram * n_action.T == n_gram
        and n_action**2 == sp.eye(12)
        and len((n_action - sp.eye(12)).nullspace()) == 11
        and (anti_n * n_gram * anti_n.T)[0] == -4
        and sp.gcd_list(list(anti_n * n_gram)) == 2
        and n_half_basis.rows == 6
        and gf2_rank(matrix_payload(n_half_basis)) == 6
        and n_fixed_histogram == {"0": 8, "1": 8, "2": 8, "3": 8}
    ):
        raise AssertionError("the archived Picard discriminant action failed")

    simple_d6 = sp.zeros(6)
    for index in range(5):
        simple_d6[index, index] = 1
        simple_d6[index, index + 1] = -1
    simple_d6[5, 4] = simple_d6[5, 5] = 1
    t_gram = sp.zeros(10)
    t_gram[0, 1] = t_gram[1, 0] = 2
    t_gram[2, 3] = t_gram[3, 2] = 2
    t_gram[4:10, 4:10] = -simple_d6 * simple_d6.T
    t_smith_matrix = smith_normal_form(t_gram, domain=ZZ)
    t_smith = [abs(int(t_smith_matrix[index, index])) for index in range(10)]

    e7_cartan = sp.Matrix(
        [
            [2, -1, 0, 0, 0, 0, 0],
            [-1, 2, -1, 0, 0, 0, 0],
            [0, -1, 2, -1, 0, 0, -1],
            [0, 0, -1, 2, -1, 0, 0],
            [0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, -1, 2, 0],
            [0, 0, -1, 0, 0, 0, 2],
        ]
    )
    anti_t = sp.Matrix(
        [
            [0, 1, 0, 1, 0, -1, -1, 0, -1, 0],
            [1, 1, 1, -1, -1, -2, -3, -2, -1, -1],
            [0, 0, 0, 0, 1, 2, 3, 4, 2, 2],
            [0, 0, 0, 0, 0, -1, -1, -2, -1, -2],
            [0, 1, 0, 0, -1, 0, -1, 0, 0, 0],
            [1, 0, 0, -1, 1, 1, 0, 0, -1, 0],
            [1, 0, 0, 1, -1, -2, -3, -2, -1, -1],
        ]
    )
    anti_t_smith_matrix = smith_normal_form(anti_t, domain=ZZ)
    anti_t_smith = [
        abs(int(anti_t_smith_matrix[index, index])) for index in range(7)
    ]
    t_action = sp.simplify(
        sp.eye(10)
        - 2
        * t_gram
        * anti_t.T
        * (anti_t * t_gram * anti_t.T).inv()
        * anti_t
    )
    fixed_t = saturated_row_kernel(t_action - sp.eye(10))
    saturated_anti_t = saturated_row_kernel(t_action + sp.eye(10))
    fixed_t_gram = fixed_t * t_gram * fixed_t.T
    fixed_t_smith_matrix = smith_normal_form(fixed_t_gram, domain=ZZ)
    fixed_t_smith = [
        abs(int(fixed_t_smith_matrix[index, index])) for index in range(3)
    ]
    eigenstack_index = abs(int(fixed_t.col_join(saturated_anti_t).det()))
    if not (
        t_gram.det() == 64
        and exact_inertia(t_gram) == (2, 8, 0)
        and t_smith == [1] * 4 + [2] * 6
        and anti_t * t_gram * anti_t.T == -2 * e7_cartan
        and anti_t_smith == [1] * 7
        and all(value.q == 1 for value in t_action)
        and t_action * t_gram * t_action.T == t_gram
        and t_action**2 == sp.eye(10)
        and sp.factor(t_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 3 * (sp.Symbol("lambda") + 1) ** 7
        and same_row_lattice(anti_t, saturated_anti_t)
        and exact_inertia(fixed_t_gram) == (2, 1, 0)
        and fixed_t_gram.det() == -16
        and fixed_t_smith == [2, 2, 4]
        and all(fixed_t_gram[index, index] % 4 == 0 for index in range(3))
        and eigenstack_index == 8
    ):
        raise AssertionError("the primitive E7(-2) near miss failed")

    t_half_basis = gf2_nullspace(t_gram)
    t_discr_action = discriminant_action(t_action, t_half_basis)
    t_fixed_histogram = fixed_q2_histogram(t_discr_action, t_half_basis, t_gram)
    if not (
        t_half_basis.rows == 6
        and gf2_rank(matrix_payload(t_half_basis)) == 6
        and t_fixed_histogram == {"0": 16, "1": 0, "2": 16, "3": 0}
        and t_fixed_histogram != n_fixed_histogram
    ):
        raise AssertionError("the exact discriminant-action mismatch failed")

    # A tempting index-four glue does produce a lattice with the right coarse
    # invariants and an integral reversal, but one generator is pure K.  Thus
    # the displayed E7(-2) has saturation index two and is not primitive.
    false_fixed = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 4]])
    orthogonal_false_gram = sp.diag(false_fixed, -2 * e7_cartan)
    false_glue_a = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
    false_glue_b = [0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
    false_basis = hermite_normal_form(
        (2 * sp.eye(10)).col_join(sp.Matrix([false_glue_a, false_glue_b])).T
    ).T / 2
    false_gram = sp.simplify(false_basis * orthogonal_false_gram * false_basis.T)
    false_action = sp.simplify(
        false_basis * sp.diag(sp.eye(3), -sp.eye(7)) * false_basis.inv()
    )
    false_k_in_glued_basis = (
        sp.Matrix.hstack(sp.zeros(7, 3), sp.eye(7)) * false_basis.inv()
    )
    false_k_smith_matrix = smith_normal_form(false_k_in_glued_basis, domain=ZZ)
    false_k_smith = [
        abs(int(false_k_smith_matrix[index, index])) for index in range(7)
    ]
    false_gram_smith_matrix = smith_normal_form(false_gram, domain=ZZ)
    false_gram_smith = [
        abs(int(false_gram_smith_matrix[index, index])) for index in range(10)
    ]
    if not (
        false_glue_a[:3] == [0, 0, 0]
        and any(false_glue_a[3:])
        and all(value.q == 1 for value in false_gram)
        and all(false_gram[index, index] % 2 == 0 for index in range(10))
        and abs(false_gram.det()) == 64
        and false_gram_smith == [1] * 4 + [2] * 6
        and all(value.q == 1 for value in false_action)
        and false_action * false_gram * false_action.T == false_gram
        and false_k_smith == [1] * 6 + [2]
    ):
        raise AssertionError("the index-four nonprimitive false positive changed")

    certificate: dict[str, Any] = {
        "schema": "r12_h10_45_e7_embedding_search_boundary_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "case_id": case["id"],
        "archived_picard_action": {
            "ambient_lattice": "U+D4(-1)^2+A1(-1)^2",
            "action": "swap source basis vectors 1 and 2 (zero based)",
            "fixed_rank": 11,
            "coinvariant_lattice": "<-4>",
            "coinvariant_divisibility": 2,
            "discriminant_group": "(Z/2)^6",
            "discriminant_action": matrix_payload(n_discr_action),
            "fixed_q2_histogram": n_fixed_histogram,
        },
        "primitive_near_miss": {
            "complement_lattice": "U(2)^2+D6(-1)",
            "e7_embedding_basis": matrix_payload(anti_t),
            "e7_embedding_smith": anti_t_smith,
            "anti_lattice": "E7(-2)",
            "involution": matrix_payload(t_action),
            "involution_sha256": matrix_digest(t_action),
            "fixed_lattice_gram": matrix_payload(fixed_t_gram),
            "fixed_lattice_signature_determinant": [[2, 1], -16],
            "fixed_lattice_smith": fixed_t_smith,
            "fixed_lattice_norms_divisible_by_four": True,
            "eigensublattice_index": eigenstack_index,
            "discriminant_action": matrix_payload(t_discr_action),
            "fixed_q2_histogram": t_fixed_histogram,
            "equivariant_glue_exists": False,
            "failure_reason": (
                "fixed discriminant classes have only even q2 values, while "
                "the Picard transposition has eight fixed classes in each of "
                "q2=0,1,2,3"
            ),
        },
        "index_four_false_positive": {
            "fixed_stock": "U+<4>",
            "glue_generators_binary": [false_glue_a, false_glue_b],
            "coarse_overlattice_invariants": {
                "signature": [2, 8],
                "determinant": 64,
                "smith": false_gram_smith,
                "integral_reversal": True,
            },
            "pure_K_glue_generator": false_glue_a,
            "K_embedding_smith": false_k_smith,
            "K_saturation_index": 2,
            "rejected_by_primitive_embedding_gate": True,
        },
        "bounded_discovery_search": {
            "theorem_status": "evidence_only_not_used_for_exhaustion",
            "fixed_central_run": {
                "trials": 100000,
                "seed": 990077,
                "complete_embeddings": 57611,
                "primitive_embeddings": 18203,
                "integral_reversals": 9046,
                "fixed_q2_histograms": {"16,0,16,0": 9046},
            },
            "random_central_run": {
                "trials": 100000,
                "seed": 881122,
                "central_batch": 500,
                "complete_embeddings": 47145,
                "primitive_embeddings": 13852,
                "integral_reversals": 6598,
                "fixed_q2_histograms": {"16,0,16,0": 6598},
            },
        },
        "vinberg_handoff": {
            "dimension": 10,
            "facet_count": 45,
            "finite_volume": True,
            "compact": False,
            "elliptic_edge_count": 7332,
            "finite_vertex_count": 1048,
            "ideal_vertex_count": 93,
            "connected_affine_component_count": 302,
            "bad_edge_count": 0,
            "bad_parabolic_neighborhood_count": 0,
            "criterion": volume["criterion"],
        },
        "scope": (
            "explicit primitive near miss and bounded search boundary; the "
            "separate universal primitive-glue theorem supplies the actual "
            "Nikulin-lane no-go"
        ),
        "next_positive_criterion": (
            "an automorphism lift must leave the symplectic order-two "
            "E8(-2) lane, for example by using an order-four or higher "
            "transcendental action with the exact Picard discriminant action"
        ),
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if archived != certificate:
            raise AssertionError("archived search-boundary certificate changed")
    else:
        args.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "all_checks_pass": True,
                "payload_sha256": certificate["certificate_payload_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
