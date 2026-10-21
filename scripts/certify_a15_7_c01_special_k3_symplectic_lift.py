#!/usr/bin/env python3
"""Exact special-period K3 lift of the A15_7_1_C01 H10/27 section.

The archived involution of the even 2-elementary lattice N=(15,7,1)
extends over an explicit graph gluing with

    T = U(2) + U(2) + <-2>^3.

Its anti-invariant lattice in T is a primitive D4(-2).  An exhaustive
seven-bit discriminant-form check supplies the equivariant anti-isometry.
An explicit generic fixed period then gives a projective Picard-rank-19 K3
surface whose full anti-invariant cohomology lattice is E8(-2).  Global
Torelli therefore realizes the action as a symplectic Nikulin involution,
while its fixed Neron--Severi lattice remains the archived rank-11 lattice.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
    same_row_lattice,
    saturated_row_kernel,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
VINBERG = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_vinberg.json"
)
LINEAGE = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_lineage.json"
)
C03_C04_HELPER = (
    ROOT / "scripts" / "certify_a18_c03_c04_special_k3_symplectic_lifts.py"
)
OUTPUT = (
    ROOT / "data" / "enumerations" / "a15_7_c01_special_k3_symplectic_lift.json"
)

EXPECTED_SOURCE_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "lineage": "dea8024f9e6e58abd0483e339774e8023b1bd511a6ac208d7db6e5e709187173",
    "c03_c04_helper": "1c15dd20d9fbb2754b154570ef453423fe45961cd6fb263d4d703f40adc6b598",
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
    if not work:
        return 0
    rank = 0
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


def discriminant_action_on_unit_half_basis(
    action: sp.Matrix, coordinate_indices: Sequence[int]
) -> sp.Matrix:
    rows = []
    coordinate_set = set(coordinate_indices)
    for coordinate in coordinate_indices:
        image = [int(action[coordinate, column]) & 1 for column in range(action.cols)]
        if any(image[column] for column in range(action.cols) if column not in coordinate_set):
            raise AssertionError("discriminant action left the displayed half-basis span")
        rows.append([image[column] for column in coordinate_indices])
    return sp.Matrix(rows)


def q2_of_half_vector(bits: Sequence[int], gram: sp.Matrix) -> int:
    """Return twice the discriminant value, in Z/4."""

    numerator = sum(
        int(bits[i]) * int(gram[i, j]) * int(bits[j])
        for i in range(gram.rows)
        for j in range(gram.cols)
    )
    if numerator % 2:
        raise AssertionError("the displayed half-vector is not dual")
    return (numerator // 2) % 4


def build_certificate() -> dict[str, Any]:
    actual_hashes = {
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
        "c03_c04_helper": sha256(C03_C04_HELPER),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A15 source or exact helper changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        construction["all_checks_pass"]
        and vinberg["all_checks_pass"]
        and lineage["all_checks_pass"]
    ):
        raise AssertionError("an upstream construction/volume/lineage audit failed")

    inventory = next(
        row for row in construction["inventory"]["examples"] if row["key"] == [15, 7, 1]
    )
    case = next(row for row in construction["cases"] if row["id"] == "A15_7_1_C01")
    volume = next(row for row in vinberg["cases"] if row["id"] == "A15_7_1_C01")
    lineage_case = next(row for row in lineage["cases"] if row["id"] == "A15_7_1_C01")
    duplicate = lineage["prior_ledger_exact_duplicate"]
    if not (
        all(case["checks"].values())
        and all(case["construction_checks"].values())
        and case["case"]["subgroup_order"] == 2
        and case["case"]["fixed_rank"] == 11
        and case["case"]["source_wall_count"] == 60
        and case["exact_reduction"]["active_hyperbolic_facet_count"] == 27
        and volume["classification"] == "finite_volume_noncompact_Coxeter_polytope"
        and volume["finite_volume"] is True
        and volume["compact"] is False
        and volume["bad_edge_count"] == 0
        and volume["bad_parabolic_neighborhood_count"] == 0
        and lineage_case["id"] == "A15_7_1_C01"
        and duplicate["source"] == "A15_7_1_C01"
        and duplicate["same_signed_normalized_gram"] is True
        and duplicate["target"].startswith("D16-D4-C7-H10-27")
    ):
        raise AssertionError("the upstream A15_7_1_C01 target changed")

    n_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in inventory["ambient"]["gram"]]
    )
    n_action = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["case"]["representative_nonidentity_primal_actions"][0]
        ]
    )
    wall_rows = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in inventory["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    if not (
        n_gram.det() == 128
        and exact_inertia(n_gram) == (1, 14, 0)
        and n_action**2 == sp.eye(15)
        and n_action * n_gram * n_action.T == n_gram
        and len((n_action - sp.eye(15)).nullspace()) == 11
        and sp.factor(n_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 11 * (sp.Symbol("lambda") + 1) ** 4
        and wall_rows.rows == 60
    ):
        raise AssertionError("ambient N=(15,7,1) action failed")

    # Seven independent unit half-vectors account for the entire order-2^7
    # discriminant group of N.
    discr_coordinates = [3, 4, 7, 8, 11, 12, 13]
    n_discr_basis = sp.zeros(7, 15)
    for index, coordinate in enumerate(discr_coordinates):
        n_discr_basis[index, coordinate] = 1
    if not (
        gf2_rank(matrix_payload(n_discr_basis)) == 7
        and all(
            int(value) % 2 == 0
            for value in n_discr_basis * n_gram
        )
    ):
        raise AssertionError("displayed N discriminant basis failed")
    n_discr_action = discriminant_action_on_unit_half_basis(
        n_action, discr_coordinates
    )

    g0 = sp.Matrix(
        [
            [0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, -1, 0, 0],
            [0, 0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 0, 0, -1],
        ]
    )
    t_gram = 2 * g0
    anti_t_basis = sp.Matrix(
        [
            [0, 0, 0, 0, 1, 1, 0],
            [-2, -2, -2, 1, -2, 1, -1],
            [-2, -1, -2, 0, -2, 1, 1],
            [-2, 0, 0, 2, 0, -1, -1],
        ]
    )
    d4_negative = sp.Matrix(
        [
            [-2, 1, 1, 1],
            [1, -2, 0, 0],
            [1, 0, -2, 0],
            [1, 0, 0, -2],
        ]
    )
    anti_t_smith = smith_normal_form(anti_t_basis, domain=ZZ)
    if not (
        anti_t_basis * g0 * anti_t_basis.T == d4_negative
        and [abs(int(anti_t_smith[i, i])) for i in range(4)] == [1, 1, 1, 1]
    ):
        raise AssertionError("the primitive D4(-2) in T failed")
    t_action = sp.simplify(
        sp.eye(7)
        - 2
        * g0
        * anti_t_basis.T
        * (anti_t_basis * g0 * anti_t_basis.T).inv()
        * anti_t_basis
    )
    expected_t_action = sp.Matrix(
        [
            [25, 14, 18, -11, 12, -12, 4],
            [48, 25, 32, -24, 20, -20, 8],
            [-24, -11, -13, 14, -8, 8, -6],
            [32, 18, 24, -13, 16, -16, 4],
            [-20, -12, -16, 8, -11, 10, -2],
            [20, 12, 16, -8, 10, -11, 2],
            [-8, -4, -4, 6, -2, 2, -3],
        ]
    )
    fixed_t_basis = sp.Matrix(
        [
            [2, 0, 0, 0, 1, -1, 1],
            [0, 1, 2, 1, 1, -1, -1],
            [0, 2, 2, -3, -1, 1, -1],
        ]
    )
    fixed_t_gram = fixed_t_basis * t_gram * fixed_t_basis.T
    if not (
        t_action == expected_t_action
        and all(value.q == 1 for value in t_action)
        and t_action**2 == sp.eye(7)
        and t_action * t_gram * t_action.T == t_gram
        and sp.factor(t_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 3 * (sp.Symbol("lambda") + 1) ** 4
        and same_row_lattice(
            fixed_t_basis, saturated_row_kernel(t_action - sp.eye(7))
        )
        and same_row_lattice(
            anti_t_basis, saturated_row_kernel(t_action + sp.eye(7))
        )
        and exact_inertia(fixed_t_gram) == (2, 1, 0)
        and fixed_t_gram.det() == -32
    ):
        raise AssertionError("the T involution or its eigenspaces failed")

    identification = sp.Matrix(
        [
            [1, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 0, 1, 1],
            [0, 1, 0, 1, 0, 1, 1],
            [1, 0, 1, 0, 1, 1, 0],
            [0, 1, 0, 0, 0, 1, 1],
            [0, 1, 1, 0, 1, 0, 1],
            [0, 0, 0, 0, 1, 1, 1],
        ]
    )
    t_mod_two = t_action.applyfunc(lambda value: int(value) & 1)
    if not (
        int(identification.det()) % 2 == 1
        and (n_discr_action * identification - identification * t_mod_two).applyfunc(
            lambda value: int(value) & 1
        )
        == sp.zeros(7)
    ):
        raise AssertionError("the discriminant identification is not equivariant")

    q_table = []
    q2_histogram = {str(value): 0 for value in range(4)}
    for source in itertools.product(range(2), repeat=7):
        n_bits = [
            sum(source[index] * int(n_discr_basis[index, column]) for index in range(7))
            & 1
            for column in range(15)
        ]
        target = [
            sum(source[index] * int(identification[index, column]) for index in range(7))
            & 1
            for column in range(7)
        ]
        q2_n = q2_of_half_vector(n_bits, n_gram)
        q2_t = q2_of_half_vector(target, t_gram)
        if (q2_n + q2_t) % 4:
            raise AssertionError("the graph is not a discriminant anti-isometry")
        q2_histogram[str(q2_n)] += 1
        q_table.append(
            {
                "source": list(source),
                "target": target,
                "q2_source_mod_4": q2_n,
                "q2_target_mod_4": q2_t,
            }
        )
    if q2_histogram != {"0": 28, "1": 36, "2": 36, "3": 28}:
        raise AssertionError("the N discriminant-form histogram changed")

    # Glue N and T along the anti-isometry graph.
    glue_rows = sp.Matrix(
        [
            [*map(int, n_discr_basis.row(index)), *map(int, identification.row(index))]
            for index in range(7)
        ]
    )
    doubled_module_rows = (2 * sp.eye(22)).col_join(glue_rows)
    k3_basis = hermite_normal_form(doubled_module_rows.T).T / 2
    orthogonal_gram = sp.diag(n_gram, t_gram)
    k3_gram = sp.simplify(k3_basis * orthogonal_gram * k3_basis.T)
    orthogonal_action = sp.diag(n_action, t_action)
    k3_action = sp.simplify(k3_basis * orthogonal_action * k3_basis.inv())
    if not (
        all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(value.q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and k3_action**2 == sp.eye(22)
        and sp.factor(k3_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 14
        * (sp.Symbol("lambda") + 1) ** 8
    ):
        raise AssertionError("the glued K3 involution failed")

    # A generic explicit period in T^B.  In the displayed fixed basis let
    # u=(0,1,0), a=(-1,1,0), b=(3,0,1), and
    # v=p*a+(p+1/(4p))*b for p=2^(1/4).  Rational orthogonality to u and v
    # forces pairing with u, a+b, and b to vanish separately; the resulting
    # rational constraint matrix has determinant -32.
    fourth_root_two = sp.root(2, 4)
    u = sp.Matrix([[0, 1, 0]])
    a = sp.Matrix([[-1, 1, 0]])
    b = sp.Matrix([[3, 0, 1]])
    v = fourth_root_two * a + (
        fourth_root_two + 1 / (4 * fourth_root_two)
    ) * b
    period_fixed_coordinates = u + sp.I * v
    period_t = sp.simplify(period_fixed_coordinates * fixed_t_basis)
    rational_fixed_constraints = sp.Matrix.hstack(
        fixed_t_gram * u.T,
        fixed_t_gram * (a + b).T,
        fixed_t_gram * b.T,
    )
    if not (
        (period_t * t_action - period_t).applyfunc(sp.simplify) == sp.zeros(1, 7)
        and sp.simplify((period_t * t_gram * period_t.T)[0]) == 0
        and sp.simplify((period_t * t_gram * sp.conjugate(period_t).T)[0]) == 4
        and rational_fixed_constraints.det() == -32
    ):
        raise AssertionError("the explicit generic fixed period failed")

    # Since the period's rational orthogonal in T is exactly T^-, recover
    # NS(X) as a saturated Smith kernel in the glued K3 lattice.
    t_component = k3_basis[:, 15:22]
    ns_basis = saturated_row_kernel(t_component * (t_action + sp.eye(7)))
    ns_gram = ns_basis * k3_gram * ns_basis.T
    ns_action = induced_row_action(ns_basis, k3_action)
    ns_smith_matrix = smith_normal_form(ns_gram, domain=ZZ)
    ns_smith = [
        abs(int(ns_smith_matrix[index, index])) for index in range(ns_basis.rows)
    ]
    if not (
        ns_basis.rows == 19
        and exact_inertia(ns_gram) == (1, 18, 0)
        and ns_gram.det() == 32
        and ns_smith == [*[1] * 16, 2, 4, 4]
        and ns_action**2 == sp.eye(19)
        and sp.factor(ns_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 11 * (sp.Symbol("lambda") + 1) ** 8
    ):
        raise AssertionError("the special-period Picard lattice failed")

    n_embedding = sp.Matrix.hstack(sp.eye(15), sp.zeros(15, 7)) * k3_basis.inv()
    ns_pivots = ns_basis.rref()[1]
    n_in_ns = n_embedding[:, list(ns_pivots)] * ns_basis[:, list(ns_pivots)].inv()
    if not (
        all(value.q == 1 for value in n_embedding)
        and all(value.q == 1 for value in n_in_ns)
        and n_in_ns * ns_basis == n_embedding
    ):
        raise AssertionError("the ambient N is not integrally contained in NS(X)")

    fixed_ns_coordinates = saturated_row_kernel(ns_action - sp.eye(19))
    fixed_ns_in_k3 = fixed_ns_coordinates * ns_basis
    archived_fixed_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["saturated_fixed_lattice"]["basis_primal_rows"]
        ]
    )
    archived_fixed_in_k3 = archived_fixed_basis * n_embedding
    if not (
        fixed_ns_coordinates.rows == 11
        and same_row_lattice(fixed_ns_in_k3, archived_fixed_in_k3)
    ):
        raise AssertionError("the fixed Picard lattice enlarged")

    # The two algebraic anti-parts are individually D4(-2).
    anti_n_basis = saturated_row_kernel(n_action + sp.eye(15))
    anti_n_gram = anti_n_basis * n_gram * anti_n_basis.T
    d4_change = sp.Matrix(
        [
            [-6, 2, 1, 3],
            [1, 0, 0, -1],
            [3, -1, -1, -1],
            [3, -1, 0, -1],
        ]
    )
    if not (
        abs(int(d4_change.det())) == 1
        and d4_change * (-anti_n_gram / 2) * d4_change.T == -d4_negative
        and anti_t_basis * t_gram * anti_t_basis.T == 2 * d4_negative
    ):
        raise AssertionError("the component D4(-2) identifications failed")

    anti_basis = saturated_row_kernel(k3_action + sp.eye(22))
    anti_gram = anti_basis * k3_gram * anti_basis.T
    half_anti_gram = anti_gram / 2
    anti_smith_matrix = smith_normal_form(anti_gram, domain=ZZ)
    if not (
        anti_basis.rows == 8
        and exact_inertia(anti_gram) == (0, 8, 0)
        and anti_gram.det() == 256
        and [abs(int(anti_smith_matrix[index, index])) for index in range(8)]
        == [2] * 8
        and all(value.q == 1 for value in half_anti_gram)
        and all(half_anti_gram[index, index] % 2 == 0 for index in range(8))
        and half_anti_gram.det() == 1
        and exact_inertia(half_anti_gram) == (0, 8, 0)
    ):
        raise AssertionError("the full anti-invariant lattice is not E8(-2)")

    fixed_cohom_basis = saturated_row_kernel(k3_action - sp.eye(22))
    fixed_cohom_gram = fixed_cohom_basis * k3_gram * fixed_cohom_basis.T
    if not (
        fixed_cohom_basis.rows == 14
        and exact_inertia(fixed_cohom_gram) == (3, 11, 0)
        and fixed_cohom_gram.det() == -256
    ):
        raise AssertionError("the Nikulin fixed cohomology lattice failed")

    fixed_coordinates = sp.Matrix(
        [[sp.Rational(value) for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]]]
    )
    fixed_interior = fixed_coordinates * archived_fixed_basis
    source_pairings = list(fixed_interior * n_gram * wall_rows.T)
    wall_lookup = {tuple(wall_rows.row(index)) for index in range(60)}
    if not (
        fixed_interior * n_action == fixed_interior
        and (fixed_interior * n_gram * fixed_interior.T)[0] == 7360
        and min(source_pairings) == 16
        and all(value > 0 for value in source_pairings)
        and {tuple(wall_rows.row(index) * n_action) for index in range(60)}
        == wall_lookup
    ):
        raise AssertionError("the invariant source-chamber interior failed")

    certificate: dict[str, Any] = {
        "schema": "a15_7_c01_special_k3_symplectic_lift_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "case_id": "A15_7_1_C01",
        "ambient_picard_action": {
            "ambient_invariants": [15, 7, 1],
            "ambient_rank_signature_determinant": [15, [1, 14], 128],
            "order": 2,
            "fixed_rank": 11,
            "action_characteristic_polynomial": "(x-1)^11*(x+1)^4",
            "action_sha256": matrix_digest(n_action),
            "discriminant_half_basis_coordinate_indices_zero_based": discr_coordinates,
            "discriminant_action": [
                [int(value) for value in n_discr_action.row(index)]
                for index in range(7)
            ],
            "anti_invariant_lattice": "D4(-2)",
            "anti_invariant_basis": matrix_payload(anti_n_basis),
            "anti_invariant_gram": matrix_payload(anti_n_gram),
            "d4_unimodular_change": matrix_payload(d4_change),
        },
        "transcendental_glue": {
            "lattice": "U(2)+U(2)+<-2>^3",
            "gram": matrix_payload(t_gram),
            "involution": matrix_payload(t_action),
            "fixed_basis": matrix_payload(fixed_t_basis),
            "fixed_gram": matrix_payload(fixed_t_gram),
            "fixed_rank_signature_determinant": [3, [2, 1], -32],
            "anti_invariant_basis": matrix_payload(anti_t_basis),
            "anti_invariant_lattice": "D4(-2)",
            "anti_isometry_source_to_T_mod_two": matrix_payload(identification),
            "anti_isometry_integer_determinant": int(identification.det()),
            "actions_equivariant": True,
            "all_128_discriminant_classes_checked": True,
            "source_q2_histogram_mod_4": q2_histogram,
            "quadratic_form_table": q_table,
        },
        "k3_lattice_action": {
            "even_unimodular_signature": [3, 19],
            "action_order": 2,
            "action_characteristic_polynomial": "(x-1)^14*(x+1)^8",
            "fixed_lattice_rank_signature_determinant": [14, [3, 11], -256],
            "anti_invariant_lattice": "E8(-2)",
            "anti_invariant_rank": 8,
            "anti_invariant_determinant": 256,
            "anti_invariant_has_no_minus_two_vectors": True,
            "glued_basis_sha256": matrix_digest(k3_basis),
            "glued_gram_sha256": matrix_digest(k3_gram),
            "glued_action_sha256": matrix_digest(k3_action),
        },
        "period_and_picard_lattice": {
            "period_in_displayed_fixed_T_basis": [
                "i*(2*2^(1/4)+3/(4*2^(1/4)))",
                "1+i*2^(1/4)",
                "i*(2^(1/4)+1/(4*2^(1/4)))",
            ],
            "period_eigenvalue": "1",
            "period_square": "0",
            "period_positive_norm": "4",
            "genericity_reason": (
                "with u=(0,1,0), a=(-1,1,0), b=(3,0,1), and p=2^(1/4), "
                "the imaginary part is p*(a+b)+(p^3/8)*b; rational "
                "orthogonality forces pairing with u, a+b, and b to vanish, "
                "and the resulting constraint determinant is -32"
            ),
            "rational_T_orthogonal": "exactly the anti-invariant D4(-2)",
            "picard_rank": 19,
            "picard_signature": [1, 18],
            "picard_determinant": 32,
            "picard_smith_diagonal_absolute": ns_smith,
            "picard_gram_sha256": matrix_digest(ns_gram),
            "ambient_N_is_contained": True,
            "added_algebraic_rank": 4,
            "all_added_algebraic_directions_are_anti_invariant": True,
            "fixed_NS_is_exactly_archived_N_fixed_lattice": True,
            "fixed_NS_rank": 11,
        },
        "torelli_and_chamber": {
            "source_60_wall_chamber_is_preserved": True,
            "archived_strict_fixed_interior_square": "7360",
            "archived_strict_fixed_interior_minimum_source_pairing": "16",
            "kahler_existence_reason": (
                "the anti-invariant lattice is E8(-2), hence contains no "
                "(-2)-class; no K3 root hyperplane contains the fixed positive "
                "cone, so a generic invariant point in the nonempty strict "
                "source-section interior lies in an invariant Kahler chamber"
            ),
            "global_torelli_conclusion": (
                "a special projective Picard-rank-19 K3 surface carries a "
                "genuine symplectic Nikulin involution inducing A15_7_1_C01"
            ),
        },
        "hyperbolic_fixed_section": {
            "dimension": 10,
            "facet_count": 27,
            "finite_volume": True,
            "compact": False,
            "finite_vertex_count": volume["finite_vertex_count"],
            "ideal_vertex_count": volume["ideal_vertex_count"],
            "elliptic_edge_count": volume["elliptic_edge_count"],
            "bad_edge_count": volume["bad_edge_count"],
            "bad_parabolic_neighborhood_count": volume[
                "bad_parabolic_neighborhood_count"
            ],
            "criterion": volume["criterion"],
        },
        "lineage": {
            "same_signed_normalized_gram_as": "D16-D4-C7-H10-27",
            "facet_permutation_one_based": duplicate[
                "source_to_target_facet_permutation_one_based"
            ],
            "abstract_polytope_is_new": False,
            "dimension_is_new": False,
            "stable_lattice_or_maximal_reflection_group_is_new": False,
            "new_scope": (
                "project-level genuine special-K3 symplectic-involution "
                "realization of the already certified abstract chamber"
            ),
        },
        "theorem": (
            "A15_7_1_C01 is the H10/27 fixed divisor/reflection-chamber "
            "section of a genuine symplectic Nikulin involution on a special "
            "Picard-rank-19 K3 surface"
        ),
        "novelty_scope": (
            "new project-level automorphism realization only; the abstract "
            "H10/27 chamber is exactly the previously frozen "
            "D16-D4-C7-H10-27 and the Nikulin automorphism type is classical"
        ),
        "nef_scope_caution": (
            "the source 60-wall chamber is a complete mixed-root reflection "
            "chamber, not asserted to be the literal nef cone of the special "
            "K3 surface"
        ),
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(args.output), "all_checks_pass": True}, indent=2))


if __name__ == "__main__":
    main()
