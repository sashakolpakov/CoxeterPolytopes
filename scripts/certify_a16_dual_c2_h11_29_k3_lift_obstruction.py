#!/usr/bin/env python3
"""Exact K3-lift boundary for the A16-DUAL-C2-H11-29 section.

The transposition of Alexeev's ``N=(16,6,1)`` chamber has a nontrivial
six-bit discriminant action.  This certificate proves that it has no
finite-order lift on a K3 surface with generic Picard lattice exactly N.
It also constructs an explicit equivariant special-period order-four
Hodge lift and certifies why it is not geometric: two algebraic (-2)
classes are negated, so no invariant Kahler class can exist.

This is deliberately a boundary theorem, not a global no-go for every
possible special-period higher-order lift.
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
    / "alexeev_118_dual_h11_29_construction.json"
)
VINBERG = (
    ROOT / "data" / "enumerations" / "alexeev_118_dual_h11_29_vinberg.json"
)
LINEAGE = (
    ROOT / "data" / "enumerations" / "alexeev_118_dual_h11_29_lineage.json"
)
KERNEL_HELPER = (
    ROOT / "scripts" / "certify_a18_c03_c04_special_k3_symplectic_lifts.py"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "a16_dual_c2_h11_29_k3_lift_obstruction.json"
)

EXPECTED_SOURCE_SHA256 = {
    "construction": "53432322ecc6210e047c9b6e9ca101a2b2f86cd5cfef758532f0975bd53f9a99",
    "vinberg": "f711b1869f9106ae0c4c6d48673ff8b3de84495dc8ac72a03787e0434e9497d9",
    "lineage": "c3325df813b0c74cca4e380f21035da5bc60b054743b2de301439294f498cf31",
    "kernel_helper": "1c15dd20d9fbb2754b154570ef453423fe45961cd6fb263d4d703f40adc6b598",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(i)] for i in range(matrix.rows)]


def matrix_digest(matrix: sp.Matrix) -> str:
    raw = json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    work = [[int(value) & 1 for value in row] for row in rows]
    rank = 0
    if not work:
        return rank
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


def gf2_coordinates(
    vector: Sequence[int], basis: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    for coefficients in itertools.product((0, 1), repeat=len(basis)):
        candidate = tuple(
            sum(
                coefficients[index] * basis[index][column]
                for index in range(len(basis))
            )
            & 1
            for column in range(len(vector))
        )
        if candidate == tuple(int(value) & 1 for value in vector):
            return coefficients
    raise AssertionError("the vector left the displayed GF(2) span")


def two_elementary_basis(
    gram: sp.Matrix, expected_length: int
) -> tuple[tuple[int, ...], ...]:
    """Return bit rows v for all discriminant generators v/2."""

    size = gram.rows
    work = [[int(gram[i, j]) & 1 for j in range(size)] for i in range(size)]
    pivots: list[int] = []
    rank = 0
    for column in range(size):
        pivot = next(
            (index for index in range(rank, size) if work[index][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(size):
            if index != rank and work[index][column]:
                work[index] = [
                    left ^ right for left, right in zip(work[index], work[rank])
                ]
        pivots.append(column)
        rank += 1
    free = [column for column in range(size) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [0] * size
        vector[free_column] = 1
        for row, pivot_column in reversed(list(enumerate(pivots))):
            vector[pivot_column] = sum(
                work[row][column] * vector[column] for column in free
            ) & 1
        basis.append(tuple(vector))
    if len(basis) != expected_length or gf2_rank(basis) != expected_length:
        raise AssertionError("unexpected two-elementary discriminant length")
    return tuple(basis)


def discriminant_action(
    action: sp.Matrix, basis: Sequence[Sequence[int]]
) -> sp.Matrix:
    rows = []
    for vector in basis:
        image = tuple(
            sum(vector[index] * int(action[index, column]) for index in range(action.rows))
            & 1
            for column in range(action.cols)
        )
        rows.append(gf2_coordinates(image, basis))
    return sp.Matrix(rows)


def q_times_four_mod_eight(
    coordinates: Sequence[int],
    half_vector_basis: Sequence[Sequence[int]],
    gram: sp.Matrix,
) -> int:
    bits = tuple(
        sum(
            coordinates[index] * half_vector_basis[index][column]
            for index in range(len(half_vector_basis))
        )
        & 1
        for column in range(gram.rows)
    )
    return sum(
        bits[i] * int(gram[i, j]) * bits[j]
        for i in range(gram.rows)
        for j in range(gram.cols)
    ) % 8


def polynomial_mod_two(poly: sp.Expr, variable: sp.Symbol) -> list[int]:
    polynomial = sp.Poly(sp.expand(poly), variable)
    return [int(value) & 1 for value in polynomial.all_coeffs()]


def build_certificate() -> dict[str, Any]:
    source_hashes = {
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
        "kernel_helper": sha256(KERNEL_HELPER),
    }
    if source_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned H11/29 source or exact helper changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    case = next(row for row in construction["cases"] if row["id"] == "A16_6_1_C01")
    volume = next(row for row in vinberg["audits"] if row["id"] == "A16_6_1_C01")
    if not (
        construction["all_checks_pass"]
        and vinberg["all_checks_pass"]
        and lineage["all_checks_pass"]
        and all(case["construction_checks"].values())
        and case["case"]["subgroup_order"] == 2
        and case["case"]["fixed_rank"] == 12
        and case["exact_reduction"]["active_hyperbolic_facet_count"] == 29
        and volume["finite_volume"] is True
        and volume["compact"] is False
        and volume["bad_edge_count"] == 0
        and volume["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the upstream H11/29 target changed")

    # Alexeev's leading-root basis for the (16,6,1) lattice.  Its action is
    # read from the independently frozen complete wall-group certificate.
    n_gram = sp.Matrix(
        [
            [-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
            [0, -4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 2],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 2, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0],
            [2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -4],
        ]
    )
    n_action = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["case"]["representative_nonidentity_primal_actions"][0]
        ]
    )
    if not (
        n_gram.det() == -64
        and exact_inertia(n_gram) == (1, 15, 0)
        and n_action**2 == sp.eye(16)
        and n_action * n_gram * n_action.T == n_gram
        and sp.factor(n_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 12
        * (sp.Symbol("lambda") + 1) ** 4
    ):
        raise AssertionError("the ambient (16,6,1) transposition changed")

    n_discriminant_basis = two_elementary_basis(n_gram, 6)
    n_discriminant_action = discriminant_action(n_action, n_discriminant_basis)
    target_minus_identity = (n_discriminant_action - sp.eye(6)).applyfunc(
        lambda value: int(value) & 1
    )
    if not (
        (n_discriminant_action**2 - sp.eye(6)).applyfunc(
            lambda value: int(value) & 1
        )
        == sp.zeros(6)
        and gf2_rank(target_minus_identity.tolist()) == 2
    ):
        raise AssertionError("the six-bit Picard discriminant action changed")

    # Generic-Picard obstruction.  If NS(X)=N, every rational primary
    # summand of T not containing the period is algebraic.  Hence a finite
    # order lift must have one cyclotomic primary type on all six dimensions.
    # The elementary inverse-totient list phi(m)|6 is exactly the list below.
    x = sp.Symbol("x")
    cyclotomic_orders = [1, 2, 3, 4, 6, 7, 9, 14, 18]
    expected_inverse_totient = {
        "1": [1, 2],
        "2": [3, 4, 6],
        "6": [7, 9, 14, 18],
    }
    if {
        str(degree): [m for m in cyclotomic_orders if sp.totient(m) == degree]
        for degree in (1, 2, 6)
    } != expected_inverse_totient:
        raise AssertionError("the inverse-totient denominator changed")

    target_char_mod_two = polynomial_mod_two((x + 1) ** 6, x)
    generic_rows = []
    for order in cyclotomic_orders:
        degree = int(sp.totient(order))
        copies = 6 // degree
        characteristic = sp.expand(sp.cyclotomic_poly(order, x) ** copies)
        reduction = polynomial_mod_two(characteristic, x)
        if order in (1, 2):
            reason = "the only action is +/-identity, hence trivial on T/2T"
        elif order == 4:
            # det(B-I)=Phi_4(1)^3=8.  Rank two modulo two would give four
            # even Smith factors and therefore force divisibility by 16.
            determinant_b_minus_i = int(sp.cyclotomic_poly(order, 1) ** copies)
            if determinant_b_minus_i != 8 or determinant_b_minus_i % 16 == 0:
                raise AssertionError("the order-four Smith obstruction changed")
            reason = (
                "rank(B-I mod 2)=2 would force four even Smith factors, "
                "so 16 divides det(B-I), but det(B-I)=8"
            )
        else:
            if reduction == target_char_mod_two:
                raise AssertionError("a non-order-four cyclotomic type survived")
            reason = "the characteristic polynomial modulo two is not (x+1)^6"
        generic_rows.append(
            {
                "order": order,
                "cyclotomic_degree": degree,
                "multiplicity": copies,
                "characteristic_polynomial": str(sp.factor(characteristic)),
                "characteristic_coefficients_mod_two": reduction,
                "excluded_reason": reason,
            }
        )

    # The closest special-period lift.  T has the unique 2-elementary type
    # (signature (2,4), length 6, delta 1) in the displayed model.
    t_gram = sp.Matrix(
        [
            [0, 2, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0],
            [0, 0, 0, 2, 0, 0],
            [0, 0, 2, 0, 0, 0],
            [0, 0, 0, 0, -2, 0],
            [0, 0, 0, 0, 0, -2],
        ]
    )
    order_four_block = sp.Matrix(
        [
            [0, 0, -1, 0],
            [0, 0, 0, -1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    t_action = sp.diag(order_four_block, -sp.eye(2))
    identification = sp.Matrix(
        [
            [1, 0, 1, 0, 1, 1],
            [0, 1, 0, 1, 1, 1],
            [0, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1],
            [0, 1, 1, 0, 0, 1],
            [0, 0, 0, 1, 1, 1],
        ]
    )
    t_half_basis = tuple(
        tuple(int(index == column) for column in range(6)) for index in range(6)
    )
    t_mod_two = t_action.applyfunc(lambda value: int(value) & 1)
    if not (
        t_gram.det() == 64
        and exact_inertia(t_gram) == (2, 4, 0)
        and t_action**4 == sp.eye(6)
        and t_action**2 != sp.eye(6)
        and t_action * t_gram * t_action.T == t_gram
        and int(identification.det()) % 2 == 1
        and (n_discriminant_action * identification).applyfunc(
            lambda value: int(value) & 1
        )
        == (identification * t_mod_two).applyfunc(lambda value: int(value) & 1)
    ):
        raise AssertionError("the special order-four transcendental data failed")

    q_table = []
    q_histogram = {str(value): 0 for value in range(8)}
    for source in itertools.product((0, 1), repeat=6):
        target = tuple(
            int(value) & 1
            for value in (sp.Matrix([source]) * identification).row(0)
        )
        q_n = q_times_four_mod_eight(source, n_discriminant_basis, n_gram)
        q_t = q_times_four_mod_eight(target, t_half_basis, t_gram)
        if (q_n + q_t) % 8:
            raise AssertionError("the special graph is not an anti-isometry")
        q_histogram[str(q_n)] += 1
        q_table.append(
            {
                "source": list(source),
                "target": list(target),
                "4q_N_mod_8": q_n,
                "4q_T_mod_8": q_t,
            }
        )

    glue_rows = sp.Matrix(
        [
            [
                *n_discriminant_basis[index],
                *map(int, identification.row(index)),
            ]
            for index in range(6)
        ]
    )
    k3_basis = hermite_normal_form(
        ((2 * sp.eye(22)).col_join(glue_rows)).T
    ).T / 2
    k3_gram = sp.simplify(k3_basis * sp.diag(n_gram, t_gram) * k3_basis.T)
    k3_action = sp.simplify(
        k3_basis * sp.diag(n_action, t_action) * k3_basis.inv()
    )
    if not (
        all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(value.q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and k3_action**4 == sp.eye(22)
        and sp.factor(k3_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 12
        * (sp.Symbol("lambda") + 1) ** 6
        * (sp.Symbol("lambda") ** 2 + 1) ** 2
    ):
        raise AssertionError("the special graph glue is not a K3 lattice action")

    sqrt_two = sp.sqrt(2)
    period_t = sp.Matrix([[1, sqrt_two, sp.I, sp.I * sqrt_two, 0, 0]])
    if not (
        period_t * t_action == sp.I * period_t
        and sp.simplify((period_t * t_gram * period_t.T)[0]) == 0
        and sp.simplify(
            (period_t * t_gram * sp.conjugate(period_t).T)[0]
        )
        == 8 * sqrt_two
    ):
        raise AssertionError("the positive special eigenperiod failed")

    # Rational orthogonality to the displayed period kills the first four T
    # coordinates and leaves exactly the last two (-1)-eigenvectors.
    t_component = k3_basis[:, 16:22]
    ns_basis = saturated_row_kernel(t_component * (t_action + sp.eye(6)))
    ns_gram = ns_basis * k3_gram * ns_basis.T
    ns_action = induced_row_action(ns_basis, k3_action)
    ns_smith_matrix = smith_normal_form(ns_gram, domain=ZZ)
    ns_smith = [
        abs(int(ns_smith_matrix[index, index])) for index in range(ns_basis.rows)
    ]
    if not (
        ns_basis.rows == 18
        and exact_inertia(ns_gram) == (1, 17, 0)
        and ns_gram.det() == -16
        and ns_smith == [*[1] * 14, 2, 2, 2, 2]
        and sp.factor(ns_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 12
        * (sp.Symbol("lambda") + 1) ** 6
    ):
        raise AssertionError("the special Picard lattice changed")

    n_embedding = sp.Matrix.hstack(sp.eye(16), sp.zeros(16, 6)) * k3_basis.inv()
    ns_pivots = ns_basis.rref()[1]
    n_in_ns = (
        n_embedding[:, list(ns_pivots)] * ns_basis[:, list(ns_pivots)].inv()
    )
    if not (
        all(value.q == 1 for value in n_embedding)
        and all(value.q == 1 for value in n_in_ns)
        and n_in_ns * ns_basis == n_embedding
    ):
        raise AssertionError("N is not integrally contained in the special NS")

    fixed_ns_coordinates = saturated_row_kernel(ns_action - sp.eye(18))
    fixed_ns_in_k3 = fixed_ns_coordinates * ns_basis
    archived_fixed_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["saturated_fixed_lattice"]["basis_primal_rows"]
        ]
    )
    archived_fixed_in_k3 = archived_fixed_basis * n_embedding
    if not (
        fixed_ns_coordinates.rows == 12
        and same_row_lattice(fixed_ns_in_k3, archived_fixed_in_k3)
    ):
        raise AssertionError("the special fixed NS lattice enlarged")

    coinvariant_coordinates = saturated_row_kernel(
        ns_gram * fixed_ns_coordinates.T
    )
    coinvariant_in_k3 = coinvariant_coordinates * ns_basis
    coinvariant_gram = coinvariant_in_k3 * k3_gram * coinvariant_in_k3.T
    coinvariant_smith_matrix = smith_normal_form(coinvariant_gram, domain=ZZ)
    coinvariant_smith = [
        abs(int(coinvariant_smith_matrix[index, index]))
        for index in range(coinvariant_gram.rows)
    ]

    t_embedding = sp.Matrix.hstack(sp.zeros(6, 16), sp.eye(6)) * k3_basis.inv()
    obstruction_roots = t_embedding[4:6, :]
    root_gram = obstruction_roots * k3_gram * obstruction_roots.T
    if not (
        coinvariant_coordinates.rows == 6
        and exact_inertia(coinvariant_gram) == (0, 6, 0)
        and coinvariant_gram.det() == 256
        and coinvariant_smith == [2, 2, 2, 2, 4, 4]
        and all(value.q == 1 for value in t_embedding)
        and root_gram == -2 * sp.eye(2)
        and obstruction_roots * k3_action == -obstruction_roots
        and obstruction_roots * k3_gram * fixed_ns_in_k3.T
        == sp.zeros(2, 12)
    ):
        raise AssertionError("the exact Kahler-root obstruction changed")

    certificate: dict[str, Any] = {
        "schema": "a16_dual_c2_h11_29_k3_lift_obstruction_v1",
        "all_checks_pass": True,
        "source_file_sha256": source_hashes,
        "case_id": "A16_6_1_C01",
        "ambient_picard_action": {
            "invariants": [16, 6, 1],
            "rank_signature_determinant": [16, [1, 15], -64],
            "order": 2,
            "fixed_rank": 12,
            "characteristic_polynomial": "(x-1)^12*(x+1)^4",
            "gram_sha256": matrix_digest(n_gram),
            "action_sha256": matrix_digest(n_action),
            "discriminant_basis": [list(row) for row in n_discriminant_basis],
            "discriminant_action": [
                [int(value) for value in n_discriminant_action.row(index)]
                for index in range(6)
            ],
            "rank_of_discriminant_action_minus_identity_mod_two": 2,
        },
        "generic_picard_no_go": {
            "statement": (
                "no finite-order K3 automorphism inducing this Picard action "
                "can have NS(X) exactly N=(16,6,1)"
            ),
            "single_primary_reason": (
                "for NS(X)=N, any rational cyclotomic primary summand of T "
                "not containing the eigenperiod would be algebraic"
            ),
            "inverse_totient_lemma": expected_inverse_totient,
            "target_characteristic_coefficients_mod_two": target_char_mod_two,
            "target_rank_B_minus_I_mod_two": 2,
            "cyclotomic_exhaustion": generic_rows,
            "order_four_determinant_obstruction": {
                "determinant_B_minus_I": 8,
                "required_even_smith_factor_count": 4,
                "forced_divisor": 16,
                "contradiction": True,
            },
        },
        "explicit_special_order_four_near_miss": {
            "transcendental_lattice": "U(2)+U(2)+<-2>+<-2>",
            "transcendental_gram": matrix_payload(t_gram),
            "transcendental_action": matrix_payload(t_action),
            "transcendental_characteristic_polynomial": "(x+1)^2*(x^2+1)^2",
            "anti_isometry_source_to_T_mod_two": matrix_payload(identification),
            "anti_isometry_integer_determinant": int(identification.det()),
            "all_64_discriminant_classes_checked": True,
            "source_4q_histogram_mod_8": q_histogram,
            "quadratic_form_table": q_table,
            "k3_even_unimodular_signature": [3, 19],
            "k3_action_order": 4,
            "k3_action_characteristic_polynomial": "(x-1)^12*(x+1)^6*(x^2+1)^2",
            "k3_basis_sha256": matrix_digest(k3_basis),
            "k3_gram_sha256": matrix_digest(k3_gram),
            "k3_action_sha256": matrix_digest(k3_action),
            "period": ["1", "sqrt(2)", "i", "i*sqrt(2)", "0", "0"],
            "period_eigenvalue": "i",
            "period_square": "0",
            "period_positive_norm": "8*sqrt(2)",
            "rational_period_orthogonal_in_T": "the last two <-2> summands",
            "picard_rank_signature_determinant": [18, [1, 17], -16],
            "picard_smith_diagonal_absolute": ns_smith,
            "fixed_NS_is_exactly_archived_fixed_lattice": True,
            "fixed_NS_rank": 12,
            "coinvariant_rank_signature_determinant": [6, [0, 6], 256],
            "coinvariant_smith_diagonal_absolute": coinvariant_smith,
            "obstruction_root_gram": matrix_payload(root_gram),
            "obstruction_root_coordinates_in_k3_basis": matrix_payload(
                obstruction_roots
            ),
            "both_roots_are_negated": True,
            "kahler_gate_passes": False,
            "failure_reason": (
                "a finite-order automorphism has an invariant Kahler class "
                "obtained by averaging, but every fixed class is orthogonal "
                "to the two displayed negated (-2)-classes"
            ),
        },
        "hyperbolic_section_handoff": {
            "dimension": 11,
            "facet_count": 29,
            "finite_volume": True,
            "compact": False,
            "elliptic_edge_count": volume["elliptic_edge_count"],
            "finite_vertex_count": volume["finite_vertex_count"],
            "ideal_vertex_count": volume["ideal_vertex_count"],
            "bad_edge_count": volume["bad_edge_count"],
            "bad_parabolic_neighborhood_count": volume[
                "bad_parabolic_neighborhood_count"
            ],
        },
        "theorem_scope": (
            "exact generic-Picard finite-order nonexistence plus an exact "
            "explicit special-period arithmetic lift rejected by the Kahler "
            "gate; no global no-go for every special-period finite-order lift"
        ),
        "nef_scope_caution": (
            "the source 118-wall mixed reflection chamber is not asserted "
            "to be the literal nef cone of a special K3 surface"
        ),
    }
    raw = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(raw).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("the archived A16 K3-lift obstruction changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "all_checks_pass": True,
                "payload_sha256": payload["certificate_payload_sha256"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
