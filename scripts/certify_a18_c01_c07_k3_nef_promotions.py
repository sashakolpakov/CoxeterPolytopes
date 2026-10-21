#!/usr/bin/env python3
"""Exact K3 and nef-cone promotions of the A18 C01 and C07 sections.

The factory certificate treats the T2=C01 and T3=C07 sections only as
fixed sections of the complete 48-wall ``W_r`` chamber of the
2-elementary lattice ``N=(18,4,0)``.  This replay supplies the missing
geometric gates.

* C01 acts trivially on the discriminant form and is glued to the identity
  on ``T=U(2)+U(2)``.  Its N-coinvariant lattice is ``E8(-2)``.
* C07 has a nontrivial order-three discriminant action.  An explicit
  order-three isometry of ``T`` realizes that action.  Its N-coinvariant
  lattice is ``-(A2 tensor D4)`` and therefore has minimum norm -4.

Both glued actions preserve a generic period with exact Picard lattice N
and preserve a K3 Weyl chamber.  Because every active restricted factory
wall comes from a short (-2)-wall and every long (-4)-wall restricts
redundantly, the archived fixed cone is exactly the invariant nef section,
not merely a mixed reflection-chamber carrier.  Finite quotient pullback
then identifies it with the nef cone of the normal quotient surface.
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
from scripts.certify_a18_c06_h13_31_k3_order4_lift import (
    discriminant_action,
    gf2_rank,
    matrix_digest,
    matrix_payload,
    q_n_mod_two,
    q_t_mod_two,
    two_elementary_discriminant_basis,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
OUTPUT = (
    ROOT / "data" / "enumerations" / "a18_c01_c07_k3_nef_promotions.json"
)
C03_HELPER = ROOT / "scripts" / "certify_a18_c03_c04_special_k3_symplectic_lifts.py"
C06_HELPER = ROOT / "scripts" / "certify_a18_c06_h13_31_k3_order4_lift.py"

EXPECTED_SOURCE_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "c03_helper": "1c15dd20d9fbb2754b154570ef453423fe45961cd6fb263d4d703f40adc6b598",
    "c06_helper": "967d6e258a00df16ff44ad3f762f81817676152d882be2fb3088894e87562027",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def smith_diagonal(gram: sp.Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return [abs(int(smith[index, index])) for index in range(gram.rows)]


def gf2_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(value) % 2 for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def gf2_product(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    return (left * right).applyfunc(lambda value: int(value) % 2)


def q_table_and_equivariance(
    n_action: sp.Matrix,
    n_gram: sp.Matrix,
    discr_basis: Sequence[Sequence[int]],
    t_action: sp.Matrix,
    identification: sp.Matrix,
) -> tuple[sp.Matrix, list[dict[str, object]]]:
    n_discr_action = discriminant_action(n_action, discr_basis)
    if gf2_product(n_discr_action, identification) != gf2_product(
        identification, t_action
    ):
        raise AssertionError("the N and T discriminant actions are not equivariant")

    table = []
    for coordinates in itertools.product(range(2), repeat=4):
        target = tuple(
            int(value) % 2
            for value in (sp.Matrix([coordinates]) * identification).row(0)
        )
        q_n = q_n_mod_two(coordinates, discr_basis, n_gram)
        q_t = q_t_mod_two(target)
        if q_n != q_t:
            raise AssertionError("the displayed graph is not a discriminant anti-isometry")
        table.append(
            {"source": list(coordinates), "target": list(target), "q_mod_2": q_n}
        )
    return n_discr_action, table


def glue_k3_lattice(
    n_gram: sp.Matrix,
    n_action: sp.Matrix,
    discr_basis: Sequence[Sequence[int]],
    t_gram: sp.Matrix,
    t_action: sp.Matrix,
    identification: sp.Matrix,
    order: int,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    glue_rows = sp.Matrix(
        [
            [*discr_basis[index], *map(int, identification.row(index))]
            for index in range(4)
        ]
    )
    doubled_module_rows = (2 * sp.eye(22)).col_join(glue_rows)
    k3_basis = hermite_normal_form(doubled_module_rows.T).T / 2
    orthogonal_gram = sp.diag(n_gram, t_gram)
    k3_gram = sp.simplify(k3_basis * orthogonal_gram * k3_basis.T)
    orthogonal_action = sp.diag(n_action, t_action)
    k3_action = sp.simplify(k3_basis * orthogonal_action * k3_basis.inv())
    if not (
        k3_basis.det() == sp.Rational(1, 16)
        and all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(value.q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and k3_action**order == sp.eye(22)
    ):
        raise AssertionError("the graph glue did not produce the asserted K3 action")
    return k3_basis, k3_gram, k3_action


def primitive_ns_in_glued_lattice(
    k3_basis: sp.Matrix, k3_gram: sp.Matrix
) -> tuple[sp.Matrix, sp.Matrix]:
    """Recover L_K3 intersect N_Q and prove that it is exactly N."""

    ns_basis = saturated_row_kernel(k3_basis[:, 18:22])
    ns_gram = ns_basis * k3_gram * ns_basis.T
    if not (ns_basis.rows == 18 and exact_inertia(ns_gram) == (1, 17, 0)):
        raise AssertionError("the generic-period Picard kernel has the wrong rank")
    return ns_basis, ns_gram


def a2_tensor_d4_witness(positive_coin_gram: sp.Matrix) -> dict[str, object]:
    """Certify an integral isometry to A2 tensor D4 and its minimum four."""

    a2 = sp.Matrix([[2, -1], [-1, 2]])
    d4 = sp.Matrix(
        [
            [2, -1, 0, 0],
            [-1, 2, -1, -1],
            [0, -1, 2, 0],
            [0, -1, 0, 2],
        ]
    )
    tensor_gram = sp.kronecker_product(a2, d4)
    # PARI qfisom found this matrix; the equality below is replayed entirely
    # in exact integer arithmetic and is the certificate, not the search.
    change = sp.Matrix(
        [
            [-4, 0, 0, 0, 1, 0, 0, 1],
            [-4, 0, 0, 0, 2, -1, 0, 2],
            [-2, 0, 0, 1, -2, 0, 0, 1],
            [-2, 0, 0, 0, 1, 0, -1, 2],
            [-2, 0, 0, 0, 1, 0, 0, 0],
            [-6, 0, 1, 0, 0, 0, 0, 0],
            [-3, 0, 0, 1, -1, 0, 0, 0],
            [-3, 1, 0, 0, 0, 0, 0, 0],
        ]
    )
    if not (
        change.det() == -1
        and change.T * tensor_gram * change == positive_coin_gram
        and tensor_gram.det() == 1296
        and exact_inertia(d4) == (4, 0, 0)
        and all(d4[index, index] % 2 == 0 for index in range(4))
    ):
        raise AssertionError("the A2 tensor D4 coinvariant identification failed")

    # For x,y in D4, the tensor norm is
    #   x^2 + (x-y)^2 + y^2.
    # Every nonzero D4 vector has even norm at least 2, and at least two of
    # x, y, x-y are nonzero.  Hence the tensor minimum is at least 4; a
    # decomposable pair of roots attains 4.
    root = sp.Matrix([[1, 0, 0, 0]])
    attaining = sp.Matrix([[1, 0, 0, 0, 0, 0, 0, 0]])
    if not (
        (root * d4 * root.T)[0] == 2
        and (attaining * tensor_gram * attaining.T)[0] == 4
    ):
        raise AssertionError("the tensor minimum witness changed")
    return {
        "isometry_type": "-(A2 tensor D4)",
        "positive_tensor_gram": matrix_payload(tensor_gram),
        "unimodular_change": matrix_payload(change),
        "congruence": "change^T*(A2 tensor D4)*change = -coinvariant_gram",
        "minimum_positive_norm": 4,
        "rootlessness_proof": (
            "write a tensor vector as (x,y) in D4^2; its norm is "
            "x^2+(x-y)^2+y^2.  At least two terms are nonzero and every "
            "nonzero D4 norm is at least 2"
        ),
    }


def build_certificate() -> dict[str, Any]:
    actual_hashes = {
        "construction": sha256(CONSTRUCTION),
        "factory": sha256(FACTORY),
        "c03_helper": sha256(C03_HELPER),
        "c06_helper": sha256(C06_HELPER),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A18 factory source or exact helper changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    if not (
        construction["all_checks_pass"]
        and factory["summary"]["raw_marked_subgroup_row_count"] == 5
        and factory["summary"]["finite_volume_noncompact_count"] == 5
        and factory["summary"]["total_bad_edge_count"] == 0
        and factory["summary"]["total_bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("an upstream A18/factory certificate is not green")

    inventory = next(
        row
        for row in construction["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    n_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in inventory["ambient"]["gram"]]
    )
    wall_rows = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in inventory["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    if not (
        n_gram.det() == -16
        and exact_inertia(n_gram) == (1, 17, 0)
        and wall_rows.rows == 48
        and sorted((wall_rows * n_gram * wall_rows.T).diagonal())
        == [-4] * 24 + [-2] * 24
    ):
        raise AssertionError("the ambient (18,4,0) lattice or 48 walls changed")

    discr_basis = two_elementary_discriminant_basis(n_gram)
    if gf2_rank(discr_basis) != 4:
        raise AssertionError("the ambient discriminant basis is not independent")
    t_gram = sp.Matrix(
        [
            [0, 2, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 0, 2],
            [0, 0, 2, 0],
        ]
    )
    identification = sp.Matrix(
        [
            [0, 0, 1, 1],
            [0, 0, 0, 1],
            [1, 1, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    if not (
        t_gram.det() == 16
        and exact_inertia(t_gram) == (2, 2, 0)
        and identification.det() == 1
    ):
        raise AssertionError("the transcendental lattice or glue changed")

    zeta = (-1 + sp.I * sp.sqrt(3)) / 2
    period = sp.Matrix(
        [[sp.sqrt(2) * zeta**2, -zeta, 1, sp.sqrt(2)]]
    )
    rational_pairing_constraints = sp.Matrix(
        [
            [1, 0, 0, 2],
            [0, -1, 2, 0],
            [-1, 0, 0, 0],
            [0, -1, 0, 0],
        ]
    )
    if not (
        sp.simplify((period * t_gram * period.T)[0]) == 0
        and sp.simplify((period * t_gram * sp.conjugate(period).T)[0])
        == 6 * sp.sqrt(2)
        and rational_pairing_constraints.det() != 0
    ):
        raise AssertionError("the generic positive period failed")

    configs = {
        "A18_4_0_C01": {
            "order": 2,
            "diagram_type": "T2",
            "t_action": sp.eye(4),
            "period_eigenvalue": sp.Integer(1),
            "period_eigenvalue_text": "1",
            "automorphism_character": "symplectic Nikulin involution",
            "fixed_gram_det": -256,
            "fixed_smith": [1, 1, *([2] * 8)],
            "coin_det": 256,
            "coin_smith": [2] * 8,
            "index_fixed_plus_coin": 64,
            "folded_norm_histogram": {"-2": 4, "-4": 10},
            "finite_vertices": 90,
            "ideal_vertices": 20,
        },
        "A18_4_0_C07": {
            "order": 3,
            "diagram_type": "T3",
            "t_action": sp.Matrix(
                [
                    [0, 0, 0, -1],
                    [0, -1, -1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, -1],
                ]
            ),
            "period_eigenvalue": zeta,
            "period_eigenvalue_text": "zeta_3=(-1+i*sqrt(3))/2",
            "automorphism_character": "purely non-symplectic automorphism of order 3",
            "fixed_gram_det": -81,
            "fixed_smith": [*([1] * 6), *([3] * 4)],
            "coin_det": 1296,
            "coin_smith": [*([1] * 4), *([6] * 4)],
            "index_fixed_plus_coin": 81,
            "folded_norm_histogram": {"-2": 9, "-6": 5},
            "finite_vertices": 128,
            "ideal_vertices": 6,
        },
    }

    records = []
    for case_id, config in configs.items():
        factory_case = next(
            row for row in factory["constructions"] if row["id"] == case_id
        )
        section = factory_case["section"]
        exact_cone = factory_case["exact_cone"]
        order = int(config["order"])
        n_action = sp.Matrix(
            [
                [sp.Rational(value) for value in row]
                for row in section["case"][
                    "representative_nonidentity_primal_actions"
                ][0]
            ]
        )
        t_action: sp.Matrix = config["t_action"]
        if not (
            factory_case["diagram_type"] == config["diagram_type"]
            and section["case"]["subgroup_order"] == order
            and section["case"]["fixed_rank"] == 10
            and section["case"]["hyperbolic_dimension"] == 9
            and n_action**order == sp.eye(18)
            and n_action * n_gram * n_action.T == n_gram
            and t_action**order == sp.eye(4)
            and t_action * t_gram * t_action.T == t_gram
            and sp.simplify(period * t_action - config["period_eigenvalue"] * period)
            == sp.zeros(1, 4)
        ):
            raise AssertionError(f"the archived action or period failed for {case_id}")

        wall_lookup = {tuple(wall_rows.row(index)) for index in range(48)}
        if {
            tuple(wall_rows.row(index) * n_action) for index in range(48)
        } != wall_lookup:
            raise AssertionError(f"the source chamber is not preserved for {case_id}")

        discr_action_matrix, q_table = q_table_and_equivariance(
            n_action,
            n_gram,
            discr_basis,
            t_action,
            identification,
        )
        if not (
            gf2_product(discr_action_matrix**order, sp.eye(4)) == sp.eye(4)
            and (
                case_id != "A18_4_0_C01"
                or discr_action_matrix == sp.eye(4)
            )
            and (
                case_id != "A18_4_0_C07"
                or gf2_rank(gf2_payload(discr_action_matrix - sp.eye(4))) == 4
            )
        ):
            raise AssertionError(f"the discriminant action failed for {case_id}")

        k3_basis, k3_gram, k3_action = glue_k3_lattice(
            n_gram,
            n_action,
            discr_basis,
            t_gram,
            t_action,
            identification,
            order,
        )
        ns_basis, ns_gram = primitive_ns_in_glued_lattice(k3_basis, k3_gram)
        n_embedding = (
            sp.Matrix.hstack(sp.eye(18), sp.zeros(18, 4)) * k3_basis.inv()
        )
        if not all(value.q == 1 for value in n_embedding):
            raise AssertionError("N is not integrally embedded in the K3 lattice")
        if not (
            same_row_lattice(ns_basis, n_embedding)
            and ns_gram.det() == -16
            and smith_diagonal(ns_gram) == [*([1] * 14), *([2] * 4)]
        ):
            raise AssertionError(f"the generic period does not have NS=N for {case_id}")

        fixed_basis = saturated_row_kernel(n_action - sp.eye(18))
        norm_operator = sum((n_action**power for power in range(order)), sp.zeros(18))
        coin_basis = saturated_row_kernel(norm_operator)
        fixed_gram = fixed_basis * n_gram * fixed_basis.T
        coin_gram = coin_basis * n_gram * coin_basis.T
        coin_action = induced_row_action(coin_basis, n_action)
        archived_fixed_basis = sp.Matrix(
            [
                [sp.Rational(value) for value in row]
                for row in section["saturated_fixed_lattice"]["basis_primal_rows"]
            ]
        )
        direct_sum_index = abs(int(fixed_basis.col_join(coin_basis).det()))
        if not (
            fixed_basis.rows == 10
            and same_row_lattice(fixed_basis, archived_fixed_basis)
            and fixed_gram.det() == config["fixed_gram_det"]
            and smith_diagonal(fixed_gram) == config["fixed_smith"]
            and exact_inertia(fixed_gram) == (1, 9, 0)
            and coin_basis.rows == 8
            and coin_gram.det() == config["coin_det"]
            and smith_diagonal(coin_gram) == config["coin_smith"]
            and exact_inertia(coin_gram) == (0, 8, 0)
            and direct_sum_index == config["index_fixed_plus_coin"]
        ):
            raise AssertionError(f"the fixed/coinvariant splitting failed for {case_id}")

        if case_id == "A18_4_0_C01":
            half_positive = -coin_gram / 2
            if not (
                coin_action == -sp.eye(8)
                and all(value.q == 1 for value in half_positive)
                and all(half_positive[index, index] % 2 == 0 for index in range(8))
                and half_positive.det() == 1
                and exact_inertia(half_positive) == (8, 0, 0)
            ):
                raise AssertionError("the C01 coinvariant is not E8(-2)")
            coin_rootlessness: dict[str, object] = {
                "isometry_type": "E8(-2)",
                "minimum_negative_norm_absolute": 4,
                "rootlessness_proof": (
                    "(-coinvariant_gram)/2 is positive even unimodular of rank 8; "
                    "hence it is E8 and every nonzero norm is at least 2"
                ),
            }
        else:
            coin_rootlessness = a2_tensor_d4_witness(-coin_gram)
            if sp.factor(coin_action.charpoly().as_expr()) != (
                sp.Symbol("lambda") ** 2 + sp.Symbol("lambda") + 1
            ) ** 4:
                raise AssertionError("the C07 coinvariant action changed")

        fixed_cohom_basis = saturated_row_kernel(k3_action - sp.eye(22))
        cohom_norm_operator = sum(
            (k3_action**power for power in range(order)), sp.zeros(22)
        )
        coin_cohom_basis = saturated_row_kernel(cohom_norm_operator)
        fixed_cohom_gram = fixed_cohom_basis * k3_gram * fixed_cohom_basis.T
        coin_cohom_gram = coin_cohom_basis * k3_gram * coin_cohom_basis.T
        if case_id == "A18_4_0_C01":
            expected_cohom = {
                "charpoly": (sp.Symbol("lambda") - 1) ** 14
                * (sp.Symbol("lambda") + 1) ** 8,
                "fixed": (14, (3, 11, 0), -256, [*([1] * 6), *([2] * 8)]),
                "coin": (8, (0, 8, 0), 256, [2] * 8),
            }
        else:
            expected_cohom = {
                "charpoly": (sp.Symbol("lambda") - 1) ** 10
                * (sp.Symbol("lambda") ** 2 + sp.Symbol("lambda") + 1) ** 6,
                "fixed": (10, (1, 9, 0), -81, [*([1] * 6), *([3] * 4)]),
                "coin": (12, (2, 10, 0), 81, [*([1] * 8), *([3] * 4)]),
            }
        if not (
            sp.factor(k3_action.charpoly().as_expr()) == expected_cohom["charpoly"]
            and (
                fixed_cohom_basis.rows,
                exact_inertia(fixed_cohom_gram),
                fixed_cohom_gram.det(),
                smith_diagonal(fixed_cohom_gram),
            )
            == expected_cohom["fixed"]
            and (
                coin_cohom_basis.rows,
                exact_inertia(coin_cohom_gram),
                coin_cohom_gram.det(),
                smith_diagonal(coin_cohom_gram),
            )
            == expected_cohom["coin"]
        ):
            raise AssertionError(f"the full cohomology invariants failed for {case_id}")

        fixed_coordinates = sp.Matrix(
            [sp.Rational(value) for value in section["strict_fixed_interior"]["fixed_lattice_coordinates"]]
        )
        fixed_interior = fixed_coordinates.T * archived_fixed_basis
        source_pairings = list(fixed_interior * n_gram * wall_rows.T)
        mirror_records = section["restricted_mirrors"]["records"]
        active_records = [
            row for row in mirror_records if row["hyperbolic_status"] == "active_facet"
        ]
        active_source_norms = [
            int((wall_rows[index - 1, :] * n_gram * wall_rows[index - 1, :].T)[0])
            for row in active_records
            for index in row["source_wall_indices"]
        ]
        reduction = section["exact_reduction"]
        if not (
            fixed_interior * n_action == fixed_interior
            and (fixed_interior * n_gram * fixed_interior.T)[0] > 0
            and min(source_pairings) > 0
            and len(active_records) == 14
            and set(active_source_norms) == {-2}
            and factory_case["source_wall_audit"]["active_facet_source_norm_histogram"]
            == {"-2": 14, "-4": 0}
            and factory_case["source_wall_audit"]["folded_primitive_facet_norm_histogram"]
            == config["folded_norm_histogram"]
            and reduction["unresolved_orbit_indices"] == []
            and reduction["causally_inactive_orbit_indices"] == []
            and reduction["active_hyperbolic_facet_count"] == 14
            and reduction["conic_redundant_count"] == len(mirror_records) - 14
            and all(exact_cone["checks"].values())
            and exact_cone["decision"]["finite_volume"] is True
            and exact_cone["decision"]["compact"] is False
            and exact_cone["exact_double_description"]["finite_vertex_count"]
            == config["finite_vertices"]
            and exact_cone["exact_double_description"]["ideal_vertex_count"]
            == config["ideal_vertices"]
        ):
            raise AssertionError(f"the chamber/short-wall gate failed for {case_id}")

        records.append(
            {
                "id": case_id,
                "diagram_type": config["diagram_type"],
                "ambient_NS": {
                    "lattice": "(18,4,0)",
                    "rank_signature_determinant": [18, [1, 17], -16],
                    "action_order": order,
                    "action_sha256": matrix_digest(n_action),
                    "action_characteristic_polynomial": str(
                        sp.factor(n_action.charpoly().as_expr())
                    ),
                    "discriminant_basis_half_vectors": [list(row) for row in discr_basis],
                    "discriminant_action_mod_two": gf2_payload(discr_action_matrix),
                    "fixed_lattice": {
                        "rank_signature_determinant": [10, [1, 9], int(fixed_gram.det())],
                        "smith_diagonal_absolute": smith_diagonal(fixed_gram),
                        "gram_sha256": matrix_digest(fixed_gram),
                    },
                    "coinvariant_lattice": {
                        "rank_signature_determinant": [8, [0, 8], int(coin_gram.det())],
                        "smith_diagonal_absolute": smith_diagonal(coin_gram),
                        "basis_in_N": matrix_payload(coin_basis),
                        "gram": matrix_payload(coin_gram),
                        "action": matrix_payload(coin_action),
                        "fixed_plus_coinvariant_index_in_N": direct_sum_index,
                        "minus_two_rootless": True,
                        **coin_rootlessness,
                    },
                },
                "transcendental_and_glue": {
                    "T": "U(2)+U(2)",
                    "T_action": matrix_payload(t_action),
                    "anti_isometry_A_N_to_A_T_mod_two": matrix_payload(identification),
                    "quadratic_form_table": q_table,
                    "actions_equivariant": True,
                    "glued_even_unimodular_signature": [3, 19],
                    "glued_basis_sha256": matrix_digest(k3_basis),
                    "glued_gram_sha256": matrix_digest(k3_gram),
                    "glued_action_sha256": matrix_digest(k3_action),
                    "glued_action_characteristic_polynomial": str(
                        sp.factor(k3_action.charpoly().as_expr())
                    ),
                    "fixed_cohomology_rank_signature_determinant": [
                        fixed_cohom_basis.rows,
                        list(exact_inertia(fixed_cohom_gram)[:2]),
                        int(fixed_cohom_gram.det()),
                    ],
                    "coinvariant_cohomology_rank_signature_determinant": [
                        coin_cohom_basis.rows,
                        list(exact_inertia(coin_cohom_gram)[:2]),
                        int(coin_cohom_gram.det()),
                    ],
                },
                "period_and_torelli": {
                    "period_in_T_tensor_C": [
                        "sqrt(2)*zeta_3^2",
                        "-zeta_3",
                        "1",
                        "sqrt(2)",
                    ],
                    "period_eigenvalue": config["period_eigenvalue_text"],
                    "period_square": "0",
                    "period_positive_norm": "6*sqrt(2)",
                    "rational_T_orthogonal": "0",
                    "rational_orthogonality_equation_rank": 4,
                    "picard_lattice_is_exactly_N": True,
                    "source_48_wall_chamber_preserved": True,
                    "NS_coinvariant_has_no_minus_two_vectors": True,
                    "invariant_kahler_chamber_exists": True,
                    "global_torelli_conclusion": (
                        "there exists a projective complex K3 surface with NS=N and a genuine "
                        + str(config["automorphism_character"])
                        + " inducing the archived action"
                    ),
                },
                "nef_and_quotient_cone": {
                    "archived_dimension_facets": [9, 14],
                    "finite_volume_noncompact": True,
                    "finite_and_ideal_vertices": [
                        config["finite_vertices"],
                        config["ideal_vertices"],
                    ],
                    "active_source_wall_norm_histogram": {"-2": 14, "-4": 0},
                    "folded_primitive_facet_norm_histogram": config[
                        "folded_norm_histogram"
                    ],
                    "all_restricted_long_walls_redundant": True,
                    "invariant_nef_section_equality": (
                        "factory cone = P_r intersect Fix(g) = Nef(X) intersect Fix(g)"
                    ),
                    "equality_proof": (
                        "P_r is contained in one W_2 chamber C.  C satisfies every oriented "
                        "short P_r-wall inequality.  Exact reduction says the restriction of "
                        "all 48 P_r inequalities is generated by the 14 active short-wall "
                        "restrictions, so C intersect Fix(g) is contained in P_r intersect "
                        "Fix(g); the reverse containment follows from P_r contained in C"
                    ),
                    "quotient_nef_conclusion": (
                        "for the normal finite quotient Y=X/<g>, pullback identifies "
                        "N^1(Y)_R with N^1(X)_R^<g> and Nef(Y) with this invariant "
                        "nef section; this does not identify the nef cone of a resolution of Y"
                    ),
                },
            }
        )

    certificate: dict[str, Any] = {
        "schema": "a18_c01_c07_k3_nef_promotions_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "case_count": len(records),
        "cases": records,
        "theorem": (
            "The A18 factory cones C01=T2 and C07=T3 are respectively the full "
            "invariant nef sections, equivalently the nef cones of the normal finite "
            "quotients, for a symplectic Nikulin involution and a purely non-symplectic "
            "order-three automorphism on complex projective K3 surfaces with exact "
            "Picard lattice (18,4,0)."
        ),
        "novelty_scope": (
            "strong geometric promotions of already certified arithmetic fixed cones; "
            "not new abstract polytopes or new Coxeter-volume computations"
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
