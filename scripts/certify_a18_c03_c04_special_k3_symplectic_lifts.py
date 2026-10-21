#!/usr/bin/env python3
"""Exact special-period K3 lifts of the A18 C03 and C04 sections.

The two selected involutions do not use the generic Picard lattice alone.
Instead, an explicit discriminant-compatible involution of
``T=U(2)+U(2)`` fixes the period and makes its anti-invariant part algebraic.
The resulting K3 surfaces have Picard ranks 20 and 19.  In both cases the
full anti-invariant cohomology lattice is ``E8(-2)``, while the fixed
Neron--Severi lattice is exactly the archived fixed sublattice of
``N=(18,4,0)``.  Thus the H11/28 and H10/28 fixed sections acquire genuine
symplectic Nikulin-involution realizations without changing their dimensions.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.certify_a18_c06_h13_31_k3_order4_lift import (
    discriminant_action,
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
C06_HELPER = ROOT / "scripts" / "certify_a18_c06_h13_31_k3_order4_lift.py"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "a18_c03_c04_special_k3_symplectic_lifts.json"
)

EXPECTED_SOURCE_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "lineage": "dea8024f9e6e58abd0483e339774e8023b1bd511a6ac208d7db6e5e709187173",
    "c06_helper": "967d6e258a00df16ff44ad3f762f81817676152d882be2fb3088894e87562027",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _domain_matrix(matrix: sp.Matrix) -> DomainMatrix:
    return DomainMatrix.from_Matrix(matrix).convert_to(ZZ)


def saturated_row_kernel(equations: sp.Matrix) -> sp.Matrix:
    """Return a saturated Z-basis for row vectors x with x*equations=0."""

    transposed = equations.T
    denominator = int(
        sp.ilcm(*(sp.Rational(value).q for value in transposed))
    )
    integral = sp.Matrix(denominator * transposed)
    diagonal, _, right = smith_normal_decomp(_domain_matrix(integral))
    diagonal_matrix = diagonal.to_Matrix()
    right_matrix = right.to_Matrix()
    zero_columns = [
        column
        for column in range(integral.cols)
        if all(
            diagonal_matrix[row, column] == 0
            for row in range(diagonal_matrix.rows)
        )
    ]
    basis = right_matrix[:, zero_columns].T
    if basis * equations != sp.zeros(basis.rows, equations.cols):
        raise AssertionError("Smith kernel replay failed")
    return basis


def induced_row_action(basis: sp.Matrix, ambient_action: sp.Matrix) -> sp.Matrix:
    pivots = basis.rref()[1]
    square = basis[:, list(pivots)]
    action = (basis * ambient_action)[:, list(pivots)] * square.inv()
    if not (
        all(value.q == 1 for value in action)
        and action * basis == basis * ambient_action
    ):
        raise AssertionError("induced integral row action failed")
    return action


def same_row_lattice(left: sp.Matrix, right: sp.Matrix) -> bool:
    if left.shape != right.shape or left.rank() != left.rows:
        return False
    pivots = right.rref()[1]
    change = left[:, list(pivots)] * right[:, list(pivots)].inv()
    return (
        all(value.q == 1 for value in change)
        and abs(int(change.det())) == 1
        and change * right == left
    )


def gf2_matrix(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(value) % 2 for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def q_table_and_equivariance(
    n_action: sp.Matrix,
    n_gram: sp.Matrix,
    discr_basis: tuple[tuple[int, ...], ...],
    t_action: sp.Matrix,
    identification: sp.Matrix,
) -> tuple[sp.Matrix, list[dict[str, object]]]:
    n_discr_action = discriminant_action(n_action, discr_basis)
    t_mod_two = t_action.applyfunc(lambda value: int(value) % 2)
    if (n_discr_action * identification).applyfunc(
        lambda value: int(value) % 2
    ) != (identification * t_mod_two).applyfunc(lambda value: int(value) % 2):
        raise AssertionError("discriminant actions are not equivariant")

    table = []
    for coordinates in itertools.product(range(2), repeat=4):
        target = tuple(
            int(value) % 2
            for value in (sp.Matrix([coordinates]) * identification).row(0)
        )
        q_n = q_n_mod_two(coordinates, discr_basis, n_gram)
        q_t = q_t_mod_two(target)
        if q_n != q_t:
            raise AssertionError("discriminant graph is not an anti-isometry")
        table.append(
            {"source": list(coordinates), "target": list(target), "q_mod_2": q_n}
        )
    return n_discr_action, table


def glue_k3_lattice(
    n_gram: sp.Matrix,
    n_action: sp.Matrix,
    discr_basis: tuple[tuple[int, ...], ...],
    t_gram: sp.Matrix,
    t_action: sp.Matrix,
    identification: sp.Matrix,
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
        all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(value.q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and k3_action**2 == sp.eye(22)
        and k3_action != sp.eye(22)
        and sp.factor(k3_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 14
        * (sp.Symbol("lambda") + 1) ** 8
    ):
        raise AssertionError("glued K3 involution failed")
    return k3_basis, k3_gram, k3_action


def build_certificate() -> dict[str, Any]:
    actual_hashes = {
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
        "c06_helper": sha256(C06_HELPER),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A18 source or exact helper changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        construction["all_checks_pass"]
        and vinberg["all_checks_pass"]
        and lineage["all_checks_pass"]
    ):
        raise AssertionError("an upstream A18 construction/volume/lineage audit failed")

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
    ):
        raise AssertionError("ambient (18,4,0) Picard stock changed")
    discr_basis = two_elementary_discriminant_basis(n_gram)
    t_gram = sp.Matrix(
        [
            [0, 2, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 0, 2],
            [0, 0, 2, 0],
        ]
    )

    fourth_root_two = sp.root(2, 4)
    configs = {
        "A18_4_0_C03": {
            "t_action": sp.Matrix(
                [
                    [0, 1, 0, 0],
                    [1, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 1, 0],
                ]
            ),
            "identification": sp.Matrix(
                [
                    [0, 0, 1, 1],
                    [0, 0, 0, 1],
                    [1, 1, 0, 0],
                    [0, 1, 0, 0],
                ]
            ),
            "period": sp.Matrix([[1, 1, sp.I, sp.I]]),
            "period_text": ["1", "1", "i", "i"],
            "period_norm": sp.Integer(8),
            "period_norm_text": "8",
            "rational_orthogonal_constraints": sp.Matrix(
                [[1, 0], [1, 0], [0, 1], [0, 1]]
            ),
            "orthogonality_reason": (
                "for (a,b,c,d) in Q^4, pairing with omega is "
                "2*(a+b)+2*i*(c+d), hence a+b=c+d=0"
            ),
            "picard_rank": 20,
            "picard_determinant": -16,
            "picard_smith": [*[1] * 18, 4, 4],
            "fixed_picard_rank": 12,
            "dimension": 11,
            "facets": 28,
        },
        "A18_4_0_C04": {
            "t_action": sp.Matrix(
                [
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 1, 0],
                ]
            ),
            "identification": sp.Matrix(
                [
                    [0, 1, 1, 1],
                    [1, 0, 0, 0],
                    [0, 0, 1, 1],
                    [1, 0, 0, 1],
                ]
            ),
            "period": sp.Matrix(
                [[1, sp.sqrt(2), sp.I * fourth_root_two, sp.I * fourth_root_two]]
            ),
            "period_text": ["1", "sqrt(2)", "i*2^(1/4)", "i*2^(1/4)"],
            "period_norm": 8 * sp.sqrt(2),
            "period_norm_text": "8*sqrt(2)",
            "rational_orthogonal_constraints": sp.Matrix(
                [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]]
            ),
            "orthogonality_reason": (
                "for (a,b,c,d) in Q^4, pairing gives "
                "2*(a*sqrt(2)+b)+2*i*2^(1/4)*(c+d); "
                "thus a=b=0 and c+d=0"
            ),
            "picard_rank": 19,
            "picard_determinant": 16,
            "picard_smith": [*[1] * 16, 2, 2, 4],
            "fixed_picard_rank": 11,
            "dimension": 10,
            "facets": 28,
        },
    }

    records = []
    for case_id, config in configs.items():
        construction_case = next(
            row for row in construction["cases"] if row["id"] == case_id
        )
        volume_case = next(row for row in vinberg["cases"] if row["id"] == case_id)
        lineage_case = next(row for row in lineage["cases"] if row["id"] == case_id)
        n_action = sp.Matrix(
            [
                [sp.Rational(value) for value in row]
                for row in construction_case["case"][
                    "representative_nonidentity_primal_actions"
                ][0]
            ]
        )
        t_action: sp.Matrix = config["t_action"]
        identification: sp.Matrix = config["identification"]
        period: sp.Matrix = config["period"]
        rational_constraints: sp.Matrix = config[
            "rational_orthogonal_constraints"
        ]

        if not (
            all(construction_case["checks"].values())
            and construction_case["case"]["subgroup_order"] == 2
            and construction_case["case"]["fixed_rank"]
            == config["fixed_picard_rank"]
            and construction_case["exact_reduction"][
                "active_hyperbolic_facet_count"
            ]
            == config["facets"]
            and volume_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and volume_case["finite_volume"] is True
            and volume_case["bad_edge_count"] == 0
            and volume_case["bad_parabolic_neighborhood_count"] == 0
            and lineage_case["id"] == case_id
        ):
            raise AssertionError(f"upstream target changed for {case_id}")
        if not (
            n_action**2 == sp.eye(18)
            and n_action * n_gram * n_action.T == n_gram
            and len((n_action - sp.eye(18)).nullspace())
            == config["fixed_picard_rank"]
            and t_action**2 == sp.eye(4)
            and t_action * t_gram * t_action.T == t_gram
            and len((t_action - sp.eye(4)).nullspace())
            == 14 - config["fixed_picard_rank"]
        ):
            raise AssertionError(f"Picard/transcendental involution failed for {case_id}")

        discr_action_matrix, q_table = q_table_and_equivariance(
            n_action,
            n_gram,
            discr_basis,
            t_action,
            identification,
        )
        k3_basis, k3_gram, k3_action = glue_k3_lattice(
            n_gram,
            n_action,
            discr_basis,
            t_gram,
            t_action,
            identification,
        )

        if not (
            period * t_action == period
            and sp.simplify((period * t_gram * period.T)[0]) == 0
            and sp.simplify((period * t_gram * sp.conjugate(period).T)[0])
            == config["period_norm"]
        ):
            raise AssertionError(f"symplectic period failed for {case_id}")

        # The displayed rational constraints are exactly the rational part of
        # omega-perpendicular in T.  They are then pulled through the glued
        # basis to recover NS(X) as a saturated Smith kernel.
        t_component = k3_basis[:, 18:22]
        ns_basis = saturated_row_kernel(t_component * rational_constraints)
        ns_gram = ns_basis * k3_gram * ns_basis.T
        ns_action = induced_row_action(ns_basis, k3_action)
        ns_smith_matrix = smith_normal_form(ns_gram, domain=ZZ)
        ns_smith = [
            abs(int(ns_smith_matrix[index, index]))
            for index in range(ns_basis.rows)
        ]
        if not (
            ns_basis.rows == config["picard_rank"]
            and exact_inertia(ns_gram)
            == (1, config["picard_rank"] - 1, 0)
            and ns_gram.det() == config["picard_determinant"]
            and ns_smith == config["picard_smith"]
            and ns_action**2 == sp.eye(config["picard_rank"])
        ):
            raise AssertionError(f"special-period Picard lattice failed for {case_id}")

        n_embedding = sp.Matrix.hstack(sp.eye(18), sp.zeros(18, 4)) * k3_basis.inv()
        if not all(value.q == 1 for value in n_embedding):
            raise AssertionError("N is not integral in the glued K3 lattice")
        ns_pivots = ns_basis.rref()[1]
        n_in_ns = (
            n_embedding[:, list(ns_pivots)]
            * ns_basis[:, list(ns_pivots)].inv()
        )
        if not (
            all(value.q == 1 for value in n_in_ns)
            and n_in_ns * ns_basis == n_embedding
        ):
            raise AssertionError("the ambient Picard stock is not contained in NS(X)")

        fixed_ns_coordinates = saturated_row_kernel(ns_action - sp.eye(ns_basis.rows))
        fixed_ns_in_k3 = fixed_ns_coordinates * ns_basis
        archived_fixed_basis = sp.Matrix(
            [
                [sp.Rational(value) for value in row]
                for row in construction_case["saturated_fixed_lattice"][
                    "basis_primal_rows"
                ]
            ]
        )
        archived_fixed_in_k3 = archived_fixed_basis * n_embedding
        if not (
            fixed_ns_coordinates.rows == config["fixed_picard_rank"]
            and same_row_lattice(fixed_ns_in_k3, archived_fixed_in_k3)
        ):
            raise AssertionError(f"fixed NS lattice enlarged for {case_id}")

        anti_basis = saturated_row_kernel(k3_action + sp.eye(22))
        anti_gram = anti_basis * k3_gram * anti_basis.T
        half_anti_gram = anti_gram / 2
        anti_smith_matrix = smith_normal_form(anti_gram, domain=ZZ)
        if not (
            anti_basis.rows == 8
            and exact_inertia(anti_gram) == (0, 8, 0)
            and anti_gram.det() == 256
            and [
                abs(int(anti_smith_matrix[index, index])) for index in range(8)
            ]
            == [2] * 8
            and all(value.q == 1 for value in half_anti_gram)
            and all(half_anti_gram[index, index] % 2 == 0 for index in range(8))
            and half_anti_gram.det() == 1
            and exact_inertia(half_anti_gram) == (0, 8, 0)
        ):
            raise AssertionError(f"anti-invariant lattice is not E8(-2) for {case_id}")

        fixed_cohom_basis = saturated_row_kernel(k3_action - sp.eye(22))
        fixed_cohom_gram = fixed_cohom_basis * k3_gram * fixed_cohom_basis.T
        if not (
            fixed_cohom_basis.rows == 14
            and exact_inertia(fixed_cohom_gram) == (3, 11, 0)
            and fixed_cohom_gram.det() == -256
        ):
            raise AssertionError(f"Nikulin invariant lattice failed for {case_id}")

        fixed_coordinates = sp.Matrix(
            [
                sp.Rational(value)
                for value in construction_case["strict_fixed_interior"][
                    "fixed_lattice_coordinates"
                ]
            ]
        )
        fixed_interior = fixed_coordinates.T * archived_fixed_basis
        source_pairings = list(fixed_interior * n_gram * wall_rows.T)
        wall_lookup = {tuple(wall_rows.row(index)) for index in range(48)}
        if not (
            fixed_interior * n_action == fixed_interior
            and (fixed_interior * n_gram * fixed_interior.T)[0] > 0
            and all(value > 0 for value in source_pairings)
            and {
                tuple(wall_rows.row(index) * n_action) for index in range(48)
            }
            == wall_lookup
        ):
            raise AssertionError(f"invariant source-chamber interior failed for {case_id}")

        records.append(
            {
                "id": case_id,
                "ambient_picard_action": {
                    "ambient_invariants": [18, 4, 0],
                    "order": 2,
                    "fixed_rank": config["fixed_picard_rank"],
                    "action_sha256": matrix_digest(n_action),
                    "discriminant_action": gf2_matrix(discr_action_matrix),
                },
                "transcendental_glue": {
                    "lattice": "U(2)+U(2)",
                    "involution": matrix_payload(t_action),
                    "fixed_rank": 14 - config["fixed_picard_rank"],
                    "anti_isometry_source_to_T_mod_two": matrix_payload(
                        identification
                    ),
                    "actions_equivariant": True,
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
                    "period_in_T_tensor_C": config["period_text"],
                    "period_eigenvalue": "1",
                    "period_square": "0",
                    "period_positive_norm": config["period_norm_text"],
                    "rational_T_orthogonal": config["orthogonality_reason"],
                    "picard_rank": config["picard_rank"],
                    "picard_signature": [1, config["picard_rank"] - 1],
                    "picard_determinant": config["picard_determinant"],
                    "picard_smith_diagonal_absolute": ns_smith,
                    "picard_gram_sha256": matrix_digest(ns_gram),
                    "ambient_N_is_contained": True,
                    "fixed_NS_is_exactly_archived_N_fixed_lattice": True,
                    "fixed_NS_rank": config["fixed_picard_rank"],
                },
                "torelli_and_chamber": {
                    "source_48_wall_chamber_is_preserved": True,
                    "archived_strict_fixed_interior_minimum_source_pairing": str(
                        min(source_pairings)
                    ),
                    "kahler_existence_reason": (
                        "the anti-invariant lattice is E8(-2), hence contains no "
                        "(-2)-class; therefore no root hyperplane contains the "
                        "fixed positive cone, and a generic invariant point in the "
                        "strict source-section interior lies in an invariant Kahler "
                        "chamber"
                    ),
                    "global_torelli_conclusion": (
                        "a special complex projective K3 surface carries a genuine "
                        "symplectic Nikulin involution inducing the selected action"
                    ),
                },
                "hyperbolic_fixed_section": {
                    "dimension": config["dimension"],
                    "facet_count": config["facets"],
                    "finite_volume": True,
                    "compact": False,
                    "finite_vertex_count": volume_case["finite_vertex_count"],
                    "ideal_vertex_count": volume_case["ideal_vertex_count"],
                    "elliptic_edge_count": volume_case["elliptic_edge_count"],
                    "bad_edge_count": volume_case["bad_edge_count"],
                    "bad_parabolic_neighborhood_count": volume_case[
                        "bad_parabolic_neighborhood_count"
                    ],
                    "criterion": volume_case["criterion"],
                },
            }
        )

    certificate: dict[str, Any] = {
        "schema": "a18_c03_c04_special_k3_symplectic_lifts_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "case_count": len(records),
        "cases": records,
        "theorem": (
            "A18_4_0_C03 (H11/28) and A18_4_0_C04 (H10/28) are fixed "
            "divisor/reflection-chamber sections of genuine symplectic Nikulin "
            "involutions on special Picard-rank 20 and 19 K3 surfaces"
        ),
        "novelty_scope": (
            "new project-level algebraic lifts of already certified chambers; "
            "not new abstract polytopes, dimensions, arithmetic lattices, or "
            "automorphism types"
        ),
        "nef_scope_caution": (
            "the source 48-wall chamber is a complete mixed (-2)/(-4)-root "
            "reflection chamber, not asserted to be the literal nef cone of the "
            "special K3 surfaces"
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

