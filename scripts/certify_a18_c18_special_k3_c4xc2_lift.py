#!/usr/bin/env python3
"""Exact C4 x C2 K3 lift of the A18 C18 V4 fixed section.

The three nonidentity Picard actions of C18 lift simultaneously to commuting
isometries of the K3 lattice.  Two lifts have order four, have the same
central square, and generate C4 x C2; their restriction to N=(18,4,0) is the
original V4.  A common special period gives Picard rank 20, while the common
fixed Neron--Severi lattice remains exactly the archived rank-11 lattice.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
    matrix_digest,
    matrix_payload,
    q_table_and_equivariance,
    same_row_lattice,
    saturated_row_kernel,
    sha256,
)
from scripts.certify_a18_c06_h13_31_k3_order4_lift import (
    discriminant_action,
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
C03_C04_HELPER = (
    ROOT / "scripts" / "certify_a18_c03_c04_special_k3_symplectic_lifts.py"
)
C06_HELPER = ROOT / "scripts" / "certify_a18_c06_h13_31_k3_order4_lift.py"
OUTPUT = (
    ROOT / "data" / "enumerations" / "a18_c18_special_k3_c4xc2_lift.json"
)

EXPECTED_SOURCE_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "lineage": "dea8024f9e6e58abd0483e339774e8023b1bd511a6ac208d7db6e5e709187173",
    "c03_c04_helper": "1c15dd20d9fbb2754b154570ef453423fe45961cd6fb263d4d703f40adc6b598",
    "c06_helper": "967d6e258a00df16ff44ad3f762f81817676152d882be2fb3088894e87562027",
}


def fraction_inverse(
    matrix: tuple[tuple[Fraction, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    size = len(matrix)
    work = [
        [*matrix[row], *(Fraction(row == column) for column in range(size))]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][index] - scale * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def fraction_quadratic(
    vector: tuple[int, ...], matrix: tuple[tuple[Fraction, ...], ...]
) -> Fraction:
    return sum(
        Fraction(vector[i]) * matrix[i][j] * Fraction(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def fraction_prefix_schur(
    matrix: tuple[tuple[Fraction, ...], ...]
) -> tuple[tuple[tuple[Fraction, ...], ...], ...]:
    size = len(matrix)
    output = []
    for prefix in range(1, size + 1):
        leading = tuple(
            tuple(matrix[i][j] for j in range(prefix)) for i in range(prefix)
        )
        if prefix == size:
            output.append(leading)
            continue
        cross = tuple(
            tuple(matrix[i][j] for j in range(prefix, size))
            for i in range(prefix)
        )
        tail = tuple(
            tuple(matrix[i][j] for j in range(prefix, size))
            for i in range(prefix, size)
        )
        tail_inverse = fraction_inverse(tail)
        output.append(
            tuple(
                tuple(
                    leading[i][j]
                    - sum(
                        cross[i][a] * tail_inverse[a][b] * cross[j][b]
                        for a in range(size - prefix)
                        for b in range(size - prefix)
                    )
                    for j in range(prefix)
                )
                for i in range(prefix)
            )
        )
    return tuple(output)


def floor_sqrt_fraction(value: Fraction) -> int:
    candidate = isqrt(value.numerator // value.denominator)
    while Fraction((candidate + 1) ** 2) <= value:
        candidate += 1
    while Fraction(candidate**2) > value:
        candidate -= 1
    return candidate


def integer_short_vector_audit(
    gram: sp.Matrix, target: int = 2
) -> dict[str, object]:
    rational = tuple(
        tuple(Fraction(int(gram[i, j])) for j in range(gram.cols))
        for i in range(gram.rows)
    )
    inverse = fraction_inverse(rational)
    bounds = tuple(
        floor_sqrt_fraction(Fraction(target) * inverse[index][index])
        for index in range(len(rational))
    )
    prefixes = fraction_prefix_schur(rational)
    visited = [0] * (len(rational) + 1)
    pruned = [0] * (len(rational) + 1)
    vectors: list[tuple[int, ...]] = []

    def extend(prefix: tuple[int, ...]) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth and fraction_quadratic(prefix, prefixes[depth - 1]) > target:
            pruned[depth] += 1
            return
        if depth == len(rational):
            if fraction_quadratic(prefix, prefixes[-1]) == target:
                vectors.append(prefix)
            return
        for entry in range(-bounds[depth], bounds[depth] + 1):
            extend((*prefix, entry))

    extend(())
    return {
        "target": target,
        "coordinate_bounds": list(bounds),
        "visited_by_depth": visited,
        "schur_pruned_by_depth": pruned,
        "vector_count": len(vectors),
        "vectors": [list(vector) for vector in vectors],
    }


def build_certificate() -> dict[str, Any]:
    actual_hashes = {
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
        "c03_c04_helper": sha256(C03_C04_HELPER),
        "c06_helper": sha256(C06_HELPER),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A18 source/helper changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        construction["all_checks_pass"]
        and vinberg["all_checks_pass"]
        and lineage["all_checks_pass"]
    ):
        raise AssertionError("an upstream C18 construction/volume/lineage audit failed")

    inventory = next(
        row
        for row in construction["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    case = next(row for row in construction["cases"] if row["id"] == "A18_4_0_C18")
    volume = next(row for row in vinberg["cases"] if row["id"] == "A18_4_0_C18")
    lineage_case = next(row for row in lineage["cases"] if row["id"] == "A18_4_0_C18")
    if not (
        all(case["checks"].values())
        and case["case"]["subgroup_order"] == 4
        and case["case"]["cyclic"] is False
        and case["case"]["fixed_rank"] == 11
        and case["exact_reduction"]["active_hyperbolic_facet_count"] == 21
        and volume["classification"] == "finite_volume_noncompact_Coxeter_polytope"
        and volume["finite_volume"] is True
        and volume["bad_edge_count"] == 0
        and volume["bad_parabolic_neighborhood_count"] == 0
        and lineage_case["id"] == "A18_4_0_C18"
    ):
        raise AssertionError("the upstream C18 target changed")

    n_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in inventory["ambient"]["gram"]]
    )
    wall_rows = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in inventory["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    n_actions = [
        sp.Matrix([[sp.Rational(value) for value in row] for row in payload])
        for payload in case["case"]["representative_nonidentity_primal_actions"]
    ]
    if not (
        len(n_actions) == 3
        and all(action**2 == sp.eye(18) for action in n_actions)
        and all(action * n_gram * action.T == n_gram for action in n_actions)
        and n_actions[1] * n_actions[2] == n_actions[0]
        and all(
            n_actions[left] * n_actions[right]
            == n_actions[right] * n_actions[left]
            for left in range(3)
            for right in range(left)
        )
    ):
        raise AssertionError("the archived Picard V4 action changed")

    t_gram = sp.Matrix(
        [
            [0, 2, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 0, 2],
            [0, 0, 2, 0],
        ]
    )
    b_left = sp.Matrix(
        [
            [0, 0, 0, -1],
            [0, 0, -1, 0],
            [0, 1, 0, 0],
            [1, 0, 0, 0],
        ]
    )
    b_right = sp.Matrix(
        [
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [-1, 0, 0, 0],
            [0, -1, 0, 0],
        ]
    )
    t_actions = [b_left * b_right, b_left, b_right]
    identification = sp.Matrix(
        [
            [0, 0, 1, 1],
            [0, 0, 0, 1],
            [1, 1, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    if not (
        all(action * t_gram * action.T == t_gram for action in t_actions)
        and t_actions[0] ** 2 == sp.eye(4)
        and t_actions[1] ** 2 == -sp.eye(4)
        and t_actions[2] ** 2 == -sp.eye(4)
        and t_actions[1] * t_actions[2] == t_actions[0]
        and t_actions[1] * t_actions[2] == t_actions[2] * t_actions[1]
    ):
        raise AssertionError("the commuting transcendental lift failed")

    discr_basis = two_elementary_discriminant_basis(n_gram)
    discr_actions = []
    q_table = None
    for n_action, t_action in zip(n_actions, t_actions, strict=True):
        discr_action_matrix, current_table = q_table_and_equivariance(
            n_action,
            n_gram,
            discr_basis,
            t_action,
            identification,
        )
        discr_actions.append(discr_action_matrix)
        if q_table is None:
            q_table = current_table
        elif current_table != q_table:
            raise AssertionError("the common discriminant anti-isometry changed")

    glue_rows = sp.Matrix(
        [
            [*discr_basis[index], *map(int, identification.row(index))]
            for index in range(4)
        ]
    )
    k3_basis = hermite_normal_form(
        (2 * sp.eye(22)).col_join(glue_rows).T
    ).T / 2
    k3_gram = sp.simplify(k3_basis * sp.diag(n_gram, t_gram) * k3_basis.T)
    k3_actions = [
        sp.simplify(k3_basis * sp.diag(n_action, t_action) * k3_basis.inv())
        for n_action, t_action in zip(n_actions, t_actions, strict=True)
    ]
    lambda_symbol = sp.Symbol("lambda")
    expected_characteristic_polynomials = [
        (lambda_symbol - 1) ** 14 * (lambda_symbol + 1) ** 8,
        (lambda_symbol - 1) ** 14
        * (lambda_symbol + 1) ** 4
        * (lambda_symbol**2 + 1) ** 2,
        (lambda_symbol - 1) ** 14
        * (lambda_symbol + 1) ** 4
        * (lambda_symbol**2 + 1) ** 2,
    ]
    if not (
        all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(all(value.q == 1 for value in action) for action in k3_actions)
        and all(
            action * k3_gram * action.T == k3_gram for action in k3_actions
        )
        and k3_actions[1] * k3_actions[2] == k3_actions[0]
        and k3_actions[1] * k3_actions[2] == k3_actions[2] * k3_actions[1]
        and k3_actions[1] ** 2 == k3_actions[2] ** 2
        and k3_actions[1] ** 2 != sp.eye(22)
        and all(
            sp.factor(action.charpoly().as_expr()) == expected
            for action, expected in zip(
                k3_actions, expected_characteristic_polynomials, strict=True
            )
        )
    ):
        raise AssertionError("simultaneous K3 glue failed")

    generated = {
        tuple(k3_actions[1] ** left_power * k3_actions[2] ** right_power)
        for left_power in range(4)
        for right_power in range(4)
    }
    if len(generated) != 8:
        raise AssertionError("the glued group is not order eight")

    period = sp.Matrix([[-sp.I, -sp.I, 1, 1]])
    period_eigenvalues = [sp.Integer(1), sp.I, -sp.I]
    if not (
        all(
            period * action == eigenvalue * period
            for action, eigenvalue in zip(
                t_actions, period_eigenvalues, strict=True
            )
        )
        and (period * t_gram * period.T)[0] == 0
        and (period * t_gram * sp.conjugate(period).T)[0] == 8
    ):
        raise AssertionError("common C4 x C2 period failed")

    rational_constraints = sp.Matrix([[1, 0], [1, 0], [0, 1], [0, 1]])
    ns_basis = saturated_row_kernel(
        k3_basis[:, 18:22] * rational_constraints
    )
    ns_gram = ns_basis * k3_gram * ns_basis.T
    ns_smith_matrix = smith_normal_form(ns_gram, domain=ZZ)
    ns_smith = [
        abs(int(ns_smith_matrix[index, index])) for index in range(ns_basis.rows)
    ]
    ns_actions = [induced_row_action(ns_basis, action) for action in k3_actions]
    if not (
        ns_basis.rows == 20
        and exact_inertia(ns_gram) == (1, 19, 0)
        and ns_gram.det() == -16
        and ns_smith == [*[1] * 18, 4, 4]
    ):
        raise AssertionError("the special C18 Picard lattice failed")

    fixed_ns_coordinates = saturated_row_kernel(
        sp.Matrix.hstack(
            *(action - sp.eye(ns_basis.rows) for action in ns_actions)
        )
    )
    fixed_ns_in_k3 = fixed_ns_coordinates * ns_basis
    n_embedding = sp.Matrix.hstack(sp.eye(18), sp.zeros(18, 4)) * k3_basis.inv()
    archived_fixed_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["saturated_fixed_lattice"]["basis_primal_rows"]
        ]
    )
    archived_fixed_in_k3 = archived_fixed_basis * n_embedding
    fixed_ns_gram = fixed_ns_coordinates * ns_gram * fixed_ns_coordinates.T
    if not (
        fixed_ns_coordinates.rows == 11
        and exact_inertia(fixed_ns_gram) == (1, 10, 0)
        and fixed_ns_gram.det() == 128
        and same_row_lattice(fixed_ns_in_k3, archived_fixed_in_k3)
    ):
        raise AssertionError("the common fixed NS lattice enlarged")

    coinvariant_coordinates = saturated_row_kernel(
        ns_gram * fixed_ns_coordinates.T
    )
    coinvariant_gram = (
        coinvariant_coordinates * ns_gram * coinvariant_coordinates.T
    )
    coinvariant_smith_matrix = smith_normal_form(coinvariant_gram, domain=ZZ)
    coinvariant_smith = [
        abs(int(coinvariant_smith_matrix[index, index]))
        for index in range(coinvariant_coordinates.rows)
    ]
    root_audit = integer_short_vector_audit(-coinvariant_gram, target=2)
    if not (
        coinvariant_coordinates.rows == 9
        and exact_inertia(coinvariant_gram) == (0, 9, 0)
        and coinvariant_gram.det() == -512
        and coinvariant_smith == [1, 1, 2, 2, 2, 2, 2, 2, 8]
        and root_audit["vector_count"] == 0
    ):
        raise AssertionError("the common algebraic coinvariant is not rootless")

    fixed_cohom_coordinates = saturated_row_kernel(
        sp.Matrix.hstack(
            *(action - sp.eye(22) for action in k3_actions)
        )
    )
    fixed_cohom_gram = (
        fixed_cohom_coordinates * k3_gram * fixed_cohom_coordinates.T
    )
    if not (
        fixed_cohom_coordinates.rows == 11
        and exact_inertia(fixed_cohom_gram) == (1, 10, 0)
        and fixed_cohom_gram.det() == 128
        and same_row_lattice(fixed_cohom_coordinates, archived_fixed_in_k3)
    ):
        raise AssertionError("common fixed cohomology is not the fixed Picard lattice")

    fixed_coordinates = sp.Matrix(
        [
            sp.Rational(value)
            for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]
        ]
    )
    fixed_interior = fixed_coordinates.T * archived_fixed_basis
    source_pairings = list(fixed_interior * n_gram * wall_rows.T)
    wall_lookup = {tuple(wall_rows.row(index)) for index in range(48)}
    if not (
        all(fixed_interior * action == fixed_interior for action in n_actions)
        and (fixed_interior * n_gram * fixed_interior.T)[0] > 0
        and all(value > 0 for value in source_pairings)
        and all(
            {
                tuple(wall_rows.row(index) * action) for index in range(48)
            }
            == wall_lookup
            for action in n_actions
        )
    ):
        raise AssertionError("the C18 source-chamber interior is not preserved")

    certificate: dict[str, Any] = {
        "schema": "a18_c18_special_k3_c4xc2_lift_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "ambient_picard_action": {
            "ambient_invariants": [18, 4, 0],
            "subgroup_on_N": "V4",
            "common_fixed_rank": 11,
            "nonidentity_action_sha256": [
                matrix_digest(action) for action in n_actions
            ],
            "nonidentity_discriminant_actions": [
                [
                    [int(value) % 2 for value in action.row(index)]
                    for index in range(4)
                ]
                for action in discr_actions
            ],
        },
        "simultaneous_transcendental_glue": {
            "lattice": "U(2)+U(2)",
            "actions_in_V4_nonidentity_order": [
                matrix_payload(action) for action in t_actions
            ],
            "action_orders": [2, 4, 4],
            "two_order_four_actions_commute": True,
            "two_order_four_actions_have_common_square_minus_identity": True,
            "anti_isometry_source_to_T_mod_two": matrix_payload(identification),
            "all_three_actions_equivariant": True,
            "quadratic_form_table": q_table,
        },
        "k3_group_action": {
            "abstract_group": "C4 x C2",
            "order": 8,
            "restriction_to_N": "V4",
            "central_kernel_on_N": "the common square acts as +1 on N and -1 on T",
            "nonidentity_lift_characteristic_polynomials": [
                "(x-1)^14*(x+1)^8",
                "(x-1)^14*(x+1)^4*(x^2+1)^2",
                "(x-1)^14*(x+1)^4*(x^2+1)^2",
            ],
            "glued_basis_sha256": matrix_digest(k3_basis),
            "glued_gram_sha256": matrix_digest(k3_gram),
            "three_lift_action_sha256": [
                matrix_digest(action) for action in k3_actions
            ],
        },
        "period_and_picard_lattice": {
            "common_period_in_T_tensor_C": ["-i", "-i", "1", "1"],
            "three_eigenvalues": ["1", "i", "-i"],
            "period_square": "0",
            "period_positive_norm": "8",
            "rational_T_orthogonal": (
                "for (a,b,c,d) in Q^4, pairing gives "
                "2*(-i*(a+b)+(c+d)); hence a+b=c+d=0"
            ),
            "picard_rank": 20,
            "picard_signature": [1, 19],
            "picard_determinant": -16,
            "picard_smith_diagonal_absolute": ns_smith,
            "picard_gram_sha256": matrix_digest(ns_gram),
            "common_fixed_NS_is_exactly_archived_C18_fixed_lattice": True,
            "common_fixed_NS_rank": 11,
            "common_fixed_NS_determinant": 128,
        },
        "rootless_coinvariant_and_torelli": {
            "common_algebraic_coinvariant_rank": 9,
            "common_algebraic_coinvariant_signature": [0, 9],
            "common_algebraic_coinvariant_determinant": -512,
            "common_algebraic_coinvariant_smith_diagonal_absolute": coinvariant_smith,
            "minus_two_vector_audit_on_negative_coinvariant": root_audit,
            "rootless": True,
            "source_48_wall_chamber_is_preserved": True,
            "archived_strict_fixed_interior_minimum_source_pairing": str(
                min(source_pairings)
            ),
            "kahler_reason": (
                "no (-2)-class is orthogonal to the entire common fixed positive "
                "cone, because the common algebraic coinvariant is rootless; a "
                "generic invariant point in the strict source-section interior "
                "therefore selects a group-invariant Kahler chamber"
            ),
            "global_torelli_conclusion": (
                "a special complex projective K3 surface carries a genuine "
                "C4 x C2 automorphism group whose restriction to N is C18=V4"
            ),
        },
        "hyperbolic_fixed_section": {
            "id": "A18_4_0_C18",
            "dimension": 10,
            "facet_count": 21,
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
        "theorem": (
            "A18_4_0_C18 (H10/21) is the common fixed divisor/reflection-"
            "chamber section of a genuine C4 x C2 action on a special "
            "Picard-rank-20 K3 surface; the action restricts to C18=V4 on N"
        ),
        "novelty_scope": (
            "new project-level algebraic lift of an already certified chamber; "
            "not a new abstract polytope, arithmetic class, dimension, or finite "
            "K3-group type"
        ),
        "nef_scope_caution": (
            "the source 48-wall chamber is a complete mixed (-2)/(-4)-root "
            "reflection chamber, not asserted to be the literal nef cone"
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
