#!/usr/bin/env python3
"""Exact order-four K3 lift of the A18_4_0_C06 H13/31 chamber.

The existing construction and Vinberg certificates start with the complete
48-wall chamber of the even 2-elementary Picard lattice (18,4,0), but use its
weighted-diagram involution only as an abstract lattice isometry.  This replay
proves that the selected involution extends to an order-four Hodge isometry of
the K3 lattice and hence, after choosing a generic compatible period, to a
genuine non-symplectic K3 automorphism.

The key point is explicit discriminant gluing.  The Picard discriminant form
is identified with that of T=U(2)+U(2).  The order-two Picard action matches
the reduction of an order-four isometry B of T with B^2=-1.  Gluing produces
the even unimodular K3 lattice, and an explicit positive i-eigenperiod has no
rational orthogonal class in T.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

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
    ROOT / "data" / "enumerations" / "a18_c06_h13_31_global_lineage.json"
)
OUTPUT = (
    ROOT / "data" / "enumerations" / "a18_c06_h13_31_k3_order4_lift.json"
)

EXPECTED_SOURCE_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "lineage": "53fb9ffb1cd151ca591a1922958ed25b3b6248daf0848e5ee653c32ec5e2de1b",
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


def two_elementary_discriminant_basis(gram: sp.Matrix) -> tuple[tuple[int, ...], ...]:
    """Basis bits a such that a/2 lies in the dual lattice."""

    size = gram.rows
    work = [[int(gram[i, j]) & 1 for j in range(size)] for i in range(size)]
    pivot_columns: list[int] = []
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
        pivot_columns.append(column)
        rank += 1
    free_columns = [column for column in range(size) if column not in pivot_columns]
    basis = []
    for free in free_columns:
        vector = [0] * size
        vector[free] = 1
        for row, pivot in reversed(list(enumerate(pivot_columns))):
            vector[pivot] = sum(
                work[row][column] * vector[column] for column in free_columns
            ) & 1
        basis.append(tuple(vector))
    if len(basis) != 4:
        raise AssertionError("the (18,4,0) discriminant length is not four")
    return tuple(basis)


def gf2_coordinates(
    vector: Sequence[int], basis: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    for coefficients in itertools.product(range(2), repeat=len(basis)):
        candidate = tuple(
            sum(coefficients[index] * basis[index][column] for index in range(len(basis)))
            & 1
            for column in range(len(vector))
        )
        if candidate == tuple(int(value) & 1 for value in vector):
            return tuple(coefficients)
    raise AssertionError("vector is not in the displayed GF(2) span")


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


def q_n_mod_two(coordinates: Sequence[int], basis: Sequence[Sequence[int]], gram: sp.Matrix) -> int:
    bits = tuple(
        sum(coordinates[index] * basis[index][column] for index in range(len(basis)))
        & 1
        for column in range(gram.rows)
    )
    numerator = sum(
        bits[i] * int(gram[i, j]) * bits[j]
        for i in range(gram.rows)
        for j in range(gram.cols)
    )
    if numerator % 4:
        raise AssertionError("dual half-vector has nonintegral discriminant value")
    return (numerator // 4) % 2


def q_t_mod_two(coordinates: Sequence[int]) -> int:
    return (coordinates[0] * coordinates[1] + coordinates[2] * coordinates[3]) % 2


def build_certificate() -> dict[str, object]:
    actual_hashes = {
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A18/C06 source certificate changed")

    construction = json.loads(CONSTRUCTION.read_text())
    vinberg = json.loads(VINBERG.read_text())
    lineage = json.loads(LINEAGE.read_text())
    if not (construction["all_checks_pass"] and vinberg["all_checks_pass"] and lineage["all_checks_pass"]):
        raise AssertionError("an upstream A18/C06 certificate is not green")

    inventory = next(
        example
        for example in construction["inventory"]["examples"]
        if example["key"] == [18, 4, 0]
    )
    case = next(row for row in construction["cases"] if row["id"] == "A18_4_0_C06")
    volume = next(row for row in vinberg["cases"] if row["id"] == "A18_4_0_C06")
    if not (
        all(case["checks"].values())
        and case["case"]["subgroup_order"] == 2
        and case["case"]["fixed_rank"] == 14
        and case["exact_reduction"]["active_hyperbolic_facet_count"] == 31
        and volume["classification"] == "finite_volume_noncompact_Coxeter_polytope"
        and volume["finite_volume"] is True
        and volume["bad_edge_count"] == 0
        and volume["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the upstream H13/31 construction or volume audit changed")

    n_gram = sp.Matrix([[sp.Rational(value) for value in row] for row in inventory["ambient"]["gram"]])
    n_action = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["case"]["representative_nonidentity_primal_actions"][0]
        ]
    )
    if not (
        n_gram.det() == -16
        and exact_inertia(n_gram) == (1, 17, 0)
        and n_action**2 == sp.eye(18)
        and n_action * n_gram * n_action.T == n_gram
        and len((n_action - sp.eye(18)).nullspace()) == 14
    ):
        raise AssertionError("ambient Picard involution failed")

    discr_basis = two_elementary_discriminant_basis(n_gram)
    if gf2_rank(discr_basis) != 4:
        raise AssertionError("Picard discriminant generators are not independent")
    discr_action = discriminant_action(n_action, discr_basis)

    # T=U(2)+U(2), in row-vector convention.
    t_gram = sp.Matrix(
        [
            [0, 2, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 0, 2],
            [0, 0, 2, 0],
        ]
    )
    t_action = sp.Matrix(
        [
            [0, 0, -1, 0],
            [0, 0, 0, -1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    discriminant_identification = sp.Matrix(
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
        and t_action**2 == -sp.eye(4)
        and t_action**4 == sp.eye(4)
        and t_action * t_gram * t_action.T == t_gram
        and abs(int(discriminant_identification.det())) % 2 == 1
    ):
        raise AssertionError("transcendental order-four witness failed")

    t_action_mod_two = t_action.applyfunc(lambda value: int(value) % 2)
    if (discr_action * discriminant_identification).applyfunc(
        lambda value: int(value) % 2
    ) != (discriminant_identification * t_action_mod_two).applyfunc(
        lambda value: int(value) % 2
    ):
        raise AssertionError("discriminant actions are not equivariant")

    q_table = []
    for coordinates in itertools.product(range(2), repeat=4):
        target = tuple(
            int(value) % 2
            for value in (sp.Matrix([coordinates]) * discriminant_identification).row(0)
        )
        q_n = q_n_mod_two(coordinates, discr_basis, n_gram)
        q_t = q_t_mod_two(target)
        if q_n != q_t:
            raise AssertionError("discriminant identification is not an anti-isometry")
        q_table.append({"source": list(coordinates), "target": list(target), "q_mod_2": q_n})

    # Graph-glue all four discriminant generators.  Multiplying by two turns
    # the overlattice basis problem into an ordinary integral HNF problem.
    doubled_generators = 2 * sp.eye(22)
    glue_rows = sp.Matrix(
        [
            [*discr_basis[index], *map(int, discriminant_identification.row(index))]
            for index in range(4)
        ]
    )
    doubled_module_rows = doubled_generators.col_join(glue_rows)
    k3_basis = hermite_normal_form(doubled_module_rows.T).T / 2
    orthogonal_sum_gram = sp.diag(n_gram, t_gram)
    k3_gram = sp.simplify(k3_basis * orthogonal_sum_gram * k3_basis.T)
    if not (
        all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
    ):
        raise AssertionError("discriminant glue is not the even unimodular K3 lattice")

    orthogonal_action = sp.diag(n_action, t_action)
    k3_action = sp.simplify(k3_basis * orthogonal_action * k3_basis.inv())
    if not (
        all(value.q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and k3_action**4 == sp.eye(22)
        and k3_action**2 != sp.eye(22)
        and sp.factor(k3_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 14
        * (sp.Symbol("lambda") + 1) ** 4
        * (sp.Symbol("lambda") ** 2 + 1) ** 2
    ):
        raise AssertionError("glued K3 action failed")

    sqrt_two = sp.sqrt(2)
    period = sp.Matrix([[1, sqrt_two, sp.I, sp.I * sqrt_two]])
    if not (
        period * t_action == sp.I * period
        and (period * t_gram * period.T)[0] == 0
        and sp.simplify((period * t_gram * sp.conjugate(period).T)[0]) == 8 * sqrt_two
    ):
        raise AssertionError("explicit positive i-eigenperiod failed")

    fixed_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["saturated_fixed_lattice"]["basis_primal_rows"]
        ]
    )
    interior_coordinates = sp.Matrix(
        [sp.Rational(value) for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    interior = interior_coordinates.T * fixed_basis
    wall_rows = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in inventory["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    wall_lookup = {tuple(wall_rows.row(index)) for index in range(wall_rows.rows)}
    chamber_wall_pairings = list(interior * n_gram * wall_rows.T)
    polarization_denominator = 1
    for value in interior:
        polarization_denominator = int(
            sp.ilcm(polarization_denominator, sp.Rational(value).q)
        )
    polarization = polarization_denominator * interior
    polarization_coordinates = [int(value) for value in polarization]
    if not (
        interior * n_action == interior
        and (interior * n_gram * interior.T)[0] == 385
        and case["strict_fixed_interior"]["minimum_folded_pairing"] == "1"
        and wall_rows.rows == 48
        and all(value > 0 for value in chamber_wall_pairings)
        and {
            tuple(wall_rows.row(index) * n_action)
            for index in range(wall_rows.rows)
        }
        == wall_lookup
        and polarization_denominator == 2
        and all(value.q == 1 for value in polarization)
        and sp.gcd_list(polarization_coordinates) == 1
        and (polarization * n_gram * polarization.T)[0] == 1540
    ):
        raise AssertionError("fixed strict positive chamber class failed")

    certificate: dict[str, object] = {
        "schema": "a18_c06_h13_31_k3_order4_lift_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "ambient_picard_lattice": {
            "invariants": [18, 4, 0],
            "rank": 18,
            "signature": [1, 17],
            "determinant": -16,
            "gram_sha256": matrix_digest(n_gram),
            "selected_action_order": 2,
            "selected_action_fixed_rank": 14,
            "selected_action_sha256": matrix_digest(n_action),
        },
        "discriminant_gluing": {
            "picard_half_vector_basis_bits": [list(vector) for vector in discr_basis],
            "picard_discriminant_action": [
                [int(value) % 2 for value in discr_action.row(index)]
                for index in range(4)
            ],
            "transcendental_lattice": "U(2)+U(2)",
            "transcendental_gram": matrix_payload(t_gram),
            "order_four_transcendental_action": matrix_payload(t_action),
            "action_square_is_minus_identity": True,
            "anti_isometry_source_to_transcendental_mod_two": matrix_payload(
                discriminant_identification
            ),
            "quadratic_form_table": q_table,
            "actions_equivariant": True,
        },
        "k3_lattice": {
            "rank": 22,
            "signature": [3, 19],
            "even": True,
            "unimodular_determinant": -1,
            "glued_gram_sha256": matrix_digest(k3_gram),
            "glued_basis_sha256": matrix_digest(k3_basis),
            "glued_action_sha256": matrix_digest(k3_action),
            "action_order": 4,
            "action_characteristic_polynomial": "(x-1)^14*(x+1)^4*(x^2+1)^2",
        },
        "period_and_torelli": {
            "period_in_T_tensor_C": ["1", "sqrt(2)", "i", "i*sqrt(2)"],
            "eigenvalue": "i",
            "period_square": "0",
            "period_positive_norm": "8*sqrt(2)",
            "no_nonzero_rational_T_class_is_orthogonal": (
                "for (a,b,c,d) in Q^4, orthogonality gives a*sqrt(2)+b=0 "
                "and c*sqrt(2)+d=0"
            ),
            "resulting_picard_lattice_is_exactly_N": True,
            "fixed_rational_chamber_class_square": "385",
            "fixed_rational_chamber_minimum_source_wall_pairing": str(
                min(chamber_wall_pairings)
            ),
            "primitive_integral_polarization_scale": polarization_denominator,
            "primitive_integral_polarization_coordinates": polarization_coordinates,
            "primitive_integral_polarization_square": "1540",
            "source_48_wall_chamber_is_preserved": True,
            "kahler_chamber_reason": (
                "the fixed class is strict in the complete mixed-root chamber; "
                "that chamber lies inside one (-2)-Weyl chamber, which the action "
                "therefore preserves, so the primitive integral double can be ample"
            ),
            "global_torelli_conclusion": (
                "a marked complex projective K3 surface with NS=N carries a genuine "
                "non-symplectic automorphism of order four inducing the C06 involution on NS"
            ),
        },
        "hyperbolic_fixed_section": {
            "id": "A18_4_0_C06",
            "dimension": 13,
            "facet_count": 31,
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
        "geometric_conclusion": (
            "The 31-wall finite-volume H13 chamber is realized as the fixed divisor "
            "section of the complete 48-wall (18,4,0) K3 reflection chamber under "
            "a genuine order-four K3 automorphism."
        ),
        "novelty_scope": (
            "new project-level geometric lift of an already certified chamber; not a "
            "new abstract polytope, first dimension population, lattice, or arithmetic class"
        ),
        "nef_scope_caution": (
            "the 48-wall source is Alexeev's complete mixed-root reflection chamber, "
            "not asserted here to be the literal nef cone of the K3 surface"
        ),
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "all_checks_pass": True}, indent=2))


if __name__ == "__main__":
    main()
