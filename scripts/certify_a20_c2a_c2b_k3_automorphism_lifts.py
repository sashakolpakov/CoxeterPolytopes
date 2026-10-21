#!/usr/bin/env python3
"""Exact K3 automorphism lifts of the two high A20 fixed sections.

Alexeev's complete 50-wall chamber for the even 2-elementary lattice
``N=(20,2,1)`` has two involution classes with fixed rank above ten.  This
certificate glues their Picard actions to the positive complement
``T=<2>+<2>``.  Class C01 requires an order-four rotation of T, while C02
glues to ``-1_T``.  The resulting Hodge isometries preserve a K3 Kahler
chamber, so Strong Torelli realizes the H12/19 and H11/20 fixed sections by
genuine non-symplectic K3 automorphisms.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_a18_c06_h13_31_k3_order4_lift import (
    discriminant_action,
    gf2_rank,
    matrix_digest,
    matrix_payload,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_wall_symmetry_inventory.json"
)
CONSTRUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_fixed_section_census.json"
)
VINBERG = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_literal_vinberg_census.json"
)
LINEAGE = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_lineage_bounds.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "a20_c2a_c2b_k3_automorphism_lifts.json"
)

EXPECTED_SOURCE_SHA256 = {
    "inventory": "20a257e60b768571a7e7c8aa64e60350bdf7940f5d221cfe031fc2e6ad879e09",
    "construction": "0294897befa7357f485d45df17d5e3e203656a3605c9c3231ab362d3a30754a6",
    "vinberg": "1ca27c1287a05824fcbede3d13a0351ea9f3bc7889463bb024229faaf847a283",
    "lineage": "4b2ac553eba574168b49654bbe9c101e1129d9436d72eabe5229aa716e0bde12",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def two_elementary_basis(gram: sp.Matrix, expected_length: int) -> tuple[tuple[int, ...], ...]:
    """Return basis bits ``v`` for the discriminant half-vectors ``v/2``."""

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


def q_times_four_mod_eight(
    coordinates: Sequence[int],
    half_vector_basis: Sequence[Sequence[int]],
    gram: sp.Matrix,
) -> int:
    """Encode the discriminant value q in Q/2Z as 4q modulo eight."""

    bits = tuple(
        sum(
            coordinates[index] * half_vector_basis[index][column]
            for index in range(len(half_vector_basis))
        )
        & 1
        for column in range(gram.rows)
    )
    return (
        sum(
            bits[i] * int(gram[i, j]) * bits[j]
            for i in range(gram.rows)
            for j in range(gram.cols)
        )
        % 8
    )


def glue_k3_lattice(
    n_gram: sp.Matrix,
    n_action: sp.Matrix,
    n_discriminant_basis: Sequence[Sequence[int]],
    t_gram: sp.Matrix,
    t_action: sp.Matrix,
    identification: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    glue_rows = sp.Matrix(
        [
            [
                *n_discriminant_basis[index],
                *map(int, identification.row(index)),
            ]
            for index in range(len(n_discriminant_basis))
        ]
    )
    doubled_rows = (2 * sp.eye(22)).col_join(glue_rows)
    k3_basis = hermite_normal_form(doubled_rows.T).T / 2
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
    ):
        raise AssertionError("the equivariant graph glue is not a K3 lattice action")
    return k3_basis, k3_gram, k3_action


def build_certificate() -> dict[str, Any]:
    source_hashes = {
        "inventory": sha256(INVENTORY),
        "construction": sha256(CONSTRUCTION),
        "vinberg": sha256(VINBERG),
        "lineage": sha256(LINEAGE),
    }
    if source_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned A20 source certificate changed")

    inventory_data = json.loads(INVENTORY.read_text(encoding="utf-8"))
    construction_data = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg_data = json.loads(VINBERG.read_text(encoding="utf-8"))
    lineage_data = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        construction_data["summary"]["survivor_count_requiring_volume_audit"] == 3
        and vinberg_data["summary"]["finite_volume_noncompact_count"] == 3
        and vinberg_data["summary"]["bad_edge_count"] == 0
        and vinberg_data["summary"]["bad_parabolic_neighborhood_count"] == 0
        and lineage_data["all_checks_pass"]
    ):
        raise AssertionError("an upstream A20 certificate is not green")

    inventory = next(
        row for row in inventory_data["examples"] if row["key"] == [20, 2, 1]
    )
    n_gram = sp.Matrix(inventory["ambient"]["gram"])
    wall_rows = sp.Matrix(inventory["walls"]["coordinates_in_leading_root_basis"])
    if not (
        n_gram.rows == 20
        and n_gram.det() == -4
        and exact_inertia(n_gram) == (1, 19, 0)
        and wall_rows.rows == 50
    ):
        raise AssertionError("the ambient (20,2,1) K3 lattice changed")

    discriminant_basis = two_elementary_basis(n_gram, 2)
    t_gram = sp.diag(2, 2)
    identification = sp.Matrix([[0, 1], [1, 1]])
    if not (
        t_gram.det() == 4
        and exact_inertia(t_gram) == (2, 0, 0)
        and int(identification.det()) % 2 == 1
    ):
        raise AssertionError("the positive transcendental complement changed")

    q_table: list[dict[str, object]] = []
    t_half_basis = ((1, 0), (0, 1))
    for coordinates in itertools.product(range(2), repeat=2):
        target = tuple(
            int(value) % 2
            for value in (sp.Matrix([coordinates]) * identification).row(0)
        )
        q_n = q_times_four_mod_eight(coordinates, discriminant_basis, n_gram)
        q_t = q_times_four_mod_eight(target, t_half_basis, t_gram)
        if (q_n + q_t) % 8:
            raise AssertionError("the displayed discriminant graph is not isotropic")
        q_table.append(
            {
                "source": list(coordinates),
                "target": list(target),
                "4q_N_mod_8": q_n,
                "4q_T_mod_8": q_t,
            }
        )

    configurations = {
        "A20_2_1_C01": {
            "ledger_id": "A20-C2a-H12-19",
            "dimension": 12,
            "facets": 19,
            "fixed_rank": 13,
            "automorphism_order": 4,
            "t_action": sp.Matrix([[0, -1], [1, 0]]),
            "period_eigenvalue": "i",
            "charpoly": "(x-1)^13*(x+1)^7*(x^2+1)",
        },
        "A20_2_1_C02": {
            "ledger_id": "A20-C2b-H11-20",
            "dimension": 11,
            "facets": 20,
            "fixed_rank": 12,
            "automorphism_order": 2,
            "t_action": -sp.eye(2),
            "period_eigenvalue": "-1",
            "charpoly": "(x-1)^12*(x+1)^10",
        },
    }

    cases: list[dict[str, object]] = []
    common_k3_gram_digest: str | None = None
    common_k3_basis_digest: str | None = None
    wall_lookup = {tuple(wall_rows.row(index)) for index in range(wall_rows.rows)}
    x = sp.Symbol("lambda")
    for case_id, config in configurations.items():
        construction = next(
            row for row in construction_data["cases"] if row["id"] == case_id
        )
        volume = next(row for row in vinberg_data["cases"] if row["id"] == case_id)
        if not (
            all(construction["checks"].values())
            and construction["case"]["subgroup_order"] == 2
            and construction["case"]["fixed_rank"] == config["fixed_rank"]
            and construction["exact_reduction"]["active_hyperbolic_facet_count"]
            == config["facets"]
            and volume["decision"]["classification"]
            == "finite_volume_noncompact"
            and volume["decision"]["finite_volume"] is True
            and volume["decision"]["compact"] is False
            and all(volume["checks"].values())
            and volume["literal_endpoint_audit"]["bad_edge_count"] == 0
            and volume["literal_endpoint_audit"][
                "bad_parabolic_neighborhood_count"
            ]
            == 0
        ):
            raise AssertionError(f"the upstream {case_id} chamber audit changed")

        n_action = sp.Matrix(
            construction["case"]["representative_nonidentity_primal_actions"][0]
        )
        t_action = config["t_action"]
        n_discriminant_action = discriminant_action(n_action, discriminant_basis)
        if not (
            n_action**2 == sp.eye(20)
            and n_action * n_gram * n_action.T == n_gram
            and len((n_action - sp.eye(20)).nullspace()) == config["fixed_rank"]
            and t_action * t_gram * t_action.T == t_gram
            and (
                n_discriminant_action * identification
            ).applyfunc(lambda value: int(value) % 2)
            == (identification * t_action).applyfunc(lambda value: int(value) % 2)
        ):
            raise AssertionError(f"the {case_id} equivariant action data failed")

        k3_basis, k3_gram, k3_action = glue_k3_lattice(
            n_gram,
            n_action,
            discriminant_basis,
            t_gram,
            t_action,
            identification,
        )
        if config["automorphism_order"] == 4:
            expected_charpoly = (x - 1) ** 13 * (x + 1) ** 7 * (x**2 + 1)
            order_check = k3_action**4 == sp.eye(22) and k3_action**2 != sp.eye(22)
        else:
            expected_charpoly = (x - 1) ** 12 * (x + 1) ** 10
            order_check = k3_action**2 == sp.eye(22) and k3_action != sp.eye(22)
        if not (
            order_check
            and sp.factor(k3_action.charpoly().as_expr()) == expected_charpoly
        ):
            raise AssertionError(f"the glued {case_id} action order changed")

        k3_gram_digest = matrix_digest(k3_gram)
        k3_basis_digest = matrix_digest(k3_basis)
        if common_k3_gram_digest is None:
            common_k3_gram_digest = k3_gram_digest
            common_k3_basis_digest = k3_basis_digest
        elif not (
            k3_gram_digest == common_k3_gram_digest
            and k3_basis_digest == common_k3_basis_digest
        ):
            raise AssertionError("the two actions did not use the same marked K3 lattice")

        period = sp.Matrix([[1, sp.I]])
        expected_eigenvalue = sp.I if config["automorphism_order"] == 4 else -1
        if not (
            period * t_action == expected_eigenvalue * period
            and (period * t_gram * period.T)[0] == 0
            and (period * t_gram * sp.conjugate(period).T)[0] == 4
        ):
            raise AssertionError(f"the {case_id} period is not a positive eigenperiod")

        fixed_basis = sp.Matrix(
            construction["saturated_fixed_lattice"]["basis_primal_rows"]
        )
        fixed_coordinates = sp.Matrix(
            [
                sp.Rational(value)
                for value in construction["strict_fixed_interior"][
                    "fixed_lattice_coordinates"
                ]
            ]
        )
        interior = fixed_coordinates.T * fixed_basis
        pairings = list(interior * n_gram * wall_rows.T)
        denominator = int(sp.ilcm(*(sp.Rational(value).q for value in interior)))
        polarization = denominator * interior
        polarization_coordinates = [int(value) for value in polarization]
        if not (
            interior * n_action == interior
            and all(value > 0 for value in pairings)
            and min(pairings) == 1
            and {
                tuple(wall_rows.row(index) * n_action)
                for index in range(wall_rows.rows)
            }
            == wall_lookup
            and all(value.q == 1 for value in polarization)
            and sp.gcd_list(polarization_coordinates) == 1
        ):
            raise AssertionError(f"the {case_id} invariant ample chamber gate failed")

        fixed_square = sp.simplify(interior * n_gram * interior.T)[0]
        polarization_square = sp.simplify(
            polarization * n_gram * polarization.T
        )[0]
        cases.append(
            {
                "id": case_id,
                "ledger_id": config["ledger_id"],
                "picard_action_order": 2,
                "picard_action_fixed_rank": config["fixed_rank"],
                "picard_action_sha256": matrix_digest(n_action),
                "picard_discriminant_action_mod_two": [
                    [int(value) % 2 for value in n_discriminant_action.row(index)]
                    for index in range(2)
                ],
                "transcendental_action": matrix_payload(t_action),
                "transcendental_action_order": config["automorphism_order"],
                "period": ["1", "i"],
                "period_eigenvalue": config["period_eigenvalue"],
                "period_square": "0",
                "period_positive_norm": "4",
                "resulting_picard_lattice_is_exactly_N": True,
                "k3_action_order": config["automorphism_order"],
                "k3_action_characteristic_polynomial": config["charpoly"],
                "k3_action_sha256": matrix_digest(k3_action),
                "fixed_rational_chamber_class_square": str(fixed_square),
                "primitive_invariant_polarization_scale": denominator,
                "primitive_invariant_polarization_coordinates": polarization_coordinates,
                "primitive_invariant_polarization_square": str(polarization_square),
                "minimum_source_wall_pairing": str(min(pairings)),
                "source_50_wall_chamber_preserved": True,
                "fixed_neron_severi_rank": config["fixed_rank"],
                "fixed_neron_severi_is_archived_fixed_lattice": True,
                "hyperbolic_fixed_section": {
                    "dimension": config["dimension"],
                    "facet_count": config["facets"],
                    "finite_volume": True,
                    "compact": False,
                    "finite_vertex_count": volume["literal_endpoint_audit"][
                        "finite_vertex_count"
                    ],
                    "ideal_vertex_count": volume["literal_endpoint_audit"][
                        "ideal_vertex_count"
                    ],
                    "elliptic_edge_count": volume["literal_endpoint_audit"][
                        "edge_count"
                    ],
                    "bad_edge_count": 0,
                    "bad_parabolic_neighborhood_count": 0,
                },
            }
        )

    certificate: dict[str, Any] = {
        "schema": "a20_c2a_c2b_k3_automorphism_lifts_v1",
        "all_checks_pass": True,
        "source_file_sha256": source_hashes,
        "ambient_picard_lattice": {
            "invariants": [20, 2, 1],
            "rank": 20,
            "signature": [1, 19],
            "determinant": -4,
            "gram_sha256": matrix_digest(n_gram),
            "source_wall_count": 50,
        },
        "discriminant_gluing": {
            "picard_half_vector_basis_bits": [
                list(vector) for vector in discriminant_basis
            ],
            "transcendental_lattice": "<2>+<2>",
            "transcendental_gram": matrix_payload(t_gram),
            "anti_isometry_source_to_transcendental_mod_two": matrix_payload(
                identification
            ),
            "quadratic_form_table": q_table,
            "both_actions_equivariant": True,
        },
        "common_marked_k3_lattice": {
            "rank": 22,
            "signature": [3, 19],
            "even": True,
            "unimodular_determinant": -1,
            "glued_gram_sha256": common_k3_gram_digest,
            "glued_basis_sha256": common_k3_basis_digest,
            "common_period": ["1", "i"],
            "no_nonzero_rational_T_class_is_orthogonal": (
                "for (a,b) in Q^2, orthogonality is 2a+2ib=0"
            ),
            "resulting_picard_lattice_is_exactly_N": True,
        },
        "cases": cases,
        "torelli_conclusion": (
            "the same marked singular projective K3 surface with NS=(20,2,1) "
            "admits a non-symplectic order-four automorphism realizing C01 and "
            "a non-symplectic involution realizing C02"
        ),
        "kahler_chamber_reason": (
            "each action preserves the complete 50-wall mixed reflection chamber "
            "and a strict positive class; that chamber lies in one (-2)-Weyl "
            "chamber, so the action preserves an ample chamber"
        ),
        "geometric_conclusion": (
            "the finite-volume noncompact H12/19 and H11/20 chambers are genuine "
            "automorphism-fixed divisor/reflection-chamber sections on a K3 surface"
        ),
        "novelty_scope": (
            "new project-level geometric realizations of already certified chambers; "
            "not new abstract polytopes, lattices, dimensions, or arithmetic classes"
        ),
        "nef_scope_caution": (
            "the 50-wall source is Alexeev's complete mixed-root reflection chamber, "
            "not asserted to be the literal nef cone"
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
            raise AssertionError("archived A20 K3-lift certificate changed")
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
                "certificate_payload_sha256": certificate[
                    "certificate_payload_sha256"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
