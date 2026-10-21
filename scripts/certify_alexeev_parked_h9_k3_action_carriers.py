#!/usr/bin/env python3
"""Promote four parked H9 sections to genuine symplectic K3 carriers.

The five-seed certificate starts with involutions of abstract hyperbolic
lattices.  This replay keeps the original actions visible.  For the four
rank-four algebraic coinvariants it identifies the coinvariant with
``D4(-2)``, completes it primitively to ``E8(-2)``, and constructs an exact
rank-18 Picard pushout.  A displayed isotropic discriminant subgroup then
gives an even determinant-16 Picard lattice whose coinvariant is still
``E8(-2)``.  Gluing its ``(Z/4)^2`` discriminant form to an explicit even
signature-``(2,2)`` determinant-16 lattice produces the K3 lattice and a
Nikulin involution.

The rank-six fifth row fails this particular symplectic completion: half of
its coinvariant is odd, so it cannot embed in the even E8 lattice.  This is a
route-specific obstruction, not a no-surface-action theorem.

The script also records the independent pointwise-invariant alternative:
the five maximal reflection-stable fixed-lattice overclosures are classical
2-elementary invariant lattices of nonsymplectic K3 involutions.  That
alternative does not realize the original ambient involution.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    same_row_lattice,
    saturated_row_kernel,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _overlattice_basis,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
PARKED = ROOT / "data/enumerations/alexeev_parked_rank10_h9_examples.json"
DENOMINATOR = (
    ROOT / "data/enumerations/alexeev_ten_small_zero_high_denominator.json"
)
LINEAGE = (
    ROOT
    / "data/enumerations/alexeev_parked_rank10_h9_lineage_overclosures.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/alexeev_parked_h9_k3_action_carriers.json"
)

EXPECTED_SOURCE_SHA256 = {
    "parked": "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430",
    "denominator": "9c518c02083588cffaa15f38f5887b6f9c2bc109eed7d93611b065a0e0c6fdff",
    "lineage": "972f1eb50940046c90c3770ff356a5a923477dec2e87d1d7acbffd85bf785bcd",
}


# PARI/GP qfisom witnesses, replayed below in exact integer arithmetic.
# Each P satisfies P^T D4 P = -(coinvariant Gram)/2.
SPECS: dict[str, dict[str, Any]] = {
    "A14_8_0_C01": {
        "d4_isometry": [
            [0, 0, -1, 1],
            [0, -2, 0, 1],
            [0, -1, 1, 0],
            [1, -1, 0, 0],
        ],
        "raw_to_picard_index": 8,
        "isotropic_generators_reduced_smith": [
            [1, 1, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0, 1, 0, 0],
        ],
    },
    "A14_8_1_C02": {
        "d4_isometry": [
            [0, -1, -1, 1],
            [0, -1, -2, 2],
            [0, -2, 0, 1],
            [1, 0, -1, 0],
        ],
        "raw_to_picard_index": 8,
        "isotropic_generators_reduced_smith": [
            [1, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 1, 0, 0],
        ],
    },
    "A14_10_1_C01": {
        "d4_isometry": [
            [0, -3, 4, 0],
            [0, 4, -6, -1],
            [0, 1, -1, 0],
            [1, -2, 1, 0],
        ],
        "raw_to_picard_index": 16,
        "isotropic_generators_reduced_smith": [
            [1, 0, 1, 0, 1, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0, 1, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
        ],
    },
    "A14_10_1_C02": {
        "d4_isometry": [
            [-1, -1, 0, 1],
            [-3, -3, 0, 2],
            [-1, -1, -1, 2],
            [-2, -1, -1, 1],
        ],
        "raw_to_picard_index": 16,
        "isotropic_generators_reduced_smith": [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
        ],
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [str(matrix[row, col]) for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _digest(matrix: sp.Matrix) -> str:
    return hashlib.sha256(
        json.dumps(_payload(matrix), separators=(",", ":")).encode()
    ).hexdigest()


def _smith(gram: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(gram, domain=ZZ)
    return [abs(int(diagonal[index, index])) for index in range(gram.rows)]


def _row_lattice_basis(generators: sp.Matrix) -> sp.Matrix:
    denominator = int(
        sp.ilcm(*(sp.Rational(value).q for value in generators))
    )
    integral = (denominator * generators).applyfunc(int)
    return hermite_normal_form(integral.T).T / denominator


def _binary_subgroup(
    generators: Sequence[Sequence[int]], length: int
) -> tuple[tuple[int, ...], ...]:
    result = {tuple([0] * length)}
    for generator in generators:
        if len(generator) != length:
            raise AssertionError("a frozen discriminant generator changed length")
        old = tuple(result)
        result.update(
            tuple(left[index] ^ int(generator[index]) for index in range(length))
            for left in old
        )
    return tuple(sorted(result))


def _full_smith_cosets(
    reduced: Iterable[Sequence[int]], moduli: Sequence[int]
) -> list[list[int]]:
    nontrivial = [index for index, modulus in enumerate(moduli) if modulus > 1]
    rows = []
    for vector in reduced:
        if len(vector) != len(nontrivial):
            raise AssertionError("reduced Smith-coordinate length changed")
        row = [0] * len(moduli)
        for position, index in enumerate(nontrivial):
            row[index] = int(vector[position])
        rows.append(row)
    return rows


def _e8_data() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    e8 = 2 * sp.eye(8)
    for index in range(6):
        e8[index, index + 1] = -1
        e8[index + 1, index] = -1
    e8[2, 7] = e8[7, 2] = -1
    embedding = sp.zeros(8, 4)
    for column, row in enumerate((1, 2, 3, 7)):
        embedding[row, column] = 1
    d4 = embedding.T * e8 * embedding
    expected_d4 = sp.Matrix(
        [
            [2, -1, 0, 0],
            [-1, 2, -1, -1],
            [0, -1, 2, 0],
            [0, -1, 0, 2],
        ]
    )
    if not (
        e8.det() == 1
        and exact_inertia(e8) == (8, 0, 0)
        and d4 == expected_d4
        and all(e8[index, index] % 2 == 0 for index in range(8))
    ):
        raise AssertionError("the standard E8/D4 model changed")
    return e8, embedding, d4


def _isotropic_audit(
    gram: sp.Matrix,
    moduli: Sequence[int],
    left_inverse: sp.Matrix,
    cosets: Sequence[Sequence[int]],
) -> dict[str, Any]:
    _, _, discriminant_gram = _discriminant_data(gram)
    vectors = [sp.Matrix(row) for row in cosets]
    closed = {
        tuple(
            (int(left[index]) + int(right[index])) % moduli[index]
            for index in range(len(moduli))
        )
        for left in cosets
        for right in cosets
    } == {tuple(map(int, row)) for row in cosets}
    bilinear_integral = all(
        (left.T * discriminant_gram * right)[0].q == 1
        for left in vectors
        for right in vectors
    )
    even = all(
        sp.factor((vector.T * discriminant_gram * vector)[0]) % 2 == 0
        for vector in vectors
    )
    if not (closed and bilinear_integral and even):
        raise AssertionError("the frozen Picard overclosure subgroup is not even isotropic")
    return {
        "order": len(cosets),
        "closed_under_addition": closed,
        "bilinear_integral": bilinear_integral,
        "quadratic_values_even": even,
        "cosets_in_smith_coordinates": [list(map(int, row)) for row in cosets],
        "left_smith_inverse_sha256": _digest(left_inverse),
    }


def _discriminant_action_is_identity(
    gram: sp.Matrix, row_action: sp.Matrix
) -> bool:
    moduli, left_inverse, _ = _discriminant_data(gram)
    column_action = row_action.T
    action = (
        left_inverse.inv()
        * gram
        * column_action
        * gram.inv()
        * left_inverse
    )
    if not all(value.q == 1 for value in action):
        return False
    return all(
        (int(action[row, col]) - int(row == col)) % moduli[row] == 0
        for row in range(gram.rows)
        for col in range(gram.cols)
    )


def _find_anti_isometry_z4_squared(
    n_gram: sp.Matrix, t_gram: sp.Matrix
) -> tuple[sp.Matrix, dict[str, Any]]:
    n_moduli, n_left_inverse, n_discriminant = _discriminant_data(n_gram)
    t_moduli, t_left_inverse, t_discriminant = _discriminant_data(t_gram)
    n_indices = [index for index, modulus in enumerate(n_moduli) if modulus > 1]
    t_indices = [index for index, modulus in enumerate(t_moduli) if modulus > 1]
    if not (
        [n_moduli[index] for index in n_indices] == [4, 4]
        and [t_moduli[index] for index in t_indices] == [4, 4]
    ):
        raise AssertionError("the final NS/T discriminant is not (Z/4)^2")
    n_form = n_discriminant.extract(n_indices, n_indices)
    t_form = t_discriminant.extract(t_indices, t_indices)
    elements = tuple(itertools.product(range(4), repeat=2))
    witness = None
    for entries in itertools.product(range(4), repeat=4):
        matrix = sp.Matrix(2, 2, entries)
        if int(matrix.det()) % 2 == 0:
            continue
        good = True
        for element in elements:
            source = sp.Matrix(element)
            target = (matrix * source).applyfunc(lambda value: int(value) % 4)
            if sp.factor(
                (source.T * n_form * source)[0]
                + (target.T * t_form * target)[0]
            ) % 2 != 0:
                good = False
                break
        if good:
            witness = matrix
            break
    if witness is None:
        n_histogram: dict[str, int] = {}
        t_histogram: dict[str, int] = {}
        for element in elements:
            vector = sp.Matrix(element)
            q_n = str(sp.factor((vector.T * n_form * vector)[0]) % 2)
            q_t = str(sp.factor((vector.T * t_form * vector)[0]) % 2)
            n_histogram[q_n] = n_histogram.get(q_n, 0) + 1
            t_histogram[q_t] = t_histogram.get(q_t, 0) + 1
        raise AssertionError(
            "no (Z/4)^2 discriminant anti-isometry was found: "
            f"NS={n_histogram}, T={t_histogram}"
        )
    table = []
    for element in elements:
        source = sp.Matrix(element)
        target = (witness * source).applyfunc(lambda value: int(value) % 4)
        q_n = sp.factor((source.T * n_form * source)[0]) % 2
        q_t = sp.factor((target.T * t_form * target)[0]) % 2
        if (q_n + q_t) % 2 != 0:
            raise AssertionError("the frozen anti-isometry table failed")
        table.append(
            {
                "source": list(map(int, element)),
                "target": [int(value) for value in target],
                "q_NS_mod_2": str(q_n),
                "q_T_mod_2": str(q_t),
            }
        )
    return witness, {
        "NS_nontrivial_smith_indices": n_indices,
        "T_nontrivial_smith_indices": t_indices,
        "NS_left_smith_inverse": n_left_inverse,
        "T_left_smith_inverse": t_left_inverse,
        "table": table,
    }


def _glue_k3(
    n_gram: sp.Matrix,
    n_action: sp.Matrix,
    t_gram: sp.Matrix,
    anti_isometry: sp.Matrix,
    anti_data: dict[str, Any],
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    n_moduli, n_left_inverse, _ = _discriminant_data(n_gram)
    t_moduli, t_left_inverse, _ = _discriminant_data(t_gram)
    n_indices = anti_data["NS_nontrivial_smith_indices"]
    t_indices = anti_data["T_nontrivial_smith_indices"]
    glue_rows = []
    for column in range(2):
        n_coordinate = sp.zeros(18, 1)
        n_coordinate[n_indices[column]] = 1
        t_coordinate_small = anti_isometry * sp.eye(2).col(column)
        t_coordinate = sp.zeros(4, 1)
        for position, index in enumerate(t_indices):
            t_coordinate[index] = int(t_coordinate_small[position]) % 4
        n_vector = n_gram.inv() * n_left_inverse * n_coordinate
        t_vector = t_gram.inv() * t_left_inverse * t_coordinate
        glue_rows.append(list(n_vector) + list(t_vector))
    generators = sp.eye(22).col_join(sp.Matrix(glue_rows))
    k3_basis = _row_lattice_basis(generators)
    orthogonal_gram = sp.diag(n_gram, t_gram)
    k3_gram = sp.simplify(k3_basis * orthogonal_gram * k3_basis.T)
    orthogonal_action = sp.diag(n_action, sp.eye(4))
    k3_action = sp.simplify(k3_basis * orthogonal_action * k3_basis.inv())
    if not (
        abs(k3_basis.det()) == sp.Rational(1, 16)
        and all(value.q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(value.q == 1 for value in k3_action)
        and k3_action**2 == sp.eye(22)
        and k3_action * k3_gram * k3_action.T == k3_gram
    ):
        raise AssertionError("the determinant-16 graph glue did not give a K3 action")
    return k3_basis, k3_gram, k3_action, orthogonal_action


def build_certificate() -> dict[str, Any]:
    actual_hashes = {
        "parked": _sha256(PARKED),
        "denominator": _sha256(DENOMINATOR),
        "lineage": _sha256(LINEAGE),
    }
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned parked-H9 input changed")
    parked = json.loads(PARKED.read_text(encoding="utf-8"))
    denominator = json.loads(DENOMINATOR.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        parked["summary"]["finite_volume_noncompact_count"] == 5
        and lineage["all_checks_pass"] is True
    ):
        raise AssertionError("an upstream parked-H9 certificate is not green")

    ambient_lookup = {
        tuple(record["key"]): _matrix(record["ambient"]["gram"])
        for record in denominator["examples"]
    }
    case_lookup = {record["id"]: record for record in parked["cases"]}
    e8, d4_embedding, d4 = _e8_data()
    sqrt_two = sp.sqrt(2)
    exotic_t = sp.Matrix(
        [
            [2, 1, 1, -2],
            [1, 6, 1, 1],
            [1, 1, -6, 1],
            [-2, 1, 1, 2],
        ]
    )
    exotic_witt_basis = sp.Matrix.hstack(
        sp.Matrix([-1, 0, 0, -1]),
        sp.Matrix([-5, -1, -5, 2]),
        sp.Matrix([4, -1, 1, 5]),
        sp.Matrix([8, -3, 3, 5]),
    )
    split_t = sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 4],
            [0, 0, 4, 0],
        ]
    )
    transcendental_candidates = (
        {
            "name": "T16a",
            "description": "explicit nonsplit-discriminant determinant-16 lattice",
            "gram": exotic_t,
            "witt_basis": exotic_witt_basis,
            "witt_gram": sp.Matrix(
                [
                    [0, 12, 0, 0],
                    [12, 0, 0, 0],
                    [0, 0, 0, -12],
                    [0, 0, -12, 0],
                ]
            ),
            "period_witt": sp.Matrix(
                [1, sqrt_two, sp.I, -sp.I * sqrt_two]
            ),
            "period_text": ["1", "sqrt(2)", "i", "-i*sqrt(2)"],
            "positive_norm": 48 * sqrt_two,
        },
        {
            "name": "U+U(4)",
            "description": "U + U(4)",
            "gram": split_t,
            "witt_basis": sp.eye(4),
            "witt_gram": split_t,
            "period_witt": sp.Matrix(
                [1, 4 * sqrt_two, sp.I, sp.I * sqrt_two]
            ),
            "period_text": ["1", "4*sqrt(2)", "i", "i*sqrt(2)"],
            "positive_norm": 16 * sqrt_two,
        },
    )
    if not all(
        candidate["gram"].det() == 16
        and exact_inertia(candidate["gram"]) == (2, 2, 0)
        and _smith(candidate["gram"]) == [1, 1, 4, 4]
        and candidate["witt_basis"].T
        * candidate["gram"]
        * candidate["witt_basis"]
        == candidate["witt_gram"]
        for candidate in transcendental_candidates
    ):
        raise AssertionError("a determinant-16 transcendental lattice changed")

    strong_records = []
    for identifier, spec in SPECS.items():
        case = case_lookup[identifier]
        ambient_gram = ambient_lookup[tuple(case["ambient_key"])]
        ambient_action = _matrix(
            case["section"]["case"]["representative_nonidentity_primal_actions"][0]
        )
        rank = ambient_gram.rows
        if not (
            rank == 14
            and ambient_action**2 == sp.eye(rank)
            and ambient_action * ambient_gram * ambient_action.T == ambient_gram
            and case["section"]["case"]["fixed_rank"] == 10
        ):
            raise AssertionError(f"the archived ambient involution changed for {identifier}")
        fixed_basis = saturated_row_kernel(ambient_action - sp.eye(rank))
        coin_basis = saturated_row_kernel(ambient_action + sp.eye(rank))
        fixed_gram = fixed_basis * ambient_gram * fixed_basis.T
        coin_gram = coin_basis * ambient_gram * coin_basis.T
        d4_isometry = sp.Matrix(spec["d4_isometry"])
        if not (
            coin_basis.rows == 4
            and coin_gram.det() == 64
            and exact_inertia(coin_gram) == (0, 4, 0)
            and d4_isometry.det() in (-1, 1)
            and d4_isometry.T * d4 * d4_isometry == -coin_gram / 2
        ):
            raise AssertionError(f"the D4(-2) coinvariant gate failed for {identifier}")

        # Work first in F_Q + E8(-2)_Q.  The original ambient lattice is
        # inserted through its exact F+K glue, with K mapped to D4(-2) in E8(-2).
        fixed_coin_basis = fixed_basis.col_join(coin_basis)
        target_map = sp.zeros(rank, 18)
        target_map[:10, :10] = sp.eye(10)
        target_map[10:, 10:] = d4_isometry.T * d4_embedding.T
        ambient_in_orthogonal_coordinates = sp.simplify(
            fixed_coin_basis.inv() * target_map
        )
        orthogonal_gram = sp.diag(fixed_gram, -2 * e8)
        if not (
            ambient_in_orthogonal_coordinates
            * orthogonal_gram
            * ambient_in_orthogonal_coordinates.T
            == ambient_gram
        ):
            raise AssertionError("the D4-to-E8 ambient embedding changed")
        raw_generators = ambient_in_orthogonal_coordinates.col_join(
            sp.zeros(8, 10).row_join(sp.eye(8))
        )
        raw_basis = _row_lattice_basis(raw_generators)
        raw_gram = sp.simplify(raw_basis * orthogonal_gram * raw_basis.T)
        orthogonal_action = sp.diag(sp.eye(10), -sp.eye(8))
        raw_action = sp.simplify(
            raw_basis * orthogonal_action * raw_basis.inv()
        )
        if not (
            all(value.q == 1 for value in raw_gram)
            and all(raw_gram[index, index] % 2 == 0 for index in range(18))
            and raw_gram.det() == fixed_gram.det()
            and all(value.q == 1 for value in raw_action)
            and _discriminant_action_is_identity(raw_gram, raw_action)
        ):
            raise AssertionError("the raw Nikulin pushout failed")

        moduli, left_inverse, _ = _discriminant_data(raw_gram)
        nontrivial_moduli = [modulus for modulus in moduli if modulus > 1]
        reduced_generators = spec["isotropic_generators_reduced_smith"]
        reduced_group = _binary_subgroup(
            reduced_generators, len(nontrivial_moduli)
        )
        cosets = _full_smith_cosets(reduced_group, moduli)
        isotropic_audit = _isotropic_audit(
            raw_gram, moduli, left_inverse, cosets
        )
        expected_index = int(spec["raw_to_picard_index"])
        if isotropic_audit["order"] != expected_index:
            raise AssertionError("the frozen Picard subgroup has wrong order")
        picard_overbasis = _overlattice_basis(
            raw_gram, left_inverse, cosets
        )
        picard_gram = sp.simplify(
            picard_overbasis.T * raw_gram * picard_overbasis
        )
        picard_action = sp.simplify(
            picard_overbasis.inv() * raw_action.T * picard_overbasis
        ).T
        picard_basis_in_orthogonal_coordinates = (
            picard_overbasis.T * raw_basis
        )
        ambient_embedding = sp.simplify(
            ambient_in_orthogonal_coordinates
            * picard_basis_in_orthogonal_coordinates.inv()
        )
        if not (
            all(value.q == 1 for value in picard_gram)
            and all(picard_gram[index, index] % 2 == 0 for index in range(18))
            and picard_gram.det() == -16
            and _smith(picard_gram) == [1] * 16 + [4, 4]
            and exact_inertia(picard_gram) == (1, 17, 0)
            and all(value.q == 1 for value in picard_action)
            and picard_action**2 == sp.eye(18)
            and picard_action * picard_gram * picard_action.T == picard_gram
            and _discriminant_action_is_identity(picard_gram, picard_action)
            and all(value.q == 1 for value in ambient_embedding)
            and ambient_embedding * picard_gram * ambient_embedding.T
            == ambient_gram
            and ambient_action * ambient_embedding
            == ambient_embedding * picard_action
        ):
            raise AssertionError(f"the final rank-18 Picard carrier failed for {identifier}")

        picard_fixed_basis = saturated_row_kernel(
            picard_action - sp.eye(18)
        )
        picard_coin_basis = saturated_row_kernel(
            picard_action + sp.eye(18)
        )
        picard_fixed_gram = (
            picard_fixed_basis * picard_gram * picard_fixed_basis.T
        )
        picard_coin_gram = (
            picard_coin_basis * picard_gram * picard_coin_basis.T
        )
        half_coin = -picard_coin_gram / 2
        if not (
            picard_fixed_basis.rows == 10
            and exact_inertia(picard_fixed_gram) == (1, 9, 0)
            and picard_coin_basis.rows == 8
            and picard_coin_gram.det() == 256
            and exact_inertia(picard_coin_gram) == (0, 8, 0)
            and all(value.q == 1 for value in half_coin)
            and all(half_coin[index, index] % 2 == 0 for index in range(8))
            and half_coin.det() == 1
            and exact_inertia(half_coin) == (8, 0, 0)
        ):
            raise AssertionError("the E8(-2) Kähler gate failed")

        selected_transcendental = None
        for candidate in transcendental_candidates:
            try:
                anti_isometry, anti_data = _find_anti_isometry_z4_squared(
                    picard_gram, candidate["gram"]
                )
            except AssertionError:
                continue
            selected_transcendental = candidate
            break
        if selected_transcendental is None:
            raise AssertionError("neither determinant-16 T discriminant form matched")
        t_gram = selected_transcendental["gram"]
        k3_basis, k3_gram, k3_action, _ = _glue_k3(
            picard_gram,
            picard_action,
            t_gram,
            anti_isometry,
            anti_data,
        )
        picard_embedding = (
            sp.eye(18).row_join(sp.zeros(18, 4)) * k3_basis.inv()
        )
        recovered_picard = saturated_row_kernel(k3_basis[:, 18:22])
        primitive_checks = {
            "embedding_integral": all(value.q == 1 for value in picard_embedding),
            "same_row_lattice": same_row_lattice(
                recovered_picard, picard_embedding
            ),
            "recovered_determinant": (
                recovered_picard * k3_gram * recovered_picard.T
            ).det()
            == -16,
        }
        if not all(primitive_checks.values()):
            raise AssertionError(
                "the rank-18 Picard lattice is not primitive in K3: "
                f"{primitive_checks}"
            )

        witt_basis = selected_transcendental["witt_basis"]
        witt_gram = selected_transcendental["witt_gram"]
        period_witt = selected_transcendental["period_witt"]
        period = (witt_basis * period_witt).T
        period_square = sp.simplify((period * t_gram * period.T)[0])
        period_positive = sp.simplify(
            (period * t_gram * sp.conjugate(period).T)[0]
        )
        rational_pairing_constraints = witt_basis.T * t_gram
        # In Witt coordinates the period is (1,sqrt(2),i,-i*sqrt(2)).
        # A rational orthogonal vector therefore has four zero coefficients
        # after splitting real/imaginary and Q+sqrt(2)Q.  The displayed
        # rational constraint matrix is nonsingular, so that vector is zero.
        if not (
            period_square == 0
            and period_positive == selected_transcendental["positive_norm"]
            and rational_pairing_constraints.det() != 0
            and sp.simplify(period * sp.eye(4) - period) == sp.zeros(1, 4)
        ):
            raise AssertionError("the generic symplectic period failed")

        cohom_fixed = saturated_row_kernel(k3_action - sp.eye(22))
        cohom_coin = saturated_row_kernel(k3_action + sp.eye(22))
        cohom_fixed_gram = cohom_fixed * k3_gram * cohom_fixed.T
        cohom_coin_gram = cohom_coin * k3_gram * cohom_coin.T
        if not (
            sp.factor(k3_action.charpoly().as_expr())
            == (sp.Symbol("lambda") - 1) ** 14
            * (sp.Symbol("lambda") + 1) ** 8
            and cohom_fixed.rows == 14
            and exact_inertia(cohom_fixed_gram) == (3, 11, 0)
            and cohom_coin.rows == 8
            and exact_inertia(cohom_coin_gram) == (0, 8, 0)
            and cohom_coin_gram.det() == 256
        ):
            raise AssertionError("the full Nikulin cohomology type failed")

        strong_records.append(
            {
                "id": identifier,
                "original_ambient_action": {
                    "rank_signature_determinant": [
                        rank,
                        [1, rank - 1],
                        int(ambient_gram.det()),
                    ],
                    "action_sha256": _digest(ambient_action),
                    "fixed_rank": fixed_basis.rows,
                    "coinvariant_rank_determinant": [
                        coin_basis.rows,
                        int(coin_gram.det()),
                    ],
                    "coinvariant_isometry": "D4(-2)",
                    "D4_isometry": _payload(d4_isometry),
                    "embedded_integrally_in_final_Picard": True,
                    "action_intertwining_exact": True,
                    "ambient_embedding_sha256": _digest(ambient_embedding),
                },
                "nikulin_completion": {
                    "primitive_coinvariant_completion": "D4(-2) subset E8(-2)",
                    "D4_coordinate_embedding_in_E8": _payload(d4_embedding),
                    "raw_pushout_determinant": int(raw_gram.det()),
                    "raw_pushout_smith": _smith(raw_gram),
                    "raw_action_trivial_on_discriminant": True,
                    "picard_overclosure_index": expected_index,
                    "isotropic_subgroup": isotropic_audit,
                },
                "final_Picard_lattice": {
                    "rank_signature_determinant": [18, [1, 17], -16],
                    "smith_diagonal_absolute": [1] * 16 + [4, 4],
                    "gram_sha256": _digest(picard_gram),
                    "action_sha256": _digest(picard_action),
                    "action_trivial_on_discriminant": True,
                    "fixed_lattice_rank_determinant": [
                        10,
                        int(picard_fixed_gram.det()),
                    ],
                    "fixed_lattice_smith_diagonal_absolute": _smith(
                        picard_fixed_gram
                    ),
                    "coinvariant_lattice": {
                        "isometry_type": "E8(-2)",
                        "rank_determinant": [8, 256],
                        "gram": _payload(picard_coin_gram),
                        "half_positive_gram_sha256": _digest(half_coin),
                        "minus_two_rootless": True,
                        "minimum_negative_norm_absolute": 4,
                        "proof": (
                            "(-coinvariant Gram)/2 is positive even unimodular "
                            "of rank 8, hence E8 and has minimum 2"
                        ),
                    },
                },
                "transcendental_and_K3_glue": {
                    "T": selected_transcendental["description"],
                    "T_gram": _payload(t_gram),
                    "anti_isometry_A_NS_to_A_T": _payload(anti_isometry),
                    "anti_isometry_table": anti_data["table"],
                    "glued_even_unimodular_signature": [3, 19],
                    "k3_basis_sha256": _digest(k3_basis),
                    "k3_gram_sha256": _digest(k3_gram),
                    "k3_action_sha256": _digest(k3_action),
                    "Picard_primitive_and_exact_for_period": True,
                    "period_in_Witt_coordinates": selected_transcendental[
                        "period_text"
                    ],
                    "Witt_basis_columns_in_T": _payload(witt_basis),
                    "period_square": "0",
                    "period_positive_norm": str(
                        selected_transcendental["positive_norm"]
                    ),
                    "period_eigenvalue": "1",
                },
                "torelli_and_geometry": {
                    "genuine_K3_automorphism": True,
                    "automorphism_type": "symplectic Nikulin involution",
                    "invariant_Kahler_chamber_exists": True,
                    "Kahler_gate": "algebraic coinvariant E8(-2) has no (-2)-roots",
                    "carrier_tier": "genuine invariant-divisor carrier",
                    "original_Coxeter_section_preserved_in_real_fixed_space": True,
                    "effective_wall_claim": False,
                    "imported_cone_equals_nef_claim": False,
                    "normal_quotient_nef_claim": False,
                },
            }
        )

    fifth = case_lookup["A16_14_1_C01"]
    fifth_gram = ambient_lookup[tuple(fifth["ambient_key"])]
    fifth_action = _matrix(
        fifth["section"]["case"]["representative_nonidentity_primal_actions"][0]
    )
    fifth_coin_basis = saturated_row_kernel(
        fifth_action + sp.eye(fifth_gram.rows)
    )
    fifth_coin_gram = fifth_coin_basis * fifth_gram * fifth_coin_basis.T
    fifth_half = -fifth_coin_gram / 2
    fifth_obstruction = {
        "id": "A16_14_1_C01",
        "coinvariant_rank_determinant": [6, int(fifth_coin_gram.det())],
        "half_coinvariant_integral": all(value.q == 1 for value in fifth_half),
        "half_coinvariant_determinant": int(fifth_half.det()),
        "half_coinvariant_even": all(
            fifth_half[index, index] % 2 == 0
            for index in range(fifth_half.rows)
        ),
        "odd_norm_witness": {
            "coordinates": [1, 0, 0, 0, 0, 0],
            "norm": int(fifth_half[0, 0]),
        },
        "symplectic_completion_to_E8_minus_2_impossible": True,
        "proof": (
            "an isometric embedding of the coinvariant in E8(-2) would, "
            "after division by 2, embed this odd lattice in the even E8 lattice"
        ),
        "scope": (
            "obstructs only the universal Nikulin symplectic-completion route; "
            "it is not a no-surface-action theorem"
        ),
    }
    if not (
        fifth_coin_basis.rows == 6
        and fifth_coin_gram.det() == 1024
        and exact_inertia(fifth_coin_gram) == (0, 6, 0)
        and fifth_obstruction["half_coinvariant_integral"] is True
        and fifth_half.det() == 16
        and fifth_obstruction["half_coinvariant_even"] is False
        and fifth_half[0, 0] % 2 == 1
    ):
        raise AssertionError("the fifth-row parity obstruction changed")

    expected_types = {
        "A14_8_0_C01": [10, 6, 0],
        "A14_8_1_C02": [10, 4, 0],
        "A14_10_1_C01": [10, 8, 0],
        "A14_10_1_C02": [10, 6, 0],
        "A16_14_1_C01": [10, 4, 0],
    }
    alternatives = []
    for record in lineage["cases"]:
        identifier = record["id"]
        carrier = record["unique_maximal_reflection_stable_overclosure"]
        if not (
            carrier["two_elementary_invariants_r_a_delta"]
            == expected_types[identifier]
            and carrier["unique"] is True
            and all(carrier["checks"].values())
        ):
            raise AssertionError("a maximal pointwise carrier changed")
        alternatives.append(
            {
                "id": identifier,
                "maximal_carrier_invariants_r_a_delta": expected_types[identifier],
                "maximal_carrier_isometry_type": carrier["isometry_type"],
                "index_over_original_fixed_lattice": carrier[
                    "index_over_original_lattice"
                ],
                "all_Coxeter_reflections_extend_integrally": True,
                "surface_realization": (
                    "very general K3 with this 2-elementary invariant lattice "
                    "and its classical nonsymplectic involution"
                ),
                "involution_fixes_full_Picard_pointwise": True,
                "realizes_original_ambient_action": False,
                "carrier_tier": "genuine pointwise-invariant divisor carrier",
                "effective_wall_claim": False,
                "imported_cone_equals_nef_claim": False,
            }
        )

    return {
        "schema": "coxeter.alexeev_parked_h9.k3_action_carriers.v1",
        "source": {
            "paths": {
                "parked": str(PARKED.relative_to(ROOT)),
                "denominator": str(DENOMINATOR.relative_to(ROOT)),
                "lineage": str(LINEAGE.relative_to(ROOT)),
            },
            "sha256": actual_hashes,
        },
        "strong_original_action_promotions": strong_records,
        "fifth_route_obstruction": fifth_obstruction,
        "pointwise_invariant_alternative_carriers": alternatives,
        "surjectivity_Phi_C": {
            "weak_lattice_carrier_input_count": 5,
            "original_action_genuine_invariant_divisor_carrier_count": 4,
            "alternative_genuine_pointwise_invariant_divisor_carrier_count": 5,
            "effective_or_nef_promotion_count": 0,
            "surface_free_Coxeter_diagram_forces_carrier_type": False,
            "additional_marking_needed": (
                "the ambient divisor lattice, integral involution, and its "
                "coinvariant embedding data; the facet diagram alone does not "
                "recover the K3 carrier"
            ),
            "precise_converse_statement": (
                "four of the five marked Alexeev sections lie in the image of "
                "the invariant-divisor carrier map via their original actions; "
                "all five lie in the image after forgetting the original action "
                "and using their maximal reflection-stable carrier"
            ),
        },
        "checks": {
            "all_four_strong_promotions_pass": len(strong_records) == 4,
            "all_four_glued_actions_are_Nikulin_involutions": all(
                record["torelli_and_geometry"]["genuine_K3_automorphism"]
                for record in strong_records
            ),
            "all_four_algebraic_coinvariants_are_E8_minus_2": all(
                record["final_Picard_lattice"]["coinvariant_lattice"][
                    "isometry_type"
                ]
                == "E8(-2)"
                for record in strong_records
            ),
            "fifth_parity_obstruction_passes": fifth_obstruction[
                "symplectic_completion_to_E8_minus_2_impossible"
            ],
            "all_five_alternative_carriers_pass": len(alternatives) == 5,
            "no_effective_or_nef_overclaim": True,
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = build_certificate()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise FileNotFoundError(args.output)
        if args.output.read_text(encoding="utf-8") != text:
            raise AssertionError("the archived K3-action carrier certificate changed")
        print(json.dumps({"status": "ok", "path": str(args.output)}))
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(json.dumps({"status": "wrote", "path": str(args.output)}))


if __name__ == "__main__":
    main()
