#!/usr/bin/env python3
"""Certify two genuine K3-automorphism lifts of the U+D8+D4 folds.

Two rows in Brandhorst--Hofmann's complete database have generic Picard
lattice ``U+D8(-1)+D4(-1)``.  Their restrictions to the Picard lattice are
conjugate to Roulleau wall-diagram classes 6 and 14.  Consequently the
finite-volume H11/16 and H10/14 fixed sections are not merely divisor-lattice
folds: after choosing the invariant nef chamber they are induced by genuine
purely non-symplectic K3 automorphisms of orders eight and sixteen.

The offline certificate pins the fixed and coinvariant lattices extracted
from the two database rows.  The optional ``--brandhorst-hofmann-zip`` flag
replays both complete 22-dimensional K3-lattice rows from the published
CC-BY archive.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence
import zipfile

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
)
from scripts.certify_a18_c06_h13_31_k3_order4_lift import matrix_digest
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    _parse_bh_row,
    _rational,
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
LINEAGE = ROOT / "data/enumerations/roulleau_ud8d4_fold_lineage.json"
OUTPUT = ROOT / "data/enumerations/roulleau_ud8d4_bh_k3_lifts.json"

EXPECTED_SOURCE_SHA256 = {
    "construction": "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
    "lineage": "654cbb85d6bf4956f513d9912463b43aa7679a166ee4fde67ccfff285f0a1538",
}

BH_ARTICLE_DOI = "10.1017/fms.2023.50"
BH_DATA_DOI = "10.5281/zenodo.10410278"
BH_DATA_VERSION = "2.1"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"


CASES: tuple[dict[str, Any], ...] = (
    {
        "name": "order8_type2",
        "bh_id": "0.8.1.3",
        "automorphism_order": 8,
        "member": "K3Groups/results/purely_ns/order08.txt",
        "member_sha256": "719eb569ec3f4d399a9948499a795e273160dfb54bb607e6e18eb31e188138d7",
        "member_row_count": 38,
        "row_one_based": 10,
        "row_sha256_with_newline": "b543f8d219c913c3d0834188aea0e62dfe724b562569a017e887b6d199d59bd9",
        "ns_gram_digest": "e07d9b83a4379fabaf6c562d03dd5153f5cf3065c39b1519c7523a09721ec01e",
        "ns_action_digest": "eb2ec4afe454ff935c628956e218ba4d039e7e6b1d9ecade5bb9df990c38140e",
        "archive_class": 6,
        "type_id": "R14_UD8D4_TYPE_2",
        "ns_action_order": 2,
        "fixed_gram": [
            [-2, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -2, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, -1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, -1],
            [0, 0, 0, 0, -1, -2, 0, 0, 0, 0, 0, -1],
            [0, 0, 0, 0, 0, 0, -2, -1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, -1, -2, 1, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, 1, 1, -2, -1, 0, -1],
            [0, 0, 0, 0, 0, 0, 1, 0, -1, -2, -1, -1],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -2, -1],
            [0, 0, 0, 0, -1, -1, 1, 1, -1, -1, -1, -2],
        ],
        "coinvariant_gram": [[-8, 4], [4, -4]],
        "coinvariant_action": [[-1, 0], [0, -1]],
        # If P*G_archive*P^t=G_BH, then P^{-1} is the coordinate isometry
        # from the archived coinvariant to the BH coinvariant.
        "coinvariant_basis_intertwiner": [[1, -1], [0, 1]],
        "bh_glue_pairs": [
            [[0] * 12, [0, 0]],
            [[0] * 10 + [2, 2], [2, 0]],
            [[0] * 10 + [2, 0], [2, 2]],
            [[0] * 11 + [2], [0, 2]],
        ],
        "fixed_discriminant_isometry": [[1, 1], [1, 0]],
        "expected": {
            "fixed_rank": 12,
            "fixed_determinant": -16,
            "fixed_smith": [*[1] * 10, 4, 4],
            "coinvariant_determinant": 16,
            "coinvariant_smith": [4, 4],
            "characteristic_polynomial": "(lambda - 1)**12*(lambda + 1)**2",
            "dimension": 11,
            "facets": 16,
            "finite_vertices": 112,
            "cusps": 9,
            "elliptic_edges": 750,
            "parabolic_audits": 18,
            "reflection_index": None,
        },
    },
    {
        "name": "order16_type4",
        "bh_id": "0.16.0.3",
        "automorphism_order": 16,
        "member": "K3Groups/results/purely_ns/order16.txt",
        "member_sha256": "2c8a58d846eb6b563a62b958bd665680790598b829d99723da463d9fa6e315ea",
        "member_row_count": 7,
        "row_one_based": 3,
        "row_sha256_with_newline": "2c78a6da28d3672826467dce6a7594d0adefec5f8d4fa2df869795751a8ffebb",
        "ns_gram_digest": "28410fa0543209b198c3aa27ba4a3a1d8cc7af264dfd8884f17b8b1e89ffb496",
        "ns_action_digest": "d3d9c4e9a4d979a48fd2b5cd49ac58aec84b38a61eb032fb403670d3e42dc7bc",
        "archive_class": 14,
        "type_id": "R14_UD8D4_TYPE_4",
        "ns_action_order": 4,
        "fixed_gram": [
            [-2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -2, -1, 0, 0, 0, 0, 0, -1],
            [0, 0, 0, -1, -2, 0, 0, 0, 0, 0, -1],
            [0, 0, 0, 0, 0, -2, -1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, -1, -2, 1, 0, 0, 1],
            [0, 0, 0, 0, 0, 1, 1, -2, -1, 0, -1],
            [0, 0, 0, 0, 0, 1, 0, -1, -2, -1, -1],
            [0, 0, 0, 0, 0, 1, 0, 0, -1, -2, -1],
            [0, 0, 0, -1, -1, 1, 1, -1, -1, -1, -2],
        ],
        "coinvariant_gram": [[-4, 0, 2], [0, -4, 2], [2, 2, -4]],
        "coinvariant_action": [[0, 0, 1], [-1, -1, -1], [0, 1, 0]],
        "coinvariant_basis_intertwiner": [[0, 1, 0], [1, 0, -1], [-1, 0, 0]],
        "bh_glue_pairs": [
            [[0] * 11, [0, 0, 0]],
            [[0] * 10 + [4], [0, 0, 4]],
            [[0] * 10 + [6], [0, 0, 6]],
            [[0] * 10 + [2], [0, 0, 2]],
        ],
        "fixed_discriminant_isometry": [[7]],
        "expected": {
            "fixed_rank": 11,
            "fixed_determinant": 8,
            "fixed_smith": [*[1] * 10, 8],
            "coinvariant_determinant": -32,
            "coinvariant_smith": [2, 2, 8],
            "characteristic_polynomial": "(lambda - 1)**11*(lambda + 1)*(lambda**2 + 1)",
            "dimension": 10,
            "facets": 14,
            "finite_vertices": 44,
            "cusps": 3,
            "elliptic_edges": 254,
            "parabolic_audits": 7,
            "reflection_index": 992,
        },
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if sp.Rational(value).q == 1 else str(value) for value in matrix.row(i)]
        for i in range(matrix.rows)
    ]


def exact_order(action: sp.Matrix, order: int) -> bool:
    if action**order != sp.eye(action.rows):
        return False
    return all(
        action**divisor != sp.eye(action.rows)
        for divisor in range(1, order)
        if order % divisor == 0
    )


def cyclotomic_evaluation(action: sp.Matrix, order: int) -> sp.Matrix:
    variable = sp.Symbol("x")
    polynomial = sp.Poly(sp.cyclotomic_poly(order, variable), variable)
    result = sp.zeros(action.rows)
    for (power,), coefficient in polynomial.terms():
        result += coefficient * action**power
    return result


def discriminant_coordinate(coordinates: sp.Matrix, gram: sp.Matrix) -> tuple[int, ...]:
    moduli, left_inverse, _ = _discriminant_data(gram)
    value = left_inverse.inv() * gram * coordinates
    if not all(sp.Rational(entry).q == 1 for entry in value):
        raise AssertionError("a glue representative did not land in the discriminant group")
    return tuple(int(value[i]) % moduli[i] for i in range(len(moduli)))


def glue_pairs(
    fixed_basis: sp.Matrix, coinvariant_basis: sp.Matrix, ambient_gram: sp.Matrix
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Return the graph of N/(F+C) in A_F+A_C in Smith coordinates."""

    split_basis = fixed_basis.col_join(coinvariant_basis)
    if split_basis.shape != ambient_gram.shape:
        raise AssertionError("the invariant and coinvariant ranks do not fill NS")
    inverse = split_basis.inv()
    generators = []
    for index in range(ambient_gram.rows):
        coordinates = sp.eye(ambient_gram.rows).row(index) * inverse
        residue = tuple(sp.Rational(value) % 1 for value in coordinates)
        if any(residue):
            generators.append(residue)
    zero = (sp.Rational(0),) * ambient_gram.rows
    subgroup = {zero}
    changed = True
    while changed:
        changed = False
        for left in tuple(subgroup):
            for right in generators:
                value = tuple((left[i] + right[i]) % 1 for i in range(len(left)))
                if value not in subgroup:
                    subgroup.add(value)
                    changed = True

    fixed_gram = fixed_basis * ambient_gram * fixed_basis.T
    coinvariant_gram = coinvariant_basis * ambient_gram * coinvariant_basis.T
    pairs = []
    for coordinates in subgroup:
        fixed = sp.Matrix(coordinates[: fixed_basis.rows])
        coinvariant = sp.Matrix(coordinates[fixed_basis.rows :])
        pairs.append(
            (
                discriminant_coordinate(fixed, fixed_gram),
                discriminant_coordinate(coinvariant, coinvariant_gram),
            )
        )
    return sorted(pairs)


def discriminant_map(
    coordinate_isometry: sp.Matrix, source_gram: sp.Matrix, target_gram: sp.Matrix
) -> sp.Matrix:
    """Induced map A_source -> A_target in Smith coordinates."""

    source_moduli, source_left_inverse, _ = _discriminant_data(source_gram)
    target_moduli, target_left_inverse, _ = _discriminant_data(target_gram)
    columns = []
    for index in range(len(source_moduli)):
        element = sp.zeros(len(source_moduli), 1)
        element[index] = 1
        source_coordinates = source_gram.inv() * source_left_inverse * element
        target_coordinates = coordinate_isometry.T * source_coordinates
        target_element = target_left_inverse.inv() * target_gram * target_coordinates
        if not all(sp.Rational(value).q == 1 for value in target_element):
            raise AssertionError("an integral isometry gave a nonintegral discriminant map")
        columns.append(
            [int(target_element[i]) % target_moduli[i] for i in range(len(target_moduli))]
        )
    return sp.Matrix(len(target_moduli), len(source_moduli), lambda i, j: columns[j][i])


def apply_discriminant_map(
    action: sp.Matrix, element: Sequence[int], target_moduli: Sequence[int]
) -> tuple[int, ...]:
    value = action * sp.Matrix(element)
    return tuple(int(value[i]) % target_moduli[i] for i in range(len(target_moduli)))


def quadratic_value(
    discriminant_gram: sp.Matrix, element: Sequence[int]
) -> sp.Rational:
    vector = sp.Matrix(element)
    return sp.Rational((vector.T * discriminant_gram * vector)[0]) % 2


def verify_fixed_discriminant_isometry(
    source_gram: sp.Matrix, target_gram: sp.Matrix, reduced_map: sp.Matrix
) -> dict[str, Any]:
    source_moduli, _, source_q = _discriminant_data(source_gram)
    target_moduli, _, target_q = _discriminant_data(target_gram)
    source_indices = [i for i, modulus in enumerate(source_moduli) if modulus > 1]
    target_indices = [i for i, modulus in enumerate(target_moduli) if modulus > 1]
    reduced_moduli = [source_moduli[i] for i in source_indices]
    if not (
        reduced_moduli == [target_moduli[i] for i in target_indices]
        and reduced_map.shape == (len(target_indices), len(source_indices))
    ):
        raise AssertionError("the fixed discriminant groups have different invariant factors")

    checked = 0
    images: set[tuple[int, ...]] = set()
    for reduced_source in itertools.product(*(range(modulus) for modulus in reduced_moduli)):
        reduced_target = tuple(
            int(value) % reduced_moduli[i]
            for i, value in enumerate(reduced_map * sp.Matrix(reduced_source))
        )
        images.add(reduced_target)
        source = [0] * len(source_moduli)
        target = [0] * len(target_moduli)
        for position, index in enumerate(source_indices):
            source[index] = reduced_source[position]
        for position, index in enumerate(target_indices):
            target[index] = reduced_target[position]
        if quadratic_value(source_q, source) != quadratic_value(target_q, target):
            raise AssertionError("the proposed fixed discriminant map does not preserve q")
        checked += 1
    if len(images) != math.prod(reduced_moduli):
        raise AssertionError("the proposed fixed discriminant map is not bijective")
    return {
        "source_invariant_factors": reduced_moduli,
        "target_invariant_factors": reduced_moduli,
        "matrix_on_nontrivial_smith_coordinates": matrix_payload(reduced_map),
        "elements_checked": checked,
        "bijective": True,
        "quadratic_form_preserved_mod_2Z": True,
    }


def apply_reduced_fixed_map(
    element: Sequence[int], source_gram: sp.Matrix, target_gram: sp.Matrix, action: sp.Matrix
) -> tuple[int, ...]:
    source_moduli, _, _ = _discriminant_data(source_gram)
    target_moduli, _, _ = _discriminant_data(target_gram)
    source_indices = [i for i, modulus in enumerate(source_moduli) if modulus > 1]
    target_indices = [i for i, modulus in enumerate(target_moduli) if modulus > 1]
    reduced = sp.Matrix([element[index] for index in source_indices])
    image = action * reduced
    result = [0] * len(target_moduli)
    for position, index in enumerate(target_indices):
        result[index] = int(image[position]) % target_moduli[index]
    return tuple(result)


def glue_overlattice(
    fixed_gram: sp.Matrix,
    coinvariant_gram: sp.Matrix,
    coinvariant_action: sp.Matrix,
    pairs: Iterable[tuple[Sequence[int], Sequence[int]]],
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Reconstruct the ambient lattice and action from F, C, and their glue."""

    fixed_moduli, fixed_left_inverse, _ = _discriminant_data(fixed_gram)
    coin_moduli, coin_left_inverse, _ = _discriminant_data(coinvariant_gram)
    split_gram = sp.diag(fixed_gram, coinvariant_gram)
    representatives = []
    for fixed_element, coin_element in pairs:
        if not (any(fixed_element) or any(coin_element)):
            continue
        fixed_coordinates = fixed_gram.inv() * fixed_left_inverse * sp.Matrix(fixed_element)
        coin_coordinates = coinvariant_gram.inv() * coin_left_inverse * sp.Matrix(coin_element)
        representatives.append(fixed_coordinates.col_join(coin_coordinates).T)
    generators = sp.eye(split_gram.rows)
    if representatives:
        generators = sp.Matrix.vstack(generators, *representatives)
    denominator = math.lcm(*(sp.Rational(value).q for value in generators))
    integral_generators = (denominator * generators).applyfunc(int)
    basis = hermite_normal_form(integral_generators.T).T / denominator
    if basis.shape != split_gram.shape:
        raise AssertionError("gluing did not produce a full-rank lattice")
    gram = sp.simplify(basis * split_gram * basis.T)
    split_action = sp.diag(sp.eye(fixed_gram.rows), coinvariant_action)
    action = sp.simplify(basis * split_action * basis.inv())
    if not (
        all(sp.Rational(value).q == 1 for value in gram)
        and all(sp.Rational(value).q == 1 for value in action)
    ):
        raise AssertionError("the pinned glue is not stable and integral")
    return basis, gram, action


def derive_bh_row(line: str, order: int) -> dict[str, Any]:
    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    source_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    source_action = sp.Matrix(22, 22, data["distinguished_generator"])
    k3_gram = sp.simplify(basis * source_gram * basis.T)
    k3_action = sp.simplify(basis * source_action * basis.inv())
    if not (
        k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in k3_gram)
        and all(k3_gram[i, i] % 2 == 0 for i in range(22))
        and all(sp.Rational(value).q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and exact_order(k3_action, order)
    ):
        raise AssertionError("the Brandhorst--Hofmann K3-lattice row failed")

    transcendental_basis = saturated_row_kernel(cyclotomic_evaluation(k3_action, order))
    ns_basis = saturated_row_kernel(k3_gram * transcendental_basis.T)
    transcendental_gram = sp.simplify(
        transcendental_basis * k3_gram * transcendental_basis.T
    )
    ns_gram = sp.simplify(ns_basis * k3_gram * ns_basis.T)
    ns_action = induced_row_action(ns_basis, k3_action)
    fixed_basis = saturated_row_kernel(ns_action - sp.eye(ns_action.rows))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * ns_gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, ns_action)
    return {
        "id": data["id"],
        "k3_gram": k3_gram,
        "k3_action": k3_action,
        "transcendental_gram": transcendental_gram,
        "ns_gram": ns_gram,
        "ns_action": ns_action,
        "fixed_gram": fixed_gram,
        "coinvariant_gram": coinvariant_gram,
        "coinvariant_action": coinvariant_action,
        "glue_pairs": glue_pairs(fixed_basis, coinvariant_basis, ns_gram),
    }


def replay_external_zip(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(payload) == BH_ZIP_SHA256
    ):
        raise AssertionError("the Brandhorst--Hofmann ZIP changed")
    replayed = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        for case in CASES:
            member = archive.read(case["member"])
            if sha256_bytes(member) != case["member_sha256"]:
                raise AssertionError("a pinned Brandhorst--Hofmann member changed")
            lines = member.decode("utf-8").splitlines()
            if len(lines) != case["member_row_count"]:
                raise AssertionError("a Brandhorst--Hofmann row count changed")
            line = lines[case["row_one_based"] - 1]
            if sha256_bytes((line + "\n").encode()) != case["row_sha256_with_newline"]:
                raise AssertionError("a pinned Brandhorst--Hofmann row changed")
            row = derive_bh_row(line, case["automorphism_order"])
            expected_pairs = [
                (tuple(fixed), tuple(coin)) for fixed, coin in case["bh_glue_pairs"]
            ]
            if not (
                row["id"] == case["bh_id"]
                and matrix_digest(row["ns_gram"]) == case["ns_gram_digest"]
                and matrix_digest(row["ns_action"]) == case["ns_action_digest"]
                and row["fixed_gram"] == sp.Matrix(case["fixed_gram"])
                and row["coinvariant_gram"] == sp.Matrix(case["coinvariant_gram"])
                and row["coinvariant_action"] == sp.Matrix(case["coinvariant_action"])
                and row["glue_pairs"] == sorted(expected_pairs)
                and row["transcendental_gram"].shape == (8, 8)
                and exact_inertia(row["transcendental_gram"]) == (2, 6, 0)
                and abs(int(row["transcendental_gram"].det())) == 16
            ):
                raise AssertionError("the exact Brandhorst--Hofmann row replay changed")
            replayed.append(case["bh_id"])
    return {"zip_sha256": BH_ZIP_SHA256, "rows_replayed": replayed}


def _archived_case(
    construction: dict[str, Any], class_index: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    case = next(
        row for row in construction["high_fixed_sections"] if row["class_index"] == class_index
    )
    geometry = next(
        row
        for row in construction["abstract_diagram_types"]
        if class_index in row["class_indices"]
    )
    return case, geometry


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    hashes = {"construction": sha256(CONSTRUCTION), "lineage": sha256(LINEAGE)}
    if hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned U+D8+D4 source certificate changed")
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (construction["all_checks_pass"] and lineage["all_checks_pass"]):
        raise AssertionError("an upstream U+D8+D4 certificate no longer passes")
    if external_zip is not None:
        replay_external_zip(external_zip)

    ns_gram = sp.Matrix(construction["source"]["basis_gram"])
    if not (
        ns_gram.det() == -16
        and exact_inertia(ns_gram) == (1, 13, 0)
        and smith_diagonal(ns_gram) == [*[1] * 10, *[2] * 4]
        and _q_histogram(ns_gram) == {"0": 6, "1": 10}
    ):
        raise AssertionError("the ambient (14,4,0) Picard lattice changed")

    results = []
    for specification in CASES:
        expected = specification["expected"]
        archived, geometry = _archived_case(construction, specification["archive_class"])
        lineage_type = next(
            row for row in lineage["types"] if row["type_id"] == specification["type_id"]
        )
        archived_action = sp.Matrix(
            archived["ns_action"]["representative_generator_action_matrices"][0]
        )
        fixed_basis = saturated_row_kernel(archived_action - sp.eye(14))
        coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
        archived_fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
        archived_coinvariant_gram = sp.simplify(
            coinvariant_basis * ns_gram * coinvariant_basis.T
        )
        archived_coinvariant_action = induced_row_action(
            coinvariant_basis, archived_action
        )
        archived_glue = glue_pairs(fixed_basis, coinvariant_basis, ns_gram)

        bh_fixed_gram = sp.Matrix(specification["fixed_gram"])
        bh_coinvariant_gram = sp.Matrix(specification["coinvariant_gram"])
        bh_coinvariant_action = sp.Matrix(specification["coinvariant_action"])
        bh_glue = sorted(
            (tuple(fixed), tuple(coin))
            for fixed, coin in specification["bh_glue_pairs"]
        )
        split_basis, reconstructed_ns_gram, reconstructed_ns_action = glue_overlattice(
            bh_fixed_gram, bh_coinvariant_gram, bh_coinvariant_action, bh_glue
        )

        characteristic = str(sp.factor(archived_action.charpoly().as_expr()))
        picard_checks = {
            "archived_action_order": exact_order(
                archived_action, specification["ns_action_order"]
            ),
            "archived_action_isometry": (
                archived_action * ns_gram * archived_action.T == ns_gram
            ),
            "characteristic_polynomial": (
                characteristic == expected["characteristic_polynomial"]
            ),
            "fixed_rank": archived_fixed_gram.rows == expected["fixed_rank"],
            "fixed_determinant": (
                int(archived_fixed_gram.det()) == expected["fixed_determinant"]
            ),
            "fixed_smith": smith_diagonal(archived_fixed_gram)
            == expected["fixed_smith"],
            "coinvariant_determinant": int(archived_coinvariant_gram.det())
            == expected["coinvariant_determinant"],
            "coinvariant_smith": smith_diagonal(archived_coinvariant_gram)
            == expected["coinvariant_smith"],
            "glue_order_four": len(archived_glue) == len(bh_glue) == 4,
            "glue_index_four": abs(split_basis.det()) == sp.Rational(1, 4),
            "reconstructed_determinant": reconstructed_ns_gram.det() == -16,
            "reconstructed_signature": exact_inertia(reconstructed_ns_gram)
            == (1, 13, 0),
            "reconstructed_even": all(
                reconstructed_ns_gram[i, i] % 2 == 0 for i in range(14)
            ),
            "reconstructed_smith": smith_diagonal(reconstructed_ns_gram)
            == [*[1] * 10, *[2] * 4],
            "reconstructed_q": _q_histogram(reconstructed_ns_gram)
            == _q_histogram(ns_gram),
            "reconstructed_action_order": exact_order(
                reconstructed_ns_action, specification["ns_action_order"]
            ),
        }
        if not all(picard_checks.values()):
            raise AssertionError(
                f"a pinned Picard action or reconstructed BH lattice changed: {picard_checks}"
            )

        # The coinvariant integral isometry and its discriminant action.
        coin_intertwiner = sp.Matrix(specification["coinvariant_basis_intertwiner"])
        if not (
            abs(int(coin_intertwiner.det())) == 1
            and coin_intertwiner
            * archived_coinvariant_gram
            * coin_intertwiner.T
            == bh_coinvariant_gram
            and coin_intertwiner * archived_coinvariant_action
            == bh_coinvariant_action * coin_intertwiner
        ):
            raise AssertionError("the coinvariant action intertwiner changed")
        coin_coordinate_isometry = coin_intertwiner.inv()
        coin_discriminant_map = discriminant_map(
            coin_coordinate_isometry,
            archived_coinvariant_gram,
            bh_coinvariant_gram,
        )

        fixed_discriminant_map = sp.Matrix(
            specification["fixed_discriminant_isometry"]
        )
        fixed_discriminant_certificate = verify_fixed_discriminant_isometry(
            archived_fixed_gram, bh_fixed_gram, fixed_discriminant_map
        )
        target_coin_moduli, _, _ = _discriminant_data(bh_coinvariant_gram)
        bh_glue_graph = {coin: fixed for fixed, coin in bh_glue}
        for archived_fixed, archived_coin in archived_glue:
            target_fixed = apply_reduced_fixed_map(
                archived_fixed,
                archived_fixed_gram,
                bh_fixed_gram,
                fixed_discriminant_map,
            )
            target_coin = apply_discriminant_map(
                coin_discriminant_map, archived_coin, target_coin_moduli
            )
            if bh_glue_graph.get(target_coin) != target_fixed:
                raise AssertionError("the fixed and coinvariant maps do not carry the glue graph")

        # The algebraic coinvariant has no (-2)-vectors: half its Gram is an
        # even negative-definite integral lattice, so all nonzero norms are
        # negative multiples of four.
        half_coinvariant = bh_coinvariant_gram / 2
        coinvariant_rootless = (
            exact_inertia(bh_coinvariant_gram)
            == (0, bh_coinvariant_gram.rows, 0)
            and all(sp.Rational(value).q == 1 for value in half_coinvariant)
            and all(half_coinvariant[i, i] % 2 == 0 for i in range(half_coinvariant.rows))
        )
        if not coinvariant_rootless:
            raise AssertionError("the BH algebraic coinvariant rootlessness changed")

        criterion = geometry["geometry"]["vinberg_volume_criterion"]
        if not (
            archived["geometry_type_id"] == geometry["type_id"] == specification["type_id"]
            and geometry["all_checks_pass"]
            and geometry["volume_type"] == "finite_volume_noncompact"
            and geometry["hyperbolic_dimension"] == expected["dimension"]
            and geometry["facet_count"] == expected["facets"]
            and geometry["geometry"]["finite_vertex_count"] == expected["finite_vertices"]
            and geometry["geometry"]["ideal_vertex_count"] == expected["cusps"]
            and criterion["decision"] == "finite_volume_noncompact"
            and criterion["elliptic_edge_count"] == expected["elliptic_edges"]
            and len(criterion["finite_volume_bad_edges"]) == 0
            and criterion["parabolic_neighborhood_audit_count"]
            == expected["parabolic_audits"]
            and lineage_type["reflection_subgroup_index"] == expected["reflection_index"]
        ):
            raise AssertionError("the Vinberg or reflection-lineage handoff changed")

        fixed_length = len(
            [value for value in smith_diagonal(archived_fixed_gram) if value > 1]
        )
        if not (
            exact_inertia(archived_fixed_gram)
            == exact_inertia(bh_fixed_gram)
            == (1, expected["fixed_rank"] - 1, 0)
            and _q_histogram(archived_fixed_gram) == _q_histogram(bh_fixed_gram)
            and expected["fixed_rank"] >= fixed_length + 2
        ):
            raise AssertionError("the Nikulin fixed-lattice gate changed")

        results.append(
            {
                "name": specification["name"],
                "theorem": (
                    f"A purely non-symplectic K3 automorphism of order "
                    f"{specification['automorphism_order']} realizes archived class "
                    f"{specification['archive_class']} and its finite-volume "
                    f"H{expected['dimension']}/{expected['facets']} fixed section"
                ),
                "brandhorst_hofmann_row": {
                    "id": specification["bh_id"],
                    "automorphism_order": specification["automorphism_order"],
                    "member": specification["member"],
                    "member_sha256": specification["member_sha256"],
                    "member_row_count": specification["member_row_count"],
                    "row_one_based": specification["row_one_based"],
                    "row_sha256_with_newline": specification[
                        "row_sha256_with_newline"
                    ],
                    "full_k3_action_replay_available": True,
                    "generic_transcendental_lattice": {
                        "rank": 8,
                        "signature": [2, 6, 0],
                        "absolute_determinant": 16,
                        "action_characteristic_polynomial": (
                            "(lambda**4 + 1)**2"
                            if specification["automorphism_order"] == 8
                            else "lambda**8 + 1"
                        ),
                    },
                    "generic_picard_lattice": {
                        "isometry_class": "U+D8(-1)+D4(-1)",
                        "nikulin_invariants": [14, 4, 0],
                        "rank_signature_determinant": [14, [1, 13, 0], -16],
                        "smith": [*[1] * 10, *[2] * 4],
                        "ns_gram_digest": specification["ns_gram_digest"],
                        "ns_action_digest": specification["ns_action_digest"],
                    },
                },
                "archived_divisor_action": {
                    "class_one_based": specification["archive_class"],
                    "type_id": specification["type_id"],
                    "group_structure": archived["structure"],
                    "ns_action_order": specification["ns_action_order"],
                    "characteristic_polynomial": characteristic,
                    "fixed_rank": expected["fixed_rank"],
                },
                "exact_action_conjugacy": {
                    "offline_bh_picard_reconstruction": {
                        "from_pinned_fixed_coinvariant_and_glue": True,
                        "rank_signature_determinant": [14, [1, 13, 0], -16],
                        "smith": [*[1] * 10, *[2] * 4],
                        "discriminant_q_histogram_mod_2Z": _q_histogram(
                            reconstructed_ns_gram
                        ),
                        "action_order": specification["ns_action_order"],
                        "action_characteristic_polynomial": str(
                            sp.factor(reconstructed_ns_action.charpoly().as_expr())
                        ),
                        "gram_digest": matrix_digest(reconstructed_ns_gram),
                        "action_digest": matrix_digest(reconstructed_ns_action),
                    },
                    "fixed_lattices": {
                        "signature": [1, expected["fixed_rank"] - 1, 0],
                        "determinant": expected["fixed_determinant"],
                        "smith": expected["fixed_smith"],
                        "archived_q_histogram_mod_2Z": _q_histogram(
                            archived_fixed_gram
                        ),
                        "bh_q_histogram_mod_2Z": _q_histogram(bh_fixed_gram),
                        "discriminant_isometry": fixed_discriminant_certificate,
                        "nikulin_lifting": {
                            "reference": "Nikulin, Integral symmetric bilinear forms, Theorem 1.14.2",
                            "doi": "10.1070/IM1980v014n01ABEH001060",
                            "indefinite_even": True,
                            "rank": expected["fixed_rank"],
                            "discriminant_length": fixed_length,
                            "rank_at_least_length_plus_two": True,
                            "conclusion": (
                                "the prescribed discriminant-form isometry lifts "
                                "to an integral isometry of the fixed lattices"
                            ),
                        },
                    },
                    "coinvariant_lattices": {
                        "archived_gram": matrix_payload(archived_coinvariant_gram),
                        "bh_gram": matrix_payload(bh_coinvariant_gram),
                        "smith": expected["coinvariant_smith"],
                        "basis_intertwiner_P_with_PGarchivePt_equal_Gbh": matrix_payload(
                            coin_intertwiner
                        ),
                        "action_intertwined": True,
                        "induced_discriminant_map": matrix_payload(
                            coin_discriminant_map
                        ),
                        "negative_definite": True,
                        "rootless": True,
                        "rootless_reason": (
                            "one half of the Gram is integral, even and negative "
                            "definite, so every nonzero norm is at most -4"
                        ),
                    },
                    "primitive_glue": {
                        "index_of_fixed_plus_coinvariant": 4,
                        "archived_pairs_in_smith_coordinates": [
                            [list(fixed), list(coin)] for fixed, coin in archived_glue
                        ],
                        "bh_pairs_in_smith_coordinates": [
                            [list(fixed), list(coin)] for fixed, coin in bh_glue
                        ],
                        "all_four_pairs_carried_exactly": True,
                    },
                    "conclusion": (
                        "the fixed and coinvariant isometries carry the primitive "
                        "glue graph and intertwine the actions, hence extend to a "
                        "conjugacy of the full Picard-lattice actions"
                    ),
                },
                "kahler_torelli_handoff": {
                    "database_row_is_a_genuine_k3_automorphism_class": True,
                    "algebraic_coinvariant_is_rootless": True,
                    "archived_complete_nef_chamber_is_action_stable": True,
                    "refined_period_map_choice": (
                        "use the BH period and choose the marked Kahler chamber to "
                        "be the archived action-stable Roulleau chamber"
                    ),
                    "strong_torelli_conclusion": (
                        "the transported Hodge isometry preserving that chamber is "
                        "induced by a K3 automorphism"
                    ),
                    "global_order_vs_picard_order": [
                        specification["automorphism_order"],
                        specification["ns_action_order"],
                    ],
                },
                "hyperbolic_polytope": {
                    "dimension": expected["dimension"],
                    "facets": expected["facets"],
                    "classification": "finite_volume_noncompact",
                    "finite_vertices": expected["finite_vertices"],
                    "cusps": expected["cusps"],
                    "elliptic_edges": expected["elliptic_edges"],
                    "bad_edges": 0,
                    "parabolic_neighborhood_audits": expected[
                        "parabolic_audits"
                    ],
                    "containing_full_reflection_group_index": expected[
                        "reflection_index"
                    ],
                },
                "scope": {
                    "new_positive_surface_realization_in_this_project": True,
                    "new_abstract_coxeter_diagram_claimed": False,
                    "worldwide_surface_realization_priority_claimed": False,
                    "general_roulleau_surface_action_claimed": False,
                    "special_k3_in_bh_deformation_class_used": True,
                },
            }
        )

    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau_ud8d4.bh_k3_lifts.v1",
        "all_checks_pass": True,
        "theorem": (
            "Two previously diagram-only Roulleau U+D8+D4 folds admit genuine "
            "K3-automorphism realizations: TYPE_2 H11/16 via BH 0.8.1.3 and "
            "TYPE_4 H10/14 via BH 0.16.0.3"
        ),
        "source_sha256": hashes,
        "brandhorst_hofmann": {
            "article": "Finite subgroups of automorphisms of K3 surfaces",
            "article_doi": BH_ARTICLE_DOI,
            "database_doi": BH_DATA_DOI,
            "database_version": BH_DATA_VERSION,
            "license": "CC-BY-4.0",
            "zip_md5": BH_ZIP_MD5,
            "zip_sha256": BH_ZIP_SHA256,
            "external_replay_command": (
                "python3 scripts/certify_roulleau_ud8d4_bh_k3_lifts.py "
                "--brandhorst-hofmann-zip /path/to/K3Groups.zip --check-archive"
            ),
        },
        "cases": results,
        "summary": {
            "genuine_k3_automorphism_lifts": 2,
            "dimensions": [11, 10],
            "finite_volume_noncompact": 2,
            "compact": 0,
            "new_abstract_diagram_claims": 0,
        },
    }
    canonical = json.dumps(certificate, sort_keys=True, separators=(",", ":"))
    certificate["payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_zip)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise AssertionError("the archived U+D8+D4 K3-lift certificate changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
