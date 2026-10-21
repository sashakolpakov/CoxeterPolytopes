#!/usr/bin/env python3
"""Certify two genuine K3-automorphism realizations of the U+E8+D4 folds.

Brandhorst--Hofmann rows 0.6.3.1 and 0.16.0.1 have generic Picard lattice
``U+E8(-1)+D4(-1)``.  Their Picard actions are respectively the order-three
and order-two outer diagram symmetries of D4.  They therefore realize the
already certified Roulleau H11/13 and H12/14 finite-volume Coxeter sections
on genuine K3 surfaces.

The default certificate is offline: it verifies frozen exact Gram/action
matrices derived from the published rows and the pinned Roulleau/Vinberg
artifact.  ``--brandhorst-hofmann-zip`` additionally replays both complete
22-dimensional rows from the published CC-BY database archive.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import itertools
import json
from pathlib import Path
from typing import Any, Sequence
import zipfile

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    _parse_bh_row,
    _rational,
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json"
OUTPUT = (
    ROOT
    / "data/enumerations/roulleau_ue8d4_generic_ns_k3_automorphisms.json"
)
EXPECTED_SOURCE_SHA256 = (
    "de41f13f9b47fa2c1e82bbd9255da74667694658e5acac2461054bb5d2996f04"
)

BH_ARTICLE_DOI = "10.1017/fms.2023.50"
BH_DATA_DOI = "10.5281/zenodo.10410278"
BH_DATA_VERSION = "2.1"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_README_MEMBER = "K3Groups/readme.txt"
BH_README_SHA256 = "588cfbaed63ade81f8c34ef256240ca26ece3f69074417b5f301f31c418bb9e9"


D4_GRAM = [
    [-2, 1, 0, 0],
    [1, -2, 1, 1],
    [0, 1, -2, 0],
    [0, 1, 0, -2],
]


CASES: tuple[dict[str, Any], ...] = (
    {
        "name": "order6_triality",
        "bh_id": "0.6.3.1",
        "automorphism_order": 6,
        "member": "K3Groups/results/purely_ns/order06.txt",
        "member_sha256": "45074265c680c21efafe549b03e34cfde10467f10e2a6cd73ef3a9eb19401f63",
        "member_row_count": 150,
        "row_one_based": 57,
        "row_sha256_with_newline": "aacb4c5b6aa9e9f622921803a076c4b0baa0ff40a19d3f4cfd377325bd56b344",
        "type_id": "R14_UE8D4_C3_S3",
        "picard_action_order": 3,
        "ns_gram_digest": "1e97b94f50ac272a897006ffd419af6d2bb2d4b8909e8364041e410f3431782e",
        "ns_action_digest": "2adeb56426c8bddffbfccf5c22a3413d92c6f37fe229e9694818ee280f548cf9",
        "transcendental_gram_digest": "c4c0fdd31752d2907ed9a66b370f7b01f1fc19eced6caaa7d82675c32a657f3b",
        "ns_gram": [
            [-2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0],
            [-2, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, -2, -1, 1, 1, 1, -1, 1, 1, 0, 0, 0],
            [0, 0, 0, -1, -2, 1, 1, 1, -1, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, -2, 0, -1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, -2, -1, 1, -1, -1, 0, 0, 0],
            [0, 0, 0, 1, 1, -1, -1, -2, 1, -1, -1, 0, 0, 0],
            [0, 0, 0, -1, -1, 0, 1, 1, -2, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, -1, -1, 1, -2, -1, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, -1, -1, 1, -1, -2, 0, 0, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
        "ns_action": [
            [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
            [-3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        "transcendental_gram": [
            [-12, 10, 2, 2, 2, 0, 2, 4],
            [10, -12, -2, -2, 0, 2, -2, -4],
            [2, -2, 0, 0, 0, 0, -1, -1],
            [2, -2, 0, 0, 0, 0, 0, -1],
            [2, 0, 0, 0, -2, -1, 0, 0],
            [0, 2, 0, 0, -1, -2, 0, 0],
            [2, -2, -1, 0, 0, 0, -2, -1],
            [4, -4, -1, -1, 0, 0, -1, -2],
        ],
        "complement_indices": [2, 13, 3, 4, 5, 6, 7, 8, 9, 10],
        "d4_indices": [0, 1, 11, 12],
        "d4_simple_basis_in_d4_coordinates": [
            [-2, 1, 1, 1],
            [0, 0, -1, 0],
            [1, -1, 0, 0],
            [1, 0, 0, 0],
        ],
        "diagram_action": [
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [1, 0, 0, 0],
        ],
        "fixed_rank": 12,
        "fixed_determinant": -3,
        "fixed_smith": [*[1] * 11, 3],
        "coinvariant_gram": [[-12, 6], [6, -4]],
        "coinvariant_action": [[-2, -3], [1, 1]],
        "coinvariant_canonical_change": [[-1, -2], [0, 1]],
        "coinvariant_class": "A2(-2)",
        "full_action_characteristic": "(lambda - 1)**12*(lambda**2 - lambda + 1)**4*(lambda**2 + lambda + 1)",
        "transcendental_characteristic": "(lambda**2 - lambda + 1)**4",
        "picard_characteristic": "(lambda - 1)**12*(lambda**2 + lambda + 1)",
        "dimension": 11,
        "facets": 13,
        "elliptic_edges": 162,
        "finite_vertices": 27,
        "cusps": 1,
        "parabolic_audits": 2,
        "vinberg_sha256": "52818c7f1c2e8813d49b0b7f9cdc3af956d823b2fddd2bb302b177e59632d80f",
    },
    {
        "name": "order16_transposition",
        "bh_id": "0.16.0.1",
        "automorphism_order": 16,
        "member": "K3Groups/results/purely_ns/order16.txt",
        "member_sha256": "2c8a58d846eb6b563a62b958bd665680790598b829d99723da463d9fa6e315ea",
        "member_row_count": 7,
        "row_one_based": 1,
        "row_sha256_with_newline": "22f66c944bc0bdaa0cfdfc0ba1e3f935c2ba6593eff91cec36fc57bd41199c18",
        "type_id": "R14_UE8D4_C2",
        "picard_action_order": 2,
        "ns_gram_digest": "3a814211e13c52c671cfa34867ee8d3fe3c209246788a223df671c6ea6dd0746",
        "ns_action_digest": "9258ac7cf8f7e40d11e1368cd7f288b065aaab17b84a2961f6d795b56a824754",
        "transcendental_gram_digest": "6a293432179a50a1c8e66c5fb23ade0ba63bd1d8dea6ae8129c2e3750694b4d5",
        "ns_gram": [
            [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, -2, -1, 1, 1, 1, -1, 1, 1, 0, 0, 0, 0],
            [0, 0, -1, -2, 1, 1, 1, -1, 1, 1, 0, 0, 0, 0],
            [0, 0, 1, 1, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, -2, -1, 1, -1, -1, 0, 0, 0, 0],
            [0, 0, 1, 1, -1, -1, -2, 1, -1, -1, 0, 0, 0, 0],
            [0, 0, -1, -1, 0, 1, 1, -2, 1, 1, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, -1, -1, 1, -2, -1, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, -1, -1, 1, -1, -2, 0, 0, 0, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, -1, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, 0, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -2, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
        "ns_action": [
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        "transcendental_gram": [
            [-28, -14, 0, 4, 4, 10, 6, 10],
            [-14, -28, 4, 4, 10, 6, 10, 4],
            [0, 4, -4, -2, -3, -1, 0, 1],
            [4, 4, -2, -4, -2, -3, -1, 0],
            [4, 10, -3, -2, -4, -2, -3, -1],
            [10, 6, -1, -3, -2, -4, -2, -3],
            [6, 10, 0, -1, -3, -2, -4, -2],
            [10, 4, 1, 0, -1, -3, -2, -4],
        ],
        "complement_indices": [1, 13, 2, 3, 4, 5, 6, 7, 8, 9],
        "d4_indices": [0, 10, 11, 12],
        "d4_simple_basis_in_d4_coordinates": [
            [0, 1, -1, 0],
            [0, -1, 0, 0],
            [-1, 0, 1, 1],
            [1, 0, 0, 0],
        ],
        "diagram_action": [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
        ],
        "fixed_rank": 13,
        "fixed_determinant": 4,
        "fixed_smith": [*[1] * 12, 4],
        "coinvariant_gram": [[-4]],
        "coinvariant_action": [[-1]],
        "coinvariant_canonical_change": [[1]],
        "coinvariant_class": "<-4>",
        "full_action_characteristic": "(lambda - 1)**13*(lambda + 1)*(lambda**8 + 1)",
        "transcendental_characteristic": "lambda**8 + 1",
        "picard_characteristic": "(lambda - 1)**13*(lambda + 1)",
        "dimension": 12,
        "facets": 14,
        "elliptic_edges": 234,
        "finite_vertices": 35,
        "cusps": 2,
        "parabolic_audits": 3,
        "vinberg_sha256": "9a50c14e82ac368380003ce60ac6fff706fec5a02125d71de69d32eb8c1af9c8",
    },
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(x) if sp.Rational(x).q == 1 else str(x) for x in matrix.row(i)]
        for i in range(matrix.rows)
    ]


def matrix_digest(matrix: sp.Matrix) -> str:
    canonical = json.dumps(
        [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)],
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode()).hexdigest()


def exact_order(action: sp.Matrix, order: int) -> bool:
    return action**order == sp.eye(action.rows) and all(
        action**divisor != sp.eye(action.rows)
        for divisor in range(1, order)
        if order % divisor == 0
    )


def induced_row_action(basis: sp.Matrix, ambient_action: sp.Matrix) -> sp.Matrix:
    pivots = basis.rref()[1]
    square = basis[:, list(pivots)]
    action = (basis * ambient_action)[:, list(pivots)] * square.inv()
    if not (
        action * basis == basis * ambient_action
        and all(sp.Rational(value).q == 1 for value in action)
    ):
        raise AssertionError("an invariant sublattice action failed")
    return action


def cyclotomic_evaluation(action: sp.Matrix, order: int) -> sp.Matrix:
    x = sp.Symbol("x")
    polynomial = sp.Poly(sp.cyclotomic_poly(order, x), x)
    result = sp.zeros(action.rows)
    for (power,), coefficient in polynomial.terms():
        result += coefficient * action**power
    return result


def _gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    work = [[int(value) & 1 for value in row] for row in rows]
    rank = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and work[index][column]:
                work[index] = [a ^ b for a, b in zip(work[index], work[rank])]
        rank += 1
    return rank


def two_elementary_discriminant_basis(gram: sp.Matrix) -> sp.Matrix:
    """Rows v modulo two for which v/2 lies in the dual lattice."""

    work = [[int(gram[i, j]) & 1 for j in range(gram.cols)] for i in range(gram.rows)]
    pivots: list[int] = []
    rank = 0
    for column in range(gram.cols):
        pivot = next(
            (index for index in range(rank, gram.rows) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(gram.rows):
            if index != rank and work[index][column]:
                work[index] = [a ^ b for a, b in zip(work[index], work[rank])]
        pivots.append(column)
        rank += 1
    free = [column for column in range(gram.cols) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [0] * gram.cols
        vector[free_column] = 1
        for row, pivot in reversed(list(enumerate(pivots))):
            vector[pivot] = sum(
                work[row][column] * vector[column] for column in free
            ) & 1
        basis.append(vector)
    return sp.Matrix(basis)


def _gf2_coordinates(vector: Sequence[int], basis: sp.Matrix) -> tuple[int, ...]:
    target = tuple(int(value) & 1 for value in vector)
    for coefficients in itertools.product(range(2), repeat=basis.rows):
        image = sp.Matrix([coefficients]) * basis
        if tuple(int(value) & 1 for value in image) == target:
            return tuple(coefficients)
    raise AssertionError("a discriminant vector left its GF(2) basis")


def _gf2_cycle_lengths(action: sp.Matrix) -> list[int]:
    unseen = set(itertools.product(range(2), repeat=action.rows))
    unseen.discard((0,) * action.rows)
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        orbit = []
        while current not in orbit:
            orbit.append(current)
            image = sp.Matrix([current]) * action
            current = tuple(int(value) & 1 for value in image)
        unseen.difference_update(orbit)
        lengths.append(len(orbit))
    return sorted(lengths)


def discriminant_audit(gram: sp.Matrix, action: sp.Matrix) -> dict[str, Any]:
    basis = two_elementary_discriminant_basis(gram)
    if not (basis.rows == 2 and _gf2_rank(basis.tolist()) == 2):
        raise AssertionError("the ambient discriminant is not (Z/2)^2")
    action_rows = []
    for row in range(basis.rows):
        action_rows.append(_gf2_coordinates(basis.row(row) * action, basis))
    discr_action = sp.Matrix(action_rows)
    table = []
    histogram: Counter[int] = Counter()
    for coordinates in itertools.product(range(2), repeat=2):
        vector = sp.Matrix([coordinates]) * basis
        numerator = int((vector * gram * vector.T)[0])
        if numerator % 4:
            raise AssertionError("a half-vector has a nonintegral q value")
        q_value = (numerator // 4) % 2
        image = sp.Matrix([coordinates]) * discr_action
        image_coordinates = [int(value) & 1 for value in image]
        histogram[q_value] += 1
        table.append(
            {
                "coordinates": list(coordinates),
                "q_mod_2Z": q_value,
                "image_coordinates": image_coordinates,
            }
        )
    return {
        "group": "(Z/2)^2",
        "half_vector_basis_rows": matrix_payload(basis),
        "q_histogram_mod_2Z": {str(key): value for key, value in sorted(histogram.items())},
        "action_matrix_mod_2": matrix_payload(discr_action),
        "nonzero_action_cycle_lengths": _gf2_cycle_lengths(discr_action),
        "table": table,
        "_action": discr_action,
    }


def gf2_conjugator(source: sp.Matrix, target: sp.Matrix) -> sp.Matrix:
    for values in itertools.product(range(2), repeat=4):
        candidate = sp.Matrix(2, 2, values)
        if int(candidate.det()) % 2 and (source * candidate - candidate * target).applyfunc(
            lambda value: int(value) & 1
        ) == sp.zeros(2):
            return candidate
    raise AssertionError("the two discriminant actions are not conjugate")


def derive_bh_row(line: str, order: int) -> dict[str, Any]:
    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    ambient_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    ambient_action = sp.Matrix(22, 22, data["distinguished_generator"])
    k3_gram = sp.simplify(basis * ambient_gram * basis.T)
    k3_action = sp.simplify(basis * ambient_action * basis.inv())
    if not (
        data["transcendental_value"] == order
        and data["connected_components"] == 1
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and all(sp.Rational(value).q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and exact_order(k3_action, order)
    ):
        raise AssertionError("a Brandhorst--Hofmann K3 row failed")

    transcendental_basis = saturated_row_kernel(
        cyclotomic_evaluation(k3_action, order)
    )
    ns_basis = saturated_row_kernel(k3_gram * transcendental_basis.T)
    transcendental_gram = sp.simplify(
        transcendental_basis * k3_gram * transcendental_basis.T
    )
    ns_gram = sp.simplify(ns_basis * k3_gram * ns_basis.T)
    ns_action = induced_row_action(ns_basis, k3_action)
    return {
        "id": data["id"],
        "k3_action": k3_action,
        "transcendental_gram": transcendental_gram,
        "transcendental_action": induced_row_action(
            transcendental_basis, k3_action
        ),
        "ns_gram": ns_gram,
        "ns_action": ns_action,
    }


def replay_external_zip(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(payload) == BH_ZIP_SHA256
    ):
        raise AssertionError("the Brandhorst--Hofmann ZIP changed")
    rows_replayed = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        readme = archive.read(BH_README_MEMBER)
        if not (
            sha256_bytes(readme) == BH_README_SHA256
            and b"each line represents a moduli space parametrizing pairs (X,G)" in readme
            and b"effective character" in readme
        ):
            raise AssertionError("the Brandhorst--Hofmann README changed")
        for case in CASES:
            member = archive.read(case["member"])
            if sha256_bytes(member) != case["member_sha256"]:
                raise AssertionError("a Brandhorst--Hofmann member changed")
            lines = member.decode("utf-8").splitlines()
            if len(lines) != case["member_row_count"]:
                raise AssertionError("a Brandhorst--Hofmann row count changed")
            line = lines[case["row_one_based"] - 1]
            if sha256_bytes((line + "\n").encode()) != case["row_sha256_with_newline"]:
                raise AssertionError("a Brandhorst--Hofmann target row changed")
            row = derive_bh_row(line, case["automorphism_order"])
            if not (
                row["id"] == case["bh_id"]
                and row["ns_gram"] == sp.Matrix(case["ns_gram"])
                and row["ns_action"] == sp.Matrix(case["ns_action"])
                and row["transcendental_gram"]
                == sp.Matrix(case["transcendental_gram"])
                and str(sp.factor(row["k3_action"].charpoly().as_expr()))
                == case["full_action_characteristic"]
                and str(
                    sp.factor(row["transcendental_action"].charpoly().as_expr())
                )
                == case["transcendental_characteristic"]
            ):
                raise AssertionError("an exact Brandhorst--Hofmann replay changed")
            rows_replayed.append(case["bh_id"])
    return {"zip_sha256": BH_ZIP_SHA256, "rows_replayed": rows_replayed}


def _coordinate_rows(indices: Sequence[int], rank: int = 14) -> sp.Matrix:
    identity = sp.eye(rank)
    return identity[list(indices), :]


def _source_d4_decomposition(
    source: dict[str, Any], archived_action: sp.Matrix
) -> dict[str, Any]:
    gram = sp.Matrix(source["source"]["basis_gram"])
    walls = sp.Matrix(source["source"]["wall_coordinate_rows"])
    # Outer A13, central A12, outer A14, outer A15.
    d4_basis = walls[[12, 11, 13, 14], :]
    complement = saturated_row_kernel(gram * d4_basis.T)
    block_basis = complement.col_join(d4_basis)
    block_gram = sp.simplify(block_basis * gram * block_basis.T)
    block_action = induced_row_action(block_basis, archived_action)
    if not (
        abs(int(block_basis.det())) == 1
        and block_gram[:10, 10:] == sp.zeros(10, 4)
        and block_gram[10:, :10] == sp.zeros(4, 10)
        and block_gram[10:, 10:] == sp.Matrix(D4_GRAM)
        and exact_inertia(block_gram[:10, :10]) == (1, 9, 0)
        and block_gram[:10, :10].det() == -1
        and all(block_gram[index, index] % 2 == 0 for index in range(10))
        and block_action[:10, :10] == sp.eye(10)
        and block_action[:10, 10:] == sp.zeros(10, 4)
        and block_action[10:, :10] == sp.zeros(4, 10)
    ):
        raise AssertionError("the archived II(1,9)+D4 decomposition changed")
    return {
        "complement_gram": block_gram[:10, :10],
        "d4_action": block_action[10:, 10:],
        "block_basis_determinant": int(block_basis.det()),
    }


def _bh_d4_decomposition(case: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(case["ns_gram"])
    action = sp.Matrix(case["ns_action"])
    complement = _coordinate_rows(case["complement_indices"])
    d4_coordinates = _coordinate_rows(case["d4_indices"])
    d4_basis = sp.Matrix(case["d4_simple_basis_in_d4_coordinates"]) * d4_coordinates
    block_basis = complement.col_join(d4_basis)
    block_gram = sp.simplify(block_basis * gram * block_basis.T)
    block_action = induced_row_action(block_basis, action)
    complement_gram = block_gram[:10, :10]
    e8_gram = complement_gram[2:, 2:]
    if not (
        abs(int(block_basis.det())) == 1
        and complement_gram[:2, :2] == sp.Matrix([[0, 1], [1, 0]])
        and complement_gram[:2, 2:] == sp.zeros(2, 8)
        and complement_gram[2:, :2] == sp.zeros(8, 2)
        and exact_inertia(e8_gram) == (0, 8, 0)
        and e8_gram.det() == 1
        and all(e8_gram[index, index] % 2 == 0 for index in range(8))
        and block_gram[10:, 10:] == sp.Matrix(D4_GRAM)
        and block_gram[:10, 10:] == sp.zeros(10, 4)
        and block_action[:10, :10] == sp.eye(10)
        and block_action[:10, 10:] == sp.zeros(10, 4)
        and block_action[10:, :10] == sp.zeros(4, 10)
        and block_action[10:, 10:] == sp.Matrix(case["diagram_action"])
    ):
        raise AssertionError("a BH U+E8+D4 block decomposition changed")
    return {
        "complement_gram": complement_gram,
        "d4_action": block_action[10:, 10:],
        "block_basis_determinant": int(block_basis.det()),
        "e8_gram_digest": matrix_digest(e8_gram),
    }


def _archive_case(source: dict[str, Any], type_id: str) -> dict[str, Any]:
    return next(row for row in source["types"] if row["type_id"] == type_id)


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise AssertionError("the pinned Roulleau U+E8+D4 certificate changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if not (
        source["all_checks_pass"]
        and source["weighted_diagram_group"]["structure"] == "S3"
        and source["weighted_diagram_group"]["order"] == 6
    ):
        raise AssertionError("the archived Roulleau chamber is no longer green")
    if external_zip is not None:
        replay_external_zip(external_zip)

    archived_ns_gram = sp.Matrix(source["source"]["basis_gram"])
    if not (
        archived_ns_gram.det() == -4
        and exact_inertia(archived_ns_gram) == (1, 13, 0)
        and smith_diagonal(archived_ns_gram) == [*[1] * 12, 2, 2]
        and all(archived_ns_gram[i, i] % 2 == 0 for i in range(14))
    ):
        raise AssertionError("the archived Picard lattice changed")

    results = []
    source_discriminant_actions = []
    for case in CASES:
        archived = _archive_case(source, case["type_id"])
        geometry = archived["geometry"]
        criterion = geometry["vinberg_volume_criterion"]
        archived_action = sp.Matrix(
            archived["ns_action"]["generator_action_matrices"][0]
        )
        ns_gram = sp.Matrix(case["ns_gram"])
        ns_action = sp.Matrix(case["ns_action"])
        transcendental_gram = sp.Matrix(case["transcendental_gram"])

        fixed_basis = saturated_row_kernel(ns_action - sp.eye(14))
        fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
        coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
        coinvariant_gram = sp.simplify(
            coinvariant_basis * ns_gram * coinvariant_basis.T
        )
        coinvariant_action = induced_row_action(coinvariant_basis, ns_action)
        canonical_change = sp.Matrix(case["coinvariant_canonical_change"])
        canonical_coinvariant = sp.simplify(
            canonical_change * coinvariant_gram * canonical_change.T
        )

        expected_canonical = (
            sp.Matrix([[-4, 2], [2, -4]])
            if case["coinvariant_class"] == "A2(-2)"
            else sp.Matrix([[-4]])
        )
        algebraic_checks = {
            "ns_signature": exact_inertia(ns_gram) == (1, 13, 0),
            "ns_determinant": ns_gram.det() == -4,
            "ns_smith": smith_diagonal(ns_gram) == [*[1] * 12, 2, 2],
            "ns_even": all(ns_gram[i, i] % 2 == 0 for i in range(14)),
            "action_isometry": ns_action * ns_gram * ns_action.T == ns_gram,
            "action_order": exact_order(ns_action, case["picard_action_order"]),
            "action_characteristic": str(sp.factor(ns_action.charpoly().as_expr()))
            == case["picard_characteristic"],
            "fixed_rank": fixed_gram.rows == case["fixed_rank"],
            "fixed_signature": exact_inertia(fixed_gram)
            == (1, case["fixed_rank"] - 1, 0),
            "fixed_determinant": fixed_gram.det() == case["fixed_determinant"],
            "fixed_smith": smith_diagonal(fixed_gram) == case["fixed_smith"],
            "coinvariant_gram": coinvariant_gram
            == sp.Matrix(case["coinvariant_gram"]),
            "coinvariant_action": coinvariant_action
            == sp.Matrix(case["coinvariant_action"]),
            "coinvariant_canonical": canonical_coinvariant == expected_canonical,
            "coinvariant_rootless": (
                exact_inertia(canonical_coinvariant)
                == (0, canonical_coinvariant.rows, 0)
                and all(value % 2 == 0 for value in canonical_coinvariant)
            ),
            "transcendental_signature": exact_inertia(transcendental_gram)
            == (2, 6, 0),
            "transcendental_determinant": transcendental_gram.det() == 4,
            "transcendental_smith": smith_diagonal(transcendental_gram)
            == [*[1] * 6, 2, 2],
            "frozen_matrix_digests": (
                matrix_digest(ns_gram) == case["ns_gram_digest"]
                and matrix_digest(ns_action) == case["ns_action_digest"]
                and matrix_digest(transcendental_gram)
                == case["transcendental_gram_digest"]
            ),
        }
        if not all(algebraic_checks.values()):
            raise AssertionError(
                f"a frozen BH generic-NS calculation changed: {algebraic_checks}"
            )

        source_block = _source_d4_decomposition(source, archived_action)
        bh_block = _bh_d4_decomposition(case)
        if not (
            source_block["d4_action"] == bh_block["d4_action"]
            == sp.Matrix(case["diagram_action"])
        ):
            raise AssertionError("the exact D4 diagram actions do not match")

        source_discriminant = discriminant_audit(
            archived_ns_gram, archived_action
        )
        bh_discriminant = discriminant_audit(ns_gram, ns_action)
        conjugator = gf2_conjugator(
            source_discriminant["_action"], bh_discriminant["_action"]
        )
        source_discriminant_actions.append(source_discriminant["_action"])
        for audit in (source_discriminant, bh_discriminant):
            del audit["_action"]
        if not (
            source_discriminant["q_histogram_mod_2Z"]
            == bh_discriminant["q_histogram_mod_2Z"]
            == {"0": 1, "1": 3}
            and source_discriminant["nonzero_action_cycle_lengths"]
            == bh_discriminant["nonzero_action_cycle_lengths"]
        ):
            raise AssertionError("the discriminant q/action class did not match")

        geometry_checks = {
            "upstream_type_green": archived["all_checks_pass"],
            "complete_source_roots": criterion["simple_roots_complete"],
            "finite_volume": geometry["finite_volume_certified"],
            "noncompact": geometry["compact"] is False,
            "decision": criterion["decision"] == "finite_volume_noncompact",
            "dimension": geometry["hyperbolic_dimension"] == case["dimension"],
            "facets": geometry["facet_count"] == case["facets"],
            "elliptic_edges": criterion["elliptic_edge_count"]
            == case["elliptic_edges"],
            "finite_vertices": geometry["finite_vertex_count"]
            == case["finite_vertices"],
            "cusps": geometry["ideal_vertex_count"] == case["cusps"],
            "no_bad_edges": len(criterion["finite_volume_bad_edges"]) == 0,
            "parabolic_neighborhoods": criterion[
                "parabolic_neighborhood_audit_count"
            ]
            == case["parabolic_audits"],
            "all_parabolic_neighborhoods_pass": geometry["checks"][
                "all_full_parabolic_neighborhoods_pass"
            ],
            "canonical_vinberg_hash": criterion["canonical_full_audit_sha256"]
            == case["vinberg_sha256"],
            "full_reflection_chamber": archived["maximal_reflection_lineage"][
                "reflection_subgroup_index"
            ]
            == 1,
        }
        if not all(geometry_checks.values()):
            raise AssertionError(f"a Vinberg handoff changed: {geometry_checks}")

        results.append(
            {
                "name": case["name"],
                "theorem": (
                    f"Brandhorst--Hofmann {case['bh_id']} gives a genuine "
                    f"order-{case['automorphism_order']} K3 automorphism whose "
                    f"Picard action realizes {case['type_id']} and its finite-volume "
                    f"H{case['dimension']}/{case['facets']} Coxeter section"
                ),
                "brandhorst_hofmann_row": {
                    "id": case["bh_id"],
                    "automorphism_order": case["automorphism_order"],
                    "member": case["member"],
                    "member_sha256": case["member_sha256"],
                    "member_row_count": case["member_row_count"],
                    "row_one_based": case["row_one_based"],
                    "row_sha256_with_newline": case["row_sha256_with_newline"],
                    "database_semantics": (
                        "one moduli space of effective characters parametrizing "
                        "pairs (X,G), as stated in the pinned database README"
                    ),
                    "generic_transcendental_lattice": {
                        "rank_signature_determinant": [8, [2, 6, 0], 4],
                        "smith": [*[1] * 6, 2, 2],
                        "characteristic_polynomial": case[
                            "transcendental_characteristic"
                        ],
                        "gram_digest": case["transcendental_gram_digest"],
                    },
                    "generic_picard_lattice": {
                        "isometry_class": "U+E8(-1)+D4(-1)",
                        "rank_signature_determinant": [14, [1, 13, 0], -4],
                        "smith": [*[1] * 12, 2, 2],
                        "discriminant_q_histogram_mod_2Z": {"0": 1, "1": 3},
                        "gram_digest": case["ns_gram_digest"],
                        "action_digest": case["ns_action_digest"],
                    },
                },
                "exact_picard_action_match": {
                    "archived_type_id": case["type_id"],
                    "picard_action_order": case["picard_action_order"],
                    "picard_characteristic_polynomial": case[
                        "picard_characteristic"
                    ],
                    "fixed_lattice": {
                        "rank": case["fixed_rank"],
                        "signature": [1, case["fixed_rank"] - 1, 0],
                        "determinant": case["fixed_determinant"],
                        "smith": case["fixed_smith"],
                    },
                    "coinvariant_lattice": {
                        "class": case["coinvariant_class"],
                        "gram_in_bh_basis": matrix_payload(coinvariant_gram),
                        "action_in_bh_basis": matrix_payload(coinvariant_action),
                        "canonical_change_of_basis": matrix_payload(
                            canonical_change
                        ),
                        "rootless": True,
                        "rootless_reason": (
                            "the canonical Gram is twice an even negative-definite "
                            "integral Gram, so every nonzero norm is at most -4"
                        ),
                    },
                    "equivariant_orthogonal_decomposition": {
                        "archived": "II_(1,9) [identity] + D4(-1) [diagram action]",
                        "brandhorst_hofmann": "U+E8(-1) [identity] + D4(-1) [same diagram action]",
                        "archived_block_basis_determinant": source_block[
                            "block_basis_determinant"
                        ],
                        "bh_block_basis_determinant": bh_block[
                            "block_basis_determinant"
                        ],
                        "d4_gram": D4_GRAM,
                        "d4_action": matrix_payload(source_block["d4_action"]),
                        "classification_used": (
                            "the even unimodular lattice of signature (1,9) is "
                            "unique, namely U+E8(-1)"
                        ),
                        "actions_conjugate_integrally": True,
                    },
                    "discriminant_form_and_action": {
                        "archived": source_discriminant,
                        "brandhorst_hofmann": bh_discriminant,
                        "basis_conjugator_archived_to_bh": matrix_payload(conjugator),
                        "same_q_action_class": True,
                    },
                    "all_algebraic_checks": algebraic_checks,
                },
                "kahler_and_torelli": {
                    "database_row_is_effective_genuine_automorphism_moduli": True,
                    "generic_picard_is_exactly_rank_14": True,
                    "algebraic_coinvariant_is_rootless": True,
                    "chamber_transport": (
                        "transport the BH invariant Kahler chamber to the complete "
                        "15-wall chamber by a (-2)-Weyl element; Weyl reflections "
                        "act trivially on the discriminant, and the chamber S3 acts "
                        "faithfully on its three nonzero discriminant classes, so a "
                        "final chamber symmetry conjugates to the archived generator"
                    ),
                    "strong_torelli_conclusion": (
                        "after this marking the genuine BH automorphism preserves "
                        "the Roulleau nef chamber and induces the archived fold"
                    ),
                    "global_order_vs_picard_order": [
                        case["automorphism_order"],
                        case["picard_action_order"],
                    ],
                },
                "hyperbolic_polytope": {
                    "dimension": case["dimension"],
                    "facets": case["facets"],
                    "classification": "finite_volume_noncompact",
                    "finite_vertices": case["finite_vertices"],
                    "cusps": case["cusps"],
                    "elliptic_edges": case["elliptic_edges"],
                    "bad_edges": 0,
                    "parabolic_neighborhood_audits": case[
                        "parabolic_audits"
                    ],
                    "canonical_vinberg_audit_sha256": case["vinberg_sha256"],
                    "full_reflection_chamber_index": 1,
                },
                "scope": {
                    "new_positive_k3_automorphism_realization_in_this_project": True,
                    "new_abstract_coxeter_polytope_claimed": False,
                    "worldwide_priority_claimed": False,
                    "original_general_roulleau_surface_action_claimed": False,
                    "generic_ns_on_the_bh_automorphism_locus": True,
                },
            }
        )

    # The two archived discriminant images generate GL(2,2), hence the S3
    # chamber action on the three nonzero discriminant classes is faithful.
    identity = sp.eye(2)
    generated = {tuple(int(x) & 1 for x in identity)}
    changed = True
    while changed:
        changed = False
        for flat in tuple(generated):
            left = sp.Matrix(2, 2, flat)
            for generator in source_discriminant_actions:
                product = (left * generator).applyfunc(lambda value: int(value) & 1)
                key = tuple(int(x) for x in product)
                if key not in generated:
                    generated.add(key)
                    changed = True
    if len(generated) != 6:
        raise AssertionError("the archived S3 discriminant action is not faithful")

    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau_ue8d4.generic_ns_k3_automorphisms.v1",
        "all_checks_pass": True,
        "theorem": (
            "Two Roulleau U+E8+D4 folds previously recorded only as divisor-"
            "lattice symmetries have genuine K3-automorphism realizations: "
            "H11/13 via BH 0.6.3.1 and H12/14 via BH 0.16.0.1"
        ),
        "source_file_sha256": source_hash,
        "brandhorst_hofmann": {
            "article": "Finite subgroups of automorphisms of K3 surfaces",
            "article_doi": BH_ARTICLE_DOI,
            "database_doi": BH_DATA_DOI,
            "database_version": BH_DATA_VERSION,
            "license": "CC-BY-4.0",
            "zip_md5": BH_ZIP_MD5,
            "zip_sha256": BH_ZIP_SHA256,
            "readme_member": BH_README_MEMBER,
            "readme_sha256": BH_README_SHA256,
            "readme_semantics_pinned": (
                "each purely_ns line is a moduli space of pairs (X,G) and an "
                "effective character with respect to an effective subgroup"
            ),
            "external_replay_command": (
                "python3 scripts/certify_roulleau_ue8d4_generic_ns_k3_automorphisms.py "
                "--brandhorst-hofmann-zip /path/to/K3Groups.zip --check-archive"
            ),
        },
        "chamber_transport_criterion": {
            "complete_15_wall_neron_severi_chamber": True,
            "weighted_diagram_group": "S3",
            "discriminant_orthogonal_group": "GL(2,2)=S3",
            "diagram_action_on_discriminant_is_faithful": True,
            "weyl_reflections_act_trivially_on_discriminant": True,
            "conclusion": (
                "every effective BH action of either nontrivial conjugacy class "
                "can be marked so that its Kahler chamber is the archived chamber "
                "and its wall action is the archived C2 or C3 generator"
            ),
        },
        "cases": results,
        "summary": {
            "genuine_k3_automorphism_realizations": 2,
            "dimensions": [11, 12],
            "finite_volume_noncompact": 2,
            "compact": 0,
            "new_abstract_polytope_claims": 0,
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
            raise AssertionError("the archived U+E8+D4 K3 realization changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
