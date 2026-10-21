#!/usr/bin/env python3
"""Certify three full-group enrichments of Roulleau K3 chamber folds.

The previously certified Brandhorst--Hofmann rows were purely
non-symplectic cyclic actions.  Three rows in the *full* effective-group
database add a commuting symplectic involution.  The strongest row has
group C16 x C2 and, on one K3 surface with generic Picard lattice
U+D8(-1)+D4(-1), simultaneously realizes the archived TYPE_2 and TYPE_4
fixed-section polytopes.

The default (offline) certificate checks the pinned chamber certificates
and the exact cyclic-action conjugacies already proved there.  Supplying
``--brandhorst-hofmann-zip`` replays the complete full and purely
non-symplectic database rows, proves the subgroup lineage by exhaustive
table comparison, and directly rechecks the integral fixed/coinvariant
glue conjugacies.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import io
import itertools
import json
from pathlib import Path
from typing import Any, Iterable, Sequence
import zipfile

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_roulleau_r14_c2d_special_k3_order4_lift import (
    _parse_bh_row,
    _rational,
    induced_row_action,
    matrix_digest,
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.certify_roulleau_ud8d4_bh_k3_lifts import (
    apply_discriminant_map,
    apply_reduced_fixed_map,
    cyclotomic_evaluation,
    discriminant_map,
    exact_order,
    glue_pairs,
    verify_fixed_discriminant_isometry,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
P14_SOURCE = (
    ROOT / "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json"
)
UD8D4_SOURCE = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
P14_PURE_CERT = (
    ROOT / "data/enumerations/roulleau_r14_c2d_special_k3_order4_lift.json"
)
UD8D4_PURE_CERT = ROOT / "data/enumerations/roulleau_ud8d4_bh_k3_lifts.json"
OUTPUT = ROOT / "data/enumerations/roulleau_bh_full_group_enrichments.json"

EXPECTED_SOURCE_SHA256 = {
    "p14_source": "5c39aa26fc4e2c7f2af5bb08815532ec9cce7b2222279c65f257e383f356f4ef",
    "ud8d4_source": "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
    "p14_pure_certificate": "eb58b165bff1918dfa7596e557f93c386c138808060b1016667c1fc1ac60b041",
    "ud8d4_pure_certificate": "b5d480d30ba0b8c39bf48d8021127ae114f4919529b620759cd20c061555c4c8",
}

BH_ARTICLE_DOI = "10.1017/fms.2023.50"
BH_DATA_DOI = "10.5281/zenodo.10410278"
BH_DATA_VERSION = "2.1"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_README_SHA256 = "588cfbaed63ade81f8c34ef256240ca26ece3f69074417b5f301f31c418bb9e9"


MEMBERS: dict[str, dict[str, Any]] = {
    "full04": {
        "path": "K3Groups/results/order04.txt",
        "sha256": "828763a2fc6609e7406e22bda835bd32701da8a667134908c22c8c6a6c30a6ca",
        "rows": 670,
    },
    "full08": {
        "path": "K3Groups/results/order08.txt",
        "sha256": "b7fd264843c27f0bac1803f3d5c82f30122531a5504c3d0e18141c5c5220e789",
        "rows": 106,
    },
    "full16": {
        "path": "K3Groups/results/order16.txt",
        "sha256": "89fdf926233d8391093e8d5c386e5507f7d5fb82c34cb7eede78348f11f696ba",
        "rows": 3,
    },
    "pure04": {
        "path": "K3Groups/results/purely_ns/order04.txt",
        "sha256": "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8",
        "rows": 79,
    },
    "pure08": {
        "path": "K3Groups/results/purely_ns/order08.txt",
        "sha256": "719eb569ec3f4d399a9948499a795e273160dfb54bb607e6e18eb31e188138d7",
        "rows": 38,
    },
    "pure16": {
        "path": "K3Groups/results/purely_ns/order16.txt",
        "sha256": "2c8a58d846eb6b563a62b958bd665680790598b829d99723da463d9fa6e315ea",
        "rows": 7,
    },
}


FULL_ROWS: dict[str, dict[str, Any]] = {
    "1.4.3.1": {
        "member": "full04",
        "row": 120,
        "row_sha256": "c3283896d9ce8253b1f058196bfe396b3aeb05f015ac6083fe3350e2885ac398",
        "n": 4,
        "small_group_ids": [[2, 1], [8, 2]],
        "group": "C4 x C2",
        "generic_ns": "U+E8(-1)+A1(-1)^4",
        "ns_q": {"0": 2, "1": 6, "1/2": 4, "3/2": 4},
        "ns_gram_digest": "83d0b16cb4443a2f925796ff53060a29ba58763f198c30226ea7222fa6380355",
        "ns_s_digest": "1795970721125033b45a42e9143f5477019f3ab63cc02c9d6968c769c6492ef5",
        "ns_g_digest": "ce7afd5deeeff58084fcb08bfac81a39c781316c9232c8d325a7c864f8ab74ba",
        "ns_image": "C2 x C2",
        "ns_image_order": 4,
        "ns_kernel": "<g^2> = C2",
        "common_fixed": [6, [1, 5, 0], -64],
    },
    "1.8.1.1": {
        "member": "full08",
        "row": 18,
        "row_sha256": "bf5e3a8aebb370c0ea74ee16fc65e9d52384a4b6627b747033ff986be6290b65",
        "n": 8,
        "small_group_ids": [[2, 1], [16, 5]],
        "group": "C8 x C2",
        "generic_ns": "U+D8(-1)+D4(-1)",
        "ns_q": {"0": 6, "1": 10},
        "ns_gram_digest": "9e3602106412256d72f7e5ce220de118a8ca69f43bff9837d2b6c002646c2d88",
        "ns_s_digest": "c6b8f5c23ade0b087183eed59deeb5862c203790580984d726705fe3d47cb49a",
        "ns_g_digest": "7c8af6a93d9d7b7c57696f7548e233f2ea975120143cc24e1c3a86e32792f30f",
        "ns_image": "C2 x C2",
        "ns_image_order": 4,
        "ns_kernel": "<g^2> = C4",
        "common_fixed": [6, [1, 5, 0], -16],
    },
    "1.16.0.1": {
        "member": "full16",
        "row": 1,
        "row_sha256": "5bd04356c93bf66261e6e6303567da9174f97addcaf310d6b0f46f496837a509",
        "n": 16,
        "small_group_ids": [[2, 1], [32, 16]],
        "group": "C16 x C2",
        "generic_ns": "U+D8(-1)+D4(-1)",
        "ns_q": {"0": 6, "1": 10},
        "ns_gram_digest": "a21d4ec45fbb2a5e886e96eece885ec848c5555d8cbfde2085b18fcb26a20382",
        "ns_s_digest": "d8ecdd7abf1c8a4aaa7fcd05b26f696981fcc92e00bdbf3d3f5f8843a154a34e",
        "ns_g_digest": "9683b059142f7d1b2df4239d6b3c2d37678bd5fab8de4a140831d8b4f9251081",
        "ns_image": "C4 x C2",
        "ns_image_order": 8,
        "ns_kernel": "<g^4> = C4",
        "common_fixed": [6, [1, 5, 0], -16],
    },
}


# The coordinate convention is P*G_archive*P^t=G_BH and
# P*A_archive=A_BH*P on the negative-definite coinvariant lattice.
ENRICHMENTS: tuple[dict[str, Any], ...] = (
    {
        "name": "p14_c2d_c4xc2",
        "full_id": "1.4.3.1",
        "element": "g",
        "element_order": 4,
        "element_ns_order": 2,
        "pure_id": "0.4.3.4",
        "pure_member": "pure04",
        "pure_row": 28,
        "pure_row_sha256": "303fd69804142bc57a27799551fa96de8dc10bfe764449e48e1ffe5d55cd7cf3",
        "source": "p14",
        "case_id": "R14-UE8A1x4-C2D-H11-19",
        "archive_class": None,
        "type_id": "R14-UE8A1x4-C2D-H11-19",
        "dimension": 11,
        "facets": 19,
        "fixed_rank": 12,
        "fixed_det": -16,
        "fixed_smith": [*[1] * 10, 4, 4],
        "fixed_q": {"0": 2, "1": 2, "3/2": 4, "3/4": 4, "7/4": 4},
        "coinvariant_gram": [[-4, 4], [4, -8]],
        "coinvariant_action": [[-1, 0], [0, -1]],
        "coin_intertwiner": [[-1, 0], [0, 1]],
        "fixed_discriminant_map": [[0, 1], [1, 2]],
    },
    {
        "name": "ud8d4_type2_c8xc2",
        "full_id": "1.8.1.1",
        "element": "g",
        "element_order": 8,
        "element_ns_order": 2,
        "pure_id": "0.8.1.3",
        "pure_member": "pure08",
        "pure_row": 10,
        "pure_row_sha256": "b543f8d219c913c3d0834188aea0e62dfe724b562569a017e887b6d199d59bd9",
        "source": "ud8d4",
        "case_id": "R14-D8D4-T2-H11-16",
        "archive_class": 6,
        "type_id": "R14_UD8D4_TYPE_2",
        "dimension": 11,
        "facets": 16,
        "fixed_rank": 12,
        "fixed_det": -16,
        "fixed_smith": [*[1] * 10, 4, 4],
        "fixed_q": {"0": 2, "1": 2, "1/4": 4, "3/2": 4, "5/4": 4},
        "coinvariant_gram": [[-4, 0], [0, -4]],
        "coinvariant_action": [[-1, 0], [0, -1]],
        "coin_intertwiner": [[-1, 0], [0, -1]],
        "fixed_discriminant_map": [[0, 1], [1, 0]],
    },
    {
        "name": "ud8d4_type2_c16xc2",
        "full_id": "1.16.0.1",
        "element": "g2",
        "element_order": 8,
        "element_ns_order": 2,
        "pure_id": "0.8.1.3",
        "pure_member": "pure08",
        "pure_row": 10,
        "pure_row_sha256": "b543f8d219c913c3d0834188aea0e62dfe724b562569a017e887b6d199d59bd9",
        "source": "ud8d4",
        "case_id": "R14-D8D4-T2-H11-16",
        "archive_class": 6,
        "type_id": "R14_UD8D4_TYPE_2",
        "dimension": 11,
        "facets": 16,
        "fixed_rank": 12,
        "fixed_det": -16,
        "fixed_smith": [*[1] * 10, 4, 4],
        "fixed_q": {"0": 2, "1": 2, "1/4": 4, "3/2": 4, "5/4": 4},
        "coinvariant_gram": [[-4, 4], [4, -8]],
        "coinvariant_action": [[-1, 0], [0, -1]],
        "coin_intertwiner": [[-1, 0], [1, -1]],
        "fixed_discriminant_map": [[0, 1], [3, 0]],
    },
    {
        "name": "ud8d4_type4_c16xc2",
        "full_id": "1.16.0.1",
        "element": "sg",
        "element_order": 16,
        "element_ns_order": 4,
        "pure_id": "0.16.0.3",
        "pure_member": "pure16",
        "pure_row": 3,
        "pure_row_sha256": "2c78a6da28d3672826467dce6a7594d0adefec5f8d4fa2df869795751a8ffebb",
        "source": "ud8d4",
        "case_id": "R14-D8D4-T4-H10-14",
        "archive_class": 14,
        "type_id": "R14_UD8D4_TYPE_4",
        "dimension": 10,
        "facets": 14,
        "fixed_rank": 11,
        "fixed_det": 8,
        "fixed_smith": [*[1] * 10, 8],
        "fixed_q": {"0": 2, "11/8": 2, "3/2": 2, "3/8": 2},
        "coinvariant_gram": [[-4, 2, 4], [2, -4, 0], [4, 0, -8]],
        "coinvariant_action": [[0, 1, 0], [1, 0, 1], [-2, -2, -1]],
        "coin_intertwiner": [[-1, 0, 1], [0, 1, -1], [1, -1, -1]],
        "fixed_discriminant_map": [[7]],
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(value)
            for value in matrix.row(index)
        ]
        for index in range(matrix.rows)
    ]


def _ambient_row(line: str) -> dict[str, Any]:
    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    source_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    k3_gram = sp.simplify(basis * source_gram * basis.T)
    symplectic = [
        sp.simplify(basis * sp.Matrix(22, 22, action) * basis.inv())
        for action in data["symplectic"]
    ]
    generator = sp.simplify(
        basis * sp.Matrix(22, 22, data["distinguished_generator"]) * basis.inv()
    )
    if not (
        k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and all(
            action * k3_gram * action.T == k3_gram
            for action in [*symplectic, generator]
        )
    ):
        raise AssertionError("a full Brandhorst--Hofmann K3 row changed")
    return {
        "data": data,
        "gram": k3_gram,
        "symplectic": symplectic,
        "generator": generator,
    }


def _derive_full_row(line: str, specification: dict[str, Any]) -> dict[str, Any]:
    row = _ambient_row(line)
    data = row["data"]
    gram = row["gram"]
    generator = row["generator"]
    symplectic = row["symplectic"]
    n = specification["n"]
    if len(symplectic) != 1:
        raise AssertionError("the pinned full row no longer has one symplectic generator")
    s = symplectic[0]
    equations = cyclotomic_evaluation(generator, n).row_join(s - sp.eye(22))
    transcendental_basis = saturated_row_kernel(equations)
    ns_basis = saturated_row_kernel(gram * transcendental_basis.T)
    ns_gram = sp.simplify(ns_basis * gram * ns_basis.T)
    ns_s = induced_row_action(ns_basis, s)
    ns_g = induced_row_action(ns_basis, generator)

    full_elements = [s**a * generator**b for a in range(2) for b in range(n)]
    ns_elements = [ns_s**a * ns_g**b for a in range(2) for b in range(n)]
    distinct_full = {matrix_digest(action) for action in full_elements}
    distinct_ns = {matrix_digest(action) for action in ns_elements}
    common_fixed = saturated_row_kernel(
        (s - sp.eye(22)).row_join(generator - sp.eye(22))
    )
    common_fixed_gram = sp.simplify(common_fixed * gram * common_fixed.T)
    symplectic_fixed = saturated_row_kernel(s - sp.eye(22))
    symplectic_coinvariant = saturated_row_kernel(gram * symplectic_fixed.T)
    symplectic_coinvariant_gram = sp.simplify(
        symplectic_coinvariant * gram * symplectic_coinvariant.T
    )
    small_ids = [list(pair) for pair in dict(data["properties"])["small group ids"]]
    if not (
        data["id"] in FULL_ROWS
        and data["transcendental_value"] == n
        and small_ids == specification["small_group_ids"]
        and exact_order(s, 2)
        and exact_order(generator, n)
        and s * generator == generator * s
        and len(distinct_full) == 2 * n
        and len(distinct_ns) == specification["ns_image_order"]
        and transcendental_basis.rows == 8
        and exact_inertia(transcendental_basis * gram * transcendental_basis.T)
        == (2, 6, 0)
        and abs(int((transcendental_basis * gram * transcendental_basis.T).det()))
        == 16
        and ns_gram.shape == (14, 14)
        and ns_gram.det() == -16
        and exact_inertia(ns_gram) == (1, 13, 0)
        and smith_diagonal(ns_gram) == [*[1] * 10, *[2] * 4]
        and _q_histogram(ns_gram) == specification["ns_q"]
        and matrix_digest(ns_gram) == specification["ns_gram_digest"]
        and matrix_digest(ns_s) == specification["ns_s_digest"]
        and matrix_digest(ns_g) == specification["ns_g_digest"]
        and [
            common_fixed.rows,
            list(exact_inertia(common_fixed_gram)),
            int(common_fixed_gram.det()),
        ]
        == specification["common_fixed"]
        and symplectic_coinvariant.rows == 8
        and exact_inertia(symplectic_coinvariant_gram) == (0, 8, 0)
        and symplectic_coinvariant_gram.det() == 256
        and all(
            sp.Rational(value).q == 1 for value in symplectic_coinvariant_gram / 2
        )
        and all(
            symplectic_coinvariant_gram[index, index] % 4 == 0
            for index in range(8)
        )
    ):
        raise AssertionError("a pinned full-group structure or generic NS changed")
    row.update(
        {
            "transcendental_basis": transcendental_basis,
            "ns_basis": ns_basis,
            "ns_gram": ns_gram,
            "ns_s": ns_s,
            "ns_g": ns_g,
            "common_fixed_gram": common_fixed_gram,
            "symplectic_coinvariant_gram": symplectic_coinvariant_gram,
        }
    )
    return row


def _selected_actions(row: dict[str, Any], expression: str) -> tuple[sp.Matrix, sp.Matrix]:
    s = row["symplectic"][0]
    g = row["generator"]
    ns_s = row["ns_s"]
    ns_g = row["ns_g"]
    if expression == "g":
        return g, ns_g
    if expression == "g2":
        return g**2, ns_g**2
    if expression == "sg":
        return s * g, ns_s * ns_g
    raise AssertionError("unknown selected group element")


def _fixed_type(gram: sp.Matrix, action: sp.Matrix) -> tuple[Any, ...]:
    fixed_basis = saturated_row_kernel(action - sp.eye(action.rows))
    fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
    return (
        fixed_basis.rows,
        int(fixed_gram.det()) if fixed_basis.rows else 1,
        tuple(smith_diagonal(fixed_gram)),
    )


def _cyclic_local_type(
    gram: sp.Matrix, action: sp.Matrix, order: int
) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        (power, *_fixed_type(gram, action**power))
        for power in range(1, order + 1)
        if order % power == 0
    )


def _full_character_fingerprint(
    gram: sp.Matrix, s: sp.Matrix, g: sp.Matrix, order: int
) -> tuple[tuple[tuple[Any, ...], ...], ...]:
    # The unordered pair treats g^b and s*g^b as the same H_s-coset.
    return tuple(
        tuple(sorted((_fixed_type(gram, g**power), _fixed_type(gram, s * g**power))))
        for power in range(order)
    )


def _archive_pair(
    specification: dict[str, Any], p14: dict[str, Any], ud8d4: dict[str, Any]
) -> dict[str, Any]:
    if specification["source"] == "p14":
        case = next(
            row for row in p14["cases"] if row["id"] == specification["case_id"]
        )
        lane = p14["source_reconstruction"]["authoritative_lane"]
        ns_gram = sp.Matrix(lane["basis_gram"])
        walls = sp.Matrix(lane["wall_coordinate_columns"])
        permutation = tuple(
            value - 1 for value in case["subgroup"]["elements_one_based"][1]
        )
        action = sp.Matrix.hstack(
            *(walls[:, permutation[index]] for index in range(14))
        ).T
        vinberg = case["geometry"]["vinberg"]
        geometry = {
            "classification": "finite_volume_noncompact",
            "dimension": vinberg["dimension"],
            "facets": vinberg["facet_count"],
            "finite_vertices": vinberg["finite_vertex_count"],
            "cusps": vinberg["ideal_vertex_count"],
            "elliptic_edges": vinberg["elliptic_edge_count"],
            "bad_edges": vinberg["bad_edge_count"],
            "bad_parabolic_neighborhoods": vinberg[
                "bad_parabolic_neighborhood_count"
            ],
        }
    else:
        archived = next(
            row
            for row in ud8d4["high_fixed_sections"]
            if row["class_index"] == specification["archive_class"]
        )
        diagram = next(
            row
            for row in ud8d4["abstract_diagram_types"]
            if row["type_id"] == specification["type_id"]
        )
        ns_gram = sp.Matrix(ud8d4["source"]["basis_gram"])
        action = sp.Matrix(
            archived["ns_action"]["representative_generator_action_matrices"][0]
        )
        criterion = diagram["geometry"]["vinberg_volume_criterion"]
        geometry = {
            "classification": diagram["volume_type"],
            "dimension": diagram["hyperbolic_dimension"],
            "facets": diagram["facet_count"],
            "finite_vertices": criterion["finite_vertex_count"],
            "cusps": criterion["ideal_vertex_count"],
            "elliptic_edges": criterion["elliptic_edge_count"],
            "bad_edges": len(criterion["finite_volume_bad_edges"]),
            "bad_parabolic_neighborhoods": 0,
        }
    fixed_basis = saturated_row_kernel(action - sp.eye(14))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    return {
        "ns_gram": ns_gram,
        "action": action,
        "fixed_basis": fixed_basis,
        "fixed_gram": sp.simplify(fixed_basis * ns_gram * fixed_basis.T),
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": sp.simplify(
            coinvariant_basis * ns_gram * coinvariant_basis.T
        ),
        "coinvariant_action": induced_row_action(coinvariant_basis, action),
        "glue": glue_pairs(fixed_basis, coinvariant_basis, ns_gram),
        "geometry": geometry,
    }


def _direct_pair_audit(
    row: dict[str, Any], specification: dict[str, Any], archive: dict[str, Any]
) -> dict[str, Any]:
    global_action, ns_action = _selected_actions(row, specification["element"])
    if not (
        exact_order(global_action, specification["element_order"])
        and exact_order(ns_action, specification["element_ns_order"])
    ):
        raise AssertionError("a selected full-group element changed order")
    ns_gram = row["ns_gram"]
    fixed_basis = saturated_row_kernel(ns_action - sp.eye(14))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * ns_gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, ns_action)
    bh_glue = glue_pairs(fixed_basis, coinvariant_basis, ns_gram)
    if not (
        fixed_basis.rows == specification["fixed_rank"]
        and fixed_gram.det() == specification["fixed_det"]
        and smith_diagonal(fixed_gram) == specification["fixed_smith"]
        and _q_histogram(fixed_gram) == specification["fixed_q"]
        and coinvariant_gram == sp.Matrix(specification["coinvariant_gram"])
        and coinvariant_action == sp.Matrix(specification["coinvariant_action"])
    ):
        raise AssertionError("a full-row selected Picard action changed")

    archived_fixed = archive["fixed_gram"]
    archived_coin = archive["coinvariant_gram"]
    archived_coin_action = archive["coinvariant_action"]
    coin_intertwiner = sp.Matrix(specification["coin_intertwiner"])
    if not (
        abs(int(coin_intertwiner.det())) == 1
        and coin_intertwiner * archived_coin * coin_intertwiner.T
        == coinvariant_gram
        and coin_intertwiner * archived_coin_action
        == coinvariant_action * coin_intertwiner
    ):
        raise AssertionError("the direct coinvariant intertwiner changed")
    fixed_map = sp.Matrix(specification["fixed_discriminant_map"])
    fixed_certificate = verify_fixed_discriminant_isometry(
        archived_fixed, fixed_gram, fixed_map
    )
    coin_map = discriminant_map(
        coin_intertwiner.inv(), archived_coin, coinvariant_gram
    )
    target_coin_moduli, _, _ = _discriminant_data(coinvariant_gram)
    bh_graph = {coin: fixed for fixed, coin in bh_glue}
    for archived_fixed_element, archived_coin_element in archive["glue"]:
        target_fixed = apply_reduced_fixed_map(
            archived_fixed_element, archived_fixed, fixed_gram, fixed_map
        )
        target_coin = apply_discriminant_map(
            coin_map, archived_coin_element, target_coin_moduli
        )
        if bh_graph.get(target_coin) != target_fixed:
            raise AssertionError("the direct pair maps do not carry the glue graph")

    fixed_length = len(
        [value for value in smith_diagonal(fixed_gram) if value > 1]
    )
    if not (
        exact_inertia(fixed_gram) == (1, fixed_gram.rows - 1, 0)
        and fixed_gram.rows >= fixed_length + 2
        and exact_inertia(coinvariant_gram)
        == (0, coinvariant_gram.rows, 0)
        and all(sp.Rational(value).q == 1 for value in coinvariant_gram / 2)
        and all(
            coinvariant_gram[index, index] % 4 == 0
            for index in range(coinvariant_gram.rows)
        )
    ):
        raise AssertionError("the Nikulin lift/rootless gate changed")
    return {
        "name": specification["name"],
        "global_element_order": specification["element_order"],
        "picard_element_order": specification["element_ns_order"],
        "fixed_discriminant_isometry": fixed_certificate,
        "coinvariant_intertwiner": matrix_payload(coin_intertwiner),
        "coinvariant_discriminant_map": matrix_payload(coin_map),
        "primitive_glue_order": len(bh_glue),
        "all_glue_pairs_carried": True,
        "nikulin_fixed_isometry_lifts": True,
        "full_picard_pairs_integrally_conjugate": True,
        "algebraic_coinvariant_rootless": True,
    }


def _read_member(archive: zipfile.ZipFile, key: str) -> list[str]:
    specification = MEMBERS[key]
    payload = archive.read(specification["path"])
    if sha256_bytes(payload) != specification["sha256"]:
        raise AssertionError(f"Brandhorst--Hofmann member {key} changed")
    lines = payload.decode("utf-8").splitlines()
    if len(lines) != specification["rows"]:
        raise AssertionError(f"Brandhorst--Hofmann member {key} row count changed")
    return lines


def replay_external_zip(
    path: Path, p14: dict[str, Any], ud8d4: dict[str, Any]
) -> dict[str, Any]:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(payload) == BH_ZIP_SHA256
    ):
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        if sha256_bytes(archive.read("K3Groups/readme.txt")) != BH_README_SHA256:
            raise AssertionError("the Brandhorst--Hofmann database readme changed")
        lines = {key: _read_member(archive, key) for key in MEMBERS}
        derived: dict[str, dict[str, Any]] = {}
        for identifier, specification in FULL_ROWS.items():
            line = lines[specification["member"]][specification["row"] - 1]
            if sha256_bytes((line + "\n").encode()) != specification["row_sha256"]:
                raise AssertionError("a pinned full-group row changed")
            derived[identifier] = _derive_full_row(line, specification)

        direct_audits = []
        pure_identifications = []
        for specification in ENRICHMENTS:
            row = derived[specification["full_id"]]
            archive_pair = _archive_pair(specification, p14, ud8d4)
            direct_audits.append(_direct_pair_audit(row, specification, archive_pair))
            global_action, _ = _selected_actions(row, specification["element"])
            target_type = _cyclic_local_type(
                row["gram"], global_action, specification["element_order"]
            )
            pure_lines = lines[specification["pure_member"]]
            pinned_line = pure_lines[specification["pure_row"] - 1]
            if sha256_bytes((pinned_line + "\n").encode()) != specification[
                "pure_row_sha256"
            ]:
                raise AssertionError("a pinned purely non-symplectic row changed")
            matches = []
            for row_number, line in enumerate(pure_lines, 1):
                candidate = _ambient_row(line)
                if _cyclic_local_type(
                    candidate["gram"],
                    candidate["generator"],
                    specification["element_order"],
                ) == target_type:
                    matches.append([row_number, candidate["data"]["id"]])
            if matches != [[specification["pure_row"], specification["pure_id"]]]:
                raise AssertionError("the exhaustive pure-row identification changed")
            pure_identifications.append(
                {
                    "enrichment": specification["name"],
                    "complete_member": MEMBERS[specification["pure_member"]]["path"],
                    "complete_row_count": len(pure_lines),
                    "unique_matching_effective_character": specification["pure_id"],
                    "row_one_based": specification["pure_row"],
                    "small_local_type": [
                        [power, rank, determinant, list(smith)]
                        for power, rank, determinant, smith in target_type
                    ],
                }
            )

        # The C16 x C2 row contains, as <s,g^2>, the unique full C8 x C2
        # effective character 1.8.1.1.  This closes the nested enrichment.
        row16 = derived["1.16.0.1"]
        hierarchy_target = _full_character_fingerprint(
            row16["gram"], row16["symplectic"][0], row16["generator"] ** 2, 8
        )
        hierarchy_hits = []
        for row_number, line in enumerate(lines["full08"], 1):
            candidate = _ambient_row(line)
            properties = dict(candidate["data"]["properties"])
            if not (
                len(candidate["symplectic"]) == 1
                and properties.get("small group ids") == ((2, 1), (16, 5))
            ):
                continue
            if _full_character_fingerprint(
                candidate["gram"],
                candidate["symplectic"][0],
                candidate["generator"],
                8,
            ) == hierarchy_target:
                hierarchy_hits.append([row_number, candidate["data"]["id"]])
        if hierarchy_hits != [[18, "1.8.1.1"]]:
            raise AssertionError("the nested C8 x C2 subgroup lineage changed")

    return {
        "zip_sha256": BH_ZIP_SHA256,
        "full_rows_replayed": list(FULL_ROWS),
        "direct_integral_pair_audits": direct_audits,
        "exhaustive_pure_subgroup_identifications": pure_identifications,
        "nested_full_group_identification": {
            "subgroup": "<s,g^2> in BH 1.16.0.1",
            "group": "C8 x C2",
            "complete_order08_rows_checked": MEMBERS["full08"]["rows"],
            "unique_matching_full_effective_character": "1.8.1.1",
            "row_one_based": 18,
        },
    }


def _offline_sources() -> tuple[dict[str, str], dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = {
        "p14_source": P14_SOURCE,
        "ud8d4_source": UD8D4_SOURCE,
        "p14_pure_certificate": P14_PURE_CERT,
        "ud8d4_pure_certificate": UD8D4_PURE_CERT,
    }
    hashes = {key: sha256(path) for key, path in paths.items()}
    if hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned chamber or pure-lift certificate changed")
    return (
        hashes,
        json.loads(P14_SOURCE.read_text(encoding="utf-8")),
        json.loads(UD8D4_SOURCE.read_text(encoding="utf-8")),
        json.loads(P14_PURE_CERT.read_text(encoding="utf-8")),
        json.loads(UD8D4_PURE_CERT.read_text(encoding="utf-8")),
    )


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    hashes, p14, ud8d4, p14_pure, ud8d4_pure = _offline_sources()
    ud_pure_by_id = {
        row["brandhorst_hofmann_row"]["id"]: row
        for row in ud8d4_pure["cases"]
    }
    if not (
        p14["all_checks_pass"]
        and ud8d4["all_checks_pass"]
        and p14_pure["all_checks_pass"]
        and ud8d4_pure["all_checks_pass"]
        and p14_pure["brandhorst_hofmann_cross_identification"][
            "target_NS_pair_exactly_U_E8_A1x4_with_double_transposition"
        ]
        and p14_pure["brandhorst_hofmann_cross_identification"][
            "unique_row_with_fixed_determinant_minus_16"
        ]
        == "0.4.3.4"
        and ud_pure_by_id["0.8.1.3"]["exact_action_conjugacy"]["conclusion"].startswith(
            "the fixed and coinvariant isometries"
        )
        and ud_pure_by_id["0.16.0.3"]["exact_action_conjugacy"]["conclusion"].startswith(
            "the fixed and coinvariant isometries"
        )
    ):
        raise AssertionError("an exact pure-action lineage certificate changed")

    results = []
    for specification in ENRICHMENTS:
        archive = _archive_pair(specification, p14, ud8d4)
        geometry = archive["geometry"]
        if not (
            geometry["classification"] == "finite_volume_noncompact"
            and geometry["dimension"] == specification["dimension"]
            and geometry["facets"] == specification["facets"]
            and geometry["bad_edges"] == 0
            and geometry["bad_parabolic_neighborhoods"] == 0
        ):
            raise AssertionError("an upstream Vinberg finite-volume theorem changed")
        full = FULL_ROWS[specification["full_id"]]
        results.append(
            {
                "name": specification["name"],
                "full_effective_character": {
                    "brandhorst_hofmann_id": specification["full_id"],
                    "group": full["group"],
                    "symplectic_kernel_Hs": "C2 generated by s",
                    "distinguished_quotient_order": full["n"],
                    "generic_picard_lattice": full["generic_ns"],
                    "picard_image": full["ns_image"],
                    "picard_image_order": full["ns_image_order"],
                    "kernel_of_H_to_O_NS": full["ns_kernel"],
                    "common_fixed_lattice_rank_signature_determinant": full[
                        "common_fixed"
                    ],
                    "effective_group_row": True,
                    "common_invariant_kahler_class_exists": True,
                },
                "selected_pure_element": {
                    "word": specification["element"],
                    "global_order": specification["element_order"],
                    "picard_order": specification["element_ns_order"],
                    "unique_pure_effective_character": specification["pure_id"],
                    "exhaustive_identification_in_external_replay": True,
                    "exact_picard_pair_conjugacy_in_pure_certificate": True,
                    "direct_fixed_coinvariant_glue_conjugacy_available": True,
                },
                "archived_polytope": {
                    "id": specification["case_id"],
                    "type_id": specification["type_id"],
                    **geometry,
                },
                "geometric_handoff": {
                    "generic_NS_equals_complete_Roulleau_divisor_lattice": True,
                    "database_effectivity_supplies_actual_K3_automorphisms": True,
                    "common_group_invariant_kahler_chamber_exists": True,
                    "weyl_marking_to_archived_complete_nef_chamber": True,
                    "selected_diagram_class_identified_exactly": True,
                    "conclusion": (
                        "the selected element of the full K3 automorphism group "
                        "realizes the archived finite-volume fixed-section polytope"
                    ),
                },
            }
        )

    external = replay_external_zip(external_zip, p14, ud8d4) if external_zip else None
    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau.bh_full_group_enrichments.v1",
        "all_checks_pass": True,
        "theorem": (
            "Three effective Brandhorst--Hofmann full-group classes enrich the "
            "previously certified cyclic Roulleau realizations.  In particular, "
            "one K3 surface with generic NS=U+D8(-1)+D4(-1) and an effective "
            "C16 x C2 action simultaneously realizes TYPE_2 H11/16 through g^2 "
            "and TYPE_4 H10/14 through s*g."
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
            "readme_semantics": (
                "each row is an effective character parametrizing K3 pairs (X,G); "
                "the listed symplectic matrices generate H_s and, together with "
                "the distinguished generator, the actual automorphism group H"
            ),
            "full_rows": FULL_ROWS,
            "member_provenance": MEMBERS,
            "external_replay_performed": external is not None,
            "external_replay": external,
        },
        "realizations": results,
        "nested_positive_lineage": {
            "p14": "0.4.3.4 < 1.4.3.1 (C4 < C4 x C2)",
            "ud8d4": (
                "0.8.1.3 < 1.8.1.1 = <s,g^2> in 1.16.0.1; "
                "also <s*g> in 1.16.0.1 is 0.16.0.3"
            ),
            "simultaneous_surface_result": (
                "BH 1.16.0.1 realizes both R14-D8D4-T2-H11-16 and "
                "R14-D8D4-T4-H10-14 on one generic-NS K3 surface"
            ),
        },
        "novelty_scope": {
            "new_project_level_full_group_geometric_realizations": 3,
            "selected_element_realizations_certified": 4,
            "distinct_archived_polytopes_enriched": 3,
            "simultaneous_two_polytope_K3_realization": True,
            "new_abstract_coxeter_polytope_claimed": False,
            "new_dimension_claimed": False,
            "new_arithmetic_class_claimed": False,
            "worldwide_priority_claimed": False,
            "interpretation": (
                "a stronger group-enriched realization of known/certified "
                "polytopes, extracted from a published exhaustive K3 database"
            ),
        },
        "summary": {
            "full_effective_group_rows": 3,
            "cyclic_subgroup_realizations": 4,
            "finite_volume_noncompact": 4,
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
            raise AssertionError("the archived full-group enrichment certificate changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
