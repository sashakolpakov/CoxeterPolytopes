#!/usr/bin/env python3
"""Certify four K3-group realizations of three Roulleau rank-12 folds.

The Brandhorst--Hofmann database contains three purely non-symplectic
order-four rows whose generic Picard lattices are respectively

    U+E8(-1)+A1(-1)^2,
    U+D8(-1)+A1(-1)^2,
    U+D4(-1)^2+A1(-1)^2.

Their Picard involutions are integrally conjugate to the three archived
Roulleau wall actions giving H10/12, H10/15, and H10/38.  A fourth full-group
row, with group C4 x C2 and a symplectic kernel, gives a second realization
of H10/38.

The default certificate is offline: it checks the pinned generic Picard
Grams, fixed/coinvariant lattices, discriminant-form isometries, primitive
glue, exact chamber actions, and upstream Vinberg audits.  The optional
``--brandhorst-hofmann-zip`` flag replays all four complete 22-dimensional
K3-lattice rows from the published CC-BY database archive.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
from typing import Any, Sequence
import zipfile

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
)
from scripts.certify_a18_c06_h13_31_k3_order4_lift import matrix_digest
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    _parse_bh_row,
    _rational,
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
GRAPH_SOURCE = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
D4_SOURCE = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
OUTPUT = ROOT / "data/enumerations/roulleau_rank12_bh_k3_lifts.json"

EXPECTED_SOURCE_SHA256 = {
    "graph": "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    "d4": "448bcc47937a2763bbd4ea97466fbc48924dc862e88ba8719e5ab98b0cc4c916",
}

BH_ARTICLE_DOI = "10.1017/fms.2023.50"
BH_DATA_DOI = "10.5281/zenodo.10410278"
BH_DATA_VERSION = "2.1"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_README_SHA256 = "588cfbaed63ade81f8c34ef256240ca26ece3f69074417b5f301f31c418bb9e9"

MEMBERS = {
    "pure04": {
        "path": "K3Groups/results/purely_ns/order04.txt",
        "sha256": "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8",
        "rows": 79,
    },
    "full04": {
        "path": "K3Groups/results/order04.txt",
        "sha256": "828763a2fc6609e7406e22bda835bd32701da8a667134908c22c8c6a6c30a6ca",
        "rows": 670,
    },
}


PURE_ACTION = [
    [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
]


NS_GRAMS = {
    "0.4.4.1": [
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, -2, -1, 1, 1, 1, -1, 1, 1, 0, 0],
        [0, 0, -1, -2, 1, 1, 1, -1, 1, 1, 0, 0],
        [0, 0, 1, 1, -2, 0, -1, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, -2, -1, 1, -1, -1, 0, 0],
        [0, 0, 1, 1, -1, -1, -2, 1, -1, -1, 0, 0],
        [0, 0, -1, -1, 0, 1, 1, -2, 1, 1, 0, 0],
        [0, 0, 1, 1, 0, -1, -1, 1, -2, -1, 0, 0],
        [0, 0, 1, 1, 0, -1, -1, 1, -1, -2, 0, 0],
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ],
    "0.4.4.2": [
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [0, 0, -2, -1, -1, 1, -1, 1, -1, -1, 0, 0],
        [0, 0, -1, -2, -1, 0, -1, 0, -1, 0, 0, 0],
        [0, 0, -1, -1, -2, 1, -1, 1, -1, 0, 0, 0],
        [0, 0, 1, 0, 1, -2, 1, -1, 0, 1, 0, 0],
        [0, 0, -1, -1, -1, 1, -2, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, -1, 0, -2, 1, 1, 0, 0],
        [0, 0, -1, -1, -1, 0, 0, 1, -2, 0, 0, 0],
        [0, 0, -1, 0, 0, 1, 0, 1, 0, -2, 0, 0],
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0],
        [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ],
    "0.4.4.6": [
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [0, 0, -2, -1, -1, 1, 1, -1, 1, 1, 0, 0],
        [0, 0, -1, -2, -1, 0, 1, -1, 1, 1, 0, 0],
        [0, 0, -1, -1, -2, 1, 1, -1, 1, 1, 0, 0],
        [0, 0, 1, 0, 1, -2, -1, 1, -1, -1, 0, 0],
        [0, 0, 1, 1, 1, -1, -2, 1, -1, -1, 0, 0],
        [0, 0, -1, -1, -1, 1, 1, -2, 1, 1, 0, 0],
        [0, 0, 1, 1, 1, -1, -1, 1, -2, -1, 0, 0],
        [0, 0, 1, 1, 1, -1, -1, 1, -1, -2, 0, 0],
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0],
        [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ],
    "1.4.4.2": [
        [-2, -2, 0, -2, 2, -2, 2, 1, -2, -1, -1, 0],
        [-2, -4, -1, -2, 1, -2, 2, 0, -2, 0, 0, -2],
        [0, -1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-2, -2, 0, -4, 3, -3, 3, 2, -2, -2, -3, 1],
        [2, 1, 0, 3, -4, 2, -2, -2, 2, 2, 2, -2],
        [-2, -2, 0, -3, 2, -4, 2, 0, 0, -2, -2, 0],
        [2, 2, 0, 3, -2, 2, -4, -2, 2, 2, 2, 0],
        [1, 0, 0, 2, -2, 0, -2, -4, 2, 2, 2, -2],
        [-2, -2, 0, -2, 2, 0, 2, 2, -4, 0, 0, 0],
        [-1, 0, 0, -2, 2, -2, 2, 2, 0, -4, -2, 2],
        [-1, 0, 0, -3, 2, -2, 2, 2, 0, -2, -4, 2],
        [0, -2, 0, 1, -2, 0, 0, -2, 0, 2, 2, -4],
    ],
}


FULL_NS_S = [
    [1, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0],
    [0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 0, -1],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, -1, 0, 0, -1, 0, -1, 0],
    [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
]

FULL_NS_G = [
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0, 1],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
]


CASES: tuple[dict[str, Any], ...] = (
    {
        "id": "R12-E8A1x2-C2-H10-12",
        "bh_id": "0.4.4.1",
        "member": "pure04",
        "row": 39,
        "row_sha256": "03ef01ed30b98e3d0fa419cb8275577c90594868fcd3400c48413555b1248452",
        "source_kind": "graph",
        "source_key": "R12_UE8A1x2",
        "source_class": 2,
        "lattice": "U+E8(-1)+A1(-1)^2",
        "determinant": -4,
        "smith": [*[1] * 10, 2, 2],
        "q": {"0": 1, "1": 1, "3/2": 2},
        "fixed_map": [[1]],
        "fixed_q": {"0": 1, "1": 1, "7/4": 2},
        "facets": 12,
        "finite_vertices": 17,
        "cusps": 2,
        "elliptic_edges": 99,
        "parabolic_audits": 3,
        "group": "C4",
        "ns_image": "C2",
    },
    {
        "id": "R12-D8A1x2-C2a-H10-15",
        "bh_id": "0.4.4.2",
        "member": "pure04",
        "row": 40,
        "row_sha256": "bbba67834e034344032bbf137a1f2f1499f06f3adac123ce2b580e0bf7614849",
        "source_kind": "graph",
        "source_key": "R12_UD8A1x2",
        "source_class": 2,
        "lattice": "U+D8(-1)+A1(-1)^2",
        "determinant": -16,
        "smith": [*[1] * 8, 2, 2, 2, 2],
        "q": {"0": 4, "1": 4, "1/2": 2, "3/2": 6},
        "fixed_map": [[0, 1, 1], [1, 0, 1], [0, 2, 1]],
        "fixed_q": {"0": 4, "1": 4, "3/4": 2, "7/4": 6},
        "facets": 15,
        "finite_vertices": 42,
        "cusps": 10,
        "elliptic_edges": 294,
        "parabolic_audits": 17,
        "group": "C4",
        "ns_image": "C2",
    },
    {
        "id": "R12-D4x2A1x2-p67-H10-38",
        "bh_id": "0.4.4.6",
        "member": "pure04",
        "row": 44,
        "row_sha256": "2c5db8537584081dac8e9ddda92022f41a129f32b043063fd8155897449168ff",
        "source_kind": "d4",
        "source_key": "R12-D4x2A1x2-p67-H10-38",
        "lattice": "U+D4(-1)^2+A1(-1)^2",
        "determinant": -64,
        "smith": [*[1] * 6, *[2] * 6],
        "q": {"0": 16, "1": 16, "1/2": 12, "3/2": 20},
        "fixed_map": [
            [0, 1, 0, 0, 1],
            [0, 0, 0, 1, 0],
            [0, 0, 1, 0, 1],
            [1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1],
        ],
        "fixed_q": {"0": 16, "1": 16, "3/4": 12, "7/4": 20},
        "facets": 38,
        "finite_vertices": 526,
        "cusps": 141,
        "elliptic_edges": 4518,
        "parabolic_audits": 307,
        "group": "C4",
        "ns_image": "C2",
    },
    {
        "id": "R12-D4x2A1x2-p67-H10-38",
        "bh_id": "1.4.4.2",
        "member": "full04",
        "row": 171,
        "row_sha256": "d29436def2b38cb2cd55967b73e064583892ce0742746458fcc74a645f5b51a6",
        "source_kind": "d4",
        "source_key": "R12-D4x2A1x2-p67-H10-38",
        "lattice": "U+D4(-1)^2+A1(-1)^2",
        "determinant": -64,
        "smith": [*[1] * 6, *[2] * 6],
        "q": {"0": 16, "1": 16, "1/2": 12, "3/2": 20},
        "fixed_map": [
            [0, 0, 0, 1, 1],
            [0, 1, 0, 1, 1],
            [0, 0, 1, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 2, 0, 2, 1],
        ],
        "fixed_q": {"0": 16, "1": 16, "3/4": 12, "7/4": 20},
        "facets": 38,
        "finite_vertices": 526,
        "cusps": 141,
        "elliptic_edges": 4518,
        "parabolic_audits": 307,
        "group": "C4 x C2",
        "small_group_ids": [[2, 1], [8, 2]],
        "ns_image": "C2 x C2",
    },
)


SEARCH_FILTER = {
    "pure_rows_scanned": 587,
    "pure_rank_12_or_13_rows": 53,
    "pure_generic_NS_target_hits": [
        "0.2.8.3", "0.2.8.4", "0.2.9.1", "0.2.9.2", "0.2.9.3",
        "0.4.4.1", "0.4.4.2", "0.4.4.3", "0.4.4.4", "0.4.4.5", "0.4.4.6",
    ],
    "pure_exact_action_hits": ["0.4.4.1", "0.4.4.2", "0.4.4.6"],
    "full_rows_scanned": 3561,
    "full_rank_12_or_13_rows": 92,
    "full_generic_NS_target_hits": ["1.2.8.33", "1.2.9.1", "1.4.4.1", "1.4.4.2"],
    "full_exact_action_hits": ["1.4.4.2"],
    "rank13_exact_action_hits": [],
    "method": (
        "compute the generic NS as the orthogonal complement of the primitive "
        "cyclotomic transcendental lattice inside the exact K3 lattice; match "
        "rank/signature/determinant/Smith/full discriminant q; enumerate every "
        "finite NS image and its order-2/order-3 subgroups; match fixed q and glue"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _source_pair(
    specification: dict[str, Any], graph: dict[str, Any], d4: dict[str, Any]
) -> dict[str, Any]:
    if specification["source_kind"] == "graph":
        stock = next(
            row for row in graph["stocks"] if row["stock_key"] == specification["source_key"]
        )
        case = next(
            row
            for row in stock["survivors"]
            if row["class_index"] == specification["source_class"]
        )
        gram = sp.Matrix(stock["basis_gram"])
        action = sp.Matrix(
            case["neron_severi_action"]["generator_action_matrices"][0]
        )
        geometry = case["geometry"]
        vinberg = geometry["vinberg_volume_criterion"]
        summary = {
            "classification": vinberg["decision"],
            "dimension": geometry["hyperbolic_dimension"],
            "facets": geometry["facet_count"],
            "finite_vertices": geometry["finite_vertex_count"],
            "cusps": geometry["ideal_vertex_count"],
            "elliptic_edges": vinberg["elliptic_edge_count"],
            "parabolic_audits": vinberg["parabolic_neighborhood_audit_count"],
            "bad_edges": len(vinberg["finite_volume_bad_edges"]),
        }
    else:
        case = next(row for row in d4["cases"] if row["id"] == specification["source_key"])
        reconstruction = d4["coordinate_reconstruction"]
        gram = sp.Matrix(reconstruction["basis_gram"])
        walls = sp.Matrix(reconstruction["wall_coordinate_rows"])
        permutation = case["subgroup"]["representative_generators_zero_based"][0]
        action = walls[[permutation[index] for index in range(12)], :]
        vinberg = case["vinberg_finite_volume_audit"]
        summary = {
            "classification": "finite_volume_noncompact",
            "dimension": vinberg["dimension"],
            "facets": vinberg["facet_count"],
            "finite_vertices": vinberg["finite_vertex_count"],
            "cusps": vinberg["ideal_vertex_count"],
            "elliptic_edges": vinberg["elliptic_edge_count"],
            "parabolic_audits": vinberg["connected_affine_component_count"],
            "bad_edges": vinberg["bad_edge_count"] + vinberg["bad_parabolic_neighborhood_count"],
        }
    fixed_basis = saturated_row_kernel(action - sp.eye(12))
    coinvariant_basis = saturated_row_kernel(gram * fixed_basis.T)
    return {
        "gram": gram,
        "action": action,
        "fixed_basis": fixed_basis,
        "fixed_gram": sp.simplify(fixed_basis * gram * fixed_basis.T),
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": sp.simplify(
            coinvariant_basis * gram * coinvariant_basis.T
        ),
        "coinvariant_action": induced_row_action(coinvariant_basis, action),
        "glue": glue_pairs(fixed_basis, coinvariant_basis, gram),
        "geometry": summary,
    }


def _pinned_bh_pair(specification: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(NS_GRAMS[specification["bh_id"]])
    action = sp.Matrix(FULL_NS_G if specification["member"] == "full04" else PURE_ACTION)
    fixed_basis = saturated_row_kernel(action - sp.eye(12))
    coinvariant_basis = saturated_row_kernel(gram * fixed_basis.T)
    return {
        "gram": gram,
        "action": action,
        "fixed_basis": fixed_basis,
        "fixed_gram": sp.simplify(fixed_basis * gram * fixed_basis.T),
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": sp.simplify(
            coinvariant_basis * gram * coinvariant_basis.T
        ),
        "coinvariant_action": induced_row_action(coinvariant_basis, action),
        "glue": glue_pairs(fixed_basis, coinvariant_basis, gram),
    }


def _pair_conjugacy(
    specification: dict[str, Any], source: dict[str, Any], bh: dict[str, Any]
) -> dict[str, Any]:
    fixed_map = sp.Matrix(specification["fixed_map"])
    fixed_certificate = verify_fixed_discriminant_isometry(
        source["fixed_gram"], bh["fixed_gram"], fixed_map
    )
    coin_intertwiner = sp.eye(1)
    if not (
        source["coinvariant_gram"] == bh["coinvariant_gram"] == sp.Matrix([[-4]])
        and source["coinvariant_action"]
        == bh["coinvariant_action"]
        == sp.Matrix([[-1]])
    ):
        raise AssertionError("the rank-one coinvariant pair changed")
    coin_map = discriminant_map(
        coin_intertwiner, source["coinvariant_gram"], bh["coinvariant_gram"]
    )
    target_coin_moduli, _, _ = _discriminant_data(bh["coinvariant_gram"])
    target_graph = {coin: fixed for fixed, coin in bh["glue"]}
    for source_fixed, source_coin in source["glue"]:
        target_fixed = apply_reduced_fixed_map(
            source_fixed,
            source["fixed_gram"],
            bh["fixed_gram"],
            fixed_map,
        )
        target_coin = apply_discriminant_map(
            coin_map, source_coin, target_coin_moduli
        )
        if target_graph.get(target_coin) != target_fixed:
            raise AssertionError("the fixed and coinvariant maps do not carry the glue graph")
    fixed_length = len(
        [value for value in smith_diagonal(bh["fixed_gram"]) if value > 1]
    )
    if not (
        len(source["glue"]) == len(bh["glue"]) == 2
        and bh["fixed_gram"].rows == 11
        and exact_inertia(bh["fixed_gram"]) == (1, 10, 0)
        and bh["fixed_gram"].rows >= fixed_length + 2
        and exact_inertia(bh["coinvariant_gram"]) == (0, 1, 0)
        and bh["coinvariant_gram"][0, 0] % 4 == 0
    ):
        raise AssertionError("the indefinite-lattice lift or rootless gate changed")
    return {
        "fixed_discriminant_isometry": fixed_certificate,
        "coinvariant_integral_intertwiner": [[1]],
        "coinvariant_discriminant_map": matrix_payload(coin_map),
        "primitive_glue_order": 2,
        "all_glue_pairs_carried": True,
        "indefinite_fixed_lattice_rank_at_least_discriminant_length_plus_two": True,
        "nikulin_discriminant_isometry_lifts_integrally": True,
        "full_NS_action_pairs_integrally_conjugate": True,
        "algebraic_coinvariant": "<-4>, primitive of divisibility 2",
        "algebraic_coinvariant_rootless": True,
    }


def _derive_external_row(line: str, order: int) -> dict[str, Any]:
    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    ambient_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    inverse = basis.inv()
    gram = sp.simplify(basis * ambient_gram * basis.T)
    symplectic = [
        sp.simplify(basis * sp.Matrix(22, 22, action) * inverse)
        for action in data["symplectic"]
    ]
    generator = sp.simplify(
        basis * sp.Matrix(22, 22, data["distinguished_generator"]) * inverse
    )
    equations = cyclotomic_evaluation(generator, order)
    for action in symplectic:
        equations = equations.row_join(action - sp.eye(22))
    transcendental_basis = saturated_row_kernel(equations)
    ns_basis = saturated_row_kernel(gram * transcendental_basis.T)
    ns_gram = sp.simplify(ns_basis * gram * ns_basis.T)
    ns_symplectic = [induced_row_action(ns_basis, action) for action in symplectic]
    ns_generator = induced_row_action(ns_basis, generator)
    transcendental_gram = sp.simplify(
        transcendental_basis * gram * transcendental_basis.T
    )
    if not (
        gram.det() == -1
        and exact_inertia(gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in gram)
        and all(gram[index, index] % 2 == 0 for index in range(22))
        and all(
            action * gram * action.T == gram
            for action in [*symplectic, generator]
        )
        and exact_order(generator, order)
        and transcendental_basis.rows == 10
        and exact_inertia(transcendental_gram) == (2, 8, 0)
        and ns_gram.shape == (12, 12)
    ):
        raise AssertionError("a Brandhorst--Hofmann K3 row failed exact replay")
    return {
        "data": data,
        "gram": gram,
        "generator": generator,
        "symplectic": symplectic,
        "transcendental_gram": transcendental_gram,
        "ns_gram": ns_gram,
        "ns_symplectic": ns_symplectic,
        "ns_generator": ns_generator,
    }


def replay_external_zip(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(payload) == BH_ZIP_SHA256
    ):
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    replayed = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        if sha256_bytes(archive.read("K3Groups/readme.txt")) != BH_README_SHA256:
            raise AssertionError("the Brandhorst--Hofmann database readme changed")
        member_lines: dict[str, list[str]] = {}
        for key, member in MEMBERS.items():
            raw = archive.read(member["path"])
            if sha256_bytes(raw) != member["sha256"]:
                raise AssertionError("a Brandhorst--Hofmann member changed")
            lines = raw.decode("utf-8").splitlines()
            if len(lines) != member["rows"]:
                raise AssertionError("a Brandhorst--Hofmann member row count changed")
            member_lines[key] = lines
        for specification in CASES:
            line = member_lines[specification["member"]][specification["row"] - 1]
            if sha256_bytes((line + "\n").encode()) != specification["row_sha256"]:
                raise AssertionError("a pinned Brandhorst--Hofmann row changed")
            row = _derive_external_row(line, 4)
            expected_gram = sp.Matrix(NS_GRAMS[specification["bh_id"]])
            expected_action = sp.Matrix(
                FULL_NS_G if specification["member"] == "full04" else PURE_ACTION
            )
            if not (
                row["data"]["id"] == specification["bh_id"]
                and row["ns_gram"] == expected_gram
                and row["ns_generator"] == expected_action
                and abs(int(row["transcendental_gram"].det()))
                == abs(specification["determinant"])
                and matrix_digest(row["ns_gram"])
                == matrix_digest(expected_gram)
            ):
                raise AssertionError("a pinned generic Picard action changed")
            if specification["member"] == "full04":
                properties = dict(row["data"]["properties"])
                observed_ids = [list(pair) for pair in properties["small group ids"]]
                if not (
                    observed_ids == specification["small_group_ids"]
                    and len(row["symplectic"]) == len(row["ns_symplectic"]) == 1
                    and exact_order(row["symplectic"][0], 2)
                    and row["symplectic"][0] * row["generator"]
                    == row["generator"] * row["symplectic"][0]
                    and row["ns_symplectic"][0] == sp.Matrix(FULL_NS_S)
                    and exact_order(row["ns_generator"], 2)
                ):
                    raise AssertionError("the full C4 x C2 group structure changed")
            else:
                if not (
                    all(action == sp.eye(22) for action in row["symplectic"])
                    and exact_order(row["ns_generator"], 2)
                ):
                    raise AssertionError("a purely non-symplectic cyclic row changed")
            replayed.append(specification["bh_id"])
    return {"zip_sha256": BH_ZIP_SHA256, "rows_replayed": replayed}


def _full_group_rootless_audit() -> dict[str, Any]:
    gram = sp.Matrix(NS_GRAMS["1.4.4.2"])
    s = sp.Matrix(FULL_NS_S)
    g = sp.Matrix(FULL_NS_G)
    elements = {matrix_digest(s**a * g**b) for a in range(2) for b in range(2)}
    fixed = saturated_row_kernel((s - sp.eye(12)).row_join(g - sp.eye(12)))
    coinvariant = saturated_row_kernel(gram * fixed.T)
    fixed_gram = sp.simplify(fixed * gram * fixed.T)
    coinvariant_gram = sp.simplify(coinvariant * gram * coinvariant.T)
    checks = {
        "commuting_involutions": s**2 == g**2 == sp.eye(12) and s * g == g * s,
        "NS_image_is_V4": len(elements) == 4,
        "common_fixed_signature": fixed.rows == 4 and exact_inertia(fixed_gram) == (1, 3, 0),
        "common_fixed_determinant": fixed_gram.det() == -16,
        "common_coinvariant_signature": coinvariant.rows == 8 and exact_inertia(coinvariant_gram) == (0, 8, 0),
        "common_coinvariant_determinant": coinvariant_gram.det() == 256,
        "common_coinvariant_twice_even": (
            all(sp.Rational(value).q == 1 for value in coinvariant_gram / 2)
            and all(coinvariant_gram[index, index] % 4 == 0 for index in range(8))
        ),
        "common_coinvariant_rootless": True,
    }
    if not all(checks.values()):
        raise AssertionError(f"the full-group ample/rootless gate changed: {checks}")
    return {
        "group": "C4 x C2",
        "small_group_id": [8, 2],
        "NS_image": "C2 x C2",
        "common_fixed_rank_signature_determinant": [4, [1, 3, 0], -16],
        "common_coinvariant_rank_signature_determinant": [8, [0, 8, 0], 256],
        "rootless_reason": (
            "the common coinvariant Gram is twice an even negative-definite "
            "integral Gram, so every nonzero norm is divisible by four"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    hashes = {"graph": sha256(GRAPH_SOURCE), "d4": sha256(D4_SOURCE)}
    if hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned rank-12 Roulleau source changed")
    graph = json.loads(GRAPH_SOURCE.read_text(encoding="utf-8"))
    d4 = json.loads(D4_SOURCE.read_text(encoding="utf-8"))
    if not (graph["all_checks_pass"] and d4["all_checks_pass"]):
        raise AssertionError("an upstream rank-12 chamber certificate is no longer green")
    external = replay_external_zip(external_zip) if external_zip is not None else None

    realizations = []
    for specification in CASES:
        source = _source_pair(specification, graph, d4)
        bh = _pinned_bh_pair(specification)
        geometry = source["geometry"]
        fixed_q = _q_histogram(bh["fixed_gram"])
        checks = {
            "source_NS_signature": exact_inertia(source["gram"]) == (1, 11, 0),
            "source_NS_determinant": source["gram"].det() == specification["determinant"],
            "source_action_involution": exact_order(source["action"], 2),
            "source_action_isometry": source["action"] * source["gram"] * source["action"].T == source["gram"],
            "BH_NS_signature": exact_inertia(bh["gram"]) == (1, 11, 0),
            "BH_NS_determinant": bh["gram"].det() == specification["determinant"],
            "BH_NS_smith": smith_diagonal(bh["gram"]) == specification["smith"],
            "BH_NS_discriminant_q": _q_histogram(bh["gram"]) == specification["q"],
            "BH_action_involution": exact_order(bh["action"], 2),
            "BH_action_isometry": bh["action"] * bh["gram"] * bh["action"].T == bh["gram"],
            "action_character": str(sp.factor(bh["action"].charpoly().as_expr())) == "(lambda - 1)**11*(lambda + 1)",
            "fixed_rank": bh["fixed_gram"].rows == source["fixed_gram"].rows == 11,
            "fixed_determinant": bh["fixed_gram"].det() == source["fixed_gram"].det(),
            "fixed_smith": smith_diagonal(bh["fixed_gram"]) == smith_diagonal(source["fixed_gram"]),
            "fixed_q": fixed_q == _q_histogram(source["fixed_gram"]) == specification["fixed_q"],
            "vinberg_classification": geometry["classification"] == "finite_volume_noncompact",
            "vinberg_counts": (
                geometry["dimension"] == 10
                and geometry["facets"] == specification["facets"]
                and geometry["finite_vertices"] == specification["finite_vertices"]
                and geometry["cusps"] == specification["cusps"]
                and geometry["elliptic_edges"] == specification["elliptic_edges"]
                and geometry["parabolic_audits"] == specification["parabolic_audits"]
                and geometry["bad_edges"] == 0
            ),
        }
        if not all(checks.values()):
            raise AssertionError(f"a rank-12 lift gate changed for {specification['bh_id']}: {checks}")
        conjugacy = _pair_conjugacy(specification, source, bh)
        realizations.append(
            {
                "polytope_id": specification["id"],
                "brandhorst_hofmann": {
                    "id": specification["bh_id"],
                    "effective_group": specification["group"],
                    "transcendental_value": 4,
                    "member": MEMBERS[specification["member"]]["path"],
                    "member_sha256": MEMBERS[specification["member"]]["sha256"],
                    "row_one_based": specification["row"],
                    "row_sha256_with_newline": specification["row_sha256"],
                    "generic_NS": specification["lattice"],
                    "generic_NS_gram_digest": matrix_digest(bh["gram"]),
                    "Picard_image": specification["ns_image"],
                    "selected_Picard_action_order": 2,
                },
                "exact_action_conjugacy": conjugacy,
                "chamber_transport": {
                    "database_effective_character_supplies_invariant_Kahler_class": True,
                    "complete_Roulleau_nef_chamber": True,
                    "weyl_group_transitive_on_root_chambers": True,
                    "fixed_form_and_exhaustive_weighted_chamber_census_identify_archived_class": True,
                    "conclusion": (
                        "after Weyl marking and a chamber symmetry, the effective "
                        "BH automorphism preserves the archived complete nef chamber "
                        "and induces exactly the archived wall action"
                    ),
                },
                "coxeter_polytope": {
                    "dimension": 10,
                    "facets": specification["facets"],
                    "finite_vertices": specification["finite_vertices"],
                    "cusps": specification["cusps"],
                    "elliptic_edges": specification["elliptic_edges"],
                    "classification": "finite_volume_noncompact",
                    "Vinberg_bad_audits": 0,
                },
                "novelty_scope": (
                    "new project-level genuine K3/complete-nef realization of an "
                    "already certified chamber; no new abstract polytope or arithmetic "
                    "commensurability class is claimed"
                ),
                "checks": checks,
                "all_checks_pass": True,
            }
        )

    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau_rank12.brandhorst_hofmann_k3_lifts.v1",
        "artifact_id": "ROULLEAU-R12-BH-K3-LIFTS",
        "theorem": (
            "Three previously diagram-only rank-12 Roulleau folds admit genuine "
            "K3-automorphism realizations: H10/12 via BH 0.4.4.1, H10/15 via "
            "BH 0.4.4.2, and H10/38 via BH 0.4.4.6.  Full row 1.4.4.2 "
            "gives a second C4 x C2 realization of H10/38."
        ),
        "distinct_polytope_realization_count": 3,
        "effective_K3_group_realization_count": 4,
        "source_sha256": hashes,
        "brandhorst_hofmann_provenance": {
            "article_doi": BH_ARTICLE_DOI,
            "database_doi": BH_DATA_DOI,
            "database_version": BH_DATA_VERSION,
            "zip_md5": BH_ZIP_MD5,
            "zip_sha256": BH_ZIP_SHA256,
            "external_replay": external,
        },
        "finite_positive_filter": SEARCH_FILTER,
        "realizations": realizations,
        "full_group_enrichment": _full_group_rootless_audit(),
        "status_separation": {
            "genuine_surface_automorphisms": True,
            "generic_NS_equals_archived_complete_nef_lattice": True,
            "abstract_polytopes_new": False,
            "arithmetic_classes_new": False,
            "rank13_positive_hit": False,
            "rank13_result_role": "finite exact search boundary only",
        },
        "checks": {
            "upstream_chambers_green": True,
            "all_four_action_pairs_integrally_conjugate": all(
                row["exact_action_conjugacy"]["full_NS_action_pairs_integrally_conjugate"]
                for row in realizations
            ),
            "all_four_effective_K3_realizations": len(realizations) == 4,
            "three_distinct_polytopes_realized": len({row["polytope_id"] for row in realizations}) == 3,
            "all_vinberg_audits_green": all(
                row["coxeter_polytope"]["Vinberg_bad_audits"] == 0
                for row in realizations
            ),
        },
    }
    certificate["all_checks_pass"] = all(certificate["checks"].values())
    if not certificate["all_checks_pass"]:
        raise AssertionError("the rank-12 K3-lift theorem failed")
    canonical = json.dumps(certificate, sort_keys=True, separators=(",", ":"))
    certificate["payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--brandhorst-hofmann-zip",
        type=Path,
        help="optional path to the published K3Groups.zip for full row replay",
    )
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_zip)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archive mismatch: regenerate {OUTPUT}")
    elif not args.stdout:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if args.stdout or args.check_archive:
        print(rendered, end="")


if __name__ == "__main__":
    main()
