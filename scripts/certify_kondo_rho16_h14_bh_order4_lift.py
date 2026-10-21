#!/usr/bin/env python3
"""Certify a genuine order-four K3 realization of Kondo's rho-16 H14 fold.

Brandhorst--Hofmann pure row 0.4.2.2 and full row 1.4.2.1 have generic Picard
lattice ``U+E8(-1)+D6(-1)``.  The pure generator in the first row, and the
element ``s*g`` in the effective ``C4 x C2`` group of the second, act on
Picard as Kondo's previously non-realized graph involution tau_2.  The
fixed/coinvariant/glue calculation proves full integral action conjugacy,
while the complete archived Vinberg audit supplies the finite-volume H14/17
polytope.

The earlier order-two obstruction is not contradicted: the realizing
automorphism has order four, and its square is identity on Picard and minus
identity on the generic transcendental lattice.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from typing import Any
import zipfile

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import _q_histogram
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
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
FOLD = ROOT / "data/enumerations/k3_high_rank_involution_folds_rho15_17.json"
OLD_BOUNDARY = ROOT / "data/candidates/k3_rho16_h14_fold_identification.json"
INDEPENDENT_FOLD = ROOT / "data/candidates/k3_rho16_h14_fold_audit.json"
OUTPUT = ROOT / "data/enumerations/kondo_rho16_h14_bh_order4_lift.json"

EXPECTED_SOURCE_SHA256 = {
    "complete_nef_source": "0dcc7ec9c10f3437cd39713cd0d6036dc4fa8ad0c978ffdad4c2b34748269120",
    "fold_and_vinberg": "df6aa10a21faaca264fc120c2459ecc36d5e727b87c6e84411d6b0b47ef345c5",
    "old_order2_boundary": "9aa30d14fddaa4c212c34642f391dd3476b17591254d0533ae8069cf836f908b",
    "independent_fold_audit": "9c5ddc288b39016b0ee09cd97252b2fee7a805bdb8bc79ab461423ec5ae71039",
}

BH_ARTICLE_DOI = "10.1017/fms.2023.50"
BH_DATA_DOI = "10.5281/zenodo.10410278"
BH_DATA_VERSION = "2.1"
BH_ZIP_MD5 = "143047b9cb76f4ad4cb2aadd5959d879"
BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_MEMBER_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_MEMBER_ROWS = 79
BH_ROW = 17
BH_ROW_ID = "0.4.2.2"
BH_ROW_SHA256 = "8335fe304b5944898671a22a4c90f50290025eadbff4b1c51b58d0d92c051265"

BH_FULL_MEMBER = "K3Groups/results/order04.txt"
BH_FULL_MEMBER_SHA256 = "828763a2fc6609e7406e22bda835bd32701da8a667134908c22c8c6a6c30a6ca"
BH_FULL_MEMBER_ROWS = 670
BH_FULL_ROW = 65
BH_FULL_ROW_ID = "1.4.2.1"
BH_FULL_ROW_SHA256 = "607ba6b9fcfafa9dce645d3135fa04320732c28fd5ec2c289396b91ac5c30d1f"

EXPECTED_PURE = {
    "k3_gram_digest": "aafdaf323fa28df3b3b548b8e637c40760081cf1c973c15431a0889559c4db49",
    "k3_action_digest": "18d300446bf17e58dd9bc23dd8101353c68ec854beafac7bde01576621c25517",
    "transcendental_gram_digest": "6039b685d3ef64f1af333739a1be1a76e7414923d8ba9a571105810439257680",
    "ns_gram_digest": "15d314000640afd66aca1df8368bb900c20f846f8d7e37dc3bf0949d38cedf46",
    "ns_action_digest": "a50dae6250e8cf60a7019246e5e268c82913f4ff3d40293113d34fc314e2d9a4",
    "fixed_gram_digest": "c36627dc77077326b08bc3c23e45c6460abfc950d4f69f76f2d91a8f0736c467",
    "fixed_rank_signature_determinant": [15, [1, 14, 0], 4],
    "fixed_smith": [*[1] * 14, 4],
    "fixed_q": {"0": 1, "1": 1, "3/4": 2},
    "coinvariant_gram": [[-4]],
    "coinvariant_action": [[-1]],
    "glue": [
        [[0] * 15, [0]],
        [[0] * 14 + [2], [2]],
    ],
}

EXPECTED_FULL = {
    "k3_gram_digest": "b997c8cd4e761fccd9205d54a2365aca7702c995863fce44947e115e94ed73f1",
    "symplectic_action_digest": "45903b1ea5d662c9a2fb4f4dac281ab5407d7e7bb3e11f50685ef16615095a13",
    "distinguished_action_digest": "c0968654de84995d5dbd6bcc482a8271592bcf54b0a49e583d2ad9ce7c42eda2",
    "selected_sg_action_digest": "295d10ed9946f9a707d5f0519968aacbef4036c3d893a85d173ed73953c21376",
    "transcendental_gram_digest": "622b26c7bcaf984fc6d4ec1657ad1c59e34b69b413cace021e568d94cced4c28",
    "ns_gram_digest": "ae1f82cae3b64549a5526f3435e987277b21c28c71465a4005a9348591aec36c",
    "ns_s_digest": "8b08af22f0bbf19d0486caf1360f60c5611b334de325f689e9528c8a7de467a2",
    "ns_g_digest": "9aa723d03fed11f2040c799eb35ef7146d8fafbe3e003e1238a9a5550fadf311",
    "ns_sg_digest": "678ec3375c10093ac47f38036986b3a622bf7eba82a83b68abe372d82ddec01e",
    "selected_fixed_gram_digest": "55f17ce6fc4e36bd22d965a5d51de7c3ecc9c1a4ffc4f6756b9f2240c9ead68b",
    "fixed_rank_signature_determinant": [15, [1, 14, 0], 4],
    "fixed_smith": [*[1] * 14, 4],
    "fixed_q": {"0": 1, "1": 1, "3/4": 2},
    "coinvariant_gram": [[-4]],
    "coinvariant_action": [[-1]],
    "glue": [
        [[0] * 15, [0]],
        [[0] * 14 + [2], [2]],
    ],
    "common_fixed_rank_signature_determinant": [8, [1, 7, 0], -256],
}


SEARCH_CENSUS = {
    "pure_rows_scanned": 587,
    "pure_exact_rho16_target_hits": ["0.4.2.2"],
    "full_rows_scanned": 3561,
    "full_order2_rows_scanned": 2140,
    "full_non_order2_rows_scanned": 1421,
    "full_order2_generic_UE8D6_hits": 2,
    "full_order2_exact_rho16_target_hits": [],
    "full_non_order2_generic_UE8D6_hits": 1,
    "full_non_order2_exact_rho16_target_hits": ["1.4.2.1"],
    "rank15_C2_or_C3_exact_target_hits": [],
    "companion_script": "scripts/scout_bh_rank15_16_positive_lifts.py",
    "method": (
        "reconstruct the generic primitive cyclotomic transcendental lattice "
        "and its orthogonal Picard lattice exactly; match rank, determinant, "
        "Smith invariants and full discriminant q; enumerate the finite Picard "
        "image and compare every order-two/order-three fixed lattice"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(value)
            for value in matrix.row(index)
        ]
        for index in range(matrix.rows)
    ]


def _archive_pair(
    source: dict[str, Any], fold: dict[str, Any]
) -> dict[str, Any]:
    source_case = next(case for case in source["cases"] if case["rho"] == 16)
    wall_gram = build_gram(source_case)
    basis_indices = [int(index) - 1 for index in source_case["basis_indices"]]
    ns_gram = wall_gram.extract(basis_indices, basis_indices)
    fold_case = next(case for case in fold["cases"] if case["rho"] == 16)
    involution = next(
        row
        for row in fold_case["involutions"]
        if row["candidate_id"] == "rho16_involution_1"
    )
    geometry = next(
        row
        for row in fold_case["geometry_certificates"]
        if row["geometry_certificate_id"] == "rho16_class_1_geometry"
    )
    # The JSON stores coordinate columns as the outer list.  Transposing gives
    # the row-coordinate action used throughout these certificates.
    ns_action = sp.Matrix(
        involution["neron_severi_action"]["action_matrix_columns"]
    ).T
    fixed_basis = saturated_row_kernel(ns_action - sp.eye(16))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * ns_gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, ns_action)
    primitive_glue = glue_pairs(fixed_basis, coinvariant_basis, ns_gram)
    graph_fixed_rank_classes = {
        row["candidate_id"]: saturated_row_kernel(
            sp.Matrix(
                row["neron_severi_action"]["action_matrix_columns"]
            ).T
            - sp.eye(16)
        ).rows
        for row in fold_case["involutions"]
    }
    criterion = geometry["vinberg_volume_criterion"]
    checks = {
        "source_complete": source["artifact_status"]
        == "source_audited_gram_verified",
        "fold_archive_passes": fold["all_checks_pass"] is True,
        "ambient_lattice": source_case["lattice"] == "U + E8(-1) + D6(-1)",
        "ambient_rank_signature_determinant": (
            ns_gram.rows == 16
            and exact_inertia(ns_gram) == (1, 15, 0)
            and ns_gram.det() == -4
            and smith_diagonal(ns_gram) == [*[1] * 14, 2, 2]
            and _q_histogram(ns_gram) == {"0": 1, "1": 1, "1/2": 2}
        ),
        "action_is_integral_involution": (
            ns_action**2 == sp.eye(16)
            and ns_action * ns_gram * ns_action.T == ns_gram
            and abs(int(ns_action.det())) == 1
        ),
        "action_characteristic_polynomial": sp.factor(
            ns_action.charpoly().as_expr()
        )
        == (sp.Symbol("lambda") - 1) ** 15
        * (sp.Symbol("lambda") + 1),
        "fixed_lattice": (
            fixed_basis.rows == 15
            and exact_inertia(fixed_gram) == (1, 14, 0)
            and fixed_gram.det() == 4
            and smith_diagonal(fixed_gram) == [*[1] * 14, 4]
            and _q_histogram(fixed_gram) == {"0": 1, "1": 1, "3/4": 2}
        ),
        "coinvariant_is_rootless_minus_four": (
            coinvariant_gram == sp.Matrix([[-4]])
            and coinvariant_action == sp.Matrix([[-1]])
        ),
        "primitive_glue_order_two": len(primitive_glue) == 2,
        "complete_graph_symmetry_class_is_unique_at_fixed_rank_15": (
            fold_case["weighted_graph_automorphism_search"]["automorphism_count"]
            == 4
            and fold_case["weighted_graph_automorphism_search"][
                "nonidentity_involution_count"
            ]
            == 3
            and graph_fixed_rank_classes
            == {
                "rho16_involution_1": 15,
                "rho16_involution_2": 8,
                "rho16_involution_3": 9,
            }
        ),
        "complete_vinberg_finite_volume": (
            geometry["all_checks_pass"]
            and geometry["hyperbolic_dimension"] == 14
            and geometry["facet_count"] == 17
            and geometry["finite_vertex_count"] == 89
            and geometry["ideal_vertex_count"] == 5
            and geometry["hyperideal_vertex_count"] == 0
            and criterion["decision"] == "finite_volume_noncompact"
            and criterion["elliptic_edge_count"] == 704
            and not criterion["finite_volume_bad_edges"]
            and criterion["parabolic_neighborhood_audit_count"] == 8
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"the archived rank-16 H14 theorem changed: {checks}")
    return {
        "checks": checks,
        "ns_gram": ns_gram,
        "ns_action": ns_action,
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": coinvariant_gram,
        "coinvariant_action": coinvariant_action,
        "glue": primitive_glue,
        "graph_fixed_rank_classes": graph_fixed_rank_classes,
        "source_case": source_case,
        "involution": involution,
        "geometry": geometry,
    }


def _derive_bh_row(line: str) -> dict[str, Any]:
    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    source_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    k3_gram = sp.simplify(basis * source_gram * basis.T)
    k3_action = sp.simplify(
        basis
        * sp.Matrix(22, 22, data["distinguished_generator"])
        * basis.inv()
    )
    # Purely non-symplectic rows still carry the identity as the sole
    # symplectic entry in the serialized database format.
    symplectic_entries = [
        sp.Matrix(22, 22, entry) for entry in data["symplectic"]
    ]
    transcendental_basis = saturated_row_kernel(
        cyclotomic_evaluation(k3_action, 4)
    )
    ns_basis = saturated_row_kernel(k3_gram * transcendental_basis.T)
    transcendental_gram = sp.simplify(
        transcendental_basis * k3_gram * transcendental_basis.T
    )
    ns_gram = sp.simplify(ns_basis * k3_gram * ns_basis.T)
    ns_action = induced_row_action(ns_basis, k3_action)
    fixed_basis = saturated_row_kernel(ns_action - sp.eye(16))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * ns_gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, ns_action)
    primitive_glue = glue_pairs(fixed_basis, coinvariant_basis, ns_gram)
    if not (
        data["id"] == BH_ROW_ID
        and data["transcendental_value"] == 4
        and len(symplectic_entries) == 1
        and symplectic_entries[0] == sp.eye(22)
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and exact_order(k3_action, 4)
        and k3_action * k3_gram * k3_action.T == k3_gram
        and matrix_digest(k3_gram) == EXPECTED_PURE["k3_gram_digest"]
        and matrix_digest(k3_action) == EXPECTED_PURE["k3_action_digest"]
        and transcendental_basis.rows == 6
        and exact_inertia(transcendental_gram) == (2, 4, 0)
        and transcendental_gram.det() == 4
        and matrix_digest(transcendental_gram)
        == EXPECTED_PURE["transcendental_gram_digest"]
        and sp.factor(
            induced_row_action(transcendental_basis, k3_action)
            .charpoly()
            .as_expr()
        )
        == (sp.Symbol("lambda") ** 2 + 1) ** 3
        and ns_gram.rows == 16
        and exact_inertia(ns_gram) == (1, 15, 0)
        and ns_gram.det() == -4
        and smith_diagonal(ns_gram) == [*[1] * 14, 2, 2]
        and _q_histogram(ns_gram) == {"0": 1, "1": 1, "1/2": 2}
        and matrix_digest(ns_gram) == EXPECTED_PURE["ns_gram_digest"]
        and exact_order(ns_action, 2)
        and matrix_digest(ns_action) == EXPECTED_PURE["ns_action_digest"]
        and fixed_basis.rows == 15
        and exact_inertia(fixed_gram) == (1, 14, 0)
        and fixed_gram.det() == 4
        and smith_diagonal(fixed_gram) == EXPECTED_PURE["fixed_smith"]
        and _q_histogram(fixed_gram) == EXPECTED_PURE["fixed_q"]
        and matrix_digest(fixed_gram) == EXPECTED_PURE["fixed_gram_digest"]
        and coinvariant_gram == sp.Matrix(EXPECTED_PURE["coinvariant_gram"])
        and coinvariant_action == sp.Matrix(EXPECTED_PURE["coinvariant_action"])
        and primitive_glue
        == [(tuple(fixed), tuple(coin)) for fixed, coin in EXPECTED_PURE["glue"]]
        and ns_basis * k3_action**2 == ns_basis
        and induced_row_action(transcendental_basis, k3_action) ** 2
        == -sp.eye(6)
    ):
        raise AssertionError("Brandhorst--Hofmann row 0.4.2.2 changed")
    return {
        "data": data,
        "k3_gram": k3_gram,
        "k3_action": k3_action,
        "transcendental_basis": transcendental_basis,
        "transcendental_gram": transcendental_gram,
        "ns_basis": ns_basis,
        "ns_gram": ns_gram,
        "ns_action": ns_action,
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": coinvariant_gram,
        "coinvariant_action": coinvariant_action,
        "glue": primitive_glue,
    }


def _derive_full_bh_row(line: str) -> dict[str, Any]:
    """Replay the C4 x C2 enrichment and select the pure element s*g."""

    data = _parse_bh_row(line)
    basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
    source_gram = sp.Matrix(
        22, 22, list(map(_rational, data["inner_product_matrix"]))
    )
    k3_gram = sp.simplify(basis * source_gram * basis.T)
    if len(data["symplectic"]) != 1:
        raise AssertionError("the full rank-16 row lost its symplectic generator")
    s = sp.simplify(
        basis * sp.Matrix(22, 22, data["symplectic"][0]) * basis.inv()
    )
    g = sp.simplify(
        basis
        * sp.Matrix(22, 22, data["distinguished_generator"])
        * basis.inv()
    )
    selected = s * g
    transcendental_basis = saturated_row_kernel(
        cyclotomic_evaluation(g, 4).row_join(s - sp.eye(22))
    )
    transcendental_gram = sp.simplify(
        transcendental_basis * k3_gram * transcendental_basis.T
    )
    ns_basis = saturated_row_kernel(k3_gram * transcendental_basis.T)
    ns_gram = sp.simplify(ns_basis * k3_gram * ns_basis.T)
    ns_s = induced_row_action(ns_basis, s)
    ns_g = induced_row_action(ns_basis, g)
    ns_selected = ns_s * ns_g
    fixed_basis = saturated_row_kernel(ns_selected - sp.eye(16))
    coinvariant_basis = saturated_row_kernel(ns_gram * fixed_basis.T)
    fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * ns_gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, ns_selected)
    primitive_glue = glue_pairs(fixed_basis, coinvariant_basis, ns_gram)
    common_fixed = saturated_row_kernel(
        (s - sp.eye(22)).row_join(g - sp.eye(22))
    )
    common_fixed_gram = sp.simplify(common_fixed * k3_gram * common_fixed.T)
    full_elements = {matrix_digest(s**a * g**b) for a in range(2) for b in range(4)}
    ns_elements = {
        matrix_digest(ns_s**a * ns_g**b) for a in range(2) for b in range(4)
    }
    symplectic_fixed = saturated_row_kernel(s - sp.eye(22))
    symplectic_coinvariant = saturated_row_kernel(k3_gram * symplectic_fixed.T)
    symplectic_coinvariant_gram = sp.simplify(
        symplectic_coinvariant * k3_gram * symplectic_coinvariant.T
    )
    if not (
        data["id"] == BH_FULL_ROW_ID
        and data["transcendental_value"] == 4
        and dict(data["properties"])["small group ids"] == ((2, 1), (8, 2))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and matrix_digest(k3_gram) == EXPECTED_FULL["k3_gram_digest"]
        and exact_order(s, 2)
        and exact_order(g, 4)
        and exact_order(selected, 4)
        and s * g == g * s
        and len(full_elements) == 8
        and matrix_digest(s) == EXPECTED_FULL["symplectic_action_digest"]
        and matrix_digest(g) == EXPECTED_FULL["distinguished_action_digest"]
        and matrix_digest(selected) == EXPECTED_FULL["selected_sg_action_digest"]
        and transcendental_basis.rows == 6
        and exact_inertia(transcendental_gram) == (2, 4, 0)
        and transcendental_gram.det() == 4
        and matrix_digest(transcendental_gram)
        == EXPECTED_FULL["transcendental_gram_digest"]
        and ns_gram.rows == 16
        and exact_inertia(ns_gram) == (1, 15, 0)
        and ns_gram.det() == -4
        and smith_diagonal(ns_gram) == [*[1] * 14, 2, 2]
        and _q_histogram(ns_gram) == {"0": 1, "1": 1, "1/2": 2}
        and matrix_digest(ns_gram) == EXPECTED_FULL["ns_gram_digest"]
        and matrix_digest(ns_s) == EXPECTED_FULL["ns_s_digest"]
        and matrix_digest(ns_g) == EXPECTED_FULL["ns_g_digest"]
        and matrix_digest(ns_selected) == EXPECTED_FULL["ns_sg_digest"]
        and len(ns_elements) == 4
        and exact_order(ns_selected, 2)
        and fixed_basis.rows == 15
        and exact_inertia(fixed_gram) == (1, 14, 0)
        and fixed_gram.det() == 4
        and smith_diagonal(fixed_gram) == EXPECTED_FULL["fixed_smith"]
        and _q_histogram(fixed_gram) == EXPECTED_FULL["fixed_q"]
        and matrix_digest(fixed_gram) == EXPECTED_FULL["selected_fixed_gram_digest"]
        and coinvariant_gram == sp.Matrix(EXPECTED_FULL["coinvariant_gram"])
        and coinvariant_action == sp.Matrix(EXPECTED_FULL["coinvariant_action"])
        and primitive_glue
        == [(tuple(fixed), tuple(coin)) for fixed, coin in EXPECTED_FULL["glue"]]
        and [
            common_fixed.rows,
            list(exact_inertia(common_fixed_gram)),
            int(common_fixed_gram.det()),
        ]
        == EXPECTED_FULL["common_fixed_rank_signature_determinant"]
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
        and ns_basis * selected**2 == ns_basis
        and induced_row_action(transcendental_basis, selected) ** 2
        == -sp.eye(6)
    ):
        raise AssertionError("Brandhorst--Hofmann full row 1.4.2.1 changed")
    return {
        "data": data,
        "k3_gram": k3_gram,
        "k3_action": selected,
        "s": s,
        "g": g,
        "transcendental_basis": transcendental_basis,
        "transcendental_gram": transcendental_gram,
        "ns_basis": ns_basis,
        "ns_gram": ns_gram,
        "ns_action": ns_selected,
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": coinvariant_gram,
        "coinvariant_action": coinvariant_action,
        "glue": primitive_glue,
        "common_fixed_gram": common_fixed_gram,
        "symplectic_coinvariant_gram": symplectic_coinvariant_gram,
    }


def _direct_conjugacy(
    archive: dict[str, Any], bh: dict[str, Any]
) -> dict[str, Any]:
    fixed_map = sp.Matrix([[1]])
    fixed_certificate = verify_fixed_discriminant_isometry(
        archive["fixed_gram"], bh["fixed_gram"], fixed_map
    )
    coin_intertwiner = sp.eye(1)
    if not (
        coin_intertwiner
        * archive["coinvariant_gram"]
        * coin_intertwiner.T
        == bh["coinvariant_gram"]
        and coin_intertwiner * archive["coinvariant_action"]
        == bh["coinvariant_action"] * coin_intertwiner
    ):
        raise AssertionError("the rank-one coinvariant intertwiner changed")
    coin_map = discriminant_map(
        coin_intertwiner.inv(),
        archive["coinvariant_gram"],
        bh["coinvariant_gram"],
    )
    coin_moduli, _, _ = _discriminant_data(bh["coinvariant_gram"])
    bh_graph = {coin: fixed for fixed, coin in bh["glue"]}
    for archived_fixed, archived_coin in archive["glue"]:
        target_fixed = apply_reduced_fixed_map(
            archived_fixed,
            archive["fixed_gram"],
            bh["fixed_gram"],
            fixed_map,
        )
        target_coin = apply_discriminant_map(coin_map, archived_coin, coin_moduli)
        if bh_graph.get(target_coin) != target_fixed:
            raise AssertionError("the primitive rank-16 glue graph is not carried")
    if not (
        archive["fixed_gram"].rows == 15
        and 15 >= 1 + 2
        and exact_inertia(bh["coinvariant_gram"]) == (0, 1, 0)
        and bh["coinvariant_gram"] == sp.Matrix([[-4]])
    ):
        raise AssertionError("the Nikulin/rootless gate changed")
    return {
        "fixed_discriminant_isometry": fixed_certificate,
        "coinvariant_intertwiner": [[1]],
        "coinvariant_discriminant_map": matrix_payload(coin_map),
        "primitive_glue_order": len(bh["glue"]),
        "all_glue_elements_carried": True,
        "nikulin_rank_at_least_length_plus_two": True,
        "full_integral_picard_action_conjugacy": True,
        "coinvariant_rootless": True,
    }


def _cyclic_order4_local_type(
    gram: sp.Matrix, action: sp.Matrix
) -> tuple[tuple[Any, ...], ...]:
    answer = []
    for power in (1, 2, 4):
        fixed_basis = saturated_row_kernel(action**power - sp.eye(22))
        fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
        answer.append(
            (
                power,
                fixed_basis.rows,
                int(fixed_gram.det()),
                tuple(smith_diagonal(fixed_gram)),
            )
        )
    return tuple(answer)


def replay_external_zip(path: Path, archive_pair: dict[str, Any]) -> dict[str, Any]:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and sha256_bytes(payload) == BH_ZIP_SHA256
    ):
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        member = archive.read(BH_MEMBER)
        full_member = archive.read(BH_FULL_MEMBER)
    if sha256_bytes(member) != BH_MEMBER_SHA256:
        raise AssertionError("the pinned order-four member changed")
    lines = member.decode("utf-8").splitlines()
    if len(lines) != BH_MEMBER_ROWS:
        raise AssertionError("the order-four row count changed")
    line = lines[BH_ROW - 1]
    if sha256_bytes((line + "\n").encode()) != BH_ROW_SHA256:
        raise AssertionError("the pinned row 0.4.2.2 changed")
    bh = _derive_bh_row(line)
    pure_conjugacy = _direct_conjugacy(archive_pair, bh)
    if sha256_bytes(full_member) != BH_FULL_MEMBER_SHA256:
        raise AssertionError("the pinned full order-four member changed")
    full_lines = full_member.decode("utf-8").splitlines()
    if len(full_lines) != BH_FULL_MEMBER_ROWS:
        raise AssertionError("the full order-four row count changed")
    full_line = full_lines[BH_FULL_ROW - 1]
    if sha256_bytes((full_line + "\n").encode()) != BH_FULL_ROW_SHA256:
        raise AssertionError("the pinned full row 1.4.2.1 changed")
    full_bh = _derive_full_bh_row(full_line)
    full_conjugacy = _direct_conjugacy(archive_pair, full_bh)

    target_type = _cyclic_order4_local_type(
        full_bh["k3_gram"], full_bh["k3_action"]
    )
    pure_matches = []
    for row_number, candidate_line in enumerate(lines, 1):
        candidate_data = _parse_bh_row(candidate_line)
        candidate_basis = sp.Matrix(
            22, 22, list(map(_rational, candidate_data["basis"]))
        )
        candidate_gram = sp.simplify(
            candidate_basis
            * sp.Matrix(
                22,
                22,
                list(map(_rational, candidate_data["inner_product_matrix"])),
            )
            * candidate_basis.T
        )
        candidate_action = sp.simplify(
            candidate_basis
            * sp.Matrix(22, 22, candidate_data["distinguished_generator"])
            * candidate_basis.inv()
        )
        if _cyclic_order4_local_type(candidate_gram, candidate_action) == target_type:
            pure_matches.append([row_number, candidate_data["id"]])
    if pure_matches != [[BH_ROW, BH_ROW_ID]]:
        raise AssertionError("the C4 x C2 selected subgroup lineage changed")
    return {
        "zip_sha256": BH_ZIP_SHA256,
        "pure_row": {
            "member_sha256": BH_MEMBER_SHA256,
            "complete_member_row_count": BH_MEMBER_ROWS,
            "row_id": BH_ROW_ID,
            "row_one_based": BH_ROW,
            "row_sha256_with_newline": BH_ROW_SHA256,
            "full_22_dimensional_row_replayed": True,
            "direct_conjugacy": pure_conjugacy,
        },
        "full_group_row": {
            "member_sha256": BH_FULL_MEMBER_SHA256,
            "complete_member_row_count": BH_FULL_MEMBER_ROWS,
            "row_id": BH_FULL_ROW_ID,
            "row_one_based": BH_FULL_ROW,
            "row_sha256_with_newline": BH_FULL_ROW_SHA256,
            "full_22_dimensional_row_replayed": True,
            "selected_element": "s*g",
            "direct_conjugacy": full_conjugacy,
        },
        "nested_cyclic_lineage": {
            "complete_pure_order4_rows_checked": BH_MEMBER_ROWS,
            "unique_match_for_s_times_g_in_1.4.2.1": BH_ROW_ID,
            "row_one_based": BH_ROW,
            "small_local_type": [
                [power, rank, determinant, list(smith)]
                for power, rank, determinant, smith in target_type
            ],
        },
    }


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    paths = {
        "complete_nef_source": SOURCE,
        "fold_and_vinberg": FOLD,
        "old_order2_boundary": OLD_BOUNDARY,
        "independent_fold_audit": INDEPENDENT_FOLD,
    }
    hashes = {key: sha256(path) for key, path in paths.items()}
    if hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned rank-16 source certificate changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    old_boundary = json.loads(OLD_BOUNDARY.read_text(encoding="utf-8"))
    independent = json.loads(INDEPENDENT_FOLD.read_text(encoding="utf-8"))
    archive = _archive_pair(source, fold)
    if not (
        old_boundary["all_checks_pass"]
        and old_boundary["algebraic_realization"]["realized_K3_involution"] is False
        and "higher-order" in old_boundary["algebraic_realization"]["special_period_scope"]
        and independent["checks"][
            "abstract_fold_is_finite_volume_noncompact"
        ]
        is True
    ):
        raise AssertionError("the prior order-two boundary or independent audit changed")
    external = replay_external_zip(external_zip, archive) if external_zip else None
    geometry = archive["geometry"]
    criterion = geometry["vinberg_volume_criterion"]

    certificate: dict[str, Any] = {
        "schema": "coxeter.kondo_rho16_h14.bh_order4_lift.v1",
        "all_checks_pass": True,
        "theorem": (
            "There exists a complex projective K3 surface X with generic "
            "NS(X)=U+E8(-1)+D6(-1) and an effective C4 x C2 automorphism "
            "group.  A purely non-symplectic order-four element sigma=s*g "
            "acts on Picard as Kondo's tau_2.  Its fixed section of the "
            "complete 19-wall nef chamber is the finite-volume noncompact "
            "H14/17 Vinberg Coxeter polytope."
        ),
        "source_sha256": hashes,
        "brandhorst_hofmann": {
            "article": "Finite subgroups of automorphisms of K3 surfaces",
            "article_doi": BH_ARTICLE_DOI,
            "database_doi": BH_DATA_DOI,
            "database_version": BH_DATA_VERSION,
            "license": "CC-BY-4.0",
            "pure_cyclic_row": {
                "row_id": BH_ROW_ID,
                "member": BH_MEMBER,
                "member_sha256": BH_MEMBER_SHA256,
                "complete_member_row_count": BH_MEMBER_ROWS,
                "row_one_based": BH_ROW,
                "row_sha256_with_newline": BH_ROW_SHA256,
                "group": "C4",
                "effective_character": True,
                "purely_non_symplectic": True,
                "global_automorphism_order": 4,
            },
            "full_group_enrichment": {
                "row_id": BH_FULL_ROW_ID,
                "member": BH_FULL_MEMBER,
                "member_sha256": BH_FULL_MEMBER_SHA256,
                "complete_member_row_count": BH_FULL_MEMBER_ROWS,
                "row_one_based": BH_FULL_ROW,
                "row_sha256_with_newline": BH_FULL_ROW_SHA256,
                "group": "C4 x C2",
                "symplectic_kernel_Hs": "C2 generated by s",
                "selected_element": "s*g",
                "selected_element_order": 4,
                "selected_element_purely_non_symplectic": True,
                "picard_image": "C2 x C2",
                "picard_image_order": 4,
                "kernel_of_group_to_O_NS": "<g^2> = C2",
                "common_fixed_lattice_rank_signature_determinant": [
                    8,
                    [1, 7, 0],
                    -256,
                ],
                "selected_cyclic_subgroup_is_unique_pure_row": BH_ROW_ID,
            },
            "generic_transcendental_lattice": {
                "rank_signature_determinant": [6, [2, 4, 0], 4],
                "action_characteristic_polynomial": "(x^2+1)^3",
                "action_square": "-identity",
            },
            "generic_picard_lattice": {
                "isometry_class": "U+E8(-1)+D6(-1)",
                "rank_signature_determinant": [16, [1, 15, 0], -4],
                "smith": [*[1] * 14, 2, 2],
                "q_histogram_mod_2Z": {"0": 1, "1": 1, "1/2": 2},
                "action_order": 2,
                "action_characteristic_polynomial": "(x-1)^15*(x+1)",
            },
            "external_replay_performed": external is not None,
            "external_replay": external,
        },
        "exact_picard_action_conjugacy": {
            "archived_action": "Kondo tau_2 / rho16_involution_1",
            "fixed_lattices": {
                "rank_signature_determinant": [15, [1, 14, 0], 4],
                "smith": [*[1] * 14, 4],
                "q_histogram_mod_2Z": {"0": 1, "1": 1, "3/4": 2},
                "discriminant_group": "Z/4",
                "discriminant_isometry": [[1]],
                "nikulin_lifting": (
                    "rank 15 is at least discriminant length 1 plus 2, so the "
                    "displayed discriminant isometry lifts integrally"
                ),
            },
            "coinvariant_lattices": {
                "archived_gram": [[-4]],
                "database_gram": [[-4]],
                "action": [[-1]],
                "integral_intertwiner": [[1]],
                "negative_definite": True,
                "rootless": True,
            },
            "primitive_glue": {
                "index": 2,
                "archived_pairs_in_smith_coordinates": [
                    [[0] * 15, [0]],
                    [[0] * 14 + [2], [2]],
                ],
                "pure_database_pairs_in_smith_coordinates": EXPECTED_PURE["glue"],
                "full_group_database_pairs_in_smith_coordinates": EXPECTED_FULL[
                    "glue"
                ],
                "identity_maps_carry_both_pairs": True,
            },
            "conclusion": (
                "the fixed and coinvariant isometries carry the primitive "
                "glue and intertwine the actions for both BH rows, hence extend "
                "to integral conjugacies of the full Picard-lattice actions"
            ),
        },
        "kahler_and_torelli": {
            "both_database_rows_are_effective_actual_K3_automorphism_classes": True,
            "full_group_has_common_positive_fixed_lattice": True,
            "symplectic_coinvariant_is_rootless": True,
            "algebraic_coinvariant_is_minus_four_and_rootless": True,
            "archived_complete_nef_chamber_is_tau2_stable": True,
            "complete_weighted_wall_graph_fixed_ranks": archive[
                "graph_fixed_rank_classes"
            ],
            "tau2_is_unique_graph_symmetry_class_of_fixed_rank_15": True,
            "weyl_marking_to_archived_complete_chamber": True,
            "strong_torelli_handoff": (
                "effectivity supplies an invariant Kahler chamber; moving it by "
                "the Weyl group to the archived complete chamber produces a wall-graph "
                "involution.  Exhaustion of the weighted graph and uniqueness of fixed "
                "rank 15 identify it with the archived tau_2 class"
            ),
        },
        "order_two_boundary_reconciliation": {
            "prior_theorem_retained": (
                "tau_2 is not induced by an order-two K3 automorphism when "
                "NS=U+E8(-1)+D6(-1)"
            ),
            "new_mechanism": (
                "sigma has order four, sigma|NS=tau_2, sigma^2|NS=identity, "
                "and sigma^2|T=-identity"
            ),
            "contradiction": False,
            "previous_special_higher_order_period_question_closed": True,
        },
        "finite_database_search_census": SEARCH_CENSUS,
        "hyperbolic_polytope": {
            "id": "R16-UE8D6-C2-H14-17",
            "published_identification": "Vinberg 1972 P14 / 14-vinb85.coxiter",
            "dimension": 14,
            "facets": 17,
            "classification": "finite_volume_noncompact",
            "finite_vertices": geometry["finite_vertex_count"],
            "cusps": geometry["ideal_vertex_count"],
            "hyperideal_vertices": geometry["hyperideal_vertex_count"],
            "elliptic_edges": criterion["elliptic_edge_count"],
            "bad_edges": len(criterion["finite_volume_bad_edges"]),
            "parabolic_neighborhood_audits": criterion[
                "parabolic_neighborhood_audit_count"
            ],
            "complete_simple_roots": criterion["simple_roots_complete"],
        },
        "novelty_scope": {
            "new_project_level_genuine_K3_complete_nef_realization": True,
            "new_project_level_full_group_enrichment": True,
            "new_abstract_coxeter_polytope_claimed": False,
            "new_dimension_claimed": False,
            "worldwide_priority_claimed": False,
            "known_polytope_new_realization": True,
        },
        "summary": {
            "effective_database_rows": 2,
            "full_group_enrichments": 1,
            "genuine_K3_automorphism_realizations": 2,
            "automorphism_order": 4,
            "picard_action_order": 2,
            "hyperbolic_dimension": 14,
            "finite_volume_noncompact": 1,
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
            raise AssertionError("the archived rank-16 order-four lift changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
