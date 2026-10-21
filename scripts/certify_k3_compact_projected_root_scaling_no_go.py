#!/usr/bin/env python3
"""Exact replay for the compact projected-root scaling counterfamily."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_projected_root_scaling_no_go.md"
SOURCE_ARCHIVE = ROOT / "data/enumerations/golden_lanner_complete_chambers.json"
ARCHIVE = ROOT / "data/enumerations/k3_compact_projected_root_scaling_no_go.json"

SOURCE_ARCHIVE_SHA256 = (
    "eeee7aae8ebe339298d57ad24133820926f57501993f9a9f56cc3b11fd1c28b5"
)


Q = sp.Matrix(
    [
        [-2, -2, 0, 0, 1, 2],
        [-2, -4, 0, 0, 2, 3],
        [0, 0, -2, -2, 1, 2],
        [0, 0, -2, -4, 2, 3],
        [1, 2, 1, 2, -2, -2],
        [2, 3, 2, 3, -2, -4],
    ]
)
ROOT_COLUMNS = (0, 2, 4)
R = sp.zeros(6, 3)
for column, row in enumerate(ROOT_COLUMNS):
    R[row, column] = 1
H = R.T * Q * R
B = Q * R


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def ambient_gram(scale: int) -> sp.Matrix:
    if scale < 2:
        raise ValueError("scale must be at least two")
    return (scale * Q).row_join(B).col_join(B.T.row_join(H))


def direct_sum_gram(scale: int) -> sp.Matrix:
    if scale < 2:
        raise ValueError("scale must be at least two")
    return sp.diag(scale * Q, scale * (scale - 1) * H)


def inertia_symmetric(matrix: sp.Matrix) -> tuple[int, int, int]:
    eigenvalues = matrix.eigenvals()
    positive = negative = zero = 0
    for eigenvalue, multiplicity in eigenvalues.items():
        numerical = complex(sp.N(eigenvalue, 80))
        if abs(numerical.imag) > 1e-40:
            raise AssertionError("symmetric matrix acquired a nonreal eigenvalue")
        if numerical.real > 1e-30:
            positive += int(multiplicity)
        elif numerical.real < -1e-30:
            negative += int(multiplicity)
        else:
            zero += int(multiplicity)
    return positive, negative, zero


def sample_row(scale: int) -> dict[str, Any]:
    gram = ambient_gram(scale)
    return {
        "scale": scale,
        "ambient_determinant": int(gram.det()),
        "ambient_signature_positive_negative_zero": list(inertia_symmetric(gram)),
        "direct_sum_to_ambient_overlattice_index": scale**3,
        "trace_determinant": int((scale * Q).det()),
        "complement_determinant": int((scale * (scale - 1) * H).det()),
        "physical_root_gram": matrix_rows(H),
        "physical_root_square": -2,
    }


def source_triangle() -> dict[str, Any]:
    archived = json.loads(SOURCE_ARCHIVE.read_text(encoding="utf-8"))
    matches = [row for row in archived["examples"] if row["seed_id"] == "triangle_255"]
    if len(matches) != 1:
        raise AssertionError("frozen triangle_255 row is not unique")
    row = matches[0]
    trace = row["trace_k3_package"]
    if trace["gram_interleaved_basis_ei_phi_ei"] != matrix_rows(Q):
        raise AssertionError("frozen triangle trace Gram changed")
    if [[int(value) for value in line] for line in trace["simple_curve_gram"]] != matrix_rows(H):
        raise AssertionError("frozen simple curve Gram changed")
    return {
        "seed_id": row["seed_id"],
        "labels": row["input_coxeter_labels"],
        "compact": row["vinberg_compactness"]["compact"],
        "trace_rank": trace["rank"],
        "trace_determinant": trace["determinant"],
        "exact_eigen_nef_return": row["k3_eigen_nef_promotion"][
            "projectivized_eigen_nef_section_equals_closed_seed"
        ],
        "complete_visible_root_arrangement": trace["visible_k3_wall_shell_collapse"][
            "restricted_hyperplanes_equal"
        ],
    }


def build_certificate() -> dict[str, Any]:
    if sha256_file(SOURCE_ARCHIVE) != SOURCE_ARCHIVE_SHA256:
        raise AssertionError("golden source archive hash changed")
    triangle = source_triangle()
    if Q.det() != -16:
        raise AssertionError("source trace determinant changed")
    if H != sp.Matrix([[-2, 0, 1], [0, -2, 1], [1, 1, -2]]):
        raise AssertionError("source A3 Gram changed")
    if [int((-H)[:i, :i].det()) for i in range(1, 4)] != [2, 4, 4]:
        raise AssertionError("negative definiteness check changed")
    if B.T * Q.inv() * B != H:
        raise AssertionError("root-column identity changed")

    scale_symbol = sp.symbols("s", positive=True, integer=True)
    symbolic_gram = (scale_symbol * Q).row_join(B).col_join(B.T.row_join(H))
    symbolic_schur = sp.simplify(H - B.T * (scale_symbol * Q).inv() * B)
    expected_schur = (scale_symbol - 1) * H / scale_symbol
    if sp.simplify(symbolic_schur - expected_schur) != sp.zeros(3):
        raise AssertionError("symbolic Schur complement changed")
    determinant = sp.factor(symbolic_gram.det())
    expected_determinant = 64 * scale_symbol**3 * (scale_symbol - 1) ** 3
    if determinant != expected_determinant:
        raise AssertionError("symbolic determinant formula changed")

    projection_rows = []
    for j in range(1, 7):
        expression_at_s = sp.expand(2 * ((j - 1) * scale_symbol - j))
        minimum_at_s_ge_2 = -2 if j == 1 else int(expression_at_s.subs(scale_symbol, 2))
        projection_rows.append(
            {
                "A3_half_norm_j": j,
                "trace_square_formula": str(expression_at_s),
                "minimum_for_scale_at_least_two": minimum_at_s_ge_2,
                "visible_negative_possible": j == 1,
            }
        )

    samples = [sample_row(scale) for scale in (2, 3, 5, 11)]
    for row in samples:
        scale = row["scale"]
        if row["ambient_determinant"] != 64 * scale**3 * (scale - 1) ** 3:
            raise AssertionError("sample determinant formula failed")
        if row["ambient_signature_positive_negative_zero"] != [1, 8, 0]:
            raise AssertionError("sample signature failed")

    return {
        "schema": "coxeter.k3-compact-projected-root-scaling-no-go.v1",
        "report_sha256": sha256_file(REPORT),
        "source": {
            "golden_archive_sha256": SOURCE_ARCHIVE_SHA256,
            "triangle": triangle,
            "trace_gram_Q": matrix_rows(Q),
            "simple_root_gram_H": matrix_rows(H),
            "minus_H_leading_principal_minors": [2, 4, 4],
        },
        "all_parameter_family": {
            "integer_scale_domain": "s>=2",
            "trace_gram": "s*Q",
            "complement_gram": "s*(s-1)*H",
            "glue_vectors": "delta_i=(r_i+c_i)/s",
            "ambient_basis": "trace basis followed by delta_1,delta_2,delta_3",
            "ambient_gram_block_formula": "[[s*Q,Q*R],[R^T*Q,H]]",
            "schur_complement": "((s-1)/s)*H",
            "signature_positive_negative_zero": [1, 8, 0],
            "ambient_determinant": "64*s^3*(s-1)^3",
            "direct_sum_to_ambient_overlattice_index": "s^3",
            "trace_embedding_primitive": True,
            "orthogonal_complement_exactly_Ks": True,
            "ambient_even_integral": True,
            "trace_contains_minus_two_roots": False,
            "displayed_delta_roots_square_minus_two": True,
            "samples": samples,
        },
        "root_projection_exhaustion": {
            "root_equation": "Q(y,y)=-2*s-(s-1)*H(a,a)",
            "A3_norm_parameter": "H(a,a)=-2*j with integer j>=1",
            "simplified_trace_square": "2*((j-1)*s-j)",
            "rows": projection_rows,
            "visible_mirror_forces_negative_full_trace_square": True,
            "only_visible_case": "j=1 and Q(y,y)=-2",
            "all_visible_roots_project_to_old_golden_root_mirrors": True,
            "three_seed_facets_present_via_delta_i": True,
            "restricted_chamber_exactly_triangle_255": True,
        },
        "K3_gate": {
            "ambient_rank": 9,
            "ambient_signature": [1, 8],
            "K3_complement_rank": 13,
            "K3_complement_signature": [2, 11],
            "discriminant_length_upper_bound": 9,
            "stable_embedding_inequality": "13>=9+2",
            "primitive_K3_embedding_for_every_scale": True,
            "generic_period_exact_Picard": True,
            "Weyl_framing_exact_eigen_nef_return": True,
            "every_restricted_facet_has_physical_minus_two_wall": True,
            "every_physical_facet_root_lies_outside_trace": True,
        },
        "fixed_data": {
            "field": "Q(sqrt(5))",
            "order": "Z[(1+sqrt(5))/2]",
            "compact_polytope": "Coxeter triangle (2,5,5)",
            "normalized_Cartan_cycles": True,
            "native_field_reflections": True,
        },
        "unbounded_data": {
            "common_trace_scale": True,
            "glue_index": True,
            "ambient_discriminant": True,
            "negative_definite_NS_complement": True,
        },
        "scope": {
            "contradicts_literal_in_trace_physical_star_theorem": False,
            "violates_same_vector_in_trace_hypothesis": True,
            "produces_infinitely_many_abstract_chambers": False,
            "proves_projected_physical_walls_insufficient_for_finite_hosts": True,
            "integral_order_action_on_ambient_complement_claimed": False,
            "rational_divisor_product_correspondence_on_trace_block": True,
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archived projected-root scaling certificate differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
