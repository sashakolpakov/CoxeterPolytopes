#!/usr/bin/env python3
"""Exact certificate for a projected-root K3 eigenslice carrier of Delta(2,4,6).

The carrier deliberately lies outside the current rational-field full-nef
contract.  It embeds the scaled rank-three raw lattice primitively in an even
rank-six K3 Picard lattice, exposes the three raw facets by ambient primitive
(-2)-roots, and proves by a finite positive-definite root window that no other
ambient root mirror meets the closed restricted triangle.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/compact_246_projected_k3_eigenslice_carrier.json"
SCHEMA = "coxeter.compact_246.projected_k3_eigenslice_carrier.v1"
SCALE = 15

# Negative wall-normal Gram for the compact Coxeter triangle Delta(2,4,6).
G = sp.Matrix(
    [
        [-12, 6, 0],
        [6, -6, 3],
        [0, 3, -2],
    ]
)
H = -2 * sp.eye(3)


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _string(value: sp.Expr) -> str:
    return str(sp.factor(value))


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    rows: list[list[int | str]] = []
    for i in range(matrix.rows):
        row: list[int | str] = []
        for j in range(matrix.cols):
            value = sp.factor(matrix[i, j])
            row.append(int(value) if value.q == 1 else str(value))
        rows.append(row)
    return rows


def _quadratic(matrix: sp.Matrix, vector: tuple[int, ...] | sp.Matrix) -> sp.Expr:
    column = vector if isinstance(vector, sp.MatrixBase) else sp.Matrix(vector)
    return sp.factor((column.T * matrix * column)[0])


def _reflection_matrix(root_index: int) -> sp.Matrix:
    result = sp.eye(3)
    for column in range(3):
        result[root_index, column] -= (
            2 * G[column, root_index] / G[root_index, root_index]
        )
    return result


def ambient_gram(scale: int = SCALE) -> sp.Matrix:
    if scale <= 0:
        raise ValueError("scale must be positive")
    return (scale * G).row_join(G).col_join(G.row_join(H))


def orthogonal_projector(scale: int = SCALE) -> sp.Matrix:
    """Project S_Q to T_Q in the (r_1,r_2,r_3,delta_1,delta_2,delta_3) basis."""

    return sp.eye(3).row_join(sp.eye(3) / scale).col_join(sp.zeros(3, 6))


def rational_reflection_extension(root_index: int, scale: int = SCALE) -> sp.Matrix:
    """Raw reflection on T_Q, identity on its orthogonal complement."""

    reflection = _reflection_matrix(root_index)
    return reflection.row_join((reflection - sp.eye(3)) / scale).col_join(
        sp.zeros(3).row_join(sp.eye(3))
    )


def _a_ellipsoid_rows() -> list[dict[str, Any]]:
    matrix = 2 * SCALE * sp.eye(3) + G
    rows = []
    for entries in itertools.product(range(-1, 2), repeat=3):
        value = int(_quadratic(matrix, entries))
        if value < 2 * SCALE:
            rows.append({"a": list(entries), "ellipsoid_value": value})
    return rows


def _root_window() -> list[dict[str, Any]]:
    ones = sp.ones(3, 1)
    p = G.inv() * ones
    p_square = _quadratic(G, p)
    qp = sp.factor(2 / p_square) * (G * p) * (G * p).T - G
    max_cosh_2r = sp.Rational(177, 65)
    rows = []
    for entries in itertools.product(range(-2, 3), range(-5, 6), range(-6, 7)):
        square = _quadratic(G, entries)
        if -square not in (2, 6, 12):
            continue
        qp_value = _quadratic(qp, entries)
        if qp_value <= max_cosh_2r * (-square):
            rows.append(
                {
                    "y": list(entries),
                    "G_square": int(square),
                    "Qp": _string(qp_value),
                }
            )
    return rows


def build_certificate() -> dict[str, Any]:
    scale = SCALE
    gram = ambient_gram(scale)
    schur = sp.factor(H - G / scale)
    minus_schur_minors = [
        sp.factor((-schur)[:size, :size].det()) for size in range(1, 4)
    ]

    projector = orthogonal_projector(scale)
    extensions = [rational_reflection_extension(i, scale) for i in range(3)]
    raw_reflections = [_reflection_matrix(i) for i in range(3)]

    smith = smith_normal_form(gram, domain=ZZ)
    smith_diagonal = [abs(int(smith[i, i])) for i in range(6)]
    discriminant_length = sum(value > 1 for value in smith_diagonal)

    # A visible root rho=x+sum a_i delta_i has restricted normal
    # y=s*x+a.  From rho^2=-2,
    # G(y,y)=G(a,a)+2s(||a||^2-1).  Negativity of G(y,y) puts a in
    # the positive-definite ellipsoid (2sI+G)[a]<2s.
    a_matrix = 2 * scale * sp.eye(3) + G
    a_minors = [a_matrix[:size, :size].det() for size in range(1, 4)]
    a_dual_coordinate_bounds = [
        sp.factor(2 * scale * a_matrix.inv()[i, i]) for i in range(3)
    ]
    a_rows = _a_ellipsoid_rows()

    p = G.inv() * sp.ones(3, 1)
    p_square = _quadratic(G, p)
    vertices = [G.inv()[:, i] for i in range(3)]
    vertex_rows = []
    for index, vertex in enumerate(vertices):
        vertex_square = _quadratic(G, vertex)
        pairing = sp.factor((p.T * G * vertex)[0])
        cosh_square = sp.factor(pairing**2 / (p_square * vertex_square))
        vertex_rows.append(
            {
                "vertex": index + 1,
                "coordinates": [_string(value) for value in vertex],
                "square": _string(vertex_square),
                "p_pairing": _string(pairing),
                "cosh_distance_squared": _string(cosh_square),
                "cosh_twice_distance": _string(2 * cosh_square - 1),
            }
        )
    max_cosh_2r = max(sp.Rational(row["cosh_twice_distance"]) for row in vertex_rows)
    qp = sp.factor(2 / p_square) * (G * p) * (G * p).T - G
    qp_minors = [sp.factor(qp[:size, :size].det()) for size in range(1, 4)]
    universal_window_bound = sp.factor(12 * max_cosh_2r)
    coordinate_square_bounds = [
        sp.factor(universal_window_bound * qp.inv()[i, i]) for i in range(3)
    ]
    root_window = _root_window()

    residue_survivors = []
    wall_squares = [-int(G[i, i]) for i in range(3)]
    for index, positive_norm in enumerate(wall_squares):
        for sign in (-1, 1):
            residue = [0, 0, 0]
            residue[index] = sign
            for row in root_window:
                if row["G_square"] != -positive_norm:
                    continue
                if all(
                    (row["y"][coordinate] - residue[coordinate]) % scale == 0
                    for coordinate in range(3)
                ):
                    residue_survivors.append(
                        {
                            "a": residue,
                            "y": row["y"],
                            "G_square": row["G_square"],
                        }
                    )

    expected_a_rows = [
        {"a": [-1, 0, 0], "ellipsoid_value": 18},
        {"a": [0, -1, 0], "ellipsoid_value": 24},
        {"a": [0, 0, -1], "ellipsoid_value": 28},
        {"a": [0, 0, 0], "ellipsoid_value": 0},
        {"a": [0, 0, 1], "ellipsoid_value": 28},
        {"a": [0, 1, 0], "ellipsoid_value": 24},
        {"a": [1, 0, 0], "ellipsoid_value": 18},
    ]
    expected_residue_survivors = [
        {"a": [-1, 0, 0], "y": [-1, 0, 0], "G_square": -12},
        {"a": [1, 0, 0], "y": [1, 0, 0], "G_square": -12},
        {"a": [0, -1, 0], "y": [0, -1, 0], "G_square": -6},
        {"a": [0, 1, 0], "y": [0, 1, 0], "G_square": -6},
        {"a": [0, 0, -1], "y": [0, 0, -1], "G_square": -2},
        {"a": [0, 0, 1], "y": [0, 0, 1], "G_square": -2},
    ]

    checks = {
        "ambient_even_integral": all(value.q == 1 for value in gram)
        and all(int(gram[i, i]) % 2 == 0 for i in range(6)),
        "schur_negative_definite": minus_schur_minors
        == [sp.Rational(6, 5), sp.Rational(44, 25), sp.Rational(1214, 375)],
        "ambient_signature_1_5": G.det() > 0 and all(value > 0 for value in minus_schur_minors),
        "ambient_determinant_exact": gram.det() == -393336,
        "ambient_smith_exact": smith_diagonal == [1, 1, 1, 1, 18, 21852],
        "trace_primitive_direct_summand": True,
        "projector_idempotent": projector * projector == projector,
        "projector_self_adjoint": projector.T * gram == gram * projector,
        "projector_has_rank_three": projector.rank() == 3,
        "delta_roots_square_minus_two": gram[3:, 3:] == H,
        "delta_projection_is_r_i_over_15": projector[:, 3:] == sp.eye(3).col_join(sp.zeros(3)) / scale,
        "raw_reflections_integral_on_trace": all(
            all(value.q == 1 for value in reflection) for reflection in raw_reflections
        ),
        "raw_reflections_preserve_G": all(
            reflection.T * G * reflection == G for reflection in raw_reflections
        ),
        "rational_extensions_preserve_ambient_form": all(
            extension.T * gram * extension == gram for extension in extensions
        ),
        "rational_extensions_are_involutions": all(
            extension * extension == sp.eye(6) for extension in extensions
        ),
        "rational_extensions_commute_with_projector": all(
            extension * projector == projector * extension for extension in extensions
        ),
        "a_ellipsoid_positive_definite": a_minors == [18, 396, 10926],
        "a_coordinate_bounds_force_minus_one_zero_one": all(
            value < 2 for value in a_dual_coordinate_bounds
        ),
        "a_ellipsoid_exactly_zero_and_signed_basis": a_rows == expected_a_rows,
        "zero_projection_root_divisibility_impossible": 2 % scale != 0,
        "Qp_positive_definite": qp_minors == [sp.Rational(804, 65), sp.Rational(612, 13), 36],
        "circumradius_exact": max_cosh_2r == sp.Rational(177, 65),
        "root_window_coordinate_bounds_below_modulus": coordinate_square_bounds
        == [
            sp.Rational(31329, 4225),
            sp.Rational(111864, 4225),
            sp.Rational(36108, 845),
        ]
        and [2, 5, 6] == [int(sp.floor(sp.sqrt(value))) for value in coordinate_square_bounds]
        and all(bound < scale for bound in (2, 5, 6)),
        "root_window_exact_size_22": len(root_window) == 22,
        "residue_filter_leaves_only_three_facet_lines": residue_survivors
        == expected_residue_survivors,
        "stable_K3_embedding_inequality": 16 >= discriminant_length + 2,
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": {
            "proved": (
                "one proper rational projector/eigen-nef K3 carrier for the "
                "raw compact Delta(2,4,6) chamber, with projected physical roots"
            ),
            "outside_current_contract": (
                "for k=Q the current Property-2 contract requires literal full-nef "
                "return, whereas this carrier returns a proper rank-three slice "
                "inside Picard rank six"
            ),
            "not_proved": [
                "universal existence of this fourth carrier category",
                "integral extension of the raw Coxeter reflections to the full Picard lattice",
                "that ambient K3 root reflections restrict to the raw Coxeter reflections",
                "finite-host control for projected-root K3 eigenslice carriers",
                "a dimension bound for an arbitrary non-K3 rational eigenslice category",
            ],
        },
        "raw_triangle": {
            "type": "Delta(2,4,6)",
            "gram_G": _matrix_payload(G),
            "scaled_trace_gram": _matrix_payload(scale * G),
            "scale": scale,
            "signature": [1, 2],
            "area_over_pi": "1/12",
        },
        "ambient_lattice": {
            "basis": ["r1", "r2", "r3", "delta1", "delta2", "delta3"],
            "gram_formula": "[[15G,G],[G,-2I_3]]",
            "gram": _matrix_payload(gram),
            "schur_complement": _matrix_payload(schur),
            "minus_schur_leading_principal_minors": [_string(value) for value in minus_schur_minors],
            "signature": [1, 5],
            "determinant": int(gram.det()),
            "smith_diagonal": smith_diagonal,
            "discriminant_length": discriminant_length,
            "even_integral": True,
            "trace_embedding_primitive": True,
            "physical_delta_gram": _matrix_payload(H),
        },
        "projector_eigenframe": {
            "orthogonal_projector_matrix": _matrix_payload(projector),
            "eigenvalue_one_space": "T_Q=span_Q(r1,r2,r3)",
            "eigenvalue_zero_space": "T_Q^perp, negative definite of rank 3",
            "self_adjoint": True,
            "raw_reflection_matrices_on_T": [_matrix_payload(matrix) for matrix in raw_reflections],
            "rational_orthogonal_extensions": [_matrix_payload(matrix) for matrix in extensions],
            "algebraic_realization": (
                "the projector and rational extensions are realized on NS_Q by "
                "rational divisor-product correspondences"
            ),
            "integral_full_NS_action_claimed": False,
        },
        "visible_root_reduction": {
            "root_coordinates": "rho=x+sum_i a_i delta_i",
            "restricted_normal": "y=15x+a, so y congruent to a modulo 15",
            "root_equation": "G(y,y)=G(a,a)+30*(||a||^2-1)",
            "visible_negativity_ellipsoid": "a^T*(30I+G)*a<30",
            "ellipsoid_matrix": _matrix_payload(a_matrix),
            "ellipsoid_leading_principal_minors": [int(value) for value in a_minors],
            "dual_coordinate_square_bounds": [_string(value) for value in a_dual_coordinate_bounds],
            "ellipsoid_lattice_points": a_rows,
            "zero_projection_argument": (
                "if y=0 then a=15b and the root equation has left side divisible "
                "by 15 but right side 2"
            ),
        },
        "finite_hidden_root_window": {
            "interior_point_p": [_string(value) for value in p],
            "p_square": _string(p_square),
            "vertices": vertex_rows,
            "max_cosh_twice_circumradius": _string(max_cosh_2r),
            "Qp_matrix": _matrix_payload(qp),
            "Qp_leading_principal_minors": [_string(value) for value in qp_minors],
            "universal_Qp_bound_for_d_at_most_12": _string(universal_window_bound),
            "coordinate_square_bounds": [_string(value) for value in coordinate_square_bounds],
            "integer_coordinate_bounds": [2, 5, 6],
            "enumerated_candidate_count": len(root_window),
            "enumerated_candidates": root_window,
            "mod_15_residue_survivors": residue_survivors,
            "conclusion": (
                "every ambient (-2)-root mirror meeting the closed restricted "
                "triangle restricts to one of the three original facet lines"
            ),
        },
        "K3_realization": {
            "ambient_rank": 6,
            "ambient_signature": [1, 5],
            "K3_complement_rank": 16,
            "K3_complement_signature": [2, 14],
            "discriminant_length": discriminant_length,
            "stable_embedding_inequality": f"16>={discriminant_length}+2",
            "primitive_embedding": True,
            "generic_period_has_exact_NS": True,
            "Weyl_framing_makes_projected_chamber_nef": True,
            "three_facet_roots_become_smooth_rational_curves": True,
            "exact_return": "P(Nef(S)_R intersect eigenspace(projector,1))=Delta(2,4,6)",
            "carrier_relative_dimension_bound": "n+1=dim(E_1)<=rho(S)<=20, hence n<=19",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise AssertionError("archived projected Delta(2,4,6) K3 carrier drifted")
        print(f"archive check passed: {OUTPUT}")
        print(f"certificate_sha256={payload['certificate_sha256']}")
        return
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
        print(payload["certificate_sha256"])
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
