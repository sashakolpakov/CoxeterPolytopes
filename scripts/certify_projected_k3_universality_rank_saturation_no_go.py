"""Certify the projected-K3 rank-saturation obstruction and block formulas.

This certificate has two logically separate parts.

* It reconstructs the exact rank-20 mixed-root data of Vinberg's X3
  chamber and checks the rational square-class obstruction to an exact-module
  K3 eigen-nef return.
* It calibrates the general one-root-per-facet block construction on the
  Delta(2,4,6) Gram matrix.  The accompanying report proves the construction
  in arbitrary rank and states the hidden-root and K3-embedding gates.

The expensive finite-volume enumeration is not repeated here.  Its frozen
owner archive and source artifacts are hash-pinned, and the wall Gram is
independently rebuilt from the owner's combinatorial constructor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_vinberg_x3_h19 import build_wall_gram
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "projected_k3_universality_rank_saturation_no_go.json"
)

OWNER_PINS = {
    "results/vinberg_x3_h19.md": (
        "332a38ed8161ead76ef2f7e1755537ec22290c8929912ef4b022cf2bb0fa8833"
    ),
    "scripts/certify_vinberg_x3_h19.py": (
        "aff3ab46aa984c84131a5a6303f52de4038538d072420d5928d347e099f1fd41"
    ),
    "data/enumerations/vinberg_x3_h19.json": (
        "95284cd5808370d2ca4c182c52aafaba3f8e90e2f5a3a18843d53caaa48b1953"
    ),
    "tests/test_vinberg_x3_h19.py": (
        "55fce21ca5cdae1fc1d1ba299ab68a96309d164e168855b2e1a35c4a612e7d55"
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _nonzero_smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    ]


def direct_block(
    gram: sp.Matrix, facet_columns: sp.Matrix, scale: int
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return (ambient Gram, projector, negative Schur complement).

    If ``G`` is the trace Gram and the columns of ``V`` are facet normals,
    the ambient Gram is

        [[sG, GV], [V^tG, -2I]].

    The projector sends the i-th new root to v_i/s.
    """

    rank = gram.rows
    facets = facet_columns.cols
    if gram.cols != rank or facet_columns.rows != rank:
        raise ValueError("incompatible Gram and facet matrices")
    if scale <= 0:
        raise ValueError("scale must be positive")

    cross = gram * facet_columns
    ambient = (scale * gram).row_join(cross).col_join(
        cross.T.row_join(-2 * sp.eye(facets))
    )
    projector = sp.eye(rank).row_join(facet_columns / scale).col_join(
        sp.zeros(facets, rank + facets)
    )
    schur = -2 * sp.eye(facets) - facet_columns.T * gram * facet_columns / scale
    return ambient, projector, schur


def _block_calibration() -> dict[str, Any]:
    gram = sp.Matrix(
        [
            [-12, 6, 0],
            [6, -6, 3],
            [0, 3, -2],
        ]
    )
    facets = sp.eye(3)
    scale = 24
    ambient, projector, schur = direct_block(gram, facets, scale)

    root_identity_holds = True
    for x1 in range(-1, 2):
        for x2 in range(-1, 2):
            for x3 in range(-1, 2):
                x = sp.Matrix([x1, x2, x3])
                for a1 in range(-1, 2):
                    for a2 in range(-1, 2):
                        for a3 in range(-1, 2):
                            a = sp.Matrix([a1, a2, a3])
                            rho = x.col_join(a)
                            y = scale * x + facets * a
                            left = int((rho.T * ambient * rho)[0])
                            right = sp.Rational(
                                int((y.T * gram * y)[0])
                                - int((a.T * facets.T * gram * facets * a)[0]),
                                scale,
                            ) - 2 * int((a.T * a)[0])
                            if sp.Rational(left) != right:
                                root_identity_holds = False

    c_matrix = facets.T * gram * facets
    shell_matrix = 2 * scale * sp.eye(3) + c_matrix
    row_sum_bound = max(
        sum(abs(int(c_matrix[row, column])) for column in range(c_matrix.cols))
        for row in range(c_matrix.rows)
    )
    smith = _nonzero_smith_diagonal(ambient)
    discriminant_length = sum(value != 1 for value in smith)
    complement_rank = 22 - ambient.rank()

    checks = {
        "ambient_is_integral": all(value.q == 1 for value in ambient),
        "ambient_is_even": all(int(ambient[i, i]) % 2 == 0 for i in range(6)),
        "ambient_signature_is_1_5": exact_inertia(ambient) == (1, 5, 0),
        "schur_is_negative_definite": exact_inertia(schur) == (0, 3, 0),
        "projector_is_idempotent": projector * projector == projector,
        "projector_is_self_adjoint": projector.T * ambient == ambient * projector,
        "new_roots_have_square_minus_two": all(
            ambient[3 + i, 3 + i] == -2 for i in range(3)
        ),
        "new_roots_project_to_facets_over_scale": all(
            projector[:3, 3 + i] == facets[:, i] / scale for i in range(3)
        ),
        "root_identity_checked_on_729_pairs": root_identity_holds,
        "scale_exceeds_exact_row_sum_shell_bound": scale > row_sum_bound,
        "scale_is_even_and_divisible_by_facet_norms": (
            scale % 2 == 0
            and all(scale % abs(int(gram[i, i])) == 0 for i in range(3))
        ),
        "stable_k3_embedding_inequality": (
            complement_rank >= discriminant_length + 2
        ),
    }
    return {
        "raw_gram": [[int(value) for value in gram.row(i)] for i in range(3)],
        "facet_matrix": [[int(value) for value in facets.row(i)] for i in range(3)],
        "scale": scale,
        "row_sum_bound": row_sum_bound,
        "ambient_rank": ambient.rank(),
        "ambient_determinant": int(ambient.det()),
        "ambient_smith_nonzero": smith,
        "ambient_discriminant_length": discriminant_length,
        "k3_complement_rank": complement_rank,
        "ambient_gram": [
            [int(ambient[row, column]) for column in range(ambient.cols)]
            for row in range(ambient.rows)
        ],
        "projector": [
            [str(projector[row, column]) for column in range(projector.cols)]
            for row in range(projector.rows)
        ],
        "shell_matrix": [
            [int(shell_matrix[row, column]) for column in range(shell_matrix.cols)]
            for row in range(shell_matrix.rows)
        ],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _weak_mixed_divisor_calibration() -> dict[str, Any]:
    """Calibrate the root-free scaled-trace weak divisor construction."""

    gram = sp.Matrix(
        [
            [-12, 6, 0],
            [6, -6, 3],
            [0, 3, -2],
        ]
    )
    scale = 4
    scaled = scale * gram
    smith = _nonzero_smith_diagonal(scaled)
    discriminant_length = sum(value != 1 for value in smith)
    complement_rank = 22 - scaled.rank()
    checks = {
        "scaled_trace_is_integral_even": (
            all(value.q == 1 for value in scaled)
            and all(int(scaled[i, i]) % 2 == 0 for i in range(scaled.rows))
        ),
        "scaled_trace_signature_is_1_2": exact_inertia(scaled) == (1, 2, 0),
        "scale_greater_than_two_excludes_minus_two_norm": scale > 2,
        "rank_at_most_ten_makes_stable_embedding_pass": (
            scaled.rank() <= 10
            and complement_rank >= discriminant_length + 2
        ),
    }
    return {
        "scale": scale,
        "rank": scaled.rank(),
        "smith_nonzero": smith,
        "discriminant_length": discriminant_length,
        "k3_complement_rank": complement_rank,
        "root_free_reason": (
            "every scaled norm is divisible by s=4, hence cannot equal -2"
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _h19_obstruction() -> dict[str, Any]:
    owner_archive_path = ROOT / "data/enumerations/vinberg_x3_h19.json"
    owner_archive = json.loads(owner_archive_path.read_text())
    rebuilt = build_wall_gram()
    archived = sp.Matrix(owner_archive["wall_gram"])
    smith = _nonzero_smith_diagonal(rebuilt)

    diagonal = [int(rebuilt[index, index]) for index in range(rebuilt.rows)]
    norm_counts = {
        "minus_two": diagonal.count(-2),
        "minus_six": diagonal.count(-6),
    }

    # Under one rational similarity c, rational line multipliers a,b would
    # satisfy c*a^2=1 and 3*c*b^2=1.  Their quotient would square to 3.
    square_class_ratio = sp.Rational(6, 2)
    numerator_root, numerator_exact = sp.integer_nthroot(
        int(square_class_ratio.p), 2
    )
    denominator_root, denominator_exact = sp.integer_nthroot(
        int(square_class_ratio.q), 2
    )
    ratio_is_rational_square = numerator_exact and denominator_exact

    hyperbolic_dimension = int(owner_archive["hyperbolic_dimension"])
    eigen_rank = hyperbolic_dimension + 1
    k3_picard_ceiling = 20
    forced_kernel_rank = k3_picard_ceiling - eigen_rank

    checks = {
        "owner_pins_match": all(
            _sha256(ROOT / relative) == expected
            for relative, expected in OWNER_PINS.items()
        ),
        "rebuilt_wall_gram_matches_archive": rebuilt == archived,
        "wall_gram_rank_is_20": rebuilt.rank() == 20,
        "wall_gram_inertia_is_1_19_16": exact_inertia(rebuilt) == (1, 19, 16),
        "wall_gram_smith_is_1_19_then_3": smith == [1] * 19 + [3],
        "mixed_facet_norm_counts_are_24_and_12": norm_counts
        == {"minus_two": 24, "minus_six": 12},
        "owner_certifies_finite_volume": bool(owner_archive["volume"]["finite_volume"]),
        "owner_certifies_complete_strict_chamber": "complete walls"
        in owner_archive["completeness"]["mathematical_statement"],
        "h19_forces_picard_rank_20": eigen_rank == k3_picard_ceiling,
        "h19_forces_zero_projector_complement": forced_kernel_rank == 0,
        "facet_norm_ratio_three_is_not_a_rational_square": (
            square_class_ratio == 3 and not ratio_is_rational_square
        ),
    }
    return {
        "input": {
            "name": "Vinberg X3 mixed-root chamber",
            "lattice": "U + 2E8(-1) + A2(-1)",
            "hyperbolic_dimension": hyperbolic_dimension,
            "trace_rank": eigen_rank,
            "facet_count": rebuilt.rows,
            "facet_norm_counts": norm_counts,
            "wall_gram_rank": rebuilt.rank(),
            "wall_gram_inertia": list(exact_inertia(rebuilt)),
            "wall_gram_smith_nonzero": smith,
        },
        "rank_saturation": {
            "complex_k3_picard_ceiling": k3_picard_ceiling,
            "forced_ambient_picard_rank": eigen_rank,
            "forced_projector_kernel_rank": forced_kernel_rank,
        },
        "square_class_equations": {
            "raw_minus_two_facet": "c*a^2 = 1",
            "raw_minus_six_facet": "3*c*b^2 = 1",
            "forced_quotient": "(a/b)^2 = 3",
            "three_is_a_rational_square": ratio_is_rational_square,
            "three_adic_valuation": 1,
        },
        "owner_pins": OWNER_PINS,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def certify() -> dict[str, Any]:
    block = _block_calibration()
    weak_mixed = _weak_mixed_divisor_calibration()
    h19 = _h19_obstruction()
    checks = {
        "direct_block_calibration_passes": block["all_checks_pass"],
        "weak_mixed_divisor_calibration_passes": weak_mixed["all_checks_pass"],
        "h19_rank_saturation_obstruction_passes": h19["all_checks_pass"],
    }
    return {
        "artifact_status": "exact projected-K3 universality no-go and block calibration",
        "scope": {
            "disproved": (
                "Every finite-volume arithmetic Coxeter chamber with trace/eigen "
                "rank at most 20 has an exact-module projected K3 eigen-nef carrier."
            ),
            "counterexample_maximality": (
                "complete for every reflection of the pinned discriminant-three "
                "host; absolute maximality among all commensurable reflection "
                "overgroups is not asserted"
            ),
            "not_disproved": (
                "universality after forgetting the rational integral module, or a "
                "maximal-arithmetic-only projected-carrier conjecture"
            ),
        },
        "direct_block_calibration": block,
        "weak_mixed_divisor_calibration": weak_mixed,
        "h19_obstruction": h19,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    args = parser.parse_args()

    payload = certify()
    if args.write_archive:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        if not ARCHIVE.exists():
            raise FileNotFoundError(ARCHIVE)
        archived = json.loads(ARCHIVE.read_text())
        if archived != payload:
            raise AssertionError("projected-K3 universality archive has drifted")
    if not payload["all_checks_pass"]:
        failed = [key for key, value in payload["checks"].items() if not value]
        raise AssertionError(f"failed projected-K3 checks: {failed}")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
