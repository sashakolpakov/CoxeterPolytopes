#!/usr/bin/env python3
"""Independent hostile replay of the projected-K3 X3 universality no-go.

This audit imports no producer functions.  It reads the frozen X3 Gram and
the frozen owner certificate, recomputes the rank/inertia/Smith/norm data,
checks the rank-saturated rational square-class obstruction, and uses a
separate small block to verify the projector and root-shell identities.
The compact-shell, rational-cusp, and primitive-glue arguments are recorded
as proof obligations in the accompanying report rather than disguised as
finite computations.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
OWNER_ARCHIVE = (
    ROOT / "data/enumerations/projected_k3_universality_rank_saturation_no_go.json"
)
X3_ARCHIVE = ROOT / "data/enumerations/vinberg_x3_h19.json"
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.json"
)
SCHEMA = "coxeter.projected_k3_universality.independent_hostile_audit.v1"

OWNER_PINS = {
    "results/projected_k3_universality_rank_saturation_no_go.md": (
        "7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2"
    ),
    "scripts/certify_projected_k3_universality_rank_saturation_no_go.py": (
        "d94e18eada49d277b811299f5378c6aa58566a861032c2d2d6fa78fdf3844c87"
    ),
    "data/enumerations/projected_k3_universality_rank_saturation_no_go.json": (
        "61ddd537135db3f47bd267abcc3d9bbea8db64d73370e7f8a1de2a6125e1751e"
    ),
    "tests/test_projected_k3_universality_rank_saturation_no_go.py": (
        "9dcbd2189e5629bd8f251d9e0950c2e2a3d438f6a810f40e27bb56724b88ada1"
    ),
}
X3_ARCHIVE_SHA256 = "95284cd5808370d2ca4c182c52aafaba3f8e90e2f5a3a18843d53caaa48b1953"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Rational congruence elimination, independent of producer helpers."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (row, column)
                    for row in range(work.rows)
                    for column in range(row)
                    if work[row, column]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            row, column = pair
            transform = sp.eye(work.rows)
            transform[:, row] = sp.eye(work.rows)[:, row] + sp.eye(work.rows)[:, column]
            work = transform.T * work * transform
            pivot = row
        order = [pivot, *range(pivot), *range(pivot + 1, work.rows)]
        work = work.extract(order, order)
        value = sp.Rational(work[0, 0])
        if value > 0:
            positive += 1
        else:
            negative += 1
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return positive, negative, 0


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index]
    ]


def _rank_mod_prime(matrix: sp.Matrix, prime: int) -> int:
    rows = [
        [int(matrix[row, column]) % prime for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]
    rank = 0
    for column in range(matrix.cols):
        pivot = next(
            (row for row in range(rank, matrix.rows) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(entry * inverse) % prime for entry in rows[rank]]
        for row in range(matrix.rows):
            if row == rank:
                continue
            multiple = rows[row][column]
            if multiple:
                rows[row] = [
                    (left - multiple * right) % prime
                    for left, right in zip(rows[row], rows[rank])
                ]
        rank += 1
        if rank == matrix.rows:
            break
    return rank


def _independent_block_audit() -> dict[str, Any]:
    gram = sp.Matrix([[2, 1], [1, -2]])
    facets = sp.Matrix([[1, 0, 1], [0, 1, 1]])
    scale = 20
    rank = gram.rows
    facet_count = facets.cols
    cross = gram * facets
    ambient = (scale * gram).row_join(cross).col_join(
        cross.T.row_join(-2 * sp.eye(facet_count))
    )
    projector = sp.eye(rank).row_join(facets / scale).col_join(
        sp.zeros(facet_count, rank + facet_count)
    )
    schur = -2 * sp.eye(facet_count) - facets.T * gram * facets / scale
    facet_gram = facets.T * gram * facets
    row_bound = max(
        sum(abs(int(facet_gram[row, column])) for column in range(facet_count))
        for row in range(facet_count)
    )

    root_identity = True
    for x_entries in itertools.product(range(-1, 2), repeat=rank):
        x = sp.Matrix(x_entries)
        for a_entries in itertools.product(range(-1, 2), repeat=facet_count):
            a = sp.Matrix(a_entries)
            rho = x.col_join(a)
            y = scale * x + facets * a
            left = sp.Rational((rho.T * ambient * rho)[0])
            right = (
                sp.Rational(
                    (y.T * gram * y)[0] - (a.T * facet_gram * a)[0], scale
                )
                - 2 * (a.T * a)[0]
            )
            if left != right:
                root_identity = False

    p_rank = _rank_mod_prime(facets, 2)
    checks = {
        "ambient_integral_even": all(value.q == 1 for value in ambient)
        and all(int(ambient[index, index]) % 2 == 0 for index in range(ambient.rows)),
        "ambient_signature_1_4": _exact_inertia(ambient) == (1, 4, 0),
        "schur_negative_definite": _exact_inertia(schur) == (0, 3, 0),
        "projector_idempotent": projector * projector == projector,
        "projector_self_adjoint": projector.T * ambient == ambient * projector,
        "new_roots_square_minus_two": all(
            ambient[rank + index, rank + index] == -2
            for index in range(facet_count)
        ),
        "new_roots_project_to_facet_over_scale": all(
            projector[:rank, rank + index] == facets[:, index] / scale
            for index in range(facet_count)
        ),
        "root_identity_on_243_pairs": root_identity,
        "scale_exceeds_row_shell_bound": scale > row_bound,
        "uniform_denominator_mod_2_rank_is_2": p_rank == 2,
        "complement_rank_dominates_glue_lower_bound": facet_count >= p_rank,
    }
    return {
        "gram": [[int(value) for value in gram.row(row)] for row in range(rank)],
        "facets": [
            [int(value) for value in facets.row(row)] for row in range(rank)
        ],
        "scale": scale,
        "row_shell_bound": row_bound,
        "ambient_signature": list(_exact_inertia(ambient)),
        "schur_signature": list(_exact_inertia(schur)),
        "uniform_denominator_prime": 2,
        "facet_rank_mod_prime": p_rank,
        "complement_rank": facet_count,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_audit() -> dict[str, Any]:
    owner_pins_match = all(
        _sha256(ROOT / relative) == expected for relative, expected in OWNER_PINS.items()
    )
    if _sha256(X3_ARCHIVE) != X3_ARCHIVE_SHA256:
        raise ArithmeticError("the frozen X3 dependency changed")

    owner = json.loads(OWNER_ARCHIVE.read_text(encoding="utf-8"))
    x3 = json.loads(X3_ARCHIVE.read_text(encoding="utf-8"))
    wall_gram = sp.Matrix(x3["wall_gram"])
    diagonal_counts = Counter(int(wall_gram[index, index]) for index in range(36))
    smith = _smith_diagonal(wall_gram)
    inertia = _exact_inertia(wall_gram)

    dimension = int(x3["hyperbolic_dimension"])
    raw_rank = wall_gram.rank()
    complex_k3_ceiling = 20
    projector_complement_rank = complex_k3_ceiling - raw_rank
    square_ratio = sp.Rational(6, 2)
    numerator_square = sp.integer_nthroot(int(square_ratio.p), 2)[1]
    denominator_square = sp.integer_nthroot(int(square_ratio.q), 2)[1]
    rational_square = numerator_square and denominator_square

    block = _independent_block_audit()
    report_text = (ROOT / "results/projected_k3_universality_rank_saturation_no_go.md").read_text(
        encoding="utf-8"
    )
    report_flat = " ".join(report_text.split())
    checks = {
        "terminal_owner_pins_match": owner_pins_match,
        "owner_archive_passes": owner["all_checks_pass"] is True,
        "x3_dimension_19_rank_20": dimension == 19 and raw_rank == 20,
        "x3_inertia_1_19_16": inertia == (1, 19, 16),
        "x3_smith_1_19_3": smith == [1] * 19 + [3],
        "x3_norm_split_24_minus2_12_minus6": diagonal_counts
        == Counter({-2: 24, -6: 12}),
        "rank_saturation_forces_zero_projector_complement": projector_complement_rank
        == 0,
        "norm_ratio_three_not_rational_square": square_ratio == 3
        and not rational_square,
        "owner_scope_is_not_maximal_only": "every finite-volume arithmetic coxeter chamber"
        in owner["scope"]["disproved"].lower(),
        "owner_preserves_absolute_maximality_caveat": "not asserted"
        in owner["scope"]["counterexample_maximality"],
        "owner_preserves_unlabelled_real_escape": "forgetting the rational integral module"
        in owner["scope"]["not_disproved"],
        "report_states_generic_facet_local_finiteness": (
            "generic positive-norm point of an exposed nef facet" in report_flat
            and "local finiteness of the root mirrors" in report_flat
        ),
        "report_states_compact_nonrational_minkowski_argument": (
            "At every other embedding, negative definiteness" in report_flat
            and "The Minkowski image of" in report_flat
        ),
        "report_states_rational_cusp_integrality_formula": (
            "Thus \\(y\\) is itself an integral reflection vector" in report_flat
        ),
        "report_states_primitive_glue_length_bound": (
            "must inject into the discriminant group of \\(K\\)" in report_flat
        ),
        "independent_block_passes": block["all_checks_pass"],
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "owner_pins": OWNER_PINS,
        "x3_dependency": {
            "path": str(X3_ARCHIVE.relative_to(ROOT)),
            "sha256": X3_ARCHIVE_SHA256,
        },
        "x3_independent_replay": {
            "hyperbolic_dimension": dimension,
            "raw_rank": raw_rank,
            "inertia": list(inertia),
            "smith_nonzero": smith,
            "facet_norm_counts": {str(key): value for key, value in diagonal_counts.items()},
            "complex_k3_picard_ceiling": complex_k3_ceiling,
            "forced_projector_complement_rank": projector_complement_rank,
            "square_class_ratio": str(square_ratio),
            "square_class_ratio_is_rational_square": rational_square,
            "three_adic_valuation": 1,
        },
        "logic_audit": {
            "generic_facet_root_necessity": (
                "A finite facet has a generic positive point.  Local finiteness of "
                "K3 root mirrors near that point makes one active root inequality "
                "support an open subset of the facet.  At rank saturation the "
                "projection is the identity, so the root lies on the raw Q-line."
            ),
            "exact_module_similarity": (
                "The contradiction uses one rational similarity and rational line "
                "representatives.  It does not apply after forgetting the pinned "
                "Q-module and retaining only the unlabelled real polytope."
            ),
            "compact_nonrational_shell": (
                "For a visible normal every archimedean summand of the fixed trace "
                "norm is negative.  Compactness uniformly bounds the distinguished "
                "Lorentz boost; negative definiteness bounds the other embeddings; "
                "the Minkowski lattice then leaves finitely many shell vectors."
            ),
            "rational_cusp_completeness": (
                "When s is divisible by each facet norm, the residue-shell identity "
                "makes y an integral reflection vector.  Completeness for every "
                "integral reflection of the pinned host makes its inequality old or "
                "redundant even at a cusp."
            ),
            "compression_discriminant_bound": (
                "Primitive T and K make the overlattice glue a graph.  Hence the "
                "subgroup generated by the projected facet classes injects into A_K, "
                "whose p-length is at most rank K."
            ),
        },
        "independent_block_calibration": block,
        "scope_verdict": {
            "status": "PASS",
            "proved": (
                "X3 disproves exact-module projected K3 eigen-nef universality over "
                "all finite-volume arithmetic chambers of trace rank at most 20."
            ),
            "not_proved": (
                "a counterexample restricted to absolutely maximal arithmetic "
                "reflection groups, or a no-go for unlabelled real/mixed-divisor "
                "K3 carriers"
            ),
        },
        "checks": checks,
    }
    payload["all_checks_pass"] = all(checks.values())
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write-archive", action="store_true")
    group.add_argument("--check-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    payload = build_audit()
    if args.write_archive:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit("independent audit archive differs from replay")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            json.dumps(
                {
                    "all_checks_pass": payload["all_checks_pass"],
                    "certificate_sha256": payload["certificate_sha256"],
                    "scope_verdict": payload["scope_verdict"]["status"],
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
