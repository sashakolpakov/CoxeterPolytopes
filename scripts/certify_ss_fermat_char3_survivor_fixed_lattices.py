#!/usr/bin/env python3
"""Certify the saturated fixed lattices of the four char-3 volume survivors.

The inputs are the exact full-stabilizer cyclic screen and the exact V4 fold
screen for the complete 5,944-wall Fermat-quartic chamber.  Smith kernels are
computed in the published 22-dimensional Neron--Severi basis.  In particular,
the V4-H10 lattice is identified with ``U+E8(-2)+<-4>`` by its full
discriminant form and Nikulin's indefinite uniqueness criterion.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.certify_jacobian_kummer_h10_class0007_fixed_lattice import (
    _discriminant_basis,
    _finite_form_split,
    _target_gram,
)
from scripts.search_ss_fermat_quartic_char3_folds import SOURCE, parse_source
from scripts.screen_ss_fermat_quartic_char3_c2_survivors import wall_cast
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CYCLIC = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_element_class_screen.json"
)
V4 = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_v4_fold_screen.json"
)
C2_REDUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_survivor_screen.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)

# Columns give the target standard basis in the archived short V4-H10 fixed
# basis.  It was found by splitting a primitive U and applying exact
# positive-definite isometry on the rank-nine complement.  The certificate
# replays it entry by entry, so no external isometry program is required.
V4_H10_TO_TARGET_ISOMETRY = sp.Matrix(
    [
        [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, -2],
        [1, 2, -2, -1, 0, 1, 1, -1, 0, 2, 1],
        [0, 0, 1, 0, 0, -1, 1, 1, -2, 0, 1],
        [0, 0, -1, -1, -1, 2, 0, -1, 1, -1, -2],
        [0, 0, -1, 0, 0, 1, 0, -1, 0, 0, -1],
        [0, 1, 0, -1, 0, 1, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, -1, 0, 1, -1, 1, 0, 1],
        [0, 0, -1, 1, 1, -1, 0, 1, -1, 0, 0],
        [0, 0, 1, 1, 0, -1, 0, 1, -1, 0, 1],
        [0, 0, 1, -2, 0, 2, -3, 1, 1, 0, -2],
        [0, 0, 0, -1, -1, 2, -1, 0, 0, 0, -2],
    ]
)

# Columns give the standard basis of McLeod's determinant-three lattice in
# the archived short C2-H10 fixed basis.  The fixed Gram is twice an odd
# integral form; this matrix certifies
#
#     P^T * (G_C2 / 2) * P = diag(3,-1,...,-1).
#
# Multiplying a quadratic form by 2 changes neither its integral orthogonal
# group nor the crystallographic reflection formula, so this identifies the
# ambient reflection group with the one studied by McLeod for f_3^10.
C2_H10_TO_MCLEOD_ISOMETRY = sp.Matrix(
    [
        [-3, 0, 0, -1, 0, -1, 1, -1, -1, -1, 0],
        [-1, 0, 0, 0, 0, -1, 0, 1, 0, -1, -1],
        [8, -1, 1, 2, 0, 2, -1, 1, 1, 3, 3],
        [-3, 0, -1, 0, 1, 0, 0, 0, -1, 0, -2],
        [0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0],
        [-2, 0, 0, -1, 0, -1, 0, 0, 0, -1, -1],
        [4, -1, 1, 1, 0, 1, 0, 0, 0, 2, 2],
        [-1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1],
        [-1, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0],
        [0, 0, 1, 0, -1, -1, 0, 0, 0, -1, 0],
        [-4, 1, 0, -1, 0, -1, 1, -1, -1, -1, -1],
    ]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _primitive_integer_column(vector: sp.Matrix) -> sp.Matrix:
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in vector]))
    entries = [int(sp.Rational(value) * denominator) for value in vector]
    content = math.gcd(*(abs(entry) for entry in entries))
    if content == 0:
        raise ArithmeticError("zero mirror vector")
    entries = [entry // content for entry in entries]
    if next(entry for entry in entries if entry) < 0:
        entries = [-entry for entry in entries]
    return sp.Matrix(entries)


def _smith_decomposition(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    diagonal, left, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(ZZ)
    )
    return diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()


def _saturated_fixed_basis(
    generators: Sequence[sp.Matrix], ambient_rank: int
) -> tuple[sp.Matrix, dict[str, Any]]:
    if not generators:
        raise ValueError("at least one generator is required")
    identity = sp.eye(ambient_rank)
    equations = generators[0].T - identity
    for generator in generators[1:]:
        equations = equations.col_join(generator.T - identity)
    diagonal, _, right = _smith_decomposition(equations)
    equation_rank = equations.rank()
    smith_basis = right[:, equation_rank:].T
    reduced_basis = smith_basis.lll(delta=sp.Rational(3, 4))
    _, pivot_columns = smith_basis.rref()
    pivot_columns = tuple(pivot_columns)
    basis_change = (
        reduced_basis[:, list(pivot_columns)]
        * smith_basis[:, list(pivot_columns)].inv()
    )
    if not (
        basis_change * smith_basis == reduced_basis
        and abs(int(basis_change.det())) == 1
        and all(reduced_basis * generator == reduced_basis for generator in generators)
    ):
        raise ArithmeticError("LLL basis is not the same saturated fixed lattice")
    smith_diagonal = [
        int(diagonal[index, index])
        for index in range(min(diagonal.shape))
        if diagonal[index, index] != 0
    ]
    return reduced_basis, {
        "equation_matrix_rank": equation_rank,
        "nonzero_smith_diagonal": smith_diagonal,
        "kernel_rank": reduced_basis.rows,
        "lll_basis_change_determinant": int(basis_change.det()),
        "saturated_smith_kernel": True,
    }


def _lattice_record(
    case_id: str,
    subgroup_type: str,
    hyperbolic_dimension: int,
    generators: Sequence[sp.Matrix],
    ambient_gram: sp.Matrix,
    source_metadata: dict[str, Any],
) -> dict[str, Any]:
    ambient_rank = ambient_gram.rows
    identity = sp.eye(ambient_rank)
    if any(generator.shape != (ambient_rank, ambient_rank) for generator in generators):
        raise ArithmeticError("generator has the wrong shape")
    if any(generator * ambient_gram * generator.T != ambient_gram for generator in generators):
        raise ArithmeticError("generator does not preserve the Neron--Severi Gram")
    if any(generator == identity for generator in generators):
        raise ArithmeticError("a survivor generator is trivial")

    basis, kernel = _saturated_fixed_basis(generators, ambient_rank)
    gram = basis * ambient_gram * basis.T
    rank = gram.rows
    smith = smith_normal_form(gram, domain=sp.ZZ)
    smith_diagonal = [abs(int(smith[index, index])) for index in range(rank)]
    determinant = abs(int(gram.det()))
    inertia = exact_inertia(gram)
    even = all(int(gram[index, index]) % 2 == 0 for index in range(rank))
    if not (
        rank == hyperbolic_dimension + 1
        and inertia == (1, hyperbolic_dimension, 0)
        and even
    ):
        raise ArithmeticError(f"unexpected fixed-lattice data for {case_id}")
    return {
        "id": case_id,
        "subgroup_type": subgroup_type,
        "hyperbolic_dimension": hyperbolic_dimension,
        "source_metadata": source_metadata,
        "generator_count": len(generators),
        "generator_orders": [
            next(
                order
                for order in range(1, 13)
                if generator**order == identity
            )
            for generator in generators
        ],
        "smith_kernel": kernel,
        "basis_rows_in_published_NS_basis": _matrix_payload(basis),
        "gram": _matrix_payload(gram),
        "rank": rank,
        "signature": [1, hyperbolic_dimension],
        "determinant_absolute": determinant,
        "smith_diagonal_absolute": smith_diagonal,
        "discriminant_group_order": determinant,
        "discriminant_group_length": sum(value > 1 for value in smith_diagonal),
        "even": even,
    }


def build_certificate() -> dict[str, Any]:
    source = parse_source(SOURCE)
    ambient_gram: sp.Matrix = source["gram"]
    if ambient_gram.shape != (22, 22) or exact_inertia(ambient_gram) != (1, 21, 0):
        raise ArithmeticError("unexpected Fermat-quartic Neron--Severi lattice")

    cyclic = json.loads(CYCLIC.read_text(encoding="utf-8"))
    cyclic_survivors = [
        record
        for record in cyclic["high_fixed_rank_casts"]
        if record["bad_pair_count"] == 0
    ]
    if sorted(record["hyperbolic_dimension"] for record in cyclic_survivors) != [10, 11]:
        raise ArithmeticError("the cyclic survivor queue changed")

    v4 = json.loads(V4.read_text(encoding="utf-8"))
    if sorted(cast["subgroup"]["hyperbolic_dimension"] for cast in v4["casts"]) != [10, 11]:
        raise ArithmeticError("the V4 survivor queue changed")

    cases: list[dict[str, Any]] = []
    for record in sorted(cyclic_survivors, key=lambda row: -row["hyperbolic_dimension"]):
        dimension = int(record["hyperbolic_dimension"])
        cases.append(
            _lattice_record(
                f"char3-C2-H{dimension}",
                "C2",
                dimension,
                [sp.Matrix(record["representative_matrix"])],
                ambient_gram,
                {
                    "gap_element_class_index": int(record["gap_class_index"]),
                    "element_class_size": int(record["class_size"]),
                    "strict_facet_count_before_volume_audit": int(record["strict_facet_count"]),
                },
            )
        )
    for cast in sorted(v4["casts"], key=lambda row: -row["subgroup"]["hyperbolic_dimension"]):
        dimension = int(cast["subgroup"]["hyperbolic_dimension"])
        cases.append(
            _lattice_record(
                f"char3-V4-H{dimension}",
                "C2xC2",
                dimension,
                [sp.Matrix(matrix) for matrix in cast["subgroup"]["generator_ns_matrices"]],
                ambient_gram,
                {
                    "subgroup_conjugacy_class_size": int(
                        cast["subgroup"]["conjugacy_class_size"]
                    ),
                    "normalizer_order": int(cast["subgroup"]["normalizer_order"]),
                    "irredundant_facet_count": int(
                        cast["global_facet_reduction"]["irredundant_facet_count"]
                    ),
                },
            )
        )

    by_id = {case["id"]: case for case in cases}
    expected = {
        "char3-C2-H11": (1_024, [1, 1, *([2] * 10)]),
        "char3-C2-H10": (6_144, [*([2] * 10), 6]),
        "char3-V4-H11": (768, [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 12]),
        "char3-V4-H10": (1_024, [1, 1, *([2] * 8), 4]),
    }
    fingerprints_match = all(
        by_id[case_id]["determinant_absolute"] == determinant
        and by_id[case_id]["smith_diagonal_absolute"] == smith
        for case_id, (determinant, smith) in expected.items()
    )

    v4_h10_gram = sp.Matrix(by_id["char3-V4-H10"]["gram"])
    orders, _, discriminant_form = _discriminant_basis(v4_h10_gram)
    v4_h10_split = _finite_form_split(discriminant_form, orders)
    target = _target_gram()
    target_orders, _, target_discriminant_form = _discriminant_basis(target)
    target_split = _finite_form_split(target_discriminant_form, target_orders)
    discriminant_forms_match = (
        v4_h10_split["all_checks_pass"]
        and target_split["all_checks_pass"]
        and v4_h10_split["canonical_form"] == target_split["canonical_form"]
    )
    nikulin_applies = (
        by_id["char3-V4-H10"]["rank"]
        >= by_id["char3-V4-H10"]["discriminant_group_length"] + 2
    )
    explicit_isometry = V4_H10_TO_TARGET_ISOMETRY
    explicit_isometry_replays = (
        abs(int(explicit_isometry.det())) == 1
        and explicit_isometry.T * v4_h10_gram * explicit_isometry == target
    )

    c2_h10_gram = sp.Matrix(by_id["char3-C2-H10"]["gram"])
    c2_h10_half_gram = c2_h10_gram / 2
    mcleod_target = sp.diag(3, *([-1] * 10))
    c2_h10_isometry = C2_H10_TO_MCLEOD_ISOMETRY
    c2_h10_half_is_integral = all(
        value.is_Integer is True for value in c2_h10_half_gram
    )
    c2_h10_explicit_isometry_replays = (
        c2_h10_half_is_integral
        and abs(int(c2_h10_isometry.det())) == 1
        and c2_h10_isometry.T * c2_h10_half_gram * c2_h10_isometry
        == mcleod_target
    )

    c2_h10_record = next(
        record
        for record in cyclic_survivors
        if record["hyperbolic_dimension"] == 10
    )
    c2_h10_cast = wall_cast(
        np.asarray(c2_h10_record["representative_matrix"], dtype=np.int64)
    )
    c2_reduction = json.loads(C2_REDUCTION.read_text(encoding="utf-8"))
    c2_h10_reduction = next(
        cast
        for cast in c2_reduction["casts"]
        if cast["hyperbolic_dimension"] == 10
    )
    c2_h10_basis = sp.Matrix(
        by_id["char3-C2-H10"]["basis_rows_in_published_NS_basis"]
    )
    c2_h10_mirror_records: list[dict[str, Any]] = []
    c2_h10_mirror_types: Counter[tuple[int, int]] = Counter()
    c2_h10_all_reflections_integral = True
    for restricted_ray_index in c2_h10_reduction["facet_restricted_ray_indices"]:
        ambient_dual_normal = sp.Matrix(
            c2_h10_cast["rays"][int(restricted_ray_index) - 1]
        )
        fixed_coordinates = _primitive_integer_column(
            c2_h10_half_gram.inv() * c2_h10_basis * ambient_dual_normal
        )
        square = int(
            (fixed_coordinates.T * c2_h10_half_gram * fixed_coordinates)[0]
        )
        divisibility = math.gcd(
            *(
                abs(int(value))
                for value in c2_h10_half_gram * fixed_coordinates
            )
        )
        reflection = sp.eye(11) - sp.Rational(2, square) * fixed_coordinates * (
            fixed_coordinates.T * c2_h10_half_gram
        )
        reflection_valid = (
            all(value.is_Integer is True for value in reflection)
            and reflection.T * c2_h10_half_gram * reflection == c2_h10_half_gram
            and reflection * reflection == sp.eye(11)
            and reflection * fixed_coordinates == -fixed_coordinates
        )
        c2_h10_all_reflections_integral = (
            c2_h10_all_reflections_integral and reflection_valid
        )
        target_coordinates = c2_h10_isometry.inv() * fixed_coordinates
        if any(value.is_Integer is not True for value in target_coordinates):
            raise ArithmeticError("C2-H10 mirror left the McLeod lattice")
        c2_h10_mirror_types[(square, divisibility)] += 1
        c2_h10_mirror_records.append(
            {
                "restricted_ray_index": int(restricted_ray_index),
                "fixed_basis_coordinates": [
                    int(value) for value in fixed_coordinates
                ],
                "mcleod_lattice_coordinates": [
                    int(value) for value in target_coordinates
                ],
                "half_scaled_square": square,
                "half_scaled_divisibility": divisibility,
                "original_even_square": 2 * square,
                "original_even_divisibility": 2 * divisibility,
                "reflection_integral_involutive_isometry": reflection_valid,
            }
        )

    v4_h10_case = by_id["char3-V4-H10"]
    v4_h10_basis = sp.Matrix(v4_h10_case["basis_rows_in_published_NS_basis"])
    v4_h10_cast = next(
        cast for cast in v4["casts"] if cast["subgroup"]["hyperbolic_dimension"] == 10
    )
    mirror_records: list[dict[str, Any]] = []
    mirror_types: Counter[tuple[int, int]] = Counter()
    all_reflections_integral = True
    for facet in v4_h10_cast["strict_facet_gate"]["certificates"]:
        ambient_dual_normal = sp.Matrix(facet["primitive_normal_covector"])
        fixed_coordinates = _primitive_integer_column(
            v4_h10_gram.inv() * v4_h10_basis * ambient_dual_normal
        )
        square = int((fixed_coordinates.T * v4_h10_gram * fixed_coordinates)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in v4_h10_gram * fixed_coordinates)
        )
        reflection = sp.eye(11) - sp.Rational(2, square) * fixed_coordinates * (
            fixed_coordinates.T * v4_h10_gram
        )
        reflection_integral = all(value.is_Integer is True for value in reflection)
        reflection_valid = (
            reflection_integral
            and reflection.T * v4_h10_gram * reflection == v4_h10_gram
            and reflection * reflection == sp.eye(11)
            and reflection * fixed_coordinates == -fixed_coordinates
        )
        all_reflections_integral = all_reflections_integral and reflection_valid
        target_coordinates = explicit_isometry.inv() * fixed_coordinates
        if any(value.is_Integer is not True for value in target_coordinates):
            raise ArithmeticError("explicit isometry sent a mirror outside the target lattice")
        mirror_types[(square, divisibility)] += 1
        mirror_records.append(
            {
                "restricted_ray_index": int(facet["restricted_ray_index"]),
                "fixed_basis_coordinates": [int(value) for value in fixed_coordinates],
                "target_lattice_coordinates": [int(value) for value in target_coordinates],
                "square": square,
                "divisibility": divisibility,
                "reflection_integral_involutive_isometry": reflection_valid,
            }
        )
    checks = {
        "exactly_four_angle_gate_survivors": len(cases) == 4,
        "all_fixed_bases_are_saturated_smith_kernels": all(
            case["smith_kernel"]["saturated_smith_kernel"] for case in cases
        ),
        "all_fixed_lattices_are_even_lorentzian": all(
            case["even"] and case["signature"] == [1, case["hyperbolic_dimension"]]
            for case in cases
        ),
        "four_fingerprints_match": fingerprints_match,
        "v4_h10_discriminant_form_matches_target": discriminant_forms_match,
        "nikulin_rank_length_criterion_applies": nikulin_applies,
        "explicit_unimodular_isometry_replays": explicit_isometry_replays,
        "c2_h10_half_scaled_lattice_is_mcleod_f3_h10": (
            c2_h10_explicit_isometry_replays
        ),
        "all_196_facets_are_primitive_integral_reflecting_roots": (
            len(c2_h10_mirror_records) == 196
            and c2_h10_all_reflections_integral
            and len(
                {
                    tuple(record["mcleod_lattice_coordinates"])
                    for record in c2_h10_mirror_records
                }
            )
            == 196
            and c2_h10_mirror_types
            == Counter({(-1, 1): 80, (-2, 1): 36, (-6, 3): 80})
        ),
        "all_123_facets_are_primitive_integral_reflecting_roots": (
            len(mirror_records) == 123
            and all_reflections_integral
            and mirror_types
            == Counter({(-4, 4): 82, (-4, 2): 31, (-2, 1): 10})
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"fixed-lattice checks failed: {checks}")

    return {
        "schema": "ss_fermat_quartic_char3_survivor_fixed_lattices_v1",
        "source": {
            "surface": "supersingular Fermat quartic in characteristic 3",
            "paper": "Kondo--Shimada, The automorphism group of a supersingular K3 surface with Artin invariant 1 in characteristic 3",
            "paper_arxiv": "https://arxiv.org/abs/1205.6520",
            "published_ns_source": str(SOURCE.relative_to(ROOT)),
            "published_ns_source_sha256": _sha256(SOURCE),
            "cyclic_screen": str(CYCLIC.relative_to(ROOT)),
            "cyclic_screen_sha256": _sha256(CYCLIC),
            "v4_screen": str(V4.relative_to(ROOT)),
            "v4_screen_sha256": _sha256(V4),
            "c2_global_reduction": str(C2_REDUCTION.relative_to(ROOT)),
            "c2_global_reduction_sha256": _sha256(C2_REDUCTION),
        },
        "cases": cases,
        "v4_h10_identification": {
            "candidate_id": "char3-V4-H10",
            "candidate_discriminant_form": v4_h10_split,
            "target_lattice": "U + E8(-2) + <-4>",
            "target_gram": _matrix_payload(target),
            "target_discriminant_form": target_split,
            "rank": 11,
            "discriminant_group_length": 9,
            "nikulin_criterion": "even indefinite, same signature and discriminant form, rank >= length(A)+2",
            "explicit_isometry_convention": "P^T * G_candidate * P = G_target",
            "explicit_isometry_matrix": _matrix_payload(explicit_isometry),
            "explicit_isometry_determinant": int(explicit_isometry.det()),
            "integrally_isometric": True,
            "maximal_reflection_chamber_lineage": "Vinberg f_2^10, 13 walls",
        },
        "c2_h10_identification": {
            "candidate_id": "char3-C2-H10",
            "candidate_gram_scalar": 2,
            "half_scaled_gram": _matrix_payload(c2_h10_half_gram),
            "target_lattice": "<3> + <-1>^10",
            "target_quadratic_form": "3*x0^2 - x1^2 - ... - x10^2",
            "target_gram": _matrix_payload(mcleod_target),
            "explicit_isometry_convention": "P^T * (G_candidate/2) * P = G_McLeod",
            "explicit_isometry_matrix": _matrix_payload(c2_h10_isometry),
            "explicit_isometry_determinant": int(c2_h10_isometry.det()),
            "integrally_isometric_after_dividing_form_by_2": True,
            "orthogonal_group_unchanged_by_scalar_2": True,
            "maximal_reflection_chamber_lineage": "McLeod f_3^10, 14 walls",
        },
        "c2_h10_196_mirrors": {
            "facet_count": len(c2_h10_mirror_records),
            "half_scaled_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(
                    c2_h10_mirror_types.items()
                )
            },
            "original_even_square_divisibility_histogram": {
                f"{2 * square},{2 * divisibility}": count
                for (square, divisibility), count in sorted(
                    c2_h10_mirror_types.items()
                )
            },
            "all_reflections_integral_involutive_isometries": (
                c2_h10_all_reflections_integral
            ),
            "mirrors": c2_h10_mirror_records,
            "subgroup_conclusion": (
                "W_196 is a literal reflection subgroup of McLeod W(f_3^10)"
            ),
        },
        "v4_h10_123_mirrors": {
            "facet_count": len(mirror_records),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(mirror_types.items())
            },
            "all_reflections_integral_involutive_isometries": all_reflections_integral,
            "mirrors": mirror_records,
            "subgroup_conclusion": "W_123 is a literal reflection subgroup of W(U+E8(-2)+<-4>)",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
