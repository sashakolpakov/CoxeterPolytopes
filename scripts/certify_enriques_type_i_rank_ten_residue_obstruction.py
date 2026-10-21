#!/usr/bin/env python3
"""Close the unique rank-ten Type-I anti-Gram residue.

The raw generated anti-lattice is negative definite and rootless, but the
half-sum glue forced by the two final downstairs walls enlarges it by index
two and creates anti-invariant (-2)-roots.  Hence the residue cannot be the
anti-algebraic lattice of a free Enriques involution.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

try:
    from scripts.certify_enriques_variable_lift_package import parse_gap_records
except ModuleNotFoundError:  # direct execution
    from certify_enriques_variable_lift_package import parse_gap_records
from coxeter_polytopes.exact.matrix import inertia


ROOT = Path(__file__).resolve().parents[1]
CENSUS = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_downstairs_anti_gram_census.json"
)
SHIMADA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_rank_ten_residue_obstruction.json"
)

CENSUS_SHA256 = "df864bed625d1411c7db46d7bbb500606ccb011297e420de8272e512ab2e22ea"
SHIMADA_SHA256 = "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("matrix is not integral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    ]


def _gcd_maximal_minors(matrix: sp.Matrix) -> int:
    rank = matrix.cols
    minors = [
        abs(int(matrix.extract(rows, range(rank)).det()))
        for rows in itertools.combinations(range(matrix.rows), rank)
    ]
    return int(sp.gcd_list(minors))


def _rank_mod_two(matrix: sp.Matrix) -> int:
    rows = [[int(value) % 2 for value in row] for row in matrix.tolist()]
    row_count = len(rows)
    column_count = len(rows[0]) if rows else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for row in range(row_count):
            if row != pivot_row and rows[row][column]:
                rows[row] = [left ^ right for left, right in zip(rows[row], rows[pivot_row])]
        pivot_row += 1
    return pivot_row


def _kernel_mod_two(matrix: sp.Matrix) -> list[list[int]]:
    rows = [[int(value) % 2 for value in row] for row in matrix.tolist()]
    row_count = len(rows)
    column_count = len(rows[0]) if rows else 0
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for row in range(row_count):
            if row != pivot_row and rows[row][column]:
                rows[row] = [left ^ right for left, right in zip(rows[row], rows[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
    free_columns = [column for column in range(column_count) if column not in pivot_columns]
    kernel = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = rows[row][free]
        kernel.append(vector)
    return kernel


def _quotient_coordinates(
    gram: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, list[int]]:
    """Return the saturated radical, generator coordinates, and quotient Gram."""

    rational_kernel = gram.nullspace()
    if len(rational_kernel) != 2:
        raise ArithmeticError("rank-ten residue no longer has radical rank two")
    primitive = []
    for vector in rational_kernel:
        denominator = sp.ilcm(*[entry.q for entry in vector])
        integral = sp.Matrix([int(entry * denominator) for entry in vector])
        primitive.append(integral / sp.gcd_list(list(integral)))
    radical = sp.Matrix.hstack(primitive[0], (primitive[0] + primitive[1]) / 2)
    if any(value.q != 1 for value in radical):
        raise ArithmeticError("saturated radical basis is not integral")
    if gram * radical != sp.zeros(12, 2) or _gcd_maximal_minors(radical) != 1:
        raise ArithmeticError("radical basis is not saturated")

    # The rows for walls 2 and 9 form a unit minor in the saturated radical.
    eliminated = [1, 8]
    retained = [index for index in range(12) if index not in eliminated]
    pivot = radical.extract(eliminated, [0, 1])
    if abs(int(pivot.det())) != 1:
        raise ArithmeticError("quotient pivot changed")
    coordinates = sp.zeros(12, 10)
    for coordinate, generator in enumerate(retained):
        coordinates[generator, coordinate] = 1
    for generator in eliminated:
        target = sp.Matrix([1 if index == generator else 0 for index in eliminated])
        coefficients = pivot.inv() * target
        quotient_row = -radical.extract(retained, [0, 1]) * coefficients
        for coordinate, value in enumerate(quotient_row):
            coordinates[generator, coordinate] = value
    quotient_gram = gram.extract(retained, retained)
    if coordinates * quotient_gram * coordinates.T != gram:
        raise ArithmeticError("quotient coordinates do not reconstruct the anti-Gram")
    return radical, coordinates, quotient_gram, retained


def analyze() -> dict[str, Any]:
    if _sha256(CENSUS) != CENSUS_SHA256:
        raise AssertionError("Type-I anti-Gram census changed")
    if _sha256(SHIMADA) != SHIMADA_SHA256:
        raise AssertionError("Shimada input changed")

    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    gram = sp.Matrix(census["rank_ten_residue"]["canonical_gram"])
    exact_inertia = inertia(_rows(gram))
    if asdict(exact_inertia) != {"positive": 0, "negative": 10, "zero": 2}:
        raise ArithmeticError("rank-ten residue inertia changed")

    record = parse_gap_records(SHIMADA)["12A"]
    walls = sp.Matrix(record["walls"])
    gram_l10 = sp.Matrix(record["GramL10"])
    if walls.shape != (12, 10) or walls * gram_l10 * walls.T != sp.Matrix(record["configuration"]):
        raise ArithmeticError("Shimada Type-I wall marking changed")

    radical, anti_coordinates, raw_anti_gram, retained = _quotient_coordinates(gram)
    first_block = raw_anti_gram[:8, :8]
    scaled_e8 = -first_block / 2
    if (
        raw_anti_gram[:8, 8:] != sp.zeros(8, 2)
        or raw_anti_gram[8:, :8] != sp.zeros(2, 8)
        or raw_anti_gram[8:, 8:] != sp.diag(-4, -4)
        or any(value.q != 1 for value in scaled_e8)
        or int(scaled_e8.det()) != 1
        or any(int(scaled_e8[index, index]) % 2 for index in range(8))
        or asdict(inertia(_rows(scaled_e8))) != {"positive": 8, "negative": 0, "zero": 0}
    ):
        raise ArithmeticError("raw anti-lattice is not E8(-2) plus two <-4> factors")
    raw_smith = _smith_diagonal(raw_anti_gram)
    if raw_smith != [2] * 8 + [4, 4] or abs(int(raw_anti_gram.det())) != 2**12:
        raise ArithmeticError("raw anti-lattice discriminant changed")
    if any(int(raw_anti_gram[index, index]) % 4 for index in range(10)) or any(
        int(2 * raw_anti_gram[left, right]) % 4
        for left in range(10)
        for right in range(left + 1, 10)
    ):
        raise ArithmeticError("raw anti-lattice norm divisibility changed")

    # The last two roots are congruent downstairs and orthogonal upstairs.
    downstairs_difference = walls.row(10) - walls.row(11)
    if any(int(value) % 2 for value in downstairs_difference):
        raise ArithmeticError("walls 11 and 12 are no longer congruent modulo two")
    if gram[10, 10] != -4 or gram[11, 11] != -4 or gram[10, 11] != 0:
        raise ArithmeticError("last two proposed anti-vectors changed")
    half_difference = (anti_coordinates.row(10) - anti_coordinates.row(11)) / 2
    half_sum = (anti_coordinates.row(10) + anti_coordinates.row(11)) / 2
    for forced_root in (half_difference, half_sum):
        if (
            (forced_root * raw_anti_gram * forced_root.T)[0] != -2
            or any(value.q not in (1, 2) for value in forced_root)
            or any(value.q != 1 for value in forced_root * raw_anti_gram)
        ):
            raise ArithmeticError("forced anti-invariant root identity changed")

    # Adjoining either forced root is the exact index-two anti-saturation.
    saturation_numerator = sp.Matrix([int(2 * value) for value in half_difference]).reshape(10, 1)
    anti_hnf = hermite_normal_form(sp.Matrix.hstack(2 * sp.eye(10), saturation_numerator))
    anti_saturated_gram = anti_hnf.T * raw_anti_gram * anti_hnf / 4
    if (
        any(value.q != 1 for value in anti_saturated_gram)
        or abs(int(anti_hnf.det())) != 2**9
        or abs(int(anti_saturated_gram.det())) != 2**10
        or _smith_diagonal(anti_saturated_gram) != [2] * 10
    ):
        raise ArithmeticError("forced anti-saturation changed")

    # Construct the full half-sum Picard overmodule as a consistency check.
    invariant_gram = 2 * gram_l10
    direct_gram = sp.diag(invariant_gram, raw_anti_gram)
    half_sum_numerators = sp.Matrix.hstack(walls, anti_coordinates)
    glue_rank = _rank_mod_two(half_sum_numerators)
    downstairs_glue_rank = _rank_mod_two(walls)
    if glue_rank != 10 or downstairs_glue_rank != 9:
        raise ArithmeticError("half-sum glue ranks changed")
    picard_hnf = hermite_normal_form(
        sp.Matrix.hstack(2 * sp.eye(20), half_sum_numerators.T)
    )
    picard_gram = picard_hnf.T * direct_gram * picard_hnf / 4
    if (
        any(value.q != 1 for value in picard_gram)
        or any(int(picard_gram[index, index]) % 2 for index in range(20))
        or int(picard_gram.det()) != -4
        or _smith_diagonal(picard_gram) != [1] * 18 + [2, 2]
        or asdict(inertia(_rows(picard_gram))) != {"positive": 1, "negative": 19, "zero": 0}
    ):
        raise ArithmeticError("rank-twenty half-sum Picard lattice changed")

    discriminant_kernel = _kernel_mod_two(picard_gram)
    if len(discriminant_kernel) != 2:
        raise ArithmeticError("Picard discriminant kernel changed")
    discriminant_values = []
    for coefficients in ((1, 0), (0, 1), (1, 1)):
        numerator = sp.Matrix(
            [
                sum(coefficients[index] * discriminant_kernel[index][column] for index in range(2)) % 2
                for column in range(20)
            ]
        )
        value = sp.Rational((numerator.T * picard_gram * numerator)[0], 4) % 2
        discriminant_values.append(str(value))
    if sorted(discriminant_values) != ["1", "3/2", "3/2"]:
        raise ArithmeticError("Picard discriminant form is not -q_(A1+A1)")

    # A primitive abstract embedding of the forced saturation in
    # U + U(2) + E8(-2): the two A1 roots live in U + U(2).
    u_u2_gram = sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 2],
            [0, 0, 2, 0],
        ]
    )
    a1_embedding = sp.Matrix(
        [
            [1, 1],
            [-1, 1],
            [0, 1],
            [0, -1],
        ]
    )
    if (
        a1_embedding.T * u_u2_gram * a1_embedding != sp.diag(-2, -2)
        or _gcd_maximal_minors(a1_embedding) != 1
    ):
        raise ArithmeticError("primitive A1^2 embedding in U+U(2) changed")

    checks = {
        "unique_residue_replayed_with_inertia_0_10_2": True,
        "saturated_integral_radical_and_rank_ten_quotient_reconstructed": True,
        "raw_generated_lattice_is_E8_minus2_plus_two_minus4_and_rootless": True,
        "walls_11_12_are_congruent_mod_2_and_their_anti_vectors_are_orthogonal": True,
        "compatible_half_sum_glue_forces_an_anti_invariant_minus2_root": True,
        "forced_anti_saturation_is_index_2_E8_minus2_plus_A1_minus1_squared": True,
        "rank_20_half_sum_Picard_lattice_has_complement_A1_squared": True,
        "primitive_embedding_exists_but_rootless_Enriques_period_condition_fails": True,
        "L3_root_regular_invariant_ample_class_is_impossible_so_L4_is_not_reached": True,
    }

    return {
        "schema": "enriques_type_i_rank_ten_residue_obstruction_v1",
        "result": {
            "decision": "IMPOSSIBLE_GOOD_ENRIQUES_LIFT_PACKAGE",
            "reason": (
                "half-sum glue for the congruent final two downstairs walls forces "
                "an anti-invariant square-minus-two vector"
            ),
            "type_i_downstairs_projection_consequence": (
                "the already-certified rank-nine geometric orbit is the only "
                "negative-semidefinite sheet-switch orbit that can pass rootlessness"
            ),
            "new_abstract_polytope_claimed": False,
        },
        "pinned_inputs": {
            "census_path": str(CENSUS.relative_to(ROOT)),
            "census_sha256": CENSUS_SHA256,
            "shimada_path": str(SHIMADA.relative_to(ROOT)),
            "shimada_sha256": SHIMADA_SHA256,
        },
        "raw_generated_anti_lattice": {
            "rank": 10,
            "inertia": {"positive": 0, "negative": 10, "zero": 0},
            "isometry_type": "E8(-2) + <-4> + <-4>",
            "determinant_absolute": 2**12,
            "smith_diagonal": raw_smith,
            "discriminant_group": "(Z/2)^8 + (Z/4)^2",
            "discriminant_form": "q_E8(-2) + <-1/4> + <-1/4>",
            "all_norms_divisible_by_4": True,
            "rootless": True,
            "saturated_radical_basis_in_Z12_columns": _rows(radical),
            "quotient_retained_wall_indices_one_based": [index + 1 for index in retained],
            "quotient_gram": _rows(raw_anti_gram),
        },
        "forced_glue_root": {
            "downstairs_relation": "r_11-r_12 is in 2*L10",
            "downstairs_half_difference_L10_coordinates": [
                int(value / 2) for value in downstairs_difference
            ],
            "anti_pairings": {
                "delta_11_square": -4,
                "delta_12_square": -4,
                "delta_11_delta_12": 0,
            },
            "identity_for_arbitrary_sheet_signs": (
                "(eps11*delta11-eps12*delta12)/2 = "
                "b11-b12-(u11-u12)/2 lies in the anti-invariant lattice N"
            ),
            "forced_vector_square": -2,
            "both_possible_sign_classes_are_roots": True,
        },
        "compatible_anti_saturation": {
            "index_over_raw_generated_lattice": 2,
            "isometry_type": "E8(-2) + A1(-1)^2",
            "determinant_absolute": 2**10,
            "smith_diagonal": [2] * 10,
            "discriminant_form": "q_E8(-2) + <-1/2> + <-1/2>",
            "contains_minus_two_roots": True,
            "primitive_embedding_in_N_exists": True,
            "explicit_A1_squared_embedding_columns_in_U_plus_U2": _rows(a1_embedding),
        },
        "rank_twenty_period_edge": {
            "half_sum_glue_rank_over_F2": glue_rank,
            "downstairs_glue_image_rank_over_F2": downstairs_glue_rank,
            "anti_saturation_kernel_dimension": glue_rank - downstairs_glue_rank,
            "picard_rank": 20,
            "picard_signature": [1, 19],
            "picard_determinant": -4,
            "picard_smith_diagonal": [1] * 18 + [2, 2],
            "picard_discriminant_nonzero_q_values_mod_2Z": sorted(discriminant_values),
            "positive_transcendental_complement": "A1 + A1",
            "ordinary_K3_period_lattice_gate": "passes",
            "Enriques_period_root_hyperplane_avoidance": "fails identically",
            "explanation": (
                "the rank-two positive period plane is orthogonal to the entire "
                "anti-algebraic saturation, including its forced (-2)-roots"
            ),
        },
        "good_package_gates": {
            "L1_pairwise_anti_gram_and_local_half_sum_congruences": "passes",
            "L2_negative_definiteness": "passes after quotient",
            "L2_rootlessness_of_saturated_anti_lattice": "FAILS",
            "L3_invariant_root_regular_ample_class": "FAILS_FOR_EVERY_INVARIANT_CLASS",
            "L4_irreducibility": "not_reached_and_irrelevant_after_L2_L3_failure",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived certificate differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
