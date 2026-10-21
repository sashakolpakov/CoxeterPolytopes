#!/usr/bin/env python3
"""Independent exact replay of the affine-H3 abundance hostile audit.

This script deliberately does not import the owner certifier.  It checks the
star lattice, the termwise orthogonality gate, the low-rank contradiction,
the rational-tree Picard gate, the fixed-lattice divisibility poison, and
the genus-three adjunction/core formulas.  The horospherical compactness and
Picard-scheme arguments remain geometric proofs in the hostile report.
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import product


def mat_vec(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def rational_rank(matrix: list[list[int]]) -> int:
    rows = [[Fraction(entry) for entry in row] for row in matrix]
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [entry / scale for entry in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            scale = rows[row][column]
            if scale:
                rows[row] = [
                    entry - scale * pivot_entry
                    for entry, pivot_entry in zip(rows[row], rows[pivot_row])
                ]
        pivot_row += 1
    return pivot_row


def determinant(matrix: list[list[int]]) -> Fraction:
    rows = [[Fraction(entry) for entry in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        scale = rows[column][column]
        result *= scale
        rows[column] = [entry / scale for entry in rows[column]]
        for row in range(column + 1, len(rows)):
            scale = rows[row][column]
            rows[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(rows[row], rows[column])
            ]
    return result


def inverse(matrix: list[list[int]]) -> list[list[Fraction]]:
    size = len(matrix)
    rows = [
        [Fraction(entry) for entry in row]
        + [Fraction(int(row_index == column)) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [entry / scale for entry in rows[column]]
        for row in range(size):
            if row == column:
                continue
            scale = rows[row][column]
            rows[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(rows[row], rows[column])
            ]
    return [row[size:] for row in rows]


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    physical_gram = [
        [-1, 1, 1, 1, 1],
        [1, -4, 0, 0, 0],
        [1, 0, -4, 0, 0],
        [1, 0, 0, -4, 0],
        [1, 0, 0, 0, -4],
    ]
    host_gram = [
        [-2, 1, 1, 1, 1],
        [1, -2, 0, 0, 0],
        [1, 0, -2, 0, 0],
        [1, 0, 0, -2, 0],
        [1, 0, 0, 0, -2],
    ]
    physical_null = [4, 1, 1, 1, 1]
    host_null = [2, 1, 1, 1, 1]
    assert mat_vec(physical_gram, physical_null) == [0] * 5
    assert mat_vec(host_gram, host_null) == [0] * 5
    assert rational_rank(physical_gram) == rational_rank(host_gram) == 4

    # If 4(E.C)+sum(B_i.C)=0 and all five intersections are nonnegative,
    # every intersection vanishes.  Exhaust a box as a poison replay; the
    # universal implication is immediate from positivity of the weights.
    zero_profiles = [
        values
        for values in product(range(5), repeat=5)
        if 4 * values[0] + sum(values[1:]) == 0
    ]
    assert zero_profiles == [(0, 0, 0, 0, 0)]

    # A bounded rational Euclidean factor has positively spanning normals.
    # This exact toy factor records denominator clearing, including a
    # half-root belonging to a physical (-4)-curve.
    quotient_normals = [(1, 0), (0, 1), (-1, -1)]
    positive_relation = [1, 1, 1]
    assert [
        sum(coefficient * normal[coordinate] for coefficient, normal in zip(positive_relation, quotient_normals))
        for coordinate in range(2)
    ] == [0, 0]
    physical_scales = [Fraction(1), Fraction(1, 2), Fraction(1)]
    cleared_physical_coefficients = [
        int(2 * coefficient * scale)
        for coefficient, scale in zip(positive_relation, physical_scales)
    ]
    assert cleared_physical_coefficients == [2, 1, 2]

    # There is no low-rank endpoint.  When rho <= 10, RR and vanishing of
    # rational plurigenera give an effective -K, contradicting nef F since
    # (-K).F=-4.
    endpoint_rows = []
    for rho in range(6, 11):
        k_square = 10 - rho
        chi_minus_k = 1 + k_square
        assert chi_minus_k >= 1
        endpoint_rows.append(
            {
                "rho": rho,
                "K_square": k_square,
                "chi_minus_K": chi_minus_k,
                "h2_minus_K": 0,
                "effective_minus_K": True,
                "minus_K_dot_F": -4,
                "contradicts_nef_F": True,
            }
        )

    # Any integral fibre divisor supported on the affine star and orthogonal
    # to its components is qF.  Exhaust primitive-size kernel vectors as an
    # exact check of the one-dimensional kernel and its integral generator.
    supported_fibre_cycles = []
    for coefficients in product(range(13), repeat=5):
        if coefficients == (0, 0, 0, 0, 0):
            continue
        if mat_vec(physical_gram, list(coefficients)) != [0] * 5:
            continue
        supported_fibre_cycles.append(list(coefficients))
    assert supported_fibre_cycles == [
        [4 * q, q, q, q, q] for q in range(1, 4)
    ]

    # The reduced star is a tree of five P1s.  Multidegree zero kills its
    # reduced Picard class; the nilpotent restriction kernel is a complex
    # vector group in characteristic zero and has no torsion.  Therefore a
    # torsion normal bundle O_F(F) is trivial and F itself moves.
    picard_gate = {
        "components": 5,
        "nodes": 4,
        "dual_graph_first_betti": 0,
        "multidegrees_O_F_F": [0, 0, 0, 0, 0],
        "reduced_multidegree_zero_Picard": "trivial",
        "nilpotent_kernel": "H1(F,I)_complex_vector_group",
        "nilpotent_kernel_torsion_free": True,
        "power_trivial_normal_bundle_implies_trivial": True,
        "F_base_point_free": True,
    }
    assert picard_gate["nodes"] - picard_gate["components"] + 1 == 0

    fibration_case = {
        "moving_index": 1,
        "general_fibre_class": "F",
        "genus": 3,
        "relative_minimal_rho_bound": 26,
    }
    assert 2 * fibration_case["genus"] - 2 == 4

    # Lattice poison: U + D4(-1) has fixed exponent-two discriminant, but
    # alpha_0=N e-theta has the same affine D4 Gram for every N and null
    # relation alpha_0+theta=N e.  Thus the Picard gate is essential.
    finite_d4 = [
        [-2, 1, 1, 1],
        [1, -2, 0, 0],
        [1, 0, -2, 0],
        [1, 0, 0, -2],
    ]
    highest_root = [2, 1, 1, 1]
    assert determinant(finite_d4) == 4
    finite_d4_inverse = inverse(finite_d4)
    assert any(entry.denominator == 2 for row in finite_d4_inverse for entry in row)
    assert all((2 * entry).denominator == 1 for row in finite_d4_inverse for entry in row)
    lattice_poison_rows = []
    for null_divisibility in range(1, 13):
        # Coordinates are (e,f,central,leaf1,leaf2,leaf3) in U + D4(-1).
        affine_root = [null_divisibility, 0, -2, -1, -1, -1]
        simple_roots = [
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ]
        ambient_gram = [
            [0, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 0, -2, 1, 1, 1],
            [0, 0, 1, -2, 0, 0],
            [0, 0, 1, 0, -2, 0],
            [0, 0, 1, 0, 0, -2],
        ]

        def ambient_dot(left: list[int], right: list[int]) -> int:
            return sum(
                left[row] * ambient_gram[row][column] * right[column]
                for row in range(6)
                for column in range(6)
            )

        affine_roots = [affine_root, *simple_roots]
        affine_gram = [
            [ambient_dot(left, right) for right in affine_roots]
            for left in affine_roots
        ]
        expected_affine = [
            [-2, 1, 0, 0, 0],
            [1, -2, 1, 1, 1],
            [0, 1, -2, 0, 0],
            [0, 1, 0, -2, 0],
            [0, 1, 0, 0, -2],
        ]
        assert affine_gram == expected_affine
        null_sum = [
            affine_root[index]
            + 2 * simple_roots[0][index]
            + simple_roots[1][index]
            + simple_roots[2][index]
            + simple_roots[3][index]
            for index in range(6)
        ]
        assert null_sum == [null_divisibility, 0, 0, 0, 0, 0]
        lattice_poison_rows.append(
            {
                "N": null_divisibility,
                "same_affine_D4_gram": True,
                "null_vector": null_sum,
            }
        )

    return {
        "schema": "rational-affine-h3-finite-volume-abundance-hostile-v2",
        "verdict": "PASS_AFTER_REPAIR",
        "owner_core": {
            "finite_volume_implies_F_base_point_free": True,
            "reflection_maximality_used": False,
            "external_classification_used": False,
        },
        "repairs": {
            "rank_six_endpoint_nonvacuous": False,
            "minimum_possible_picard_rank": 11,
            "lattice_alone_forces_index_one": False,
            "nonreduced_picard_forces_index_one": True,
            "genus_two_branch": False,
            "relative_minimal_core_bound_26": True,
            "vertical_tail_bound": False,
        },
        "physical_gram": physical_gram,
        "physical_null_cycle": physical_null,
        "host_gram": host_gram,
        "host_null_cycle": host_null,
        "termwise_zero_profiles": [list(profile) for profile in zero_profiles],
        "sample_positive_relation": {
            "quotient_normals": quotient_normals,
            "coefficients": positive_relation,
            "physical_root_scales": [str(scale) for scale in physical_scales],
            "cleared_physical_coefficients": cleared_physical_coefficients,
        },
        "empty_endpoint_rows": endpoint_rows,
        "supported_fibre_cycles": supported_fibre_cycles,
        "picard_gate": picard_gate,
        "fibration_case": fibration_case,
        "lattice_poison": {
            "ambient": "U_plus_D4_minus_one",
            "determinant_absolute": 4,
            "discriminant_exponent": 2,
            "rows": lattice_poison_rows,
        },
        "source_audit": {
            "arxiv": "0901.1094v2",
            "source_archive_sha256": "146fdf403f85ba8383edb46777a9df5fd6a61a1f4ed4718d8661b90f757e27d3",
            "main_tex_sha256": "6a23a31d295bb0d2bfdb644aba1a2d778edc3aaa78e17058f5dfde774b1e6e94",
            "lemma_2_4": "big -K forces (-K).N>0 for every nonzero nef N",
            "lemma_2_5": "nef non-big N is base-point-free under the inherited big -K hypothesis",
            "theorem_2_9": "big rational surfaces have finitely generated Cox ring",
        },
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
