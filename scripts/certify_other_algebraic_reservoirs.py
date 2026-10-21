#!/usr/bin/env python3
"""Exact lightweight certificates for the non-K3 algebraic reservoir audit.

This script deliberately does not enumerate Vinberg roots.  It freezes three
small exact calculations used by the accompanying theorem note:

* the curve classes on four flex-cubic blow-ups are Vinberg's
  ``P10`` through ``P13`` simple roots in the odd unimodular lattices; the
  last two leave the cubic after point eleven and follow a free exceptional
  tail;
* the Reed--Muller description of Borcherds' ``H17/960`` chamber has 28
  anisotropic discriminant classes, with 32 long roots in every class; and
* the one-place-at-infinity cone theorem has a sharp integral-root ceiling
  in dimension eleven (the smooth cubic with a two-point fibre chain); and
* the standard plane-blow-up Cremona simplex acquires a Lorentzian proper
  vertex subdiagram as soon as its dimension exceeds nine.

The algebro-geometric completeness arguments (genus/intersection inequality,
nef-cone equality, and primitive-embedding glue obstruction) are proofs in
``results/other_algebraic_reservoirs_audit.md``.  The output here certifies
their finite combinatorial inputs without a floating tolerance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from functools import reduce
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data" / "enumerations" / "other_algebraic_reservoirs_certificate.json"


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _inertia_symmetric(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return exact inertia by an LDL-style symmetric congruence recursion."""

    matrix = sp.Matrix(matrix)
    if matrix != matrix.T:
        raise ValueError("inertia requires a symmetric matrix")

    positive = negative = zero = 0
    while matrix.rows:
        size = matrix.rows
        pivot = next((i for i in range(size) if matrix[i, i] != 0), None)
        if pivot is not None:
            if pivot:
                permutation = list(range(size))
                permutation[0], permutation[pivot] = permutation[pivot], permutation[0]
                matrix = matrix.extract(permutation, permutation)
            value = sp.Rational(matrix[0, 0])
            if value > 0:
                positive += 1
            else:
                negative += 1
            if size == 1:
                break
            column = matrix[1:, 0]
            matrix = sp.simplify(matrix[1:, 1:] - column * column.T / value)
            continue

        pair = next(
            ((i, j) for i in range(size) for j in range(i + 1, size) if matrix[i, j] != 0),
            None,
        )
        if pair is None:
            zero += size
            break
        i, j = pair
        permutation = [i, j] + [k for k in range(size) if k not in (i, j)]
        matrix = matrix.extract(permutation, permutation)
        block = matrix[:2, :2]
        # A zero-diagonal nonzero symmetric 2-by-2 block has inertia (1,1).
        positive += 1
        negative += 1
        if size == 2:
            break
        cross = matrix[2:, :2]
        matrix = sp.simplify(matrix[2:, 2:] - cross * block.inv() * cross.T)
    return positive, negative, zero


def _extreme_ray_census(
    intersection: sp.Matrix, root_columns: sp.Matrix
) -> dict[str, object]:
    """Enumerate a small rational polyhedral cone's rays exactly."""

    ambient_rank = intersection.rows
    inequalities = root_columns.T * intersection
    rays: dict[tuple[int, ...], int] = {}
    for active in combinations(range(root_columns.cols), ambient_rank - 1):
        active_matrix = inequalities[list(active), :]
        if active_matrix.rank() != ambient_rank - 1:
            continue
        vector = active_matrix.nullspace()[0]
        values = inequalities * vector
        if all(value >= 0 for value in values):
            pass
        elif all(value <= 0 for value in values):
            vector = -vector
        else:
            continue

        denominator = sp.ilcm(*(sp.denom(value) for value in vector))
        integral = [int(value * denominator) for value in vector]
        content = reduce(gcd, (abs(value) for value in integral if value), 0)
        primitive = tuple(value // content for value in integral)
        if next(value for value in primitive if value) < 0:
            primitive = tuple(-value for value in primitive)
        square = int((sp.Matrix(primitive).T * intersection * sp.Matrix(primitive))[0])
        rays[primitive] = square

    return {
        "total": len(rays),
        "finite": sum(square > 0 for square in rays.values()),
        "ideal": sum(square == 0 for square in rays.values()),
        "hyperideal": sum(square < 0 for square in rays.values()),
        "primitive_rays": [
            {"vector": list(vector), "square": square}
            for vector, square in sorted(rays.items())
        ],
    }


def _pn_surface_certificate(number_of_points: int) -> dict[str, object]:
    if number_of_points not in (10, 11, 12, 13):
        raise ValueError("the certified flex-cubic construction is used for 10 <= n <= 13")
    intersection = sp.diag(1, *([-1] * number_of_points))
    roots: list[list[int]] = []
    labels: list[str] = []

    for i in range(1, number_of_points):
        vector = [0] * (number_of_points + 1)
        vector[i] = 1
        vector[i + 1] = -1
        roots.append(vector)
        labels.append(f"E{i}-E{i + 1}")

    last_exceptional = [0] * (number_of_points + 1)
    last_exceptional[number_of_points] = 1
    roots.append(last_exceptional)
    labels.append(f"E{number_of_points}")

    flex_line = [1, -1, -1, -1] + [0] * (number_of_points - 3)
    roots.append(flex_line)
    labels.append("H-E1-E2-E3")

    cubic_points = min(number_of_points, 11)
    cubic = [3] + [-1] * cubic_points + [0] * (number_of_points - cubic_points)
    roots.append(cubic)
    labels.append(
        f"3H-E1-...-E{cubic_points}"
        + ("=-K_X" if number_of_points <= 11 else "=C_tilde")
    )

    columns = sp.Matrix(roots).T
    gram = columns.T * intersection * columns
    norms = [int(gram[i, i]) for i in range(gram.rows)]
    positive_intersections = {
        f"{labels[i]} :: {labels[j]}": int(gram[i, j])
        for i in range(gram.rows)
        for j in range(i + 1, gram.cols)
        if gram[i, j] > 0
    }
    unexpected = [
        (i, j, int(gram[i, j]))
        for i in range(gram.rows)
        for j in range(i + 1, gram.cols)
        if gram[i, j] < 0
    ]

    # Deleting the anticanonical root gives a lattice basis of determinant one.
    basis = columns[:, : number_of_points + 1]
    expected_norms = (
        Counter({-2: 10, -1: 2})
        if number_of_points == 10
        else Counter({-2: number_of_points + 1, -1: 1})
    )
    ray_census = _extreme_ray_census(intersection, columns)
    expected_ray_census = {
        10: {"total": 19, "finite": 17, "ideal": 2, "hyperideal": 0},
        11: {"total": 28, "finite": 26, "ideal": 2, "hyperideal": 0},
        12: {"total": 37, "finite": 35, "ideal": 2, "hyperideal": 0},
        13: {"total": 46, "finite": 43, "ideal": 3, "hyperideal": 0},
    }[number_of_points]
    checks = {
        "ambient_signature": _inertia_symmetric(intersection)
        == (1, number_of_points, 0),
        "ambient_determinant_abs_one": abs(int(intersection.det())) == 1,
        "root_rank_is_full": columns.rank() == number_of_points + 1,
        "roots_before_cubic_are_integral_basis": abs(int(basis.det())) == 1,
        "expected_norm_histogram": Counter(norms) == expected_norms,
        "all_distinct_curve_intersections_nonnegative": not unexpected,
        "cubic_class_exact": roots[-1]
        == [3] + [-1] * cubic_points + [0] * (number_of_points - cubic_points),
        "exact_extreme_ray_census": all(
            ray_census[key] == value for key, value in expected_ray_census.items()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"P{number_of_points} surface certificate changed: {checks}")

    known_face_data = {
        10: {"finite_vertices": 17, "cusps": 2},
        11: {"finite_vertices": 26, "cusps": 2},
        12: {"finite_vertices": 35, "cusps": 2},
        13: {"finite_vertices": 43, "cusps": 3},
    }[number_of_points]

    if number_of_points <= 11:
        surface = (
            f"blow up a smooth plane cubic successively at {number_of_points} "
            "infinitely-near points along the cubic above a flex"
        )
        completeness = (
            "D_tilde^2 = [(d-1)(d-2)-sum m_i(m_i-1)] "
            "+ [3d-sum m_i] - 2 >= -2 for D != cubic"
        )
    else:
        tail_length = number_of_points - 11
        surface = (
            "blow up eleven successive points along a smooth cubic above a flex, "
            f"then a free exceptional tail of length {tail_length} away from the cubic"
        )
        completeness = (
            "for A_j=sum_(i<=j)m_i and y=max(0,(A_n-3d)/2), every other "
            "plane curve is xL+yC+sum mu_i R_i with x=d-3y and "
            "mu_j=x min(j,3)+y min(j,11)-A_j >= 0"
        )

    return {
        "surface": surface,
        "lattice": (
            f"I_(1,{number_of_points}) with basis H,E1,...,E{number_of_points}"
        ),
        "hyperbolic_dimension": number_of_points,
        "simple_curve_labels": labels,
        "simple_curve_columns": _matrix_payload(columns),
        "intersection_gram": _matrix_payload(gram),
        "norm_histogram": {str(key): value for key, value in sorted(Counter(norms).items())},
        "nonzero_off_diagonal_intersections": positive_intersections,
        "known_chamber": (
            f"Vinberg P{number_of_points}, {number_of_points + 2} facets, "
            "finite-volume noncompact"
        ),
        "known_exact_face_data": known_face_data,
        "exact_cone_extreme_rays": ray_census,
        "completeness_identity": completeness,
        "checks": checks,
        "all_checks_pass": True,
    }


def _rm_first_order_code() -> set[int]:
    words: set[int] = set()
    for linear_mask in range(16):
        for constant in range(2):
            word = 0
            for point in range(16):
                bit = ((linear_mask & point).bit_count() & 1) ^ constant
                word |= bit << point
            words.add(word)
    return words


def _dual_code(code: Iterable[int]) -> set[int]:
    # Constant one and the four coordinate functions generate RM(1,4).
    code = set(code)
    basis = [next(word for word in code if word.bit_count() == 16)]
    coordinate_words = []
    for coordinate in range(4):
        word = sum(((point >> coordinate) & 1) << point for point in range(16))
        coordinate_words.append(word)
    basis.extend(coordinate_words)
    return {
        word
        for word in range(1 << 16)
        if all(((word & generator).bit_count() & 1) == 0 for generator in basis)
    }


def _reed_muller_certificate() -> dict[str, object]:
    code = _rm_first_order_code()
    dual = _dual_code(code)
    code_weights = Counter(word.bit_count() for word in code)
    dual_weights = Counter(word.bit_count() for word in dual)

    unseen = set(dual)
    cosets: list[set[int]] = []
    while unseen:
        representative = min(unseen)
        coset = {representative ^ codeword for codeword in code}
        if not coset <= dual:
            raise ArithmeticError("RM(1,4) is no longer self-orthogonal")
        cosets.append(coset)
        unseen -= coset

    q_histogram: Counter[int] = Counter()
    coset_size_histogram: Counter[int] = Counter()
    anisotropic_weight_histograms: Counter[tuple[tuple[int, int], ...]] = Counter()
    for coset in cosets:
        q_values = {(word.bit_count() // 2) & 1 for word in coset}
        if len(q_values) != 1:
            raise ArithmeticError("weight/2 modulo 2 is not constant on an RM coset")
        q_value = next(iter(q_values))
        q_histogram[q_value] += 1
        coset_size_histogram[len(coset)] += 1
        if q_value:
            weight_histogram = tuple(sorted(Counter(word.bit_count() for word in coset).items()))
            anisotropic_weight_histograms[weight_histogram] += 1

    long_words = {word for word in dual if word.bit_count() % 4 == 2}
    checks = {
        "rm_dimension_five": len(code) == 32,
        "rm_weight_enumerator": code_weights == Counter({8: 30, 0: 1, 16: 1}),
        "dual_dimension_eleven": len(dual) == 2048,
        "dual_weight_enumerator": dual_weights
        == Counter({8: 870, 6: 448, 10: 448, 4: 140, 12: 140, 0: 1, 16: 1}),
        "discriminant_has_sixty_four_classes": len(cosets) == 64,
        "all_cosets_have_thirty_two_words": coset_size_histogram == Counter({32: 64}),
        "plus_type_q_histogram": q_histogram == Counter({0: 36, 1: 28}),
        "long_word_count": len(long_words) == 896,
        "each_anisotropic_class_has_thirty_two_long_words": (
            anisotropic_weight_histograms
            == Counter({((6, 16), (10, 16)): 28})
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Reed--Muller certificate changed: {checks}")

    return {
        "lattice": "II_(1,17)(2^{+6}) = (18,6,0)",
        "signature": [1, 17],
        "determinant_absolute": 64,
        "simple_roots": {"norm_-2": 64, "norm_-4": 896, "total": 960},
        "rm_weight_enumerator": {str(k): v for k, v in sorted(code_weights.items())},
        "dual_weight_enumerator": {str(k): v for k, v in sorted(dual_weights.items())},
        "discriminant_quadratic_histogram": {str(k): v for k, v in sorted(q_histogram.items())},
        "long_roots_per_anisotropic_discriminant_class": 32,
        "primitive_K3_3_embedding_glue": {
            "ambient_lattice": "Lambda_K3[3] = 3U + 2E8(-1) + <-4>",
            "ambient_discriminant_group": "Z/4",
            "orthogonal_complement_rank": 5,
            "forced_glue_projection_dimension_in_A_L": 5,
            "orthogonal_discriminant_classes_available_for_divisibility_two": 1,
            "maximum_compatible_norm_minus_four_walls": 32,
            "minimum_incompatible_norm_minus_four_walls": 864,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def _cremona_roots(number_of_points: int) -> tuple[sp.Matrix, sp.Matrix]:
    intersection = sp.diag(1, *([-1] * number_of_points))
    roots: list[list[int]] = []
    alpha_zero = [1, -1, -1, -1] + [0] * (number_of_points - 3)
    roots.append(alpha_zero)
    for i in range(1, number_of_points):
        vector = [0] * (number_of_points + 1)
        vector[i] = 1
        vector[i + 1] = -1
        roots.append(vector)
    columns = sp.Matrix(roots).T
    return columns, columns.T * intersection * columns


def _cremona_ceiling_certificate() -> dict[str, object]:
    roots10, gram10 = _cremona_roots(10)
    roots11, gram11 = _cremona_roots(11)
    # In the n=11 diagram, delete the terminal alpha_10.  The remaining
    # principal block is literally the n=10 Lorentzian T_(2,3,7) diagram.
    inherited = gram11[:-1, :-1]
    checks = {
        "n10_root_rank": roots10.rank() == 10,
        "n10_lorentzian_inertia": _inertia_symmetric(gram10) == (1, 9, 0),
        "n11_root_rank": roots11.rank() == 11,
        "n11_lorentzian_inertia": _inertia_symmetric(gram11) == (1, 10, 0),
        "n11_has_proper_lorentzian_vertex_subdiagram": inherited == gram10,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Cremona ceiling certificate changed: {checks}")
    return {
        "simple_diagram": "T_(2,3,n-3) on K_X^perp for the blow-up of P2 at n points",
        "n_equals_10": "T_(2,3,7), finite-volume noncompact simplex in H9",
        "n_at_least_11": (
            "deleting the long-arm terminal leaves a Lorentzian T_(2,3,7) "
            "proper vertex subdiagram, hence the canonical simplex has a hyperideal vertex"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def _pencil_at_infinity_root_ceiling_certificate() -> dict[str, object]:
    """Freeze the elementary numerical part of the Galindo et al. ceiling.

    For a smooth degree-d curve with one place at infinity, the pencil
    <C,dL> has d^2 free base points and L follows the first d.  If k fibre
    branches of lengths s_j are then attached at the terminal exceptional,
    the potentially obstructing negative curves have squares

        L^2 = 1-d, E_terminal^2 = -1-k, C_j^2 = -s_j.

    In the unimodular plane-blow-up lattice these primitive classes define
    integral reflections only when their negative square is -1 or -2.
    """

    admissible_smooth_cases: list[dict[str, object]] = []
    for degree in range(1, 8):
        base_points = degree * degree
        for branches in range(0, 5):
            length_choices = [()] if branches == 0 else []
            if branches:
                # Equal lengths suffice for the bounded audit; the formula
                # below records the componentwise condition as well.
                length_choices = [(length,) * branches for length in range(1, 6)]
            for lengths in length_choices:
                negative_squares = []
                if 1 - degree < 0:
                    negative_squares.append(1 - degree)
                negative_squares.append(-1 - branches)
                negative_squares.extend(-length for length in lengths)
                integral_root_facets = all(square in (-1, -2) for square in negative_squares)
                if integral_root_facets:
                    admissible_smooth_cases.append({
                        "degree": degree,
                        "base_points": base_points,
                        "branches": branches,
                        "branch_lengths": list(lengths),
                        "hyperbolic_dimension": base_points + sum(lengths),
                    })

    checks = {
        "smooth_cubic_has_nine_base_points": 3 * 3 == 9,
        "infinity_line_on_smooth_cubic_has_square_minus_two": 1 - 3 == -2,
        "at_most_one_added_fibre_branch": all(
            case["branches"] <= 1 for case in admissible_smooth_cases
        ),
        "each_added_fibre_chain_has_length_at_most_two": all(
            all(length <= 2 for length in case["branch_lengths"])
            for case in admissible_smooth_cases
        ),
        "maximum_dimension_is_eleven": max(
            case["hyperbolic_dimension"] for case in admissible_smooth_cases
        ) == 11,
        "unique_dimension_eleven_case_is_cubic_plus_two": [
            case
            for case in admissible_smooth_cases
            if case["hyperbolic_dimension"] == 11
        ] == [{
            "degree": 3,
            "base_points": 9,
            "branches": 1,
            "branch_lengths": [2],
            "hyperbolic_dimension": 11,
        }],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"pencil-at-infinity ceiling changed: {checks}")
    return {
        "source": (
            "Galindo--Monserrat--Moreno-Avila--Perez-Callejo, "
            "arXiv:2602.15703, Theorems 3 and 7"
        ),
        "smooth_base_cluster": {
            "number_of_base_points": "d^2",
            "all_base_points_free": True,
            "line_at_infinity_last_base_point": "d",
        },
        "negative_generator_squares": {
            "line_at_infinity": "1-d",
            "terminal_base_exceptional_after_k_branches": "-1-k",
            "selected_fibre_with_chain_length_s_j": "-s_j",
            "nonterminal_free_chain_exceptionals": "-2",
            "terminal_branch_exceptionals": "-1",
        },
        "integral_root_conditions": {
            "degree": "d <= 3",
            "number_of_added_branches": "k <= 1",
            "length_of_the_unique_added_branch": "s_1 <= 2",
        },
        "sharp_maximum": {
            "degree": 3,
            "base_points": 9,
            "added_branch_length": 2,
            "hyperbolic_dimension": 11,
            "chamber": "P11",
        },
        "singular_branch_obstruction": (
            "a satellite base point makes an older strict exceptional curve "
            "have square at most -3"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema": "coxeter-polytopes.other-algebraic-reservoirs.v2",
        "p10_rational_surface": _pn_surface_certificate(10),
        "p11_rational_surface": _pn_surface_certificate(11),
        "p12_rational_surface": _pn_surface_certificate(12),
        "p13_rational_surface": _pn_surface_certificate(13),
        "borcherds_h17_reed_muller": _reed_muller_certificate(),
        "pencil_at_infinity_root_ceiling": _pencil_at_infinity_root_ceiling_certificate(),
        "plane_blowup_cremona_ceiling": _cremona_ceiling_certificate(),
        "toric_surface_obstruction": {
            "rho_at_least_11_implies_boundary_components_at_least": 13,
            "identity": "sum D_i^2 = 12 - 3r",
            "if_all_negative_squares_at_least_minus_two": "sum D_i^2 >= -2r",
            "contradiction_for": "r >= 13",
            "conclusion": (
                "a smooth complete toric surface of Picard rank at least 11 has an invariant "
                "negative curve of square at most -3, whose reflection is nonintegral"
            ),
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "payload_sha256": payload["payload_sha256"],
        "all_checks_pass": all(
            payload[key]["all_checks_pass"]
            for key in (
                "p10_rational_surface",
                "p11_rational_surface",
                "p12_rational_surface",
                "p13_rational_surface",
                "borcherds_h17_reed_muller",
                "pencil_at_infinity_root_ceiling",
                "plane_blowup_cremona_ceiling",
            )
        ),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
