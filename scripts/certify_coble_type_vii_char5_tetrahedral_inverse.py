#!/usr/bin/env python3
"""Certify a strict P-first plane realization of the Type-VII H10/72 datum.

The certifier starts from the frozen 72-wall lattice, finds all five maximal
ten-curve blowdown frames, and derives their 4+6 line/conic/quartic design.
It then realizes that design over F5 without importing the known existence of
a Type-VII surface.  Exact finite-field, Singular, cddlib, and lattice checks
return the constructed strong Coble surface to the same marked H10/72 cone.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp


PRIME = 5
ROOT = Path(__file__).resolve().parents[1]
RAW_DESIGN = ROOT / "data/enumerations/coble_remaining_seven_contraction_designs.json"
RAW_POLYTOPE = (
    ROOT / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json"
)
OUTPUT = (
    ROOT / "data/enumerations/coble_type_vii_char5_tetrahedral_inverse.json"
)
SCHEMA = "coxeter.coble_type_vii_char5_tetrahedral_inverse.v1"

Point = tuple[int, int, int]
Pair = tuple[int, int]
PointLabel = tuple[str, int]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inverse(value: int) -> int:
    value %= PRIME
    if value == 0:
        raise ZeroDivisionError("zero is not invertible")
    return pow(value, PRIME - 2, PRIME)


def normalize(values: Iterable[int]) -> tuple[int, ...]:
    vector = tuple(value % PRIME for value in values)
    leading = next((value for value in vector if value), None)
    if leading is None:
        raise ArithmeticError("zero has no projective normalization")
    scale = inverse(leading)
    return tuple(value * scale % PRIME for value in vector)


def rref(rows: Sequence[Sequence[int]]) -> tuple[list[list[int]], list[int]]:
    matrix = [[value % PRIME for value in row] for row in rows]
    if not matrix:
        return [], []
    row_index = 0
    pivots: list[int] = []
    for column in range(len(matrix[0])):
        pivot = next(
            (
                candidate
                for candidate in range(row_index, len(matrix))
                if matrix[candidate][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[row_index], matrix[pivot] = matrix[pivot], matrix[row_index]
        scale = inverse(matrix[row_index][column])
        matrix[row_index] = [scale * value % PRIME for value in matrix[row_index]]
        for other in range(len(matrix)):
            if other == row_index:
                continue
            scale = matrix[other][column]
            if scale:
                matrix[other] = [
                    (left - scale * right) % PRIME
                    for left, right in zip(matrix[other], matrix[row_index])
                ]
        pivots.append(column)
        row_index += 1
        if row_index == len(matrix):
            break
    return matrix, pivots


def nullspace(rows: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    reduced, pivots = rref(rows)
    column_count = len(rows[0])
    free = [column for column in range(column_count) if column not in pivots]
    result: list[tuple[int, ...]] = []
    for free_column in free:
        vector = [0] * column_count
        vector[free_column] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row_index][free_column] % PRIME
        result.append(tuple(vector))
    return result


def monomials(degree: int) -> list[tuple[int, int, int]]:
    return [
        (x_degree, y_degree, degree - x_degree - y_degree)
        for x_degree in range(degree, -1, -1)
        for y_degree in range(degree - x_degree, -1, -1)
    ]


def monomial_value(exponent: tuple[int, int, int], point: Point) -> int:
    value = 1
    for coordinate, power in zip(point, exponent):
        value = value * pow(coordinate, power, PRIME) % PRIME
    return value


def evaluation_row(degree: int, point: Point) -> list[int]:
    return [monomial_value(exponent, point) for exponent in monomials(degree)]


def derivative_row(degree: int, variable: int, point: Point) -> list[int]:
    row: list[int] = []
    for exponent in monomials(degree):
        power = exponent[variable]
        if power == 0:
            row.append(0)
            continue
        remaining = list(exponent)
        remaining[variable] -= 1
        row.append(power * monomial_value(tuple(remaining), point) % PRIME)
    return row


def derivative_value(
    coefficients: Sequence[int],
    degree: int,
    orders: tuple[int, int, int],
    point: Point,
) -> int:
    total = 0
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        multiplier = coefficient % PRIME
        remaining: list[int] = []
        for power, order in zip(exponent, orders):
            if order > power:
                multiplier = 0
                break
            for offset in range(order):
                multiplier = multiplier * (power - offset) % PRIME
            remaining.append(power - order)
        if multiplier:
            total += multiplier * monomial_value(tuple(remaining), point)
    return total % PRIME


def evaluate(coefficients: Sequence[int], degree: int, point: Point) -> int:
    return sum(
        coefficient * monomial_value(exponent, point)
        for coefficient, exponent in zip(coefficients, monomials(degree))
    ) % PRIME


def polynomial_string(coefficients: Sequence[int], degree: int) -> str:
    terms: list[str] = []
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        coefficient %= PRIME
        if coefficient == 0:
            continue
        factors = []
        for variable, power in zip(("x", "y", "z"), exponent):
            if power == 1:
                factors.append(variable)
            elif power > 1:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors) or "1"
        terms.append(monomial if coefficient == 1 else f"{coefficient}*{monomial}")
    return "+".join(terms)


def cross(left: Point, right: Point) -> tuple[int, int, int]:
    return normalize(
        (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )
    )  # type: ignore[return-value]


def line_value(line: Sequence[int], point: Point) -> int:
    return sum(left * right for left, right in zip(line, point)) % PRIME


def conic_determinant(coefficients: Sequence[int]) -> int:
    xx, xy, xz, yy, yz, zz = coefficients
    half = inverse(2)
    matrix = sp.Matrix(
        [
            [xx, half * xy, half * xz],
            [half * xy, yy, half * yz],
            [half * xz, half * yz, zz],
        ]
    )
    return int(matrix.det()) % PRIME


def perfect_matchings() -> list[tuple[Pair, Pair, Pair]]:
    edges = set(itertools.combinations(range(6), 2))
    result = []
    for triple in itertools.combinations(sorted(edges), 3):
        if sorted(vertex for edge in triple for vertex in edge) == list(range(6)):
            result.append(tuple(sorted(triple)))
    return sorted(result)


def one_factorizations() -> list[tuple[tuple[Pair, Pair, Pair], ...]]:
    edges = set(itertools.combinations(range(6), 2))
    matchings = perfect_matchings()
    result: set[tuple[tuple[Pair, Pair, Pair], ...]] = set()

    def search(
        remaining: set[Pair], blocks: list[tuple[Pair, Pair, Pair]]
    ) -> None:
        if not remaining:
            result.add(tuple(sorted(blocks)))
            return
        first = min(remaining)
        for matching in matchings:
            matching_set = set(matching)
            if first in matching_set and matching_set <= remaining:
                search(remaining - matching_set, [*blocks, matching])

    search(edges, [])
    return sorted(result)


def all_projective_points() -> list[Point]:
    return sorted(
        {
            normalize(vector)
            for vector in itertools.product(range(PRIME), repeat=3)
            if any(vector)
        }
    )  # type: ignore[return-value]


def coordinate_configuration() -> dict[str, object]:
    ordinary: list[Point] = [
        normalize((parameter * parameter, parameter, 1))  # type: ignore[arg-type]
        for parameter in range(PRIME)
    ]
    ordinary.append((1, 0, 0))
    if len(set(ordinary)) != 6 or any(
        (point[0] * point[2] - point[1] * point[1]) % PRIME
        for point in ordinary
    ):
        raise ArithmeticError("the six ordinary points are not the F5 conic")

    chord_lines = {
        edge: cross(ordinary[edge[0]], ordinary[edge[1]])
        for edge in itertools.combinations(range(6), 2)
    }
    matching_poles: dict[tuple[Pair, Pair, Pair], Point] = {}
    for matching in perfect_matchings():
        pole = cross(chord_lines[matching[0]], chord_lines[matching[1]])
        if line_value(chord_lines[matching[2]], pole) == 0:
            matching_poles[matching] = pole
    if len(matching_poles) != 10 or len(set(matching_poles.values())) != 10:
        raise ArithmeticError("the conic does not have the ten expected internal poles")

    factorization_records = []
    usable = []
    concurrency_histogram: Counter[int] = Counter()
    for factorization in one_factorizations():
        concurrent = tuple(
            matching for matching in factorization if matching in matching_poles
        )
        residual = tuple(
            matching for matching in factorization if matching not in matching_poles
        )
        concurrency_histogram[len(concurrent)] += 1
        record = {
            "matchings": [
                [list(edge) for edge in matching] for matching in factorization
            ],
            "concurrent_matching_count": len(concurrent),
            "concurrent_poles": [list(matching_poles[m]) for m in concurrent],
            "residual_matchings": [
                [list(edge) for edge in matching] for matching in residual
            ],
        }
        factorization_records.append(record)
        if len(concurrent) == 4 and len(residual) == 1:
            usable.append((concurrent, residual[0]))
    if (
        len(one_factorizations()) != 6
        or concurrency_histogram != Counter({4: 5, 0: 1})
        or len(usable) != 5
    ):
        raise ArithmeticError("the five tetrahedral factorization frames changed")

    selected_matchings = (
        ((0, 1), (2, 3), (4, 5)),
        ((0, 2), (1, 4), (3, 5)),
        ((0, 3), (1, 5), (2, 4)),
        ((0, 5), (1, 2), (3, 4)),
    )
    residual_matching = ((0, 4), (1, 3), (2, 5))
    if not all(matching in matching_poles for matching in selected_matchings):
        raise ArithmeticError("the selected four matchings are not concurrent")
    if tuple(sorted((*selected_matchings, residual_matching))) not in one_factorizations():
        raise ArithmeticError("the selected five matchings do not factor K6")
    special = [matching_poles[matching] for matching in selected_matchings]
    expected_special = [(1, 1, 4), (1, 3, 1), (1, 2, 2), (1, 0, 2)]
    if special != expected_special:
        raise ArithmeticError("the canonical special coordinates changed")

    centers: dict[PointLabel, Point] = {
        **{("s", index): point for index, point in enumerate(special)},
        **{("o", index): point for index, point in enumerate(ordinary)},
    }
    if len(set(centers.values())) != 10:
        raise ArithmeticError("the special and ordinary centers collided")

    lines = []
    for special_index, matching in enumerate(selected_matchings):
        for pair_index, pair in enumerate(matching):
            coefficients = chord_lines[pair]
            support = {
                label
                for label, point in centers.items()
                if line_value(coefficients, point) == 0
            }
            expected_support = {
                ("s", special_index),
                ("o", pair[0]),
                ("o", pair[1]),
            }
            if support != expected_support:
                raise ArithmeticError("a chord line has accidental center incidence")
            lines.append(
                {
                    "node": ("line", special_index, pair_index),
                    "special_index": special_index,
                    "ordinary_pair": pair,
                    "degree": 1,
                    "coefficients": coefficients,
                    "support_multiplicities": {
                        label: 1 for label in sorted(expected_support)
                    },
                }
            )
    if len({tuple(record["coefficients"]) for record in lines}) != 12:
        raise ArithmeticError("the twelve chord lines are not distinct")

    conic_supports: list[tuple[tuple[str, object], dict[PointLabel, int]]] = [
        (("ordinary",), {("o", index): 1 for index in range(6)})
    ]
    for pair_index, pair in enumerate(residual_matching):
        conic_supports.append(
            (
                ("residual", pair_index),
                {
                    **{("s", index): 1 for index in range(4)},
                    ("o", pair[0]): 1,
                    ("o", pair[1]): 1,
                },
            )
        )
    conics = []
    for label, support in conic_supports:
        matrix = [evaluation_row(2, centers[point]) for point in support]
        kernel = nullspace(matrix)
        if len(rref(matrix)[1]) != 5 or len(kernel) != 1:
            raise ArithmeticError("six centers do not determine a unique conic")
        coefficients = normalize(kernel[0])
        exact_support = {
            point
            for point, coordinates in centers.items()
            if evaluate(coefficients, 2, coordinates) == 0
        }
        if exact_support != set(support):
            raise ArithmeticError("a conic has accidental center incidence")
        determinant = conic_determinant(coefficients)
        if determinant == 0:
            raise ArithmeticError("a prescribed conic is singular")
        conics.append(
            {
                "node": ("conic", *label),
                "label": label,
                "degree": 2,
                "coefficients": coefficients,
                "support_multiplicities": support,
                "symmetric_matrix_determinant_mod_5": determinant,
            }
        )

    quartics = []
    for omitted_special in range(4):
        conditions = [evaluation_row(4, centers[("o", index)]) for index in range(6)]
        for special_index in range(4):
            if special_index == omitted_special:
                continue
            conditions.extend(
                derivative_row(4, variable, centers[("s", special_index)])
                for variable in range(3)
            )
        kernel = nullspace(conditions)
        if len(rref(conditions)[1]) != 14 or len(kernel) != 1:
            raise ArithmeticError("a prescribed quartic system is not a unique vector")
        coefficients = normalize(kernel[0])
        support = {
            **{("o", index): 1 for index in range(6)},
            **{
                ("s", index): 2
                for index in range(4)
                if index != omitted_special
            },
        }
        center_zero_set = {
            label
            for label, point in centers.items()
            if evaluate(coefficients, 4, point) == 0
        }
        if center_zero_set != set(support):
            raise ArithmeticError("quartic center support is not exact")
        quartics.append(
            {
                "node": ("quartic", omitted_special),
                "omitted_special": omitted_special,
                "degree": 4,
                "condition_row_count": len(conditions),
                "condition_rank": len(rref(conditions)[1]),
                "coefficients": coefficients,
                "support_multiplicities": support,
            }
        )
    if len({tuple(record["coefficients"]) for record in quartics}) != 4:
        raise ArithmeticError("the four quartics are not distinct")

    cubic_matrix = [evaluation_row(3, point) for point in centers.values()]
    cubic_rank = len(rref(cubic_matrix)[1])
    if cubic_rank != 10 or nullspace(cubic_matrix):
        raise ArithmeticError("an anticanonical cubic survives")

    sextic_conditions = []
    for point in centers.values():
        sextic_conditions.extend(
            derivative_row(6, variable, point) for variable in range(3)
        )
    sextic_kernel = nullspace(sextic_conditions)
    sextic_rank = len(rref(sextic_conditions)[1])
    if sextic_rank != 27 or len(sextic_kernel) != 1:
        raise ArithmeticError("the double-center sextic is not unique")
    boundary = normalize(sextic_kernel[0])
    value_augmented = [
        *sextic_conditions,
        *[evaluation_row(6, point) for point in centers.values()],
    ]
    if len(rref(value_augmented)[1]) != 27:
        raise ArithmeticError("Euler value consequences failed for the sextic")

    curves = [*lines, *conics, *quartics]
    return {
        "ordinary_points": ordinary,
        "special_points": special,
        "centers": centers,
        "all_projective_point_count": len(all_projective_points()),
        "perfect_matching_count": len(perfect_matchings()),
        "concurrent_matching_count": len(matching_poles),
        "one_factorization_count": len(one_factorizations()),
        "one_factorization_concurrency_histogram": concurrency_histogram,
        "five_coordinate_frames": factorization_records,
        "selected_matchings": selected_matchings,
        "residual_matching": residual_matching,
        "lines": lines,
        "conics": conics,
        "quartics": quartics,
        "curves": curves,
        "anticanonical": {
            "ambient_vector_dimension": 10,
            "evaluation_rank": cubic_rank,
            "kernel_dimension": 0,
        },
        "antibicanonical": {
            "ambient_vector_dimension": 28,
            "gradient_condition_row_count": len(sextic_conditions),
            "gradient_condition_rank": sextic_rank,
            "kernel_dimension": len(sextic_kernel),
            "value_augmented_rank": len(rref(value_augmented)[1]),
            "coefficients": boundary,
            "polynomial": polynomial_string(boundary, 6),
        },
    }


def quadratic_tangent_cone(
    coefficients: Sequence[int], degree: int, point: Point
) -> tuple[int, tuple[int, int], Point, tuple[int, int, int]]:
    chart = next(index for index, coordinate in enumerate(point) if coordinate)
    normalized_point = tuple(
        coordinate * inverse(point[chart]) % PRIME for coordinate in point
    )
    free = tuple(index for index in range(3) if index != chart)
    first2 = tuple(2 if index == free[0] else 0 for index in range(3))
    mixed = tuple(1 if index in free else 0 for index in range(3))
    second2 = tuple(2 if index == free[1] else 0 for index in range(3))
    half = inverse(2)
    quadratic = (
        half * derivative_value(coefficients, degree, first2, normalized_point)
        % PRIME,
        derivative_value(coefficients, degree, mixed, normalized_point),
        half * derivative_value(coefficients, degree, second2, normalized_point)
        % PRIME,
    )
    return chart, free, normalized_point, quadratic


def binary_quadratic_value(
    quadratic: tuple[int, int, int], direction: tuple[int, int]
) -> int:
    first, mixed, second = quadratic
    left, right = direction
    return (
        first * left * left + mixed * left * right + second * right * right
    ) % PRIME


def binary_quadratic_discriminant(quadratic: tuple[int, int, int]) -> int:
    first, mixed, second = quadratic
    return (mixed * mixed - 4 * first * second) % PRIME


def binary_quadratic_resultant(
    left: tuple[int, int, int], right: tuple[int, int, int]
) -> int:
    first, mixed, second = left
    other_first, other_mixed, other_second = right
    sylvester = sp.Matrix(
        [
            [first, mixed, second, 0],
            [0, first, mixed, second],
            [other_first, other_mixed, other_second, 0],
            [0, other_first, other_mixed, other_second],
        ]
    )
    return int(sylvester.det()) % PRIME


def univariate_restriction(
    coefficients: Sequence[int],
    degree: int,
    normalized_point: Point,
    chart: int,
    free: tuple[int, int],
    direction: tuple[int, int],
) -> list[int]:
    """Coefficients of F(p+t*d) in the affine chart, modulo five."""

    total = [0] * (degree + 1)
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        expansion = [coefficient % PRIME]
        for variable in range(3):
            if variable == chart:
                factor = [1]
            else:
                free_index = free.index(variable)
                center = normalized_point[variable]
                delta = direction[free_index]
                factor = [
                    math.comb(exponent[variable], power)
                    * pow(center, exponent[variable] - power, PRIME)
                    * pow(delta, power, PRIME)
                    % PRIME
                    for power in range(exponent[variable] + 1)
                ]
            product = [0] * (len(expansion) + len(factor) - 1)
            for left_index, left in enumerate(expansion):
                for right_index, right in enumerate(factor):
                    product[left_index + right_index] = (
                        product[left_index + right_index] + left * right
                    ) % PRIME
            expansion = product
        for power, value in enumerate(expansion):
            total[power] = (total[power] + value) % PRIME
    return total


def point_label_string(label: PointLabel) -> str:
    return f"{label[0]}{label[1]}"


def local_geometry(configuration: dict[str, object]) -> dict[str, object]:
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    boundary: Sequence[int] = configuration["antibicanonical"]["coefficients"]  # type: ignore[index]
    boundary_records = []
    boundary_quadratics: dict[PointLabel, tuple[int, int, int]] = {}
    node_discriminants: list[int] = []
    cusp_cubic_coefficients: list[int] = []
    directions = [(1, value) for value in range(PRIME)] + [(0, 1)]
    for label, point in centers.items():
        chart, free, normalized_point, quadratic = quadratic_tangent_cone(
            boundary, 6, point
        )
        boundary_quadratics[label] = quadratic
        discriminant = binary_quadratic_discriminant(quadratic)
        tangent_directions = [
            direction
            for direction in directions
            if binary_quadratic_value(quadratic, direction) == 0
        ]
        if label[0] == "s":
            if discriminant == 0:
                raise ArithmeticError("a special boundary singularity is not a node")
            node_discriminants.append(discriminant)
            singularity = "A1 ordinary node"
            cubic_tangent_coefficient = None
        else:
            if discriminant != 0 or len(tangent_directions) != 1:
                raise ArithmeticError("an ordinary-center boundary point is not unibranch tangent")
            restriction = univariate_restriction(
                boundary,
                6,
                normalized_point,
                chart,
                free,
                tangent_directions[0],
            )
            if restriction[0] or restriction[1] or restriction[2] or not restriction[3]:
                raise ArithmeticError("the unibranch double point is not an A2 cusp")
            cubic_tangent_coefficient = restriction[3]
            cusp_cubic_coefficients.append(restriction[3])
            singularity = "A2 ordinary cusp"
        boundary_records.append(
            {
                "point": point_label_string(label),
                "coordinates": list(point),
                "chart": chart,
                "quadratic_tangent_cone": list(quadratic),
                "discriminant_mod_5": discriminant,
                "F5_tangent_directions": [list(direction) for direction in tangent_directions],
                "cubic_coefficient_along_unique_tangent": cubic_tangent_coefficient,
                "singularity": singularity,
            }
        )

    quartic_node_records = []
    quartic_discriminants: list[int] = []
    quartic_boundary_resultants: list[int] = []
    smooth_root_values: list[tuple[str, int]] = []
    smooth_root_count = 0
    for curve in configuration["curves"]:  # type: ignore[union-attr]
        kind = curve["node"][0]
        degree = curve["degree"]
        coefficients = curve["coefficients"]
        support: dict[PointLabel, int] = curve["support_multiplicities"]
        for label, multiplicity in support.items():
            point = centers[label]
            if kind == "quartic" and multiplicity == 2:
                _, _, _, quartic_quadratic = quadratic_tangent_cone(
                    coefficients, degree, point
                )
                discriminant = binary_quadratic_discriminant(quartic_quadratic)
                if discriminant == 0:
                    raise ArithmeticError("a quartic prescribed double point is not a node")
                resultant = binary_quadratic_resultant(
                    boundary_quadratics[label], quartic_quadratic
                )
                if resultant == 0:
                    raise ArithmeticError("quartic and boundary node tangent cones share a branch")
                quartic_discriminants.append(discriminant)
                quartic_boundary_resultants.append(resultant)
                quartic_node_records.append(
                    {
                        "quartic_omits_special": curve["omitted_special"],
                        "point": point_label_string(label),
                        "quartic_tangent_cone": list(quartic_quadratic),
                        "quartic_discriminant_mod_5": discriminant,
                        "boundary_quartic_tangent_cone_resultant_mod_5": resultant,
                    }
                )
                continue

            chart, free, normalized_point, boundary_quadratic = quadratic_tangent_cone(
                boundary, 6, point
            )
            gradient = []
            for variable in free:
                orders = tuple(1 if index == variable else 0 for index in range(3))
                gradient.append(
                    derivative_value(
                        coefficients, degree, orders, normalized_point
                    )
                )
            if gradient == [0, 0]:
                raise ArithmeticError("a root curve is unexpectedly singular at a simple center")
            direction = (gradient[1], -gradient[0] % PRIME)
            tangent_value = binary_quadratic_value(boundary_quadratic, direction)
            if tangent_value == 0:
                raise ArithmeticError("a smooth root tangent equals a boundary tangent")
            smooth_root_count += 1
            smooth_root_values.append((kind, tangent_value))

    if smooth_root_count != 84 or len(quartic_node_records) != 12:
        raise ArithmeticError("the root/boundary tangent census is incomplete")
    return {
        "boundary_singularities": boundary_records,
        "boundary_node_count": 4,
        "boundary_cusp_count": 6,
        "boundary_node_discriminant_histogram": {
            str(value): count
            for value, count in sorted(Counter(node_discriminants).items())
        },
        "boundary_cusp_tangent_cubic_coefficient_histogram": {
            str(value): count
            for value, count in sorted(Counter(cusp_cubic_coefficients).items())
        },
        "boundary_delta_sum": 10,
        "boundary_resolved_by_one_blowup_at_each_center": True,
        "quartic_node_records": quartic_node_records,
        "quartic_node_count": len(quartic_node_records),
        "quartic_node_discriminant_histogram": {
            str(value): count
            for value, count in sorted(Counter(quartic_discriminants).items())
        },
        "smooth_root_center_incidence_count": smooth_root_count,
        "smooth_root_boundary_tangent_value_histogram": {
            str(value): count
            for value, count in sorted(
                Counter(value for _, value in smooth_root_values).items()
            )
        },
        "smooth_root_boundary_tangent_value_histogram_by_kind": {
            kind: {
                str(value): count
                for value, count in sorted(
                    Counter(
                        value
                        for current_kind, value in smooth_root_values
                        if current_kind == kind
                    ).items()
                )
            }
            for kind in ("line", "conic", "quartic")
        },
        "quartic_node_boundary_node_resultant_count": len(
            quartic_boundary_resultants
        ),
        "quartic_node_boundary_node_resultant_histogram": {
            str(value): count
            for value, count in sorted(
                Counter(quartic_boundary_resultants).items()
            )
        },
        "all_root_branches_avoid_all_boundary_branches": True,
        "bezout_exhaustion": {
            "line": "3 centers * local 2 = 6 = 1*6",
            "conic": "6 centers * local 2 = 12 = 2*6",
            "quartic": "6 simple centers *2 + 3 double centers *4 = 24 = 4*6",
        },
    }


def point_ideal_declaration(name: str, point: Point) -> str:
    variables = ("x", "y", "z")
    chart = next(index for index, coordinate in enumerate(point) if coordinate)
    scale = inverse(point[chart])
    normalized = tuple(value * scale % PRIME for value in point)
    relations = [
        f"{variables[index]}-{normalized[index]}*{variables[chart]}"
        for index in range(3)
        if index != chart
    ]
    return f"ideal {name}={','.join(relations)};"


def singular_factor_audit(configuration: dict[str, object]) -> dict[str, object]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the exact curve audit")
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    declarations = []
    point_names: dict[PointLabel, str] = {}
    for label, point in centers.items():
        name = point_label_string(label).upper()
        point_names[label] = name
        declarations.append(point_ideal_declaration(name, point))

    quartics = configuration["quartics"]  # type: ignore[assignment]
    boundary = configuration["antibicanonical"]  # type: ignore[assignment]
    polynomial_declarations = [
        f"poly q{index}={polynomial_string(record['coefficients'], 4)};"
        for index, record in enumerate(quartics)
    ]
    polynomial_declarations.append(f"poly b={boundary['polynomial']};")

    blocks = []
    for index, record in enumerate(quartics):
        singular_labels = [
            ("s", special)
            for special in range(4)
            if special != record["omitted_special"]
        ]
        point_intersection = ",".join(point_names[label] for label in singular_labels)
        prefix = f"Q{index}"
        blocks.append(
            f'''
list {prefix}F=factorize(q{index});
ideal {prefix}J=std(sat(jacob(q{index}),maxideal(1)));
ideal {prefix}R=std(radical({prefix}J));
ideal {prefix}T=std(intersect({point_intersection}));
print("{prefix}.FACTOR_COUNT_WITH_UNIT"); print(size({prefix}F[1]));
print("{prefix}.NONUNIT_EXPONENT"); print({prefix}F[2][2]);
print("{prefix}.AFFINE_CONE_DIMENSION"); print(dim({prefix}J));
print("{prefix}.SINGULAR_DEGREE"); print(mult({prefix}J));
print("{prefix}.RADICAL_DEGREE"); print(mult({prefix}R));
print("{prefix}.J_MOD_R"); print(size(reduce({prefix}J,{prefix}R)));
print("{prefix}.R_MOD_J"); print(size(reduce({prefix}R,{prefix}J)));
print("{prefix}.J_MOD_POINTS"); print(size(reduce({prefix}J,{prefix}T)));
print("{prefix}.POINTS_MOD_J"); print(size(reduce({prefix}T,{prefix}J)));
print("{prefix}.MINIMAL_PRIMES"); print(size(minAssGTZ({prefix}J)));
'''
        )

    all_points = ",".join(point_names[label] for label in centers)
    blocks.append(
        f'''
list BF=factorize(b);
ideal BJ=std(sat(jacob(b),maxideal(1)));
ideal BR=std(radical(BJ));
ideal BT=std(intersect({all_points}));
print("B.FACTOR_COUNT_WITH_UNIT"); print(size(BF[1]));
print("B.NONUNIT_EXPONENT"); print(BF[2][2]);
print("B.AFFINE_CONE_DIMENSION"); print(dim(BJ));
print("B.SINGULAR_DEGREE"); print(mult(BJ));
print("B.RADICAL_DEGREE"); print(mult(BR));
print("B.J_MOD_R"); print(size(reduce(BJ,BR)));
print("B.R_MOD_J"); print(size(reduce(BR,BJ)));
print("B.R_MOD_POINTS"); print(size(reduce(BR,BT)));
print("B.POINTS_MOD_R"); print(size(reduce(BT,BR)));
print("B.MINIMAL_PRIMES"); print(size(minAssGTZ(BJ)));
'''
    )
    program = "\n".join(
        [
            'LIB "primdec.lib";',
            "ring r=5,(x,y,z),dp;",
            *polynomial_declarations,
            *declarations,
            *blocks,
            "exit;",
        ]
    )
    completed = subprocess.run(
        [executable, "-q", "-c", program],
        check=True,
        capture_output=True,
        text=True,
    )
    labels = {}
    for prefix in [*(f"Q{index}" for index in range(4)), "B"]:
        for field in (
            "FACTOR_COUNT_WITH_UNIT",
            "NONUNIT_EXPONENT",
            "AFFINE_CONE_DIMENSION",
            "SINGULAR_DEGREE",
            "RADICAL_DEGREE",
            "J_MOD_R",
            "R_MOD_J",
            "MINIMAL_PRIMES",
        ):
            labels[f"{prefix}.{field}"] = (prefix, field.lower())
        if prefix == "B":
            labels[f"{prefix}.R_MOD_POINTS"] = (prefix, "radical_mod_points")
            labels[f"{prefix}.POINTS_MOD_R"] = (prefix, "points_mod_radical")
        else:
            labels[f"{prefix}.J_MOD_POINTS"] = (prefix, "jacobian_mod_points")
            labels[f"{prefix}.POINTS_MOD_J"] = (prefix, "points_mod_jacobian")
    parsed: dict[str, dict[str, int]] = {
        prefix: {} for prefix in [*(f"Q{index}" for index in range(4)), "B"]
    }
    output_lines = [
        line.strip() for line in completed.stdout.splitlines() if line.strip()
    ]
    for index, line in enumerate(output_lines[:-1]):
        if line in labels:
            prefix, field = labels[line]
            parsed[prefix][field] = int(output_lines[index + 1])

    quartic_expected = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "singular_degree": 3,
        "radical_degree": 3,
        "j_mod_r": 0,
        "r_mod_j": 0,
        "jacobian_mod_points": 0,
        "points_mod_jacobian": 0,
        "minimal_primes": 3,
    }
    for index in range(4):
        if parsed[f"Q{index}"] != quartic_expected:
            raise ArithmeticError(
                f"quartic {index} Singular certificate changed: {parsed[f'Q{index}']}"
            )
    boundary_expected_fixed = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "singular_degree": 16,
        "radical_degree": 10,
        "j_mod_r": 0,
        "radical_mod_points": 0,
        "points_mod_radical": 0,
        "minimal_primes": 10,
    }
    if any(parsed["B"].get(key) != value for key, value in boundary_expected_fixed.items()):
        raise ArithmeticError(f"boundary Singular certificate changed: {parsed['B']}")
    if parsed["B"].get("r_mod_j", 0) == 0:
        raise ArithmeticError("the cusp Jacobian scheme was falsely reduced")

    return {
        "engine": "Singular exact factorization, saturated Jacobians, radicals, point-ideal equalities, and minimal primes",
        "quartics": [
            {
                **parsed[f"Q{index}"],
                "irreducible_over_F5": True,
                "singular_scheme_is_exactly_three_reduced_special_points": True,
                "geometrically_integral": True,
            }
            for index in range(4)
        ],
        "quartic_geometric_integrality_argument": (
            "F5-irreducibility makes geometric components a Frobenius orbit. "
            "Ordinary nodes force at most two components. Two conjugate conics "
            "would meet transversely at the three rational nodes, so Bezout forces "
            "a fourth intersection and hence another singular point, contrary to "
            "the exact three-point singular scheme."
        ),
        "boundary": {
            **parsed["B"],
            "irreducible_over_F5": True,
            "geometrically_integral": True,
            "radical_is_exactly_the_ten_center_ideal": True,
            "jacobian_scheme_is_nonreduced": True,
            "jacobian_length_decomposition": "4 A1 nodes *1 + 6 A2 cusps *2 = 16",
        },
        "boundary_geometric_integrality_argument": (
            "If the F5-irreducible sextic split into r>=2 equal-degree Frobenius "
            "components, every F5-rational singular center would lie on all of them. "
            "Two components would share ten points, contradicting Bezout because "
            "their common degree 6/r is at most three."
        ),
        "hostile_warning": (
            "The boundary does not have ten nodes. Its saturated Jacobian scheme "
            "has length 16 and is nonreduced at the six A2 cusps; only its radical "
            "is the reduced ten-center scheme."
        ),
    }


def integer_vector(vector: sp.Matrix) -> list[int]:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("a marked divisor coordinate is not integral")
    return [int(value) for value in vector]


def primitive_integer(values: Sequence[Fraction | int | sp.Rational]) -> tuple[int, ...]:
    fractions = [Fraction(value) for value in values]
    denominator = math.lcm(*(value.denominator for value in fractions))
    integers = [
        value.numerator * (denominator // value.denominator) for value in fractions
    ]
    divisor = reduce(math.gcd, (abs(value) for value in integers), 0)
    if divisor == 0:
        raise ArithmeticError("cannot primitive-normalize zero")
    return tuple(value // divisor for value in integers)


def cdd_input(inequalities: sp.Matrix) -> str:
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(inequalities[row, column]) for column in range(inequalities.cols))
        for row in range(inequalities.rows)
    )
    lines.append("end")
    return "\n".join(lines) + "\n"


def root_cone_audit(
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    boundary_index: int,
    root_indices: Sequence[int],
    nonboundary_minus_one_indices: Sequence[int],
) -> dict[str, object]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the independent root-cone audit")
    boundary = walls[boundary_index]
    if (
        boundary != sp.Matrix([*([0] * 10), -1])
        or any(gram[index, 10] for index in range(10))
        or gram[10, 10] != -1
        or any(walls[index][10] for index in root_indices)
    ):
        raise ArithmeticError("the frozen CM complement is no longer the first ten coordinates")
    cm_gram = gram[:10, :10]
    cm_roots = [walls[index][:10, :] for index in root_indices]
    inequalities = sp.Matrix.vstack(
        *(root.T * cm_gram for root in cm_roots)
    )
    completed = subprocess.run(
        [executable, "--rep"],
        input=cdd_input(inequalities),
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    lines = completed.stdout.splitlines()
    try:
        representation = lines.index("V-representation")
        begin = next(
            index
            for index in range(representation + 1, len(lines))
            if lines[index].strip() == "begin"
        )
    except (ValueError, StopIteration) as error:
        raise ArithmeticError("cddlib omitted the root-cone V-representation") from error
    row_count, width, kind = lines[begin + 1].split()[:3]
    if kind != "rational" or int(width) != 11:
        raise ArithmeticError("cddlib returned the wrong root-cone generator width")
    rays = []
    for line in lines[begin + 2 : begin + 2 + int(row_count)]:
        values = [Fraction(token) for token in line.split()]
        if len(values) != 11 or values[0] != 0:
            raise ArithmeticError("the homogeneous root cone acquired a vertex")
        rays.append(primitive_integer(values[1:]))
    if lines[begin + 2 + int(row_count)].strip() != "end" or len(rays) != len(set(rays)):
        raise ArithmeticError("the cdd root ray block is malformed or duplicated")

    def square(ray: Sequence[int]) -> int:
        vector = sp.Matrix(ray)
        return int((vector.T * cm_gram * vector)[0])

    classifications = Counter(
        "finite" if square(ray) > 0 else "ideal" if square(ray) == 0 else "hyperideal"
        for ray in rays
    )
    if classifications != Counter({"finite": 450, "ideal": 51}):
        raise ArithmeticError(f"the direct root-cone census changed: {classifications}")
    ideal_rays = sorted(ray for ray in rays if square(ray) == 0)

    cusp_records = []
    raw_cusp_set = set()
    for wall_index in nonboundary_minus_one_indices:
        isotropic = walls[wall_index] + boundary
        if isotropic[10] != 0:
            raise ArithmeticError("E+b left the Coble-Mukai complement")
        ray = primitive_integer([int(value) for value in isotropic[:10, :]])
        if square(ray) != 0:
            raise ArithmeticError("E+b is not isotropic")
        raw_cusp_set.add(ray)
        cusp_records.append(
            {
                "raw_minus_one_wall_index_zero_based": wall_index,
                "primitive_f_equals_E_plus_b_in_CM_basis": list(ray),
            }
        )
    if len(raw_cusp_set) != 51 or raw_cusp_set != set(ideal_rays):
        raise ArithmeticError("the 51 raw E+b classes are not exactly the root-cone cusps")

    ray_digest = hashlib.sha256(
        json.dumps(sorted(rays), separators=(",", ":")).encode()
    ).hexdigest()
    cusp_digest = hashlib.sha256(
        json.dumps(ideal_rays, separators=(",", ":")).encode()
    ).hexdigest()
    interior = sum((sp.Matrix(ray) for ray in rays), sp.zeros(10, 1))
    interior_square = int((interior.T * cm_gram * interior)[0])
    slacks = [int(value) for value in inequalities * interior]
    time_pairings = [
        int((interior.T * cm_gram * sp.Matrix(ray))[0]) for ray in rays
    ]
    if interior_square <= 0 or min(slacks) <= 0 or min(time_pairings) <= 0:
        raise ArithmeticError("the root cone did not stay on one Lorentz time sheet")
    active_histogram = Counter(
        sum(
            int((inequalities * sp.Matrix(ray))[row]) == 0
            for row in range(inequalities.rows)
        )
        for ray in rays
    )
    return {
        "engine": "direct cddexec_gmp rational H-to-V conversion of the twenty raw root inequalities",
        "ray_count": len(rays),
        "finite_ray_count": classifications["finite"],
        "ideal_ray_count": classifications["ideal"],
        "hyperideal_ray_count": classifications["hyperideal"],
        "active_root_count_histogram": {
            str(value): count for value, count in sorted(active_histogram.items())
        },
        "primitive_ray_list_sha256": ray_digest,
        "primitive_ideal_ray_list_sha256": cusp_digest,
        "strict_interior_vector": integer_vector(interior),
        "strict_interior_square": interior_square,
        "minimum_root_inequality_slack": min(slacks),
        "minimum_time_sheet_pairing": min(time_pairings),
        "all_rays_on_one_closed_positive_cone_sheet": True,
        "finite_volume_root_chamber": True,
        "cusp_records": sorted(
            cusp_records,
            key=lambda record: record["raw_minus_one_wall_index_zero_based"],
        ),
        "all_51_ideal_rays_equal_the_51_primitive_E_plus_b_classes": True,
    }


def raw_graph_for_frame(
    walls: Sequence[sp.Matrix],
    inverse_basis: sp.Matrix,
    root_indices: Sequence[int],
) -> tuple[nx.Graph, list[dict[str, object]]]:
    graph = nx.Graph()
    for position in range(10):
        graph.add_node(("point", position), kind="point")
    root_records = []
    for root_index in root_indices:
        coordinates = integer_vector(inverse_basis * walls[root_index])
        degree = coordinates[0]
        kind = {1: "line", 2: "conic", 4: "quartic"}.get(degree)
        if kind is None:
            raise ArithmeticError("a raw root has a forbidden plane degree")
        root_node = (kind, root_index)
        graph.add_node(root_node, kind=kind)
        multiplicities = {}
        for position, coefficient in enumerate(coordinates[1:]):
            if coefficient == 0:
                continue
            multiplicity = -coefficient
            if multiplicity not in (1, 2):
                raise ArithmeticError("a raw root has a forbidden center multiplicity")
            multiplicities[position] = multiplicity
            incidence_node = ("incidence", kind, root_index, position)
            graph.add_node(incidence_node, kind=f"multiplicity_{multiplicity}")
            graph.add_edge(root_node, incidence_node)
            graph.add_edge(incidence_node, ("point", position))
        expected_profile = {
            "line": Counter({1: 3}),
            "conic": Counter({1: 6}),
            "quartic": Counter({1: 6, 2: 3}),
        }[kind]
        if Counter(multiplicities.values()) != expected_profile:
            raise ArithmeticError("a raw root has the wrong full multiplicity profile")
        root_records.append(
            {
                "raw_wall_index_zero_based": root_index,
                "kind": kind,
                "coordinates_in_H_E_basis": coordinates,
                "frame_multiplicities": multiplicities,
            }
        )
    return graph, root_records


def geometric_incidence_graph(
    configuration: dict[str, object],
) -> tuple[nx.Graph, dict[tuple[object, ...], dict[str, object]]]:
    graph = nx.Graph()
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    for label in centers:
        graph.add_node(("point", label), kind="point")
    curve_lookup = {}
    for curve in configuration["curves"]:  # type: ignore[union-attr]
        node = tuple(curve["node"])
        kind = node[0]
        graph.add_node(node, kind=kind)
        curve_lookup[node] = curve
        support: dict[PointLabel, int] = curve["support_multiplicities"]
        for incidence_index, (point, multiplicity) in enumerate(sorted(support.items())):
            incidence_node = ("incidence", *node, incidence_index)
            graph.add_node(
                incidence_node, kind=f"multiplicity_{multiplicity}"
            )
            graph.add_edge(node, incidence_node)
            graph.add_edge(incidence_node, ("point", point))
    return graph, curve_lookup


def frame_marking_audit(
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    boundary: sp.Matrix,
    root_indices: Sequence[int],
    frame_set: Sequence[int],
    configuration: dict[str, object],
) -> tuple[dict[str, object], list[list[int]]]:
    sorted_frame = sorted(frame_set)
    plane_class = (
        boundary + sum((walls[index] for index in sorted_frame), sp.zeros(11, 1))
    ) / 3
    oriented_frame = list(sorted_frame)
    trial_basis = sp.Matrix.hstack(
        plane_class, *[walls[index] for index in oriented_frame]
    )
    if trial_basis.det() == -1:
        oriented_frame[-1], oriented_frame[-2] = oriented_frame[-2], oriented_frame[-1]
    basis = sp.Matrix.hstack(
        plane_class, *[walls[index] for index in oriented_frame]
    )
    expected_gram = sp.diag(1, *([-1] * 10))
    if basis.det() != 1 or basis.T * gram * basis != expected_gram:
        raise ArithmeticError("a ten-frame does not yield a determinant-one I(1,10) basis")
    inverse_basis = basis.inv()
    boundary_coordinates = integer_vector(inverse_basis * boundary)
    if boundary_coordinates != [3, *([-1] * 10)]:
        raise ArithmeticError("the raw boundary is not 3H-sum(E) in a frame")

    raw_graph, root_records = raw_graph_for_frame(
        walls, inverse_basis, root_indices
    )
    geometric_graph, curve_lookup = geometric_incidence_graph(configuration)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        geometric_graph,
        node_match=lambda left, right: left["kind"] == right["kind"],
    )
    markings = []
    every_vector_exact = True
    for mapping in matcher.isomorphisms_iter():
        encoding = ",".join(
            point_label_string(mapping[("point", position)][1])
            for position in range(10)
        )
        for record in root_records:
            raw_node = (record["kind"], record["raw_wall_index_zero_based"])
            geometric_node = mapping[raw_node]
            curve = curve_lookup[geometric_node]
            expected_coordinates = [curve["degree"]]
            support: dict[PointLabel, int] = curve["support_multiplicities"]
            expected_coordinates.extend(
                -support.get(mapping[("point", position)][1], 0)
                for position in range(10)
            )
            if record["coordinates_in_H_E_basis"] != expected_coordinates:
                every_vector_exact = False
                break
        markings.append((encoding, mapping))
    if len(markings) != 24 or len({encoding for encoding, _ in markings}) != 24:
        raise ArithmeticError("a frame does not have the expected 24 colored markings")
    if not every_vector_exact:
        raise ArithmeticError("an incidence marking failed exact root-vector equality")
    canonical_encoding, canonical_mapping = min(markings, key=lambda item: item[0])

    canonical_points = [
        {
            "raw_frame_wall_index_zero_based": oriented_frame[position],
            "geometric_center": point_label_string(
                canonical_mapping[("point", position)][1]
            ),
        }
        for position in range(10)
    ]
    canonical_roots = []
    for record in root_records:
        raw_node = (record["kind"], record["raw_wall_index_zero_based"])
        geometric_node = canonical_mapping[raw_node]
        canonical_roots.append(
            {
                "raw_wall_index_zero_based": record["raw_wall_index_zero_based"],
                "kind": record["kind"],
                "geometric_curve_node": [
                    list(item) if isinstance(item, tuple) else item
                    for item in geometric_node
                ],
                "coordinates_in_H_E_basis": record["coordinates_in_H_E_basis"],
            }
        )

    full_wall_coordinates = [integer_vector(inverse_basis * wall) for wall in walls]
    full_digest = hashlib.sha256(
        json.dumps(full_wall_coordinates, separators=(",", ":")).encode()
    ).hexdigest()
    torsor_digest = hashlib.sha256(
        "\n".join(sorted(encoding for encoding, _ in markings)).encode()
    ).hexdigest()
    profile = Counter(record["kind"] for record in root_records)
    special_positions = [
        position
        for position in range(10)
        if sum(
            record["frame_multiplicities"].get(position) == 2
            for record in root_records
        )
        == 3
    ]
    ordinary_positions = [
        position for position in range(10) if position not in special_positions
    ]
    if len(special_positions) != 4 or len(ordinary_positions) != 6:
        raise ArithmeticError("a raw frame does not intrinsically split 4+6")
    return (
        {
            "frame_set_wall_indices_zero_based": sorted_frame,
            "determinant_one_ordered_frame_wall_indices_zero_based": oriented_frame,
            "plane_class_in_raw_basis": integer_vector(plane_class),
            "basis_determinant": int(basis.det()),
            "basis_gram": [
                [int(expected_gram[row, column]) for column in range(11)]
                for row in range(11)
            ],
            "boundary_coordinates_in_H_E_basis": boundary_coordinates,
            "intrinsic_special_frame_positions": special_positions,
            "intrinsic_ordinary_frame_positions": ordinary_positions,
            "root_kind_histogram": dict(sorted(profile.items())),
            "colored_incidence_marking_count": len(markings),
            "all_24_markings_give_exact_twenty_root_vectors": True,
            "marking_torsor_sha256": torsor_digest,
            "canonical_marking_encoding": canonical_encoding,
            "canonical_point_marking": canonical_points,
            "canonical_root_marking": canonical_roots,
            "all_72_wall_coordinates_integral": True,
            "transported_full_wall_coordinates_sha256": full_digest,
        },
        full_wall_coordinates,
    )


def raw_polytope_audit(
    raw: dict[str, object], configuration: dict[str, object]
) -> dict[str, object]:
    polytope = raw["polytope"]  # type: ignore[index]
    gram = sp.Matrix(polytope["lattice_gram"])  # type: ignore[index]
    walls = [sp.Matrix(row) for row in polytope["wall_vectors"]]  # type: ignore[index]

    def pairing(left: sp.Matrix, right: sp.Matrix) -> int:
        return int((left.T * gram * right)[0])

    norm_histogram = Counter(pairing(wall, wall) for wall in walls)
    root_indices = [
        index for index, wall in enumerate(walls) if pairing(wall, wall) == -2
    ]
    minus_one_indices = [
        index for index, wall in enumerate(walls) if pairing(wall, wall) == -1
    ]
    boundary_candidates = [
        index
        for index in minus_one_indices
        if all(pairing(walls[index], walls[root]) == 0 for root in root_indices)
    ]
    if len(walls) != 72 or norm_histogram != Counter({-1: 52, -2: 20}):
        raise ArithmeticError("the raw Type-VII wall inventory changed")
    if len(boundary_candidates) != 1:
        raise ArithmeticError("the raw boundary is not intrinsically unique")
    boundary_index = boundary_candidates[0]
    boundary = walls[boundary_index]
    nonboundary = [
        index for index in minus_one_indices if index != boundary_index
    ]
    if any(pairing(boundary, walls[index]) != 1 for index in nonboundary):
        raise ArithmeticError("a raw minus-one wall violates b.E=1")

    disjointness = nx.Graph()
    disjointness.add_nodes_from(nonboundary)
    disjointness.add_edges_from(
        (left, right)
        for left, right in itertools.combinations(nonboundary, 2)
        if pairing(walls[left], walls[right]) == 0
    )
    maximal_cliques = [sorted(clique) for clique in nx.find_cliques(disjointness)]
    clique_histogram = Counter(len(clique) for clique in maximal_cliques)
    frames = sorted(clique for clique in maximal_cliques if len(clique) == 10)
    if clique_histogram != Counter({7: 80, 9: 25, 10: 5}) or len(frames) != 5:
        raise ArithmeticError("the intrinsic five-frame census changed")

    frame_records = []
    full_coordinate_archives = []
    for frame in frames:
        record, coordinates = frame_marking_audit(
            gram, walls, boundary, root_indices, frame, configuration
        )
        frame_records.append(record)
        full_coordinate_archives.append(coordinates)
    aggregate_digest = hashlib.sha256(
        json.dumps(full_coordinate_archives, separators=(",", ":")).encode()
    ).hexdigest()

    root_cone = root_cone_audit(
        gram,
        walls,
        boundary_index,
        root_indices,
        nonboundary,
    )
    exact = polytope["exact_double_description"]  # type: ignore[index]
    vinberg = polytope["literal_vinberg"]  # type: ignore[index]
    checks = raw["checks"]  # type: ignore[index]
    if not (
        raw["all_checks_pass"]  # type: ignore[index]
        and all(checks.values())  # type: ignore[union-attr]
        and polytope["finite_volume"]  # type: ignore[index]
        and not polytope["compact"]  # type: ignore[index]
        and exact["hyperideal_vertex_count"] == 0  # type: ignore[index]
        and vinberg["bad_edge_count"] == 0  # type: ignore[index]
        and vinberg["bad_parabolic_neighbourhood_count"] == 0  # type: ignore[index]
    ):
        raise ArithmeticError("the frozen raw P-side Vinberg certificate is not green")
    return {
        "raw_wall_count": len(walls),
        "primitive_wall_norm_histogram": {
            str(value): count for value, count in sorted(norm_histogram.items())
        },
        "intrinsic_boundary_wall_index_zero_based": boundary_index,
        "boundary_characterization": "unique norm-minus-one wall orthogonal to all twenty norm-minus-two roots",
        "root_wall_indices_zero_based": root_indices,
        "nonboundary_minus_one_wall_count": len(nonboundary),
        "maximal_disjoint_clique_size_histogram": {
            str(value): count for value, count in sorted(clique_histogram.items())
        },
        "intrinsic_plane_frame_count": len(frames),
        "all_five_frames": frame_records,
        "all_five_frames_have_determinant_one_I_1_10_bases": True,
        "all_120_colored_markings_give_exact_root_vectors": True,
        "aggregate_all_five_full_wall_coordinate_sha256": aggregate_digest,
        "direct_twenty_root_cone": root_cone,
        "general_one_boundary_return": {
            "input": "the constructed strong Coble surface plus its exact twenty effective root marking",
            "finite_index": "the direct 450+51+0 root-cone census makes the twenty-root reflection subgroup finite index",
            "root_exhaustion": "the general Coble Vinberg finite-index proposition excludes hidden effective irreducible roots; no classification/existence theorem is used",
            "cusp_lift": "for every primitive nef cusp f, E=f-b=f+K is the one-boundary (-1)-curve; conversely every (-1)-curve gives f=E+b",
            "negative_curve_exhaustion": "boundary, the twenty (-2)-curves, and the 51 cusp lifts exhaust all negative curves",
            "same_full_wall_set": True,
        },
        "upstream_full_P_vinberg_return": {
            "polytope_id": polytope["id"],  # type: ignore[index]
            "dimension": polytope["dimension"],  # type: ignore[index]
            "facet_count": polytope["facet_count"],  # type: ignore[index]
            "finite_vertices": exact["finite_vertex_count"],  # type: ignore[index]
            "ideal_vertices": exact["ideal_vertex_count"],  # type: ignore[index]
            "hyperideal_extreme_rays": exact["hyperideal_vertex_count"],  # type: ignore[index]
            "elliptic_edges": vinberg["exact_geometric_edge_count"],  # type: ignore[index]
            "finite_volume": polytope["finite_volume"],  # type: ignore[index]
            "compact": polytope["compact"],  # type: ignore[index]
            "all_facets_irredundant": exact["all_facets_irredundant"],  # type: ignore[index]
            "literal_vinberg_zero_failures": True,
        },
        "exact_same_marked_H10_72": True,
        "import_boundary": (
            "Only the already-frozen raw 72-wall face/Vinberg enumeration and "
            "general Coble cone/cusp theorems are imported. The existence, "
            "coordinates, singularities, roots, root-cone finite volume, and "
            "marked lattice return are reconstructed here without importing a "
            "known Type-VII surface."
        ),
    }


def curve_payload(curve: dict[str, object]) -> dict[str, object]:
    support: dict[PointLabel, int] = curve["support_multiplicities"]  # type: ignore[assignment]
    degree = int(curve["degree"])
    payload: dict[str, object] = {
        "node": list(curve["node"]),  # type: ignore[arg-type]
        "degree": degree,
        "coefficients": list(curve["coefficients"]),  # type: ignore[arg-type]
        "polynomial": polynomial_string(curve["coefficients"], degree),  # type: ignore[arg-type]
        "support_multiplicities": [
            {"point": point_label_string(label), "multiplicity": multiplicity}
            for label, multiplicity in sorted(support.items())
        ],
    }
    for field in (
        "special_index",
        "ordinary_pair",
        "label",
        "symmetric_matrix_determinant_mod_5",
        "omitted_special",
        "condition_row_count",
        "condition_rank",
    ):
        if field in curve:
            value = curve[field]
            if isinstance(value, tuple):
                value = list(value)
            payload[field] = value
    return payload


def build_certificate() -> dict[str, object]:
    raw_design = json.loads(RAW_DESIGN.read_text(encoding="utf-8"))
    raw_polytope = json.loads(RAW_POLYTOPE.read_text(encoding="utf-8"))
    design_row = raw_design["seven_row_outcome"]["one_boundary_72"]
    design_candidate = raw_design["different_finite_design_candidates"][
        "one_boundary_72"
    ]
    if not (
        raw_design["all_checks_pass"]
        and design_row["wall_count"] == 72
        and design_row["simultaneous_full_P2_frame_count"] == 5
        and design_row["root_count"] == 20
        and design_candidate["plane_frame_count"] == 5
    ):
        raise ArithmeticError("the raw finite-design starting gate changed")

    configuration = coordinate_configuration()
    local = local_geometry(configuration)
    singular = singular_factor_audit(configuration)
    raw = raw_polytope_audit(raw_polytope, configuration)

    finite_geometry = {
        "field": "F5",
        "ordinary_conic": "x*z-y^2=0",
        "ordinary_points": [
            {"point": f"o{index}", "coordinates": list(point)}
            for index, point in enumerate(configuration["ordinary_points"])
        ],
        "special_points": [
            {"point": f"s{index}", "coordinates": list(point)}
            for index, point in enumerate(configuration["special_points"])
        ],
        "projective_plane_point_count": configuration["all_projective_point_count"],
        "K6_perfect_matching_count": configuration["perfect_matching_count"],
        "concurrent_chord_matching_count": configuration[
            "concurrent_matching_count"
        ],
        "K6_one_factorization_count": configuration["one_factorization_count"],
        "one_factorization_concurrency_histogram": {
            str(value): count
            for value, count in sorted(
                configuration["one_factorization_concurrency_histogram"].items()
            )
        },
        "usable_tetrahedral_factorization_count": 5,
        "all_six_factorizations": configuration["five_coordinate_frames"],
        "selected_four_concurrent_matchings": [
            [list(edge) for edge in matching]
            for matching in configuration["selected_matchings"]
        ],
        "selected_residual_matching": [
            list(edge) for edge in configuration["residual_matching"]
        ],
        "triangle_line_count": len(configuration["lines"]),
        "conic_count": len(configuration["conics"]),
        "quartic_count": len(configuration["quartics"]),
        "lines": [curve_payload(curve) for curve in configuration["lines"]],
        "conics": [curve_payload(curve) for curve in configuration["conics"]],
        "quartics": [curve_payload(curve) for curve in configuration["quartics"]],
        "all_four_conics_smooth_geometrically_integral": True,
        "all_curve_center_supports_exact": True,
    }
    linear_systems = {
        "anticanonical_cubics": configuration["anticanonical"],
        "antibicanonical_sextics": {
            **configuration["antibicanonical"],
            "coefficients": list(
                configuration["antibicanonical"]["coefficients"]
            ),
        },
    }
    strong_coble = {
        "surface": "the blowup over algebraic closure of F5 of the ten explicit centers in P2",
        "rational_surface": True,
        "canonical_class": "K=-3H+sum(E_i); K^2=-1",
        "boundary_class": "B=6H-2sum(E_i)=-2K; B^2=-4",
        "minus_K_empty": configuration["anticanonical"]["kernel_dimension"] == 0,
        "minus_2K_unique": configuration["antibicanonical"]["kernel_dimension"] == 1,
        "boundary_geometrically_integral": singular["boundary"][
            "geometrically_integral"
        ],
        "boundary_plane_singularity_profile": "4 A1 nodes plus 6 A2 cusps",
        "boundary_arithmetic_genus": 10,
        "boundary_delta_sum": local["boundary_delta_sum"],
        "boundary_strict_transform_smooth_rational": True,
        "one_blowup_resolves_each_node_and_each_cusp": True,
        "root_class_inventory": {
            "12 lines": "H-E_s-E_o-E_o; square -2; K-pairing 0",
            "4 conics": "2H-sum(six simple E); square -2; K-pairing 0",
            "4 quartics": "4H-sum(six ordinary E)-2sum(three special E); square -2; K-pairing 0",
        },
        "all_twenty_root_strict_transforms_smooth_geometrically_integral_rational": True,
        "all_twenty_roots_disjoint_from_boundary": local[
            "all_root_branches_avoid_all_boundary_branches"
        ],
        "strong_Coble_verdict": "PASS",
    }

    checks = {
        "raw_P_intrinsically_has_exactly_five_ten_curve_frames": raw[
            "intrinsic_plane_frame_count"
        ]
        == 5,
        "every_raw_frame_has_determinant_one_I_1_10_basis": raw[
            "all_five_frames_have_determinant_one_I_1_10_bases"
        ],
        "finite_geometry_has_five_usable_factorizations": finite_geometry[
            "usable_tetrahedral_factorization_count"
        ]
        == 5,
        "twelve_lines_four_conics_four_quartics_constructed": (
            finite_geometry["triangle_line_count"],
            finite_geometry["conic_count"],
            finite_geometry["quartic_count"],
        )
        == (12, 4, 4),
        "all_four_quartics_geometrically_integral_with_exact_three_nodes": all(
            record["geometrically_integral"]
            and record["singular_scheme_is_exactly_three_reduced_special_points"]
            for record in singular["quartics"]
        ),
        "no_anticanonical_cubic": strong_coble["minus_K_empty"],
        "unique_antibicanonical_sextic": strong_coble["minus_2K_unique"],
        "boundary_geometrically_integral": strong_coble[
            "boundary_geometrically_integral"
        ],
        "boundary_exactly_four_nodes_six_cusps": local["boundary_node_count"] == 4
        and local["boundary_cusp_count"] == 6,
        "boundary_radical_exactly_ten_centers_but_Jacobian_length_16": singular[
            "boundary"
        ]["radical_is_exactly_the_ten_center_ideal"]
        and singular["boundary"]["singular_degree"] == 16,
        "all_root_boundary_tangent_gates": local[
            "all_root_branches_avoid_all_boundary_branches"
        ],
        "strong_Coble": strong_coble["strong_Coble_verdict"] == "PASS",
        "all_120_markings_match_exact_root_vectors": raw[
            "all_120_colored_markings_give_exact_root_vectors"
        ],
        "all_five_integral_full_72_wall_transports": all(
            frame["all_72_wall_coordinates_integral"]
            for frame in raw["all_five_frames"]
        ),
        "direct_root_cone_450_51_0": (
            raw["direct_twenty_root_cone"]["finite_ray_count"],
            raw["direct_twenty_root_cone"]["ideal_ray_count"],
            raw["direct_twenty_root_cone"]["hyperideal_ray_count"],
        )
        == (450, 51, 0),
        "exact_51_cusp_lifts": raw["direct_twenty_root_cone"][
            "all_51_ideal_rays_equal_the_51_primitive_E_plus_b_classes"
        ],
        "same_marked_finite_volume_noncompact_H10_72": raw[
            "exact_same_marked_H10_72"
        ]
        and raw["upstream_full_P_vinberg_return"]["finite_volume"]
        and not raw["upstream_full_P_vinberg_return"]["compact"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"a Type-VII inverse gate failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "status": "PROVED_EXPLICIT_P_FIRST_STRONG_COBLE_AND_SAME_MARKED_H10_72",
        "direction": "RAW_POLYTOPE_TO_INTRINSIC_DESIGN_TO_EXPLICIT_SURFACE_TO_EXACT_SAME_POLYTOPE",
        "sources": {
            "raw_design_archive": {
                "path": str(RAW_DESIGN.relative_to(ROOT)),
                "sha256": sha256(RAW_DESIGN),
            },
            "raw_H10_72_archive": {
                "path": str(RAW_POLYTOPE.relative_to(ROOT)),
                "sha256": sha256(RAW_POLYTOPE),
                "surface_existence_or_classification_imported": False,
                "used_fields": "lattice Gram, 72 wall vectors, and independently frozen P-side Vinberg census",
            },
        },
        "raw_starting_design_gate": {
            "lattice_rank": design_row["lattice_rank"],
            "wall_count": design_row["wall_count"],
            "root_count": design_row["root_count"],
            "plane_frame_count": design_row["simultaneous_full_P2_frame_count"],
            "raw_profile": design_candidate["root_plane_profile"],
        },
        "finite_geometry": finite_geometry,
        "linear_systems": linear_systems,
        "local_geometry_and_transversality": local,
        "exact_factor_and_singular_scheme_audit": singular,
        "blowup_and_strong_Coble": strong_coble,
        "marked_same_polytope_return": raw,
        "inverse_level": {
            "rating": "L3 structural P-first calibration, not class-wide L4",
            "positive_content": (
                "The raw P intrinsically forces five frames and the 4+6 design; one "
                "is explicitly realized, all AG gates close, and the exact marked "
                "72-wall P returns without importing known surface existence."
            ),
            "why_not_L4": (
                "This proves one rigid datum. It does not classify a variable family "
                "of raw arithmetic Coxeter polytopes or give necessary-and-sufficient "
                "carrier criteria for that family."
            ),
        },
        "lineage": {
            "known_after_the_fact": (
                "The returned root graph and chamber are the known characteristic-five "
                "Type-VII / H10/72 datum. This label is not used for existence."
            ),
            "known_printed_models": "Kondo prints a P4 Cremona-quotient model, not this plane certificate",
            "project_new": (
                "the explicit F5 plane equation package and strict raw-P-to-surface-to-"
                "same-marked-P proof"
            ),
            "worldwide_priority": "open; no worldwide-new surface, equation, or abstract-polytope claim",
            "new_abstract_hyperbolic_polytope": False,
            "new_dimension": False,
            "new_general_dimension_bound": False,
        },
        "hostile_warning": singular["hostile_warning"],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archive mismatch: {arguments.output}")
        print("Type-VII characteristic-five tetrahedral inverse archive replay passed")
        return
    if arguments.write:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
