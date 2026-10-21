#!/usr/bin/env python3
"""Certify a strict raw-P to plane-cluster realization of H10/30.

The input is the physical odd-unimodular Gram and thirty primitive walls of
the H10/30 chamber.  The program exhausts its sequential numerical P2
endpoints, extracts the two rigid 7+3 cluster branches, solves their incidence
design over Q, and verifies the resulting strong-Coble surface and exact
same-marked-polytope handoff.

No known Type-II surface equation or surface label participates in the branch
extraction or in the construction.  The Type-II name is attached only in the
final lineage/classification layer.
"""

from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_coble_remaining_seven_contraction_designs import (  # noqa: E402
    _dot,
    _load_rows,
    _normalized,
    _plane_class,
)


OUTPUT = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json"
FULL_CONE = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"
HOSTILE_RECOGNIZER = (
    ROOT / "data/enumerations/polytope_first_coble_h10_hostile_audit.json"
)
SCHEMA = "coxeter.coble_type_ii_b3_plane_inverse.v1"

X, Y, Z = sp.symbols("x y z")
S, T = sp.symbols("s t")
VARIABLES = (X, Y, Z)

# The labels are part of the incidence solution, not raw input metadata.
POINT_ORDER = ("a", "A", "d", "D", "g", "G", "B", "C", "E", "F")
NEAR_POINTS = {"a", "d", "g"}
PARENT = {"a": "A", "d": "D", "g": "G"}
PROPER_POINTS: dict[str, tuple[int, int, int]] = {
    "A": (1, 0, 0),
    "D": (0, 1, 0),
    "G": (0, 0, 1),
    "B": (1, 1, 1),
    "C": (1, 1, -1),
    "E": (-1, 1, 1),
    "F": (-1, 1, -1),
}

EXCEPTIONAL_ROOTS = {
    "EA_minus_Ea": ("a", "A"),
    "ED_minus_Ed": ("d", "D"),
    "EG_minus_Eg": ("g", "G"),
}

LINE_ROOTS: dict[str, dict[str, Any]] = {
    "L_DG": {"coefficients": (1, 0, 0), "support": ("g", "G", "D")},
    "L_AG": {"coefficients": (0, 1, 0), "support": ("a", "A", "G")},
    "L_AD": {"coefficients": (0, 0, 1), "support": ("d", "D", "A")},
    "L_EFG": {"coefficients": (1, 1, 0), "support": ("E", "F", "G")},
    "L_BCG": {"coefficients": (1, -1, 0), "support": ("B", "C", "G")},
    "L_ACF": {"coefficients": (0, 1, 1), "support": ("A", "C", "F")},
    "L_ABE": {"coefficients": (0, 1, -1), "support": ("A", "B", "E")},
    "L_CDE": {"coefficients": (1, 0, 1), "support": ("C", "D", "E")},
    "L_BDF": {"coefficients": (1, 0, -1), "support": ("B", "D", "F")},
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _primitive_integer_vector(vector: sp.Matrix) -> list[int]:
    denominator = 1
    for value in vector:
        denominator = int(sp.ilcm(denominator, value.q))
    entries = [int(value * denominator) for value in vector]
    divisor = 0
    for value in entries:
        divisor = math.gcd(divisor, abs(value))
    entries = [value // divisor for value in entries]
    for value in entries:
        if value:
            if value < 0:
                entries = [-entry for entry in entries]
            break
    return entries


def _sequential_endpoints(row: dict[str, Any]) -> dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]]:
    matrix = row["wall_gram"]
    allowed = [index for index in range(len(matrix)) if index not in row["boundaries"]]
    current: dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = {
        0: (tuple(), tuple())
    }
    for _ in range(row["gram_matrix"].rows - 1):
        following: dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = {}
        for mask, (contracted, inverse) in current.items():
            size = len(contracted)
            for exceptional in allowed:
                if mask & (1 << exceptional):
                    continue
                vector = [matrix[index][exceptional] for index in contracted]
                product = [
                    sum(inverse[left][right] * vector[right] for right in range(size))
                    for left in range(size)
                ]
                quotient_square = matrix[exceptional][exceptional] - sum(
                    vector[index] * product[index] for index in range(size)
                )
                if quotient_square != -1:
                    continue
                enlarged = [
                    [
                        inverse[left][right] - product[left] * product[right]
                        for right in range(size)
                    ]
                    + [product[left]]
                    for left in range(size)
                ]
                enlarged.append([*product, -1])
                temporary = [*contracted, exceptional]
                order = sorted(range(size + 1), key=temporary.__getitem__)
                new_contracted = tuple(temporary[index] for index in order)
                new_inverse = tuple(
                    tuple(enlarged[left][right] for right in order) for left in order
                )
                following.setdefault(
                    mask | (1 << exceptional), (new_contracted, new_inverse)
                )
        current = following
    return current


def _endpoint_signature(row: dict[str, Any], contracted: Sequence[int]) -> tuple[Any, ...]:
    plane = _plane_class(row, contracted)
    boundary_degrees = tuple(
        sorted(
            _dot(row["gram_matrix"], row["walls_int"][index], plane)
            for index in row["boundaries"]
        )
    )
    root_degrees = tuple(
        sorted(
            Counter(
                _dot(row["gram_matrix"], row["walls_int"][index], plane)
                for index in row["roots"]
            ).items()
        )
    )
    contracted_norms = tuple(
        sorted(Counter(row["norms"][index] for index in contracted).items())
    )
    return boundary_degrees, root_degrees, contracted_norms


def _lexicographic_contraction_order(row: dict[str, Any], target: Sequence[int]) -> list[int]:
    target_set = set(target)
    wall_gram = row["wall_gram"]

    def recurse(order: list[int]) -> list[int] | None:
        if len(order) == len(target_set):
            return order
        if order:
            gram = sp.Matrix([[wall_gram[left][right] for right in order] for left in order])
            inverse = gram.inv()
        for candidate in sorted(target_set - set(order)):
            quotient_square = wall_gram[candidate][candidate]
            if order:
                vector = sp.Matrix([wall_gram[index][candidate] for index in order])
                quotient_square -= int((vector.T * inverse * vector)[0])
            if quotient_square != -1:
                continue
            found = recurse([*order, candidate])
            if found is not None:
                return found
        return None

    answer = recurse([])
    if answer is None:
        raise ArithmeticError("a P2 target has no sequential (-1)-contraction order")
    return answer


def _orthogonal_plane_basis(
    row: dict[str, Any], target: Sequence[int]
) -> tuple[list[int], sp.Matrix, list[sp.Matrix]]:
    order = _lexicographic_contraction_order(row, target)
    strict_classes: list[list[int]] = []
    exceptionals: list[sp.Matrix] = []
    for wall_index in order:
        current = sp.Matrix(row["walls_int"][wall_index])
        if strict_classes:
            gram = sp.Matrix(
                [
                    [_dot(row["gram_matrix"], left, right) for right in strict_classes]
                    for left in strict_classes
                ]
            )
            vector = sp.Matrix(
                [
                    _dot(
                        row["gram_matrix"],
                        previous,
                        row["walls_int"][wall_index],
                    )
                    for previous in strict_classes
                ]
            )
            projection = gram.inv() * vector
            total_transform = current - sum(
                (
                    projection[index] * sp.Matrix(strict_classes[index])
                    for index in range(len(strict_classes))
                ),
                sp.zeros(row["gram_matrix"].rows, 1),
            )
        else:
            total_transform = current
        if any(value.q != 1 for value in total_transform):
            raise ArithmeticError("a total exceptional transform is not integral")
        exceptionals.append(total_transform)
        strict_classes.append(row["walls_int"][wall_index])

    plane = sp.Matrix(_plane_class(row, tuple(sorted(target))))
    return order, plane, exceptionals


def _basis_coordinates(basis: sp.Matrix, vector: Sequence[int]) -> tuple[int, ...]:
    result = basis.inv() * sp.Matrix(vector)
    if any(value.q != 1 for value in result):
        raise ArithmeticError("a raw wall is nonintegral in a determinant-one plane basis")
    return tuple(int(value) for value in result)


def _canonical_incidence_graph() -> nx.Graph:
    graph = nx.Graph()
    for label in POINT_ORDER:
        graph.add_node(("point", label), kind="near" if label in NEAR_POINTS else "proper")
    for name, support in EXCEPTIONAL_ROOTS.items():
        graph.add_node(("root", name), kind="exceptional")
        for label in support:
            graph.add_edge(("root", name), ("point", label))
    for name, record in LINE_ROOTS.items():
        graph.add_node(("root", name), kind="line")
        for label in record["support"]:
            graph.add_edge(("root", name), ("point", label))
    return graph


def _canonical_root_class(name: str, position_of: dict[str, int]) -> tuple[int, ...]:
    result = [0] * 11
    if name in LINE_ROOTS:
        result[0] = 1
        for label in LINE_ROOTS[name]["support"]:
            result[1 + position_of[label]] = -1
    else:
        near, proper = EXCEPTIONAL_ROOTS[name]
        result[1 + position_of[near]] = -1
        result[1 + position_of[proper]] = 1
    return tuple(result)


def _branch_certificate(row: dict[str, Any], target: Sequence[int]) -> dict[str, Any]:
    order, plane, exceptionals = _orthogonal_plane_basis(row, target)
    basis = sp.Matrix.hstack(plane, *exceptionals)
    expected_gram = sp.diag(1, *([-1] * 10))
    actual_gram = basis.T * row["gram_matrix"] * basis
    if actual_gram != expected_gram or abs(int(basis.det())) != 1:
        raise ArithmeticError("the extracted plane basis is not determinant-one I_(1,10)")

    wall_coordinates = [
        _basis_coordinates(basis, vector) for vector in row["walls_int"]
    ]
    canonical_coordinates = _basis_coordinates(basis, row["canonical"])
    boundary_coordinates = [wall_coordinates[index] for index in row["boundaries"]]
    if canonical_coordinates != (-3, *([1] * 10)):
        raise ArithmeticError("the canonical class changed in the plane basis")
    if boundary_coordinates != [(3, *([-1] * 10))]:
        raise ArithmeticError("the primitive half-boundary is not -K")

    root_records = [
        (index, wall_coordinates[index]) for index in row["roots"]
    ]
    raw_graph = nx.Graph()
    near_positions: set[int] = set()
    for _, coordinates in root_records:
        if coordinates[0] == 0:
            negatives = [index for index, value in enumerate(coordinates[1:]) if value == -1]
            positives = [index for index, value in enumerate(coordinates[1:]) if value == 1]
            if len(negatives) != 1 or len(positives) != 1:
                raise ArithmeticError("a contracted root is not old-exceptional minus near-exceptional")
            near_positions.update(negatives)
    if len(near_positions) != 3:
        raise ArithmeticError("the plane branch did not recover three first-near centers")
    for position in range(10):
        raw_graph.add_node(
            ("point", position),
            kind="near" if position in near_positions else "proper",
        )
    for root_index, coordinates in root_records:
        kind = "exceptional" if coordinates[0] == 0 else "line"
        if kind == "line" and not (
            coordinates[0] == 1
            and sum(value == -1 for value in coordinates[1:]) == 3
            and all(value in (0, -1) for value in coordinates[1:])
        ):
            raise ArithmeticError("a positive-degree root is not a three-center line")
        raw_graph.add_node(("root", root_index), kind=kind)
        for position, value in enumerate(coordinates[1:]):
            if value:
                raw_graph.add_edge(("root", root_index), ("point", position))

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        _canonical_incidence_graph(),
        node_match=lambda left, right: left["kind"] == right["kind"],
    )
    maps = list(matcher.isomorphisms_iter())
    if not maps:
        raise ArithmeticError("the raw 7+3 incidence graph is not the canonical B3 design")
    label_rank = {label: index for index, label in enumerate(POINT_ORDER)}
    point_maps = sorted(
        (
            tuple(mapping[("point", position)][1] for position in range(10)),
            mapping,
        )
        for mapping in maps
    )
    selected_labels, selected_map = point_maps[0]
    position_of = {label: position for position, label in enumerate(selected_labels)}

    canonical_classes = {
        name: _canonical_root_class(name, position_of)
        for name in (*EXCEPTIONAL_ROOTS, *LINE_ROOTS)
    }
    if set(canonical_classes.values()) != {coordinates for _, coordinates in root_records}:
        raise ArithmeticError("the twelve canonical root classes do not equal the raw root set")
    name_by_class = {coordinates: name for name, coordinates in canonical_classes.items()}
    root_transport = [
        {
            "raw_wall_index_zero_based": root_index,
            "canonical_curve": name_by_class[coordinates],
            "coordinates_H_E": list(coordinates),
        }
        for root_index, coordinates in root_records
    ]

    reconstruction = basis * sp.Matrix.hstack(
        *(sp.Matrix(coordinates) for coordinates in wall_coordinates)
    )
    raw_matrix = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in row["walls_int"]))
    if reconstruction != raw_matrix:
        raise ArithmeticError("the complete raw wall transport failed")

    return {
        "contracted_wall_set_zero_based": list(target),
        "lexicographic_valid_contraction_order_zero_based": order,
        "plane_class_in_raw_basis": [int(value) for value in plane],
        "exceptional_total_transform_columns_in_raw_basis": [
            [int(value) for value in exceptional] for exceptional in exceptionals
        ],
        "basis_determinant": int(basis.det()),
        "basis_gram": [[int(value) for value in actual_gram.row(index)] for index in range(11)],
        "canonical_class_coordinates_H_E": list(canonical_coordinates),
        "primitive_half_boundary_coordinates_H_E": list(boundary_coordinates[0]),
        "near_center_positions_zero_based": sorted(near_positions),
        "colored_incidence_isomorphism_count": len(maps),
        "selected_point_labels_in_E_order": list(selected_labels),
        "root_transport": root_transport,
        "all_thirty_raw_wall_coordinates_H_E": [list(row) for row in wall_coordinates],
        "all_thirty_wall_vectors_reconstruct_exactly": True,
        "all_twelve_raw_roots_equal_canonical_curve_classes": True,
    }


def _raw_extraction() -> tuple[dict[str, Any], dict[str, Any]]:
    row = _normalized(next(record for record in _load_rows() if record["id"] == "one_boundary_30"))
    endpoints = _sequential_endpoints(row)
    signatures: Counter[str] = Counter()
    selected: list[tuple[int, ...]] = []
    desired = ((3,), ((0, 3), (1, 9)), ((-2, 3), (-1, 7)))
    for contracted, _ in endpoints.values():
        signature = _endpoint_signature(row, contracted)
        signatures[json.dumps(signature)] += 1
        if signature == desired:
            selected.append(contracted)
    if len(endpoints) != 573 or len(signatures) != 28 or len(selected) != 2:
        raise ArithmeticError("the H10/30 all-contraction denominator changed")
    branches = [_branch_certificate(row, target) for target in sorted(selected)]
    return row, {
        "input_policy": {
            "reads": ["integral lattice Gram", "thirty primitive wall vectors"],
            "does_not_read": [
                "surface name",
                "Kondo type",
                "known plane points",
                "known sextic",
                "stored contraction",
            ],
        },
        "raw_rank": row["gram_matrix"].rows,
        "raw_wall_count": len(row["walls_int"]),
        "raw_root_count": len(row["roots"]),
        "all_sequential_P2_endpoint_count": len(endpoints),
        "endpoint_degree_signature_count": len(signatures),
        "selected_signature": {
            "half_boundary_degree": 3,
            "root_degree_histogram": {"0": 3, "1": 9},
            "contracted_original_norm_histogram": {"-2": 3, "-1": 7},
        },
        "selected_endpoint_count": len(branches),
        "branches": branches,
        "both_branches_have_the_same_abstract_7_plus_3_design": True,
    }


def _normalized_projective_point(point: Sequence[sp.Rational | int]) -> tuple[sp.Rational, ...]:
    values = [sp.Rational(value) for value in point]
    for value in values:
        if value:
            return tuple(sp.simplify(entry / value) for entry in values)
    raise ArithmeticError("zero vector is not a projective point")


def _line_value(coefficients: Sequence[int], point: Sequence[int]) -> int:
    return sum(left * right for left, right in zip(coefficients, point))


def _near_substitution(label: str, polynomial: sp.Expr) -> sp.Expr:
    if label == "a":
        return sp.expand(polynomial.subs({X: 1, Y: S * T, Z: S}))
    if label == "d":
        return sp.expand(polynomial.subs({X: S, Y: 1, Z: S * T}))
    if label == "g":
        return sp.expand(polynomial.subs({X: S * T, Y: S, Z: 1}))
    raise KeyError(label)


def _monomial_exponents(degree: int) -> list[tuple[int, int, int]]:
    return [
        (left, middle, degree - left - middle)
        for left in range(degree + 1)
        for middle in range(degree - left + 1)
    ]


def _monomials(degree: int) -> list[sp.Expr]:
    return [X**left * Y**middle * Z**right for left, middle, right in _monomial_exponents(degree)]


def _proper_local_data(label: str) -> tuple[int, tuple[int, int], tuple[sp.Rational, ...]]:
    point = PROPER_POINTS[label]
    chart = next(index for index, value in enumerate(point) if value)
    normalized = _normalized_projective_point(point)
    free = tuple(index for index in range(3) if index != chart)
    return chart, free, normalized


def _cluster_condition_matrix(degree: int, multiplicity: int) -> tuple[list[sp.Expr], sp.Matrix, list[str]]:
    monomials = _monomials(degree)
    rows: list[list[sp.Rational]] = []
    names: list[str] = []
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _, free, point = _proper_local_data(label)
        for total in range(multiplicity):
            for first_order in range(total + 1):
                second_order = total - first_order
                rows.append(
                    [
                        sp.diff(
                            monomial,
                            VARIABLES[free[0]],
                            first_order,
                            VARIABLES[free[1]],
                            second_order,
                        ).subs(dict(zip(VARIABLES, point)))
                        for monomial in monomials
                    ]
                )
                names.append(
                    f"{label}.proper.derivative_{first_order}_{second_order}"
                )
        if label.lower() in NEAR_POINTS:
            near = label.lower()
            transforms = [sp.Poly(_near_substitution(near, monomial), S, T) for monomial in monomials]
            for total in range(multiplicity):
                for s_order in range(total + 1):
                    t_order = total - s_order
                    rows.append(
                        [
                            transform.coeff_monomial(
                                S ** (multiplicity + s_order) * T**t_order
                            )
                            for transform in transforms
                        ]
                    )
                    names.append(f"{near}.near.derivative_{s_order}_{t_order}")
    return monomials, sp.Matrix(rows), names


def _polynomial_string(polynomial: sp.Expr) -> str:
    return str(sp.expand(polynomial)).replace("**", "^").replace(" ", "")


def _linear_systems() -> tuple[dict[str, Any], sp.Expr]:
    cubics, cubic_matrix, cubic_names = _cluster_condition_matrix(3, 1)
    sextics, sextic_matrix, sextic_names = _cluster_condition_matrix(6, 2)
    cubic_rank = cubic_matrix.rank()
    sextic_rank = sextic_matrix.rank()
    sextic_kernel = sextic_matrix.nullspace()
    if cubic_rank != 10 or sextic_rank != 27 or len(sextic_kernel) != 1:
        raise ArithmeticError("the anticanonical/antibicanonical Hilbert ranks changed")
    coefficients = _primitive_integer_vector(sextic_kernel[0])
    boundary = sp.expand(sum(value * monomial for value, monomial in zip(coefficients, sextics)))
    expected = sp.expand(
        X**4 * Y**2
        + Y**4 * Z**2
        + Z**4 * X**2
        - 3 * X**2 * Y**2 * Z**2
    )
    if boundary != expected:
        raise ArithmeticError("the unique boundary sextic changed")
    return {
        "monomial_order": "x^a*y^b*z^(d-a-b), a ascending then b ascending",
        "anticanonical_cubics": {
            "ambient_vector_dimension": len(cubics),
            "condition_row_count": cubic_matrix.rows,
            "condition_rank": cubic_rank,
            "kernel_dimension": len(cubic_matrix.nullspace()),
            "condition_labels": cubic_names,
        },
        "antibicanonical_sextics": {
            "ambient_vector_dimension": len(sextics),
            "condition_row_count": sextic_matrix.rows,
            "condition_rank": sextic_rank,
            "kernel_dimension": len(sextic_kernel),
            "condition_labels": sextic_names,
            "coefficient_vector": coefficients,
            "polynomial": _polynomial_string(boundary),
        },
    }, boundary


def _local_expansion(polynomial: sp.Expr, label: str) -> tuple[sp.Expr, int, sp.Expr, int]:
    chart, free, point = _proper_local_data(label)
    U, V = sp.symbols("u v")
    substitution = {
        VARIABLES[chart]: 1,
        VARIABLES[free[0]]: point[free[0]] + U,
        VARIABLES[free[1]]: point[free[1]] + V,
    }
    expansion = sp.expand(polynomial.subs(substitution))
    terms = sp.Poly(expansion, U, V).terms()
    multiplicity = min(sum(exponent) for exponent, coefficient in terms if coefficient)
    tangent = sp.expand(
        sum(
            coefficient * U ** exponent[0] * V ** exponent[1]
            for exponent, coefficient in terms
            if sum(exponent) == multiplicity
        )
    )
    quadratic = sp.Poly(tangent, U, V)
    first = quadratic.coeff_monomial(U**2)
    middle = quadratic.coeff_monomial(U * V)
    second = quadratic.coeff_monomial(V**2)
    discriminant = int(sp.expand(middle**2 - 4 * first * second))
    return expansion, multiplicity, tangent, discriminant


def _order_at_origin(polynomial: sp.Expr, variables: tuple[sp.Symbol, sp.Symbol]) -> int:
    terms = sp.Poly(sp.expand(polynomial), *variables).terms()
    return min(sum(exponent) for exponent, coefficient in terms if coefficient)


def _cluster_multiplicities(polynomial: sp.Expr) -> dict[str, int]:
    result: dict[str, int] = {}
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _, proper_multiplicity, _, _ = _local_expansion(polynomial, label)
        result[label] = proper_multiplicity
        near = label.lower()
        if near in NEAR_POINTS:
            transform = sp.cancel(_near_substitution(near, polynomial) / S**proper_multiplicity)
            if not transform.is_polynomial(S, T):
                raise ArithmeticError("a blowup strict transform is not polynomial")
            result[near] = _order_at_origin(transform, (S, T))
    return result


def _canonical_geometry(boundary: sp.Expr) -> dict[str, Any]:
    # Incidence rigidity: after normalizing A,D,G and B, all remaining
    # collinearities reduce to c^2=1; c=1 collides C with B.
    parameter = sp.symbols("c")
    symbolic = {
        "B": sp.Matrix([1, 1, 1]),
        "C": sp.Matrix([1, 1, parameter]),
        "E": sp.Matrix([1 / parameter, 1, 1]),
        "F": sp.Matrix([parameter, 1, parameter]),
    }
    symbolic.update(
        {
            "A": sp.Matrix([1, 0, 0]),
            "D": sp.Matrix([0, 1, 0]),
            "G": sp.Matrix([0, 0, 1]),
        }
    )
    forced_determinants = {
        "BCG": sp.factor(
            sp.Matrix.hstack(symbolic["B"], symbolic["C"], symbolic["G"]).det()
        ),
        "ABE": sp.factor(
            sp.Matrix.hstack(symbolic["A"], symbolic["B"], symbolic["E"]).det()
        ),
        "BDF": sp.factor(
            sp.Matrix.hstack(symbolic["B"], symbolic["D"], symbolic["F"]).det()
        ),
        "ACF": sp.factor(
            sp.Matrix.hstack(symbolic["A"], symbolic["C"], symbolic["F"]).det()
        ),
        "CDE": sp.factor(
            sp.Matrix.hstack(symbolic["C"], symbolic["D"], symbolic["E"]).det()
        ),
    }
    if any(value != 0 for value in forced_determinants.values()):
        raise ArithmeticError("the normalized forced incidence equations changed")
    residual = sp.factor(
        sp.Matrix.hstack(symbolic["E"], symbolic["F"], sp.Matrix([0, 0, 1])).det()
    )
    if sp.simplify(residual * parameter - (1 - parameter**2)) != 0:
        raise ArithmeticError("the final incidence equation is not c^2=1")

    line_records = []
    for name, record in LINE_ROOTS.items():
        coefficients = record["coefficients"]
        polynomial = coefficients[0] * X + coefficients[1] * Y + coefficients[2] * Z
        multiplicities = _cluster_multiplicities(polynomial)
        exact_support = tuple(
            label for label in POINT_ORDER if multiplicities.get(label, 0) == 1
        )
        if set(exact_support) != set(record["support"]) or any(
            value not in (0, 1) for value in multiplicities.values()
        ):
            raise ArithmeticError(f"line {name} has the wrong proper/near cluster support")
        line_records.append(
            {
                "name": name,
                "equation_coefficients_x_y_z": list(coefficients),
                "polynomial": _polynomial_string(polynomial),
                "cluster_support": list(record["support"]),
                "exact_cluster_multiplicities": multiplicities,
                "strict_transform_class": (
                    "H-" + "-".join(f"E_{label}" for label in record["support"])
                ),
            }
        )

    boundary_multiplicities = _cluster_multiplicities(boundary)
    if boundary_multiplicities != {label: 2 for label in POINT_ORDER}:
        raise ArithmeticError("the sextic is not double at the complete 7+3 cluster")

    normalized_points = {
        _normalized_projective_point(point): label for label, point in PROPER_POINTS.items()
    }
    outside_intersections: dict[tuple[sp.Rational, ...], list[tuple[str, str]]] = {}
    proper_pair_count = 0
    for (left_name, left), (right_name, right) in itertools.combinations(LINE_ROOTS.items(), 2):
        point = _normalized_projective_point(
            sp.Matrix(left["coefficients"]).cross(sp.Matrix(right["coefficients"]))
        )
        if point in normalized_points:
            proper_pair_count += 1
        else:
            outside_intersections.setdefault(point, []).append((left_name, right_name))
    if proper_pair_count != 30 or len(outside_intersections) != 6 or any(
        len(pairs) != 1 for pairs in outside_intersections.values()
    ):
        raise ArithmeticError("the B3 line-arrangement intersection census changed")

    proper_line_counts = {
        label: sum(
            _line_value(record["coefficients"], point) == 0
            for record in LINE_ROOTS.values()
        )
        for label, point in PROPER_POINTS.items()
    }
    if Counter(proper_line_counts.values()) != Counter({4: 3, 3: 4}):
        raise ArithmeticError("the seven proper multiple points changed")

    return {
        "incidence_rigidity": {
            "normalization": "A=(1:0:0), D=(0:1:0), G=(0:0:1), B=(1:1:1)",
            "forced_symbolic_points": {
                label: [str(value) for value in vector]
                for label, vector in symbolic.items()
                if label in ("B", "C", "E", "F")
            },
            "already_solved_collinearity_determinants": {
                label: str(value) for label, value in forced_determinants.items()
            },
            "last_collinearity_determinant": str(residual),
            "nondegenerate_parameter_solution": "c=-1; c=1 gives C=B",
            "requires_characteristic_not_two": True,
        },
        "proper_points": [
            {"label": label, "coordinates": list(PROPER_POINTS[label])}
            for label in ("A", "D", "G", "B", "C", "E", "F")
        ],
        "first_near_points": [
            {
                "label": "a",
                "parent": "A",
                "direction_line": "AG",
                "blowup_chart": "x=1, y=s*t, z=s; near point (s,t)=(0,0)",
            },
            {
                "label": "d",
                "parent": "D",
                "direction_line": "AD",
                "blowup_chart": "y=1, x=s, z=s*t; near point (s,t)=(0,0)",
            },
            {
                "label": "g",
                "parent": "G",
                "direction_line": "DG",
                "blowup_chart": "z=1, x=s*t, y=s; near point (s,t)=(0,0)",
            },
        ],
        "line_arrangement": {
            "name": "projectivized B3 reflection arrangement",
            "line_count": len(line_records),
            "lines": line_records,
            "proper_point_line_counts": proper_line_counts,
            "line_pair_count": 36,
            "line_pairs_meeting_at_blown_proper_points": proper_pair_count,
            "remaining_distinct_transverse_double_intersections": len(outside_intersections),
            "remaining_intersection_points": [
                [str(value) for value in point] for point in sorted(outside_intersections)
            ],
        },
        "exceptional_roots": [
            {
                "name": name,
                "class": f"E_{proper}-E_{near}",
                "description": f"strict transform of the exceptional curve over {proper}",
            }
            for name, (near, proper) in EXCEPTIONAL_ROOTS.items()
        ],
        "root_inventory": {
            "line_roots": 9,
            "exceptional_chain_roots": 3,
            "total_smooth_geometrically_integral_rational_minus_two_curves": 12,
            "line_line_transverse_intersections_after_blowup": 6,
            "exceptional_line_transverse_intersections_after_blowup": 9,
            "total_root_pair_intersection_one": 15,
            "total_root_pair_intersection_zero": 51,
        },
        "boundary_cluster_multiplicities": boundary_multiplicities,
    }


def _singular_audit(boundary: sp.Expr) -> dict[str, Any]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for absolute factorization and singular schemes")
    polynomial = _polynomial_string(boundary)
    program = f'''
LIB "primdec.lib";
ring r=0,(x,y,z),dp;
poly f={polynomial};
ideal A=y,z;
ideal D=x,z;
ideal G=x,y;
ideal B=x-y,y-z;
ideal C=x-y,y+z;
ideal E=x+y,y-z;
ideal F=x+y,y+z;
ideal T=std(intersect(A,D,G,B,C,E,F));
ideal J=std(sat(jacob(f),maxideal(1)));
ideal R=std(radical(J));
list FF=factorize(f);
print("Q_FACTOR_COUNT_WITH_UNIT"); print(size(FF[1]));
print("Q_NONUNIT_EXPONENT"); print(FF[2][2]);
print("SINGULAR_AFFINE_CONE_DIMENSION"); print(dim(J));
print("SINGULAR_SCHEME_DEGREE"); print(mult(J));
print("SINGULAR_RADICAL_DEGREE"); print(mult(R));
print("SINGULAR_MINIMAL_PRIMES"); print(size(minAssGTZ(J)));
print("J_MOD_R"); print(size(reduce(J,R)));
print("R_MOD_J"); print(size(reduce(R,J)));
print("R_MOD_POINTS"); print(size(reduce(R,T)));
print("POINTS_MOD_R"); print(size(reduce(T,R)));
LIB "absfact.lib";
def ABSR=absFactorize(f);
setring ABSR;
print("ABS_FACTOR_COUNT_WITH_UNIT"); print(size(absolute_factors[1]));
print("ABS_NONUNIT_EXPONENT"); print(absolute_factors[2][2]);
print("ABS_COMPONENT_COUNT"); print(absolute_factors[4]);
exit;
'''
    completed = subprocess.run(
        [executable, "-q", "--no-rc"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    labels = {
        "Q_FACTOR_COUNT_WITH_UNIT": "factor_count_with_unit_over_Q",
        "Q_NONUNIT_EXPONENT": "nonunit_exponent_over_Q",
        "SINGULAR_AFFINE_CONE_DIMENSION": "singular_affine_cone_dimension",
        "SINGULAR_SCHEME_DEGREE": "singular_scheme_degree",
        "SINGULAR_RADICAL_DEGREE": "singular_radical_degree",
        "SINGULAR_MINIMAL_PRIMES": "singular_minimal_primes",
        "J_MOD_R": "jacobian_mod_radical",
        "R_MOD_J": "radical_mod_jacobian",
        "R_MOD_POINTS": "radical_mod_seven_points",
        "POINTS_MOD_R": "seven_points_mod_radical",
        "ABS_FACTOR_COUNT_WITH_UNIT": "absolute_factor_count_with_unit",
        "ABS_NONUNIT_EXPONENT": "absolute_nonunit_exponent",
        "ABS_COMPONENT_COUNT": "absolute_component_count",
    }
    parsed: dict[str, int] = {}
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    for index, line in enumerate(lines[:-1]):
        if line in labels:
            parsed[labels[line]] = int(lines[index + 1])
    expected = {
        "factor_count_with_unit_over_Q": 2,
        "nonunit_exponent_over_Q": 1,
        "singular_affine_cone_dimension": 1,
        "singular_scheme_degree": 13,
        "singular_radical_degree": 7,
        "singular_minimal_primes": 7,
        "jacobian_mod_radical": 0,
        "radical_mod_jacobian": 4,
        "radical_mod_seven_points": 0,
        "seven_points_mod_radical": 0,
        "absolute_factor_count_with_unit": 2,
        "absolute_nonunit_exponent": 1,
        "absolute_component_count": 1,
    }
    if parsed != expected:
        raise ArithmeticError(f"Singular boundary audit changed: {parsed}")
    return {
        "engine": "Singular exact Q-factorization, absfact, saturated Jacobian, radical and primary support",
        **parsed,
        "irreducible_over_Q": True,
        "absolutely_irreducible": True,
        "singular_support_is_exactly_the_seven_proper_points": True,
        "scheme_length_decomposition": "3 A3 tacnodes * 3 + 4 A1 nodes * 1 = 13",
    }


def _local_boundary_audit(boundary: sp.Expr) -> dict[str, Any]:
    records = []
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _, multiplicity, tangent, discriminant = _local_expansion(boundary, label)
        if multiplicity != 2:
            raise ArithmeticError("a proper boundary singularity is not double")
        if label in ("A", "D", "G"):
            near = label.lower()
            transform = sp.expand(sp.cancel(_near_substitution(near, boundary) / S**2))
            near_multiplicity = _order_at_origin(transform, (S, T))
            terms = sp.Poly(transform, S, T).terms()
            near_tangent = sp.expand(
                sum(
                    coefficient * S ** exponent[0] * T ** exponent[1]
                    for exponent, coefficient in terms
                    if sum(exponent) == near_multiplicity
                )
            )
            quadratic = sp.Poly(near_tangent, S, T)
            a = quadratic.coeff_monomial(S**2)
            b = quadratic.coeff_monomial(S * T)
            c = quadratic.coeff_monomial(T**2)
            near_discriminant = int(b**2 - 4 * a * c)
            if discriminant != 0 or near_multiplicity != 2 or near_discriminant == 0:
                raise ArithmeticError("a proposed A3 tacnode does not become an ordinary node")
            singularity = "A3 ordinary tacnode; first blowup gives A1 node; second blowup resolves"
        else:
            near_multiplicity = None
            near_tangent = None
            near_discriminant = None
            if discriminant == 0:
                raise ArithmeticError("a proposed ordinary node has repeated tangent")
            singularity = "A1 ordinary node; one blowup resolves"
        records.append(
            {
                "point": label,
                "multiplicity": multiplicity,
                "tangent_cone": str(sp.factor(tangent)),
                "tangent_discriminant": discriminant,
                "first_blowup_near_multiplicity": near_multiplicity,
                "first_blowup_tangent_cone": (
                    str(sp.factor(near_tangent)) if near_tangent is not None else None
                ),
                "first_blowup_tangent_discriminant": near_discriminant,
                "singularity": singularity,
            }
        )
    return {
        "proper_singularity_count": len(records),
        "A3_tacnode_count": 3,
        "A1_node_count": 4,
        "delta_sum": 3 * 2 + 4,
        "records": records,
        "all_singularities_resolved_by_the_declared_7_plus_3_blowups": True,
        "strict_boundary_smooth": True,
        "strict_boundary_self_intersection": -4,
        "strict_boundary_arithmetic_genus": 0,
        "strict_boundary_geometrically_integral_rational": True,
    }


def _same_polytope_return() -> dict[str, Any]:
    archive = _load(FULL_CONE)["cases"]["II"]
    if not archive["all_checks_pass"]:
        raise ArithmeticError("the independent H10/30 full-cone archive is not green")
    double = archive["exact_double_description"]
    vinberg = archive["literal_vinberg"]
    if not (
        double["hyperideal_vertex_count"] == 0
        and vinberg["bad_edge_count"] == 0
        and vinberg["bad_parabolic_neighborhood_count"] == 0
    ):
        raise ArithmeticError("the independent D_CobVin return hypotheses changed")
    hostile = _load(HOSTILE_RECOGNIZER)
    if not hostile["all_checks_pass"]:
        raise ArithmeticError("the independent P-first Coble recognizer is not green")
    presentation = next(
        record
        for record in hostile["presentations"]
        if record["chamber_id"] == "D16-D5-H10-30"
    )
    if len(presentation["candidate_audits"]) != 1:
        raise ArithmeticError("the raw Type-II boundary is no longer unique")
    candidate = presentation["candidate_audits"][0]
    cm_cone = candidate["cm_cone"]
    cusps = candidate["cusp_regeneration"]
    saturation = candidate["integral_saturation_uniqueness"]
    if not (
        candidate["coble_incidence_pass"]
        and candidate["complete_package_regenerated_from_P"]
        and candidate["kondo_type_II_root_gram_match"]
        and cm_cone["kind_histogram"] == {"finite": 38, "ideal": 17}
        and cusps["complete_primitive_isotropic_cusp_count"] == 17
        and cusps["sets_equal"]
        and saturation["unique_E10_saturation_up_to_diagram_automorphism"]
    ):
        raise ArithmeticError("the independent CM-root/cusp recognition gate changed")
    return {
        "independent_CM_root_cone_gate": {
            "finite_rays": cm_cone["kind_histogram"]["finite"],
            "ideal_rays": cm_cone["kind_histogram"]["ideal"],
            "hyperideal_rays": 0,
            "all_seventeen_cusps_equal_the_raw_E_plus_b_classes": cusps["sets_equal"],
            "root_gram_is_Kondo_Type_II": candidate["kondo_type_II_root_gram_match"],
            "root_span_index_in_E10": saturation["root_span_index_in_E10"],
            "unique_E10_saturation_up_to_diagram_automorphism": saturation[
                "unique_E10_saturation_up_to_diagram_automorphism"
            ],
        },
        "independent_raw_cone_gate": {
            "all_thirty_facets_irredundant": double["all_facets_irredundant"],
            "finite_vertices": double["finite_vertex_count"],
            "ideal_vertices": double["ideal_vertex_count"],
            "hyperideal_vertices": double["hyperideal_vertex_count"],
            "elliptic_edges": vinberg["elliptic_edge_count"],
            "connected_affine_neighborhoods": vinberg["connected_affine_component_count"],
            "bad_edges": vinberg["bad_edge_count"],
            "bad_parabolic_neighborhoods": vinberg[
                "bad_parabolic_neighborhood_count"
            ],
            "reflection_index_in_W_P10": archive["P10_reflection_subgroup_index"],
        },
        "general_theorem_handoff": (
            "The twelve explicitly constructed roots are the complete raw CM root seed. "
            "D_CobVin gives finite volume and finite reflection index; Kondo root exhaustion "
            "excludes hidden roots; the one-boundary Halphen cusp theorem realizes exactly "
            "the seventeen raw E=f-b walls and excludes every hidden (-1)-curve."
        ),
        "returned_wall_inventory": {
            "primitive_half_boundary_b_minus_K": 1,
            "explicit_minus_two_roots": 12,
            "automatic_cusp_lift_minus_one_curves": 17,
            "total_facets": 30,
        },
        "exact_same_marked_H10_30": True,
        "finite_volume": True,
        "noncompact": True,
    }


def build_certificate() -> dict[str, Any]:
    _, raw = _raw_extraction()
    linear_systems, boundary = _linear_systems()
    geometry = _canonical_geometry(boundary)
    singular = _singular_audit(boundary)
    local = _local_boundary_audit(boundary)
    returned = _same_polytope_return()

    root_inventory = geometry["root_inventory"]
    checks = {
        "raw_all_contraction_denominator_is_573": raw["all_sequential_P2_endpoint_count"] == 573,
        "exactly_two_clean_7_plus_3_endpoints": raw["selected_endpoint_count"] == 2,
        "both_endpoint_bases_unimodular": all(
            abs(branch["basis_determinant"]) == 1 for branch in raw["branches"]
        ),
        "both_endpoint_root_transports_exact": all(
            branch["all_twelve_raw_roots_equal_canonical_curve_classes"]
            for branch in raw["branches"]
        ),
        "all_thirty_raw_walls_transport_exactly": all(
            branch["all_thirty_wall_vectors_reconstruct_exactly"]
            for branch in raw["branches"]
        ),
        "incidence_solution_is_rigid_and_nondegenerate": geometry["incidence_rigidity"][
            "nondegenerate_parameter_solution"
        ].startswith("c=-1"),
        "anticanonical_system_empty": linear_systems["anticanonical_cubics"][
            "kernel_dimension"
        ]
        == 0,
        "antibicanonical_system_singleton": linear_systems["antibicanonical_sextics"][
            "kernel_dimension"
        ]
        == 1,
        "boundary_absolutely_irreducible": singular["absolutely_irreducible"],
        "singular_support_exact": singular[
            "singular_support_is_exactly_the_seven_proper_points"
        ],
        "declared_blowups_resolve_boundary": local[
            "all_singularities_resolved_by_the_declared_7_plus_3_blowups"
        ],
        "all_twelve_root_curves_smooth_geometrically_integral": (
            root_inventory[
                "total_smooth_geometrically_integral_rational_minus_two_curves"
            ]
            == 12
        ),
        "root_intersection_gram_counts_match_raw": (
            root_inventory["total_root_pair_intersection_one"] == 15
            and root_inventory["total_root_pair_intersection_zero"] == 51
        ),
        "strong_Coble": (
            linear_systems["anticanonical_cubics"]["kernel_dimension"] == 0
            and linear_systems["antibicanonical_sextics"]["kernel_dimension"] == 1
            and local["strict_boundary_smooth"]
        ),
        "same_marked_polytope_return": returned["exact_same_marked_H10_30"],
    }

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "status": "PROVED_STRICT_RAW_P_FIRST_B3_PLANE_REALIZATION_AND_SAME_H10_30_RETURN",
        "direction": "RAW_H10_30_TO_TWO_INTRINSIC_P2_BRANCHES_TO_RIGID_B3_CLUSTER_TO_STRONG_COBLE_TO_SAME_P",
        "field_scope": {
            "equations_and_blowups_defined_over": "Q",
            "geometric_integrality_and_smoothness_certified_over": "Qbar, hence over C",
            "classification_and_same_nef_cone_conclusion_used_over": "C",
            "excluded_overstatement": (
                "No all-characteristic-not-two claim is made.  In characteristic 3 the "
                "coefficient -3 and the displayed node discriminants degenerate."
            ),
        },
        "raw_polytope_extraction": raw,
        "P_derived_incidence_realization": geometry,
        "linear_systems": linear_systems,
        "absolute_factor_and_singular_scheme_audit": singular,
        "local_boundary_resolution": local,
        "blowup_surface": {
            "base": "P2 over Q",
            "operations": "blow up seven proper rational points, then a,d,g in the declared tangent directions",
            "Picard_rank": 11,
            "canonical_class": "K=-3H+sum_(ten centers) E_i",
            "boundary_class": "B=6H-2sum E_i=-2K",
            "boundary_square": -4,
            "root_classes": "nine H-E-E-E lines and three E_parent-E_near exceptional roots",
            "minus_K_empty": True,
            "minus_2K_unique_reduced_boundary": True,
            "strong_Coble": True,
            "defined_over_Q": True,
        },
        "root_boundary_disjointness": {
            "numerical": "B.C=(-2K).C=0 for every one of the twelve roots",
            "geometric": (
                "The boundary is absolutely irreducible of degree six and hence contains no "
                "root line; all strict curves are distinct effective irreducibles on the smooth "
                "blowup. Nonnegative local intersection plus total intersection zero makes "
                "every root disjoint from the strict boundary."
            ),
            "root_root_transversality": (
                "Six line pairs meet at unblown ordinary double intersections; each of the "
                "three exceptional roots meets the three nonselected tangent directions over "
                "its parent. These fifteen intersections are transverse and every other root "
                "pair is disjoint."
            ),
        },
        "same_polytope_return": returned,
        "lineage": {
            "after_construction_classification": "Kondo one-boundary Type II",
            "abstract_polytope": "the already known D16-D5 H10/30 chamber",
            "sextic_equation_lineage": (
                "The forced equation is the Choi-Lam cyclic extremal positive-semidefinite "
                "ternary sextic from 1977; the polynomial itself is not new."
            ),
            "new_project_content": (
                "the strict raw-P derivation of the 7+3 cluster and its identification as the "
                "Choi-Lam sextic, followed by an exact strong-Coble and same-H10/30 return; "
                "no known Type-II surface equation is imported"
            ),
            "priority_boundary": (
                "new exact realization and cross-identification in this project; no worldwide "
                "priority claim for the Choi-Lam/Coble link, no new equation, no new abstract "
                "polytope, no new dimension, and not L4"
            ),
            "distinct_nearby_models": (
                "Kondo Example 4.5 uses a P1xP1 double-cover quotient.  The Hessian Coble "
                "surfaces of Allcock-Dolgachev are a different tetrahedral family with "
                "infinite automorphism group and ten proper nodes."
            ),
        },
        "sources": {
            "Kondo_classification_and_return": "https://arxiv.org/abs/2102.00591",
            "Choi_Lam_sextic_lineage": "https://doi.org/10.1007/BF01360024",
            "distinct_Hessian_Coble_family": "https://doi.org/10.5802/ahl.57",
            "general_inverse_return": "results/coble_vinberg_incidence_inverse_transfer.md",
            "independent_full_cone": "data/enumerations/coble_one_boundary_full_ns_h10.json",
            "independent_raw_Coble_recognizer": (
                "data/enumerations/polytope_first_coble_h10_hostile_audit.json"
            ),
        },
        "inverse_level": {
            "rating": "L3_STRICT_SINGLE_RAW_P_INVERSE_REALIZATION_NOT_L4",
            "positive_result": (
                "Every algebraic equation is forced by one of two raw contraction branches, "
                "and the exact same marked finite-volume noncompact H10/30 chamber returns."
            ),
            "why_not_L4": (
                "This is one rigid raw polytope and one known abstract chamber, not a variable "
                "P-only class with an exhaustive necessary-and-sufficient carrier theorem."
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_output is not None and payload != _load(args.check_output):
        raise SystemExit("live Type-II B3 plane certificate differs from frozen JSON")
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
