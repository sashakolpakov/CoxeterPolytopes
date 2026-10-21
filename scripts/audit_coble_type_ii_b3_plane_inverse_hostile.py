#!/usr/bin/env python3
"""Hostile, independent audit of the raw H10/30 -> B3 plane construction.

This file deliberately does *not* import ``certify_coble_type_ii_b3_plane_inverse``
or any of its helpers.  It starts with the archived rank-11 lattice and thirty
wall vectors, independently enumerates every sequential P2 endpoint, rebuilds
the two plane bases and the B3 cluster, and checks the surface and cone with
exact arithmetic.  The last section records precisely what this one chamber
does, and does not, prove about a general polytope-to-surface transition.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
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


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"
TARGET = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json"
CONE = ROOT / "data/enumerations/shimada_disc16_d5_h10_30_exact_cone_certificate.json"
FACES = ROOT / "data/enumerations/shimada_disc16_d5_h10_30_vinberg_face_certificate.json"
RECOGNIZER = ROOT / "data/enumerations/polytope_first_coble_h10_hostile_audit.json"
OUTPUT = ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse_hostile_audit.json"
SCHEMA = "coxeter.coble_type_ii_b3_plane_inverse.hostile_audit.v1"

X, Y, Z = sp.symbols("x y z")
S, T = sp.symbols("s t")
U, V = sp.symbols("u v")

POINT_ORDER = ("a", "A", "d", "D", "g", "G", "B", "C", "E", "F")
NEAR = {"a", "d", "g"}
PROPER = {
    "A": (1, 0, 0),
    "D": (0, 1, 0),
    "G": (0, 0, 1),
    "B": (1, 1, 1),
    "C": (1, 1, -1),
    "E": (-1, 1, 1),
    "F": (-1, 1, -1),
}
CHAINS = {
    "EA_minus_Ea": ("a", "A"),
    "ED_minus_Ed": ("d", "D"),
    "EG_minus_Eg": ("g", "G"),
}
LINES = {
    "L_DG": ((1, 0, 0), ("g", "G", "D")),
    "L_AG": ((0, 1, 0), ("a", "A", "G")),
    "L_AD": ((0, 0, 1), ("d", "D", "A")),
    "L_EFG": ((1, 1, 0), ("E", "F", "G")),
    "L_BCG": ((1, -1, 0), ("B", "C", "G")),
    "L_ACF": ((0, 1, 1), ("A", "C", "F")),
    "L_ABE": ((0, 1, -1), ("A", "B", "E")),
    "L_CDE": ((1, 0, 1), ("C", "D", "E")),
    "L_BDF": ((1, 0, -1), ("B", "D", "F")),
}


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(gram: sp.Matrix, left: Sequence[Any] | sp.Matrix, right: Sequence[Any] | sp.Matrix) -> sp.Rational:
    return sp.Rational((sp.Matrix(left).T * gram * sp.Matrix(right))[0])


def primitive(vector: sp.Matrix, orient: bool = True) -> tuple[int, ...]:
    denominator = math.lcm(*(int(sp.denom(value)) for value in vector))
    entries = [int(value * denominator) for value in vector]
    divisor = math.gcd(*(abs(value) for value in entries))
    if divisor == 0:
        raise ArithmeticError("zero vector cannot be primitive")
    entries = [value // divisor for value in entries]
    if orient:
        first = next(value for value in entries if value)
        if first < 0:
            entries = [-value for value in entries]
    return tuple(entries)


def raw_input() -> dict[str, Any]:
    archive = load(RAW)
    case = archive["cases"]["II"]
    gram = sp.Matrix(case["lattice_gram"])
    walls = tuple(sp.Matrix(row) for row in case["wall_vectors"])
    wall_gram = sp.Matrix([[dot(gram, left, right) for right in walls] for left in walls])
    if wall_gram != sp.Matrix(case["negative_wall_gram"]):
        raise ArithmeticError("wall Gram does not reconstruct from raw vectors")
    charpoly = gram.charpoly().as_poly()
    inertia = (
        int(charpoly.count_roots(0, sp.oo)),
        int(charpoly.count_roots(-sp.oo, 0)),
        11 - int(charpoly.count_roots(0, sp.oo)) - int(charpoly.count_roots(-sp.oo, 0)),
    )
    if gram.det() != 1 or inertia != (1, 10, 0) or wall_gram.rank() != 11:
        raise ArithmeticError("raw H10 lattice invariants changed")
    norms = tuple(int(wall_gram[index, index]) for index in range(30))

    def characteristic(vector: sp.Matrix) -> bool:
        return all(
            (int(dot(gram, vector, sp.eye(11).col(index))) - int(gram[index, index])) % 2 == 0
            for index in range(11)
        )

    characteristic_walls = tuple(index for index in range(30) if norms[index] == -1 and characteristic(walls[index]))
    if characteristic_walls != (0,):
        raise ArithmeticError("boundary is not the unique characteristic norm-minus-one wall")
    boundary = walls[0]
    roots = tuple(index for index, norm in enumerate(norms) if norm == -2)
    exceptional = tuple(index for index, norm in enumerate(norms) if norm == -1 and index != 0)
    cusp_classes = tuple(primitive(walls[index] + boundary, orient=False) for index in exceptional)
    if not (
        roots == tuple(range(1, 13))
        and exceptional == tuple(range(13, 30))
        and len(set(cusp_classes)) == 17
        and all(dot(gram, boundary, walls[index]) == 0 for index in roots)
        and all(dot(gram, boundary, walls[index]) == 1 for index in exceptional)
        and all(dot(gram, cusp, cusp) == 0 and dot(gram, cusp, boundary) == 0 for cusp in cusp_classes)
    ):
        raise ArithmeticError("raw b/root/cusp-lift partition changed")
    return {
        "archive": archive,
        "case": case,
        "gram": gram,
        "walls": walls,
        "wall_gram": wall_gram,
        "norms": norms,
        "boundary": boundary,
        "roots": roots,
        "exceptional": exceptional,
        "cusps": cusp_classes,
        "inertia": inertia,
    }


def sequential_endpoints(data: dict[str, Any]) -> tuple[dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]], list[int]]:
    matrix: sp.Matrix = data["wall_gram"]
    current: dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = {0: (tuple(), tuple())}
    trail = [1]
    for _depth in range(10):
        following: dict[int, tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]] = {}
        for mask, (chosen, inverse) in current.items():
            size = len(chosen)
            for candidate in range(1, 30):
                if mask & (1 << candidate):
                    continue
                vector = [int(matrix[index, candidate]) for index in chosen]
                product = [sum(inverse[left][right] * vector[right] for right in range(size)) for left in range(size)]
                schur = int(matrix[candidate, candidate]) - sum(vector[index] * product[index] for index in range(size))
                if schur != -1:
                    continue
                enlarged = [
                    [inverse[left][right] - product[left] * product[right] for right in range(size)] + [product[left]]
                    for left in range(size)
                ]
                enlarged.append([*product, -1])
                temporary = [*chosen, candidate]
                order = sorted(range(size + 1), key=temporary.__getitem__)
                sorted_chosen = tuple(temporary[index] for index in order)
                sorted_inverse = tuple(tuple(enlarged[left][right] for right in order) for left in order)
                following.setdefault(mask | (1 << candidate), (sorted_chosen, sorted_inverse))
        current = following
        trail.append(len(current))
    return current, trail


def plane_class(data: dict[str, Any], chosen: Sequence[int]) -> tuple[int, ...]:
    gram: sp.Matrix = data["gram"]
    walls: tuple[sp.Matrix, ...] = data["walls"]
    equations = sp.Matrix.hstack(*(walls[index] for index in chosen)).T * gram
    kernel = equations.nullspace()
    if len(kernel) != 1:
        raise ArithmeticError("a ten-contraction complement is not a line")
    vector = list(primitive(kernel[0]))
    if dot(gram, data["boundary"], vector) < 0:
        vector = [-value for value in vector]
    if dot(gram, vector, vector) != 1:
        raise ArithmeticError("endpoint complement is not a plane class")
    return tuple(vector)


def endpoint_signature(data: dict[str, Any], chosen: Sequence[int]) -> tuple[Any, ...]:
    h = plane_class(data, chosen)
    return (
        (int(dot(data["gram"], data["boundary"], h)),),
        tuple(sorted(Counter(int(dot(data["gram"], data["walls"][index], h)) for index in data["roots"]).items())),
        tuple(sorted(Counter(data["norms"][index] for index in chosen).items())),
    )


def contraction_order(data: dict[str, Any], target: Sequence[int]) -> list[int]:
    wall_gram: sp.Matrix = data["wall_gram"]
    wanted = set(target)

    def visit(order: list[int]) -> list[int] | None:
        if len(order) == 10:
            return order
        inverse = wall_gram.extract(order, order).inv() if order else None
        for candidate in sorted(wanted - set(order)):
            square = wall_gram[candidate, candidate]
            if order:
                pair = sp.Matrix([wall_gram[index, candidate] for index in order])
                square -= (pair.T * inverse * pair)[0]
            if square == -1:
                answer = visit([*order, candidate])
                if answer is not None:
                    return answer
        return None

    answer = visit([])
    if answer is None:
        raise ArithmeticError("selected endpoint has no valid sequential order")
    return answer


def canonical_graph() -> nx.Graph:
    graph = nx.Graph()
    for label in POINT_ORDER:
        graph.add_node(("point", label), kind="near" if label in NEAR else "proper")
    for name, support in CHAINS.items():
        graph.add_node(("root", name), kind="exceptional")
        for label in support:
            graph.add_edge(("root", name), ("point", label))
    for name, (_coefficients, support) in LINES.items():
        graph.add_node(("root", name), kind="line")
        for label in support:
            graph.add_edge(("root", name), ("point", label))
    return graph


def canonical_root(name: str, positions: dict[str, int]) -> tuple[int, ...]:
    answer = [0] * 11
    if name in LINES:
        answer[0] = 1
        for label in LINES[name][1]:
            answer[1 + positions[label]] = -1
    else:
        near, parent = CHAINS[name]
        answer[1 + positions[near]] = -1
        answer[1 + positions[parent]] = 1
    return tuple(answer)


def branch_audit(data: dict[str, Any], target: Sequence[int]) -> dict[str, Any]:
    order = contraction_order(data, target)
    strict: list[sp.Matrix] = []
    totals: list[sp.Matrix] = []
    for index in order:
        curve = data["walls"][index]
        if strict:
            local_gram = sp.Matrix([[dot(data["gram"], left, right) for right in strict] for left in strict])
            pair = sp.Matrix([dot(data["gram"], left, curve) for left in strict])
            coefficients = local_gram.inv() * pair
            total = curve - sum((coefficients[position] * strict[position] for position in range(len(strict))), sp.zeros(11, 1))
        else:
            total = curve
        if any(sp.denom(value) != 1 for value in total):
            raise ArithmeticError("nonintegral total exceptional transform")
        strict.append(curve)
        totals.append(total)
    h = sp.Matrix(plane_class(data, target))
    basis = sp.Matrix.hstack(h, *totals)
    basis_gram = basis.T * data["gram"] * basis
    if basis_gram != sp.diag(1, *([-1] * 10)) or abs(int(basis.det())) != 1:
        raise ArithmeticError("plane/exceptional basis is not unimodular I(1,10)")
    inverse = basis.inv()
    coordinates = []
    for wall in data["walls"]:
        coordinate = inverse * wall
        if any(sp.denom(value) != 1 for value in coordinate):
            raise ArithmeticError("wall transport is not integral")
        coordinates.append(tuple(int(value) for value in coordinate))
    if basis * sp.Matrix.hstack(*(sp.Matrix(row) for row in coordinates)) != sp.Matrix.hstack(*data["walls"]):
        raise ArithmeticError("thirty-wall inverse transport failed")
    if tuple(inverse * (-data["boundary"])) != (-3, *([1] * 10)):
        raise ArithmeticError("canonical class is not -3H+sum(E)")
    if coordinates[0] != (3, *([-1] * 10)):
        raise ArithmeticError("boundary is not -K in plane coordinates")

    raw_graph = nx.Graph()
    root_rows = [(index, coordinates[index]) for index in data["roots"]]
    near_positions: set[int] = set()
    for _index, row in root_rows:
        if row[0] == 0:
            negative = [position for position, value in enumerate(row[1:]) if value == -1]
            positive = [position for position, value in enumerate(row[1:]) if value == 1]
            if len(negative) != 1 or len(positive) != 1:
                raise ArithmeticError("degree-zero root is not E_parent-E_near")
            near_positions.add(negative[0])
    if len(near_positions) != 3:
        raise ArithmeticError("branch does not contain three near centers")
    for position in range(10):
        raw_graph.add_node(("point", position), kind="near" if position in near_positions else "proper")
    for index, row in root_rows:
        kind = "exceptional" if row[0] == 0 else "line"
        if kind == "line" and not (row[0] == 1 and sum(value == -1 for value in row[1:]) == 3 and all(value in (0, -1) for value in row[1:])):
            raise ArithmeticError("degree-one root is not a three-center line")
        raw_graph.add_node(("root", index), kind=kind)
        for position, value in enumerate(row[1:]):
            if value:
                raw_graph.add_edge(("root", index), ("point", position))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        canonical_graph(),
        node_match=lambda left, right: left["kind"] == right["kind"],
    )
    maps = list(matcher.isomorphisms_iter())
    if not maps:
        raise ArithmeticError("transported roots are not the B3 colored incidence graph")
    labelings = sorted((tuple(mapping[("point", position)][1] for position in range(10)), mapping) for mapping in maps)
    labels, _mapping = labelings[0]
    positions = {label: position for position, label in enumerate(labels)}
    canonical = {name: canonical_root(name, positions) for name in (*CHAINS, *LINES)}
    if set(canonical.values()) != {row for _index, row in root_rows}:
        raise ArithmeticError("canonical B3 root classes do not equal the raw root set")
    by_class = {row: name for name, row in canonical.items()}
    return {
        "contracted_wall_set_zero_based": list(target),
        "lexicographic_contraction_order_zero_based": order,
        "plane_class_in_raw_basis": list(map(int, h)),
        "basis_determinant": int(basis.det()),
        "basis_gram": [list(map(int, basis_gram.row(index))) for index in range(11)],
        "near_positions_zero_based": sorted(near_positions),
        "colored_B3_isomorphism_count": len(maps),
        "selected_labels_in_exceptional_order": list(labels),
        "root_transport": [
            {"raw_wall_zero_based": index, "curve": by_class[row], "coordinates_H_E": list(row)}
            for index, row in root_rows
        ],
        "all_thirty_wall_coordinates_H_E": [list(row) for row in coordinates],
        "all_thirty_reconstruct_exactly": True,
        "all_twelve_roots_equal_B3_classes": True,
    }


def raw_endpoint_audit(data: dict[str, Any]) -> dict[str, Any]:
    endpoints, trail = sequential_endpoints(data)
    desired = ((3,), ((0, 3), (1, 9)), ((-2, 3), (-1, 7)))
    histogram: Counter[str] = Counter()
    selected = []
    for chosen, _inverse in endpoints.values():
        signature = endpoint_signature(data, chosen)
        key = json.dumps(signature, separators=(",", ":"))
        histogram[key] += 1
        if signature == desired:
            selected.append(chosen)
    selected.sort()
    if trail != [1, 17, 87, 329, 948, 2094, 3472, 4417, 4455, 3291, 573]:
        raise ArithmeticError("sequential contraction trail changed")
    if len(histogram) != 28 or selected != [
        (1, 5, 9, 13, 15, 17, 19, 22, 26, 28),
        (3, 7, 11, 13, 15, 17, 19, 22, 26, 28),
    ]:
        raise ArithmeticError("clean 7+3 endpoint selection changed")
    branches = [branch_audit(data, target) for target in selected]
    return {
        "state_count_by_contraction_depth_0_through_10": trail,
        "all_sequential_P2_endpoint_count": len(endpoints),
        "distinct_degree_signature_count": len(histogram),
        "selected_signature": {
            "boundary_degree": 3,
            "root_degree_histogram": {"0": 3, "1": 9},
            "contracted_norm_histogram": {"-2": 3, "-1": 7},
        },
        "selected_endpoint_count": len(selected),
        "branches": branches,
    }


def monomials(degree: int) -> list[sp.Expr]:
    return [X**a * Y**b * Z ** (degree - a - b) for a in range(degree + 1) for b in range(degree - a + 1)]


def near_substitution(label: str, polynomial: sp.Expr) -> sp.Expr:
    if label == "a":
        return sp.expand(polynomial.subs({X: 1, Y: S * T, Z: S}))
    if label == "d":
        return sp.expand(polynomial.subs({X: S, Y: 1, Z: S * T}))
    if label == "g":
        return sp.expand(polynomial.subs({X: S * T, Y: S, Z: 1}))
    raise KeyError(label)


def proper_chart(label: str) -> tuple[int, tuple[int, int], tuple[sp.Rational, ...]]:
    point = tuple(sp.Rational(value) for value in PROPER[label])
    chart = next(index for index, value in enumerate(point) if value)
    normalized = tuple(value / point[chart] for value in point)
    free = tuple(index for index in range(3) if index != chart)
    return chart, free, normalized


def condition_matrix(degree: int, multiplicity: int) -> sp.Matrix:
    basis = monomials(degree)
    variables = (X, Y, Z)
    rows = []
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _chart, free, point = proper_chart(label)
        substitution = dict(zip(variables, point))
        for total in range(multiplicity):
            for first in range(total + 1):
                second = total - first
                rows.append([sp.diff(item, variables[free[0]], first, variables[free[1]], second).subs(substitution) for item in basis])
        near = label.lower()
        if near in NEAR:
            transforms = [sp.Poly(near_substitution(near, item), S, T) for item in basis]
            for total in range(multiplicity):
                for first in range(total + 1):
                    second = total - first
                    rows.append([item.coeff_monomial(S ** (multiplicity + first) * T**second) for item in transforms])
    return sp.Matrix(rows)


def modular_rank(matrix: sp.Matrix, prime: int) -> int:
    rows = [[int(value) % prime for value in matrix.row(index)] for index in range(matrix.rows)]
    rank = 0
    for column in range(matrix.cols):
        pivot = next((index for index in range(rank, matrix.rows) if rows[index][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(value * inverse) % prime for value in rows[rank]]
        for index in range(matrix.rows):
            if index != rank and rows[index][column]:
                scale = rows[index][column]
                rows[index] = [(left - scale * right) % prime for left, right in zip(rows[index], rows[rank])]
        rank += 1
    return rank


def primitive_kernel(vector: sp.Matrix) -> list[int]:
    return list(primitive(vector))


def polynomial_text(polynomial: sp.Expr) -> str:
    return str(sp.expand(polynomial)).replace("**", "^").replace(" ", "")


def local_expansion(polynomial: sp.Expr, label: str) -> tuple[sp.Expr, int, sp.Expr, int]:
    chart, free, point = proper_chart(label)
    variables = (X, Y, Z)
    expression = sp.expand(
        polynomial.subs(
            {
                variables[chart]: 1,
                variables[free[0]]: point[free[0]] + U,
                variables[free[1]]: point[free[1]] + V,
            }
        )
    )
    terms = sp.Poly(expression, U, V).terms()
    order = min(sum(exponent) for exponent, coefficient in terms if coefficient)
    tangent = sp.expand(sum(coefficient * U ** exponent[0] * V ** exponent[1] for exponent, coefficient in terms if sum(exponent) == order))
    quadratic = sp.Poly(tangent, U, V)
    a = quadratic.coeff_monomial(U**2)
    b = quadratic.coeff_monomial(U * V)
    c = quadratic.coeff_monomial(V**2)
    return expression, order, tangent, int(b**2 - 4 * a * c)


def origin_order(polynomial: sp.Expr) -> int:
    return min(sum(exponent) for exponent, coefficient in sp.Poly(sp.expand(polynomial), S, T).terms() if coefficient)


def cluster_multiplicities(polynomial: sp.Expr) -> dict[str, int]:
    answer: dict[str, int] = {}
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _expression, order, _tangent, _discriminant = local_expansion(polynomial, label)
        answer[label] = order
        near = label.lower()
        if near in NEAR:
            strict = sp.cancel(near_substitution(near, polynomial) / S**order)
            if not strict.is_polynomial(S, T):
                raise ArithmeticError("first blowup transform is not polynomial")
            answer[near] = origin_order(strict)
    return answer


def incidence_and_linear_systems() -> tuple[dict[str, Any], sp.Expr]:
    c = sp.symbols("c")
    symbolic = {
        "A": sp.Matrix([1, 0, 0]),
        "D": sp.Matrix([0, 1, 0]),
        "G": sp.Matrix([0, 0, 1]),
        "B": sp.Matrix([1, 1, 1]),
        "C": sp.Matrix([1, 1, c]),
        "E": sp.Matrix([1 / c, 1, 1]),
        "F": sp.Matrix([c, 1, c]),
    }
    forced = {
        name: sp.factor(sp.Matrix.hstack(*(symbolic[label] for label in labels)).det())
        for name, labels in {
            "BCG": ("B", "C", "G"),
            "ABE": ("A", "B", "E"),
            "BDF": ("B", "D", "F"),
            "ACF": ("A", "C", "F"),
            "CDE": ("C", "D", "E"),
        }.items()
    }
    residual = sp.factor(sp.Matrix.hstack(symbolic["E"], symbolic["F"], symbolic["G"]).det())
    if any(forced.values()) or sp.simplify(c * residual - (1 - c**2)) != 0:
        raise ArithmeticError("B3 incidence reconstruction changed")

    line_supports: dict[str, dict[str, int]] = {}
    for name, (coefficients, support) in LINES.items():
        line = coefficients[0] * X + coefficients[1] * Y + coefficients[2] * Z
        multiplicities = cluster_multiplicities(line)
        actual = {label for label, value in multiplicities.items() if value == 1}
        if actual != set(support) or any(value not in (0, 1) for value in multiplicities.values()):
            raise ArithmeticError(f"wrong cluster support for {name}")
        line_supports[name] = multiplicities

    normalized_points = {}
    for label, point in PROPER.items():
        first = next(value for value in point if value)
        normalized_points[tuple(sp.Rational(value, first) for value in point)] = label
    proper_pairs = 0
    outside = set()
    for (_left_name, (left, _)), (_right_name, (right, _)) in itertools.combinations(LINES.items(), 2):
        cross = sp.Matrix(left).cross(sp.Matrix(right))
        first = next(value for value in cross if value)
        point = tuple(sp.Rational(value, first) for value in cross)
        if point in normalized_points:
            proper_pairs += 1
        else:
            outside.add(tuple(str(value) for value in point))
    if proper_pairs != 30 or len(outside) != 6:
        raise ArithmeticError("B3 line-pair census changed")

    cubic = condition_matrix(3, 1)
    sextic = condition_matrix(6, 2)
    kernel = sextic.nullspace()
    if cubic.shape != (10, 10) or cubic.rank() != 10 or sextic.shape != (30, 28) or sextic.rank() != 27 or len(kernel) != 1:
        raise ArithmeticError("cluster Hilbert functions changed")
    coefficients = primitive_kernel(kernel[0])
    boundary = sp.expand(sum(value * item for value, item in zip(coefficients, monomials(6))))
    expected = X**4 * Y**2 + Y**4 * Z**2 + Z**4 * X**2 - 3 * X**2 * Y**2 * Z**2
    if sp.expand(boundary - expected) != 0 or cluster_multiplicities(boundary) != {label: 2 for label in POINT_ORDER}:
        raise ArithmeticError("unique boundary sextic changed")
    return {
        "incidence_normalization": "A,D,G are the coordinate points and B=(1:1:1)",
        "five_forced_determinants": {name: str(value) for name, value in forced.items()},
        "last_determinant": str(residual),
        "solutions": {"c=1": "degenerate because C=B", "c=-1": "unique nondegenerate branch"},
        "field_warning": "incidence solution itself requires characteristic not 2",
        "proper_points": {label: list(point) for label, point in PROPER.items()},
        "near_points": {
            "a": {"parent": "A", "direction": "AG", "chart": "x=1,y=s*t,z=s"},
            "d": {"parent": "D", "direction": "AD", "chart": "y=1,x=s,z=s*t"},
            "g": {"parent": "G", "direction": "DG", "chart": "z=1,x=s*t,y=s"},
        },
        "line_cluster_multiplicities": line_supports,
        "line_pair_census": {"total": 36, "at_blown_proper_points": proper_pairs, "outside_transverse": len(outside), "outside_points": [list(point) for point in sorted(outside)]},
        "linear_systems_over_Q": {
            "cubic_matrix_shape": list(cubic.shape),
            "cubic_rank": cubic.rank(),
            "h0_minus_K": cubic.cols - cubic.rank(),
            "sextic_matrix_shape": list(sextic.shape),
            "sextic_rank": sextic.rank(),
            "h0_minus_2K": sextic.cols - sextic.rank(),
            "sextic_coefficients": coefficients,
            "unique_sextic": polynomial_text(boundary),
        },
        "linear_system_ranks_mod_3": {
            "cubic_rank": modular_rank(cubic, 3),
            "sextic_rank": modular_rank(sextic, 3),
        },
    }, boundary


def parse_markers(output: str, names: dict[str, str]) -> dict[str, int]:
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    answer: dict[str, int] = {}
    for index, line in enumerate(lines[:-1]):
        if line in names:
            answer[names[line]] = int(lines[index + 1])
    return answer


def run_singular(program: str) -> str:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for hostile algebraic checks")
    completed = subprocess.run([executable, "-q", "--no-rc"], input=program, text=True, capture_output=True, check=True)
    return completed.stdout


def singular_char0(boundary: sp.Expr) -> dict[str, Any]:
    f = polynomial_text(boundary)
    program = f'''
LIB "primdec.lib";
ring r=0,(x,y,z),dp;
poly f={f};
ideal A=y,z; ideal D=x,z; ideal G=x,y;
ideal B=x-y,y-z; ideal C=x-y,y+z; ideal E=x+y,y-z; ideal F=x+y,y+z;
ideal P=std(intersect(A,D,G,B,C,E,F));
ideal J=f,jacob(f); J=std(sat(J,maxideal(1))); ideal R=std(radical(J));
list QF=factorize(f);
print("QFC"); print(size(QF[1])); print("QEXP"); print(QF[2][2]);
print("JDIM"); print(dim(J)); print("JDEG"); print(mult(J)); print("RDEG"); print(mult(R));
print("MINP"); print(size(minAssGTZ(J))); print("JMR"); print(size(reduce(J,R)));
print("RMJ"); print(size(reduce(R,J))); print("RMP"); print(size(reduce(R,P))); print("PMR"); print(size(reduce(P,R)));
LIB "absfact.lib"; def AR=absFactorize(f); setring AR;
print("AFC"); print(size(absolute_factors[1])); print("AEXP"); print(absolute_factors[2][2]); print("ACOMP"); print(absolute_factors[4]);
exit;
'''
    names = {
        "QFC": "factor_count_with_unit_Q", "QEXP": "nonunit_exponent_Q", "JDIM": "affine_cone_dimension",
        "JDEG": "singular_scheme_degree", "RDEG": "radical_degree", "MINP": "minimal_primes",
        "JMR": "J_mod_radical", "RMJ": "radical_mod_J", "RMP": "radical_mod_points", "PMR": "points_mod_radical",
        "AFC": "absolute_factor_count_with_unit", "AEXP": "absolute_nonunit_exponent", "ACOMP": "absolute_component_count",
    }
    parsed = parse_markers(run_singular(program), names)
    expected = {
        "factor_count_with_unit_Q": 2, "nonunit_exponent_Q": 1, "affine_cone_dimension": 1,
        "singular_scheme_degree": 13, "radical_degree": 7, "minimal_primes": 7,
        "J_mod_radical": 0, "radical_mod_J": 4, "radical_mod_points": 0, "points_mod_radical": 0,
        "absolute_factor_count_with_unit": 2, "absolute_nonunit_exponent": 1, "absolute_component_count": 1,
    }
    if parsed != expected:
        raise ArithmeticError(f"characteristic-zero Singular audit changed: {parsed}")
    return {**parsed, "absolutely_irreducible": True, "support_exactly_seven_points": True}


def local_char0(boundary: sp.Expr) -> dict[str, Any]:
    records = []
    for label in ("A", "D", "G", "B", "C", "E", "F"):
        _expression, order, tangent, discriminant = local_expansion(boundary, label)
        if order != 2:
            raise ArithmeticError("boundary point is not double")
        if label in ("A", "D", "G"):
            near = label.lower()
            strict = sp.expand(sp.cancel(near_substitution(near, boundary) / S**2))
            near_order = origin_order(strict)
            tangent_terms = sp.Poly(strict, S, T).terms()
            near_tangent = sp.expand(sum(coefficient * S ** exponent[0] * T ** exponent[1] for exponent, coefficient in tangent_terms if sum(exponent) == near_order))
            q = sp.Poly(near_tangent, S, T)
            near_discriminant = int(q.coeff_monomial(S * T) ** 2 - 4 * q.coeff_monomial(S**2) * q.coeff_monomial(T**2))
            if discriminant != 0 or near_order != 2 or near_discriminant == 0:
                raise ArithmeticError("A3 resolution record changed")
            kind = "A3"
        else:
            near_order = None
            near_tangent = None
            near_discriminant = None
            if discriminant == 0:
                raise ArithmeticError("A1 tangent cone became repeated in characteristic zero")
            kind = "A1"
        records.append({
            "point": label, "multiplicity": order, "tangent": str(sp.factor(tangent)), "discriminant": discriminant,
            "kind": kind, "first_blowup_near_multiplicity": near_order,
            "first_blowup_near_tangent": str(sp.factor(near_tangent)) if near_tangent is not None else None,
            "first_blowup_near_discriminant": near_discriminant,
        })
    return {
        "records": records, "A3_count": 3, "A1_count": 4, "delta_sum": 10,
        "declared_seven_proper_plus_three_near_blowups_resolve": True,
        "strict_boundary_smooth_rational": True, "strict_boundary_square": -4,
    }


def singular_char3(boundary: sp.Expr) -> dict[str, Any]:
    f = polynomial_text(boundary)
    program = f'''
LIB "primdec.lib";
ring r=3,(x,y,z),dp;
poly f={f};
ideal A=y,z; ideal D=x,z; ideal G=x,y;
ideal B=x-y,y-z; ideal C=x-y,y+z; ideal E=x+y,y-z; ideal F=x+y,y+z;
ideal P=std(intersect(A,D,G,B,C,E,F));
ideal J=f,jacob(f); J=std(sat(J,maxideal(1))); ideal R=std(radical(J)); list FF=factorize(f);
print("FC"); print(size(FF[1])); print("EXP"); print(FF[2][2]); print("JDIM"); print(dim(J));
print("JDEG"); print(mult(J)); print("RDIM"); print(dim(R)); print("RDEG"); print(mult(R));
print("RMP"); print(size(reduce(R,P))); print("PMR"); print(size(reduce(P,R))); exit;
'''
    names = {"FC": "factor_count_with_unit_F3", "EXP": "nonunit_exponent_F3", "JDIM": "affine_cone_dimension", "JDEG": "singular_scheme_degree", "RDIM": "radical_affine_cone_dimension", "RDEG": "radical_degree", "RMP": "radical_mod_points", "PMR": "points_mod_radical"}
    parsed = parse_markers(run_singular(program), names)
    expected = {
        "factor_count_with_unit_F3": 2, "nonunit_exponent_F3": 1, "affine_cone_dimension": 1,
        "singular_scheme_degree": 21, "radical_affine_cone_dimension": 1, "radical_degree": 7,
        "radical_mod_points": 0, "points_mod_radical": 0,
    }
    if parsed != expected:
        raise ArithmeticError(f"characteristic-three Singular audit changed: {parsed}")

    repeated = {}
    for label in ("B", "C", "E", "F"):
        expression, _order, tangent, discriminant = local_expansion(boundary, label)
        tangent_mod = sp.Poly(tangent, U, V, modulus=3).as_expr()
        epsilon = -1 if label in ("B", "E") else 1
        expected_tangent = sp.Poly((U - epsilon * V) ** 2, U, V, modulus=3).as_expr()
        if sp.Poly(tangent_mod - expected_tangent, U, V, modulus=3).as_expr() != 0 or discriminant % 3 != 0:
            raise ArithmeticError("characteristic-three repeated tangent check changed")
        local_mod = sp.Poly(expression, U, V, modulus=3).as_expr()
        strict = sp.Poly(sp.cancel(local_mod.subs(V, U * (epsilon + T)) / U**2), U, T, modulus=3).as_expr()
        gradient = [int(sp.Poly(sp.diff(strict, variable), U, T, modulus=3).eval({U: 0, T: 0})) % 3 for variable in (U, T)]
        if gradient == [0, 0]:
            raise ArithmeticError("first blowup did not smooth a characteristic-three repeated tangent")
        repeated[label] = {"tangent_mod_3": str(tangent_mod), "discriminant_mod_3": 0, "first_blowup_gradient_mod_3": gradient}
    return {
        **parsed,
        "repeated_tangent_points": repeated,
        "correct_rejection": "the four Q/C A1 nodes are not A1 in characteristic 3 and the singular scheme length jumps 13 to 21",
        "explicit_nonclaim": "the reduced projective singular support is not positive-dimensional; it is still the same seven points",
        "absolute_irreducibility_over_F3bar_not_claimed": True,
    }


def root_geometry(boundary_absolutely_irreducible: bool) -> dict[str, Any]:
    positions = {label: index for index, label in enumerate(POINT_ORDER)}
    roots = [canonical_root(name, positions) for name in (*CHAINS, *LINES)]
    gram = sp.diag(1, *([-1] * 10))
    canonical = sp.Matrix([-3, *([1] * 10)])
    boundary = -2 * canonical
    if not all(dot(gram, root, root) == -2 and dot(gram, canonical, root) == 0 and dot(gram, boundary, root) == 0 for root in roots):
        raise ArithmeticError("canonical root/boundary intersections changed")
    pairings = [int(dot(gram, left, right)) for left, right in itertools.combinations(roots, 2)]
    histogram = dict(sorted(Counter(pairings).items()))
    if histogram != {0: 51, 1: 15}:
        raise ArithmeticError("root intersection histogram changed")
    if not boundary_absolutely_irreducible:
        raise ArithmeticError("boundary disjointness needs geometric irreducibility")
    return {
        "root_count": 12,
        "all_roots_square_minus_two_and_K_orthogonal": True,
        "pair_intersection_histogram": {str(key): value for key, value in histogram.items()},
        "nine_line_roots_smooth_rational": True,
        "three_exceptional_chain_roots_smooth_rational": True,
        "boundary_class_H_E": list(map(int, boundary)),
        "boundary_intersection_with_every_root": 0,
        "boundary_disjointness_argument": "the strict boundary is an irreducible degree-six curve, contains no root line or exceptional component, and distinct effective irreducible curves have nonnegative local intersections; total intersection zero forces disjointness",
    }


def cdd_rays(inequalities: sp.Matrix) -> tuple[tuple[int, ...], ...]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for hostile cone replay")
    lines = ["H-representation", "begin", f"{inequalities.rows} {inequalities.cols + 1} rational"]
    lines.extend("0 " + " ".join(str(inequalities[row, column]) for column in range(inequalities.cols)) for row in range(inequalities.rows))
    lines.append("end")
    completed = subprocess.run([executable, "--rep"], input="\n".join(lines) + "\n", text=True, capture_output=True, check=True)
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cdd emitted wrong generator representation")
    answer = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = line.split()
        if Fraction(tokens[0]) != 0:
            raise ArithmeticError("homogeneous cone emitted a vertex row")
        values = [Fraction(token) for token in tokens[1:]]
        denominator = math.lcm(*(value.denominator for value in values))
        integers = [value.numerator * (denominator // value.denominator) for value in values]
        divisor = math.gcd(*(abs(value) for value in integers))
        answer.append(tuple(value // divisor for value in integers))
    if len(set(answer)) != len(answer):
        raise ArithmeticError("cdd emitted duplicate primitive rays")
    return tuple(answer)


def cone_and_return(data: dict[str, Any]) -> dict[str, Any]:
    wall_gram: sp.Matrix = data["wall_gram"]
    pivots = tuple(wall_gram.rref()[1])
    inequalities = wall_gram[:, list(pivots)]
    rays = cdd_rays(inequalities)
    basis_gram = wall_gram.extract(pivots, pivots)
    kinds: Counter[str] = Counter()
    active_histogram: Counter[str] = Counter()
    facet_bases: list[list[list[int]]] = [[] for _ in range(30)]
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("cdd ray violates a raw wall inequality")
        square = (vector.T * basis_gram * vector)[0]
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        active = [index for index, value in enumerate(slacks) if value == 0]
        kinds[kind] += 1
        active_histogram[f"{kind}_{len(active)}"] += 1
        for facet in active:
            if len(facet_bases[facet]) == 10:
                continue
            trial = [*facet_bases[facet], list(ray)]
            if sp.Matrix(trial).rank() > len(facet_bases[facet]):
                facet_bases[facet].append(list(ray))
    digest = hashlib.sha256(json.dumps(sorted(rays), separators=(",", ":")).encode()).hexdigest()
    expected_active = {"finite_10": 1333, "ideal_10": 48, "ideal_11": 98, "ideal_12": 52, "ideal_13": 9}
    if not (
        pivots == tuple(range(11))
        and len(rays) == 1540
        and kinds == Counter({"finite": 1333, "ideal": 207})
        and dict(sorted(active_histogram.items())) == expected_active
        and all(len(rows) == 10 for rows in facet_bases)
    ):
        raise ArithmeticError("independent raw cdd cone census changed")
    raw_cone = data["case"]["exact_double_description"]
    independent_cone = load(CONE)
    face_archive = load(FACES)
    if not (
        digest == raw_cone["primitive_ray_list_sha256"]
        and raw_cone["finite_vertex_count"] == 1333
        and raw_cone["ideal_vertex_count"] == 207
        and raw_cone["hyperideal_vertex_count"] == 0
        and independent_cone["exact_double_description"]["active_facet_count_histogram"] == expected_active
        and independent_cone["decision"]["classification"] == "finite_volume_noncompact"
        and face_archive["literal_endpoint_audit"]["edge_count"] == 10030
        and face_archive["literal_endpoint_audit"]["edge_endpoint_count_histogram"] == {"2": 10030}
        and face_archive["literal_endpoint_audit"]["edge_endpoint_type_histogram"] == {"ff": 4423, "fi": 4484, "ii": 1123}
        and face_archive["literal_endpoint_audit"]["finite_vertex_count"] == 1333
        and face_archive["literal_endpoint_audit"]["ideal_vertex_count"] == 207
        and face_archive["diagram_enumeration"]["connected_affine_component_count"] == 436
    ):
        raise ArithmeticError("independent cone/Vinberg archive comparison failed")

    hostile = load(RECOGNIZER)
    presentation = next(item for item in hostile["presentations"] if item["chamber_id"] == "D16-D5-H10-30")
    if len(presentation["candidate_audits"]) != 1:
        raise ArithmeticError("P-first recognizer no longer has a unique candidate")
    candidate = presentation["candidate_audits"][0]
    cm = candidate["cm_cone"]
    cusp = candidate["cusp_regeneration"]
    saturation = candidate["integral_saturation_uniqueness"]
    if not (
        hostile["all_checks_pass"]
        and candidate["coble_incidence_pass"]
        and candidate["complete_package_regenerated_from_P"]
        and candidate["kondo_type_II_root_gram_match"]
        and cm["kind_histogram"] == {"finite": 38, "ideal": 17}
        and cusp["complete_primitive_isotropic_cusp_count"] == 17
        and cusp["sets_equal"]
        and saturation["root_span_index_in_E10"] == 2
        and saturation["unique_E10_saturation_up_to_diagram_automorphism"]
    ):
        raise ArithmeticError("independent exact same-P return gate changed")
    return {
        "direct_raw_cdd": {
            "coordinate_pivots_zero_based": list(pivots), "extreme_rays": len(rays),
            "finite": kinds["finite"], "ideal": kinds["ideal"], "hyperideal": kinds["hyperideal"],
            "active_facet_histogram": dict(sorted(active_histogram.items())),
            "all_thirty_facets_irredundant": all(len(rows) == 10 for rows in facet_bases),
            "primitive_ray_digest": digest, "finite_volume_noncompact": True,
        },
        "vinberg_crosscheck": {
            "elliptic_rank9_edges": 10030, "endpoint_histogram": {"ff": 4423, "fi": 4484, "ii": 1123},
            "every_edge_has_two_endpoints": True, "connected_affine_components": 436,
            "bad_edges": 0, "bad_parabolic_neighborhoods": 0,
        },
        "same_marked_P_return": {
            "wall_inventory": {"boundary": 1, "roots": 12, "cusp_lift_minus_one": 17, "total": 30},
            "all_seventeen_f_equal_E_plus_b_recovered_from_raw_walls": True,
            "independent_recognizer_complete_package": True,
            "root_span_index_in_E10": 2,
            "unique_E10_saturation_up_to_diagram_automorphism": True,
            "exact_same_marked_H10_30": True,
        },
    }


def target_crosscheck(raw: dict[str, Any], incidence: dict[str, Any], char0: dict[str, Any], local: dict[str, Any]) -> dict[str, Any]:
    target = load(TARGET)
    target_branches = target["raw_polytope_extraction"]["branches"]
    if not (
        target["all_checks_pass"]
        and target["raw_polytope_extraction"]["all_sequential_P2_endpoint_count"] == raw["all_sequential_P2_endpoint_count"]
        and [item["contracted_wall_set_zero_based"] for item in target_branches] == [item["contracted_wall_set_zero_based"] for item in raw["branches"]]
        and [item["plane_class_in_raw_basis"] for item in target_branches] == [item["plane_class_in_raw_basis"] for item in raw["branches"]]
        and target["linear_systems"]["anticanonical_cubics"]["kernel_dimension"] == incidence["linear_systems_over_Q"]["h0_minus_K"]
        and target["linear_systems"]["antibicanonical_sextics"]["kernel_dimension"] == incidence["linear_systems_over_Q"]["h0_minus_2K"]
        and target["absolute_factor_and_singular_scheme_audit"]["singular_scheme_degree"] == char0["singular_scheme_degree"]
        and target["local_boundary_resolution"]["A3_tacnode_count"] == local["A3_count"]
        and target["local_boundary_resolution"]["A1_node_count"] == local["A1_count"]
    ):
        raise ArithmeticError("hostile reconstruction disagrees with frozen target")
    return {
        "target_json_sha256": sha256(TARGET),
        "headline_claims_agree": True,
        "target_was_read_only_after_independent_reconstruction": True,
        "target_builder_or_helpers_imported": False,
    }


def build_audit() -> dict[str, Any]:
    data = raw_input()
    endpoints = raw_endpoint_audit(data)
    incidence, boundary = incidence_and_linear_systems()
    char0 = singular_char0(boundary)
    local = local_char0(boundary)
    char3 = singular_char3(boundary)
    roots = root_geometry(char0["absolutely_irreducible"])
    returned = cone_and_return(data)
    crosscheck = target_crosscheck(endpoints, incidence, char0, local)

    checks = {
        "raw_input_only_before_crosscheck": True,
        "raw_lattice_is_odd_unimodular_signature_1_10": data["inertia"] == (1, 10, 0),
        "unique_characteristic_boundary_and_12_plus_17_partition": True,
        "all_573_sequential_endpoints_enumerated": endpoints["all_sequential_P2_endpoint_count"] == 573,
        "exactly_two_clean_7_plus_3_branches": endpoints["selected_endpoint_count"] == 2,
        "both_plane_bases_unimodular": all(abs(branch["basis_determinant"]) == 1 for branch in endpoints["branches"]),
        "all_30_walls_transport_in_both_branches": all(branch["all_thirty_reconstruct_exactly"] for branch in endpoints["branches"]),
        "B3_incidence_rigid": incidence["solutions"]["c=-1"].startswith("unique"),
        "h0_minus_K_zero_h0_minus_2K_one": incidence["linear_systems_over_Q"]["h0_minus_K"] == 0 and incidence["linear_systems_over_Q"]["h0_minus_2K"] == 1,
        "boundary_absolutely_irreducible_over_Qbar": char0["absolutely_irreducible"],
        "singular_support_and_types_exact_over_Qbar": char0["support_exactly_seven_points"] and local["A3_count"] == 3 and local["A1_count"] == 4,
        "declared_blowups_resolve_boundary": local["declared_seven_proper_plus_three_near_blowups_resolve"],
        "root_geometry_and_boundary_disjointness": roots["root_count"] == 12,
        "direct_cdd_no_hyperideal_rays": returned["direct_raw_cdd"]["hyperideal"] == 0,
        "exact_same_marked_P_returns": returned["same_marked_P_return"]["exact_same_marked_H10_30"],
        "char3_overstatement_rejected_correctly": char3["singular_scheme_degree"] == 21 and len(char3["repeated_tangent_points"]) == 4,
        "independent_target_agreement": crosscheck["headline_claims_agree"],
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "verdict": "GREEN_WITH_CORRECTED_CHARACTERISTIC_3_BOUNDARY",
        "scope": {
            "positive": "project-new strict P-first plane realization of the known abstract D16-D5 H10/30 chamber over Q/C",
            "negative": "not a new abstract polytope, not a new dimension, not compact, and not an L4 general inverse theorem",
            "independence": "no import from the target certificate or its builder; frozen target read only at the final crosscheck",
        },
        "raw_source": {"path": str(RAW.relative_to(ROOT)), "sha256": sha256(RAW), "rank": 11, "wall_count": 30, "inertia_positive_negative_zero": list(data["inertia"])},
        "raw_endpoint_reconstruction": endpoints,
        "B3_incidence_and_linear_systems": incidence,
        "characteristic_zero_absolute_and_singular_audit": char0,
        "characteristic_zero_local_resolution": local,
        "root_curve_geometry": roots,
        "direct_cone_and_same_P_return": returned,
        "characteristic_three_scope_correction": char3,
        "target_crosscheck": crosscheck,
        "lineage": {
            "polynomial": "x^4 y^2 + y^4 z^2 + z^4 x^2 - 3 x^2 y^2 z^2",
            "historical_identification": "the cyclic Choi-Lam ternary sextic; Choi and Lam, Extremal Positive Semidefinite Forms, Math. Ann. 231 (1977), 1-18",
            "doi": "https://doi.org/10.1007/BF01360024",
            "modern_primary_formula_confirmation": "https://doi.org/10.1017/fms.2026.10221",
            "priority_statement": "the equation and abstract H10/30 chamber are known; the project-new content is the strict raw-P derivation, B3/Choi-Lam identification, and exact same-marked-P return",
        },
        "polytope_to_surface_theoretical_leverage": {
            "proved_here": "For this marked H10/30 input, the lattice-visible 7+3 contraction signature plus unimodular transport, realizable B3 incidence, interpolation and effectivity/exhaustion gates are sufficient to construct a strong Coble carrier whose complete chamber is exactly P.",
            "reusable_sufficient_gate": [
                "a primitive characteristic b^2=-1 wall splitting the other walls into b-orthogonal (-2)-roots and E with f=E+b primitive isotropic",
                "a sequential rank-ten contraction endpoint with unimodular I(1,10) plane basis",
                "a realizable nondegenerate incidence scheme for the transported root classes",
                "the required Hilbert functions, reduced/irreducible boundary, and resolution conditions",
                "a complete cone plus root/(-1)-curve exhaustion theorem returning every wall",
            ],
            "necessary_parts_for_any_carrier_of_this_marked_kind": "the lattice partition, unimodular plane marking, transported effective root incidence, Hilbert functions and complete-wall return",
            "missing_for_L4": "prove these conditions are exhaustive up to allowed birational carrier types for a variable class of arithmetic chambers, then classify the resulting incidence schemes/moduli; this audit proves neither that every compact arithmetic P has a surface carrier nor a higher-dimensional bound",
            "classification_program": "enumerate lattice contraction signatures from P, form incidence realization schemes, stratify by Hilbert/effectivity conditions, and use emptiness or AG classification for nonexistence while using rational points for new examples",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", type=Path)
    args = parser.parse_args()
    payload = build_audit()
    if args.check_archive is not None and payload != load(args.check_archive):
        raise SystemExit("live hostile audit differs from frozen JSON")
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
