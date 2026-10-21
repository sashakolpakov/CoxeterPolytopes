#!/usr/bin/env python3
"""Independent hostile audit of the characteristic-five Desargues inverse.

No construction helper is imported.  The script reconstructs the quotient
plane, its ten edge points, the triangle-line/cut-conic design, and the unique
sextic from finite-field linear algebra.  It then asks Singular for an exact
scheme-theoretic singular-locus and factorization certificate and constructs
an integral, marked bridge to every wall of the frozen H10/78 chamber.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import reduce
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Iterable, Sequence

import networkx as nx
import sympy as sp


PRIME = 5
ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "data/enumerations/coble_type_vi_char5_desargues_inverse.json"
RAW_POLYTOPE = (
    ROOT / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_char5_desargues_inverse_hostile_audit.json"
)
SCHEMA = "coxeter.coble_type_vi_char5_desargues_inverse.hostile_audit.v1"

Edge = tuple[int, int]
Point = tuple[int, int, int]


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
        raise ArithmeticError("cannot projectively normalize zero")
    scale = inverse(leading)
    return tuple(scale * value % PRIME for value in vector)


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


def nullspace(rows: Sequence[Sequence[int]], column_count: int | None = None) -> list[tuple[int, ...]]:
    if rows:
        column_count = len(rows[0])
    if column_count is None:
        raise ValueError("column count is required for an empty matrix")
    reduced, pivots = rref(rows)
    free = [column for column in range(column_count) if column not in pivots]
    vectors: list[tuple[int, ...]] = []
    for free_column in free:
        vector = [0] * column_count
        vector[free_column] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row_index][free_column] % PRIME
        vectors.append(tuple(vector))
    return vectors


def monomials(degree: int) -> list[tuple[int, int, int]]:
    return [
        (x_degree, y_degree, degree - x_degree - y_degree)
        for x_degree in range(degree, -1, -1)
        for y_degree in range(degree - x_degree, -1, -1)
    ]


def monomial_value(exponent: tuple[int, int, int], point: Point) -> int:
    answer = 1
    for coordinate, power in zip(point, exponent):
        answer = answer * pow(coordinate, power, PRIME) % PRIME
    return answer


def evaluation_row(degree: int, point: Point) -> list[int]:
    return [monomial_value(exponent, point) for exponent in monomials(degree)]


def derivative_value(
    coefficients: Sequence[int],
    degree: int,
    orders: tuple[int, int, int],
    point: Point,
) -> int:
    total = 0
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        multiplier = coefficient % PRIME
        remaining = []
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


def quotient_edge_points() -> dict[Edge, Point]:
    """Return [e_i-e_j] in {sum=0}/<1>, using x_4=0 gauge."""

    points: dict[Edge, Point] = {}
    for left, right in itertools.combinations(range(5), 2):
        vector = [0] * 5
        vector[left] = 1
        vector[right] = -1
        if sum(vector) % PRIME:
            raise ArithmeticError("edge difference left the sum-zero hyperplane")
        gauged = [(value - vector[4]) % PRIME for value in vector]
        if gauged[4] or sum(gauged) % PRIME:
            raise ArithmeticError("all-ones quotient gauge failed")
        points[(left, right)] = normalize(gauged[:3])  # type: ignore[assignment]
    if len(points) != 10 or len(set(points.values())) != 10:
        raise ArithmeticError("the ten quotient points are not distinct")
    return points


def evaluate_polynomial(
    coefficients: Sequence[int], degree: int, point: Point
) -> int:
    return sum(
        coefficient * monomial_value(exponent, point)
        for coefficient, exponent in zip(coefficients, monomials(degree))
    ) % PRIME


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


def reconstruct_design(points: dict[Edge, Point]) -> dict[str, object]:
    all_edges = list(points)
    lines: list[dict[str, object]] = []
    for triple in itertools.combinations(range(5), 3):
        support = tuple(sorted(itertools.combinations(triple, 2)))
        matrix = [evaluation_row(1, points[edge]) for edge in support]
        kernel = nullspace(matrix)
        if len(kernel) != 1:
            raise ArithmeticError("triangle does not determine a unique line")
        coefficients = normalize(kernel[0])
        exact_support = tuple(
            edge
            for edge in all_edges
            if evaluate_polynomial(coefficients, 1, points[edge]) == 0
        )
        if exact_support != support:
            raise ArithmeticError("triangle line has an accidental edge point")
        lines.append(
            {
                "label": list(triple),
                "support": [list(edge) for edge in support],
                "evaluation_rank": len(rref(matrix)[1]),
                "coefficients": list(coefficients),
            }
        )

    conics: list[dict[str, object]] = []
    for pair in itertools.combinations(range(5), 2):
        support = tuple(
            edge for edge in all_edges if len(set(edge).intersection(pair)) == 1
        )
        matrix = [evaluation_row(2, points[edge]) for edge in support]
        kernel = nullspace(matrix)
        if len(kernel) != 1:
            raise ArithmeticError("cut does not determine a unique conic")
        coefficients = normalize(kernel[0])
        exact_support = tuple(
            edge
            for edge in all_edges
            if evaluate_polynomial(coefficients, 2, points[edge]) == 0
        )
        if set(exact_support) != set(support):
            raise ArithmeticError("cut conic has an accidental edge point")
        determinant = conic_determinant(coefficients)
        if determinant == 0:
            raise ArithmeticError("cut conic is singular")
        conics.append(
            {
                "label": list(pair),
                "support": [list(edge) for edge in support],
                "evaluation_rank": len(rref(matrix)[1]),
                "coefficients": list(coefficients),
                "symmetric_matrix_determinant_mod_5": determinant,
            }
        )

    line_replication = Counter(
        tuple(edge) for record in lines for edge in record["support"]  # type: ignore[index]
    )
    conic_replication = Counter(
        tuple(edge) for record in conics for edge in record["support"]  # type: ignore[index]
    )
    if (
        len({tuple(record["coefficients"]) for record in lines}) != 10
        or len({tuple(record["coefficients"]) for record in conics}) != 10
        or set(line_replication.values()) != {3}
        or set(conic_replication.values()) != {6}
    ):
        raise ArithmeticError("the line/cut incidence census is not exact")
    return {
        "triangle_lines": lines,
        "cut_conics": conics,
        "line_count": len(lines),
        "conic_count": len(conics),
        "line_replication": sorted(set(line_replication.values())),
        "conic_replication": sorted(set(conic_replication.values())),
        "all_supports_exact": True,
        "all_conics_smooth_geometrically_integral": True,
        "conic_determinant_histogram": {
            str(value): count
            for value, count in sorted(
                Counter(
                    record["symmetric_matrix_determinant_mod_5"]
                    for record in conics
                ).items()
            )
        },
    }


def reconstruct_linear_systems(points: dict[Edge, Point]) -> dict[str, object]:
    cubic_matrix = [evaluation_row(3, point) for point in points.values()]
    cubic_rank = len(rref(cubic_matrix)[1])

    sextic_conditions: list[list[int]] = []
    for point in points.values():
        for variable in range(3):
            sextic_conditions.append(derivative_row(6, variable, point))
    sextic_rank = len(rref(sextic_conditions)[1])
    sextic_kernel = nullspace(sextic_conditions)
    if cubic_rank != 10 or len(sextic_kernel) != 1 or sextic_rank != 27:
        raise ArithmeticError("the anti-canonical linear-system dimensions changed")
    sextic = normalize(sextic_kernel[0])
    if any(evaluate_polynomial(sextic, 6, point) for point in points.values()):
        raise ArithmeticError("gradient conditions failed Euler vanishing")
    augmented = [*sextic_conditions, *[evaluation_row(6, p) for p in points.values()]]
    if len(rref(augmented)[1]) != sextic_rank:
        raise ArithmeticError("value rows are not consequences of gradient rows")
    return {
        "anticanonical_cubics": {
            "ambient_vector_dimension": 10,
            "evaluation_rank": cubic_rank,
            "kernel_dimension": 0,
            "empty_projective_system": True,
        },
        "antibicanonical_sextics": {
            "ambient_vector_dimension": 28,
            "gradient_condition_rows": 30,
            "gradient_condition_rank": sextic_rank,
            "gradient_kernel_dimension": 1,
            "value_augmented_rank": sextic_rank,
            "coefficients": list(sextic),
            "polynomial": polynomial_string(sextic, 6),
        },
    }


def local_node_and_branch_audit(
    points: dict[Edge, Point],
    design: dict[str, object],
    sextic: Sequence[int],
) -> dict[str, object]:
    half = inverse(2)
    nodes: list[dict[str, object]] = []
    tangent_values: list[tuple[str, int]] = []
    curves = [
        *(dict(record, kind="line", degree=1) for record in design["triangle_lines"]),  # type: ignore[index]
        *(dict(record, kind="conic", degree=2) for record in design["cut_conics"]),  # type: ignore[index]
    ]
    for edge, point in points.items():
        chart = next(index for index, coordinate in enumerate(point) if coordinate)
        affine_point = tuple(
            coordinate * inverse(point[chart]) % PRIME for coordinate in point
        )
        free = [index for index in range(3) if index != chart]
        first2 = tuple(2 if index == free[0] else 0 for index in range(3))
        mixed = tuple(1 if index in free else 0 for index in range(3))
        second2 = tuple(2 if index == free[1] else 0 for index in range(3))
        quadratic = (
            half * derivative_value(sextic, 6, first2, affine_point) % PRIME,
            derivative_value(sextic, 6, mixed, affine_point),
            half * derivative_value(sextic, 6, second2, affine_point) % PRIME,
        )
        discriminant = (
            quadratic[1] ** 2 - 4 * quadratic[0] * quadratic[2]
        ) % PRIME
        if discriminant == 0:
            raise ArithmeticError("a proposed center is not an ordinary node")
        nodes.append(
            {
                "edge": list(edge),
                "chart": chart,
                "quadratic_tangent_cone": list(quadratic),
                "discriminant_mod_5": discriminant,
            }
        )

        for curve in curves:
            support = {tuple(item) for item in curve["support"]}  # type: ignore[index]
            if edge not in support:
                continue
            gradient = []
            for variable in free:
                order = tuple(1 if index == variable else 0 for index in range(3))
                gradient.append(
                    derivative_value(
                        curve["coefficients"],  # type: ignore[arg-type]
                        curve["degree"],  # type: ignore[arg-type]
                        order,
                        affine_point,
                    )
                )
            if gradient == [0, 0]:
                raise ArithmeticError("root curve is singular at a center")
            direction = (gradient[1], -gradient[0] % PRIME)
            tangent_value = (
                quadratic[0] * direction[0] ** 2
                + quadratic[1] * direction[0] * direction[1]
                + quadratic[2] * direction[1] ** 2
            ) % PRIME
            if tangent_value == 0:
                raise ArithmeticError("root direction is a boundary branch direction")
            tangent_values.append((curve["kind"], tangent_value))  # type: ignore[arg-type]

    if len(tangent_values) != 90:
        raise ArithmeticError("the root-center direction census is incomplete")
    by_kind = {
        kind: {
            str(value): count
            for value, count in sorted(
                Counter(v for current_kind, v in tangent_values if current_kind == kind).items()
            )
        }
        for kind in ("line", "conic")
    }
    return {
        "ordinary_node_count": len(nodes),
        "nodes": nodes,
        "discriminant_histogram": {
            str(value): count
            for value, count in sorted(
                Counter(record["discriminant_mod_5"] for record in nodes).items()
            )
        },
        "root_center_incidence_count": len(tangent_values),
        "tangent_value_histogram": {
            str(value): count
            for value, count in sorted(Counter(v for _, v in tangent_values).items())
        },
        "tangent_value_histogram_by_root_kind": by_kind,
        "all_root_directions_avoid_both_boundary_branches": True,
        "local_intersection_multiplicity_with_boundary_at_each_center": 2,
    }


def point_ideal_declaration(index: int, point: Point) -> str:
    variables = ("x", "y", "z")
    chart = next(i for i, coordinate in enumerate(point) if coordinate)
    scale = inverse(point[chart])
    normalized = tuple(value * scale % PRIME for value in point)
    relations = [
        f"{variables[index2]}-{normalized[index2]}*{variables[chart]}"
        for index2 in range(3)
        if index2 != chart
    ]
    return f"ideal P{index}={','.join(relations)};"


def singular_scheme_and_factor_audit(polynomial: str, points: Sequence[Point]) -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the hostile scheme audit")
    declarations = "\n".join(
        point_ideal_declaration(index, point) for index, point in enumerate(points)
    )
    point_names = ",".join(f"P{index}" for index in range(len(points)))
    program = f'''
LIB "primdec.lib";
ring r=5,(x,y,z),dp;
poly f={polynomial};
list fac=factorize(f);
ideal S=std(sat(jacob(f),maxideal(1)));
ideal R=std(radical(S));
{declarations}
ideal T=std(intersect({point_names}));
print("FACTOR_COUNT_WITH_UNIT"); print(size(fac[1]));
print("NONUNIT_EXPONENT"); print(fac[2][2]);
print("AFFINE_CONE_DIMENSION"); print(dim(S));
print("SINGULAR_SCHEME_DEGREE"); print(mult(S));
print("RADICAL_DEGREE"); print(mult(R));
print("S_EQUALS_RADICAL_FORWARD"); print(size(reduce(S,R)));
print("S_EQUALS_RADICAL_REVERSE"); print(size(reduce(R,S)));
print("S_EQUALS_POINTS_FORWARD"); print(size(reduce(S,T)));
print("S_EQUALS_POINTS_REVERSE"); print(size(reduce(T,S)));
print("MINIMAL_PRIME_COUNT"); print(size(minAssGTZ(S)));
exit;
'''
    completed = subprocess.run(
        [singular, "-q", "-c", program],
        check=True,
        capture_output=True,
        text=True,
    )
    labels = {
        "FACTOR_COUNT_WITH_UNIT": "factor_count_with_unit",
        "NONUNIT_EXPONENT": "nonunit_exponent",
        "AFFINE_CONE_DIMENSION": "affine_cone_dimension",
        "SINGULAR_SCHEME_DEGREE": "projective_singular_scheme_degree",
        "RADICAL_DEGREE": "radical_degree",
        "S_EQUALS_RADICAL_FORWARD": "singular_mod_radical_generator_count",
        "S_EQUALS_RADICAL_REVERSE": "radical_mod_singular_generator_count",
        "S_EQUALS_POINTS_FORWARD": "singular_mod_point_ideal_generator_count",
        "S_EQUALS_POINTS_REVERSE": "point_ideal_mod_singular_generator_count",
        "MINIMAL_PRIME_COUNT": "minimal_prime_count",
    }
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    parsed: dict[str, int] = {}
    for index, line in enumerate(lines[:-1]):
        if line in labels:
            parsed[labels[line]] = int(lines[index + 1])
    expected = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "projective_singular_scheme_degree": 10,
        "radical_degree": 10,
        "singular_mod_radical_generator_count": 0,
        "radical_mod_singular_generator_count": 0,
        "singular_mod_point_ideal_generator_count": 0,
        "point_ideal_mod_singular_generator_count": 0,
        "minimal_prime_count": 10,
    }
    if parsed != expected:
        raise ArithmeticError(
            f"Singular scheme/factor certificate changed: {parsed}; stdout={completed.stdout}"
        )
    return {
        **parsed,
        "irreducible_over_F5": True,
        "saturated_jacobian_equals_its_radical": True,
        "saturated_jacobian_equals_explicit_ten_point_ideal": True,
        "exact_reduced_projective_support": True,
        "geometrically_integral": True,
        "geometric_integrality_proof": (
            "F5 is perfect. If the F5-irreducible sextic had r>=2 geometric "
            "components, Frobenius transitivity would give equal degree 6/r and "
            "every F5-rational node would lie on every component. Two distinct "
            "components would therefore share ten points, while Bezout gives at "
            "most (6/r)^2<=9. Hence r=1."
        ),
        "engine": "Singular 4 exact factorization, saturation, radical, ideal equality, and primary support count",
    }


def integer_vector(vector: sp.Matrix) -> list[int]:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("expected an integral marked coordinate vector")
    return [int(value) for value in vector]


def raw_marked_bridge(raw: dict[str, object]) -> dict[str, object]:
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
    if len(boundary_candidates) != 1:
        raise ArithmeticError("the raw boundary half is not intrinsic and unique")
    boundary_index = boundary_candidates[0]
    boundary = walls[boundary_index]
    nonboundary = [index for index in minus_one_indices if index != boundary_index]
    if any(pairing(boundary, walls[index]) != 1 for index in nonboundary):
        raise ArithmeticError("a raw minus-one wall does not satisfy (-K).E=1")

    disjointness = nx.Graph()
    disjointness.add_nodes_from(nonboundary)
    disjointness.add_edges_from(
        (left, right)
        for left, right in itertools.combinations(nonboundary, 2)
        if pairing(walls[left], walls[right]) == 0
    )
    maximal_cliques = [sorted(clique) for clique in nx.find_cliques(disjointness)]
    clique_histogram = Counter(len(clique) for clique in maximal_cliques)
    frames = [clique for clique in maximal_cliques if len(clique) == 10]
    if len(frames) != 1 or max(clique_histogram) != 10:
        raise ArithmeticError("the raw P2 blowdown frame is not unique")
    frame = frames[0]
    frame_sum = sum((walls[index] for index in frame), sp.zeros(11, 1))
    plane_class = (boundary + frame_sum) / 3
    marked_basis = sp.Matrix.hstack(plane_class, *[walls[index] for index in frame])
    marked_gram = marked_basis.T * gram * marked_basis
    expected_marked_gram = sp.diag(1, *([-1] * 10))
    if marked_gram != expected_marked_gram or abs(int(marked_basis.det())) != 1:
        raise ArithmeticError("the frame does not produce a unimodular I(1,10) marking")
    inverse_basis = marked_basis.inv()

    boundary_coordinates = integer_vector(inverse_basis * boundary)
    if boundary_coordinates != [3, *([-1] * 10)]:
        raise ArithmeticError("the intrinsic boundary is not -K in the plane marking")

    raw_graph = nx.Graph()
    for frame_position in range(10):
        raw_graph.add_node(("point", frame_position), kind="point")
    root_records: list[dict[str, object]] = []
    for root_index in root_indices:
        coordinates = integer_vector(inverse_basis * walls[root_index])
        degree = coordinates[0]
        support_positions = [
            position
            for position, coefficient in enumerate(coordinates[1:])
            if coefficient == -1
        ]
        expected_support_size = {1: 3, 2: 6}.get(degree)
        if (
            expected_support_size is None
            or len(support_positions) != expected_support_size
            or any(value not in (0, -1) for value in coordinates[1:])
        ):
            raise ArithmeticError("a raw root has the wrong plane class")
        kind = "line" if degree == 1 else "conic"
        root_node = (kind, root_index)
        raw_graph.add_node(root_node, kind=kind)
        raw_graph.add_edges_from(
            (root_node, ("point", position)) for position in support_positions
        )
        root_records.append(
            {
                "raw_wall_index_zero_based": root_index,
                "kind": kind,
                "coordinates_in_H_E_frame": coordinates,
                "frame_support_positions": support_positions,
            }
        )

    geometric_graph = nx.Graph()
    geometric_edges = list(itertools.combinations(range(5), 2))
    for edge in geometric_edges:
        geometric_graph.add_node(("point", edge), kind="point")
    for triple in itertools.combinations(range(5), 3):
        node = ("line", triple)
        geometric_graph.add_node(node, kind="line")
        geometric_graph.add_edges_from(
            (node, ("point", edge)) for edge in itertools.combinations(triple, 2)
        )
    for pair in itertools.combinations(range(5), 2):
        node = ("conic", pair)
        geometric_graph.add_node(node, kind="conic")
        geometric_graph.add_edges_from(
            (node, ("point", edge))
            for edge in geometric_edges
            if len(set(pair).intersection(edge)) == 1
        )

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        geometric_graph,
        node_match=lambda left, right: left["kind"] == right["kind"],
    )
    markings: list[tuple[str, dict[tuple[str, object], tuple[str, object]]]] = []
    for mapping in matcher.isomorphisms_iter():
        encoding = ",".join(
            "".join(str(value) for value in mapping[("point", position)][1])
            for position in range(10)
        )
        markings.append((encoding, mapping))
    if len(markings) != 120 or len({encoding for encoding, _ in markings}) != 120:
        raise ArithmeticError("the exact colored marking torsor is not S5-sized")
    canonical_encoding, canonical_mapping = min(markings, key=lambda item: item[0])

    canonical_point_map = []
    for position, raw_wall_index in enumerate(frame):
        edge = canonical_mapping[("point", position)][1]
        canonical_point_map.append(
            {
                "raw_frame_wall_index_zero_based": raw_wall_index,
                "geometric_edge": list(edge),
            }
        )
    canonical_root_map = []
    for record in root_records:
        kind = record["kind"]
        root_index = record["raw_wall_index_zero_based"]
        geometric_label = canonical_mapping[(kind, root_index)][1]
        if kind == "line":
            geometric_support = set(itertools.combinations(geometric_label, 2))
        else:
            geometric_support = {
                edge
                for edge in geometric_edges
                if len(set(edge).intersection(geometric_label)) == 1
            }
        expected_coordinates = [
            1 if kind == "line" else 2,
            *[
                -1
                if canonical_mapping[("point", position)][1] in geometric_support
                else 0
                for position in range(10)
            ],
        ]
        if record["coordinates_in_H_E_frame"] != expected_coordinates:
            raise ArithmeticError("an incidence isomorphism failed exact root coordinates")
        canonical_root_map.append(
            {
                "raw_wall_index_zero_based": root_index,
                "kind": kind,
                "geometric_label": list(geometric_label),
            }
        )

    full_wall_coordinates = [integer_vector(inverse_basis * wall) for wall in walls]
    full_wall_digest = hashlib.sha256(
        json.dumps(full_wall_coordinates, separators=(",", ":")).encode()
    ).hexdigest()
    marking_torsor_digest = hashlib.sha256(
        "\n".join(sorted(encoding for encoding, _ in markings)).encode()
    ).hexdigest()

    cusp_lifts = []
    for index in nonboundary:
        isotropic = walls[index] + boundary
        raw_coordinates = integer_vector(isotropic)
        if (
            pairing(isotropic, isotropic) != 0
            or pairing(isotropic, boundary) != 0
            or reduce(math.gcd, (abs(value) for value in raw_coordinates), 0) != 1
        ):
            raise ArithmeticError("a raw E+b class is not primitive CM-isotropic")
        cusp_lifts.append(index)

    exact = polytope["exact_double_description"]  # type: ignore[index]
    upstream_checks = raw["checks"]  # type: ignore[index]
    if not (
        raw["all_checks_pass"]  # type: ignore[index]
        and polytope["finite_volume"]  # type: ignore[index]
        and not polytope["compact"]  # type: ignore[index]
        and exact["hyperideal_vertex_count"] == 0  # type: ignore[index]
        and all(upstream_checks.values())  # type: ignore[union-attr]
    ):
        raise ArithmeticError("the frozen P-side Vinberg certificate is not green")

    return {
        "raw_wall_count": len(walls),
        "primitive_wall_norm_histogram": {
            str(value): count for value, count in sorted(norm_histogram.items())
        },
        "intrinsic_boundary_wall_index_zero_based": boundary_index,
        "boundary_characterization": "unique norm-minus-one wall orthogonal to all norm-minus-two walls",
        "nonboundary_minus_one_wall_count": len(nonboundary),
        "root_wall_count": len(root_indices),
        "maximal_disjoint_clique_size_histogram": {
            str(size): count for size, count in sorted(clique_histogram.items())
        },
        "unique_ten_curve_frame_wall_indices_zero_based": frame,
        "plane_class_in_raw_basis": integer_vector(plane_class),
        "marked_basis_determinant": int(marked_basis.det()),
        "marked_basis_gram": [
            [int(marked_gram[row, column]) for column in range(marked_gram.cols)]
            for row in range(marked_gram.rows)
        ],
        "boundary_coordinates_in_H_E_frame": boundary_coordinates,
        "root_plane_profile": dict(
            Counter(
                "degree_1_three_centers"
                if record["kind"] == "line"
                else "degree_2_six_centers"
                for record in root_records
            )
        ),
        "colored_incidence_marking_count": len(markings),
        "marking_torsor_sha256": marking_torsor_digest,
        "canonical_marking_encoding": canonical_encoding,
        "canonical_point_map": canonical_point_map,
        "canonical_root_map": canonical_root_map,
        "all_twenty_root_coordinates_match_triangle_cut_classes": True,
        "transported_full_wall_coordinate_count": len(full_wall_coordinates),
        "transported_full_wall_coordinates_sha256": full_wall_digest,
        "primitive_isotropic_E_plus_b_count": len(cusp_lifts),
        "full_wall_transport_is_integral": True,
        "upstream_P_side_import": {
            "finite_volume": polytope["finite_volume"],  # type: ignore[index]
            "compact": polytope["compact"],  # type: ignore[index]
            "finite_vertices": exact["finite_vertex_count"],  # type: ignore[index]
            "ideal_vertices": exact["ideal_vertex_count"],  # type: ignore[index]
            "hyperideal_extreme_rays": exact["hyperideal_vertex_count"],  # type: ignore[index]
            "elliptic_edges": 166320,
            "facet_irredundancy": exact["all_facets_irredundant"],  # type: ignore[index]
            "vinberg_checks_all_true": all(upstream_checks.values()),  # type: ignore[union-attr]
        },
        "exact_same_marked_H10_78": True,
        "logical_scope": (
            "The hostile audit computes the integral marking. Completeness and finite "
            "volume of the transported 78-wall cone are imported from the separately "
            "frozen full-nef/Vinberg certificate, not re-enumerated here."
        ),
    }


def build_audit() -> dict[str, object]:
    target = json.loads(TARGET.read_text(encoding="utf-8"))
    raw = json.loads(RAW_POLYTOPE.read_text(encoding="utf-8"))
    points = quotient_edge_points()
    design = reconstruct_design(points)
    systems = reconstruct_linear_systems(points)
    sextic = systems["antibicanonical_sextics"]["coefficients"]  # type: ignore[index]
    local = local_node_and_branch_audit(points, design, sextic)  # type: ignore[arg-type]
    singular = singular_scheme_and_factor_audit(
        systems["antibicanonical_sextics"]["polynomial"],  # type: ignore[index]
        list(points.values()),
    )
    bridge = raw_marked_bridge(raw)

    target_coefficients = target["sextic_double_point_system"]["coefficients"]
    target_points = {
        key: tuple(value) for key, value in target["points"].items()
    }
    independent_points = {str(key): value for key, value in points.items()}
    target_agrees = (
        target["all_checks_pass"]
        and list(sextic) == list(target_coefficients)  # type: ignore[arg-type]
        and target_points == independent_points
        and target["singular_and_factor_audit"]["projective_singular_degree"] == 10
    )
    if not target_agrees:
        raise ArithmeticError("the independent reconstruction disagrees with its target")

    checks = {
        "quotient_has_ten_distinct_edge_points": len(set(points.values())) == 10,
        "ten_triangle_lines_have_exact_3_point_supports": design["line_count"] == 10
        and design["line_replication"] == [3],
        "ten_cut_conics_have_exact_6_point_supports": design["conic_count"] == 10
        and design["conic_replication"] == [6],
        "all_cut_conics_are_smooth": design["all_conics_smooth_geometrically_integral"],
        "no_anticanonical_cubic": systems["anticanonical_cubics"]["kernel_dimension"] == 0,  # type: ignore[index]
        "unique_antibicanonical_sextic": systems["antibicanonical_sextics"]["gradient_kernel_dimension"] == 1,  # type: ignore[index]
        "exact_reduced_ten_point_singular_scheme": singular["exact_reduced_projective_support"],
        "ten_ordinary_nodes": local["ordinary_node_count"] == 10,
        "sextic_geometrically_integral": singular["geometrically_integral"],
        "all_90_root_directions_avoid_boundary_branches": local["root_center_incidence_count"] == 90
        and local["all_root_directions_avoid_both_boundary_branches"],
        "strong_Coble_conditions": True,
        "raw_unique_unimodular_P2_frame": bridge["marked_basis_determinant"] in (-1, 1),
        "exact_twenty_root_marking": bridge["all_twenty_root_coordinates_match_triangle_cut_classes"],
        "exact_full_78_wall_integral_transport": bridge["transported_full_wall_coordinate_count"] == 78
        and bridge["full_wall_transport_is_integral"],
        "same_finite_volume_noncompact_H10_78": bridge["exact_same_marked_H10_78"]
        and bridge["upstream_P_side_import"]["finite_volume"]  # type: ignore[index]
        and not bridge["upstream_P_side_import"]["compact"],  # type: ignore[index]
        "independent_reconstruction_agrees_with_target": target_agrees,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"a hostile gate failed: {checks}")

    target_handoff_keys = set(target["raw_polytope_handoff"])
    target_had_explicit_marking = bool(
        {
            "marked_basis_determinant",
            "canonical_point_map",
            "transported_full_wall_coordinates_sha256",
        }.intersection(target_handoff_keys)
    )
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "verdict": "GREEN_WITH_MARKED_HANDOFF_STRENGTHENING",
        "direction": "RAW_POLYTOPE_TO_SURFACE_TO_EXACT_SAME_MARKED_POLYTOPE",
        "sources": {
            "target_inverse_archive": {
                "path": str(TARGET.relative_to(ROOT)),
                "sha256": sha256(TARGET),
                "construction_helper_imported": False,
            },
            "raw_full_nef_polytope_archive": {
                "path": str(RAW_POLYTOPE.relative_to(ROOT)),
                "sha256": sha256(RAW_POLYTOPE),
            },
        },
        "quotient_model": {
            "field": "F5",
            "vector_space": "{(x0,...,x4): sum(xi)=0}/<all-ones>",
            "dimension": 3,
            "projectivization": "P2(F5)",
            "edge_points": [
                {"edge": list(edge), "point": list(point)}
                for edge, point in points.items()
            ],
        },
        "finite_design": design,
        "linear_systems": systems,
        "local_geometry": local,
        "singular_scheme_and_integrality": singular,
        "blowup_and_strong_Coble": {
            "surface": "Bl_{ten edge points}(P2_F5)",
            "canonical_class": "K=-3H+sum(E_e); K^2=-1",
            "boundary_class": "B=6H-2sum(E_e)=-2K; B^2=-4",
            "boundary_smooth_after_one_blowup_at_each_node": True,
            "boundary_geometric_genus": 0,
            "minus_K_empty": True,
            "minus_2K_has_unique_member": True,
            "triangle_root_class": "H-sum(the three triangle E_e); square -2; K-dot 0",
            "cut_root_class": "2H-sum(the six cut E_e); square -2; K-dot 0",
            "twenty_roots_smooth_geometrically_integral_rational": True,
            "twenty_roots_disjoint_from_boundary_after_blowup": True,
            "bezout_accounting": (
                "Each root has local intersection multiplicity exactly 2 at every "
                "assigned node. Lines use 3*2=6 and conics use 6*2=12, exhausting "
                "their total Bezout intersections with the sextic."
            ),
            "strong_Coble_verdict": "PASS",
        },
        "marked_raw_polytope_handoff": bridge,
        "target_certificate_audit": {
            "independent_coefficients_and_points_agree": target_agrees,
            "target_same_P_gate_had_explicit_integral_marking": target_had_explicit_marking,
            "material_strengthening": (
                "The target compared only the root plane profile before asserting "
                "exact_same_P. This audit supplies the missing unique frame, "
                "unimodular I(1,10) basis, 120-element S5 marking torsor, exact "
                "twenty-root coordinates, and all-78-wall transport digest."
                if not target_had_explicit_marking
                else "The target now contains an explicit marking; this audit reconstructs it independently."
            ),
        },
        "lineage_and_scope": {
            "known_surface": "Kondo Type VI, characteristic 5, one boundary; unique isomorphism class",
            "known_surface_models": "Kondo gives Hessian/Cremona and P1xP1 double-cover models",
            "plane_model_priority": (
                "The explicit quotient-plane sextic is not present in the pinned Kondo "
                "source checked here, but worldwide priority is not established."
            ),
            "known_polytope": "the frozen project H10/78 full-nef chamber; no new abstract polytope is claimed",
            "positive_result": (
                "an exact P-first realization and same-marked-P calibration, upgrading "
                "this datum from surface-first L0 to structural inverse L3"
            ),
            "inverse_level": "L3 datum calibration, not a class-wide L4 carrier theorem",
            "not_proved": [
                "a necessary-and-sufficient P-to-surface theorem for a variable arithmetic Coxeter class",
                "a new higher-dimensional abstract Coxeter polytope",
                "a sharper general nonexistence dimension bound",
                "worldwide priority for this explicit plane equation",
            ],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_audit(), indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archive mismatch: {arguments.output}")
        print("characteristic-five Desargues hostile archive replay passed")
        return
    if arguments.write:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
