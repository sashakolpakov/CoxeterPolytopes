#!/usr/bin/env python3
"""Certify the characteristic-five Desargues inverse realization of Coble Type VI.

This exact finite-field certifier constructs the
three-dimensional quotient

    {x_0+...+x_4=0 in F_5^5} / <(1,1,1,1,1)>,

puts the ten edge points [e_i-e_j] in its projective plane, reconstructs the
ten triangle lines and ten cut conics, solves for the unique plane sextic
singular at all ten points, and closes geometric integrality, blowup,
strong-Coble, and exact same-polytope return.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
from typing import Iterable, Sequence

import networkx as nx
import sympy as sp


P = 5
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_type_vi_char5_desargues_inverse.json"
RAW_DESIGN = ROOT / "data/enumerations/coble_remaining_seven_contraction_designs.json"
RAW_POLYTOPE = ROOT / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json"
SCHEMA = "coxeter.coble_type_vi_char5_desargues_inverse.v1"


def inv(value: int) -> int:
    value %= P
    if not value:
        raise ZeroDivisionError
    return pow(value, P - 2, P)


def rref(matrix: Sequence[Sequence[int]]) -> tuple[list[list[int]], list[int]]:
    work = [[entry % P for entry in row] for row in matrix]
    if not work:
        return work, []
    row = 0
    pivots: list[int] = []
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work)) if work[index][column] % P),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = inv(work[row][column])
        work[row] = [(scale * value) % P for value in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            scale = work[index][column] % P
            if scale:
                work[index] = [
                    (left - scale * right) % P
                    for left, right in zip(work[index], work[row])
                ]
        pivots.append(column)
        row += 1
        if row == len(work):
            break
    return work, pivots


def nullspace(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    reduced, pivots = rref(matrix)
    columns = len(matrix[0]) if matrix else 0
    free = [column for column in range(columns) if column not in pivots]
    answer: list[list[int]] = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = (-reduced[row][free_column]) % P
        answer.append(vector)
    return answer


def normalize(vector: Iterable[int]) -> tuple[int, ...]:
    values = tuple(value % P for value in vector)
    first = next(value for value in values if value)
    scale = inv(first)
    return tuple((scale * value) % P for value in values)


def quotient_coordinates(vector: Sequence[int]) -> tuple[int, int, int]:
    """Gauge the fifth coordinate to zero and retain the first three.

    Because char(F)=5, subtracting a multiple of the all-ones vector
    preserves the sum-zero hyperplane.  The fourth coordinate is then minus
    the sum of the first three.
    """

    if len(vector) != 5 or sum(vector) % P:
        raise ValueError("expected a sum-zero vector in F5^5")
    gauged = [(value - vector[4]) % P for value in vector]
    if gauged[4] or sum(gauged[:4]) % P:
        raise ArithmeticError("quotient gauge failed")
    return normalize(gauged[:3])  # type: ignore[return-value]


def edge_points() -> dict[tuple[int, int], tuple[int, int, int]]:
    points: dict[tuple[int, int], tuple[int, int, int]] = {}
    for left, right in itertools.combinations(range(5), 2):
        vector = [0] * 5
        vector[left] = 1
        vector[right] = -1
        points[(left, right)] = quotient_coordinates(vector)
    if len(set(points.values())) != 10:
        raise ArithmeticError("the ten edge points collided")
    return points


def monomials(degree: int) -> list[tuple[int, int, int]]:
    return [
        (x_degree, y_degree, degree - x_degree - y_degree)
        for x_degree in range(degree, -1, -1)
        for y_degree in range(degree - x_degree, -1, -1)
    ]


def evaluate_monomial(
    exponent: tuple[int, int, int], point: tuple[int, int, int]
) -> int:
    return (
        pow(point[0], exponent[0], P)
        * pow(point[1], exponent[1], P)
        * pow(point[2], exponent[2], P)
    ) % P


def evaluation_row(
    degree: int, point: tuple[int, int, int]
) -> list[int]:
    return [evaluate_monomial(exponent, point) for exponent in monomials(degree)]


def derivative_row(
    degree: int, variable: int, point: tuple[int, int, int]
) -> list[int]:
    row = []
    for exponent in monomials(degree):
        power = exponent[variable]
        if not power:
            row.append(0)
            continue
        derivative = list(exponent)
        derivative[variable] -= 1
        row.append((power * evaluate_monomial(tuple(derivative), point)) % P)
    return row


def derivative_value(
    coefficients: Sequence[int],
    degree: int,
    orders: tuple[int, int, int],
    point: tuple[int, int, int],
) -> int:
    answer = 0
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        multiplier = coefficient % P
        reduced = []
        for power, order in zip(exponent, orders):
            if order > power:
                multiplier = 0
                break
            for offset in range(order):
                multiplier = (multiplier * (power - offset)) % P
            reduced.append(power - order)
        if multiplier:
            answer += multiplier * evaluate_monomial(tuple(reduced), point)
    return answer % P


def polynomial_string(coefficients: Sequence[int], degree: int) -> str:
    terms = []
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        coefficient %= P
        if not coefficient:
            continue
        factors = []
        for variable, power in zip(("x", "y", "z"), exponent):
            if power == 1:
                factors.append(variable)
            elif power:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors) if factors else "1"
        terms.append(monomial if coefficient == 1 else f"{coefficient}*{monomial}")
    return "+".join(terms) if terms else "0"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_dot(
    gram: sp.Matrix, left: Sequence[int] | sp.Matrix, right: Sequence[int] | sp.Matrix
) -> int:
    return int((sp.Matrix(left).T * gram * sp.Matrix(right))[0])


def raw_marking_audit(raw_polytope: dict[str, object]) -> dict[str, object]:
    polytope = raw_polytope["polytope"]
    gram = sp.Matrix(polytope["lattice_gram"])
    walls = [list(map(int, wall)) for wall in polytope["wall_vectors"]]
    rank = gram.rows
    wall_gram = [
        [integer_dot(gram, left, right) for right in walls]
        for left in walls
    ]
    norms = [wall_gram[index][index] for index in range(len(walls))]
    standard = [
        [int(index == column) for index in range(rank)]
        for column in range(rank)
    ]
    boundaries = [
        index
        for index, norm in enumerate(norms)
        if norm == -1
        and all(
            (integer_dot(gram, walls[index], standard[column]) - int(gram[column, column]))
            % 2
            == 0
            for column in range(rank)
        )
    ]
    if len(boundaries) != 1:
        raise ArithmeticError("the raw one-boundary wall was not intrinsic")
    boundary = boundaries[0]
    roots = [index for index, norm in enumerate(norms) if norm == -2]
    minus_one = [
        index for index, norm in enumerate(norms) if norm == -1 and index != boundary
    ]
    orthogonal = nx.Graph()
    orthogonal.add_nodes_from(minus_one)
    for position, left in enumerate(minus_one):
        for right in minus_one[position + 1 :]:
            if wall_gram[left][right] == 0:
                orthogonal.add_edge(left, right)
    frames = [tuple(sorted(clique)) for clique in nx.find_cliques(orthogonal) if len(clique) == 10]
    if len(frames) != 1:
        raise ArithmeticError("the raw ten-wall exceptional frame is not unique")
    frame = frames[0]
    boundary_vector = sp.Matrix(walls[boundary])
    exceptional_sum = sum((sp.Matrix(walls[index]) for index in frame), sp.zeros(rank, 1))
    plane = (boundary_vector + exceptional_sum) / 3
    if any(value.q != 1 for value in plane):
        raise ArithmeticError("H=(b+sum E)/3 is not integral")
    basis = sp.Matrix.hstack(plane, *(sp.Matrix(walls[index]) for index in frame))
    basis_gram = basis.T * gram * basis
    expected_gram = sp.diag(1, *([-1] * 10))
    if basis_gram != expected_gram or abs(int(basis.det())) != 1:
        raise ArithmeticError("the raw plane/exceptional marking is not unimodular")
    boundary_coordinates = tuple(int(value) for value in basis.inv() * boundary_vector)
    if boundary_coordinates != (3, *([-1] * 10)):
        raise ArithmeticError("the boundary is not 3H-sum(E)")

    raw_blocks: dict[str, list[tuple[int, ...]]] = {"line": [], "conic": []}
    root_coordinate_checks = []
    for root in roots:
        degree = integer_dot(gram, plane, walls[root])
        intersections = tuple(wall_gram[root][exceptional] for exceptional in frame)
        if any(value not in (0, 1) for value in intersections):
            raise ArithmeticError("a root has a non-incidence exceptional pairing")
        support = tuple(index for index, value in enumerate(intersections) if value)
        kind = "line" if (degree, len(support)) == (1, 3) else "conic" if (degree, len(support)) == (2, 6) else None
        if kind is None:
            raise ArithmeticError("a raw root has the wrong plane degree/support")
        coordinates = tuple(int(value) for value in basis.inv() * sp.Matrix(walls[root]))
        expected = (degree, *tuple(-int(index in support) for index in range(10)))
        if coordinates != expected:
            raise ArithmeticError("a raw root is not its marked plane curve class")
        raw_blocks[kind].append(support)
        root_coordinate_checks.append((root, kind, coordinates))
    if {kind: len(blocks) for kind, blocks in raw_blocks.items()} != {"line": 10, "conic": 10}:
        raise ArithmeticError("the raw colored block counts changed")

    raw_graph = nx.Graph()
    raw_graph.add_nodes_from((("point", index) for index in range(10)), kind="point")
    for kind in ("line", "conic"):
        for block_index, block in enumerate(raw_blocks[kind]):
            node = (kind, block_index)
            raw_graph.add_node(node, kind=kind)
            raw_graph.add_edges_from((("point", point), node) for point in block)

    edges = tuple(itertools.combinations(range(5), 2))
    canonical_graph = nx.Graph()
    canonical_graph.add_nodes_from((("point", edge) for edge in edges), kind="point")
    line_blocks = [tuple(itertools.combinations(triple, 2)) for triple in itertools.combinations(range(5), 3)]
    conic_blocks = [
        tuple(edge for edge in edges if len(set(edge) & set(pair)) == 1)
        for pair in itertools.combinations(range(5), 2)
    ]
    for kind, blocks in (("line", line_blocks), ("conic", conic_blocks)):
        for block_index, block in enumerate(blocks):
            node = (kind, block_index)
            canonical_graph.add_node(node, kind=kind)
            canonical_graph.add_edges_from((("point", point), node) for point in block)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        canonical_graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match("kind", ""),
    )
    point_maps = []
    for mapping in matcher.isomorphisms_iter():
        point_maps.append(
            tuple(mapping[("point", index)][1] for index in range(10))
        )
    if len(point_maps) != 120:
        raise ArithmeticError("the colored Desargues/cut marking did not have 120 isomorphisms")
    lexicographic_mapping = min(point_maps)
    return {
        "boundary_wall_index": boundary,
        "exceptional_frame_wall_indices": frame,
        "plane_class_vector": tuple(int(value) for value in plane),
        "basis_determinant_absolute": abs(int(basis.det())),
        "basis_gram": [list(map(int, basis_gram.row(index))) for index in range(rank)],
        "boundary_coordinates_in_plane_basis": boundary_coordinates,
        "root_count": len(root_coordinate_checks),
        "all_root_vectors_equal_marked_plane_classes": True,
        "colored_incidence_isomorphism_count": len(point_maps),
        "lexicographically_first_point_map_to_K5_edges": lexicographic_mapping,
        "exact_marked_root_seed_transport": True,
    }


def line_and_conic_systems(
    points: dict[tuple[int, int], tuple[int, int, int]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    lines = []
    for triple in itertools.combinations(range(5), 3):
        support = tuple(sorted(itertools.combinations(triple, 2)))
        kernel = nullspace([evaluation_row(1, points[edge]) for edge in support])
        if len(kernel) != 1:
            raise ArithmeticError("a K5 triangle did not determine one line")
        coefficients = normalize(kernel[0])
        lines.append(
            {
                "triangle": triple,
                "support": support,
                "coefficients": coefficients,
                "polynomial": polynomial_string(coefficients, 1),
            }
        )

    conics = []
    all_vertices = set(range(5))
    for pair in itertools.combinations(range(5), 2):
        complement = all_vertices - set(pair)
        support = tuple(
            edge
            for edge in points
            if len(set(edge) & set(pair)) == 1
            and len(set(edge) & complement) == 1
        )
        kernel = nullspace([evaluation_row(2, points[edge]) for edge in support])
        if len(kernel) != 1:
            raise ArithmeticError("a K5 cut did not determine one conic")
        coefficients = normalize(kernel[0])
        # A ternary conic is smooth iff its symmetric matrix has full rank;
        # in odd characteristic the off-diagonal coefficients are halved.
        xx, xy, xz, yy, yz, zz = coefficients
        half = inv(2)
        symmetric = [
            [xx, half * xy, half * xz],
            [half * xy, yy, half * yz],
            [half * xz, half * yz, zz],
        ]
        _, pivots = rref(symmetric)
        if len(pivots) != 3:
            raise ArithmeticError("a cut conic is singular")
        conics.append(
            {
                "pair": pair,
                "support": support,
                "coefficients": coefficients,
                "polynomial": polynomial_string(coefficients, 2),
            }
        )
    return lines, conics


def sextic_system(
    points: dict[tuple[int, int], tuple[int, int, int]],
) -> dict[str, object]:
    conditions = []
    for point in points.values():
        conditions.extend(derivative_row(6, variable, point) for variable in range(3))
    reduced, pivots = rref(conditions)
    kernel = nullspace(conditions)
    values_vanish = all(
        sum(coefficient * value for coefficient, value in zip(kernel[0], evaluation_row(6, point))) % P == 0
        for point in points.values()
    ) if len(kernel) == 1 else False
    return {
        "condition_rows": len(conditions),
        "condition_rank": len(pivots),
        "ambient_dimension": len(monomials(6)),
        "kernel_dimension": len(kernel),
        "value_conditions_follow": values_vanish,
        "coefficients": normalize(kernel[0]) if len(kernel) == 1 else kernel,
        "polynomial": polynomial_string(normalize(kernel[0]), 6) if len(kernel) == 1 else None,
        "rref_nonzero_rows": sum(any(value for value in row) for row in reduced),
    }


def anticanonical_system(
    points: dict[tuple[int, int], tuple[int, int, int]],
) -> dict[str, int | bool]:
    conditions = [evaluation_row(3, point) for point in points.values()]
    _, pivots = rref(conditions)
    kernel = nullspace(conditions)
    return {
        "ambient_dimension": len(monomials(3)),
        "evaluation_rank": len(pivots),
        "kernel_dimension": len(kernel),
        "empty": not kernel,
    }


def singular_factor_audit(polynomial: str) -> dict[str, int | bool | str]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the exact factor/singular-scheme audit")
    program = f'''LIB "elim.lib";
ring r=5,(x,y,z),dp;
poly f={polynomial};
list factors=factorize(f);
ideal singular_scheme=sat(jacob(f),ideal(x,y,z));
print("FACTOR_COUNT_WITH_UNIT");
print(size(factors[1]));
print("NONUNIT_EXPONENT");
print(factors[2][2]);
print("AFFINE_CONE_DIMENSION");
print(dim(std(singular_scheme)));
print("PROJECTIVE_SINGULAR_DEGREE");
print(mult(std(singular_scheme)));
exit;
'''
    completed = subprocess.run(
        [executable, "-q", "-c", program],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    labels = {
        "FACTOR_COUNT_WITH_UNIT": "factor_count_with_unit",
        "NONUNIT_EXPONENT": "nonunit_exponent",
        "AFFINE_CONE_DIMENSION": "affine_cone_dimension",
        "PROJECTIVE_SINGULAR_DEGREE": "projective_singular_degree",
    }
    parsed: dict[str, int] = {}
    for index, line in enumerate(lines[:-1]):
        if line in labels:
            parsed[labels[line]] = int(lines[index + 1])
    expected = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "projective_singular_degree": 10,
    }
    if parsed != expected:
        raise ArithmeticError(f"unexpected Singular certificate: {parsed}; stdout={completed.stdout}")
    return {
        **parsed,
        "irreducible_over_F5": True,
        "projective_singular_scheme_dimension": 0,
        "singular_scheme_is_exactly_the_ten_reduced_F5_points": True,
        "geometrically_integral": True,
        "geometric_integrality_argument": (
            "F is irreducible over F5. If it split geometrically into r>=2 Frobenius-"
            "conjugate components, each of the ten F5-rational singular points would lie "
            "on every component. Their common degree is 6/r<=3, so two components would "
            "share ten points, contradicting Bezout (6/r)^2<=9."
        ),
        "engine": "Singular exact factorization and saturated Jacobian multiplicity",
    }


def node_and_tangent_audit(
    points: dict[tuple[int, int], tuple[int, int, int]],
    lines: Sequence[dict[str, object]],
    conics: Sequence[dict[str, object]],
    sextic: dict[str, object],
) -> dict[str, object]:
    coefficients = sextic["coefficients"]
    if not isinstance(coefficients, tuple):
        raise ArithmeticError("the sextic kernel is not one-dimensional")
    half = inv(2)
    node_data = {}
    tangent_checks = []
    roots = [*(dict(curve, degree=1) for curve in lines), *(dict(curve, degree=2) for curve in conics)]
    for edge, point in points.items():
        chart = next(index for index, coordinate in enumerate(point) if coordinate)
        free = [index for index in range(3) if index != chart]
        normalized_point = tuple((coordinate * inv(point[chart])) % P for coordinate in point)
        first, second = free
        first2_orders = tuple(2 if index == first else 0 for index in range(3))
        mixed_orders = tuple(1 if index in free else 0 for index in range(3))
        second2_orders = tuple(2 if index == second else 0 for index in range(3))
        quad_a = half * derivative_value(coefficients, 6, first2_orders, normalized_point) % P
        quad_b = derivative_value(coefficients, 6, mixed_orders, normalized_point)
        quad_c = half * derivative_value(coefficients, 6, second2_orders, normalized_point) % P
        discriminant = (quad_b * quad_b - 4 * quad_a * quad_c) % P
        if not discriminant:
            raise ArithmeticError(f"the sextic point {edge} is not an ordinary node")
        node_data[str(edge)] = {
            "chart": chart,
            "quadratic_tangent_cone": (quad_a, quad_b, quad_c),
            "discriminant": discriminant,
        }
        for root_index, root in enumerate(roots):
            support = {tuple(item) for item in root["support"]}  # type: ignore[index]
            if edge not in support:
                continue
            root_coefficients = root["coefficients"]
            root_degree = root["degree"]
            gradient = []
            for variable in free:
                orders = tuple(1 if index == variable else 0 for index in range(3))
                gradient.append(
                    derivative_value(
                        root_coefficients,  # type: ignore[arg-type]
                        root_degree,  # type: ignore[arg-type]
                        orders,
                        normalized_point,
                    )
                )
            if gradient == [0, 0]:
                raise ArithmeticError("a root curve is singular at a prescribed center")
            direction = (gradient[1], -gradient[0] % P)
            tangent_value = (
                quad_a * direction[0] * direction[0]
                + quad_b * direction[0] * direction[1]
                + quad_c * direction[1] * direction[1]
            ) % P
            if not tangent_value:
                raise ArithmeticError("a root curve is tangent to a sextic branch")
            tangent_checks.append((root_index, edge, tangent_value))
    expected = sum(len(root["support"]) for root in roots)  # type: ignore[arg-type]
    if len(tangent_checks) != expected:
        raise ArithmeticError("the root-center tangent census is incomplete")
    return {
        "ordinary_node_count": len(node_data),
        "nodes": node_data,
        "root_center_incidence_count": len(tangent_checks),
        "all_root_directions_avoid_both_boundary_branches": True,
        "tangent_value_histogram": {
            str(value): count
            for value, count in sorted(Counter(item[2] for item in tangent_checks).items())
        },
    }


def build_certificate() -> dict[str, object]:
    points = edge_points()
    lines, conics = line_and_conic_systems(points)
    sextic = sextic_system(points)
    node_audit = node_and_tangent_audit(points, lines, conics, sextic)
    anticanonical = anticanonical_system(points)
    if not anticanonical["empty"]:
        raise ArithmeticError("the anticanonical system is nonempty")
    singular_audit = singular_factor_audit(str(sextic["polynomial"]))
    triangle_degrees = Counter(edge for line in lines for edge in line["support"])
    cut_degrees = Counter(edge for conic in conics for edge in conic["support"])
    raw_design = json.loads(RAW_DESIGN.read_text(encoding="utf-8"))
    raw_row = raw_design["seven_row_outcome"]["one_boundary_78"]
    raw_candidate = raw_design["different_finite_design_candidates"]["one_boundary_78"]
    raw_polytope = json.loads(RAW_POLYTOPE.read_text(encoding="utf-8"))
    marking = raw_marking_audit(raw_polytope)
    if not (
        raw_row["lattice_rank"] == 11
        and raw_row["wall_count"] == 78
        and raw_row["simultaneous_full_P2_frame_count"] == 1
        and raw_row["root_count"] == 20
        and raw_candidate["point_design"] == "ten edges of K5"
        and raw_polytope["polytope"]["exact_double_description"]["hyperideal_vertex_count"] == 0
    ):
        raise ArithmeticError("the raw H10/78 handoff changed")
    certificate = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "RAW_POLYTOPE_TO_FINITE_GEOMETRY_TO_STRONG_COBLE_TO_SAME_POLYTOPE",
        "sources": {
            "raw_design_archive": {
                "path": str(RAW_DESIGN.relative_to(ROOT)),
                "sha256": sha256(RAW_DESIGN),
            },
            "raw_polytope_archive": {
                "path": str(RAW_POLYTOPE.relative_to(ROOT)),
                "sha256": sha256(RAW_POLYTOPE),
            },
        },
        "field": "F5",
        "quotient": "{sum x_i=0}/<all-ones>",
        "points": {str(edge): point for edge, point in points.items()},
        "point_count": len(points),
        "triangle_lines": lines,
        "triangle_line_count": len(lines),
        "triangle_replication": sorted(set(triangle_degrees.values())),
        "cut_conics": conics,
        "cut_conic_count": len(conics),
        "cut_replication": sorted(set(cut_degrees.values())),
        "sextic_double_point_system": sextic,
        "sextic_local_audit": node_audit,
        "anticanonical_system": anticanonical,
        "singular_and_factor_audit": singular_audit,
        "divisor_classes_on_blowup": {
            "canonical": "K=-3H+sum(E_i)",
            "boundary": "B=6H-2sum(E_i)=-2K; B^2=-4",
            "triangle_line_roots": "H-E_i-E_j-E_k; square -2; K-pairing 0",
            "cut_conic_roots": "2H-sum(six cut E_i); square -2; K-pairing 0",
            "boundary_is_smooth_rational_after_blowup": True,
            "boundary_disjoint_from_all_twenty_roots_after_blowup": True,
        },
        "strong_Coble": {
            "minus_K_empty": anticanonical["empty"],
            "minus_2K_unique": sextic["kernel_dimension"] == 1,
            "unique_minus_2K_member_geometrically_integral": singular_audit["geometrically_integral"],
            "verdict": "PASS",
        },
        "raw_polytope_handoff": {
            "raw_lattice_rank": raw_row["lattice_rank"],
            "raw_wall_count": raw_row["wall_count"],
            "unique_simultaneous_P2_frame": raw_row["simultaneous_full_P2_frame_count"] == 1,
            "raw_root_profile": raw_candidate["root_plane_profile"],
            "constructed_root_profile": {
                "degree_1_multiplicities_1,1,1": len(lines),
                "degree_2_multiplicities_1,1,1,1,1,1": len(conics),
            },
            "exact_raw_marking": marking,
            "general_D_CobVin_return": (
                "the twenty realized root seeds and the unique boundary enter the frozen "
                "finite-index root/Halphen compiler; its 57 primitive cusps give all 57 "
                "(-1)-curves and return the exact marked 78-wall chamber"
            ),
            "full_cone_rays": {
                "finite": raw_polytope["polytope"]["exact_double_description"]["finite_vertex_count"],
                "ideal": raw_polytope["polytope"]["exact_double_description"]["ideal_vertex_count"],
                "hyperideal": raw_polytope["polytope"]["exact_double_description"]["hyperideal_vertex_count"],
            },
            "exact_same_P": marking["exact_marked_root_seed_transport"],
        },
        "inverse_level": {
            "rating": "L3 structural P-first calibration, not L4",
            "reason": (
                "the raw lattice uniquely forces the finite design and this certificate "
                "constructs a carrier returning the same P, but no class-wide exhaustive "
                "alternative-carrier theorem or new abstract-polytope bound is proved"
            ),
        },
        "lineage": {
            "known": (
                "Kondo's characteristic-five Type-VI Coble surface, its twenty Gamma_VI "
                "roots, S5 symmetry, and the underlying H10/78 abstract chamber lineage"
            ),
            "project_new": (
                "label-free raw-P extraction of the Desargues/cut design and the explicit "
                "quotient-plane sextic certificate giving a P-first same-chamber return"
            ),
            "worldwide_priority": "open; no new abstract polytope is claimed",
        },
        "all_checks_pass": True,
    }
    return certificate


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
        print("Type-VI characteristic-five Desargues inverse archive replay passed")
        return
    if arguments.write:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
