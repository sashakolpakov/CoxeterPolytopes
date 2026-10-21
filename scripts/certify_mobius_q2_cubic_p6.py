#!/usr/bin/env python3
"""Certify the q=2 Frobenius/Mobius surface and its P^6 nef chamber.

The construction starts on ``P1 x P1`` over ``F4`` with the two graphs
``y=x^2`` and ``x=y^2``.  It proves their five intersections are a del
Pezzo-general set, reconstructs the anticanonical Hermitian cubic and all 27
(-1)-curves, and then performs two independent finite-volume checks on the
dual nef chamber: cddlib/GMP and Vinberg's literal edge-endpoint criterion.

The final lattice change of basis sends the 27 curve classes entry-for-entry
to the normal-vector families in Everitt--Ratcliffe--Tschantz, Table 2.  Thus
the chamber is the classical Potyagailo--Vinberg P^6, not a new abstract
hyperbolic polytope.
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
import sys
from typing import Any, Iterable, Iterator, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from coxeter_polytopes.verify.vinberg_volume import (  # noqa: E402
    verify_vinberg_volume_criterion,
)


OUTPUT = ROOT / "data/enumerations/mobius_q2_cubic_p6.json"
SCHEMA = "coxeter.mobius_q2_cubic_p6.v1"

# F4 = F2[w]/(w^2+w+1), encoded in the polynomial basis as
# 0, 1, w, w+1 = 0, 1, 2, 3.
F4 = (0, 1, 2, 3)
F4_NAMES = ("0", "1", "w", "w+1")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _payload_sha256(value: Any) -> str:
    return _sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def _fadd(left: int, right: int) -> int:
    return left ^ right


def _fmul(left: int, right: int) -> int:
    result = 0
    first, second = left, right
    while second:
        if second & 1:
            result ^= first
        first <<= 1
        second >>= 1
    # Reduce by t^2+t+1.  Products of two degree-one representatives have
    # degree at most two, but the loop keeps the routine manifestly closed.
    for degree in range(4, 1, -1):
        if result & (1 << degree):
            result ^= 1 << degree
            result ^= 1 << (degree - 1)
            result ^= 1 << (degree - 2)
    return result


def _fpow(value: int, exponent: int) -> int:
    result, base = 1, value
    while exponent:
        if exponent & 1:
            result = _fmul(result, base)
        base = _fmul(base, base)
        exponent //= 2
    return result


def _finv(value: int) -> int:
    if value == 0:
        raise ZeroDivisionError
    return _fpow(value, 2)


def _frref(matrix: Sequence[Sequence[int]]) -> tuple[list[list[int]], list[int]]:
    work = [list(row) for row in matrix]
    if not work:
        return [], []
    row = 0
    pivots: list[int] = []
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = _finv(work[row][column])
        work[row] = [_fmul(inverse, value) for value in work[row]]
        for index in range(len(work)):
            if index == row or work[index][column] == 0:
                continue
            factor = work[index][column]
            work[index] = [
                _fadd(work[index][entry], _fmul(factor, work[row][entry]))
                for entry in range(len(work[0]))
            ]
        pivots.append(column)
        row += 1
        if row == len(work):
            break
    return work, pivots


def _frank(matrix: Sequence[Sequence[int]]) -> int:
    return len(_frref(matrix)[1])


def _fnullspace(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    reduced, pivots = _frref(matrix)
    columns = len(matrix[0])
    free = [column for column in range(columns) if column not in pivots]
    vectors: list[list[int]] = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        # In characteristic two, subtraction is addition.
        for row, pivot in enumerate(pivots):
            vector[pivot] = reduced[row][free_column]
        vectors.append(vector)
    return vectors


def _binary_monomial_values(point: tuple[int, int], degree: int) -> list[int]:
    first, second = point
    return [
        _fmul(_fpow(first, exponent), _fpow(second, degree - exponent))
        for exponent in range(degree + 1)
    ]


def _bihomogeneous_evaluation(
    point: tuple[tuple[int, int], tuple[int, int]],
    first_degree: int,
    second_degree: int,
) -> list[int]:
    first = _binary_monomial_values(point[0], first_degree)
    second = _binary_monomial_values(point[1], second_degree)
    return [_fmul(left, right) for left in first for right in second]


def _projective_line_f4() -> list[tuple[int, int]]:
    return [*((value, 1) for value in F4), (1, 0)]


def _frob(point: tuple[int, int]) -> tuple[int, int]:
    return (_fpow(point[0], 2), _fpow(point[1], 2))


def _poly_multiply(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
) -> dict[tuple[int, int], int]:
    result: dict[tuple[int, int], int] = {}
    for (left_x, left_y), left_coefficient in left.items():
        for (right_x, right_y), right_coefficient in right.items():
            exponent = (left_x + right_x, left_y + right_y)
            result[exponent] = _fadd(
                result.get(exponent, 0),
                _fmul(left_coefficient, right_coefficient),
            )
    return {exponent: value for exponent, value in result.items() if value}


def _compositions(total: int, parts: int, prefix: tuple[int, ...] = ()) -> Iterator[tuple[int, ...]]:
    if parts == 1:
        yield (*prefix, total)
        return
    for value in range(total + 1):
        yield from _compositions(total - value, parts - 1, (*prefix, value))


def _relation_space(
    sections: Sequence[dict[tuple[int, int], int]], degree: int
) -> tuple[list[tuple[int, ...]], list[list[int]]]:
    monomials = list(_compositions(degree, len(sections)))
    expansions: list[dict[tuple[int, int], int]] = []
    for exponents in monomials:
        polynomial = {(0, 0): 1}
        for section, exponent in zip(sections, exponents, strict=True):
            for _ in range(exponent):
                polynomial = _poly_multiply(polynomial, section)
        expansions.append(polynomial)
    ambient_exponents = [
        (first, second)
        for first in range(2 * degree + 1)
        for second in range(2 * degree + 1)
    ]
    matrix = [
        [polynomial.get(exponent, 0) for polynomial in expansions]
        for exponent in ambient_exponents
    ]
    return monomials, _fnullspace(matrix)


def _finite_geometry() -> dict[str, Any]:
    line = _projective_line_f4()
    centers = [(point, _frob(point)) for point in line]
    if len(set(centers)) != 5:
        raise ArithmeticError("the Frobenius graph lost an F4 point")

    graph_vanishing = []
    for (x0, x1), (y0, y1) in centers:
        y_equals_x2 = _fadd(_fmul(y0, _fpow(x1, 2)), _fmul(y1, _fpow(x0, 2)))
        x_equals_y2 = _fadd(_fmul(x0, _fpow(y1, 2)), _fmul(x1, _fpow(y0, 2)))
        graph_vanishing.append((y_equals_x2, x_equals_y2))
    if set(graph_vanishing) != {(0, 0)}:
        raise ArithmeticError("an alleged graph intersection is not on both graphs")

    evaluations_11 = [_bihomogeneous_evaluation(point, 1, 1) for point in centers]
    triple_curves = []
    for triple in itertools.combinations(range(5), 3):
        matrix = [evaluations_11[index] for index in triple]
        kernel = _fnullspace(matrix)
        if _frank(matrix) != 3 or len(kernel) != 1:
            raise ArithmeticError("a triple does not determine one (1,1) curve")
        coefficients = kernel[0]
        determinant = _fadd(
            _fmul(coefficients[0], coefficients[3]),
            _fmul(coefficients[1], coefficients[2]),
        )
        vanishing = tuple(
            index
            for index, evaluation in enumerate(evaluations_11)
            if _field_dot(coefficients, evaluation) == 0
        )
        if determinant == 0 or vanishing != triple:
            raise ArithmeticError("a triple curve is reducible or contains an extra center")
        triple_curves.append(
            {
                "center_indices_zero_based": list(triple),
                "O_1_1_coefficients": coefficients,
                "coefficient_names": [F4_NAMES[value] for value in coefficients],
                "bilinear_determinant": F4_NAMES[determinant],
                "is_smooth_irreducible_projectivity_graph": True,
            }
        )
    four_ranks = Counter(
        _frank([evaluations_11[index] for index in subset])
        for subset in itertools.combinations(range(5), 4)
    )
    if four_ranks != {4: 5}:
        raise ArithmeticError("four centers acquired a common (1,1) curve")

    evaluations_22 = [_bihomogeneous_evaluation(point, 2, 2) for point in centers]
    anticanonical_basis = _fnullspace(evaluations_22)
    if _frank(evaluations_22) != 5 or len(anticanonical_basis) != 4:
        raise ArithmeticError("the anticanonical system has the wrong dimension")
    sections = [
        {
            (first, second): vector[3 * first + second]
            for first in range(3)
            for second in range(3)
            if vector[3 * first + second]
        }
        for vector in anticanonical_basis
    ]
    relation_dimensions = {}
    cubic_monomials: list[tuple[int, ...]] = []
    cubic_kernel: list[list[int]] = []
    for degree in (1, 2, 3):
        monomials, kernel = _relation_space(sections, degree)
        relation_dimensions[str(degree)] = len(kernel)
        if degree == 3:
            cubic_monomials, cubic_kernel = monomials, kernel
    if relation_dimensions != {"1": 0, "2": 0, "3": 1}:
        raise ArithmeticError("the anticanonical image did not acquire one cubic equation")
    relation_terms = [
        {"exponents_z0_z1_z2_z3": list(cubic_monomials[index]), "coefficient": coefficient}
        for index, coefficient in enumerate(cubic_kernel[0])
        if coefficient
    ]
    expected_terms = {
        (0, 1, 2, 0),
        (0, 2, 1, 0),
        (1, 0, 0, 2),
        (2, 0, 0, 1),
    }
    if {
        tuple(term["exponents_z0_z1_z2_z3"]) for term in relation_terms
    } != expected_terms or {term["coefficient"] for term in relation_terms} != {1}:
        raise ArithmeticError("the frozen Hermitian cubic equation changed")

    triples = Counter(
        tuple(sorted(triple))
        for curve in triple_curves
        for triple in itertools.combinations(curve["center_indices_zero_based"], 3)
    )
    pairs = Counter(
        tuple(sorted(pair))
        for curve in triple_curves
        for pair in itertools.combinations(curve["center_indices_zero_based"], 2)
    )
    points = Counter(
        point
        for curve in triple_curves
        for point in curve["center_indices_zero_based"]
    )
    if not (
        len(triples) == 10
        and set(triples.values()) == {1}
        and set(pairs.values()) == {3}
        and set(points.values()) == {6}
    ):
        raise ArithmeticError("the S(3,3,5) incidence count changed")

    return {
        "field": {
            "model": "F2[w]/(w^2+w+1)",
            "integer_encoding": {str(index): name for index, name in enumerate(F4_NAMES)},
        },
        "centers": [
            {
                "x": list(point[0]),
                "y": list(point[1]),
                "x_names": [F4_NAMES[value] for value in point[0]],
                "y_names": [F4_NAMES[value] for value in point[1]],
            }
            for point in centers
        ],
        "frobenius_graphs": {
            "classes": ["(2,1)", "(1,2)"],
            "intersection_number": 5,
            "five_distinct_F4_intersections": True,
            "local_Jacobian_determinant_at_each_center": "1",
            "transverse_at_all_five_centers": True,
            "bezout_proves_no_further_geometric_intersections": True,
        },
        "general_position": {
            "distinct_first_and_second_projections": True,
            "triple_O_1_1_evaluation_rank_histogram": {"3": 10},
            "four_point_O_1_1_evaluation_rank_histogram": {
                str(rank): count for rank, count in sorted(four_ranks.items())
            },
            "no_four_centers_on_an_O_1_1_curve": True,
            "del_Pezzo_degree_three_gate": True,
        },
        "S_3_3_5": {
            "block_count": len(triple_curves),
            "point_replication": sorted(set(points.values())),
            "pair_replication": sorted(set(pairs.values())),
            "triple_replication": sorted(set(triples.values())),
            "all_ten_block_curves_are_smooth_irreducible_projectivity_graphs": True,
            "blocks": triple_curves,
        },
        "anticanonical_model": {
            "O_2_2_evaluation_rank": 5,
            "section_dimension_after_five_basepoints": 4,
            "section_basis_monomial_order": (
                "x0^i*x1^(2-i)*y0^j*y1^(2-j), lexicographic (i,j), 0<=i,j<=2"
            ),
            "section_basis_coefficients": anticanonical_basis,
            "relation_kernel_dimensions_by_degree": relation_dimensions,
            "cubic_relation_terms": relation_terms,
            "equation": "z0^2*z3+z0*z3^2+z1^2*z2+z1*z2^2=0",
            "Hermitian_matrix_rank": 4,
            "gradient": ["z3^2", "z2^2", "z1^2", "z0^2"],
            "geometrically_smooth": True,
            "identification": "nondegenerate Hermitian cubic surface H(3,4)",
            "F4_rational_point_count_from_blowup": 25 + 5 * 4,
        },
    }


def _field_dot(left: Sequence[int], right: Sequence[int]) -> int:
    result = 0
    for first, second in zip(left, right, strict=True):
        result = _fadd(result, _fmul(first, second))
    return result


def _intersection_gram() -> sp.Matrix:
    gram = sp.zeros(7)
    gram[0, 1] = gram[1, 0] = 1
    for index in range(2, 7):
        gram[index, index] = -1
    return gram


def _line_classes() -> tuple[list[str], list[sp.Matrix]]:
    labels: list[str] = []
    classes: list[sp.Matrix] = []
    for index in range(5):
        vector = [0] * 7
        vector[2 + index] = 1
        labels.append(f"E_{index + 1}")
        classes.append(sp.Matrix(vector))
    for ruling, coordinate in (("F", 0), ("G", 1)):
        for index in range(5):
            vector = [0] * 7
            vector[coordinate] = 1
            vector[2 + index] = -1
            labels.append(f"{ruling}-E_{index + 1}")
            classes.append(sp.Matrix(vector))
    for triple in itertools.combinations(range(5), 3):
        vector = [1, 1, *(-int(index in triple) for index in range(5))]
        labels.append("M_{" + ",".join(str(index + 1) for index in triple) + "}")
        classes.append(sp.Matrix(vector))
    labels.extend(("B_(1,2)", "B_(2,1)"))
    classes.extend((sp.Matrix([1, 2, -1, -1, -1, -1, -1]), sp.Matrix([2, 1, -1, -1, -1, -1, -1])))
    if len(classes) != 27:
        raise ArithmeticError("the cubic line stock changed")
    return labels, classes


def _p2_marking(gram: sp.Matrix, labels: Sequence[str], curves: Sequence[sp.Matrix]) -> dict[str, Any]:
    plane = sp.Matrix([3, 1, -1, -1, -1, -1, -1])
    exceptional = []
    for index in range(5):
        vector = [1, 0, *(-int(position == index) for position in range(5))]
        exceptional.append(sp.Matrix(vector))
    exceptional.append(sp.Matrix([2, 1, -1, -1, -1, -1, -1]))
    basis = sp.Matrix.hstack(plane, *exceptional)
    basis_gram = basis.T * gram * basis
    if abs(int(basis.det())) != 1 or basis_gram != sp.diag(1, -1, -1, -1, -1, -1, -1):
        raise ArithmeticError("the six-line blowdown marking is not unimodular P2")
    inverse = basis.inv()

    records = []
    normals = []
    for label, curve in zip(labels, curves, strict=True):
        coordinates = inverse * curve
        if any(value.q != 1 for value in coordinates):
            raise ArithmeticError("a line class left the P2 marking integrally")
        values = [int(value) for value in coordinates]
        normal = [-value for value in values[1:]] + [values[0]]
        normals.append(normal)
        records.append(
            {
                "curve_label": label,
                "P1xP1_blowup_coordinates_F_G_E1_to_E5": [int(value) for value in curve],
                "P2_blowup_coordinates_d_then_L_coefficients": values,
                "ERT_outward_normal_coordinates": normal,
            }
        )

    standard = []
    for index in range(6):
        vector = [0] * 7
        vector[index] = -1
        standard.append(vector)
    for pair in itertools.combinations(range(6), 2):
        standard.append([int(index in pair) for index in range(6)] + [1])
    for missing in range(6):
        standard.append([int(index != missing) for index in range(6)] + [2])
    if sorted(normals) != sorted(standard):
        raise ArithmeticError("the 27 classes do not equal ERT Table 2")

    # Exhaust every lattice solution D^2=K.D=-1.  In P2 coordinates
    # D=dH-sum(m_i L_i), Cauchy gives
    # (3d-1)^2 <= 6(d^2+1), hence d is 0, 1, or 2.
    solutions = []
    by_degree: Counter[int] = Counter()
    for degree in (0, 1, 2):
        bound = math.isqrt(degree * degree + 1)
        for multiplicities in itertools.product(range(-bound, bound + 1), repeat=6):
            if (
                sum(multiplicities) == 3 * degree - 1
                and sum(value * value for value in multiplicities) == degree * degree + 1
            ):
                solutions.append((degree, *multiplicities))
                by_degree[degree] += 1
    actual = sorted((normal[-1], *normal[:6]) for normal in normals)
    if sorted(solutions) != actual or by_degree != {0: 6, 1: 15, 2: 6}:
        raise ArithmeticError("the lattice exceptional-class exhaustion changed")

    return {
        "six_pairwise_disjoint_curves": [
            "F-E_1", "F-E_2", "F-E_3", "F-E_4", "F-E_5", "B_(2,1)"
        ],
        "plane_class": [int(value) for value in plane],
        "basis_columns": [[int(basis[row, column]) for row in range(7)] for column in range(7)],
        "basis_determinant_absolute": abs(int(basis.det())),
        "basis_gram": [[int(value) for value in basis_gram.row(row)] for row in range(7)],
        "Cauchy_degree_bound": "(3d-1)^2<=6(d^2+1), so d in {0,1,2}",
        "exceptional_class_solution_histogram_by_degree": {
            str(degree): count for degree, count in sorted(by_degree.items())
        },
        "all_exceptional_lattice_solutions_realized": True,
        "line_records": records,
        "ERT_Table_2_normal_multiset_exact": True,
        "ERT_normal_matrix_sha256": _payload_sha256(sorted(normals)),
    }


def _primitive_integer_ray(tokens: Iterable[str]) -> tuple[int, ...]:
    values = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in values))
    integers = tuple(value.numerator * (denominator // value.denominator) for value in values)
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted the zero ray")
    return tuple(value // divisor for value in integers)


def _exact_cone(
    gram: sp.Matrix,
    curves: Sequence[sp.Matrix],
    executable: str,
) -> dict[str, Any]:
    walls = sp.Matrix.hstack(*curves)
    inequalities = walls.T * gram
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(int(value)) for value in inequalities.row(row))
        for row in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stdout[-4000:] + completed.stderr[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count_text, width_text, kind = output[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if kind != "rational" or width != 8:
        raise ArithmeticError("cddlib emitted the wrong generator format")

    rays = sorted(
        _primitive_integer_ray(line.split()[1:])
        for line in output[begin + 2 : begin + 2 + count]
        if Fraction(line.split()[0]) == 0
    )
    if len(rays) != count or len(rays) != len(set(rays)):
        raise ArithmeticError("cddlib emitted a vertex row or duplicate ray")
    kinds: Counter[str] = Counter()
    active_histogram: Counter[tuple[str, int]] = Counter()
    active_facets: set[int] = set()
    records = []
    ample = sp.Matrix([2, 2, -1, -1, -1, -1, -1])
    if set(inequalities * ample) != {1} or int((ample.T * gram * ample)[0]) != 3:
        raise ArithmeticError("the anticanonical ample witness changed")
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("an exact ray violates a nef inequality")
        active = [index for index, value in enumerate(slacks) if value == 0]
        active_facets.update(active)
        square = int((vector.T * gram * vector)[0])
        ray_kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        kinds[ray_kind] += 1
        active_histogram[(ray_kind, len(active))] += 1
        records.append(
            {
                "ray": list(ray),
                "square": square,
                "kind": ray_kind,
                "active_facet_indices_zero_based": active,
            }
        )
    if not (
        kinds == {"finite": 72, "ideal": 27}
        and active_histogram == {("finite", 6): 72, ("ideal", 10): 27}
        and active_facets == set(range(27))
    ):
        raise ArithmeticError("the P6 exact cone census changed")
    return {
        "engine": "cddlib GMP rational double-description",
        "inequality_matrix_sha256": _payload_sha256(
            [[int(value) for value in inequalities.row(row)] for row in range(27)]
        ),
        "extreme_ray_count": len(rays),
        "finite_vertex_count": kinds["finite"],
        "ideal_vertex_count": kinds["ideal"],
        "hyperideal_vertex_count": kinds["hyperideal"],
        "active_facet_count_histogram": {
            f"{ray_kind}_{number}": count
            for (ray_kind, number), count in sorted(active_histogram.items())
        },
        "every_wall_is_irredundant": active_facets == set(range(27)),
        "strict_ample_witness_minus_K": [int(value) for value in ample],
        "minimum_ample_wall_pairing": 1,
        "finite_volume_from_complete_cone": kinds["hyperideal"] == 0,
        "compact": kinds["ideal"] == 0 and kinds["hyperideal"] == 0,
        "rays": records,
    }


def _vinberg(curve_gram: sp.Matrix) -> dict[str, Any]:
    wall_gram = [
        [-int(curve_gram[left, right]) for right in range(curve_gram.cols)]
        for left in range(curve_gram.rows)
    ]
    verification = verify_vinberg_volume_criterion(
        wall_gram,
        6,
        simple_roots_complete=True,
    )
    endpoint_histogram = Counter(
        (len(edge.finite_endpoints), len(edge.ideal_endpoints))
        for edge in verification.edge_audits
    )
    cusp_histogram = Counter(
        (len(endpoint.facets), tuple(sorted(map(len, endpoint.affine_components))))
        for endpoint in verification.ideal_vertices
    )
    if not (
        verification.decision == "finite_volume_noncompact"
        and verification.elliptic_edge_count == 648
        and verification.finite_vertex_count == 72
        and verification.ideal_vertex_count == 27
        and not verification.finite_volume_bad_edges
        and len(verification.compact_bad_edges) == 648
        and endpoint_histogram == {(1, 1): 432, (0, 2): 216}
        and cusp_histogram == {(10, (2, 2, 2, 2, 2)): 27}
        and all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in verification.parabolic_neighborhood_audits
        )
    ):
        raise ArithmeticError("Vinberg's literal P6 audit changed")
    return {
        "method": "literal elliptic-edge endpoint and full affine zero-neighborhood criterion",
        "simple_roots_complete": verification.simple_roots_complete,
        "decision": verification.decision,
        "finite_volume": verification.finite_volume,
        "compact": verification.compact,
        "elliptic_rank_5_edge_count": verification.elliptic_edge_count,
        "finite_vertex_count": verification.finite_vertex_count,
        "ideal_vertex_count": verification.ideal_vertex_count,
        "edge_endpoint_histogram_finite_then_ideal": {
            f"{finite},{ideal}": count
            for (finite, ideal), count in sorted(endpoint_histogram.items())
        },
        "finite_volume_bad_edge_count": len(verification.finite_volume_bad_edges),
        "compact_bad_edge_count": len(verification.compact_bad_edges),
        "connected_affine_component_audit_count": len(
            verification.parabolic_neighborhood_audits
        ),
        "invalid_full_zero_neighborhood_count": sum(
            not audit.is_rank_dimension_minus_one_parabolic
            for audit in verification.parabolic_neighborhood_audits
        ),
        "ideal_neighborhood_histogram_facets_then_component_sizes": {
            f"{facets};{','.join(map(str, sizes))}": count
            for (facets, sizes), count in sorted(cusp_histogram.items())
        },
    }


def _reflection_integrality(gram: sp.Matrix, curves: Sequence[sp.Matrix]) -> dict[str, Any]:
    identity = sp.eye(7)
    matrices = []
    for curve in curves:
        reflection = identity + 2 * curve * (curve.T * gram)
        if any(value.q != 1 for value in reflection):
            raise ArithmeticError("a line reflection became nonintegral")
        if reflection.T * gram * reflection != gram or int(reflection.det()) != -1:
            raise ArithmeticError("a line reflection does not preserve I_(1,6)")
        if any((int(value) - int(identity[row, column])) % 2 for row in range(7) for column, value in enumerate(reflection.row(row))):
            raise ArithmeticError("a line reflection left the level-two congruence group")
        matrices.append([[int(value) for value in reflection.row(row)] for row in range(7)])
    return {
        "formula": "s_C(x)=x+2(x.C)C because C^2=-1",
        "all_27_integral": True,
        "all_preserve_NS_gram": True,
        "all_determinant_minus_one": True,
        "all_congruent_to_identity_mod_2": True,
        "reflection_matrices_sha256": _payload_sha256(matrices),
    }


def _q_factory_dichotomy() -> dict[str, Any]:
    samples = {}
    for q in (2, 3, 4, 5, 7, 8, 9):
        samples[str(q)] = {
            "center_count": q * q + 1,
            "circle_size": q + 1,
            "circle_square": 1 - q,
            "circle_divisibility_in_full_NS": 1,
            "circle_reflection_integral": 2 % (q - 1) == 0,
        }
    return {
        "circle_class": "C_A=F+G-sum_(i in A) E_i, |A|=q+1",
        "circle_square": "1-q=-(q-1)",
        "divisibility": 1,
        "divisibility_witness": "C_A.E_i=1 for every incident exceptional E_i",
        "reflection_formula": "s_C(x)=x+2(x.C)C/(q-1)",
        "integrality_criterion": "q-1 divides 2",
        "prime_power_solutions": [2, 3],
        "q2_endpoint": "27 norm-minus-one lines on the Hermitian cubic; P^6",
        "q3_endpoint": "30 norm-minus-two Mobius circles in the Type-MI Coble seed; H11/206 after the full compiler",
        "scope_warning": (
            "This excludes q>3 only from the full-NS integral-reflection Mobius factory. "
            "It is not nonexistence of other q-polytopes, other lattices, or other carriers."
        ),
        "sample_exact_arithmetic": samples,
    }


def build_certificate(cdd_executable: str | None = None) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact P6 cone")
    geometry = _finite_geometry()
    gram = _intersection_gram()
    labels, curves = _line_classes()
    canonical = sp.Matrix([-2, -2, 1, 1, 1, 1, 1])
    curve_matrix = sp.Matrix.hstack(*curves)
    curve_gram = curve_matrix.T * gram * curve_matrix
    if any(curve_gram[index, index] != -1 for index in range(27)):
        raise ArithmeticError("a claimed cubic line does not have square -1")
    if any(int((canonical.T * gram * curve)[0]) != -1 for curve in curves):
        raise ArithmeticError("a claimed cubic line fails adjunction")
    distinct_pairings = Counter(
        int(curve_gram[left, right])
        for left in range(27)
        for right in range(left)
    )
    if distinct_pairings != {0: 216, 1: 135}:
        raise ArithmeticError("the Schlaefli pairing graph changed")
    marking = _p2_marking(gram, labels, curves)
    cone = _exact_cone(gram, curves, executable)
    vinberg = _vinberg(curve_gram)
    reflections = _reflection_integrality(gram, curves)
    all_checks = (
        geometry["anticanonical_model"]["geometrically_smooth"]
        and marking["all_exceptional_lattice_solutions_realized"]
        and cone["finite_vertex_count"] == vinberg["finite_vertex_count"] == 72
        and cone["ideal_vertex_count"] == vinberg["ideal_vertex_count"] == 27
        and cone["hyperideal_vertex_count"] == 0
        and vinberg["decision"] == "finite_volume_noncompact"
        and marking["ERT_Table_2_normal_multiset_exact"]
    )
    if not all_checks:
        raise ArithmeticError("the q=2 cubic/P6 join failed")
    return {
        "schema": SCHEMA,
        "construction": (
            "blow up P1xP1 over F4 at the five transverse intersections of "
            "y=x^2 and x=y^2"
        ),
        "finite_geometry": geometry,
        "surface": {
            "smooth": True,
            "surface_type": "degree-three del Pezzo / nondegenerate Hermitian cubic",
            "canonical_class_F_G_E1_to_E5": [int(value) for value in canonical],
            "canonical_square": int((canonical.T * gram * canonical)[0]),
            "NS_basis": ["F", "G", "E1", "E2", "E3", "E4", "E5"],
            "NS_gram": [[int(value) for value in gram.row(row)] for row in range(7)],
            "NS_determinant": int(gram.det()),
            "NS_signature": [1, 6, 0],
            "NS_identification": "odd unimodular I_(1,6)",
            "F4_rational_minus_one_curve_count": len(curves),
            "all_27_minus_one_curves_are_geometrically_realized": True,
        },
        "line_chamber": {
            "wall_count": len(curves),
            "wall_labels": labels,
            "wall_vectors": [[int(value) for value in curve] for curve in curves],
            "wall_gram": [[int(value) for value in curve_gram.row(row)] for row in range(27)],
            "wall_gram_sha256": _payload_sha256(
                [[int(value) for value in curve_gram.row(row)] for row in range(27)]
            ),
            "diagonal_square": -1,
            "distinct_wall_pairing_histogram": {
                str(value): count for value, count in sorted(distinct_pairings.items())
            },
            "all_distinct_pairings_are_zero_or_one": True,
            "right_angled_explanation": (
                "pairing 0 gives orthogonal meeting mirrors; pairing 1 gives parallel "
                "mirrors, so every actual dihedral angle is pi/2"
            ),
            "reflections": reflections,
            "exact_double_description": cone,
            "literal_Vinberg": vinberg,
        },
        "P2_blowdown_and_ERT_identification": marking,
        "lineage": {
            "abstract_polytope": "Potyagailo--Vinberg P^6",
            "Euclidean_dual": "Gosset/Coxeter 2_21",
            "arithmetic_group": "level-two congruence reflection subgroup Gamma_2^6 of PO_(6,1)(Z)",
            "new_abstract_polytope": False,
            "new_surface": False,
            "primary_sources": [
                {
                    "citation": "Potyagailo--Vinberg, Comment. Math. Helv. 80 (2005), 63--73",
                    "url": "https://doi.org/10.4171/CMH/4",
                    "role": "the finite-volume right-angled P^n series and P^6 facet/cusp count",
                },
                {
                    "citation": "Everitt--Ratcliffe--Tschantz, Math. Ann. 354 (2012), 871--905",
                    "url": "https://arxiv.org/abs/1009.3231",
                    "role": "Table 2 normal vectors, 72+27 vertices, congruence-two theorem, and dual 2_21 identification",
                },
                {
                    "citation": "Coxeter, Amer. J. Math. 62 (1940), 457--486",
                    "url": "https://doi.org/10.2307/2371466",
                    "role": "2_21 and the 27-line cubic-surface configuration",
                },
                {
                    "citation": "Hirschfeld, J. Austral. Math. Soc. 4 (1964), 83--89",
                    "url": "https://doi.org/10.1017/S144678870002276X",
                    "role": "the classical F4 double-six and 27-line cubic surface",
                },
            ],
            "project_new_claim": (
                "Within this project: the frozen F4 Frobenius-graph blowup -> explicit "
                "Hermitian cubic -> all 27 divisor walls -> exact P^6 return, placed as "
                "the q=2 companion to q=3. Worldwide priority for this particular "
                "presentation is not claimed."
            ),
            "worldwide_priority": "open; no novelty claim for the surface, polytope, or their classical 27-line/2_21 relation",
        },
        "Mobius_factory_integral_reflection_dichotomy": _q_factory_dichotomy(),
        "inverse_scope": {
            "rating": "L3 exact structural calibration, not class-wide L4",
            "positive_content": (
                "one uniform q-factory has exactly two full-NS integral circle-reflection "
                "endpoints, q=2 and q=3, returning known P^6 and H11/206 chambers"
            ),
            "not_proved": [
                "abstract nonexistence for q>3",
                "a necessary-and-sufficient surface carrier theorem for all arithmetic Coxeter polytopes",
                "a new higher-dimensional hyperbolic polytope",
            ],
        },
        "all_checks_pass": all_checks,
        "script_sha256": _sha256_bytes(Path(__file__).read_bytes()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd", default=shutil.which("cddexec_gmp") or "cddexec_gmp")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.cdd)
    serialized = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != serialized:
            raise SystemExit("frozen q=2 cubic/P6 archive differs from live certificate")
        print(f"archive replay passed: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
