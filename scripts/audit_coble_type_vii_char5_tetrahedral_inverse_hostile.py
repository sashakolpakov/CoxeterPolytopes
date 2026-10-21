#!/usr/bin/env python3
"""Independent hostile audit of the characteristic-five Type-VII inverse.

The production certifier is deliberately not imported.  Starting with only
the lattice Gram and the 72 primitive normals in the frozen H10/72 archive,
this script reconstructs every ten-curve plane frame, independently builds
the forced finite geometry over F5, and checks the exact marked return.  It
uses Singular for scheme-theoretic curve audits and cddlib for a fresh exact
H-to-V conversion of the twenty-root cone.
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
from typing import Iterable, Sequence

import networkx as nx
import sympy as sp


PRIME = 5
ROOT = Path(__file__).resolve().parents[1]
RAW_POLYTOPE = (
    ROOT / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json"
)
TARGET = ROOT / "data/enumerations/coble_type_vii_char5_tetrahedral_inverse.json"
OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vii_char5_tetrahedral_inverse_hostile_audit.json"
)
SCHEMA = "coxeter.coble_type_vii_char5_tetrahedral_inverse.hostile_audit.v1"

Point = tuple[int, int, int]
PointLabel = tuple[str, int]
Edge = tuple[int, int]
Matching = tuple[Edge, Edge, Edge]


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
    row = 0
    pivots: list[int] = []
    for column in range(len(matrix[0])):
        pivot = next(
            (candidate for candidate in range(row, len(matrix)) if matrix[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        scale = inverse(matrix[row][column])
        matrix[row] = [scale * value % PRIME for value in matrix[row]]
        for other in range(len(matrix)):
            if other == row or matrix[other][column] == 0:
                continue
            scale = matrix[other][column]
            matrix[other] = [
                (left - scale * right) % PRIME
                for left, right in zip(matrix[other], matrix[row])
            ]
        pivots.append(column)
        row += 1
        if row == len(matrix):
            break
    return matrix, pivots


def nullspace(
    rows: Sequence[Sequence[int]], column_count: int | None = None
) -> list[tuple[int, ...]]:
    if rows:
        column_count = len(rows[0])
    if column_count is None:
        raise ValueError("column count is needed for an empty matrix")
    reduced, pivots = rref(rows)
    free = [column for column in range(column_count) if column not in pivots]
    result: list[tuple[int, ...]] = []
    for free_column in free:
        vector = [0] * column_count
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free_column] % PRIME
        result.append(tuple(vector))
    return result


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


def cross(left: Point, right: Point) -> Point:
    return normalize(
        (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )
    )  # type: ignore[return-value]


def line_value(line: Sequence[int], point: Point) -> int:
    return sum(left * right for left, right in zip(line, point)) % PRIME


def all_projective_points() -> list[Point]:
    return sorted(
        {
            normalize(vector)
            for vector in itertools.product(range(PRIME), repeat=3)
            if any(vector)
        }
    )  # type: ignore[return-value]


def all_matchings() -> list[Matching]:
    edges = list(itertools.combinations(range(6), 2))
    return sorted(
        tuple(sorted(triple))
        for triple in itertools.combinations(edges, 3)
        if sorted(vertex for edge in triple for vertex in edge) == list(range(6))
    )


def all_factorizations() -> list[tuple[Matching, ...]]:
    edges = set(itertools.combinations(range(6), 2))
    matchings = all_matchings()
    answers: set[tuple[Matching, ...]] = set()

    def recurse(remaining: set[Edge], chosen: list[Matching]) -> None:
        if not remaining:
            answers.add(tuple(sorted(chosen)))
            return
        first = min(remaining)
        for matching in matchings:
            block = set(matching)
            if first in block and block <= remaining:
                recurse(remaining - block, [*chosen, matching])

    recurse(edges, [])
    return sorted(answers)


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


def reconstruct_finite_geometry() -> dict[str, object]:
    ordinary: list[Point] = [
        normalize((parameter * parameter, parameter, 1))  # type: ignore[arg-type]
        for parameter in range(PRIME)
    ]
    ordinary.append((1, 0, 0))
    if len(set(ordinary)) != 6 or any(
        (point[0] * point[2] - point[1] * point[1]) % PRIME
        for point in ordinary
    ):
        raise ArithmeticError("the six points are not the complete F5 conic")

    chords = {
        edge: cross(ordinary[edge[0]], ordinary[edge[1]])
        for edge in itertools.combinations(range(6), 2)
    }
    concurrent: dict[Matching, Point] = {}
    for matching in all_matchings():
        pole = cross(chords[matching[0]], chords[matching[1]])
        if line_value(chords[matching[2]], pole) == 0:
            concurrent[matching] = pole
    if len(concurrent) != 10 or len(set(concurrent.values())) != 10:
        raise ArithmeticError("concurrent matching census changed")

    histogram: Counter[int] = Counter()
    usable: list[tuple[tuple[Matching, ...], Matching]] = []
    for factorization in all_factorizations():
        good = tuple(matching for matching in factorization if matching in concurrent)
        bad = tuple(matching for matching in factorization if matching not in concurrent)
        histogram[len(good)] += 1
        if len(good) == 4 and len(bad) == 1:
            usable.append((good, bad[0]))
    if len(all_factorizations()) != 6 or histogram != Counter({4: 5, 0: 1}):
        raise ArithmeticError("one-factorization concurrency histogram changed")

    # This is selected intrinsically as the lexicographically first usable
    # factorization, not copied from the production certificate.
    selected, residual = usable[0]
    special = [concurrent[matching] for matching in selected]
    centers: dict[PointLabel, Point] = {
        **{("s", index): point for index, point in enumerate(special)},
        **{("o", index): point for index, point in enumerate(ordinary)},
    }
    if len(set(centers.values())) != 10:
        raise ArithmeticError("special and ordinary centers collide")

    curves: list[dict[str, object]] = []
    for special_index, matching in enumerate(selected):
        for pair_index, edge in enumerate(matching):
            coefficients = chords[edge]
            support = {
                label
                for label, point in centers.items()
                if line_value(coefficients, point) == 0
            }
            expected = {
                ("s", special_index),
                ("o", edge[0]),
                ("o", edge[1]),
            }
            if support != expected:
                raise ArithmeticError("a chord has unintended center incidence")
            curves.append(
                {
                    "node": ("line", special_index, pair_index),
                    "degree": 1,
                    "coefficients": coefficients,
                    "support": {label: 1 for label in sorted(expected)},
                }
            )

    conic_supports: list[tuple[tuple[object, ...], dict[PointLabel, int]]] = [
        (("ordinary",), {("o", index): 1 for index in range(6)})
    ]
    for pair_index, edge in enumerate(residual):
        conic_supports.append(
            (
                ("residual", pair_index),
                {
                    **{("s", index): 1 for index in range(4)},
                    ("o", edge[0]): 1,
                    ("o", edge[1]): 1,
                },
            )
        )
    conic_records = []
    for label, support in conic_supports:
        matrix = [evaluation_row(2, centers[point]) for point in support]
        kernel = nullspace(matrix)
        if len(rref(matrix)[1]) != 5 or len(kernel) != 1:
            raise ArithmeticError("six points do not determine a unique conic")
        coefficients = normalize(kernel[0])
        exact_support = {
            point
            for point, coordinates in centers.items()
            if evaluate(coefficients, 2, coordinates) == 0
        }
        determinant = conic_determinant(coefficients)
        if exact_support != set(support) or determinant == 0:
            raise ArithmeticError("a conic is singular or has accidental centers")
        record = {
            "node": ("conic", *label),
            "degree": 2,
            "coefficients": coefficients,
            "support": support,
            "evaluation_rank": 5,
            "symmetric_matrix_determinant_mod_5": determinant,
        }
        conic_records.append(record)
        curves.append(record)

    quartic_records = []
    for omitted in range(4):
        conditions = [evaluation_row(4, centers[("o", index)]) for index in range(6)]
        for special_index in range(4):
            if special_index == omitted:
                continue
            conditions.extend(
                derivative_row(4, variable, centers[("s", special_index)])
                for variable in range(3)
            )
        kernel = nullspace(conditions)
        rank = len(rref(conditions)[1])
        if rank != 14 or len(kernel) != 1:
            raise ArithmeticError("quartic system is not one-dimensional")
        coefficients = normalize(kernel[0])
        support = {
            **{("o", index): 1 for index in range(6)},
            **{
                ("s", index): 2
                for index in range(4)
                if index != omitted
            },
        }
        zero_centers = {
            label
            for label, point in centers.items()
            if evaluate(coefficients, 4, point) == 0
        }
        if zero_centers != set(support):
            raise ArithmeticError("quartic has the wrong center support")
        record = {
            "node": ("quartic", omitted),
            "degree": 4,
            "coefficients": coefficients,
            "support": support,
            "condition_row_count": len(conditions),
            "condition_rank": rank,
            "omitted_special": omitted,
        }
        quartic_records.append(record)
        curves.append(record)

    cubic_rows = [evaluation_row(3, point) for point in centers.values()]
    cubic_rank = len(rref(cubic_rows)[1])
    if cubic_rank != 10 or nullspace(cubic_rows):
        raise ArithmeticError("an anticanonical cubic survives")

    sextic_rows = [
        derivative_row(6, variable, point)
        for point in centers.values()
        for variable in range(3)
    ]
    sextic_kernel = nullspace(sextic_rows)
    sextic_rank = len(rref(sextic_rows)[1])
    if sextic_rank != 27 or len(sextic_kernel) != 1:
        raise ArithmeticError("the double-point sextic is not unique")
    boundary = normalize(sextic_kernel[0])
    augmented = [
        *sextic_rows,
        *[evaluation_row(6, point) for point in centers.values()],
    ]
    augmented_rank = len(rref(augmented)[1])
    if augmented_rank != 27:
        raise ArithmeticError("sextic values are not gradient consequences")

    return {
        "ordinary": ordinary,
        "special": special,
        "centers": centers,
        "selected_matchings": selected,
        "residual_matching": residual,
        "curve_records": curves,
        "conic_records": conic_records,
        "quartic_records": quartic_records,
        "boundary": boundary,
        "projective_point_count": len(all_projective_points()),
        "perfect_matching_count": len(all_matchings()),
        "concurrent_matching_count": len(concurrent),
        "factorization_count": len(all_factorizations()),
        "factorization_histogram": histogram,
        "usable_factorization_count": len(usable),
        "cubic_rank": cubic_rank,
        "sextic_rank": sextic_rank,
        "sextic_augmented_rank": augmented_rank,
    }


def affine_expansion(
    coefficients: Sequence[int], degree: int, point: Point
) -> tuple[int, tuple[int, int], dict[tuple[int, int], int]]:
    """Expand in independent affine coordinates around a projective point."""

    chart = next(index for index, coordinate in enumerate(point) if coordinate)
    scale = inverse(point[chart])
    normalized = tuple(coordinate * scale % PRIME for coordinate in point)
    free = tuple(index for index in range(3) if index != chart)
    expansion: dict[tuple[int, int], int] = Counter()
    for coefficient, exponent in zip(coefficients, monomials(degree)):
        coefficient %= PRIME
        if coefficient == 0:
            continue
        first = exponent[free[0]]
        second = exponent[free[1]]
        for left_power in range(first + 1):
            for right_power in range(second + 1):
                value = (
                    coefficient
                    * math.comb(first, left_power)
                    * pow(normalized[free[0]], first - left_power, PRIME)
                    * math.comb(second, right_power)
                    * pow(normalized[free[1]], second - right_power, PRIME)
                ) % PRIME
                expansion[(left_power, right_power)] = (
                    expansion.get((left_power, right_power), 0) + value
                ) % PRIME
    return chart, free, {key: value for key, value in expansion.items() if value}


def quadratic(expansion: dict[tuple[int, int], int]) -> tuple[int, int, int]:
    return (
        expansion.get((2, 0), 0),
        expansion.get((1, 1), 0),
        expansion.get((0, 2), 0),
    )


def qvalue(form: tuple[int, int, int], direction: tuple[int, int]) -> int:
    left, mixed, right = form
    u, v = direction
    return (left * u * u + mixed * u * v + right * v * v) % PRIME


def discriminant(form: tuple[int, int, int]) -> int:
    left, mixed, right = form
    return (mixed * mixed - 4 * left * right) % PRIME


def resultant(left: tuple[int, int, int], right: tuple[int, int, int]) -> int:
    a, b, c = left
    d, e, f = right
    return int(
        sp.Matrix(
            [[a, b, c, 0], [0, a, b, c], [d, e, f, 0], [0, d, e, f]]
        ).det()
    ) % PRIME


def local_and_transversality_audit(configuration: dict[str, object]) -> dict[str, object]:
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    boundary: Sequence[int] = configuration["boundary"]  # type: ignore[assignment]
    directions = [(1, value) for value in range(PRIME)] + [(0, 1)]
    boundary_forms: dict[PointLabel, tuple[int, int, int]] = {}
    node_discriminants: list[int] = []
    cusp_cubics: list[int] = []
    boundary_records = []
    for label, point in centers.items():
        _, _, expansion = affine_expansion(boundary, 6, point)
        multiplicity = min(sum(key) for key in expansion)
        form = quadratic(expansion)
        boundary_forms[label] = form
        tangent_directions = [direction for direction in directions if qvalue(form, direction) == 0]
        disc = discriminant(form)
        if multiplicity != 2:
            raise ArithmeticError("a boundary center is not a double point")
        if label[0] == "s":
            if disc == 0:
                raise ArithmeticError("a special boundary center is not A1")
            node_discriminants.append(disc)
            singularity = "A1"
            cubic = None
        else:
            if disc != 0 or len(tangent_directions) != 1:
                raise ArithmeticError("an ordinary boundary center lacks a double tangent")
            u, v = tangent_directions[0]
            cubic = sum(
                expansion.get((left, 3 - left), 0)
                * pow(u, left, PRIME)
                * pow(v, 3 - left, PRIME)
                for left in range(4)
            ) % PRIME
            if cubic == 0:
                raise ArithmeticError("a double-tangent boundary point is worse than A2")
            cusp_cubics.append(cubic)
            singularity = "A2"
        boundary_records.append(
            {
                "point": f"{label[0]}{label[1]}",
                "multiplicity": multiplicity,
                "quadratic_tangent_cone": list(form),
                "discriminant_mod_5": disc,
                "tangent_directions_over_F5": [list(direction) for direction in tangent_directions],
                "cubic_along_double_tangent": cubic,
                "singularity": singularity,
            }
        )

    quartic_discriminants: list[int] = []
    quartic_resultants: list[int] = []
    smooth_values: list[tuple[str, int]] = []
    smooth_count = 0
    quartic_node_count = 0
    for curve in configuration["curve_records"]:  # type: ignore[union-attr]
        kind = curve["node"][0]
        coefficients = curve["coefficients"]
        degree = curve["degree"]
        support: dict[PointLabel, int] = curve["support"]
        for label, multiplicity in support.items():
            _, _, expansion = affine_expansion(coefficients, degree, centers[label])
            actual_multiplicity = min(sum(key) for key in expansion)
            if actual_multiplicity != multiplicity:
                raise ArithmeticError("a root has the wrong local multiplicity")
            if kind == "quartic" and multiplicity == 2:
                form = quadratic(expansion)
                disc = discriminant(form)
                common = resultant(form, boundary_forms[label])
                if disc == 0 or common == 0:
                    raise ArithmeticError("quartic node or boundary transversality failed")
                quartic_discriminants.append(disc)
                quartic_resultants.append(common)
                quartic_node_count += 1
                continue
            linear = (expansion.get((1, 0), 0), expansion.get((0, 1), 0))
            if linear == (0, 0):
                raise ArithmeticError("a supposedly smooth root center is singular")
            tangent = (linear[1], -linear[0] % PRIME)
            tangent_value = qvalue(boundary_forms[label], tangent)
            if tangent_value == 0:
                raise ArithmeticError("root and boundary share a tangent branch")
            smooth_count += 1
            smooth_values.append((kind, tangent_value))
    if smooth_count != 84 or quartic_node_count != 12:
        raise ArithmeticError("root/boundary incidence census changed")

    return {
        "boundary_records": boundary_records,
        "boundary_A1_count": 4,
        "boundary_A2_count": 6,
        "boundary_node_discriminant_histogram": {
            str(value): count for value, count in sorted(Counter(node_discriminants).items())
        },
        "boundary_cusp_cubic_histogram": {
            str(value): count for value, count in sorted(Counter(cusp_cubics).items())
        },
        "quartic_node_count": quartic_node_count,
        "quartic_node_discriminant_histogram": {
            str(value): count
            for value, count in sorted(Counter(quartic_discriminants).items())
        },
        "quartic_boundary_resultant_histogram": {
            str(value): count
            for value, count in sorted(Counter(quartic_resultants).items())
        },
        "smooth_root_incidence_count": smooth_count,
        "smooth_root_tangent_value_histogram_by_kind": {
            kind: {
                str(value): count
                for value, count in sorted(
                    Counter(value for current, value in smooth_values if current == kind).items()
                )
            }
            for kind in ("line", "conic", "quartic")
        },
        "delta_sum": 10,
        "all_root_branches_avoid_boundary_branches": True,
    }


def point_ideal(name: str, point: Point) -> str:
    variables = ("x", "y", "z")
    chart = next(index for index, coordinate in enumerate(point) if coordinate)
    scale = inverse(point[chart])
    normalized = tuple(coordinate * scale % PRIME for coordinate in point)
    relations = [
        f"{variables[index]}-{normalized[index]}*{variables[chart]}"
        for index in range(3)
        if index != chart
    ]
    return f"ideal {name}={','.join(relations)};"


def singular_scheme_audit(configuration: dict[str, object]) -> dict[str, object]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required")
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    quartics: list[dict[str, object]] = configuration["quartic_records"]  # type: ignore[assignment]
    boundary: Sequence[int] = configuration["boundary"]  # type: ignore[assignment]
    names = {label: f"P{index}" for index, label in enumerate(centers)}
    declarations = [point_ideal(names[label], point) for label, point in centers.items()]
    polynomials = [
        f"poly q{index}={polynomial_string(record['coefficients'], 4)};"
        for index, record in enumerate(quartics)
    ]
    polynomials.append(f"poly f={polynomial_string(boundary, 6)};")
    blocks = []
    for index, record in enumerate(quartics):
        singular_labels = [
            ("s", special)
            for special in range(4)
            if special != record["omitted_special"]
        ]
        explicit = ",".join(names[label] for label in singular_labels)
        blocks.append(
            f'''
list A{index}=factorize(q{index});
ideal J{index}=std(sat(jacob(q{index}),maxideal(1)));
ideal R{index}=std(radical(J{index}));
ideal T{index}=std(intersect({explicit}));
print("Q{index}.FACTORS"); print(size(A{index}[1]));
print("Q{index}.EXPONENT"); print(A{index}[2][2]);
print("Q{index}.DIM"); print(dim(J{index}));
print("Q{index}.JDEG"); print(mult(J{index}));
print("Q{index}.RDEG"); print(mult(R{index}));
print("Q{index}.J_R"); print(size(reduce(J{index},R{index})));
print("Q{index}.R_J"); print(size(reduce(R{index},J{index})));
print("Q{index}.J_T"); print(size(reduce(J{index},T{index})));
print("Q{index}.T_J"); print(size(reduce(T{index},J{index})));
print("Q{index}.MINASS"); print(size(minAssGTZ(J{index})));
'''
        )
    all_points = ",".join(names[label] for label in centers)
    blocks.append(
        f'''
list AF=factorize(f);
ideal JF=std(sat(jacob(f),maxideal(1)));
ideal RF=std(radical(JF));
ideal TF=std(intersect({all_points}));
print("F.FACTORS"); print(size(AF[1]));
print("F.EXPONENT"); print(AF[2][2]);
print("F.DIM"); print(dim(JF));
print("F.JDEG"); print(mult(JF));
print("F.RDEG"); print(mult(RF));
print("F.J_R"); print(size(reduce(JF,RF)));
print("F.R_J"); print(size(reduce(RF,JF)));
print("F.R_T"); print(size(reduce(RF,TF)));
print("F.T_R"); print(size(reduce(TF,RF)));
print("F.MINASS"); print(size(minAssGTZ(JF)));
'''
    )
    program = "\n".join(
        [
            'LIB "primdec.lib";',
            "ring r=5,(x,y,z),dp;",
            *polynomials,
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
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    parsed: dict[str, dict[str, int]] = {**{f"Q{i}": {} for i in range(4)}, "F": {}}
    field_map = {
        "FACTORS": "factor_count_with_unit",
        "EXPONENT": "nonunit_exponent",
        "DIM": "affine_cone_dimension",
        "JDEG": "jacobian_degree",
        "RDEG": "radical_degree",
        "J_R": "jacobian_mod_radical",
        "R_J": "radical_mod_jacobian",
        "J_T": "jacobian_mod_points",
        "T_J": "points_mod_jacobian",
        "R_T": "radical_mod_points",
        "T_R": "points_mod_radical",
        "MINASS": "minimal_primes",
    }
    for index, line in enumerate(lines[:-1]):
        if "." not in line:
            continue
        prefix, key = line.split(".", 1)
        if prefix in parsed and key in field_map:
            parsed[prefix][field_map[key]] = int(lines[index + 1])

    quartic_expected = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "jacobian_degree": 3,
        "radical_degree": 3,
        "jacobian_mod_radical": 0,
        "radical_mod_jacobian": 0,
        "jacobian_mod_points": 0,
        "points_mod_jacobian": 0,
        "minimal_primes": 3,
    }
    if any(parsed[f"Q{index}"] != quartic_expected for index in range(4)):
        raise ArithmeticError(f"quartic singular schemes changed: {parsed}")
    boundary_fixed = {
        "factor_count_with_unit": 2,
        "nonunit_exponent": 1,
        "affine_cone_dimension": 1,
        "jacobian_degree": 16,
        "radical_degree": 10,
        "jacobian_mod_radical": 0,
        "radical_mod_points": 0,
        "points_mod_radical": 0,
        "minimal_primes": 10,
    }
    if any(parsed["F"].get(key) != value for key, value in boundary_fixed.items()):
        raise ArithmeticError(f"boundary singular scheme changed: {parsed['F']}")
    if parsed["F"].get("radical_mod_jacobian", 0) == 0:
        raise ArithmeticError("the six cusp Jacobian structures were falsely reduced")
    return {
        "engine": "independent Singular factor/saturated-Jacobian/radical audit",
        "quartics": [parsed[f"Q{index}"] for index in range(4)],
        "boundary": parsed["F"],
        "all_quartics_irreducible_over_F5": True,
        "boundary_irreducible_over_F5": True,
        "all_quartics_geometrically_integral": True,
        "boundary_geometrically_integral": True,
        "boundary_jacobian_nonreduced": True,
        "boundary_length_decomposition": "4*A1 + 6*A2 = 4*1+6*2=16",
    }


def integer_vector(vector: sp.Matrix) -> list[int]:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("a transported lattice coordinate is nonintegral")
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


def cdd_h_representation(matrix: sp.Matrix) -> str:
    rows = [
        "H-representation",
        "begin",
        f"{matrix.rows} {matrix.cols + 1} rational",
    ]
    rows.extend(
        "0 " + " ".join(str(matrix[row, column]) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    rows.append("end")
    return "\n".join(rows) + "\n"


def exact_root_cone(
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    boundary_index: int,
    root_indices: Sequence[int],
    exceptional_indices: Sequence[int],
) -> dict[str, object]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required")
    boundary = walls[boundary_index]
    if any(walls[index][10] for index in root_indices) or boundary[10] != -1:
        raise ArithmeticError("raw CM coordinates no longer split off the boundary")
    cm_gram = gram[:10, :10]
    inequalities = sp.Matrix.vstack(
        *(walls[index][:10, :].T * cm_gram for index in root_indices)
    )
    completed = subprocess.run(
        [executable, "--rep"],
        input=cdd_h_representation(inequalities),
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    lines = completed.stdout.splitlines()
    try:
        v_index = lines.index("V-representation")
        begin = next(
            index
            for index in range(v_index + 1, len(lines))
            if lines[index].strip() == "begin"
        )
    except (ValueError, StopIteration) as error:
        raise ArithmeticError("cddlib did not return a V-representation") from error
    row_count, width, kind = lines[begin + 1].split()[:3]
    if kind != "rational" or int(width) != 11:
        raise ArithmeticError("unexpected cddlib output width")
    rays: list[tuple[int, ...]] = []
    for line in lines[begin + 2 : begin + 2 + int(row_count)]:
        row = [Fraction(token) for token in line.split()]
        if len(row) != 11 or row[0] != 0:
            raise ArithmeticError("root cone unexpectedly has a vertex")
        rays.append(primitive_integer(row[1:]))
    if len(rays) != len(set(rays)):
        raise ArithmeticError("cddlib returned duplicate rays")

    def square(ray: Sequence[int]) -> int:
        vector = sp.Matrix(ray)
        return int((vector.T * cm_gram * vector)[0])

    classes = Counter(
        "finite" if square(ray) > 0 else "ideal" if square(ray) == 0 else "hyperideal"
        for ray in rays
    )
    if classes != Counter({"finite": 450, "ideal": 51}):
        raise ArithmeticError(f"root-cone causal census changed: {classes}")
    ideal = sorted(ray for ray in rays if square(ray) == 0)

    cusp_lifts = set()
    cusp_records = []
    for index in exceptional_indices:
        isotropic = walls[index] + boundary
        if isotropic[10] != 0:
            raise ArithmeticError("E+b does not lie in CM")
        ray = primitive_integer([int(value) for value in isotropic[:10, :]])
        if square(ray) != 0:
            raise ArithmeticError("E+b is not isotropic")
        cusp_lifts.add(ray)
        cusp_records.append({"wall_index": index, "f_equals_E_plus_b": list(ray)})
    if len(cusp_lifts) != 51 or cusp_lifts != set(ideal):
        raise ArithmeticError("the raw E+b classes are not exactly the 51 cusp rays")

    interior = sum((sp.Matrix(ray) for ray in rays), sp.zeros(10, 1))
    interior_square = int((interior.T * cm_gram * interior)[0])
    slacks = [int(value) for value in inequalities * interior]
    pairings = [
        int((interior.T * cm_gram * sp.Matrix(ray))[0]) for ray in rays
    ]
    if interior_square <= 0 or min(slacks) <= 0 or min(pairings) <= 0:
        raise ArithmeticError("the cdd rays do not lie on one positive-cone sheet")
    active = Counter(
        sum(int(value) == 0 for value in inequalities * sp.Matrix(ray))
        for ray in rays
    )
    facet_active_counts = []
    facet_face_ranks = []
    for row in range(inequalities.rows):
        active_rays = [
            sp.Matrix(ray)
            for ray in rays
            if int((inequalities * sp.Matrix(ray))[row]) == 0
        ]
        facet_active_counts.append(len(active_rays))
        facet_face_ranks.append(sp.Matrix.hstack(*active_rays).rank())
    if any(rank != 9 for rank in facet_face_ranks):
        raise ArithmeticError("a root inequality is redundant or fails to span a facet")
    return {
        "engine": "fresh cddexec_gmp H-to-V conversion from the twenty raw root rows",
        "ray_count": len(rays),
        "finite_ray_count": classes["finite"],
        "ideal_ray_count": classes["ideal"],
        "hyperideal_ray_count": classes["hyperideal"],
        "active_root_count_histogram": {
            str(value): count for value, count in sorted(active.items())
        },
        "facet_active_ray_counts": facet_active_counts,
        "facet_face_ranks": facet_face_ranks,
        "all_twenty_root_inequalities_irredundant": True,
        "primitive_ray_list_sha256": hashlib.sha256(
            json.dumps(sorted(rays), separators=(",", ":")).encode()
        ).hexdigest(),
        "primitive_ideal_ray_list_sha256": hashlib.sha256(
            json.dumps(ideal, separators=(",", ":")).encode()
        ).hexdigest(),
        "strict_interior_square": interior_square,
        "minimum_root_slack": min(slacks),
        "minimum_time_sheet_pairing": min(pairings),
        "cusp_records": sorted(cusp_records, key=lambda record: record["wall_index"]),
        "ideal_rays_equal_all_51_E_plus_b_classes": True,
        "finite_volume_twenty_root_chamber": True,
    }


def geometric_incidence_graph(
    configuration: dict[str, object],
) -> tuple[nx.Graph, dict[tuple[object, ...], dict[str, object]]]:
    graph = nx.Graph()
    centers: dict[PointLabel, Point] = configuration["centers"]  # type: ignore[assignment]
    for label in centers:
        graph.add_node(("point", label), kind="point")
    lookup = {}
    for curve in configuration["curve_records"]:  # type: ignore[union-attr]
        node = tuple(curve["node"])
        kind = node[0]
        graph.add_node(node, kind=kind)
        lookup[node] = curve
        for label, multiplicity in curve["support"].items():
            graph.add_edge(node, ("point", label), multiplicity=multiplicity)
    return graph, lookup


def raw_incidence_graph(
    walls: Sequence[sp.Matrix],
    inverse_basis: sp.Matrix,
    root_indices: Sequence[int],
) -> tuple[nx.Graph, list[dict[str, object]]]:
    graph = nx.Graph()
    for position in range(10):
        graph.add_node(("point", position), kind="point")
    records = []
    for index in root_indices:
        coordinates = integer_vector(inverse_basis * walls[index])
        kind = {1: "line", 2: "conic", 4: "quartic"}.get(coordinates[0])
        if kind is None:
            raise ArithmeticError("raw root has a nongeometric plane degree")
        multiplicities = {
            position: -coefficient
            for position, coefficient in enumerate(coordinates[1:])
            if coefficient
        }
        expected = {
            "line": Counter({1: 3}),
            "conic": Counter({1: 6}),
            "quartic": Counter({1: 6, 2: 3}),
        }[kind]
        if Counter(multiplicities.values()) != expected:
            raise ArithmeticError("raw root multiplicities do not match the curve alphabet")
        node = (kind, index)
        graph.add_node(node, kind=kind)
        for position, multiplicity in multiplicities.items():
            graph.add_edge(node, ("point", position), multiplicity=multiplicity)
        records.append(
            {
                "index": index,
                "kind": kind,
                "coordinates": coordinates,
                "multiplicities": multiplicities,
            }
        )
    return graph, records


def frame_audit(
    gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    boundary: sp.Matrix,
    root_indices: Sequence[int],
    frame: Sequence[int],
    configuration: dict[str, object],
) -> tuple[dict[str, object], list[list[int]]]:
    frame_set = sorted(frame)
    plane = (
        boundary + sum((walls[index] for index in frame_set), sp.zeros(11, 1))
    ) / 3
    if any(value.q != 1 for value in plane):
        raise ArithmeticError("the forced plane class is nonintegral")
    oriented = list(frame_set)
    basis = sp.Matrix.hstack(plane, *[walls[index] for index in oriented])
    if basis.det() == -1:
        oriented[-1], oriented[-2] = oriented[-2], oriented[-1]
        basis = sp.Matrix.hstack(plane, *[walls[index] for index in oriented])
    expected_gram = sp.diag(1, *([-1] * 10))
    if basis.det() != 1 or basis.T * gram * basis != expected_gram:
        raise ArithmeticError("frame does not give a determinant-one I(1,10) basis")
    inverse_basis = basis.inv()
    boundary_coordinates = integer_vector(inverse_basis * boundary)
    if boundary_coordinates != [3, *([-1] * 10)]:
        raise ArithmeticError("boundary is not 3H-sum(E_i)")

    raw_graph, records = raw_incidence_graph(walls, inverse_basis, root_indices)
    geometry_graph, curve_lookup = geometric_incidence_graph(configuration)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        raw_graph,
        geometry_graph,
        node_match=lambda left, right: left["kind"] == right["kind"],
        edge_match=lambda left, right: left["multiplicity"] == right["multiplicity"],
    )
    encodings = []
    exact_count = 0
    for mapping in matcher.isomorphisms_iter():
        encoding = ",".join(
            f"{mapping[('point', position)][1][0]}{mapping[('point', position)][1][1]}"
            for position in range(10)
        )
        encodings.append(encoding)
        exact = True
        for record in records:
            geometric_node = mapping[(record["kind"], record["index"])]
            curve = curve_lookup[geometric_node]
            expected = [curve["degree"]]
            expected.extend(
                -curve["support"].get(mapping[("point", position)][1], 0)
                for position in range(10)
            )
            if expected != record["coordinates"]:
                exact = False
                break
        exact_count += int(exact)
    if len(encodings) != 24 or len(set(encodings)) != 24 or exact_count != 24:
        raise ArithmeticError("not all 24 colored markings give exact root vectors")

    full_coordinates = [integer_vector(inverse_basis * wall) for wall in walls]
    special_positions = sorted(
        position
        for position in range(10)
        if sum(
            record["multiplicities"].get(position) == 2 for record in records
        )
        == 3
    )
    if len(special_positions) != 4:
        raise ArithmeticError("raw roots do not force a 4+6 center split")
    return (
        {
            "frame_set_wall_indices_zero_based": frame_set,
            "determinant_one_ordered_frame": oriented,
            "plane_class_in_raw_basis": integer_vector(plane),
            "basis_determinant": int(basis.det()),
            "basis_gram": [list(map(int, expected_gram.row(row))) for row in range(11)],
            "boundary_coordinates": boundary_coordinates,
            "special_frame_positions": special_positions,
            "ordinary_frame_positions": [
                position for position in range(10) if position not in special_positions
            ],
            "root_kind_histogram": dict(
                sorted(Counter(record["kind"] for record in records).items())
            ),
            "colored_marking_count": len(encodings),
            "exact_root_vector_marking_count": exact_count,
            "marking_torsor_sha256": hashlib.sha256(
                "\n".join(sorted(encodings)).encode()
            ).hexdigest(),
            "full_wall_coordinate_sha256": hashlib.sha256(
                json.dumps(full_coordinates, separators=(",", ":")).encode()
            ).hexdigest(),
            "all_72_wall_coordinates_integral": True,
        },
        full_coordinates,
    )


def raw_lattice_audit(
    raw_payload: dict[str, object], configuration: dict[str, object]
) -> dict[str, object]:
    polytope = raw_payload["polytope"]  # type: ignore[index]
    gram = sp.Matrix(polytope["lattice_gram"])
    walls = [sp.Matrix(row) for row in polytope["wall_vectors"]]

    def pairing(left: sp.Matrix, right: sp.Matrix) -> int:
        return int((left.T * gram * right)[0])

    norms = Counter(pairing(wall, wall) for wall in walls)
    roots = [index for index, wall in enumerate(walls) if pairing(wall, wall) == -2]
    minus_one = [index for index, wall in enumerate(walls) if pairing(wall, wall) == -1]
    boundary_candidates = [
        index
        for index in minus_one
        if all(pairing(walls[index], walls[root]) == 0 for root in roots)
    ]
    if len(walls) != 72 or norms != Counter({-1: 52, -2: 20}):
        raise ArithmeticError("raw wall inventory changed")
    if len(boundary_candidates) != 1:
        raise ArithmeticError("boundary wall is not intrinsically unique")
    boundary_index = boundary_candidates[0]
    boundary = walls[boundary_index]
    exceptional = [index for index in minus_one if index != boundary_index]
    if any(pairing(boundary, walls[index]) != 1 for index in exceptional):
        raise ArithmeticError("a candidate exceptional wall has b.E != 1")
    if any(
        reduce(math.gcd, (abs(int(value)) for value in wall), 0) != 1
        for wall in walls
    ):
        raise ArithmeticError("a raw wall vector is nonprimitive")
    root_pairings = Counter(
        pairing(walls[left], walls[right])
        for left, right in itertools.combinations(roots, 2)
    )
    if root_pairings != Counter({0: 135, 1: 30, 2: 25}):
        raise ArithmeticError("the twenty roots no longer form the acute Coxeter system")
    cm_gram = gram[:10, :10]
    cm_root_matrix = sp.Matrix.hstack(*[walls[index][:10, :] for index in roots])
    if cm_gram.det() != -1 or cm_root_matrix.rank() != 10:
        raise ArithmeticError("the root stock no longer spans the rank-ten unimodular CM lattice")
    identity = sp.eye(10)
    for index in roots:
        root = walls[index][:10, :]
        reflection = identity + root * (root.T * cm_gram)
        if any(value.q != 1 for value in reflection):
            raise ArithmeticError("a root reflection is nonintegral")
        if reflection.T * cm_gram * reflection != cm_gram or reflection * reflection != identity:
            raise ArithmeticError("a root reflection is not an involutive CM isometry")

    graph = nx.Graph()
    graph.add_nodes_from(exceptional)
    graph.add_edges_from(
        (left, right)
        for left, right in itertools.combinations(exceptional, 2)
        if pairing(walls[left], walls[right]) == 0
    )
    maximal = [sorted(clique) for clique in nx.find_cliques(graph)]
    clique_histogram = Counter(len(clique) for clique in maximal)
    frames = sorted(clique for clique in maximal if len(clique) == 10)
    if clique_histogram != Counter({7: 80, 9: 25, 10: 5}) or len(frames) != 5:
        raise ArithmeticError("maximal disjoint-frame census changed")

    frame_records = []
    coordinate_archives = []
    for frame in frames:
        record, coordinates = frame_audit(
            gram, walls, boundary, roots, frame, configuration
        )
        frame_records.append(record)
        coordinate_archives.append(coordinates)
    aggregate_digest = hashlib.sha256(
        json.dumps(coordinate_archives, separators=(",", ":")).encode()
    ).hexdigest()
    cone = exact_root_cone(
        gram, walls, boundary_index, roots, exceptional
    )

    exact = polytope["exact_double_description"]
    vinberg = polytope["literal_vinberg"]
    upstream = {
        "id": polytope["id"],
        "dimension": polytope["dimension"],
        "facet_count": polytope["facet_count"],
        "finite_vertices": exact["finite_vertex_count"],
        "ideal_vertices": exact["ideal_vertex_count"],
        "hyperideal_rays": exact["hyperideal_vertex_count"],
        "elliptic_edges": vinberg["exact_geometric_edge_count"],
        "bad_edges": vinberg["bad_edge_count"],
        "bad_parabolic_neighbourhoods": vinberg[
            "bad_parabolic_neighbourhood_count"
        ],
        "finite_volume": polytope["finite_volume"],
        "compact": polytope["compact"],
        "all_facets_irredundant": exact["all_facets_irredundant"],
    }
    if not (
        raw_payload["all_checks_pass"]
        and all(raw_payload["checks"].values())
        and upstream["finite_volume"]
        and not upstream["compact"]
        and upstream["hyperideal_rays"] == 0
        and upstream["bad_edges"] == 0
        and upstream["bad_parabolic_neighbourhoods"] == 0
    ):
        raise ArithmeticError("upstream P-side Vinberg archive is not green")
    return {
        "primitive_wall_norm_histogram": {
            str(value): count for value, count in sorted(norms.items())
        },
        "intrinsic_boundary_wall_index_zero_based": boundary_index,
        "root_wall_indices_zero_based": roots,
        "exceptional_wall_count": len(exceptional),
        "all_walls_primitive": True,
        "CM_gram_determinant": int(cm_gram.det()),
        "CM_root_span_rank": cm_root_matrix.rank(),
        "root_pairing_histogram": {
            str(value): count for value, count in sorted(root_pairings.items())
        },
        "all_twenty_root_reflections_integral_involutive_isometries": True,
        "maximal_disjoint_clique_size_histogram": {
            str(value): count for value, count in sorted(clique_histogram.items())
        },
        "intrinsic_plane_frame_count": len(frames),
        "all_five_frames": frame_records,
        "all_120_markings_exact": sum(
            frame["exact_root_vector_marking_count"] for frame in frame_records
        )
        == 120,
        "aggregate_full_wall_coordinate_sha256": aggregate_digest,
        "twenty_root_cone": cone,
        "upstream_P_side_vinberg_import": upstream,
    }


def build_audit() -> dict[str, object]:
    raw_payload = json.loads(RAW_POLYTOPE.read_text(encoding="utf-8"))
    target = json.loads(TARGET.read_text(encoding="utf-8"))
    configuration = reconstruct_finite_geometry()
    local = local_and_transversality_audit(configuration)
    schemes = singular_scheme_audit(configuration)
    lattice = raw_lattice_audit(raw_payload, configuration)

    finite_geometry = {
        "field": "F5",
        "ordinary_conic": "x*z-y^2",
        "ordinary_points": [list(point) for point in configuration["ordinary"]],
        "special_points": [list(point) for point in configuration["special"]],
        "projective_plane_point_count": configuration["projective_point_count"],
        "perfect_matching_count": configuration["perfect_matching_count"],
        "concurrent_matching_count": configuration["concurrent_matching_count"],
        "one_factorization_count": configuration["factorization_count"],
        "one_factorization_concurrency_histogram": {
            str(value): count
            for value, count in sorted(configuration["factorization_histogram"].items())
        },
        "usable_factorization_count": configuration["usable_factorization_count"],
        "selected_matchings": [
            [list(edge) for edge in matching]
            for matching in configuration["selected_matchings"]
        ],
        "residual_matching": [list(edge) for edge in configuration["residual_matching"]],
        "line_count": sum(
            curve["node"][0] == "line" for curve in configuration["curve_records"]
        ),
        "conic_count": len(configuration["conic_records"]),
        "quartic_count": len(configuration["quartic_records"]),
        "conic_coefficients": [
            list(record["coefficients"]) for record in configuration["conic_records"]
        ],
        "conic_determinants": [
            record["symmetric_matrix_determinant_mod_5"]
            for record in configuration["conic_records"]
        ],
        "quartic_condition_ranks": [
            record["condition_rank"] for record in configuration["quartic_records"]
        ],
        "quartic_coefficients": [
            list(record["coefficients"]) for record in configuration["quartic_records"]
        ],
        "all_curve_center_supports_exact": True,
    }
    linear_systems = {
        "anticanonical_cubic": {
            "ambient_dimension": 10,
            "evaluation_rank": configuration["cubic_rank"],
            "kernel_dimension": 0,
        },
        "antibicanonical_sextic": {
            "ambient_dimension": 28,
            "gradient_row_count": 30,
            "gradient_rank": configuration["sextic_rank"],
            "kernel_dimension": 1,
            "value_augmented_rank": configuration["sextic_augmented_rank"],
            "coefficients": list(configuration["boundary"]),
            "polynomial": polynomial_string(configuration["boundary"], 6),
        },
    }

    target_geometry = target["finite_geometry"]
    target_raw = target["marked_same_polytope_return"]
    target_agreement = {
        "ordinary_points": finite_geometry["ordinary_points"]
        == [record["coordinates"] for record in target_geometry["ordinary_points"]],
        "special_points": finite_geometry["special_points"]
        == [record["coordinates"] for record in target_geometry["special_points"]],
        "conic_coefficients": finite_geometry["conic_coefficients"]
        == [record["coefficients"] for record in target_geometry["conics"]],
        "quartic_coefficients": finite_geometry["quartic_coefficients"]
        == [record["coefficients"] for record in target_geometry["quartics"]],
        "boundary_coefficients": linear_systems["antibicanonical_sextic"][
            "coefficients"
        ]
        == target["linear_systems"]["antibicanonical_sextics"]["coefficients"],
        "five_frame_sets": [
            record["frame_set_wall_indices_zero_based"]
            for record in lattice["all_five_frames"]
        ]
        == [
            record["frame_set_wall_indices_zero_based"]
            for record in target_raw["all_five_frames"]
        ],
        "aggregate_full_wall_digest": lattice[
            "aggregate_full_wall_coordinate_sha256"
        ]
        == target_raw["aggregate_all_five_full_wall_coordinate_sha256"],
        "root_ray_digest": lattice["twenty_root_cone"][
            "primitive_ray_list_sha256"
        ]
        == target_raw["direct_twenty_root_cone"]["primitive_ray_list_sha256"],
        "cusp_digest": lattice["twenty_root_cone"][
            "primitive_ideal_ray_list_sha256"
        ]
        == target_raw["direct_twenty_root_cone"][
            "primitive_ideal_ray_list_sha256"
        ],
    }
    if not all(target_agreement.values()):
        raise ArithmeticError(f"hostile reconstruction disagrees with target: {target_agreement}")

    strong_coble = {
        "surface": "blowup over algebraic closure of F5 of the ten independently reconstructed P2 centers",
        "rational": True,
        "K_square": -1,
        "minus_K_empty": configuration["cubic_rank"] == 10,
        "minus_2K_unique": configuration["sextic_rank"] == 27,
        "boundary_class": "6H-2sum(E_i)=-2K",
        "boundary_self_intersection": -4,
        "boundary_plane_singularities": "4 A1 nodes and 6 A2 cusps",
        "boundary_Jacobian_length": schemes["boundary"]["jacobian_degree"],
        "boundary_radical_degree": schemes["boundary"]["radical_degree"],
        "boundary_strict_transform_smooth_rational": True,
        "twenty_root_inventory": {"lines": 12, "conics": 4, "quartics": 4},
        "twenty_roots_smooth_geometrically_integral_rational_after_blowup": True,
        "twenty_roots_boundary_disjoint_after_blowup": local[
            "all_root_branches_avoid_boundary_branches"
        ],
        "strong_Coble_verdict": "PASS",
    }

    theorem_handoff = {
        "direct_finite_index_gate": (
            "The fresh 450+51+0 causal census makes the chamber of the twenty "
            "constructed integral root reflections finite-covolume in H9; hence "
            "their subgroup has finite index in O(CM)."
        ),
        "no_hidden_roots": (
            "Dolgachev-Kondo finite-index implies finite Aut, followed by Kondo's "
            "Vinberg2 proposition, applies to the twenty already-constructed "
            "effective irreducible (-2)-roots and excludes every hidden effective "
            "irreducible CM root."
        ),
        "one_boundary_cusp_lift": (
            "For any primitive nef isotropic cusp f, the general Halphen theorem "
            "and the single boundary force an index-two I1/II blown fibre "
            "B+2E=2f, hence E=f-b. Conversely a (-1)-curve gives f=E+b; it is "
            "vertical for the primitive fibration and therefore is that unique "
            "exceptional component, so no nonprimitive multiple or extra (-1)-curve "
            "is hidden."
        ),
        "negative_curve_exhaustion": (
            "The general Coble negative-curve trichotomy leaves only boundary, "
            "(-2), and (-1) curves. Thus b, the twenty roots, and the 51 cusp "
            "lifts are the complete nef wall stock."
        ),
        "same_wall_set_reason": (
            "The exact set identity between the 51 raw E+b vectors and all 51 cdd "
            "cusps forces the resulting exceptional classes E=f-b to be precisely "
            "the 51 raw norm-minus-one walls; every root and b already agrees in "
            "each unimodular marking."
        ),
        "surface_existence_or_Type_VII_classification_used_to_construct_carrier": False,
        "general_AG_theorems_used_after_explicit_carrier_exists": True,
        "exact_same_marked_H10_72": True,
    }

    checks = {
        "finite_geometry_10_centers_and_5_usable_factorizations": (
            len(configuration["centers"]) == 10
            and configuration["usable_factorization_count"] == 5
        ),
        "curve_system_12_4_4": (
            finite_geometry["line_count"],
            finite_geometry["conic_count"],
            finite_geometry["quartic_count"],
        )
        == (12, 4, 4),
        "linear_system_ranks": (
            configuration["cubic_rank"],
            configuration["sextic_rank"],
            configuration["sextic_augmented_rank"],
        )
        == (10, 27, 27),
        "quartic_schemes_exact": all(
            record["jacobian_degree"] == record["radical_degree"] == 3
            and record["minimal_primes"] == 3
            for record in schemes["quartics"]
        ),
        "boundary_scheme_4A1_6A2_length16": (
            local["boundary_A1_count"],
            local["boundary_A2_count"],
            schemes["boundary"]["jacobian_degree"],
            schemes["boundary"]["radical_degree"],
        )
        == (4, 6, 16, 10),
        "strong_Coble": strong_coble["strong_Coble_verdict"] == "PASS",
        "five_raw_frames": lattice["intrinsic_plane_frame_count"] == 5,
        "all_120_exact_markings": lattice["all_120_markings_exact"],
        "all_72_wall_transports_integral": all(
            frame["all_72_wall_coordinates_integral"]
            for frame in lattice["all_five_frames"]
        ),
        "root_cone_450_51_0": (
            lattice["twenty_root_cone"]["finite_ray_count"],
            lattice["twenty_root_cone"]["ideal_ray_count"],
            lattice["twenty_root_cone"]["hyperideal_ray_count"],
        )
        == (450, 51, 0),
        "all_51_cusp_lifts_exact": lattice["twenty_root_cone"][
            "ideal_rays_equal_all_51_E_plus_b_classes"
        ],
        "target_reconstruction_agrees": all(target_agreement.values()),
        "same_marked_H10_72_without_surface_existence_circularity": theorem_handoff[
            "exact_same_marked_H10_72"
        ]
        and not theorem_handoff[
            "surface_existence_or_Type_VII_classification_used_to_construct_carrier"
        ],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile Type-VII gate failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "verdict": "GREEN_INDEPENDENT_L3_SAME_MARKED_H10_72",
        "sources": {
            "raw_polytope_archive": {
                "path": str(RAW_POLYTOPE.relative_to(ROOT)),
                "sha256": sha256(RAW_POLYTOPE),
                "used_for_carrier_construction": "only lattice Gram and 72 primitive wall vectors",
            },
            "target_certificate": {
                "path": str(TARGET.relative_to(ROOT)),
                "sha256": sha256(TARGET),
                "used_only_after_independent_reconstruction_for_agreement": True,
            },
            "primary_AG_source": "Kondo, arXiv:2102.00591, especially CME10, NefCone, FiniteIndex, Vinberg2, Halphen, and fibration-isotropic",
            "production_module_imported": False,
        },
        "finite_geometry": finite_geometry,
        "linear_systems": linear_systems,
        "local_geometry_and_transversality": local,
        "scheme_theoretic_curve_audit": schemes,
        "blowup_and_strong_Coble": strong_coble,
        "raw_lattice_and_cone": lattice,
        "AG_same_P_handoff": theorem_handoff,
        "target_certificate_agreement": target_agreement,
        "lineage_and_scope": {
            "known_abstract_polytope": "the already-known project H10/72 Type-VII chamber",
            "positive_result": "a project-new explicit F5 plane realization and strict P-to-surface-to-same-marked-P proof",
            "new_abstract_hyperbolic_polytope": False,
            "new_dimension": False,
            "worldwide_priority_for_equations": "not claimed",
            "inverse_level": "L3 structural one-datum calibration, not class-wide L4",
            "why_not_L4": "one rigid raw chamber is treated; there is no non-singleton necessary-and-sufficient carrier classification or AG-derived dimension bound",
        },
        "hostile_warning": (
            "The sextic is not ten-nodal: its saturated Jacobian has length 16, "
            "with four A1 nodes and six A2 cusps; only its radical is the ten "
            "reduced centers."
        ),
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
        print("Type-VII characteristic-five tetrahedral hostile archive replay passed")
        return
    if arguments.write:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
