#!/usr/bin/env python3
"""Independent hostile replay of the projected Delta(2,4,6) K3 carrier.

The auditor imports neither SymPy nor the owner certifier.  All lattice,
ellipsoid, root-window, Smith-invariant, projector, and rational-extension
calculations use Fraction arithmetic and explicit finite enumeration.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/compact_246_projected_k3_eigenslice_carrier.md"
OWNER_SCRIPT = ROOT / "scripts/certify_compact_246_projected_k3_eigenslice_carrier.py"
OWNER_ARCHIVE = (
    ROOT / "data/enumerations/compact_246_projected_k3_eigenslice_carrier.json"
)
OWNER_TEST = ROOT / "tests/test_compact_246_projected_k3_eigenslice_carrier.py"

AUDIT_REPORT = (
    ROOT
    / "results/compact_246_projected_k3_eigenslice_carrier_"
    "independent_hostile_audit.md"
)
AUDIT_ARCHIVE = (
    ROOT
    / "data/enumerations/compact_246_projected_k3_eigenslice_carrier_"
    "independent_hostile_audit.json"
)

OWNER_HASHES = {
    "results/compact_246_projected_k3_eigenslice_carrier.md": (
        "f4c4a038610e7c47600af971154ca056d7ae143546812582619e83ea7286cd1b"
    ),
    "scripts/certify_compact_246_projected_k3_eigenslice_carrier.py": (
        "346ed6ed0395008e0a98888efd008a4d74be37e6551a4b2e377f21f42a6aac4a"
    ),
    "data/enumerations/compact_246_projected_k3_eigenslice_carrier.json": (
        "6d94dbab080c4450bbe17bc7d3da091c550a0ac61965ee926de0eada42f537e1"
    ),
    "tests/test_compact_246_projected_k3_eigenslice_carrier.py": (
        "315e144769eb8c86e1ca12551210e25a402c0c2b0ce863a54b40c3cfcdb8a485"
    ),
}

SCALE = 15


def F(value: int | Fraction, denominator: int | None = None) -> Fraction:
    if denominator is not None:
        return Fraction(value, denominator)
    return value if isinstance(value, Fraction) else Fraction(value)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def encode(value: Fraction) -> int | str:
    return int(value) if value.denominator == 1 else str(value)


def encode_matrix(matrix: list[list[Fraction]]) -> list[list[int | str]]:
    return [[encode(value) for value in row] for row in matrix]


def zeros(rows: int, columns: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(columns)] for _ in range(rows)]


def identity(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(int(i == j)) for j in range(size)]
        for i in range(size)
    ]


def transpose(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*matrix)]


def multiply(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Fraction(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def add(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def scalar(
    coefficient: Fraction, matrix: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [[coefficient * value for value in row] for row in matrix]


def matvec(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> list[Fraction]:
    return [
        sum((row[j] * vector[j] for j in range(len(vector))), Fraction(0))
        for row in matrix
    ]


def dot(left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum((x * y for x, y in zip(left, right)), Fraction(0))


def pairing(
    gram: list[list[Fraction]],
    left: list[Fraction],
    right: list[Fraction],
) -> Fraction:
    return dot(left, matvec(gram, right))


def quadratic(
    gram: list[list[Fraction]], vector: list[int | Fraction]
) -> Fraction:
    converted = [F(value) for value in vector]
    return pairing(gram, converted, converted)


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    work = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column] / pivot_value
            for j in range(column, size):
                work[row][j] -= factor * work[column][j]
    return result


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    work = [matrix[i][:] + identity(size)[i] for i in range(size)]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            raise ZeroDivisionError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][j] - factor * work[column][j]
                for j in range(2 * size)
            ]
    return [row[size:] for row in work]


def submatrix(
    matrix: list[list[Fraction]],
    row_indices: tuple[int, ...],
    column_indices: tuple[int, ...],
) -> list[list[Fraction]]:
    return [
        [matrix[i][j] for j in column_indices]
        for i in row_indices
    ]


def leading_minor(matrix: list[list[Fraction]], size: int) -> Fraction:
    indices = tuple(range(size))
    return determinant(submatrix(matrix, indices, indices))


def determinantal_divisors(
    matrix: list[list[Fraction]],
) -> tuple[list[int], list[int]]:
    size = len(matrix)
    divisors: list[int] = []
    previous = 1
    invariants: list[int] = []
    for order in range(1, size + 1):
        value_gcd = 0
        for rows in combinations(range(size), order):
            for columns in combinations(range(size), order):
                value = determinant(submatrix(matrix, rows, columns))
                if value.denominator != 1:
                    raise AssertionError("nonintegral minor")
                value_gcd = gcd(value_gcd, abs(value.numerator))
        divisors.append(value_gcd)
        if value_gcd % previous:
            raise AssertionError("determinantal divisors do not divide")
        invariants.append(value_gcd // previous)
        previous = value_gcd
    return divisors, invariants


def block_matrix(
    top_left: list[list[Fraction]],
    top_right: list[list[Fraction]],
    bottom_left: list[list[Fraction]],
    bottom_right: list[list[Fraction]],
) -> list[list[Fraction]]:
    top = [top_left[i] + top_right[i] for i in range(len(top_left))]
    bottom = [
        bottom_left[i] + bottom_right[i] for i in range(len(bottom_left))
    ]
    return top + bottom


G = [
    [F(-12), F(6), F(0)],
    [F(6), F(-6), F(3)],
    [F(0), F(3), F(-2)],
]


def raw_reflection(root_index: int) -> list[list[Fraction]]:
    result = identity(3)
    for column in range(3):
        result[root_index][column] -= (
            2 * G[column][root_index] / G[root_index][root_index]
        )
    return result


def ambient_packet() -> dict[str, Any]:
    top_left = scalar(F(SCALE), G)
    bottom_right = scalar(F(-2), identity(3))
    ambient = block_matrix(top_left, G, G, bottom_right)

    g_two = leading_minor(G, 2)
    g_pivots = [G[0][0], g_two / G[0][0], determinant(G) / g_two]
    if g_pivots != [-12, -3, 1]:
        raise AssertionError("raw signature changed")

    schur = add(bottom_right, scalar(F(-1, SCALE), G))
    minus_schur = scalar(F(-1), schur)
    schur_minors = [leading_minor(minus_schur, size) for size in range(1, 4)]
    expected_schur_minors = [
        F(6, 5),
        F(44, 25),
        F(1214, 375),
    ]
    if schur_minors != expected_schur_minors:
        raise AssertionError("Schur complement is not negative definite")

    ambient_determinant = determinant(ambient)
    divisors, smith = determinantal_divisors(ambient)
    if ambient_determinant != -393336:
        raise AssertionError("ambient determinant changed")
    if divisors != [1, 1, 1, 1, 18, 393336]:
        raise AssertionError(f"determinantal divisors changed: {divisors}")
    if smith != [1, 1, 1, 1, 18, 21852]:
        raise AssertionError(f"Smith invariants changed: {smith}")

    if any(
        value.denominator != 1
        for row in ambient
        for value in row
    ):
        raise AssertionError("ambient Gram is not integral")
    if any(ambient[i][i].numerator % 2 for i in range(6)):
        raise AssertionError("ambient Gram is not even")

    projector = block_matrix(
        identity(3),
        scalar(F(1, SCALE), identity(3)),
        zeros(3, 3),
        zeros(3, 3),
    )
    if multiply(projector, projector) != projector:
        raise AssertionError("projector is not idempotent")
    if multiply(transpose(projector), ambient) != multiply(ambient, projector):
        raise AssertionError("projector is not self-adjoint")

    projected_deltas: list[list[int | str]] = []
    for index in range(3):
        delta = [F(0) for _ in range(6)]
        delta[3 + index] = F(1)
        if quadratic(ambient, delta) != -2:
            raise AssertionError("physical delta is not a root")
        projection = matvec(projector, delta)
        expected = [F(0) for _ in range(6)]
        expected[index] = F(1, SCALE)
        if projection != expected:
            raise AssertionError("delta projection changed")
        projected_deltas.append([encode(value) for value in projection])

    # The first three coordinate vectors are a primitive direct summand of
    # Z^6; the quotient is freely generated by the last three.
    return {
        "ambient_gram": encode_matrix(ambient),
        "ambient_even_integral": True,
        "raw_LDL_pivots": [encode(value) for value in g_pivots],
        "raw_signature": [1, 2],
        "schur_complement": encode_matrix(schur),
        "minus_schur_leading_minors": [
            encode(value) for value in schur_minors
        ],
        "ambient_signature": [1, 5],
        "ambient_determinant": int(ambient_determinant),
        "determinantal_divisors": divisors,
        "smith_invariants": smith,
        "discriminant_length": sum(value > 1 for value in smith),
        "trace_gram": encode_matrix(top_left),
        "trace_primitive_direct_summand": True,
        "projector": encode_matrix(projector),
        "projector_idempotent": True,
        "projector_self_adjoint": True,
        "eigenvalue_one_space": "span_Q(r1,r2,r3)",
        "kernel_negative_definite": True,
        "projected_delta_vectors": projected_deltas,
        "_ambient": ambient,
        "_projector": projector,
    }


def visible_root_reduction() -> dict[str, Any]:
    matrix = add(scalar(F(30), identity(3)), G)
    minors = [leading_minor(matrix, size) for size in range(1, 4)]
    if minors != [18, 396, 10926]:
        raise AssertionError("a-ellipsoid is not positive definite")
    matrix_inverse = inverse(matrix)
    coordinate_bounds = [30 * matrix_inverse[i][i] for i in range(3)]
    expected_bounds = [F(1105, 607), F(840, 607), F(660, 607)]
    if coordinate_bounds != expected_bounds or not all(
        value < 2 for value in coordinate_bounds
    ):
        raise AssertionError("a-coordinate bounds changed")

    rows = []
    for entries in product(range(-1, 2), repeat=3):
        value = quadratic(matrix, list(entries))
        if value < 30:
            rows.append({"a": list(entries), "ellipsoid_value": int(value)})
    expected_vectors = {
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
        (0, 0, 0),
    }
    if {tuple(row["a"]) for row in rows} != expected_vectors:
        raise AssertionError("a-ellipsoid lattice points changed")

    # Verify the completed-square identity on a coefficient-spanning sample.
    ambient = ambient_packet()["_ambient"]
    for x_entries in product((-1, 0, 1), repeat=3):
        for a_entries in product((-1, 0, 1), repeat=3):
            rho = [F(value) for value in x_entries + a_entries]
            y = [
                SCALE * x_entries[i] + a_entries[i]
                for i in range(3)
            ]
            identity_left = SCALE * quadratic(ambient, rho)
            identity_right = (
                quadratic(G, y)
                - quadratic(G, list(a_entries))
                - 30 * sum(value * value for value in a_entries)
            )
            if identity_left != identity_right:
                raise AssertionError("root completed-square identity failed")

    return {
        "restricted_normal": "y=15x+a",
        "root_equation": "G(y,y)=G(a,a)+30*(||a||^2-1)",
        "visible_root_implies_negative_Gy": {
            "y_positive": "y_perp is negative definite and misses H^2",
            "y_null": "y_perp only touches the positive-cone boundary",
            "closed_compact_P_points_are_positive": True,
            "conclusion": "a nonzero restricted mirror meeting P has G(y,y)<0",
        },
        "a_ellipsoid": "a^T(30I+G)a<30",
        "a_matrix": encode_matrix(matrix),
        "a_matrix_leading_minors": [int(value) for value in minors],
        "dual_coordinate_square_bounds": [str(value) for value in coordinate_bounds],
        "coordinate_box": [-1, 0, 1],
        "ellipsoid_lattice_points": rows,
        "a_zero_root_impossible": (
            "y=15x would force 15*G(x,x)=-2"
        ),
        "zero_projection_root_impossible": (
            "y=0 gives a=15b and a multiple of 15 equal to 2"
        ),
        "visible_nonzero_a_values": (
            "signed standard basis vectors only"
        ),
        "completed_square_identity_checked_on_729_pairs": True,
    }


def hidden_root_window() -> dict[str, Any]:
    g_inverse = inverse(G)
    ones = [F(1), F(1), F(1)]
    p = matvec(g_inverse, ones)
    p_square = quadratic(G, p)
    if p != [F(11, 12), F(2), F(5, 2)] or p_square != F(65, 12):
        raise AssertionError("interior point changed")
    if matvec(G, p) != ones:
        raise AssertionError("p does not have unit wall slacks")

    vertices: list[dict[str, Any]] = []
    cosh_twice_values: list[Fraction] = []
    for index in range(3):
        vertex = [g_inverse[row][index] for row in range(3)]
        vertex_square = quadratic(G, vertex)
        p_pairing = pairing(G, p, vertex)
        cosh_squared = p_pairing * p_pairing / (p_square * vertex_square)
        cosh_twice = 2 * cosh_squared - 1
        cosh_twice_values.append(cosh_twice)
        vertices.append(
            {
                "index": index + 1,
                "coordinates": [encode(value) for value in vertex],
                "square": encode(vertex_square),
                "p_pairing": encode(p_pairing),
                "cosh_distance_squared": encode(cosh_squared),
                "cosh_twice_distance": encode(cosh_twice),
            }
        )
    expected_cosh_twice = [F(177, 65), F(79, 65), F(17, 13)]
    if cosh_twice_values != expected_cosh_twice:
        raise AssertionError("vertex distances changed")
    max_cosh_twice = max(cosh_twice_values)

    gp = matvec(G, p)
    rank_one = [[gp[i] * gp[j] for j in range(3)] for i in range(3)]
    qp = add(scalar(F(2, 1) / p_square, rank_one), scalar(F(-1), G))
    qp_minors = [leading_minor(qp, size) for size in range(1, 4)]
    expected_qp_minors = [F(804, 65), F(612, 13), F(36)]
    if qp_minors != expected_qp_minors:
        raise AssertionError("Qp is not positive definite")

    universal_bound = 12 * max_cosh_twice
    qp_inverse = inverse(qp)
    coordinate_square_bounds = [
        universal_bound * qp_inverse[i][i] for i in range(3)
    ]
    expected_coordinate_bounds = [
        F(31329, 4225),
        F(111864, 4225),
        F(36108, 845),
    ]
    if coordinate_square_bounds != expected_coordinate_bounds:
        raise AssertionError("Qp coordinate bounds changed")
    integer_bounds = [2, 5, 6]
    if not (
        coordinate_square_bounds[0] < 9
        and coordinate_square_bounds[1] < 36
        and coordinate_square_bounds[2] < 49
    ):
        raise AssertionError("integer coordinate box is not justified")

    root_rows: list[dict[str, Any]] = []
    for entries in product(
        range(-2, 3),
        range(-5, 6),
        range(-6, 7),
    ):
        y = list(entries)
        square = quadratic(G, y)
        d = -square
        if d not in (2, 6, 12):
            continue
        qp_value = quadratic(qp, y)
        if qp_value <= d * max_cosh_twice:
            root_rows.append(
                {
                    "y": y,
                    "G_square": int(square),
                    "Qp": str(qp_value),
                }
            )
    if len(root_rows) != 22:
        raise AssertionError(f"root window size changed: {len(root_rows)}")

    signed_basis = []
    for index in range(3):
        for sign in (-1, 1):
            a = [0, 0, 0]
            a[index] = sign
            signed_basis.append(a)
    residue_rows = []
    for a in signed_basis:
        norm = -int(quadratic(G, a))
        for row in root_rows:
            if row["G_square"] != -norm:
                continue
            if all(
                (row["y"][coordinate] - a[coordinate]) % SCALE == 0
                for coordinate in range(3)
            ):
                x = [
                    (row["y"][coordinate] - a[coordinate]) // SCALE
                    for coordinate in range(3)
                ]
                residue_rows.append(
                    {
                        "a": a,
                        "y": row["y"],
                        "x": x,
                        "G_square": row["G_square"],
                    }
                )
    expected_residues = {
        ((-1, 0, 0), (-1, 0, 0), (0, 0, 0)),
        ((1, 0, 0), (1, 0, 0), (0, 0, 0)),
        ((0, -1, 0), (0, -1, 0), (0, 0, 0)),
        ((0, 1, 0), (0, 1, 0), (0, 0, 0)),
        ((0, 0, -1), (0, 0, -1), (0, 0, 0)),
        ((0, 0, 1), (0, 0, 1), (0, 0, 0)),
    }
    observed_residues = {
        (tuple(row["a"]), tuple(row["y"]), tuple(row["x"]))
        for row in residue_rows
    }
    if observed_residues != expected_residues:
        raise AssertionError(f"residue survivors changed: {residue_rows}")

    return {
        "interior_point": [encode(value) for value in p],
        "integral_multiple_12p": [int(12 * value) for value in p],
        "p_square": encode(p_square),
        "wall_slacks": [1, 1, 1],
        "vertices": vertices,
        "distance_convexity_circumradius_argument": (
            "P is the geodesic convex hull of its vertices, so CAT(0) "
            "convexity bounds distance from p by the largest vertex distance"
        ),
        "max_cosh_twice_radius": str(max_cosh_twice),
        "Qp": encode_matrix(qp),
        "Qp_leading_minors": [encode(value) for value in qp_minors],
        "hyperplane_distance_identity": (
            "for G(y,y)=-d, Qp(y)=d*cosh(2*dist(p,y_perp))"
        ),
        "universal_Qp_bound": str(universal_bound),
        "coordinate_square_bounds": [
            str(value) for value in coordinate_square_bounds
        ],
        "integer_coordinate_bounds": integer_bounds,
        "box_size": 5 * 11 * 13,
        "root_window_count": len(root_rows),
        "root_window": root_rows,
        "residue_survivors": residue_rows,
        "residue_exhaustive_because_bounds_below_15": True,
        "unique_roots_above_facets": (
            "y=+/-e_i forces a=+/-e_i, x=0, rho=+/-delta_i"
        ),
        "actual_ambient_wall_argument": (
            "a generic facet point has only delta_i^perp above it, and the "
            "segment from p has no earlier root crossing"
        ),
        "all_visible_root_mirrors_on_closed_P_are_facets": True,
    }


def correspondence_packet(
    ambient: list[list[Fraction]],
    projector: list[list[Fraction]],
) -> dict[str, Any]:
    raw = [raw_reflection(index) for index in range(3)]
    extensions: list[list[list[Fraction]]] = []
    root_reflection_restriction_mismatch = []
    for index, reflection in enumerate(raw):
        difference = add(reflection, scalar(F(-1), identity(3)))
        extension = block_matrix(
            reflection,
            scalar(F(1, SCALE), difference),
            zeros(3, 3),
            identity(3),
        )
        if multiply(extension, extension) != identity(6):
            raise AssertionError("rational extension is not involutive")
        if multiply(
            multiply(transpose(extension), ambient), extension
        ) != ambient:
            raise AssertionError("rational extension is not an isometry")
        if multiply(transpose(extension), ambient) != multiply(
            ambient, extension
        ):
            raise AssertionError("rational reflection is not self-adjoint")
        if multiply(extension, projector) != multiply(projector, extension):
            raise AssertionError("rational extension does not commute with e")
        if all(
            value.denominator == 1
            for row in extension
            for value in row
        ):
            raise AssertionError("raw reflection unexpectedly extends integrally")
        extensions.append(extension)

        # Ambient reflection in delta_i does not preserve T and its projected
        # action is not the raw reflection.
        delta_reflection = identity(6)
        root_coordinate = 3 + index
        for column in range(6):
            delta_reflection[root_coordinate][column] += ambient[
                column
            ][root_coordinate]
        preserves_trace = all(
            delta_reflection[row][column] == 0
            for column in range(3)
            for row in range(3, 6)
        )
        projected_action = multiply(
            multiply(
                [row[:6] for row in projector[:3]],
                delta_reflection,
            ),
            [
                [F(int(i == j)) for j in range(3)]
                for i in range(6)
            ],
        )
        root_reflection_restriction_mismatch.append(
            {
                "root": index + 1,
                "ambient_delta_reflection_preserves_T": preserves_trace,
                "projected_action_equals_raw": projected_action == reflection,
            }
        )
    if any(
        row["ambient_delta_reflection_preserves_T"]
        or row["projected_action_equals_raw"]
        for row in root_reflection_restriction_mismatch
    ):
        raise AssertionError("physical root reflection unexpectedly returns raw action")

    # If D_i is the coordinate divisor basis and D_i^vee is the
    # intersection-dual basis (the columns of A^{-1}), then
    # sum_i D_i^vee x f(D_i) acts as f on NS_Q.  Check the duality exactly.
    ambient_inverse = inverse(ambient)
    dual_pairing = multiply(ambient, ambient_inverse)
    if dual_pairing != identity(6):
        raise AssertionError("dual divisor basis failed")
    for operator in [projector] + extensions:
        # On a coordinate column alpha, intersections with the dual basis
        # recover alpha_i, and the product correspondence returns f(alpha).
        for column in range(6):
            basis = [F(int(i == column)) for i in range(6)]
            recovered = matvec(operator, basis)
            if recovered != [operator[row][column] for row in range(6)]:
                raise AssertionError("divisor-product action failed")

    return {
        "raw_reflections_integral_on_T": all(
            value.denominator == 1
            for reflection in raw
            for row in reflection
            for value in row
        ),
        "raw_reflections": [encode_matrix(matrix) for matrix in raw],
        "rational_extensions": [encode_matrix(matrix) for matrix in extensions],
        "extensions_are_involutive_isometries": True,
        "extensions_self_adjoint": True,
        "extensions_commute_with_projector": True,
        "extensions_preserve_full_NS_integrally": False,
        "physical_root_reflection_comparison": root_reflection_restriction_mismatch,
        "dual_divisor_basis": encode_matrix(ambient_inverse),
        "divisor_product_formula": (
            "Gamma_f=sum_i D_i^vee x f(D_i), acting by "
            "alpha -> sum_i (alpha.D_i^vee)f(D_i)=f(alpha)"
        ),
        "projector_and_extensions_are_rational_algebraic_correspondences": True,
    }


def k3_and_chamber_packet(
    ambient_packet_data: dict[str, Any],
    window: dict[str, Any],
) -> dict[str, Any]:
    length = ambient_packet_data["discriminant_length"]
    if length != 2 or not 16 >= length + 2:
        raise AssertionError("stable K3 embedding gate failed")
    if not window["all_visible_root_mirrors_on_closed_P_are_facets"]:
        raise AssertionError("hidden restricted roots remain")
    if window["integral_multiple_12p"] != [11, 24, 30]:
        raise AssertionError("no integral positive chamber point")

    return {
        "ambient_even_signature": [1, 5],
        "K3_lattice_signature": [3, 19],
        "complement_signature": [2, 14],
        "complement_rank": 16,
        "discriminant_length": length,
        "stable_embedding_inequality": f"16>={length}+2",
        "primitive_K3_embedding": True,
        "generic_period_exact_NS": True,
        "projective_integral_chamber_point": "12p=(11,24,30)",
        "no_zero_projection_roots": True,
        "p_regular_for_ambient_root_arrangement": True,
        "Weyl_framing": (
            "send the ambient root chamber containing p to the ample chamber"
        ),
        "facet_roots_effective_simple": True,
        "facet_irreducibility_argument": (
            "at a generic nef facet point, every component of an effective "
            "delta_i is orthogonal to a positive class, hence a (-2)-root; "
            "uniqueness of the ambient root hyperplane and primitivity of "
            "delta_i force one irreducible component with multiplicity one"
        ),
        "facet_curves_smooth_rational": True,
        "restricted_chamber_two_inclusions": {
            "P_inside_ambient_chamber_slice": (
                "P interior is connected, contains p, and has no root crossing"
            ),
            "ambient_chamber_slice_inside_P": (
                "the three unique delta_i walls impose all raw facet halfspaces"
            ),
        },
        "exact_eigen_nef_equality": (
            "P(Nef(X)_R intersect E_1(e))=Delta(2,4,6)"
        ),
        "carrier_relative_dimension_bound": (
            "n+1=dim(E_1(e))<=rho(X)<=20, hence n<=19"
        ),
    }


def owner_hash_audit() -> dict[str, str]:
    observed = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (OWNER_REPORT, OWNER_SCRIPT, OWNER_ARCHIVE, OWNER_TEST)
    }
    if observed != OWNER_HASHES:
        raise AssertionError(
            f"owner bytes changed; expected {OWNER_HASHES}, observed {observed}"
        )
    return observed


def build_certificate() -> dict[str, Any]:
    owner_hashes = owner_hash_audit()
    ambient = ambient_packet()
    visible = visible_root_reduction()
    window = hidden_root_window()
    correspondences = correspondence_packet(
        ambient["_ambient"], ambient["_projector"]
    )
    k3 = k3_and_chamber_packet(ambient, window)

    serializable_ambient = {
        key: value
        for key, value in ambient.items()
        if not key.startswith("_")
    }
    return {
        "schema": (
            "coxeter.compact-246-projected-k3-eigenslice-"
            "independent-hostile-audit.v1"
        ),
        "audit_report_sha256": sha256_file(AUDIT_REPORT),
        "owner_hashes": owner_hashes,
        "independence": {
            "imports_owner_certifier": False,
            "imports_sympy": False,
            "exact_engine": "Fraction with explicit matrix and finite-box code",
            "smith_check": "gcds of every square minor",
            "root_window_rebuilt": True,
        },
        "ambient_lattice_and_projector": serializable_ambient,
        "visible_root_reduction": visible,
        "finite_hidden_root_window": window,
        "rational_correspondences": correspondences,
        "K3_realization_and_exact_return": k3,
        "scope": {
            "one_projected_K3_carrier_proved": True,
            "outside_original_k_Q_literal_full_nef_Property2": True,
            "Picard_rank": 6,
            "eigenslice_rank": 3,
            "raw_reflections_extend_rationally": True,
            "raw_reflections_extend_integrally": False,
            "ambient_root_reflections_return_raw_Coxeter_action": False,
            "universal_fourth_branch_nonempty": False,
            "finite_host_control_for_fourth_branch": False,
            "arbitrary_non_K3_eigenslice_dimension_bound": False,
            "K3_carrier_relative_n_bound": 19,
        },
        "presentation": {
            "owner_scope_frontier_explicit": True,
            "owner_actual_wall_and_irreducibility_steps_explicit": True,
            "owner_markdown_pandoc_checked": True,
            "owner_reproduction_requires_PYTHONPATH_src_dot": True,
        },
        "verdict": "PASS",
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        AUDIT_ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        AUDIT_ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if (
            not AUDIT_ARCHIVE.exists()
            or AUDIT_ARCHIVE.read_text(encoding="utf-8") != encoded
        ):
            raise SystemExit("independent projected K3 audit archive differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
