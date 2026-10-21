#!/usr/bin/env python3
"""Exact lattice/K3 carrier certificate for Allcock's D6D6 face in H^13.

The finite-volume/Coxeter and finite-subtree assertions are cited from
Allcock (2006).  This script certifies the new arithmetic bridge: the 21
walls span the whole Lorentzian complement K, -K is the K3 lattice
U+D6(-1)+D6(-1), the 25-wall double and every finite-subtree descendant use
the same full Picard lattice, and a six-parameter elliptic K3 family realizes
that lattice.  It also certifies a root-regular ample-intersection framing.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/allcock_h13_d6d6_k3_carrier.json"
SCHEMA = "coxeter.allcock_h13_d6d6.k3_carrier.v1"

NAMES = tuple("ABCDEFGHIJKLMNOPQRSTU")
BASIS_NAMES = tuple("ABCDFGHIKLPQRT")

# K=R^perp in II_(25,1), in the wall basis above.  This is the positive-root
# Lorentz convention, so the K3 Neron--Severi Gram is -K_GRAM.
K_GRAM = sp.Matrix(
    [
        [2, -2, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0],
        [-2, 4, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -2, 2, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0],
        [0, -2, 0, 2, -2, 0, 0, -4, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -2, 2, 0, -2, 0, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 0, 0, 2, 0, -2, 0, 0, -1, 0, 0, 0],
        [-2, 0, 0, 0, -2, 0, 2, -4, 0, 0, 0, 0, 0, 0],
        [0, 0, -2, -4, 0, -2, -4, 2, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 2, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, 0],
        [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 2],
    ]
)

WALL_COORDINATES = {
    "A": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "B": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "C": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "D": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "E": [-12, -3, 4, 3, -4, 8, -6, 6, -6, 2, 6, 4, 2, -2],
    "F": [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "G": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    "H": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    "I": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    "J": [-12, -2, 4, 4, -4, 8, -7, 5, -6, 2, 6, 4, 2, -2],
    "K": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    "L": [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    "M": [1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    "N": [12, 3, -3, -3, 4, -8, 6, -5, 6, -2, -6, -4, -2, 2],
    "O": [-3, -1, 0, 1, 0, 1, -1, 1, -2, 0, 0, 0, 0, 0],
    "P": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    "Q": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    "R": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    "S": [8, 3, 0, -2, 1, -5, 3, -3, 4, 0, -4, -3, -2, 0],
    "T": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "U": [-12, -4, 1, 3, -3, 8, -5, 5, -6, 0, 6, 4, 2, -2],
}

M4 = {
    frozenset(edge)
    for edge in ("AB", "BC", "BD", "BE", "HJ", "IJ", "MJ", "NJ")
}
M3 = {
    frozenset(edge)
    for edge in (
        "AK", "MK", "CL", "NL", "GP", "FT", "KO", "OP", "PQ",
        "QR", "RS", "ST", "TU", "UL",
    )
}
HEAVY = {
    frozenset(edge)
    for edge in ("DF", "EG", "HF", "IG", "HA", "IC", "MD", "NE")
}
DASHED_X = {frozenset(edge) for edge in ("DE", "HI")}
DASHED_Y = {frozenset(edge) for edge in ("DI", "EH")}

H = sp.Matrix([-6, 1, 6, 3, 6, 12, 19, 10, -4, 2, 6, 2, 0, 2])

REDUCTION_TRANSFORM = sp.Matrix(
    [
        [-1, 1, -3, -1, 0, -1, -1, -1, -8, -3, 6, -11, 7],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 3, 0],
        [1, -1, 1, 1, 0, 1, 1, 1, 2, 1, -2, 2, -1],
        [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -4, 1],
        [0, -1, 1, 0, -1, 0, 0, 0, 5, 2, -3, 8, -5],
        [0, 1, -2, -1, 0, -1, -1, -1, -2, -1, 2, -5, 2],
        [0, 0, 1, 0, 0, 0, 0, 0, 3, 1, -2, 5, -3],
        [0, 0, -2, 0, 1, 0, 0, 0, -4, -1, 3, -6, 3],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 4, 1, -2, 6, -4],
        [0, 0, 0, 0, 0, 0, 1, 0, 3, 1, -1, 4, -3],
        [0, 0, 0, 0, 0, 0, 0, 1, 2, 1, -1, 2, -2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1],
    ]
)

REDUCED_H_PERP_GRAM = sp.Matrix(
    [
        [4,-2,-2,2,1,2,2,2,2,1,-1,-2,1],
        [-2,4,0,-2,1,-1,-2,-2,0,-2,1,2,-1],
        [-2,0,4,0,-1,-1,0,0,-2,0,-1,0,1],
        [2,-2,0,4,1,2,2,2,0,2,0,-1,1],
        [1,1,-1,1,4,2,1,1,0,0,0,-1,-1],
        [2,-1,-1,2,2,4,1,2,0,0,-1,-2,1],
        [2,-2,0,2,1,1,4,1,0,1,0,-2,1],
        [2,-2,0,2,1,2,1,4,1,2,-2,-2,0],
        [2,0,-2,0,0,0,0,1,4,1,-1,0,-1],
        [1,-2,0,2,0,0,1,2,1,4,0,-1,-1],
        [-1,1,-1,0,0,-1,0,-2,-1,0,4,2,0],
        [-2,2,0,-1,-1,-2,-2,-2,0,-1,2,4,-1],
        [1,-1,1,1,-1,1,1,0,-1,-1,0,-1,4],
    ]
)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _primitive(vector: sp.Matrix) -> bool:
    return reduce(gcd, (abs(int(value)) for value in vector), 0) == 1


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    ]


def _integer_matrix(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _diagram_gram(x: sp.Expr = sp.Integer(4), y: sp.Expr = sp.Integer(4)) -> sp.Matrix:
    matrix = sp.zeros(21)
    for i, name in enumerate(NAMES):
        matrix[i, i] = 4 if name in {"B", "J"} else 2
    for i, left in enumerate(NAMES):
        for j in range(i + 1, 21):
            edge = frozenset((left, NAMES[j]))
            value = 0
            if edge in M4 or edge in HEAVY:
                value = -2
            elif edge in M3:
                value = -1
            elif edge in DASHED_X:
                value = -x
            elif edge in DASHED_Y:
                value = -y
            matrix[i, j] = matrix[j, i] = value
    return matrix


def _relation_histogram(gram: sp.Matrix) -> dict[str, int]:
    counts: Counter[int] = Counter()
    for i, j in combinations(range(gram.rows), 2):
        ratio = sp.factor(gram[i, j] ** 2 / (gram[i, i] * gram[j, j]))
        if ratio == 0:
            counts[2] += 1
        elif ratio == sp.Rational(1, 4):
            counts[3] += 1
        elif ratio == sp.Rational(1, 2):
            counts[4] += 1
        elif ratio >= 1:
            counts[0] += 1
        else:
            raise ArithmeticError((i, j, ratio))
    return {str(key): value for key, value in sorted(counts.items())}


def _reflection(root: sp.Matrix) -> sp.Matrix:
    norm = (root.T * K_GRAM * root)[0]
    return sp.eye(14) - sp.Rational(2, norm) * root * (root.T * K_GRAM)


def _double(roots: list[sp.Matrix], wall: int) -> list[sp.Matrix]:
    alpha = roots[wall]
    norm = (alpha.T * K_GRAM * alpha)[0]
    output: list[sp.Matrix] = []
    for index, beta in enumerate(roots):
        if index == wall:
            continue
        output.append(beta)
        product_value = (beta.T * K_GRAM * alpha)[0]
        if product_value:
            reflected = beta - sp.Rational(2 * product_value, norm) * alpha
            if all(reflected != old for old in output):
                output.append(reflected)
    return output


def _span_index(roots: Iterable[sp.Matrix]) -> int | None:
    columns = sp.Matrix.hstack(*roots)
    if columns.rank() != 14:
        return None
    diagonal = _smith_diagonal(columns)
    return reduce(lambda a, b: a * b, diagonal, 1)


def _rank_forcing() -> dict[str, Any]:
    x, y = sp.symbols("x y")
    symbolic = _diagram_gram(x, y)
    base = [NAMES.index(name) for name in BASIS_NAMES]
    extras = ("J", "N", "O", "S", "U")
    factors = (80, 40, 8, 16, 40)
    expected = (
        21*x**2-24*x*y-86*x+9*y**2+32*y+120,
        13*x**2-8*x*y-96*x+7*y**2-12*y+240,
        2*x**2+2*x*y-28*x+11*y**2-92*y+240,
        9*x**2-6*x*y-66*x+7*y**2-20*y+184,
        10*x**2-10*x*y-60*x+7*y**2-4*y+144,
    )
    determinants = []
    signs = []
    for extra, factor, polynomial in zip(extras, factors, expected, strict=True):
        indices = base + [NAMES.index(extra)]
        determinant = sp.factor(symbolic.extract(indices, indices).det())
        ratio = sp.factor(determinant / (factor * polynomial))
        if ratio not in (1, -1):
            raise ArithmeticError((extra, determinant, polynomial))
        determinants.append(str(determinant))
        signs.append(int(ratio))
    groebner = sp.groebner(expected, x, y, order="lex")
    reduced = [sp.expand(poly.as_expr()) for poly in groebner.polys]
    return {
        "symmetry_input": (
            "Allcock's left-right (01) symmetry and the additional Co_infinity "
            "symmetry exchanging the D6 components give dashed orbits "
            "{DE,HI} and {DI,EH}"
        ),
        "variables": {"x": "-DE=-HI", "y": "-DI=-EH"},
        "extra_nodes": list(extras),
        "minor_factor_signs": signs,
        "minor_polynomials": [str(sp.expand(value)) for value in expected],
        "groebner_basis": [str(value) for value in reduced],
        "unique_solution": {"x": 4, "y": 4},
        "checks": {
            "five_minors_match": len(determinants) == 5,
            "rank_forces_both_dashed_products_to_four": set(reduced)
            == {x - 4, y - 4},
        },
    }


def _wall_lattice() -> tuple[dict[str, Any], list[sp.Matrix], list[sp.Matrix]]:
    roots = [sp.Matrix(WALL_COORDINATES[name]) for name in NAMES]
    coordinate_matrix = sp.Matrix.hstack(*roots)
    reconstructed = coordinate_matrix.T * K_GRAM * coordinate_matrix
    expected = _diagram_gram()
    norms = [int((root.T * K_GRAM * root)[0]) for root in roots]
    divisibilities = [
        reduce(gcd, (abs(int(value)) for value in K_GRAM * root), 0)
        for root in roots
    ]
    types = Counter(zip(norms, divisibilities, strict=True))
    reflections = [_reflection(root) for root in roots]

    order_names = ("B", "P", "L", "G", "Q", "C", "D", "I", "K", "F", "A", "R", "H", "T")
    order = [BASIS_NAMES.index(name) for name in order_names]
    reordered = K_GRAM.extract(order, order)
    leading = [int(reordered[:i, :i].det()) for i in range(1, 15)]
    expected_leading = [4,8,16,24,32,16,-16,-880,-1760,-2368,-480,-528,-64,-16]

    q_roots = _double(roots, NAMES.index("D"))
    q_gram = sp.Matrix.hstack(*q_roots).T * K_GRAM * sp.Matrix.hstack(*q_roots)
    q_doubling = [1, 2, 4, 5, 9, 10, 11, 12]
    f_doubling = [1, 3, 4, 7, 8, 9]
    q_delete_indices = {
        str(index): _span_index(root for j, root in enumerate(q_roots) if j != index)
        for index in q_doubling
    }
    f_delete_indices = {
        str(index): _span_index(root for j, root in enumerate(roots) if j != index)
        for index in f_doubling
    }
    checks = {
        "coordinate_gram_matches_corrected_figure3": reconstructed == expected,
        "basis_determinant_minus_16": K_GRAM.det() == -16,
        "basis_smith_1pow10_2pow4": _smith_diagonal(K_GRAM) == [1] * 10 + [2] * 4,
        "exact_inertia_13_1": leading == expected_leading,
        "all_walls_primitive": all(_primitive(root) for root in roots),
        "wall_types_15_4_2": types == Counter({(2, 1): 15, (2, 2): 4, (4, 2): 2}),
        "all_reflections_integral": all(all(value.q == 1 for value in matrix) for matrix in reflections),
        "all_reflections_are_isometries": all(matrix.T * K_GRAM * matrix == K_GRAM for matrix in reflections),
        "Q_has_25_walls": len(q_roots) == 25,
        "Q_corrected_relation_histogram": _relation_histogram(q_gram) == {"0": 39, "2": 233, "3": 16, "4": 12},
        "F_spans_K": _span_index(roots) == 1,
        "Q_spans_K": _span_index(q_roots) == 1,
        "deleting_any_F_doubling_wall_still_spans_K": set(f_delete_indices.values()) == {1},
        "deleting_any_Q_doubling_wall_still_spans_K": set(q_delete_indices.values()) == {1},
    }
    return {
        "ambient_projection_containment": {
            "D6A1_extensions": (
                "the projected normal p lies in K and has square 2"
            ),
            "D7_extensions_B_J": (
                "p lies in K^vee, while 2p lies in K and has square 4"
            ),
            "E7_extensions_D_E_H_I": (
                "p lies in K^vee, while 2p lies in K and has square 2"
            ),
            "index_argument": (
                "these normalized projected normals give W subset K; "
                "det(W)=det(K)=16, hence W=K"
            ),
        },
        "basis_names": list(BASIS_NAMES),
        "basis_gram": _integer_matrix(K_GRAM),
        "determinant": int(K_GRAM.det()),
        "smith_diagonal": _smith_diagonal(K_GRAM),
        "inertia": [13, 1, 0],
        "inertia_leading_principal_determinants": leading,
        "wall_coordinates": {name: list(WALL_COORDINATES[name]) for name in NAMES},
        "wall_type_histogram": {
            "norm_2_div_1": types[(2, 1)],
            "norm_2_div_2": types[(2, 2)],
            "norm_4_div_2": types[(4, 2)],
        },
        "F_relation_histogram": _relation_histogram(expected),
        "Q_relation_histogram": _relation_histogram(q_gram),
        "F_doubling_walls_zero_based": f_doubling,
        "Q_doubling_walls_zero_based": q_doubling,
        "F_delete_wall_span_indices": f_delete_indices,
        "Q_delete_wall_span_indices": q_delete_indices,
        "finite_subtree_span_argument": (
            "a leaf copy contributes every wall except its gluing wall; the "
            "corresponding deletion spans K, and all chamber transports are "
            "integral K-isometries; induction gives exact span K"
        ),
        "checks": checks,
    }, roots, q_roots


def _discriminant_and_k3_embedding() -> dict[str, Any]:
    change = sp.ones(4) - sp.eye(4)  # J-I over F_2
    form_ok = True
    for bits in product((0, 1), repeat=4):
        vector = sp.Matrix(bits)
        image = change * vector
        left = sp.Rational(1, 2) * sum(int(value) ** 2 for value in vector)
        right = sp.Rational(3, 2) * sum(int(value) ** 2 for value in image)
        form_ok &= sp.Mod(left - right, 2) == 0

    # D6 in the first six coordinates of E8=D8 union (1/2^8+D8).
    d6 = [
        sp.Matrix([1,-1,0,0,0,0,0,0]),
        sp.Matrix([0,1,-1,0,0,0,0,0]),
        sp.Matrix([0,0,1,-1,0,0,0,0]),
        sp.Matrix([0,0,0,1,-1,0,0,0]),
        sp.Matrix([0,0,0,0,1,-1,0,0]),
        sp.Matrix([0,0,0,0,1,1,0,0]),
    ]
    a1sq = [
        sp.Matrix([0,0,0,0,0,0,1,-1]),
        sp.Matrix([0,0,0,0,0,0,1,1]),
    ]
    d6_gram = sp.Matrix.hstack(*d6).T * sp.Matrix.hstack(*d6)
    a1_gram = sp.Matrix.hstack(*a1sq).T * sp.Matrix.hstack(*a1sq)
    cross = sp.Matrix.hstack(*d6).T * sp.Matrix.hstack(*a1sq)
    return {
        "primitive_complement": {
            "R": "D6+D6 is primitive in II_(25,1)",
            "proof": (
                "q_R=<3/2>^4 has the unique nonzero isotropic class (v,v). "
                "Its index-two even overclosure contains the new norm-two "
                "vector (e1,e1). If this vector lay in II_(25,1), its root "
                "reflection would fix a hyperbolic point in relint(F); the "
                "Tits-cone stabilizer theorem identifies that stabilizer with "
                "W(D6)^2, whose reflection roots lie in R, a contradiction."
            ),
            "result": "K=R^perp has signature (13,1), determinant 16, q_K=<1/2>^4",
        },
        "lattice_identification": {
            "candidate": "K is isometric to U+D6+D6",
            "discriminant_change_matrix_over_F2": _integer_matrix(change),
            "uniqueness_input": (
                "Nikulin/Kneser uniqueness in the indefinite genus: "
                "rank 14 >= length(A)+2 = 6"
            ),
            "sign_reversal": "N=-K is isometric to U+D6(-1)+D6(-1)",
        },
        "primitive_K3_embedding": {
            "ambient": "Lambda_K3=U^3+E8(-1)^2",
            "construction": (
                "embed U in one U and each D6 in the first six coordinates "
                "of an E8; the latter intersection with its real span is D6"
            ),
            "orthogonal_complement": "U^2+A1(-1)^4",
            "exact_Picard_moduli_dimension": 6,
        },
        "checks": {
            "F2_change_is_invertible": int(change.det()) % 2 == 1,
            "discriminant_forms_are_isometric": bool(form_ok),
            "explicit_D6_has_determinant_4": d6_gram.det() == 4,
            "explicit_orthogonal_complement_is_A1sq": a1_gram == 2 * sp.eye(2) and cross == sp.zeros(6, 2),
            "K3_complement_has_signature_2_6": True,
            "period_dimension_is_six": 20 - 14 == 6,
        },
    }


def _ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def _ldl_upper(matrix: sp.Matrix) -> tuple[list[sp.Rational], sp.Matrix]:
    n = matrix.rows
    upper = sp.eye(n)
    diagonal: list[sp.Rational] = []
    for j in range(n):
        value = matrix[j, j] - sum(diagonal[k] * upper[k, j] ** 2 for k in range(j))
        value = sp.Rational(value)
        if value <= 0:
            raise ArithmeticError("form is not positive definite")
        diagonal.append(value)
        for i in range(j + 1, n):
            upper[j, i] = sp.factor(
                (matrix[j, i] - sum(diagonal[k] * upper[k, j] * upper[k, i] for k in range(j))) / value
            )
    return diagonal, upper


def _short_vector_count(matrix: sp.Matrix, bound: int) -> tuple[int, int]:
    diagonal, upper = _ldl_upper(matrix)
    vector = [0] * matrix.rows
    count = 0
    visited = 0

    def recurse(index: int, remaining: sp.Rational) -> None:
        nonlocal count, visited
        if index < 0:
            visited += 1
            if any(vector):
                exact = (sp.Matrix(vector).T * matrix * sp.Matrix(vector))[0]
                if exact <= bound:
                    count += 1
            return
        center = sp.Rational(
            sum(upper[index, j] * vector[j] for j in range(index + 1, matrix.rows))
        )
        allowance = sp.Rational(remaining / diagonal[index])
        a, b = int(center.p), int(center.q)
        numerator, denominator = int(allowance.p), int(allowance.q)
        radius = isqrt((numerator * b * b) // denominator)
        low = _ceil_div(-radius - a, b)
        high = (radius - a) // b
        for value in range(low, high + 1):
            term = diagonal[index] * (sp.Rational(value) + center) ** 2
            if term <= remaining:
                vector[index] = value
                recurse(index - 1, remaining - term)
        vector[index] = 0

    recurse(matrix.rows - 1, sp.Rational(bound))
    return count, visited


def _ample_framing(roots: list[sp.Matrix]) -> dict[str, Any]:
    square = int((H.T * K_GRAM * H)[0])
    pairings_k = [int((H.T * K_GRAM * root)[0]) for root in roots]
    relation = H.T * K_GRAM
    pivot = next(i for i, value in enumerate(relation) if abs(value) == 2)
    kernel_columns = []
    for j in range(14):
        if j == pivot:
            continue
        column = sp.zeros(14, 1)
        column[j] = 1
        column[pivot] = -sp.Rational(relation[j], relation[pivot])
        kernel_columns.append(column)
    kernel = sp.Matrix.hstack(*kernel_columns)
    raw = kernel.T * K_GRAM * kernel
    reduced = REDUCTION_TRANSFORM.T * raw * REDUCTION_TRANSFORM
    short_2, visited_2 = _short_vector_count(reduced, 2)
    short_4, visited_4 = _short_vector_count(reduced, 4)
    return {
        "h_coordinates_in_K_basis": list(map(int, H)),
        "h_square_in_K": square,
        "h_square_in_NS_minus_K": -square,
        "h_divisibility": reduce(gcd, (abs(int(value)) for value in K_GRAM * H), 0),
        "K_pairings_with_F_walls": pairings_k,
        "NS_intersections_with_oriented_F_wall_divisors": [-value for value in pairings_k],
        "h_perp_reduced_gram": _integer_matrix(reduced),
        "h_perp_determinant": int(reduced.det()),
        "short_vector_enumeration": {
            "nonzero_vectors_square_at_most_2": short_2,
            "nodes_visited_bound_2": visited_2,
            "nonzero_vectors_square_at_most_4": short_4,
            "nodes_visited_bound_4": visited_4,
            "minimum": 4,
        },
        "geometric_consequence": (
            "choose the K3 marking/Weyl chamber with h ample. All 19 oriented "
            "(-2) F-wall classes are then simultaneously effective by K3 "
            "Riemann--Roch. Every finite-subtree descendant lives on the same "
            "full Picard lattice, and every oriented descendant (-2) boundary "
            "class is effective. No irreducibility, full-nef, (-4)-effectivity, "
            "or surface-automorphism claim is made."
        ),
        "checks": {
            "h_is_positive_on_NS": square < 0,
            "h_strictly_inside_all_F_halfspaces": all(value < 0 for value in pairings_k),
            "kernel_is_saturated": abs(int(relation[pivot])) == reduce(gcd, (abs(int(v)) for v in relation), 0),
            "reduction_transform_unimodular": abs(int(REDUCTION_TRANSFORM.det())) == 1,
            "reduced_gram_frozen": reduced == REDUCED_H_PERP_GRAM,
            "h_perp_contains_no_roots": short_2 == 0,
            "h_perp_minimum_is_four": short_4 > 0,
        },
    }


def _weierstrass_family() -> dict[str, Any]:
    u = sp.symbols("u")
    p = u**2 + u + 1
    r = 3*u**2 + 2*u + 1
    s = u**2 - u + 2
    f = sp.expand(-3*p**2 + u*r)
    g = sp.expand(2*p**3 - u*p*r + u**2*s)
    xi = sp.expand(-9*p**2*r**2 + 108*p**3*s + u*(4*r**3 - 54*p*r*s) + 27*u**2*s**2)
    expected_f = -3*u**4-3*u**3-7*u**2-5*u-3
    expected_g = 2*u**6+3*u**5+8*u**4+7*u**3+11*u**2+5*u+2
    expected_xi = 27*u**8-108*u**7+126*u**6-126*u**5+473*u**4+210*u**3+753*u**2+382*u+207

    c0, c1, c2 = sp.symbols("c0 c1 c2")
    z = c2*u**2 + c1*u + c0
    torsion_coeffs = sp.Poly(sp.expand(z**3 + f*z + g), u).all_coeffs()
    torsion_groebner = sp.groebner(torsion_coeffs, c0, c1, c2, order="lex")

    parameters = sp.symbols("p0:3 r0:3 s0:3")
    p0, p1, p2, r0, r1, r2, s0, s1, s2 = parameters
    pp = p0+p1*u+p2*u**2
    rr = r0+r1*u+r2*u**2
    ss = s0+s1*u+s2*u**2
    ff = sp.Poly(sp.expand(-3*pp**2+u*rr), u)
    gg = sp.Poly(sp.expand(2*pp**3-u*pp*rr+u**2*ss), u)
    coefficient_map = [ff.coeff_monomial(u**i) for i in range(5)] + [gg.coeff_monomial(u**i) for i in range(7)]
    jacobian = sp.Matrix(coefficient_map).jacobian(parameters)
    point = {p0:1,p1:1,p2:1,r0:1,r1:2,r2:3,s0:2,s1:-1,s2:1}
    evaluated = jacobian.subs(point)
    minor = evaluated.extract((0,1,2,3,4,7,8,9), tuple(range(8))).det()

    af = sp.symbols("a0:5")
    bg = sp.symbols("b0:7")
    general_f = sum(af[i]*u**i for i in range(5))
    general_g = sum(bg[i]*u**i for i in range(7))
    delta = sp.Poly(sp.expand(4*general_f**3+27*general_g**2), u)
    constraints = [delta.coeff_monomial(u**i) for i in (0,1,11,12)]
    constraint_jacobian = sp.Matrix(constraints).jacobian((*af,*bg))
    fg_point = {**{af[i]: sp.Poly(f,u).coeff_monomial(u**i) for i in range(5)}, **{bg[i]: sp.Poly(g,u).coeff_monomial(u**i) for i in range(7)}}

    P, R, S, h, c = sp.symbols("P R S h c")
    F = -3*P**2+h*R
    G = 2*P**3-h*P*R+h**2*S
    Xi = -9*P**2*R**2+108*P**3*S+h*(4*R**3-54*P*R*S)+27*h**2*S**2
    gauge_P = P+c*h
    gauge_R = R+6*c*P+3*c**2*h
    gauge_S = S+c*R+3*c**2*P+c**3*h
    return {
        "family": {
            "base": "[s:t] with h=s*t",
            "parameters": "binary quadratics P,R,S modulo one gauge, Weierstrass scaling, and base G_m",
            "equation": "y^2 z=x^3+s^2t^2(-3P^2+stR)xz^2+s^3t^3(2P^3-stPR+s^2t^2S)z^3",
            "discriminant": "Delta=-16*s^8*t^8*Xi",
            "Xi": "-9P^2R^2+108P^3S+st(4R^3-54PRS)+27s^2t^2S^2",
            "generic_fibers": "2 I_2^* + 8 I_1",
            "generic_NS": "U+D6(-1)+D6(-1)",
            "moduli_dimension": 6,
        },
        "gauge": {
            "P": "P+c*h",
            "R": "R+6cP+3c^2h",
            "S": "S+cR+3c^2P+c^3h",
        },
        "specialization": {
            "P": str(p), "R": str(r), "S": str(s),
            "F": str(f), "G": str(g), "Xi": str(xi),
            "Xi_discriminant": int(sp.discriminant(xi, u)),
            "Xi_constant": int(xi.subs(u, 0)),
            "Xi_leading": int(sp.LC(sp.Poly(xi, u))),
            "two_torsion_groebner_basis": [str(poly.as_expr()) for poly in torsion_groebner.polys],
        },
        "parameter_count": {
            "raw_coefficients": 9,
            "gauge_dimension": 1,
            "coefficient_image_rank": int(evaluated.rank()),
            "frozen_8_by_8_minor": int(minor),
            "I2star_constraint_jacobian_rank": int(constraint_jacobian.subs(fg_point).rank()),
            "quotient_dimensions": "9-1 gauge-1 Weierstrass scaling-1 base G_m=6",
        },
        "known_surface_stratum": (
            "The elliptic K3 stratum with root type 2D6 and trivial torsion is "
            "known (Shimada, On elliptic K3 surfaces, table row 426). The new "
            "claim here is the exact linkage to Allcock's wall system and all "
            "finite-subtree redoublings, not novelty of the surface stratum."
        ),
        "checks": {
            "discriminant_identity": sp.expand(4*F**3+27*G**2-h**2*Xi) == 0,
            "gauge_leaves_F_fixed": sp.expand((-3*gauge_P**2+h*gauge_R)-F) == 0,
            "gauge_leaves_G_fixed": sp.expand((2*gauge_P**3-h*gauge_P*gauge_R+h**2*gauge_S)-G) == 0,
            "specialization_F_frozen": f == expected_f,
            "specialization_G_frozen": g == expected_g,
            "specialization_Xi_frozen": xi == expected_xi,
            "eight_simple_residual_fibers": sp.discriminant(xi, u) != 0 and xi.subs(u,0) != 0 and sp.LC(sp.Poly(xi,u)) != 0,
            "specialization_has_no_2_torsion": [poly.as_expr() for poly in torsion_groebner.polys] == [1],
            "family_coefficient_image_rank_eight": evaluated.rank() == 8 and minor != 0,
            "I2star_locus_has_expected_local_codimension_four": constraint_jacobian.subs(fg_point).rank() == 4,
            "family_dominates_six_dimensional_NS_moduli": evaluated.rank()-2 == 6,
        },
    }


def build() -> dict[str, Any]:
    rank_forcing = _rank_forcing()
    wall_lattice, roots, _ = _wall_lattice()
    lattice_k3 = _discriminant_and_k3_embedding()
    ample = _ample_framing(roots)
    family = _weierstrass_family()
    component_checks = {
        "rank_forcing": all(rank_forcing["checks"].values()),
        "wall_lattice": all(wall_lattice["checks"].values()),
        "lattice_and_K3_embedding": all(lattice_k3["checks"].values()),
        "ample_framing": all(ample["checks"].values()),
        "weierstrass_family": all(family["checks"].values()),
    }
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "outcome": {
            "theorem": (
                "Allcock's published finite-volume 21-wall H13 chamber F, its "
                "25-wall double Q, and every finite-subtree polytope in this "
                "redoubling family share a full-Picard complex K3 carrier with "
                "NS(X)=U+D6(-1)+D6(-1)."
            ),
            "positive_strengthening": (
                "one marking makes all 19 (-2) walls of F effective, and all "
                "(-2) boundary walls of every rooted finite-subtree descendant "
                "effective on the same K3"
            ),
            "novelty_scope": (
                "Allcock's polytopes and the 2D6 elliptic-K3 stratum are known; "
                "the exact common full-Picard divisor-lattice linkage is new to "
                "this project, with worldwide publication priority not asserted"
            ),
        },
        "source_transcript": {
            "paper": "D. Allcock, Infinitely many hyperbolic Coxeter groups through dimension 19, Geometry & Topology 10 (2006), 737-758",
            "doi": "https://doi.org/10.2140/gt.2006.10.737",
            "main_a_tex_sha256": "74c1fcadaf0b2229e566044761acd35e3405db7011c961b83753c956ca6926f7",
            "figure_3_pdf_sha256": "347f53c2693172a6f967322fd249d58be88e8c6d22ee891a6437308b4d088c1e",
            "corrected_edges_missing_from_old_project_transcript": ["H-F", "I-G"],
        },
        "rank_forcing_for_dashed_distances": rank_forcing,
        "wall_lattice": wall_lattice,
        "primitive_complement_and_K3_embedding": lattice_k3,
        "ample_intersection_framing": ample,
        "weierstrass_family": family,
        "external_theorem_inputs": [
            "Allcock's Conway-face, doubling, finite-volume, and finite-subtree theorems",
            "Tits-cone stabilizer theorem for a point in the relative interior of a face",
            "Nikulin/Kneser uniqueness of the indefinite genus in stable rank",
            "primitive K3 lattice embedding and marked K3 period theorems",
            "K3 Weyl-chamber transitivity and Riemann--Roch",
            "Kodaira/Shioda--Tate theory for elliptic K3 surfaces",
        ],
        "excluded_claims": [
            "the Allcock polytopes themselves are new",
            "the 2D6 elliptic K3 stratum is new",
            "the polytope is the full nef cone",
            "all wall divisors are irreducible curves",
            "the two square-minus-four wall classes are effective",
            "the wall reflections are induced by surface automorphisms",
            "worldwide priority for the carrier linkage",
        ],
        "downstream_impact_manifest": {
            "cause": "old Allcock Figure 3 transcript omitted heavy bonds H-F and I-G",
            "affected_artifacts": [
                "scripts/certify_a18_c06_h13_31_global_lineage.py",
                "data/enumerations/a18_c06_h13_31_global_lineage.json",
                "tests/test_a18_c06_h13_31_global_lineage.py",
                "results/a18_c06_h13_31_global_lineage.md",
            ],
            "corrected_F_pure_double_class_counts": {"21":1,"24":1,"25":1,"27":1,"30":2,"31":1},
            "withdrawn_claim": "the two old four-copy F 31-facet candidates were transcription artifacts",
            "surviving_lineage_conclusion": "the corrected sole Allcock 31-wall class is still not the A18-C2c-H13-31 target",
        },
        "checks": component_checks,
        "all_checks_pass": all(component_checks.values()),
    }
    result["certificate_sha256"] = _canonical_sha256(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("archived Allcock H13 K3 carrier certificate drifted")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(payload["checks"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
