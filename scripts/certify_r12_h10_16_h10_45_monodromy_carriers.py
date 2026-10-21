#!/usr/bin/env python3
"""Certify explicit algebraic-monodromy carriers for H10/16 and H10/45.

The two ambient elliptic K3 families were previously constructed with exact
generic Neron--Severi lattices

    U + D8(-1) + A1(-1)^2
    U + D4(-1)^2 + A1(-1)^2.

This certificate uses meridians around elementary fiber-enhancement divisors.
For the first family, the two terminal branches in the resolved I4* fiber
exchange; this is the D8 spinor involution.  For the second family, two roots
of the Tate cubic of one I0* fiber exchange; this is a D4 outer transposition.

Explicit unimodular changes of basis identify the resulting full integral
Neron--Severi actions with the archived H10/16 and H10/45 wall actions.  The
induced maps of discriminant quadratic forms are checked on all 16 and 64
classes, respectively.  These are divisor-local-system/parallel-transport
realizations, not automorphisms of one general K3 fiber.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import itertools
import json
from math import gcd
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_ue8a1x2_order4_weierstrass_family import (
    A1X2,
    U,
    block_diag,
    negative_cartan,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
UD8_FAMILY = ROOT / "data/enumerations/ud8a1x2_order4_weierstrass_family.json"
UD4_FAMILY = ROOT / "data/enumerations/ud4x2a1x2_order4_weierstrass_family.json"
UD8_FOLDS = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
UD4_FOLDS = (
    ROOT
    / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
)
H10_16_BOUNDARY = ROOT / "data/enumerations/r12_h10_16_nikulin_obstruction.json"
H10_45_BOUNDARY = (
    ROOT / "data/enumerations/r12_h10_45_e7_embedding_search_boundary.json"
)
OUTPUT = (
    ROOT / "data/enumerations/r12_h10_16_h10_45_monodromy_carriers.json"
)

EXPECTED_SHA256 = {
    "ud8_family": "8494d89aec70059f090795dc2aa757e64a7ac686e878379cb24a021bc2f9bdc3",
    "ud4_family": "55f04a533492001858ead2528aeb4ac686ab0e07323d2076fda8fd9ba5c6302f",
    "ud8_folds": "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    "ud4_folds": "448bcc47937a2763bbd4ea97466fbc48924dc862e88ba8719e5ab98b0cc4c916",
    "h10_16_boundary": "58b5e005bc7612e060728dcc75f2eed2f52cadf826ea767a345b2b54c0be38f9",
    "h10_45_boundary": "bef0ebd838102f9f7648101b0b9f1d587593d01d4bff4e048c67f3faa71d17c9",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(row)] for row in range(matrix.rows)]


def canonical_digest(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def smith_diagonal(gram: sp.Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return [abs(int(smith[index, index])) for index in range(gram.rows)]


def _mod_one(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - sp.floor(value)


def _mod_two(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - 2 * sp.floor(value / 2)


def discriminant_q_histogram(gram: sp.Matrix) -> dict[str, int]:
    """Enumerate q(x)=x^2 mod 2Z on the full discriminant group."""

    inverse = gram.inv()
    generators = [
        tuple(_mod_one(inverse[row, column]) for row in range(gram.rows))
        for column in range(gram.cols)
    ]
    zero = tuple(sp.Rational(0) for _ in range(gram.rows))
    classes = {zero}
    queue: deque[tuple[sp.Rational, ...]] = deque([zero])
    while queue:
        current = queue.popleft()
        for generator in generators:
            candidate = tuple(
                _mod_one(left + right)
                for left, right in zip(current, generator, strict=True)
            )
            if candidate not in classes:
                classes.add(candidate)
                queue.append(candidate)
    histogram = Counter(
        str(_mod_two((sp.Matrix(vector).T * gram * sp.Matrix(vector))[0]))
        for vector in classes
    )
    return {key: histogram[key] for key in sorted(histogram, key=sp.Rational)}


def gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    work = [[int(value) & 1 for value in row] for row in rows]
    rank = 0
    if not work:
        return 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and work[index][column]:
                work[index] = [
                    left ^ right for left, right in zip(work[index], work[rank])
                ]
        rank += 1
    return rank


def gf2_nullspace(matrix: sp.Matrix) -> sp.Matrix:
    work = [[int(value) & 1 for value in matrix.row(row)] for row in range(matrix.rows)]
    row_count, column_count = matrix.shape
    rank = 0
    pivots: list[int] = []
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(row_count):
            if row != rank and work[row][column]:
                work[row] = [
                    left ^ right for left, right in zip(work[row], work[rank])
                ]
        pivots.append(column)
        rank += 1
    free = [column for column in range(column_count) if column not in pivots]
    rows = []
    for free_column in free:
        vector = [0] * column_count
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = work[row][free_column]
        rows.append(vector)
    return sp.Matrix(rows)


def gf2_coordinates(vector: Sequence[int], basis: sp.Matrix) -> list[int]:
    target = tuple(int(value) & 1 for value in vector)
    for bits in itertools.product((0, 1), repeat=basis.rows):
        image = tuple(
            sum(bits[index] * int(basis[index, column]) for index in range(basis.rows))
            & 1
            for column in range(basis.cols)
        )
        if image == target:
            return list(bits)
    raise ArithmeticError("vector is not in the displayed GF(2) span")


def discriminant_action(action: sp.Matrix, half_basis: sp.Matrix) -> sp.Matrix:
    rows = []
    for row in range(half_basis.rows):
        image = half_basis.row(row) * action
        rows.append(gf2_coordinates([int(value) & 1 for value in image], half_basis))
    return sp.Matrix(rows)


def q2(bits: Sequence[int], half_basis: sp.Matrix, gram: sp.Matrix) -> int:
    representative = sp.Matrix([bits]) * half_basis
    numerator = int((representative * gram * representative.T)[0])
    if numerator % 2:
        raise ArithmeticError("half-vector representative is not dual")
    return (numerator // 2) % 4


def fixed_basis_for_swap(rank: int, left: int, right: int) -> sp.Matrix:
    rows = []
    for index in range(rank):
        if index == right:
            continue
        row = [0] * rank
        row[index] = 1
        if index == left:
            row[right] = 1
        rows.append(row)
    return sp.Matrix(rows)


def divisibility(vector: sp.Matrix, gram: sp.Matrix) -> int:
    value = 0
    for entry in vector.T * gram:
        value = gcd(value, abs(int(entry)))
    return value


def finite_form_comparison(
    standard_gram: sp.Matrix,
    archived_gram: sp.Matrix,
    standard_action: sp.Matrix,
    archived_action: sp.Matrix,
    standard_to_archived_basis: sp.Matrix,
) -> dict[str, Any]:
    """Check the exact induced discriminant action under an integral isometry."""

    standard_half = gf2_nullspace(standard_gram)
    archived_half = gf2_nullspace(archived_gram)
    if standard_half.rows != archived_half.rows:
        raise AssertionError("discriminant lengths differ")
    finite_map_rows = []
    for row in range(standard_half.rows):
        image = standard_half.row(row) * standard_to_archived_basis
        finite_map_rows.append(
            gf2_coordinates([int(value) & 1 for value in image], archived_half)
        )
    finite_map = sp.Matrix(finite_map_rows)
    standard_finite_action = discriminant_action(standard_action, standard_half)
    archived_finite_action = discriminant_action(archived_action, archived_half)
    equivariance = (
        standard_finite_action * finite_map - finite_map * archived_finite_action
    ).applyfunc(lambda value: int(value) & 1)
    if gf2_rank(matrix_payload(finite_map)) != finite_map.rows:
        raise AssertionError("the finite discriminant map is not invertible")
    if equivariance != sp.zeros(finite_map.rows):
        raise AssertionError("the finite discriminant actions are not equivariant")

    table = []
    for bits in itertools.product((0, 1), repeat=standard_half.rows):
        coefficient_row = sp.Matrix([bits])
        target_row = (coefficient_row * finite_map).applyfunc(
            lambda value: int(value) & 1
        )
        target_bits = [int(value) for value in target_row]
        source_q2 = q2(bits, standard_half, standard_gram)
        target_q2 = q2(target_bits, archived_half, archived_gram)
        if source_q2 != target_q2:
            raise AssertionError("the displayed finite map does not preserve q")
        table.append(
            {
                "standard_coefficients": list(bits),
                "archived_coefficients": target_bits,
                "q2": source_q2,
            }
        )
    return {
        "all_checks_pass": True,
        "group": f"(Z/2)^{standard_half.rows}",
        "standard_half_vector_basis_bits": matrix_payload(standard_half),
        "archived_half_vector_basis_bits": matrix_payload(archived_half),
        "standard_action": matrix_payload(standard_finite_action),
        "archived_action": matrix_payload(archived_finite_action),
        "standard_to_archived_map": matrix_payload(finite_map),
        "map_rank": gf2_rank(matrix_payload(finite_map)),
        "actions_equivariant": True,
        "all_q_classes_checked": len(table),
        "q_class_table_sha256": canonical_digest(table),
        "q_class_table": table,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "ud8_family": UD8_FAMILY,
        "ud4_family": UD4_FAMILY,
        "ud8_folds": UD8_FOLDS,
        "ud4_folds": UD4_FOLDS,
        "h10_16_boundary": H10_16_BOUNDARY,
        "h10_45_boundary": H10_45_BOUNDARY,
    }
    actual_hashes = {key: sha256(path) for key, path in paths.items()}
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError(
            f"a pinned carrier input changed: {actual_hashes} != {EXPECTED_SHA256}"
        )

    ud8_family = json.loads(UD8_FAMILY.read_text(encoding="utf-8"))
    ud4_family = json.loads(UD4_FAMILY.read_text(encoding="utf-8"))
    ud8_folds = json.loads(UD8_FOLDS.read_text(encoding="utf-8"))
    ud4_folds = json.loads(UD4_FOLDS.read_text(encoding="utf-8"))
    boundary16 = json.loads(H10_16_BOUNDARY.read_text(encoding="utf-8"))
    boundary45 = json.loads(H10_45_BOUNDARY.read_text(encoding="utf-8"))
    if not all(
        source["all_checks_pass"]
        for source in [
            ud8_family,
            ud4_family,
            ud8_folds,
            ud4_folds,
            boundary16,
            boundary45,
        ]
    ):
        raise AssertionError("an upstream carrier input failed")

    d8 = negative_cartan(
        8,
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)],
    )
    d4 = negative_cartan(4, [(0, 1), (1, 2), (1, 3)])

    # ------------------------------------------------------------------
    # H10/16: a meridian around the terminal-slope discriminant of I4*.
    # ------------------------------------------------------------------
    stock8 = next(
        stock
        for stock in ud8_folds["stocks"]
        if stock["stock_key"] == "R12_UD8A1x2"
    )
    case16 = next(
        case for case in stock8["survivors"] if case["class_index"] == 3
    )
    class16 = next(
        case for case in stock8["classes"] if case["class_index"] == 3
    )
    archived_gram16 = sp.Matrix(stock8["basis_gram"])
    archived_action16 = sp.Matrix(
        case16["neron_severi_action"]["generator_action_matrices"][0]
    )
    standard_gram16 = block_diag(U, d8, A1X2)
    standard_action16 = sp.eye(12)
    standard_action16[8, 8] = standard_action16[9, 9] = 0
    standard_action16[8, 9] = standard_action16[9, 8] = 1

    # Rows are the standard U+D8+A1^2 basis in Roulleau's archived basis.
    # In particular rows 8 and 9 are exactly archived walls A10 and A11.
    basis16 = sp.Matrix(
        [
            [1, 3, 0, 0, -3, -6, -9, -12, -3, -6, -6, -1],
            [0, 1, 0, 0, -1, -2, -3, -4, -1, -2, -2, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 3, 0, 0, -3, -6, -9, -12, -3, -6, -6, -1],
            [2, 3, 2, 1, 0, -2, -4, -6, -2, -3, -3, 1],
        ]
    )
    fixed_basis16 = fixed_basis_for_swap(12, 8, 9)
    fixed_gram16 = fixed_basis16 * standard_gram16 * fixed_basis16.T
    coin16 = sp.zeros(12, 1)
    coin16[8], coin16[9] = 1, -1
    finite16 = finite_form_comparison(
        standard_gram16,
        archived_gram16,
        standard_action16,
        archived_action16,
        basis16,
    )

    s, z = sp.symbols("s z")
    r, alpha, beta, gamma = sp.symbols("r alpha beta gamma")
    p0 = -3 * r**2 + 3 - alpha - beta
    q1 = 9 * r - 2 * alpha * r - beta * r - r**3 - 8 + 3 * alpha - 2 * gamma
    q0 = 2 * r**3 - 2 + alpha - gamma - q1
    p8 = -3 * s**3 + alpha * s**2 + beta * s + p0
    q8 = 2 * s**4 - alpha * s**3 + gamma * s**2 + q1 * s + q0
    e8 = sp.expand(4 * p8**3 + 27 * s * q8**2)
    residual8 = sp.cancel(e8 / (s - 1) ** 2)
    witness8 = {
        r: 1,
        alpha: 2,
        beta: 1,
        gamma: sp.Rational(-2, 3),
    }
    slice8 = dict(witness8)
    slice8[gamma] = sp.Rational(-2, 3) + z
    slope_discriminant8 = sp.expand(alpha**2 - 12 * beta - 12 * gamma)
    slice_slope_discriminant8 = sp.expand(slope_discriminant8.subs(slice8))
    boundary_residual8 = sp.expand(residual8.subs(witness8))
    boundary_p8 = sp.expand(p8.subs(witness8))
    boundary_q8 = sp.expand(q8.subs(witness8))
    p8_infinity = -3 + alpha * z + beta * z**2 + p0 * z**3
    q8_infinity = 2 - alpha * z + gamma * z**2 + q1 * z**3 + q0 * z**4
    bracket8 = sp.expand(4 * p8_infinity**3 + 27 * q8_infinity**2)
    boundary_bracket8 = sp.factor(bracket8.subs(witness8))

    expected_residual8 = -112 * s**4 + 16 * s**3 + 120 * s**2 - 105 * s - 108
    expected_bracket8 = -z**3 * (z - 1) ** 2 * (
        108 * z**4 + 105 * z**3 - 120 * z**2 - 16 * z + 112
    )
    local_slope = sp.Symbol("q")
    local_slope_polynomial8 = sp.expand(
        3 * local_slope**2 + alpha * local_slope + beta + gamma
    )
    local_slice_polynomial8 = sp.expand(local_slope_polynomial8.subs(slice8))

    volume16 = case16["geometry"]["vinberg_volume_criterion"]
    target_fixed16 = case16["fixed_lattice"]
    checks16 = {
        "source_family_exact": ud8_family["neron_severi"]["generic_exactness"]
        ["conclusion"]
        == "very general NS=U+D8(-1)+A1(-1)^2",
        "source_family_four_moduli": ud8_family["neron_severi"]
        ["generic_exactness"]["parameter_moduli_dimension"]
        == 4,
        "source_family_saturated": "trivial lattice is saturated"
        in ud8_family["neron_severi"]["generic_exactness"]
        ["torsion_height_argument"],
        "source_local_resolution_exact": ud8_family["checks"]
        ["actual_d8_blowup_chart_actions"],
        "basis_unimodular": abs(int(basis16.det())) == 1,
        "basis_isometry": basis16 * archived_gram16 * basis16.T
        == standard_gram16,
        "full_action_conjugacy": basis16 * archived_action16 * basis16.inv()
        == standard_action16,
        "standard_action_isometry": standard_action16 * standard_gram16
        * standard_action16.T
        == standard_gram16,
        "standard_action_order_two": standard_action16**2 == sp.eye(12),
        "fixed_rank": fixed_gram16.rows == 11,
        "fixed_signature": exact_inertia(fixed_gram16) == (1, 10, 0),
        "fixed_determinant": fixed_gram16.det() == 16,
        "fixed_smith": smith_diagonal(fixed_gram16) == [*[1] * 8, 2, 2, 4],
        "fixed_q_form": discriminant_q_histogram(fixed_gram16)
        == target_fixed16["discriminant_q_histogram_mod_two"],
        "coinvariant_minus_four": (coin16.T * standard_gram16 * coin16)[0]
        == -4,
        "coinvariant_divisibility_two": divisibility(coin16, standard_gram16)
        == 2,
        "finite_action_equivariant": finite16["actions_equivariant"],
        "all_16_q_classes": finite16["all_q_classes_checked"] == 16,
        "archived_finite_action_exact": finite16["archived_action"]
        == boundary16["finite_discriminant_gate"]["source_action"],
        "slice_transverse": slice_slope_discriminant8 == -12 * z,
        "slice_two_terminal_roots": sp.discriminant(
            local_slice_polynomial8, local_slope
        )
        == -12 * z,
        "boundary_residual": boundary_residual8 == expected_residual8,
        "boundary_residual_squarefree": sp.discriminant(
            boundary_residual8, s
        )
        == 2220142134528,
        "boundary_no_additive_collision": sp.resultant(
            boundary_residual8, boundary_p8, s
        )
        == -665523,
        "boundary_finite_i2_stays_exact": boundary_residual8.subs(s, 1)
        == -189,
        "boundary_no_zero_collision": boundary_residual8.subs(s, 0) == -108,
        "boundary_infinity_enhancement": boundary_bracket8 == expected_bracket8,
        "archived_wall_permutation": class16[
            "representative_generators_one_based"
        ][0]
        == [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10, 12, 13, 16, 15, 14, 19, 18, 17],
        "vinberg_decision": volume16["decision"]
        == "finite_volume_noncompact",
        "vinberg_two_endpoints": volume16["edge_extension_histogram"]
        == {"2": 426},
        "vinberg_no_bad_edges": volume16["finite_volume_bad_edges"] == [],
        "vinberg_simple_roots_complete": volume16["simple_roots_complete"],
        "vinberg_digest": volume16["canonical_full_audit_sha256"]
        == "900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226",
    }
    if not all(checks16.values()):
        raise AssertionError(f"an H10/16 monodromy-carrier check failed: {checks16}")

    # ------------------------------------------------------------------
    # H10/45: a meridian around one I0* Tate-cubic discriminant.
    # ------------------------------------------------------------------
    case45 = next(
        case
        for case in ud4_folds["cases"]
        if case["id"] == "R12-D4x2A1x2-S6transposition-H10-45"
    )
    archived_gram45 = sp.Matrix(ud4_folds["coordinate_reconstruction"]["basis_gram"])
    archived_action45 = sp.eye(12)
    archived_action45[1, 1] = archived_action45[2, 2] = 0
    archived_action45[1, 2] = archived_action45[2, 1] = 1
    standard_gram45 = block_diag(U, d4, d4, A1X2)
    standard_action45 = sp.eye(12)
    standard_action45[4, 4] = standard_action45[5, 5] = 0
    standard_action45[4, 5] = standard_action45[5, 4] = 1

    # Rows are the standard U+D4+D4+A1^2 basis in Roulleau's basis.
    # Rows 4 and 5 are archived walls A2 and A3, the S6 transposition.
    basis45 = sp.Matrix(
        [
            [-1, -1, -1, 0, 0, 0, 0, 0, -1, 0, 0, -2],
            [-1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, -2],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [1, 1, 1, 0, 0, 0, 0, 0, 1, -1, -1, 2],
            [-1, -1, -1, 0, -1, 1, 1, 2, 0, 2, 1, -2],
            [-1, -1, -1, -1, 0, 1, 1, 2, 0, 2, 1, -2],
            [1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 2],
            [1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 2],
        ]
    )
    fixed_basis45 = fixed_basis_for_swap(12, 4, 5)
    fixed_gram45 = fixed_basis45 * standard_gram45 * fixed_basis45.T
    coin45 = sp.zeros(12, 1)
    coin45[4], coin45[5] = 1, -1
    finite45 = finite_form_comparison(
        standard_gram45,
        archived_gram45,
        standard_action45,
        archived_action45,
        basis45,
    )

    a, b, g, d = sp.symbols("alpha beta gamma delta")
    u = s - 1
    p4 = -3 + a * u + b * u**2
    q4 = 2 - a * u + g * u**2 + d * u**3
    e4 = sp.expand(4 * p4**3 + 27 * q4**2)
    residual4 = sp.cancel(e4 / (s - 1) ** 2)
    p4_zero = sp.expand(p4.subs(s, 0))
    q4_zero = sp.expand(q4.subs(s, 0))
    cubic_discriminant4 = sp.expand(-4 * p4_zero**3 - 27 * q4_zero**2)
    witness4 = {a: 1, b: 1, g: 1, d: 2}
    slice4 = dict(witness4)
    slice4[d] = 2 + z
    slice_discriminant4 = sp.factor(cubic_discriminant4.subs(slice4))
    boundary_residual4 = sp.expand(residual4.subs(witness4))
    boundary_other4 = sp.cancel(boundary_residual4 / s)
    boundary_p4 = sp.expand(p4.subs(witness4))
    infinity_cubic_discriminant4 = sp.expand(-4 * b**3 - 27 * d**2)
    local_root = sp.Symbol("rho")
    slice_tate_cubic4 = sp.expand(
        (local_root**3 + p4_zero * local_root + q4_zero).subs(slice4)
    )
    volume45 = case45["vinberg_finite_volume_audit"]
    target_fixed45 = case45["fixed_lattice"]
    inherited_gate45 = {
        "relation_matrix_sha256": case45["relation_matrix_sha256"],
        "vinberg_finite_volume_audit": volume45,
        "archived_action": matrix_payload(archived_action45),
    }
    inherited_gate_digest45 = canonical_digest(inherited_gate45)

    checks45 = {
        "source_family_exact": ud4_family["neron_severi"]["generic_exactness"]
        ["conclusion"]
        == "very general NS=U+D4(-1)^2+A1(-1)^2",
        "source_family_four_moduli": ud4_family["neron_severi"]
        ["generic_exactness"]["parameter_moduli_dimension"]
        == 4,
        "source_family_saturated": "trivial lattice is saturated"
        in ud4_family["neron_severi"]["generic_exactness"]
        ["torsion_height_argument"],
        "source_zero_resolution_exact": ud4_family["checks"]
        ["zero_actual_branch_blowup_identity"],
        "basis_unimodular": abs(int(basis45.det())) == 1,
        "basis_isometry": basis45 * archived_gram45 * basis45.T
        == standard_gram45,
        "full_action_conjugacy": basis45 * archived_action45 * basis45.inv()
        == standard_action45,
        "standard_action_isometry": standard_action45 * standard_gram45
        * standard_action45.T
        == standard_gram45,
        "standard_action_order_two": standard_action45**2 == sp.eye(12),
        "fixed_rank": fixed_gram45.rows == 11,
        "fixed_signature": exact_inertia(fixed_gram45) == (1, 10, 0),
        "fixed_determinant": fixed_gram45.det() == 64,
        "fixed_smith": smith_diagonal(fixed_gram45)
        == [*[1] * 6, 2, 2, 2, 2, 4],
        "fixed_q_form": discriminant_q_histogram(fixed_gram45)
        == target_fixed45["discriminant_q_histogram_mod_two"],
        "coinvariant_minus_four": (coin45.T * standard_gram45 * coin45)[0]
        == -4,
        "coinvariant_divisibility_two": divisibility(coin45, standard_gram45)
        == 2,
        "finite_action_equivariant": finite45["actions_equivariant"],
        "all_64_q_classes": finite45["all_q_classes_checked"] == 64,
        "archived_finite_action_exact": [
            [str(value) for value in row] for row in finite45["archived_action"]
        ]
        == boundary45["archived_picard_action"]["discriminant_action"],
        "slice_transverse": sp.expand(
            slice_discriminant4 - 27 * z * (4 - z)
        )
        == 0,
        "slice_cubic_discriminant": sp.expand(
            sp.discriminant(slice_tate_cubic4, local_root)
            - 27 * z * (4 - z)
        )
        == 0,
        "boundary_double_root": sp.factor(slice_tate_cubic4.subs(z, 0))
        == (local_root - 1) ** 2 * (local_root + 2),
        "boundary_residual": sp.factor(boundary_residual4)
        == s * (112 * s**3 - 328 * s**2 + 207 * s + 216),
        "boundary_other_squarefree": sp.discriminant(boundary_other4, s)
        == -14242881024,
        "boundary_no_additive_collision": sp.resultant(
            boundary_other4, boundary_p4, s
        )
        == -275427,
        "boundary_finite_i2_stays_exact": boundary_other4.subs(s, 1) == 207,
        "boundary_zero_enhancement_simple_in_s": sp.diff(e4, s)
        .subs(witness4)
        .subs(s, 0)
        == 216,
        "other_d4_stays_split": infinity_cubic_discriminant4.subs(witness4)
        == -112,
        "archived_wall_transposition": case45["subgroup"]
        ["representative_generators_zero_based"][0][:3]
        == [0, 2, 1],
        "vinberg_decision": volume45["classification"]
        == "finite_volume_noncompact_hyperbolic_Coxeter_polytope",
        "vinberg_finite_volume": volume45["finite_volume"],
        "vinberg_noncompact": not volume45["compact"],
        "vinberg_no_bad_edges": volume45["bad_edge_count"] == 0,
        "vinberg_no_bad_parabolics": volume45[
            "bad_parabolic_neighborhood_count"
        ]
        == 0,
        "vinberg_all_literal_checks": volume45["all_checks_pass"],
    }
    if not all(checks45.values()):
        raise AssertionError(f"an H10/45 monodromy-carrier check failed: {checks45}")

    certificate: dict[str, Any] = {
        "schema": "r12_h10_16_h10_45_monodromy_carriers_v1",
        "all_checks_pass": True,
        "source_file_sha256": actual_hashes,
        "theorem": (
            "The explicit rank-12 elliptic K3 families with generic NS "
            "U+D8(-1)+A1(-1)^2 and U+D4(-1)^2+A1(-1)^2 carry algebraic "
            "divisor-local-system meridians whose integral Picard actions "
            "are exactly the archived H10/16 D8-spinor and H10/45 D4-outer "
            "involutions. Their fixed nef sections are the certified "
            "finite-volume noncompact 16- and 45-facet Coxeter polytopes."
        ),
        "scope": {
            "realization_type": (
                "algebraic family plus divisor local system and parallel "
                "transport around a fiber-enhancement meridian"
            ),
            "single_fiber_automorphism_claimed": False,
            "generic_fiber_picard_exact": True,
            "new_abstract_polytope_claimed": False,
            "novelty": (
                "new explicit equation-to-monodromy-to-Coxeter carrier in "
                "the project; complementary to finite-index ambient carriers"
            ),
        },
        "H10_16": {
            "case_id": "R12-D8A1x2-C2b-H10-16",
            "family": {
                "equation": ud8_family["family"]["affine_equation"],
                "P": ud8_family["family"]["P"],
                "Q": ud8_family["family"]["Q"],
                "generic_NS": "U+D8(-1)+A1(-1)^2",
                "generic_Picard_rank": 12,
                "generic_moduli_dimension": 4,
                "Mordell_Weil_rank": 0,
                "Mordell_Weil_torsion": 0,
                "fiber_configuration": "I4*+2I2+10I1",
                "primitive_saturation_proof": ud8_family["neron_severi"]
                ["generic_exactness"]["torsion_height_argument"],
            },
            "transverse_meridian": {
                "boundary_point": {
                    "r": "1",
                    "alpha": "2",
                    "beta": "1",
                    "gamma": "-2/3",
                },
                "slice": "r=1, alpha=2, beta=1, gamma=-2/3+z",
                "loop": "z=epsilon*exp(i*theta), 0<=theta<=2*pi",
                "normal_discriminant": "alpha^2-12*beta-12*gamma=-12*z",
                "third_blowup_polynomial": str(local_slice_polynomial8),
                "root_labels": "q_+(z)=(-2+sqrt(-12z))/6, q_-(z)=(-2-sqrt(-12z))/6",
                "root_monodromy": "q_+ <-> q_- after one positive meridian",
                "central_fiber_enhancement": "I4*/D8 -> I6*/D10 at z=0",
                "boundary_P": str(boundary_p8),
                "boundary_Q": str(boundary_q8),
                "boundary_residual": str(boundary_residual8),
                "boundary_residual_discriminant": 2220142134528,
                "boundary_residual_P_resultant": -665523,
                "generic_nearby_fiber": (
                    "choose the base point in an arbitrarily small generic "
                    "deformation of this meridian, outside the countable "
                    "Noether--Lefschetz divisors; the local braid class is "
                    "unchanged and the certified generic NS is exact"
                ),
            },
            "component_transport": {
                "local_resolution": (
                    "the first two blowup centers and the transverse branch "
                    "are single-valued; the two roots of the exact third-"
                    "blowup polynomial label the two terminal/spinor curves"
                ),
                "exchanged": ["D8 spinor component S_+", "D8 spinor component S_-"],
                "fixed": [
                    "the six common D8 chain components",
                    "fiber and zero section U",
                    "the I2 component over t=+1",
                    "the I2 component over t=-1",
                ],
                "standard_NS_action": "id_U + outer(D8) + id_A1^2",
            },
            "integral_action_identification": {
                "standard_gram": matrix_payload(standard_gram16),
                "archived_gram": matrix_payload(archived_gram16),
                "standard_action": matrix_payload(standard_action16),
                "archived_action": matrix_payload(archived_action16),
                "standard_basis_rows_in_archived_basis": matrix_payload(basis16),
                "basis_determinant": int(basis16.det()),
                "basis_carries_spinors_to_archived_walls": [10, 11],
                "full_integral_action_pairs_conjugate": True,
                "fixed_gram": matrix_payload(fixed_gram16),
                "fixed_rank_signature_determinant": [11, [1, 10], 16],
                "fixed_smith": smith_diagonal(fixed_gram16),
                "fixed_discriminant_q_histogram": discriminant_q_histogram(
                    fixed_gram16
                ),
                "coinvariant": "<-4>, divisibility 2",
                "finite_discriminant_form": finite16,
            },
            "vinberg": {
                "classification": "finite_volume_noncompact",
                "dimension": 10,
                "facets": 16,
                "finite_vertices": volume16["finite_vertex_count"],
                "cusps": volume16["ideal_vertex_count"],
                "elliptic_edges": volume16["elliptic_edge_count"],
                "edge_extension_histogram": volume16[
                    "edge_extension_histogram"
                ],
                "parabolic_neighborhood_audits": volume16[
                    "parabolic_neighborhood_audit_count"
                ],
                "bad_edges": len(volume16["finite_volume_bad_edges"]),
                "canonical_full_audit_sha256": volume16[
                    "canonical_full_audit_sha256"
                ],
            },
            "checks": checks16,
            "all_checks_pass": True,
        },
        "H10_45": {
            "case_id": "R12-D4x2A1x2-S6transposition-H10-45",
            "family": {
                "equation": ud4_family["family"]["affine_equation"],
                "P": ud4_family["family"]["P"],
                "Q": ud4_family["family"]["Q"],
                "generic_NS": "U+D4(-1)^2+A1(-1)^2",
                "generic_Picard_rank": 12,
                "generic_moduli_dimension": 4,
                "Mordell_Weil_rank": 0,
                "Mordell_Weil_torsion": 0,
                "fiber_configuration": "2I0*+2I2+8I1",
                "primitive_saturation_proof": ud4_family["neron_severi"]
                ["generic_exactness"]["torsion_height_argument"],
            },
            "transverse_meridian": {
                "boundary_point": {
                    "alpha": "1",
                    "beta": "1",
                    "gamma": "1",
                    "delta": "2",
                },
                "slice": "alpha=beta=gamma=1, delta=2+z",
                "loop": "z=epsilon*exp(i*theta), 0<=theta<=2*pi",
                "zero_Tate_cubic": str(slice_tate_cubic4),
                "normal_discriminant": "27*z*(4-z)",
                "boundary_factorization": "(rho-1)^2*(rho+2)",
                "root_monodromy": (
                    "the two roots specializing to rho=1 exchange after one "
                    "positive meridian; the root specializing to -2 is fixed"
                ),
                "central_fiber_enhancement": "I0*/D4 -> IV*/E6 at z=0",
                "boundary_residual": str(boundary_residual4),
                "other_residual_cubic_discriminant": -14242881024,
                "other_residual_P_resultant": -275427,
                "infinity_Tate_cubic_discriminant": -112,
                "generic_nearby_fiber": (
                    "choose the base point in an arbitrarily small generic "
                    "deformation of this meridian, outside the countable "
                    "Noether--Lefschetz divisors; the local braid class is "
                    "unchanged and the certified generic NS is exact"
                ),
            },
            "component_transport": {
                "local_resolution": (
                    "the three roots of xi^3+P(0)xi+Q(0) label the three "
                    "outer I0*/D4 components in the literal first blowup chart"
                ),
                "exchanged": ["first D4 outer component E_+", "first D4 outer component E_-"],
                "fixed": [
                    "the first D4 central component and third outer component",
                    "all four components of the infinity D4",
                    "fiber and zero section U",
                    "the I2 component over t=+1",
                    "the I2 component over t=-1",
                ],
                "standard_NS_action": "id_U + outer_transposition(D4) + id_D4 + id_A1^2",
            },
            "integral_action_identification": {
                "standard_gram": matrix_payload(standard_gram45),
                "archived_gram": matrix_payload(archived_gram45),
                "standard_action": matrix_payload(standard_action45),
                "archived_action": matrix_payload(archived_action45),
                "standard_basis_rows_in_archived_basis": matrix_payload(basis45),
                "basis_determinant": int(basis45.det()),
                "basis_carries_outer_pair_to_archived_walls": [2, 3],
                "full_integral_action_pairs_conjugate": True,
                "fixed_gram": matrix_payload(fixed_gram45),
                "fixed_rank_signature_determinant": [11, [1, 10], 64],
                "fixed_smith": smith_diagonal(fixed_gram45),
                "fixed_discriminant_q_histogram": discriminant_q_histogram(
                    fixed_gram45
                ),
                "coinvariant": "<-4>, divisibility 2",
                "finite_discriminant_form": finite45,
            },
            "vinberg": {
                "classification": "finite_volume_noncompact",
                "dimension": volume45["dimension"],
                "facets": volume45["facet_count"],
                "finite_vertices": volume45["finite_vertex_count"],
                "cusps": volume45["ideal_vertex_count"],
                "elliptic_edges": volume45["elliptic_edge_count"],
                "endpoint_type_histogram": volume45[
                    "edge_endpoint_type_histogram"
                ],
                "parabolic_neighborhood_audits": volume45[
                    "connected_affine_component_count"
                ],
                "bad_edges": volume45["bad_edge_count"],
                "bad_parabolic_neighborhoods": volume45[
                    "bad_parabolic_neighborhood_count"
                ],
                "relation_matrix_sha256": case45["relation_matrix_sha256"],
                "inherited_gate_digest_schema": (
                    "canonical_json_v1(relation_matrix_sha256, "
                    "vinberg_finite_volume_audit, archived_action)"
                ),
                "inherited_gate_sha256": inherited_gate_digest45,
                "upstream_full_certificate_sha256": actual_hashes["ud4_folds"],
            },
            "checks": checks45,
            "all_checks_pass": True,
        },
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--stdout", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if certificate != archived:
            raise SystemExit("recomputed certificate differs from archive")
    elif not args.stdout:
        args.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.stdout:
        print(json.dumps(certificate, indent=2, sort_keys=True))
    else:
        print(
            json.dumps(
                {
                    "output": str(args.output),
                    "all_checks_pass": certificate["all_checks_pass"],
                    "payload_sha256": certificate["certificate_payload_sha256"],
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
