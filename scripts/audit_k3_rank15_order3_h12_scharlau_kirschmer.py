#!/usr/bin/env python3
"""Exact lineage and full-reflection-group audit for the H12 triality fold.

The audit has three independent parts.

1. It freezes and optionally re-parses Scharlau--Kirschmer's ``Lorentz.zip``
   source, selecting the determinant-six entry of ``List12`` and the matching
   row of ``Gram12``.
2. It gives an explicit integral unimodular isometry from the archived Gram
   matrix to the negative of the K3 triality fixed lattice and compares every
   face number and cusp count.
3. It classifies every primitive reflective root of the lattice and replays
   Vinberg's height algorithm.  PARI/GP exhausts the height-zero stabilizer;
   bounded exact Z3 searches prove that no root is missing at any positive
   height through the last wall.  The already independent finite-volume
   certificate then invokes Vinberg's stopping criterion.

No surface automorphism is used or asserted.
"""

from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Iterable, Sequence
import zipfile

import sympy as sp
import z3

from scripts.independent_certify_k3_rho16_h14_fold import (
    _canonical_sha256,
    _exact_inertia,
    _matrix_payload,
    _matrix_sha256,
    _sha256_file,
    _smith_diagonal,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_PATH = ROOT / "data/candidates/k3_rho15_h12_cyclic_fold_audit.json"
CERTIFICATE_PATH = (
    ROOT
    / "data/enumerations/k3_rank15_order3_h12_scharlau_kirschmer_audit.json"
)
DEFAULT_ARCHIVE_PATH = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")

SOURCE_PAGE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
)
SOURCE_ARCHIVE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/Lorentz.zip"
)
SOURCE_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
LIST12_SHA256 = "5a1a4972fffc83ec072e38573bb3848a84b4b42d1747abbc3c708eaa6260ab88"
GRAM12_SHA256 = "acc8f92bf41c533ea041440eccbdc333dfd17eb68c8448db6d4b694215877c0d"

# The sixth entries in List12 and Gram12.  They are frozen here so the exact
# mathematical audit remains reproducible without an implicit network fetch;
# --verify-source-archive re-parses the original ZIP byte for byte.
FROZEN_LIST_ENTRY = (
    6,
    "II_{12,1}([2^{1}]_{1} 3^{-1})",
    [17, 130, 596, 1835, 4017, 6444, 7680, 6807, 4433, 2054, 623, 93],
    3,
)
FROZEN_GRAM = [
    [4, -2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [-2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2],
]
CONTROL_VECTOR = [0, 2, 4, 2, 4, 4, 4, 4, 4, 4, 2, 2, 1]


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _parse_list12(payload: bytes) -> list[tuple[Any, ...]]:
    text = payload.decode("ascii").replace("<", "(").replace(">", ")")
    parsed = ast.literal_eval(text)
    return [tuple(entry) for entry in parsed]


def _parse_gram12(payload: bytes) -> list[list[int]]:
    parsed = ast.literal_eval(payload.decode("ascii"))
    return [[int(value) for value in entry] for entry in parsed]


def verify_source_archive(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    if _sha256_bytes(raw) != SOURCE_ARCHIVE_SHA256:
        raise AssertionError("Lorentz.zip SHA-256 mismatch")
    with zipfile.ZipFile(path) as archive:
        list_raw = archive.read("List12")
        gram_raw = archive.read("Gram12")
        member_names = sorted(archive.namelist())
    if _sha256_bytes(list_raw) != LIST12_SHA256:
        raise AssertionError("List12 SHA-256 mismatch")
    if _sha256_bytes(gram_raw) != GRAM12_SHA256:
        raise AssertionError("Gram12 SHA-256 mismatch")
    lists = _parse_list12(list_raw)
    grams = _parse_gram12(gram_raw)
    if lists[5] != FROZEN_LIST_ENTRY:
        raise AssertionError("the sixth List12 entry differs from the freeze")
    if grams[5] != [value for row in FROZEN_GRAM for value in row]:
        raise AssertionError("the sixth Gram12 entry differs from the freeze")
    return {
        "archive_path": str(path),
        "archive_sha256": SOURCE_ARCHIVE_SHA256,
        "member_count": len(member_names),
        "list12_sha256": LIST12_SHA256,
        "gram12_sha256": GRAM12_SHA256,
        "selected_entry_one_based": 6,
        "selected_entry_equal_to_frozen_literal": True,
    }


def _standard_e8() -> sp.Matrix:
    matrix = 2 * sp.eye(8)
    for left, right in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (6, 7)):
        matrix[left, right] = matrix[right, left] = -1
    return matrix


def _explicit_lattice_isometry(
    candidate_gram: sp.Matrix, source_gram: sp.Matrix
) -> dict[str, Any]:
    """Return P with P^T(-candidate_gram)P = source_gram."""

    identity = sp.eye(13)
    e = [identity[:, index] for index in range(13)]
    u = sp.Matrix([0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 3, 4, 2])
    v = u + e[2]
    projections = [
        vector
        - int((vector.T * candidate_gram * v)[0]) * u
        - int((vector.T * candidate_gram * u)[0]) * v
        for vector in e
    ]
    complement_indices = (0, 1, 3, 5, 6, 7, 8, 9, 10, 11, 12)
    candidate_decomposition = sp.Matrix.hstack(
        u, -v, *[projections[index] for index in complement_indices]
    )

    # qfisom(E8_standard, E8_source) supplied this unimodular matrix X;
    # the equality X^T E8_source X = E8_standard is checked below.
    e8_change = sp.Matrix(
        [
            [0, 1, 0, 0, 0, 0, -2, 4],
            [0, 1, 0, 0, 0, 0, -3, 7],
            [0, 0, 1, 0, 0, 0, -3, 6],
            [0, 0, 0, 1, 0, 0, -3, 5],
            [0, 0, 0, 0, 1, 0, -3, 4],
            [0, 0, 0, 0, 0, 1, -2, 3],
            [0, 0, 0, 0, 0, 0, -1, 2],
            [1, -1, 0, 0, 0, 0, 1, -2],
        ]
    )
    source_decomposition = sp.zeros(13, 13)
    source_decomposition[8, 0] = 1
    source_decomposition[9, 1] = 1
    source_decomposition[11, 2] = 1
    source_decomposition[12, 3] = 1
    source_decomposition[10, 4] = 1
    source_decomposition[:8, 5:13] = e8_change

    common_form = source_decomposition.T * source_gram * source_decomposition
    isometry = candidate_decomposition * source_decomposition.inv()
    expected_common = sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[2, -1], [-1, 2]]),
        sp.Matrix([[2]]),
        _standard_e8(),
    )
    checks = {
        "candidate_decomposition_unimodular": abs(candidate_decomposition.det()) == 1,
        "source_decomposition_unimodular": abs(source_decomposition.det()) == 1,
        "source_E8_change_unimodular": abs(e8_change.det()) == 1,
        "both_decompose_as_U_A2_A1_E8": (
            common_form == expected_common
            and candidate_decomposition.T
            * (-candidate_gram)
            * candidate_decomposition
            == expected_common
        ),
        "explicit_integral_isometry": (
            abs(isometry.det()) == 1
            and isometry.T * (-candidate_gram) * isometry == source_gram
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"explicit lattice isometry failed: {checks}")
    return {
        "sign_convention": (
            "Scharlau--Kirschmer use signature (12,1); the K3 audit uses "
            "signature (1,12), hence the comparison is to -G_K3"
        ),
        "common_orthogonal_decomposition": "U + A2(+1) + A1(+1) + E8(+1)",
        "candidate_decomposition_matrix": _matrix_payload(candidate_decomposition),
        "source_decomposition_matrix": _matrix_payload(source_decomposition),
        "integral_isometry_P": _matrix_payload(isometry),
        "identity": "P^T * (-G_K3) * P = Gram12[6]",
        "integral_isometry_P_sha256": _matrix_sha256(isometry),
        "checks": checks,
    }


def _mod_two(value: Fraction) -> Fraction:
    return value - 2 * (value // 2)


def _reflective_root_classification() -> dict[str, Any]:
    discriminant_values = {
        1: {Fraction(0)},
        2: {Fraction(3, 2)},
        3: {Fraction(4, 3)},
        6: {Fraction(5, 6)},
    }
    trials: list[dict[str, Any]] = []
    allowed: list[tuple[int, int]] = []
    for divisibility in (1, 2, 3, 6):
        for half_norm in range(1, divisibility + 1):
            if divisibility % half_norm:
                continue
            value = _mod_two(Fraction(-2 * half_norm, divisibility**2))
            admitted = value in discriminant_values[divisibility]
            trials.append(
                {
                    "half_negative_norm": half_norm,
                    "divisibility": divisibility,
                    "class_square_mod_2Z": str(value),
                    "admitted": admitted,
                }
            )
            if admitted:
                allowed.append((half_norm, divisibility))
    if allowed != [(1, 1), (1, 2), (3, 3)]:
        raise AssertionError(f"unexpected reflective root types: {allowed}")
    return {
        "argument": (
            "For a primitive r with r^2=-2k and divisibility d, reflection "
            "integrality gives k|d, while r/d has exact order d in A_L. "
            "Here A_L is the orthogonal sum of A2(-1) and A1(-1), so its "
            "six class squares determine the displayed exhaustive list."
        ),
        "discriminant_group": "Z/3 + Z/2 = Z/6",
        "discriminant_form_class_counts": {
            "order_1_q_0": 1,
            "order_2_q_3/2": 1,
            "order_3_q_4/3": 2,
            "order_6_q_5/6": 2,
        },
        "finite_trial_table": trials,
        "allowed_half_norm_and_divisibility": [list(pair) for pair in allowed],
        "conclusion": (
            "Every primitive reflective root has norm -2 and divisibility "
            "1 or 2, or norm -6 and divisibility 3."
        ),
    }


def _gp_matrix(matrix: sp.Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    ) + "]"


def _canonical_sign(vector: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(value) for value in vector)
    first = next(value for value in result if value)
    return result if first > 0 else tuple(-value for value in result)


def _height_zero_root_audit(
    gram: sp.Matrix, normals: sp.Matrix, control: sp.Matrix
) -> dict[str, Any]:
    pairing = control.T * gram
    if list(pairing) != [2] + [0] * 10 + [1, 0]:
        raise AssertionError("unexpected primitive control pairing")

    kernel_columns: list[sp.Matrix] = []
    for index in range(13):
        if index == 11:
            continue
        vector = sp.zeros(13, 1)
        vector[index] = 1
        if index == 0:
            vector[11] = -2
        kernel_columns.append(vector)
    kernel = sp.Matrix.hstack(*kernel_columns)
    kernel_form = -(kernel.T * gram * kernel)
    active_indices = [
        index
        for index in range(normals.cols)
        if int((control.T * gram * normals[:, index])[0]) == 0
    ]
    active = sp.Matrix.hstack(*[normals[:, index] for index in active_indices])
    active_gram = -(active.T * gram * active)
    kernel_coordinates = (kernel.T * kernel).inv() * kernel.T * active
    if any(value.q != 1 for value in kernel_coordinates):
        raise AssertionError("active roots are not in the saturated kernel")
    active_inverse = kernel_coordinates.inv()

    gp_program = (
        f"Q={_gp_matrix(kernel_form)};"
        f"K={_gp_matrix(kernel)};"
        f"G={_gp_matrix(gram)};"
        "R=qfminim(Q,6,,2);print(\"TOTAL|\",R[1]);V=R[3];"
        "for(j=1,matsize(V)[2],z=V[,j];n=z~*Q*z;"
        "if(n==2,print(n,\"|\",Vec(z)),"
        "if(n==6,x=K*z;y=G*x;ok=1;"
        "for(k=1,13,if(y[k]%3,ok=0));"
        "if(ok,print(n,\"|\",Vec(z))))));quit()"
    )
    completed = subprocess.run(
        ["gp", "-fq", "-s", "500000000"],
        input=gp_program,
        text=True,
        capture_output=True,
        timeout=30,
        check=True,
    )
    if completed.stderr.strip():
        raise AssertionError(f"PARI/GP error: {completed.stderr.strip()}")
    lines = completed.stdout.splitlines()
    total_vectors = int(lines[0].split("|", 1)[1])
    records: list[tuple[int, tuple[int, ...], tuple[str, ...]]] = []
    bad_sign_roots: list[tuple[int, ...]] = []
    counts: dict[int, int] = {}
    for line in lines[1:]:
        norm_text, vector_text = line.split("|", 1)
        norm = int(norm_text)
        coordinates = sp.Matrix(
            [int(value) for value in vector_text.strip()[1:-1].split(",")]
        )
        root = kernel * coordinates
        root_tuple = _canonical_sign([int(value) for value in root])
        root = sp.Matrix(root_tuple)
        kernel_coordinate = (kernel.T * kernel).inv() * kernel.T * root
        coefficients = active_inverse * kernel_coordinate
        same_sign = all(value >= 0 for value in coefficients) or all(
            value <= 0 for value in coefficients
        )
        if not same_sign:
            bad_sign_roots.append(root_tuple)
        counts[norm] = counts.get(norm, 0) + 1
        records.append(
            (
                norm,
                root_tuple,
                tuple(str(value) for value in coefficients),
            )
        )
    records.sort()
    expected_active_gram = sp.diag(
        _standard_d10(), sp.Matrix([[2]]), sp.Matrix([[6]])
    )
    checks = {
        "kernel_is_saturated_rank_12": kernel.rank() == 12,
        "kernel_form_positive_definite_determinant_12": (
            kernel_form.is_positive_definite and kernel_form.det() == 12
        ),
        "active_roots_are_D10_A1_A1": active_gram == expected_active_gram,
        "active_root_span_has_index_2_in_kernel": abs(kernel_coordinates.det()) == 2,
        "all_reflective_kernel_roots_enumerated": counts == {2: 91, 6: 1},
        "every_reflective_root_is_positive_or_negative": not bad_sign_roots,
    }
    if not all(checks.values()):
        raise AssertionError(f"height-zero root audit failed: {checks}")
    version = subprocess.run(
        ["gp", "--version-short"],
        text=True,
        capture_output=True,
        timeout=5,
        check=True,
    ).stdout.strip()
    return {
        "method": (
            "PARI/GP qfminim enumerates all vectors of the positive-definite "
            "saturated lattice p-perp through norm 6. The reflective-root "
            "classification filters this to every possible stabilizer root."
        ),
        "pari_gp_version": version,
        "control_pairing_row": [int(value) for value in pairing],
        "saturated_kernel_basis": _matrix_payload(kernel),
        "positive_kernel_gram": _matrix_payload(kernel_form),
        "qfminim_total_signed_vectors_norm_at_most_6": total_vectors,
        "reflective_root_pairs_by_negative_norm": {
            str(-norm): count for norm, count in sorted(counts.items())
        },
        "canonical_reflective_root_and_simple_coefficient_sha256": (
            _canonical_sha256(records)
        ),
        "active_candidate_facets": [index + 1 for index in active_indices],
        "active_positive_gram": _matrix_payload(active_gram),
        "active_type": "D10 + A1(norm 2) + A1(norm 6)",
        "bad_mixed_sign_root_count": len(bad_sign_roots),
        "checks": checks,
    }


def _standard_d10() -> sp.Matrix:
    matrix = 2 * sp.eye(10)
    for left, right in (
        (0, 1),
        (1, 2),
        (1, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
    ):
        matrix[left, right] = matrix[right, left] = -1
    return matrix


def _z3_quadratic(matrix: sp.Matrix, variables: Sequence[z3.ArithRef]) -> z3.ArithRef:
    diagonal = sum(
        int(matrix[index, index]) * variables[index] * variables[index]
        for index in range(matrix.rows)
    )
    off_diagonal = sum(
        2 * int(matrix[left, right]) * variables[left] * variables[right]
        for left in range(matrix.rows)
        for right in range(left + 1, matrix.cols)
    )
    return diagonal + off_diagonal


def _z3_linear(row: sp.Matrix, variables: Sequence[z3.ArithRef]) -> z3.ArithRef:
    return sum(int(row[index]) * variables[index] for index in range(len(variables)))


def _coordinate_bounds(
    positive_form_inverse: sp.Matrix, negative_norm: int, pairing: int
) -> list[int]:
    # Q=2(Gp)(Gp)^T-2G is positive definite and x^TQx=2a^2+2d
    # on x^2=-d, (p,x)=a.  Cauchy in Q gives the exact box below.
    q_value = 2 * pairing**2 + 2 * negative_norm
    return [
        math.isqrt(int(sp.floor(q_value * positive_form_inverse[index, index])))
        for index in range(positive_form_inverse.rows)
    ]


def _positive_height_root_audit(
    gram: sp.Matrix, normals: sp.Matrix, control: sp.Matrix
) -> dict[str, Any]:
    gp = gram * control
    bounding_form = 2 * gp * gp.T - 2 * gram
    if not bounding_form.is_positive_definite:
        raise AssertionError("root-search bounding form is not positive definite")
    bounding_inverse = bounding_form.inv()
    roots = [normals[:, index] for index in range(normals.cols)]
    active = [root for root in roots if int((control.T * gram * root)[0]) == 0]

    # All possible height levels through the final candidate wall.  A norm-6
    # reflective root has divisibility 3, hence its pairing with p is a
    # multiple of 3.
    levels = ((2, 1), (6, 3), (2, 2), (2, 3), (6, 6), (2, 4))
    accepted = list(active)
    audits: list[dict[str, Any]] = []
    for level_index, (negative_norm, pairing) in enumerate(levels):
        expected = [
            root
            for root in roots
            if -int((root.T * gram * root)[0]) == negative_norm
            and int((control.T * gram * root)[0]) == pairing
        ]
        variables = [z3.Int(f"root_{level_index}_{index}") for index in range(13)]
        solver = z3.Solver()
        solver.set(timeout=120_000)
        bounds = _coordinate_bounds(bounding_inverse, negative_norm, pairing)
        for variable, bound in zip(variables, bounds):
            solver.add(variable >= -bound, variable <= bound)
        solver.add(_z3_quadratic(gram, variables) == -negative_norm)
        solver.add(_z3_linear(control.T * gram, variables) == pairing)
        if negative_norm == 6:
            for row in range(13):
                solver.add(_z3_linear(gram.row(row), variables) % 3 == 0)
        for earlier in accepted:
            solver.add(_z3_linear(earlier.T * gram, variables) >= 0)
        for root in expected:
            solver.add(
                z3.Or(
                    *[
                        variables[index] != int(root[index])
                        for index in range(13)
                    ]
                )
            )
        result = solver.check()
        if result != z3.unsat:
            model = (
                [int(solver.model()[variable].as_long()) for variable in variables]
                if result == z3.sat
                else None
            )
            raise AssertionError(
                f"extra/unknown Vinberg root at {(negative_norm, pairing)}: "
                f"{result}, {model}"
            )
        if any(
            int((left.T * gram * right)[0]) < 0
            for left_index, left in enumerate(expected)
            for right in expected[left_index + 1 :]
        ):
            raise AssertionError("same-height candidate roots are incompatible")
        audits.append(
            {
                "negative_norm": negative_norm,
                "control_pairing": pairing,
                "vinberg_height_pairing_squared_over_negative_norm": str(
                    Fraction(pairing**2, negative_norm)
                ),
                "coordinate_bounds": bounds,
                "expected_candidate_facets": [
                    roots.index(root) + 1 for root in expected
                ],
                "extra_root_query": "unsat",
            }
        )
        accepted.extend(expected)

    if len(accepted) != 17:
        raise AssertionError("Vinberg replay did not recover all 17 roots")
    return {
        "method": (
            "At each exact height, bounded integer equations enumerate every "
            "root satisfying Vinberg's nonpositive/acute-wall inequalities. "
            "The positive-definite Q bound makes each Z3 UNSAT query finite."
        ),
        "z3_version": z3.get_version_string(),
        "control_vector": [int(value) for value in control],
        "control_square": int((control.T * gram * control)[0]),
        "positive_definite_bounding_form": _matrix_payload(bounding_form),
        "bounding_form_determinant": int(bounding_form.det()),
        "level_audits": audits,
        "recovered_simple_root_count": len(accepted),
        "all_positive_height_extra_root_queries_unsat": True,
    }


def _all_reflections_integral(gram: sp.Matrix, normals: sp.Matrix) -> bool:
    for column in range(normals.cols):
        root = normals[:, column]
        square = int((root.T * gram * root)[0])
        numerator = 2 * root.T * gram
        if any(int(value) % (-square) for value in numerator):
            return False
        if math.gcd(*[abs(int(value)) for value in root]) != 1:
            return False
    return True


def build_certificate() -> dict[str, Any]:
    candidate = json.loads(CANDIDATE_PATH.read_text())
    gram = sp.Matrix(candidate["matrices"]["fixed_lattice_gram_13"])
    normals = sp.Matrix(candidate["matrices"]["folded_normal_coordinates_13x17"])
    source_gram = sp.Matrix(FROZEN_GRAM)
    control = sp.Matrix(CONTROL_VECTOR)

    lattice_isometry = _explicit_lattice_isometry(gram, source_gram)
    root_classification = _reflective_root_classification()
    height_zero = _height_zero_root_audit(gram, normals, control)
    positive_heights = _positive_height_root_audit(gram, normals, control)

    determinant, genus_symbol, source_faces, source_cusps = FROZEN_LIST_ENTRY
    candidate_f_vector = candidate["complete_face_lattice"][
        "f_vector_vertices_through_polytope"
    ]
    candidate_finite_vertices = candidate["order_three_H12_fold"][
        "finite_vertex_count"
    ]
    candidate_cusps = candidate["order_three_H12_fold"]["ideal_vertex_count"]
    candidate_faces_in_source_order = (
        list(reversed(candidate_f_vector[1:-1])) + [candidate_finite_vertices]
    )
    source_full_f_vector = (
        [source_faces[-1] + source_cusps]
        + list(reversed(source_faces[:-1]))
        + [1]
    )
    candidate_volume_gate = candidate["literal_vinberg_coxiter_prop_1_10"]

    checks = {
        "source_lattice_signature_determinant_snf": (
            _exact_inertia(source_gram) == (12, 1, 0)
            and int(source_gram.det()) == -6
            and _smith_diagonal(source_gram) == [1] * 12 + [6]
        ),
        "candidate_lattice_signature_determinant_snf": (
            _exact_inertia(gram) == (1, 12, 0)
            and int(gram.det()) == 6
            and _smith_diagonal(gram) == [1] * 12 + [6]
        ),
        "explicit_integral_lattice_identity": all(lattice_isometry["checks"].values()),
        "list12_determinant_and_genus_symbol": (
            determinant == 6
            and genus_symbol == "II_{12,1}([2^{1}]_{1} 3^{-1})"
        ),
        "all_source_face_numbers_match": candidate_faces_in_source_order == source_faces,
        "source_cusp_count_matches": candidate_cusps == source_cusps == 3,
        "recovered_full_f_vector_matches": source_full_f_vector == candidate_f_vector,
        "all_17_candidate_normals_are_primitive_reflective_roots": (
            normals.cols == 17 and _all_reflections_integral(gram, normals)
        ),
        "height_zero_simple_system_complete": all(height_zero["checks"].values()),
        "positive_vinberg_layers_complete": positive_heights[
            "all_positive_height_extra_root_queries_unsat"
        ],
        "candidate_polyhedron_has_finite_volume": (
            candidate_volume_gate["finite_volume_noncompact"] is True
            and candidate["order_three_H12_fold"]["hyperideal_vertex_count"] == 0
            and candidate_volume_gate["bad_edges"] == []
        ),
        "vinberg_finite_volume_stopping_criterion_applies": True,
    }
    if not all(checks.values()):
        raise AssertionError(f"Scharlau--Kirschmer H12 audit failed: {checks}")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "artifact_id": "k3_rank15_order3_h12_scharlau_kirschmer_audit",
        "certificate_kind": (
            "exact_source_lineage_lattice_isometry_face_cusp_and_full_root_audit"
        ),
        "source_freeze": {
            "authors": "R. Scharlau and M. Kirschmer",
            "title": "Reflective Lorentzian lattices",
            "official_index_url": SOURCE_PAGE_URL,
            "official_archive_url": SOURCE_ARCHIVE_URL,
            "retrieved_utc_date": "2026-08-06",
            "archive_sha256": SOURCE_ARCHIVE_SHA256,
            "archive_member_count": 36,
            "list12_sha256": LIST12_SHA256,
            "gram12_sha256": GRAM12_SHA256,
            "selected_entry_one_based": 6,
            "list12_entry": [
                determinant,
                genus_symbol,
                source_faces,
                source_cusps,
            ],
            "list_semantics": (
                "List12 records <-det, genus symbol, #faces, #cusps>; "
                "faces[k] counts (12-k)-faces, so its final 93 counts finite "
                "vertices and the 3 ideal vertices are recorded separately."
            ),
            "gram12_entry_13x13": _matrix_payload(source_gram),
            "source_scope_warning": (
                "The archive supplies lattice Grams and face/cusp counts but "
                "does not distribute chamber roots or a weighted Coxeter Gram."
            ),
        },
        "exact_lattice_identity": lattice_isometry,
        "lattice_invariants": {
            "candidate_G_signature": list(_exact_inertia(gram)),
            "candidate_G_determinant": int(gram.det()),
            "candidate_G_smith_diagonal": _smith_diagonal(gram),
            "source_Gram12_signature": list(_exact_inertia(source_gram)),
            "source_Gram12_determinant": int(source_gram.det()),
            "source_Gram12_smith_diagonal": _smith_diagonal(source_gram),
            "genus_symbol": genus_symbol,
            "class_number_statement_from_source": (
                "The source index states that every displayed genus has class number one."
            ),
        },
        "complete_face_and_cusp_comparison": {
            "source_faces_facets_through_finite_vertices": source_faces,
            "candidate_faces_in_source_order": candidate_faces_in_source_order,
            "source_cusps": source_cusps,
            "candidate_cusps": candidate_cusps,
            "source_recovered_full_f_vector_including_ideal_vertices": (
                source_full_f_vector
            ),
            "candidate_full_f_vector": candidate_f_vector,
            "all_face_dimensions_and_cusps_equal": True,
        },
        "full_reflective_root_audit": {
            "classification": root_classification,
            "height_zero": height_zero,
            "positive_heights": positive_heights,
            "finite_volume_stopping_step": {
                "criterion": (
                    "Vinberg's algorithm terminates once the accumulated "
                    "acute-angled root polyhedron has finite volume."
                ),
                "independent_candidate_finite_volume": True,
                "finite_vertices": candidate["order_three_H12_fold"][
                    "finite_vertex_count"
                ],
                "ideal_vertices": candidate_cusps,
                "hyperideal_extreme_rays": candidate["order_three_H12_fold"][
                    "hyperideal_vertex_count"
                ],
                "conclusion": (
                    "The 17 candidate roots are the complete simple-root "
                    "system of the full reflection group W(S)."
                ),
            },
        },
        "identity_hierarchy": {
            "rung_1_same_genus_and_integral_lattice": {
                "proved": True,
                "strength": (
                    "exact integral lattice, stronger than genus: an explicit "
                    "unimodular P satisfies P^T(-G_K3)P=Gram12[6]"
                ),
            },
            "rung_2_same_complete_face_and_cusp_invariants": {
                "proved": True,
                "strength": "all twelve proper-face counts and all three cusps agree",
            },
            "rung_3_explicit_weighted_diagram_facet_permutation": {
                "proved": False,
                "status": "not testable from the distributed archive",
                "reason": (
                    "Lorentz.zip contains no root list or weighted chamber "
                    "Gram, so there is no source facet ordering to permute."
                ),
            },
            "rung_4_full_reflection_group_identity": {
                "proved": True,
                "route": (
                    "independent exhaustive classification of all reflective "
                    "roots plus Vinberg height replay and finite-volume stopping"
                ),
                "qualification": (
                    "This proves equality with the full W(S) underlying the "
                    "census entry and chamber equality up to W(S)-conjugacy; "
                    "it does not manufacture an absent source facet ordering."
                ),
            },
            "highest_defensible_statement": (
                "The K3 triality chamber is the full reflection chamber of "
                "the exact determinant-six Scharlau--Kirschmer lattice, hence "
                "the same abstract chamber up to reflection-group conjugacy."
            ),
        },
        "novelty_conclusion": {
            "new_abstract_hyperbolic_coxeter_polytope": False,
            "new_lattice_or_reflection_chamber": False,
            "possible_new_algebro_geometric_realization_or_derivation": (
                "not decided by this audit"
            ),
            "surface_automorphism_used": False,
            "summary": (
                "The H12 positive example is real and finite-volume, but its "
                "abstract lattice reflection chamber was already present in "
                "the Scharlau--Kirschmer census."
            ),
        },
        "candidate_source": {
            "path": str(CANDIDATE_PATH.relative_to(ROOT)),
            "sha256": _sha256_file(CANDIDATE_PATH),
            "positive_weighted_root_gram_sha256": candidate["matrix_sha256"][
                "positive_normal_coxeter_gram_17"
            ],
        },
        "source_sha256": {
            "certifier": _sha256_file(Path(__file__).resolve()),
            "candidate_certificate": _sha256_file(CANDIDATE_PATH),
            "official_archive": SOURCE_ARCHIVE_SHA256,
            "official_List12_member": LIST12_SHA256,
            "official_Gram12_member": GRAM12_SHA256,
        },
        "checks": checks,
        "commands": [
            "PYTHONPATH=src:. python3 scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py --verify-source-archive /path/to/Lorentz.zip",
            "PYTHONPATH=src:. python3 scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py --verify-certificate",
            "PYTHONPATH=src:. pytest -q tests/test_k3_rank15_order3_h12_scharlau_kirschmer.py",
        ],
    }
    payload["certificate_payload_sha256"] = _canonical_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, Any]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "payload_digest_matches": digest == _canonical_sha256(body),
        "candidate_source_hash_matches": payload["candidate_source"]["sha256"]
        == _sha256_file(CANDIDATE_PATH),
        "full_exact_recomputation_matches": payload == fresh,
        "all_checks_true": all(payload["checks"].values()),
    }
    if not all(checks.values()):
        raise AssertionError(f"lineage certificate verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    parser.add_argument("--verify-source-archive", type=Path)
    arguments = parser.parse_args()
    if arguments.verify_source_archive is not None:
        print(
            json.dumps(
                verify_source_archive(arguments.verify_source_archive),
                indent=2,
                sort_keys=True,
            )
        )
        return
    if arguments.write_certificate:
        payload = build_certificate()
        CERTIFICATE_PATH.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"wrote {CERTIFICATE_PATH}")
        print(payload["certificate_payload_sha256"])
        return
    payload = json.loads(CERTIFICATE_PATH.read_text())
    print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
