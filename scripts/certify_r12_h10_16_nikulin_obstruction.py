#!/usr/bin/env python3
"""Exact obstruction to the proposed Nikulin lift of R12 H10/16.

The archived H10/16 action on

    M = U + D8(-1) + A1(-1)^2

has coinvariant lattice <-4>.  A symplectic involution extending it would
have full anti-invariant lattice E8(-2), so the complement

    T = U + U(2) + D6(-1)

would have to contain a *primitive* E7(-2).  This script proves that no such
embedding exists.  The obstruction is the quadratic space

    qbar_T(x mod 2T) = x^2/2 mod 2:

its maximal totally singular dimension is six, whereas a primitive
E7(-2) would inject a seven-dimensional totally singular subspace.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "enumerations"
    / "roulleau_rank12_graph_all_subgroup_folds.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "r12_h10_16_nikulin_obstruction.json"
)
EXPECTED_SOURCE_SHA256 = (
    "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def matrix_digest(matrix: sp.Matrix) -> str:
    payload = json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _rank_bits(rows: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = int(row)
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _vector_bits(vector: Sequence[int]) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(vector))


def _bits_vector(bits: int, dimension: int) -> sp.Matrix:
    return sp.Matrix([(bits >> index) & 1 for index in range(dimension)])


def _combine_bits(basis: Sequence[int], coefficients: int) -> int:
    answer = 0
    for index, vector in enumerate(basis):
        if (coefficients >> index) & 1:
            answer ^= vector
    return answer


def _mod_two_kernel_basis(gram: sp.Matrix) -> list[int]:
    row_bits = [
        sum((int(gram[row, column]) & 1) << column for column in range(gram.cols))
        for row in range(gram.rows)
    ]
    target_dimension = gram.cols - _rank_bits(row_bits)
    basis: list[int] = []
    for vector in range(1, 1 << gram.cols):
        if any((row & vector).bit_count() & 1 for row in row_bits):
            continue
        if _rank_bits([*basis, vector]) == len(basis) + 1:
            basis.append(vector)
        if len(basis) == target_dimension:
            return basis
    raise AssertionError("failed to construct the mod-two kernel")


def _qbar(gram: sp.Matrix, vector_bits: int) -> int:
    vector = _bits_vector(vector_bits, gram.rows)
    square = int((vector.T * gram * vector)[0])
    if square % 2:
        raise AssertionError("qbar was evaluated on a non-even lattice")
    return (square // 2) & 1


def _bilinear_mod_two(gram: sp.Matrix, left_bits: int, right_bits: int) -> int:
    left = _bits_vector(left_bits, gram.rows)
    right = _bits_vector(right_bits, gram.rows)
    return int((left.T * gram * right)[0]) & 1


def _discriminant_numerator_mod_eight(gram: sp.Matrix, vector_bits: int) -> int:
    """Return b^T G b mod 8 for the dual half-vector b/2."""

    vector = _bits_vector(vector_bits, gram.rows)
    return int((vector.T * gram * vector)[0]) % 8


def _coefficient_of(basis: Sequence[int], vector: int) -> int:
    for coefficients in range(1 << len(basis)):
        if _combine_bits(basis, coefficients) == vector:
            return coefficients
    raise AssertionError("vector is not in the displayed binary basis")


def _discriminant_action(basis: Sequence[int], action: sp.Matrix) -> sp.Matrix:
    result = sp.zeros(len(basis))
    for source, vector_bits in enumerate(basis):
        image = action * _bits_vector(vector_bits, action.cols)
        image_bits = _vector_bits([int(value) for value in image])
        coefficients = _coefficient_of(basis, image_bits)
        for target in range(len(basis)):
            result[target, source] = (coefficients >> target) & 1
    return result


def _inverse_mod_two(matrix: sp.Matrix) -> sp.Matrix:
    size = matrix.rows
    rows = [
        [int(matrix[row, column]) & 1 for column in range(size)]
        + [int(row == column) for column in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(size):
            if row != column and rows[row][column]:
                rows[row] = [left ^ right for left, right in zip(rows[row], rows[column])]
    return sp.Matrix([row[size:] for row in rows])


def _find_discriminant_anti_isometry(
    source_gram: sp.Matrix,
    source_basis: Sequence[int],
    target_gram: sp.Matrix,
    target_basis: Sequence[int],
) -> sp.Matrix:
    """Find a deterministic four-bit anti-isometry by exhaustive GL(4,2)."""

    dimension = len(source_basis)
    for images in itertools.permutations(range(1, 1 << dimension), dimension):
        if _rank_bits(images) != dimension:
            continue
        valid = True
        for coefficients in range(1 << dimension):
            target_coefficients = _combine_bits(images, coefficients)
            source_vector = _combine_bits(source_basis, coefficients)
            target_vector = _combine_bits(target_basis, target_coefficients)
            if (
                _discriminant_numerator_mod_eight(source_gram, source_vector)
                + _discriminant_numerator_mod_eight(target_gram, target_vector)
            ) % 8:
                valid = False
                break
        if not valid:
            continue
        answer = sp.zeros(dimension)
        for source, image in enumerate(images):
            for target in range(dimension):
                answer[target, source] = (image >> target) & 1
        return answer
    raise AssertionError("the two discriminant forms are not anti-isometric")


def _standard_complement_gram() -> sp.Matrix:
    """Gram of U + U(2) + D6(-1)."""

    gram = sp.zeros(10)
    gram[0, 1] = gram[1, 0] = 1
    gram[2, 3] = gram[3, 2] = 2
    for index in range(6):
        gram[4 + index, 4 + index] = -2
    for left, right in [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)]:
        gram[4 + left, 4 + right] = gram[4 + right, 4 + left] = 1
    return gram


def build_certificate() -> dict[str, Any]:
    actual_source_hash = sha256(SOURCE)
    if actual_source_hash != EXPECTED_SOURCE_SHA256:
        raise AssertionError("the pinned rank-12 fold certificate changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if not source["all_checks_pass"]:
        raise AssertionError("the upstream rank-12 fold certificate failed")

    stock = next(
        row for row in source["stocks"] if row["stock_key"] == "R12_UD8A1x2"
    )
    survivor = next(row for row in stock["survivors"] if row["class_index"] == 3)
    class_record = next(row for row in stock["classes"] if row["class_index"] == 3)
    geometry = survivor["geometry"]
    volume = geometry["vinberg_volume_criterion"]
    if not (
        survivor["all_checks_pass"]
        and all(geometry["checks"].values())
        and geometry["facet_count"] == 16
        and volume["decision"] == "finite_volume_noncompact"
        and volume["finite_volume_bad_edges"] == []
        and volume["edge_extension_histogram"] == {"2": 426}
        and volume["parabolic_neighborhood_audit_count"] == 14
    ):
        raise AssertionError("the archived H10/16 Vinberg gates changed")

    m_gram = sp.Matrix(stock["basis_gram"])
    wall_rows = sp.Matrix(stock["wall_coordinate_rows"])
    permutation = class_record["representative_generators_one_based"][0]
    m_action = sp.zeros(12)
    for source_index in range(12):
        m_action[:, source_index] = wall_rows[permutation[source_index] - 1, :].T
    archived_action = sp.Matrix(
        survivor["neron_severi_action"]["generator_action_matrices"][0]
    )
    fixed_basis = sp.Matrix(survivor["fixed_lattice"]["basis_rows_in_source_basis"])
    fixed_gram = fixed_basis * m_gram * fixed_basis.T
    fixed_smith_matrix = smith_normal_form(fixed_gram, domain=ZZ)
    fixed_smith = [
        abs(int(fixed_smith_matrix[index, index])) for index in range(fixed_gram.rows)
    ]
    anti_vectors = (m_action + sp.eye(12)).nullspace()
    if len(anti_vectors) != 1:
        raise AssertionError("the archived action no longer has coinrank one")
    coinvariant = anti_vectors[0]
    denominator_lcm = sp.ilcm(*[value.q for value in coinvariant])
    coinvariant = denominator_lcm * coinvariant
    coordinate_gcd = 0
    for value in coinvariant:
        coordinate_gcd = gcd(coordinate_gcd, abs(int(value)))
    coinvariant /= coordinate_gcd
    coinvariant_square = int((coinvariant.T * m_gram * coinvariant)[0])
    coinvariant_divisibility = 0
    for value in m_gram * coinvariant:
        coinvariant_divisibility = gcd(coinvariant_divisibility, abs(int(value)))
    if not (
        m_action == archived_action
        and m_action**2 == sp.eye(12)
        and m_action.T * m_gram * m_action == m_gram
        and sp.factor(m_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 11 * (sp.Symbol("lambda") + 1)
        and m_gram.det() == -16
        and exact_inertia(m_gram) == (1, 11, 0)
        and fixed_gram == sp.Matrix(survivor["fixed_lattice"]["gram"])
        and fixed_smith == [*[1] * 8, 2, 2, 4]
        and coinvariant_square == -4
        and coinvariant_divisibility == 2
    ):
        raise AssertionError("the H10/16 Picard action reconstruction failed")

    t_gram = _standard_complement_gram()
    t_smith_matrix = smith_normal_form(t_gram, domain=ZZ)
    t_smith = [
        abs(int(t_smith_matrix[index, index])) for index in range(t_gram.rows)
    ]
    if not (
        t_gram.det() == 16
        and exact_inertia(t_gram) == (2, 8, 0)
        and t_smith == [*[1] * 6, 2, 2, 2, 2]
        and all(t_gram[index, index] % 2 == 0 for index in range(10))
    ):
        raise AssertionError("the standard (10,4,1) complement failed")

    # The discriminant anti-isometry exists: the obstruction is later and
    # genuinely stronger than a mismatch of the four-bit finite forms.
    m_discriminant_basis = _mod_two_kernel_basis(m_gram)
    t_discriminant_basis = _mod_two_kernel_basis(t_gram)
    anti_isometry = _find_discriminant_anti_isometry(
        m_gram, m_discriminant_basis, t_gram, t_discriminant_basis
    )
    m_discriminant_action = _discriminant_action(m_discriminant_basis, m_action)
    t_finite_action = (
        anti_isometry
        * m_discriminant_action
        * _inverse_mod_two(anti_isometry)
    ).applyfunc(lambda value: int(value) & 1)
    equivariance = (
        t_finite_action * anti_isometry
        - anti_isometry * m_discriminant_action
    ).applyfunc(lambda value: int(value) & 1)
    q_table = []
    for coefficients in range(16):
        source_vector = _combine_bits(m_discriminant_basis, coefficients)
        target_coefficients_vector = anti_isometry * _bits_vector(coefficients, 4)
        target_coefficients = _vector_bits(
            [int(value) for value in target_coefficients_vector]
        )
        target_vector = _combine_bits(t_discriminant_basis, target_coefficients)
        source_numerator = _discriminant_numerator_mod_eight(
            m_gram, source_vector
        )
        target_numerator = _discriminant_numerator_mod_eight(
            t_gram, target_vector
        )
        if (source_numerator + target_numerator) % 8:
            raise AssertionError("the displayed discriminant graph is not anti-isometric")
        q_table.append(
            {
                "source_coefficients": [
                    (coefficients >> index) & 1 for index in range(4)
                ],
                "target_coefficients": [
                    (target_coefficients >> index) & 1 for index in range(4)
                ],
                "source_bGb_mod_8": source_numerator,
                "target_bGb_mod_8": target_numerator,
            }
        )
    if equivariance != sp.zeros(4):
        raise AssertionError("the transported four-bit actions are not equivariant")

    # Exact mod-two obstruction.  The polar radical has dimension four, but
    # qbar is a nonzero linear form on it, leaving a three-dimensional
    # singular radical.  The quotient polar space has dimension six, hence
    # Witt index at most three.  Therefore every totally singular subspace
    # has dimension at most 3+3=6.
    polar_row_bits = [
        sum((int(t_gram[row, column]) & 1) << column for column in range(10))
        for row in range(10)
    ]
    polar_rank = _rank_bits(polar_row_bits)
    radical_basis = _mod_two_kernel_basis(t_gram)
    radical_vectors = [
        _combine_bits(radical_basis, coefficients) for coefficients in range(16)
    ]
    singular_radical = [vector for vector in radical_vectors if _qbar(t_gram, vector) == 0]
    qbar_one_radical = [vector for vector in radical_vectors if _qbar(t_gram, vector) == 1]
    sharp_singular_basis_rows = [
        [0, 1, 1, 1, 1, 0, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 1, 0, 1, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        [0, 1, 1, 1, 0, 0, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 0, 0, 0, 1, 1],
    ]
    sharp_singular_basis = [_vector_bits(row) for row in sharp_singular_basis_rows]
    if not (
        polar_rank == 6
        and len(radical_basis) == 4
        and len(singular_radical) == 8
        and _rank_bits(singular_radical) == 3
        and len(qbar_one_radical) == 8
        and _rank_bits(sharp_singular_basis) == 6
        and all(_qbar(t_gram, vector) == 0 for vector in sharp_singular_basis)
        and all(
            _bilinear_mod_two(t_gram, left, right) == 0
            for left in sharp_singular_basis
            for right in sharp_singular_basis
        )
    ):
        raise AssertionError("the sharp mod-two Witt obstruction failed")

    certificate: dict[str, Any] = {
        "schema": "r12_h10_16_nikulin_obstruction_v1",
        "all_checks_pass": True,
        "source_file_sha256": {str(SOURCE.relative_to(ROOT)): actual_source_hash},
        "case_id": "R12-D8A1x2-C2b-H10-16",
        "archived_positive_polytope": {
            "dimension": 10,
            "facet_count": 16,
            "finite_volume": True,
            "compact": False,
            "finite_vertex_count": geometry["finite_vertex_count"],
            "ideal_vertex_count": geometry["ideal_vertex_count"],
            "elliptic_edge_count": volume["elliptic_edge_count"],
            "edge_extension_histogram": volume["edge_extension_histogram"],
            "endpoint_type_histogram": volume["endpoint_type_histogram"],
            "bad_edge_count": len(volume["finite_volume_bad_edges"]),
            "parabolic_neighborhood_audit_count": volume[
                "parabolic_neighborhood_audit_count"
            ],
            "canonical_full_audit_sha256": volume[
                "canonical_full_audit_sha256"
            ],
            "simple_roots_complete": volume["simple_roots_complete"],
        },
        "picard_action": {
            "lattice": "U+D8(-1)+A1(-1)^2",
            "rank_signature_determinant": [12, [1, 11], -16],
            "action_order": 2,
            "action_characteristic_polynomial": "(x-1)^11*(x+1)",
            "action_sha256": matrix_digest(m_action),
            "fixed_lattice_rank_determinant": [11, 16],
            "fixed_lattice_smith_diagonal_absolute": fixed_smith,
            "fixed_lattice_is_two_elementary": False,
            "fixed_lattice_Z4_factor_count": fixed_smith.count(4),
            "coinvariant_basis": matrix_payload(coinvariant.T),
            "coinvariant_lattice": "<-4>",
            "coinvariant_divisibility": coinvariant_divisibility,
        },
        "k3_complement": {
            "lattice": "U+U(2)+D6(-1)",
            "rank_signature_determinant": [10, [2, 8], 16],
            "smith_diagonal_absolute": t_smith,
            "two_elementary_invariants": [10, 4, 1],
            "uniqueness_scope": (
                "the indefinite even 2-elementary classification identifies "
                "this as the complement for every primitive K3-lattice "
                "embedding of M"
            ),
            "gram": matrix_payload(t_gram),
        },
        "finite_discriminant_gate": {
            "status": "passes",
            "source_half_vector_basis_bits": m_discriminant_basis,
            "target_half_vector_basis_bits": t_discriminant_basis,
            "anti_isometry_source_to_target_mod_two": matrix_payload(anti_isometry),
            "source_action": matrix_payload(m_discriminant_action),
            "transported_target_action": matrix_payload(t_finite_action),
            "actions_equivariant": True,
            "all_16_classes_checked": True,
            "quadratic_form_table": q_table,
        },
        "primitive_E7_minus_2_obstruction": {
            "quadratic_space": "qbar(x mod 2T)=x^2/2 mod 2",
            "ambient_mod_two_dimension": 10,
            "polar_rank": polar_rank,
            "polar_radical_dimension": len(radical_basis),
            "polar_radical_basis_bits": radical_basis,
            "qbar_is_nonzero_on_polar_radical": True,
            "singular_radical_dimension": _rank_bits(singular_radical),
            "nondegenerate_quotient_dimension": 6,
            "nondegenerate_quotient_maximal_isotropic_dimension": 3,
            "maximal_totally_singular_dimension": 6,
            "bound_is_sharp": True,
            "sharp_six_dimensional_basis": sharp_singular_basis_rows,
            "E7_minus_2_mod_two_dimension": 7,
            "E7_minus_2_qbar_is_identically_zero": True,
            "primitive_embedding_exists": False,
            "proof": (
                "primitivity would inject E7(-2)/2E7(-2) into T/2T; "
                "the image would be a seven-dimensional totally qbar-singular "
                "subspace, contradicting the sharp dimension-six bound"
            ),
        },
        "k3_consequence": {
            "symplectic_Nikulin_lift_exists": False,
            "reason": (
                "a primitive square-minus-four vector in E8(-2) is a scaled "
                "E8 root, and its orthogonal is E7(-2); the complement "
                "obstruction therefore excludes every symplectic involution "
                "extending the archived action"
            ),
            "nonsymplectic_involution_with_exact_fixed_NS_exists": False,
            "nonsymplectic_reason": (
                "the invariant cohomology lattice of a nonsymplectic K3 "
                "involution is 2-elementary, whereas the required fixed "
                "Picard lattice has Smith factors 2,2,4"
            ),
            "minimum_remaining_automorphism_order": 4,
            "remaining_positive_lane": (
                "seek a special-period order-four or higher action with no "
                "new fixed Picard classes, a rootless enlarged Picard "
                "coinvariant, the displayed transported four-bit action, and "
                "an invariant Kahler class in the archived strict section"
            ),
        },
        "theorem": (
            "The finite-volume noncompact R12 H10/16 Coxeter polytope remains "
            "valid, but its proposed special-period symplectic Nikulin K3 "
            "realization is impossible; any exact fixed-rank-11 K3 "
            "automorphism realization must have order at least four"
        ),
        "novelty_scope": (
            "new sharp obstruction/criterion for the algebraic carrier only; "
            "no negative statement about the already certified abstract "
            "H10/16 polytope and no global no-go for higher-order K3 lifts"
        ),
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(args.output.read_text(encoding="utf-8")):
            raise SystemExit("recomputed certificate differs from archive")
    else:
        args.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"output": str(args.output), "all_checks_pass": True}, indent=2))


if __name__ == "__main__":
    main()
