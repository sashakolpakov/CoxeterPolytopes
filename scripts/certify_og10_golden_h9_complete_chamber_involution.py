#!/usr/bin/env python3
"""Exact certificate for the intrinsic Pell-fold symmetry of the OG10 H9 slice.

The complete fixed chamber lives in

    L9 = <-2*phi> + E8 + <2>

with coordinates ``(x,E8^8,u)``.  The Pell translation on the ``(x,u)``
plane folds back across one norm-two wall to an involution of the complete
chamber.  This script also audits whether that involution can preserve the
ambient golden H10 carrier.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Iterable, Sequence

from sympy import Matrix as SympyMatrix, ZZ, eye
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from scripts import certify_og10_golden_h10_carrier as carrier
from scripts import certify_og10_golden_h9_fixed_slice as fixed_slice
from scripts.probe_og10_golden_chamber_folding import _interior_point


ROOT = Path(__file__).resolve().parents[1]
EXTENDED_PREFIX = (
    ROOT / "data/candidates/og10_golden_vinberg_through_4_plus_5phi.json"
)
FIXED_COORDINATES = (0, *range(2, 11))

Vector = tuple[Golden, ...]
RowMatrix = tuple[Vector, ...]


def _sum(values: Iterable[Golden]) -> Golden:
    return sum(values, start=ZERO)


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return _sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def _identity(size: int) -> RowMatrix:
    return tuple(
        tuple(ONE if i == j else ZERO for j in range(size))
        for i in range(size)
    )


def _apply(matrix: RowMatrix, vector: Sequence[Golden]) -> Vector:
    return tuple(
        _sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def _compose(left: RowMatrix, right: RowMatrix) -> RowMatrix:
    return tuple(
        tuple(
            _sum(left[i][k] * right[k][j] for k in range(len(left)))
            for j in range(len(right))
        )
        for i in range(len(left))
    )


def _matrix_isometry(matrix: RowMatrix, gram: RowMatrix) -> bool:
    columns = tuple(
        tuple(matrix[i][j] for i in range(len(matrix)))
        for j in range(len(matrix))
    )
    return all(
        _pair(columns[i], gram, columns[j]) == gram[i][j]
        for i in range(len(matrix))
        for j in range(len(matrix))
    )


def _reflection(normal: Vector, gram: RowMatrix) -> RowMatrix:
    square = _pair(normal, gram, normal)
    covector = tuple(
        _sum(gram[i][j] * normal[j] for j in range(len(normal)))
        for i in range(len(normal))
    )
    return tuple(
        tuple(
            (ONE if i == j else ZERO)
            - Golden(2) * normal[i] * covector[j] / square
            for j in range(len(normal))
        )
        for i in range(len(normal))
    )


def _block_matrix(entries: tuple[Golden, Golden, Golden, Golden]) -> RowMatrix:
    result = [list(row) for row in _identity(10)]
    result[0][0], result[0][9], result[9][0], result[9][9] = entries
    return tuple(tuple(row) for row in result)


def _trace_action(field_matrix: RowMatrix) -> SympyMatrix:
    """Restriction of scalars for an O-linear matrix, in (e,phi*e) order."""

    rank = len(field_matrix)
    result = SympyMatrix.zeros(2 * rank)
    for j in range(rank):
        for i in range(rank):
            entry = field_matrix[i][j]
            if not entry.is_integral:
                raise ArithmeticError("field matrix is not O-integral")
            a, b = int(entry.a), int(entry.b)
            result[2 * i, 2 * j] = a
            result[2 * i + 1, 2 * j] = b
            result[2 * i, 2 * j + 1] = b
            result[2 * i + 1, 2 * j + 1] = a + b
    return result


def _full_h2_fixed_and_anti_bases() -> tuple[SympyMatrix, SympyMatrix, SympyMatrix]:
    """Return ambient Gram and saturated h-fixed/anti-fixed bases."""

    ambient_h_columns: list[Vector] = [
        (ONE, ZERO, *(ZERO for _ in range(9))),
        (ONE - PHI, Golden(-1), *(ZERO for _ in range(9))),
    ]
    ambient_h_columns.extend(
        tuple(ONE if i == j else ZERO for i in range(11))
        for j in range(2, 11)
    )
    ambient_h_rows = tuple(
        tuple(ambient_h_columns[j][i] for j in range(11))
        for i in range(11)
    )
    original_h = _trace_action(ambient_h_rows)
    full_change = (
        SympyMatrix(carrier.trace_basis_change())
        * carrier.t11_reduction()
        * carrier.core_small_reordering()
    )
    trace_h = full_change.inv() * original_h * full_change
    direct_h = SympyMatrix.diag(trace_h, -eye(2))
    overlattice_change = SympyMatrix.diag(eye(18), carrier.overlattice_basis())
    ambient_h = overlattice_change.inv() * direct_h * overlattice_change
    direct_gram = SympyMatrix.diag(
        carrier.core_gram(), carrier.small_n_gram(), carrier.complement_gram()
    )
    ambient_gram = overlattice_change.T * direct_gram * overlattice_change

    # SymPy's rational nullspace basis is integral and saturated in this case.
    fixed_basis = SympyMatrix.hstack(*(ambient_h - eye(24)).nullspace())
    anti_equations = DomainMatrix.from_Matrix(
        ambient_h + eye(24), fmt="dense"
    ).convert_to(ZZ)
    diagonal, _, right_change = smith_normal_decomp(anti_equations)
    diagonal_matrix = diagonal.to_Matrix()
    zero_columns = [
        j
        for j in range(24)
        if all(diagonal_matrix[i, j] == 0 for i in range(24))
    ]
    anti_basis = right_change.to_Matrix()[:, zero_columns]
    if not all(entry.is_Integer is True for entry in fixed_basis):
        raise ArithmeticError("fixed basis is unexpectedly nonintegral")
    if not all(entry.is_Integer is True for entry in anti_basis):
        raise ArithmeticError("anti basis is unexpectedly nonintegral")
    return ambient_gram, fixed_basis, anti_basis


def _fixed_field_embedding_in_full_h2() -> tuple[SympyMatrix, SympyMatrix]:
    """Embed the obvious O-basis of y=0 and return its g action."""

    fixed_indices = (0, *range(2, 11))
    original_inclusion = SympyMatrix.zeros(22, 20)
    for target, source in enumerate(fixed_indices):
        original_inclusion[2 * source, 2 * target] = 1
        original_inclusion[2 * source + 1, 2 * target + 1] = 1
    full_change = (
        SympyMatrix(carrier.trace_basis_change())
        * carrier.t11_reduction()
        * carrier.core_small_reordering()
    )
    n_inclusion = full_change.inv() * original_inclusion
    overlattice_change = SympyMatrix.diag(eye(18), carrier.overlattice_basis())
    ambient_inclusion = overlattice_change.inv() * SympyMatrix.vstack(
        n_inclusion, SympyMatrix.zeros(2, 20)
    )
    folded = _block_matrix((PHI, Golden(-1), PHI, -PHI))
    return ambient_inclusion, _trace_action(folded)


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _encoded_matrix(matrix: RowMatrix) -> list[list[dict[str, str]]]:
    return [[entry.to_json() for entry in row] for row in matrix]


def _sympy_digest(matrix: SympyMatrix) -> str:
    payload = [
        [str(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]
    encoded = json.dumps(payload, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _ambient_field_lifts(folded: RowMatrix) -> tuple[RowMatrix, RowMatrix]:
    """The only K-linear lifts preserving y=0, indexed by nu -> +/-nu."""

    lifts = []
    for epsilon in (1, -1):
        result = [list(row) for row in _identity(11)]
        coefficient = (Golden(epsilon) * (PHI - ONE) - ONE) / 2
        result[0][0] = PHI
        result[0][1] = coefficient
        result[0][10] = Golden(-1)
        result[1][0] = ZERO
        result[1][1] = Golden(epsilon)
        result[1][10] = ZERO
        result[10][0] = PHI
        result[10][1] = Golden(Fraction(-1, 2))
        result[10][10] = -PHI
        lift = tuple(tuple(row) for row in result)
        # Its restriction to y=0 is the fixed matrix, with E8 coordinates
        # shifted by one place in the ambient basis.
        fixed_columns = (0, *range(2, 11))
        restriction = tuple(
            tuple(lift[i][j] for j in fixed_columns)
            for i in fixed_columns
        )
        if restriction != folded:
            raise ArithmeticError("ambient lift has the wrong fixed restriction")
        lifts.append(lift)
    return tuple(lifts)  # type: ignore[return-value]


def _full_og10_extension_audit() -> dict[str, object]:
    """Build one positive infinite-order full-lattice extension of g.

    This is an isometry of the abstract Beauville lattice.  It deliberately
    is not claimed to preserve the original period, Neron--Severi lattice, or
    the golden H10 carrier.
    """

    ambient_gram, fixed_basis, anti_basis = _full_h2_fixed_and_anti_bases()
    fixed_embedding, fixed_action = _fixed_field_embedding_in_full_h2()
    basis_change = fixed_basis.gauss_jordan_solve(fixed_embedding)[0]
    action_on_fixed = basis_change * fixed_action * basis_change.inv()
    fixed_gram = fixed_basis.T * ambient_gram * fixed_basis
    anti_gram = anti_basis.T * ambient_gram * anti_basis

    # An explicit unimodular reduction of the anti-invariant lattice to
    # <2> + <2> + <-2> + <-6>.
    reduction = SympyMatrix.hstack(
        SympyMatrix([-6, 0, -2, 1]),
        SympyMatrix([0, 1, 0, 0]),
        SympyMatrix([4, 0, 1, 1]),
        SympyMatrix([7, 0, 2, 0]),
    )
    diagonal_gram = SympyMatrix.diag(2, 2, -2, -6)

    # C0 has the discriminant-gluing action required by g.  Multiplication
    # by the positive reflection diag(-1,1,1,1), which is trivial modulo two,
    # changes it to a positive-orientation, infinite-order lift.
    complement_involution = SympyMatrix(
        [
            [-1, -1, 1, 0],
            [-1, 2, -1, 3],
            [-1, 1, 0, 3],
            [0, -1, 1, -1],
        ]
    )
    positive_reflection = SympyMatrix.diag(-1, 1, 1, 1)
    positive_complement_action = positive_reflection * complement_involution
    anti_action = (
        reduction * positive_complement_action * reduction.inv()
    )
    split_basis = SympyMatrix.hstack(fixed_basis, anti_basis)
    ambient_action = (
        split_basis
        * SympyMatrix.diag(action_on_fixed, anti_action)
        * split_basis.inv()
    )

    # The original period complement is the final rank-two direct summand
    # before gluing.  Its ambient columns provide an exact non-Hodge check.
    overlattice_change = SympyMatrix.diag(eye(18), carrier.overlattice_basis())
    original_period = overlattice_change.inv() * SympyMatrix.vstack(
        SympyMatrix.zeros(22, 2), eye(2)
    )
    augmented_period_rank = SympyMatrix.hstack(
        original_period, ambient_action * original_period
    ).rank()

    characteristic_coefficients = [
        int(value) for value in positive_complement_action.charpoly().all_coeffs()
    ]
    one_eigenspace_dimension = len(
        (positive_complement_action - eye(4)).nullspace()
    )
    checks = {
        "fixed_field_basis_is_the_saturated_h_fixed_lattice": (
            all(entry.is_Integer is True for entry in fixed_embedding)
            and all(entry.is_Integer is True for entry in basis_change)
            and abs(int(basis_change.det())) == 1
        ),
        "anti_lattice_reduces_to_2_2_minus2_minus6": (
            abs(int(reduction.det())) == 1
            and reduction.T * anti_gram * reduction == diagonal_gram
        ),
        "complement_seed_is_an_involution_isometry": (
            complement_involution * complement_involution == eye(4)
            and complement_involution.T
            * (diagonal_gram / 2)
            * complement_involution
            == diagonal_gram / 2
        ),
        "positive_complement_action_is_isometric": (
            positive_complement_action.T
            * (diagonal_gram / 2)
            * positive_complement_action
            == diagonal_gram / 2
        ),
        "full_extension_is_integral": all(
            entry.is_Integer is True for entry in ambient_action
        ),
        "full_extension_is_an_isometry": (
            ambient_action.T * ambient_gram * ambient_action == ambient_gram
        ),
        "full_extension_restricts_to_g": (
            ambient_action * fixed_embedding
            == fixed_embedding * fixed_action
        ),
        "positive_complement_action_has_infinite_order": (
            characteristic_coefficients == [1, -2, 0, 2, -1]
            and one_eigenspace_dimension == 1
        ),
        "full_extension_preserves_positive_orientation": (
            int(positive_complement_action[:2, :2].det()) == 3
        ),
        "full_extension_does_not_preserve_original_period_plane": (
            augmented_period_rank == 4
        ),
    }
    return {
        "fixed_lattice_rank": 20,
        "anti_invariant_lattice_rank": 4,
        "anti_invariant_diagonal_gram": [2, 2, -2, -6],
        "anti_invariant_reduction": [
            [int(reduction[i, j]) for j in range(4)] for i in range(4)
        ],
        "positive_complement_action_in_diagonal_basis": [
            [int(positive_complement_action[i, j]) for j in range(4)]
            for i in range(4)
        ],
        "positive_complement_characteristic_coefficients": characteristic_coefficients,
        "positive_complement_one_eigenspace_dimension": one_eigenspace_dimension,
        "ambient_action_sha256": _sympy_digest(ambient_action),
        "scope": (
            "This is a positive infinite-order isometry of the abstract OG10 "
            "Beauville lattice extending g on H^2(X,Z)^h.  It does not preserve "
            "the original period plane, Neron--Severi lattice, golden H10 "
            "carrier, or H10 Coxeter chamber; hence it is not an automorphism "
            "or birational map of the archived OG10 period model."
        ),
        "checks": checks,
        "verified": all(checks.values()),
    }


def build() -> dict[str, object]:
    ambient_gram = og10_field_gram()
    fixed_gram = tuple(
        tuple(ambient_gram[i][j] for j in FIXED_COORDINATES)
        for i in FIXED_COORDINATES
    )
    archived = fixed_slice._load_roots(EXTENDED_PREFIX, 1144)
    ambient_fixed_roots = tuple(
        root for root in archived if fixed_slice._apply_h(root) == root
    )
    roots = tuple(
        tuple(root[i] for i in FIXED_COORDINATES)
        for root in ambient_fixed_roots
    )

    _, _, ambient_interior = _interior_point(ambient_gram)
    h_interior = fixed_slice._apply_h(ambient_interior)
    invariant_interior = tuple(
        (ambient_interior[i] + h_interior[i]) / 2 for i in range(11)
    )
    interior = tuple(invariant_interior[i] for i in FIXED_COORDINATES)

    pell = _block_matrix((PHI, ONE, PHI, PHI))
    pell_inverse = _block_matrix((PHI, Golden(-1), -PHI, PHI))
    folded = _block_matrix((PHI, Golden(-1), PHI, -PHI))
    eta = (ONE, *(ZERO for _ in range(8)), PHI * PHI)

    pell_point = _apply(pell, interior)
    inverse_point = _apply(pell_inverse, interior)
    pell_violations = tuple(
        i for i, root in enumerate(roots) if _pair(pell_point, fixed_gram, root).sign() > 0
    )
    inverse_violations = tuple(
        i for i, root in enumerate(roots) if _pair(inverse_point, fixed_gram, root).sign() > 0
    )
    pell_wall = roots[pell_violations[0]]
    inverse_wall = roots[inverse_violations[0]]
    pell_fold = _compose(_reflection(pell_wall, fixed_gram), pell)
    inverse_fold = _compose(_reflection(inverse_wall, fixed_gram), pell_inverse)
    endpoint = _apply(folded, interior)
    endpoint_pairings = tuple(_pair(endpoint, fixed_gram, root) for root in roots)

    orthogonal = (ZERO, *endpoint[1:])
    orthogonal_square = _pair(orthogonal, fixed_gram, orthogonal)
    a = endpoint[0]
    violator_bound = orthogonal_square / (
        PHI * (Golden(2) * PHI * a * a - orthogonal_square)
    )
    next_coordinate = Golden(3, 6)
    next_coordinate_square = next_coordinate * next_coordinate

    eta_reflection = _reflection(eta, fixed_gram)
    ambient_lifts = _ambient_field_lifts(folded)
    nu = (PHI - ONE, Golden(2), *(ZERO for _ in range(9)))
    ambient_lift_integrality = tuple(
        all(entry.is_integral for row in lift for entry in row)
        for lift in ambient_lifts
    )
    ambient_lift_normal_signs = tuple(
        _apply(lift, nu) == tuple(Golden(epsilon) * entry for entry in nu)
        for lift, epsilon in zip(ambient_lifts, (1, -1), strict=True)
    )
    full_og10_extension = _full_og10_extension_audit()
    checks = {
        "archive_has_61_fixed_roots": len(roots) == 61,
        "all_archived_fixed_roots_have_norm_two": all(
            _pair(root, fixed_gram, root) == Golden(2) for root in roots
        ),
        "pell_and_inverse_are_isometries": (
            _matrix_isometry(pell, fixed_gram)
            and _matrix_isometry(pell_inverse, fixed_gram)
            and _compose(pell, pell_inverse) == _identity(10)
        ),
        "each_pell_endpoint_violates_exactly_one_known_wall": (
            pell_violations == (9,) and inverse_violations == (8,)
        ),
        "both_one_wall_folds_equal_g": pell_fold == inverse_fold == folded,
        "g_is_integral_isometric_involution": (
            _matrix_isometry(folded, fixed_gram)
            and _compose(folded, folded) == _identity(10)
        ),
        "g_is_reflection_in_eta": (
            _pair(eta, fixed_gram, eta) == Golden(4, 4)
            and eta_reflection == folded
        ),
        "folded_endpoint_is_strict_for_all_61_walls": all(
            value.sign() < 0 for value in endpoint_pairings
        ),
        "violator_bound_matches_claim": violator_bound
        == Golden(
            Fraction(202727309617, 200773372141),
            Fraction(-1021288932, 200773372141),
        ),
        "violator_bound_below_next_uncomputed_level_square": (
            (next_coordinate_square - violator_bound).sign() > 0
        ),
        "only_two_H10_field_lifts_are_both_nonintegral": (
            ambient_lift_normal_signs == (True, True)
            and ambient_lift_integrality == (False, False)
            and ambient_lifts[0][10][1] == Golden(Fraction(-1, 2))
            and ambient_lifts[1][10][1] == Golden(Fraction(-1, 2))
        ),
        "positive_infinite_order_full_OG10_lattice_extension_exists": (
            full_og10_extension["verified"] is True
        ),
    }

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_OG10_golden_H9_complete_chamber_involution",
        "candidate_id": "og10_golden_h9_complete_chamber_involution",
        "verified": all(checks.values()),
        "checks": checks,
        "fixed_lattice": {
            "decomposition": "<-2*phi> + E8 + <2>",
            "coordinate_order": ["x", "E8^8", "u"],
            "rank_over_Z_phi": 10,
            "hyperbolic_dimension": 9,
        },
        "pell_fold": {
            "pell_row_major_x_u_block": [
                [PHI.to_json(), ONE.to_json()],
                [PHI.to_json(), PHI.to_json()],
            ],
            "pell_inverse_row_major_x_u_block": [
                [PHI.to_json(), Golden(-1).to_json()],
                [(-PHI).to_json(), PHI.to_json()],
            ],
            "pell_violated_wall_index_in_61": pell_violations[0],
            "pell_violated_wall": _encoded_vector(pell_wall),
            "inverse_violated_wall_index_in_61": inverse_violations[0],
            "inverse_violated_wall": _encoded_vector(inverse_wall),
            "folded_g_row_major_x_u_block": [
                [PHI.to_json(), Golden(-1).to_json()],
                [PHI.to_json(), (-PHI).to_json()],
            ],
            "reflection_normal_eta": _encoded_vector(eta),
            "reflection_normal_square": _pair(eta, fixed_gram, eta).to_json(),
        },
        "complete_chamber_membership": {
            "archived_fixed_wall_count": len(roots),
            "folded_endpoint": _encoded_vector(endpoint),
            "largest_archived_wall_pairing": max(
                endpoint_pairings,
                key=lambda value: float(value.a)
                + ((1 + 5**0.5) / 2) * float(value.b),
            ).to_json(),
            "possible_violator_x_square_bound": violator_bound.to_json(),
            "next_uncomputed_control_coordinate": next_coordinate.to_json(),
            "next_uncomputed_control_coordinate_square": (
                next_coordinate_square.to_json()
            ),
            "argument": (
                "Every positive-height simple root violated by the endpoint "
                "has x^2 at most the displayed bound.  The ambient Vinberg "
                "enumeration is exhaustive before 3+6*phi, and every fixed "
                "simple root is an ambient y=0 simple root.  Strictness on the "
                "61 archived fixed roots therefore proves complete membership."
            ),
        },
        "golden_H10_extension_obstruction": {
            "ambient_normal_nu": _encoded_vector(nu),
            "normal_square": _pair(nu, ambient_gram, nu).to_json(),
            "unique_lifts_indexed_by_nu_sign": [1, -1],
            "both_lifts_have_y_to_u_coefficient": "-1/2",
            "both_lifts_are_O_nonintegral": [
                not value for value in ambient_lift_integrality
            ],
            "conclusion": (
                "No K-linear integral isometry of the golden H10 lattice "
                "preserves y=0 and restricts to g.  Hence g is not a symmetry "
                "of the complete golden H10 Coxeter chamber."
            ),
        },
        "full_OG10_lattice_extension": full_og10_extension,
        "conclusion": (
            "The Pell translation and its inverse each fold across one wall to "
            "the same reflection g.  The exact unseen-wall bound lies below the "
            "next uncomputed Vinberg coordinate, so g stabilizes the complete H9 "
            "chamber, not merely its 61-wall prefix.  The symmetry has no integral "
            "golden-H10 lift.  It does have a positive infinite-order extension "
            "to the abstract OG10 Beauville lattice, but that transverse extension "
            "does not preserve the archived period or H10 chamber."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded)


if __name__ == "__main__":
    main()
