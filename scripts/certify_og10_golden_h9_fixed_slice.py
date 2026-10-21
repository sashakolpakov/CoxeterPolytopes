#!/usr/bin/env python3
"""Certify the h-fixed H9 slice of the golden OG10 norm-two chamber.

The complete-chamber involution is

    h(x,y,z,u) = (x + (1-phi)y, -y, z, u).

Its fixed field lattice is ``<-2 phi> + E8 + <2>``.  This script performs
the finite exact part of the fixed-slice argument:

* reconstruct the 640-root h-closed wall set from the archived exhaustive
  prefix;
* extract its 44 h-fixed roots and classify the 298 nontrivial h-orbits;
* extract all 61 fixed roots in the later 1144-wall exhaustive prefix;
* verify the arithmetic inequalities used by the global ultraparallel-orbit
  lemma; and
* certify one spacelike extreme ray of the current 61-wall outer cone; and
* identify the ``u=0`` facet with Bugaenko's compact H8 chamber.

With ``--full-census`` it also recomputes every extreme ray of the historical
44-wall outer cone, freezing the checkpoint from which the 61-wall witness
was obtained.  Wall completeness is deliberately not asserted: later fixed
roots may truncate the current hyperideal escape direction.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
from typing import Iterable, Sequence

from sympy import Matrix as SympyMatrix, ZZ, eye
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import Inertia, inertia
from coxeter_polytopes.known_examples.bugaenko import (
    bugaenko_simple_roots,
    bugaenko_wall_gram,
)
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from coxeter_polytopes.verify.chamber import (
    diagnose_chamber_cone,
    verify_hyperideal_extreme_ray,
)
from scripts.probe_og10_golden_chamber_folding import _interior_point
from scripts import certify_og10_golden_h10_carrier as carrier


ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT / "data/candidates/og10_golden_vinberg_through_3_plus_4phi.json"
EXTENDED_PREFIX = ROOT / "data/candidates/og10_golden_vinberg_through_4_plus_5phi.json"
FIXED_COORDINATES = (0, *range(2, 11))
HYPERIDEAL_WITNESS = (0, 1, 2, 4, 5, 6, 12, 48, 49)
HISTORICAL_HYPERIDEAL_WITNESS = (0, 1, 4, 5, 6, 7, 10, 12, 26)
BUGAENKO_REORDERING = (8, 6, 5, 4, 3, 2, 1, 0, 9, 10, 7)
HISTORICAL_FULL_CENSUS = {
    "finite_ray_count": 784,
    "cusp_count": 0,
    "hyperideal_ray_count": 3113,
    "opposite_sheet_ray_count": 0,
}

Vector = tuple[Golden, ...]
Matrix = tuple[tuple[Golden, ...], ...]


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


def _determinant(matrix: Sequence[Sequence[Golden]]) -> Golden:
    """Compute a small determinant exactly over Q(phi)."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [list(row) for row in matrix]
    result = ONE
    for column in range(size):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return ZERO
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            result = -result
        pivot = work[column][column]
        result *= pivot
        for row in range(column + 1, size):
            factor = work[row][column] / pivot
            for entry in range(column + 1, size):
                work[row][entry] -= factor * work[column][entry]
    return result


def _apply_h(vector: Sequence[Golden]) -> Vector:
    x, y, *tail = vector
    return (x + (ONE - PHI) * y, -y, *tail)


def _load_roots(path: Path, expected_count: int) -> tuple[Vector, ...]:
    payload = json.loads(path.read_text())
    roots = tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )
    if len(roots) != expected_count:
        raise ArithmeticError(
            f"expected {expected_count} archived roots in {path.name}, found {len(roots)}"
        )
    return roots


def _fixed_gram(ambient_gram: Sequence[Sequence[Golden]]) -> Matrix:
    return tuple(
        tuple(ambient_gram[i][j] for j in FIXED_COORDINATES)
        for i in FIXED_COORDINATES
    )


def _restrict(vector: Sequence[Golden]) -> Vector:
    return tuple(vector[i] for i in FIXED_COORDINATES)


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode(
        "ascii"
    )
    return hashlib.sha256(encoded).hexdigest()


def _counter_records(counter: Counter[Golden]) -> list[dict[str, object]]:
    return [
        {"value": value.to_json(), "count": counter[value]}
        for value in sorted(counter, key=lambda item: (item.a, item.b))
    ]


def _global_ultraparallel_lemma_checks() -> dict[str, object]:
    """Return exact scalar checks in the all-simple-root argument.

    If a nonfixed simple pair were not ultraparallel, acuteness and
    ``(r,h r)=2-(3+phi)y^2`` would put ``y^2`` between
    ``2/(3+phi)`` and ``4/(3+phi)``.  Positive definiteness at the conjugate
    embedding gives ``(y')^2 <= 4/(4-phi)``.  Writing ``y=a+b phi``, the
    triangle inequality then forces ``|b|<1``.  Thus b=0, while the first
    upper bound forces the nonzero integer a to have |a|<1, a contradiction.

    The only square roots in the paper proof are eliminated below by
    squaring after checking that both sides are positive.
    """

    identity_lower = Golden(2) / (Golden(3) + PHI)
    identity_upper = Golden(4) / (Golden(3) + PHI)
    conjugate_upper = Golden(4) / (Golden(4) - PHI)
    square_root_comparison_left = Golden(5) - identity_upper - conjugate_upper
    square_root_comparison_margin = (
        square_root_comparison_left * square_root_comparison_left
        - Golden(4) * identity_upper * conjugate_upper
    )
    checks = {
        "identity_lower_positive": identity_lower.sign() > 0,
        "identity_upper_below_one": (ONE - identity_upper).sign() > 0,
        "conjugate_upper_positive": conjugate_upper.sign() > 0,
        "five_minus_upper_bounds_positive": square_root_comparison_left.sign() > 0,
        "sqrt_identity_upper_plus_sqrt_conjugate_upper_below_sqrt5": (
            square_root_comparison_margin.sign() > 0
        ),
    }
    return {
        "pairing_formula": "(r,h(r)) = 2 - (3+phi)*y^2",
        "non_ultraparallel_identity_y_square_lower": identity_lower.to_json(),
        "non_ultraparallel_identity_y_square_upper": identity_upper.to_json(),
        "conjugate_y_square_upper": conjugate_upper.to_json(),
        "sqrt_comparison_positive_left": square_root_comparison_left.to_json(),
        "sqrt_comparison_squared_margin": square_root_comparison_margin.to_json(),
        "checks": checks,
        "verified": all(checks.values()),
    }


def _is_rational_square(value: Fraction) -> bool:
    if value < 0:
        return False
    return (
        isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def _smith_diagonal(matrix: SympyMatrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [int(normal[i, i]) for i in range(min(normal.shape)) if normal[i, i]]


def _sympy_inertia(matrix: SympyMatrix) -> Inertia:
    return inertia(
        [[Golden(int(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]
    )


def _ambient_eigenspace_audit() -> dict[str, object]:
    """Reconstruct the full OG10 h eigensublattices after discriminant gluing."""

    field_columns: list[Vector] = [
        (ONE, ZERO, *(ZERO for _ in range(9))),
        (ONE - PHI, Golden(-1), *(ZERO for _ in range(9))),
    ]
    field_columns.extend(
        tuple(ONE if i == j else ZERO for i in range(11))
        for j in range(2, 11)
    )
    original_trace_action = SympyMatrix.zeros(22)
    for j, column in enumerate(field_columns):
        for i, entry in enumerate(column):
            if not entry.is_integral:
                raise ArithmeticError("h is unexpectedly nonintegral")
            a, b = int(entry.a), int(entry.b)
            original_trace_action[2 * i, 2 * j] = a
            original_trace_action[2 * i + 1, 2 * j] = b
            original_trace_action[2 * i, 2 * j + 1] = b
            original_trace_action[2 * i + 1, 2 * j + 1] = a + b

    full_change = (
        SympyMatrix(carrier.trace_basis_change())
        * carrier.t11_reduction()
        * carrier.core_small_reordering()
    )
    trace_action = full_change.inv() * original_trace_action * full_change
    direct_action = SympyMatrix.diag(trace_action, -eye(2))
    overlattice_change = SympyMatrix.diag(eye(18), carrier.overlattice_basis())
    ambient_action = overlattice_change.inv() * direct_action * overlattice_change
    direct_gram = SympyMatrix.diag(
        carrier.core_gram(), carrier.small_n_gram(), carrier.complement_gram()
    )
    ambient_gram = overlattice_change.T * direct_gram * overlattice_change

    fixed_basis = SympyMatrix.hstack(*(ambient_action - eye(24)).nullspace())
    if not all(entry.is_Integer is True for entry in fixed_basis):
        raise ArithmeticError("the fixed-kernel basis is unexpectedly nonintegral")
    fixed_basis_smith = smith_normal_form(fixed_basis, domain=ZZ)
    fixed_basis_invariants = [
        int(fixed_basis_smith[i, i]) for i in range(fixed_basis.cols)
    ]
    fixed_gram = fixed_basis.T * ambient_gram * fixed_basis

    anti_equations = ambient_action + eye(24)
    domain_equations = DomainMatrix.from_Matrix(
        anti_equations, fmt="dense"
    ).convert_to(ZZ)
    diagonal, _, right_change = smith_normal_decomp(domain_equations)
    diagonal_matrix = diagonal.to_Matrix()
    zero_columns = [
        j
        for j in range(24)
        if all(diagonal_matrix[i, j] == 0 for i in range(24))
    ]
    anti_basis = right_change.to_Matrix()[:, zero_columns]
    anti_gram = anti_basis.T * ambient_gram * anti_basis

    fixed_smith = _smith_diagonal(fixed_gram)
    anti_smith = _smith_diagonal(anti_gram)
    combined_index = abs(int(SympyMatrix.hstack(fixed_basis, anti_basis).det()))
    checks = {
        "ambient_action_integral": all(
            entry.is_Integer is True for entry in ambient_action
        ),
        "ambient_action_is_involution": ambient_action * ambient_action == eye(24),
        "fixed_basis_is_saturated": fixed_basis_invariants == [1] * 20,
        "fixed_and_anti_ranks_are_20_and_4": (
            fixed_basis.cols == 20 and anti_basis.cols == 4
        ),
        "fixed_and_anti_are_orthogonal": (
            fixed_basis.T * ambient_gram * anti_basis
            == SympyMatrix.zeros(20, 4)
        ),
        "fixed_invariants_match_2_elementary_type_1_19_4_1": (
            _sympy_inertia(fixed_gram) == Inertia(1, 19, 0)
            and int(fixed_gram.det()) == -16
            and fixed_smith == [1] * 16 + [2] * 4
        ),
        "anti_invariants_match_2_2_discriminant_48": (
            _sympy_inertia(anti_gram) == Inertia(2, 2, 0)
            and int(anti_gram.det()) == 48
            and anti_smith == [2, 2, 2, 6]
        ),
        "fixed_plus_anti_has_index_16": combined_index == 16,
    }
    return {
        "fixed_lattice": {
            "rank": 20,
            "inertia": [1, 19, 0],
            "determinant": -16,
            "smith_diagonal": fixed_smith,
            "two_elementary_invariants": {"r": 20, "a": 4, "delta": 1},
            "classification_isometry": "<2> + <-2>^3 + E8(-1)^2",
        },
        "anti_invariant_lattice": {
            "rank": 4,
            "inertia": [2, 2, 0],
            "determinant": 48,
            "smith_diagonal": anti_smith,
            "classification_isometry": "<2>^2 + <-2> + <-6>",
        },
        "fixed_plus_anti_index_in_full_H2": combined_index,
        "checks": checks,
        "verified": all(checks.values()),
    }


def _kahler_regularization_audit() -> dict[str, object]:
    """Exclude a Kähler wall whose hyperplane contains the h-fixed space."""

    involution_normal_square = Golden(6, 2)
    wall_norms = (Golden(2), Golden(4), Golden(6, -2))
    ratios = tuple(value / involution_normal_square for value in wall_norms)
    rational_norms = tuple(ratio.norm() for ratio in ratios)
    checks = {
        "wall_ratio_norms_are_1_11_4_11_5_11": rational_norms
        == (Fraction(1, 11), Fraction(4, 11), Fraction(5, 11)),
        "no_wall_ratio_norm_is_a_rational_square": all(
            not _is_rational_square(value) for value in rational_norms
        ),
    }
    return {
        "visible_kahler_wall_norm_representatives": [
            value.to_json() for value in wall_norms
        ],
        "ratios_to_h_reflection_normal_square": [
            value.to_json() for value in ratios
        ],
        "rational_norms_of_ratios": [str(value) for value in rational_norms],
        "conclusion": (
            "No visible OG10 Kahler-wall normal is proportional over Q(sqrt(5)) "
            "to the h-anti-invariant line, so no such wall contains Fix(h)."
        ),
        "checks": checks,
        "verified": all(checks.values()),
    }
def build(*, full_census: bool = False) -> dict[str, object]:
    ambient_gram = og10_field_gram()
    fixed_gram = _fixed_gram(ambient_gram)
    archived = _load_roots(PREFIX, 457)
    extended = _load_roots(EXTENDED_PREFIX, 1144)
    closed = tuple(dict.fromkeys((*archived, *(_apply_h(root) for root in archived))))
    historical_fixed = tuple(root for root in closed if _apply_h(root) == root)
    historical_restricted = tuple(_restrict(root) for root in historical_fixed)
    fixed = tuple(root for root in extended if _apply_h(root) == root)
    restricted = tuple(_restrict(root) for root in fixed)

    visited: set[Vector] = set()
    paired_orbits: list[tuple[Vector, Vector]] = []
    for root in closed:
        if root in visited:
            continue
        image = _apply_h(root)
        visited.add(root)
        visited.add(image)
        if image != root:
            paired_orbits.append((root, image))

    pairings = Counter(
        _pair(root, ambient_gram, image) for root, image in paired_orbits
    )
    sum_squares = Counter(
        _pair(
            tuple(root[i] + image[i] for i in range(11)),
            ambient_gram,
            tuple(root[i] + image[i] for i in range(11)),
        )
        for root, image in paired_orbits
    )
    levels = Counter(root[0] for root in fixed)

    _, _, interior = _interior_point(ambient_gram)
    image_interior = _apply_h(interior)
    invariant_interior = tuple(
        (interior[i] + image_interior[i]) / 2 for i in range(11)
    )
    restricted_interior = _restrict(invariant_interior)

    witness = verify_hyperideal_extreme_ray(
        restricted,
        fixed_gram,
        HYPERIDEAL_WITNESS,
        dimension=9,
    )
    historical_witness = verify_hyperideal_extreme_ray(
        historical_restricted,
        fixed_gram,
        HISTORICAL_HYPERIDEAL_WITNESS,
        dimension=9,
    )

    u_zero_roots = tuple(root for root in restricted if not root[-1])
    reordered_u_zero_gram = [
        [
            _pair(
                u_zero_roots[BUGAENKO_REORDERING[i]],
                fixed_gram,
                u_zero_roots[BUGAENKO_REORDERING[j]],
            )
            for j in range(11)
        ]
        for i in range(11)
    ]
    reordered_u_zero_roots = tuple(
        tuple(u_zero_roots[index][:-1]) for index in BUGAENKO_REORDERING
    )
    u_zero_basis_determinant = _determinant(reordered_u_zero_roots[:9])
    archived_bugaenko_basis_determinant = _determinant(
        tuple(tuple(root) for root in bugaenko_simple_roots()[:9])
    )
    ambient_eigenspaces = _ambient_eigenspace_audit()
    regularization = _kahler_regularization_audit()
    nonfixed_formula_checks = []
    nu = (PHI - ONE, Golden(2), *(ZERO for _ in range(9)))
    for root, image in paired_orbits:
        y = root[1]
        difference = tuple(image[i] - root[i] for i in range(11))
        expected_difference = tuple(-y * entry for entry in nu)
        pairing = _pair(root, ambient_gram, image)
        projected_square = _pair(
            tuple(root[i] + image[i] for i in range(11)),
            ambient_gram,
            tuple(root[i] + image[i] for i in range(11)),
        )
        nonfixed_formula_checks.append(
            difference == expected_difference
            and pairing == Golden(2) - (Golden(3) + PHI) * y * y
            and projected_square == Golden(8) - Golden(2) * (Golden(3) + PHI) * y * y
        )

    checks = {
        "ambient_archive_has_457_roots": len(archived) == 457,
        "ambient_extended_prefix_has_1144_roots": len(extended) == 1144,
        "h_closed_wall_set_has_640_roots": len(closed) == 640,
        "old_h_closed_orbits_are_44_fixed_plus_298_pairs": (
            sum(_apply_h(root) == root for root in closed) == 44
            and len(paired_orbits) == 298
        ),
        "extended_prefix_has_61_fixed_roots": len(fixed) == 61,
        "fixed_roots_are_exactly_y_zero": all(not root[1] for root in fixed),
        "fixed_lattice_identity_signature_is_9_1": (
            inertia(fixed_gram, embedding=0) == Inertia(9, 1, 0)
        ),
        "fixed_lattice_conjugate_signature_is_10_0": (
            inertia(fixed_gram, embedding=1) == Inertia(10, 0, 0)
        ),
        "invariant_interior_is_fixed": _apply_h(invariant_interior) == invariant_interior,
        "invariant_interior_is_timelike": (
            _pair(invariant_interior, ambient_gram, invariant_interior).sign() < 0
        ),
        "invariant_interior_is_strict_for_61_walls": all(
            _pair(invariant_interior, ambient_gram, root).sign() < 0
            for root in fixed
        ),
        "all_298_archived_pair_formulas_hold": all(nonfixed_formula_checks),
        "all_298_archived_projected_normals_are_timelike": all(
            value.sign() < 0 for value in sum_squares
        ),
        "global_ultraparallel_scalar_inequalities_hold": (
            _global_ultraparallel_lemma_checks()["verified"] is True
        ),
        "61_wall_outer_cone_has_certified_hyperideal_extreme_ray": (
            witness.norm.sign() > 0
        ),
        "historical_44_wall_witness_recomputed": historical_witness.norm
        == Golden(178, -110),
        "u_zero_fixed_walls_have_bugaenko_11_wall_gram": (
            len(u_zero_roots) == 11
            and reordered_u_zero_gram == bugaenko_wall_gram()
        ),
        "u_zero_bugaenko_first_nine_roots_are_an_integral_lattice_basis": (
            u_zero_basis_determinant == Golden(-3, 2)
            and archived_bugaenko_basis_determinant == Golden(-3, 2)
            and u_zero_basis_determinant.norm() == -1
        ),
        "full_H2_eigensublattices_verified": ambient_eigenspaces["verified"]
        is True,
        "fixed_slice_not_contained_in_a_visible_kahler_wall": regularization[
            "verified"
        ]
        is True,
    }

    result: dict[str, object] = {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_OG10_golden_h9_fixed_slice_prefix_diagnosis",
        "candidate_id": "og10_golden_h9_fixed_lattice",
        "verified": all(checks.values()),
        "checks": checks,
        "involution": {
            "field_matrix_first_binary_block_row_major": [
                [ONE.to_json(), (ONE - PHI).to_json()],
                [ZERO.to_json(), Golden(-1).to_json()],
            ],
            "reflection_normal": _encoded_vector(nu),
            "reflection_normal_square": _pair(nu, ambient_gram, nu).to_json(),
            "fixed_vector_condition": "y=0",
        },
        "fixed_lattice": {
            "decomposition": "<-2*phi> + E8 + <2>",
            "rank_over_Z_phi": 10,
            "hyperbolic_dimension": 9,
            "gram": [[entry.to_json() for entry in row] for row in fixed_gram],
            "identity_inertia": [9, 1, 0],
            "conjugate_inertia": [10, 0, 0],
        },
        "global_fixed_slice_theorem": {
            "difference_formula": "h(r)-r = -y*(phi-1,2,0^9)",
            "projected_normal_square_formula": "(r+h(r))^2 = 8-2*(3+phi)*y^2",
            "nonfixed_simple_orbits": "globally_ultraparallel_and_redundant_on_H9",
            "fixed_simple_roots": "exactly_the_simple_norm_two_walls_of_the_fixed_lattice_chamber",
            "arithmetic_lemma": _global_ultraparallel_lemma_checks(),
        },
        "right_angled_bugaenko_facet": {
            "facet_normal_in_fixed_coordinates": _encoded_vector(
                (ZERO, *(ZERO for _ in range(8)), Golden(-1))
            ),
            "facet_normal_square": Golden(2).to_json(),
            "global_simple_root_dichotomy": (
                "For r=(x,z,u), acuteness gives u>=0 and conjugate "
                "positivity gives |u'|<=1.  Integrality excludes 0<u<1; "
                "hence u=0 or (s_u,r)=-2u<=-2."
            ),
            "complete_facet": "Bugaenko compact H8 chamber",
            "complete_facet_wall_count": 11,
            "complete_facet_f_vector": [41, 164, 316, 374, 294, 156, 54, 11, 1],
            "u_zero_wall_indices_in_61_wall_prefix": [
                i for i, root in enumerate(restricted) if not root[-1]
            ],
            "bugaenko_diagram_reordering": list(BUGAENKO_REORDERING),
            "verified_diagram_gram": reordered_u_zero_gram
            == bugaenko_wall_gram(),
            "first_nine_root_coordinate_determinant": (
                u_zero_basis_determinant.to_json()
            ),
            "first_nine_root_determinant_norm": str(
                u_zero_basis_determinant.norm()
            ),
            "verified_integral_lattice_identification": (
                u_zero_basis_determinant == archived_bugaenko_basis_determinant
                and u_zero_basis_determinant.norm() == -1
            ),
            "scope": (
                "This proves a compact codimension-one facet of the complete H9 "
                "chamber, not finite volume of the H9 chamber."
            ),
        },
        "algebro_geometric_involution": {
            "full_H2_eigensublattices": ambient_eigenspaces,
            "kahler_wall_containment_audit": regularization,
            "external_chamber_theorem_consequence": (
                "Choose a generic point of Fix(h) in the eigen-movable chamber, "
                "away from the locally finite Kahler-wall restrictions.  Its "
                "Kahler-type chamber is h-stable and is realized by a same-period "
                "birational model.  Since the positive Hodge-monodromy extension "
                "acts by -1 on H^(2,0) and preserves that model's Kahler cone, "
                "global Torelli realizes h as a biregular anti-symplectic involution."
            ),
            "scope": (
                "The H9 Coxeter chamber is the golden identity-embedding slice "
                "inside the rank-20 invariant Neron-Severi lattice; it is not the "
                "full H19 positive cone of that invariant lattice."
            ),
        },
        "known_wall_orbits": {
            "archived_roots": len(archived),
            "extended_prefix_roots": len(extended),
            "h_closed_roots": len(closed),
            "old_h_closed_orbits": 44 + len(paired_orbits),
            "old_h_closed_fixed_roots": 44,
            "old_h_closed_paired_orbits": len(paired_orbits),
            "extended_prefix_fixed_roots": len(fixed),
            "fixed_root_original_1144_indices": [
                i for i, root in enumerate(extended) if _apply_h(root) == root
            ],
            "fixed_roots_sha256": _digest(
                [_encoded_vector(root) for root in restricted]
            ),
            "fixed_root_control_level_counts": _counter_records(levels),
            "paired_root_inner_product_counts": _counter_records(pairings),
            "paired_projected_sum_square_counts": _counter_records(sum_squares),
        },
        "invariant_interior_point": _encoded_vector(restricted_interior),
        "known_61_wall_outer_cone": {
            "status": "not_finite_volume_before_any_later_fixed_walls_are_added",
            "facet_count": len(restricted),
            "hyperideal_witness_active_facets": list(witness.active_facets),
            "hyperideal_witness_ray": _encoded_vector(witness.ray),
            "hyperideal_witness_norm": witness.norm.to_json(),
            "hyperideal_witness_pairings_sha256": _digest(
                [entry.to_json() for entry in witness.wall_pairings]
            ),
            "scope": (
                "This obstructs finite volume for the present 61-wall outer cone only; "
                "unseen fixed simple roots may cut off the ray."
            ),
        },
        "historical_44_wall_outer_cone": {
            "facet_count": len(historical_restricted),
            "hyperideal_witness_active_facets": list(
                historical_witness.active_facets
            ),
            "hyperideal_witness_ray": _encoded_vector(historical_witness.ray),
            "hyperideal_witness_norm": historical_witness.norm.to_json(),
            "full_exact_census": HISTORICAL_FULL_CENSUS,
            "full_exact_census_recomputed_in_this_run": False,
        },
        "conclusion": (
            "The complete h-fixed slice is a norm-two Coxeter chamber for "
            "<-2*phi>+E8+<2>.  Among the currently known walls it has 61 fixed "
            "facets, whose outer cone has a certified hyperideal escape.  Its "
            "u=0 facet is exactly Bugaenko's compact H8 polytope.  Completeness "
            "and finite covolume of the fixed-lattice chamber remain open."
        ),
    }

    if full_census:
        geometry = diagnose_chamber_cone(
            historical_restricted,
            fixed_gram,
            restricted_interior,
            dimension=9,
        )
        recomputed = {
            "finite_ray_count": len(geometry.finite_vertices),
            "cusp_count": len(geometry.cusps),
            "hyperideal_ray_count": len(geometry.hyperideal_rays),
            "opposite_sheet_ray_count": len(geometry.opposite_sheet_rays),
        }
        if recomputed != HISTORICAL_FULL_CENSUS:
            raise ArithmeticError("historical 44-wall full census changed")
        result["historical_44_wall_outer_cone"]["full_exact_census"] = recomputed
        result["historical_44_wall_outer_cone"][
            "full_exact_census_recomputed_in_this_run"
        ] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full-census",
        action="store_true",
        help="recompute all extreme rays of the historical 44-wall outer cone",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write the canonical JSON certificate to this path instead of stdout",
    )
    args = parser.parse_args()
    encoded = json.dumps(
        build(full_census=args.full_census), indent=2, sort_keys=True
    ) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded)


if __name__ == "__main__":
    main()
