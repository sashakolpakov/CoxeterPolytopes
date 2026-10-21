#!/usr/bin/env python3
r"""Certify an infinite stabilizer of the complete silver norm-two chamber.

This is deliberately not a Vinberg stopping argument.  Two integral
involutions of the rank-ten ``Z[sqrt(2)]`` carrier are shown to fix regular
timelike points of the *complete* arrangement chamber.  Completeness at those
points follows from an exact Lorentzian separator bound which lies inside the
already exhaustive part of the archived Vinberg run.  Their product has
trace greater than ten and therefore infinite order.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    _silver_dot,
    _silver_inverse,
    _silver_mat_vec,
    _silver_quadratic,
    admissible_silver_control_coordinates,
    audit_silver_k3_height_zero,
    silver_k3_alexeev_field_gram,
    silver_k3_control_vector,
    silver_k3_height_zero_walls,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.json"
)
ARTIFACT = ROOT / "data/candidates/silver_k3_norm2_infinite_stabilizer.json"
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "01c8c722b0db00f757d67789d82674f2e9724ec2e81e4b7e720870e1700a67f5"
)
EXPECTED_INPUT_ROOTS_SHA256 = (
    "26605ff267008faf5895f216acb8cd7c7a420ab136021b97ba476e9a1572ad33"
)

Vector = tuple[Silver, ...]
Matrix = tuple[Vector, ...]


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _identity(size: int) -> Matrix:
    return tuple(
        tuple(Silver(1 if row == column else 0) for column in range(size))
        for row in range(size)
    )


def _compose(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(
                (left[row][middle] * right[middle][column] for middle in range(10)),
                SILVER_ZERO,
            )
            for column in range(10)
        )
        for row in range(10)
    )


def _transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][column] for row in range(10)) for column in range(10))


def _matrix_isometry(matrix: Matrix, gram: Matrix) -> bool:
    return _compose(_transpose(matrix), _compose(gram, matrix)) == gram


def _matrix_digest(matrix: Matrix) -> str:
    return _canonical_digest([[entry.to_json() for entry in row] for row in matrix])


def _ldlt_pivots(matrix: Matrix) -> tuple[Silver, ...]:
    """Return exact diagonal pivots for an unpivoted ``LDL^T`` decomposition."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("the field form must be square")
    if any(matrix[i][j] != matrix[j][i] for i in range(size) for j in range(i)):
        raise ValueError("the field form must be symmetric")
    work = [list(row) for row in matrix]
    pivots: list[Silver] = []
    for index in range(size):
        pivot = work[index][index]
        if not pivot:
            raise ValueError("the advertised exact LDL decomposition has a zero pivot")
        pivots.append(pivot)
        tail = range(index + 1, size)
        updated = {
            (row, column): (
                work[row][column]
                - work[row][index] * work[index][column] / pivot
            )
            for row in tail
            for column in tail
        }
        for (row, column), value in updated.items():
            work[row][column] = value
    return tuple(pivots)


def _signature_from_pivots(
    pivots: Sequence[Silver], *, embedding: int
) -> tuple[int, int, int]:
    signs = tuple(pivot.sign(embedding) for pivot in pivots)
    return (signs.count(1), signs.count(-1), signs.count(0))


def _decode_roots(document: Mapping[str, object]) -> tuple[Vector, ...]:
    rows = document.get("ordered_prefix_roots")
    if not isinstance(rows, list):
        raise ValueError("the silver input has no literal ordered roots")
    roots = []
    for row in rows:
        if not isinstance(row, list) or len(row) != 10:
            raise ValueError("a silver input root has the wrong rank")
        roots.append(tuple(Silver.from_json(entry) for entry in row))
    return tuple(roots)


def _columns(vectors: Sequence[Vector]) -> Matrix:
    if len(vectors) != 10:
        raise ValueError("ten columns are required")
    return tuple(tuple(vectors[column][row] for column in range(10)) for row in range(10))


def _seed_reversal(roots: Sequence[Vector]) -> Matrix:
    permutation = (1, 0, 9, 3, 4, 5, 6, 7, 8, 2)
    basis = _columns(roots[:10])
    images = _columns(tuple(roots[permutation[index]] for index in range(10)))
    return _compose(images, _silver_inverse(basis))


def _wall_swap(roots: Sequence[Vector], gram: Matrix) -> tuple[Matrix, Vector, Silver]:
    difference = tuple(roots[4][index] - roots[18][index] for index in range(10))
    gram_difference = _silver_mat_vec(gram, difference)
    difference_norm = _silver_dot(difference, gram_difference)
    coefficient = Silver(2) / difference_norm
    identity = _identity(10)
    matrix = tuple(
        tuple(
            identity[row][column]
            - coefficient * difference[row] * gram_difference[column]
            for column in range(10)
        )
        for row in range(10)
    )
    return matrix, difference, coefficient


def _seed_fixed_witness() -> Vector:
    s = SQRT2
    return (
        Silver(1),
        Silver(147) / 71 - Silver(30) / 71 * s,
        Silver(5) / 71 + Silver(112) / 71 * s,
        Silver(1) + 2 * s,
        Silver(141) / 71 + Silver(77) / 71 * s,
        Silver(54) / 71 + Silver(31) / 71 * s,
        Silver(96) / 71 + Silver(63) / 71 * s,
        Silver(96) / 71 + Silver(63) / 71 * s,
        Silver(54) / 71 + Silver(31) / 71 * s,
        Silver(6) / 71 - s / 142,
    )


def _swap_boundary_minimizer() -> Vector:
    s = SQRT2
    return (
        Silver(1),
        Silver(1),
        s,
        Silver(1) + s,
        Silver(1) + s,
        (Silver(2) + 2 * s) / 3,
        Silver(1) + s,
        (Silver(4) + 2 * s) / 3,
        (Silver(2) + s) / 3,
        (Silver(2) + s) / 6,
    )


def _swap_strict_point() -> Vector:
    s = SQRT2
    return (
        Silver(1),
        Silver(856) / 713 + Silver(45) / 713 * s,
        Silver(166) / 713 + Silver(850) / 713 * s,
        Silver(41) / 31 + Silver(40) / 31 * s,
        Silver(955) / 713 + Silver(844) / 713 * s,
        Silver(486) / 713 + Silver(487) / 713 * s,
        Silver(882) / 713 + Silver(831) / 713 * s,
        Silver(1140) / 713 + Silver(623) / 713 * s,
        Silver(586) / 713 + Silver(329) / 713 * s,
        Silver(226) / 713 + Silver(227) / 1426 * s,
    )


def _separator_bound_squared(
    control: Vector, witness: Vector, gram: Matrix
) -> Silver:
    """Bound ``x^2`` for a norm-two wall separating control and witness.

    Put ``P=control^2<0``, ``V=witness^2<0``, and ``C=(control,witness)``.
    This includes the equality case ``(witness,r)=0``.  Orthogonal
    decomposition relative to ``control`` and Cauchy--Schwarz in its
    positive-definite orthogonal complement give

        x^2 <= 2 * (C^2 + (-P)*V) / (-V).
    """

    control_norm = _silver_quadratic(control, gram)
    witness_norm = _silver_quadratic(witness, gram)
    control_witness = _silver_dot(control, _silver_mat_vec(gram, witness))
    return (
        2
        * (control_witness * control_witness + (-control_norm) * witness_norm)
        / (-witness_norm)
    )


def _charpoly_coefficients(matrix: Matrix) -> tuple[Silver, ...]:
    """Faddeev--LeVerrier coefficients, leading coefficient first."""

    identity = _identity(10)
    work = identity
    coefficients = [Silver(1)]
    for degree in range(1, 11):
        product = _compose(matrix, work)
        coefficient = -sum((product[i][i] for i in range(10)), SILVER_ZERO) / degree
        coefficients.append(coefficient)
        work = tuple(
            tuple(
                product[row][column]
                + (coefficient if row == column else SILVER_ZERO)
                for column in range(10)
            )
            for row in range(10)
        )
    return tuple(coefficients)


def _encoded_vector(vector: Sequence[Silver]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def build() -> dict[str, object]:
    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    if input_document.get("certificate_sha256") != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise ValueError("the frozen 1,291-wall input certificate changed")
    if input_document.get("ordered_prefix_roots_sha256") != EXPECTED_INPUT_ROOTS_SHA256:
        raise ValueError("the frozen 1,291-wall root ordering changed")
    roots = _decode_roots(input_document)
    if len(roots) != 1291:
        raise ValueError("the frozen silver prefix no longer has 1,291 walls")

    gram = silver_k3_alexeev_field_gram()
    control = silver_k3_control_vector()
    identity = _identity(10)
    root_images = tuple(_silver_mat_vec(gram, root) for root in roots)
    form_pivots = _ldlt_pivots(gram)
    identity_signature = _signature_from_pivots(form_pivots, embedding=0)
    conjugate_signature = _signature_from_pivots(form_pivots, embedding=1)
    form_determinant = Silver(1)
    for pivot in form_pivots:
        form_determinant *= pivot
    field_anisotropic = conjugate_signature == (10, 0, 0)
    height_zero = audit_silver_k3_height_zero()
    height_zero_walls = silver_k3_height_zero_walls()

    seed = _seed_reversal(roots)
    seed_witness = _seed_fixed_witness()
    seed_pairings = tuple(_silver_dot(seed_witness, image) for image in root_images)
    seed_bound = _separator_bound_squared(control, seed_witness, gram)

    swap, difference, swap_coefficient = _wall_swap(roots, gram)
    swap_boundary = _swap_boundary_minimizer()
    swap_strict = _swap_strict_point()
    swap_witness = tuple(
        (swap_boundary[index] + swap_strict[index]) / 2 for index in range(10)
    )
    swap_pairings = tuple(_silver_dot(swap_witness, image) for image in root_images)
    swap_boundary_pairings = tuple(
        _silver_dot(swap_boundary, image) for image in root_images
    )
    swap_strict_pairings = tuple(
        _silver_dot(swap_strict, image) for image in root_images
    )
    swap_bound = _separator_bound_squared(control, swap_witness, gram)

    product = _compose(seed, swap)
    product_trace = sum((product[index][index] for index in range(10)), SILVER_ZERO)
    charpoly = _charpoly_coefficients(product)
    expected_charpoly = (
        Silver(1),
        Silver(-10, -4),
        Silver(39, 24),
        Silver(-76, -56),
        Silver(70, 56),
        Silver(0),
        Silver(-70, -56),
        Silver(76, 56),
        Silver(-39, -24),
        Silver(10, 4),
        Silver(-1),
    )

    expected_short_coordinates = tuple(admissible_silver_control_coordinates(4))
    expected_swap_coordinates = tuple(admissible_silver_control_coordinates(3))
    archived_coordinates = tuple(
        Silver.from_json(entry) for entry in input_document["candidate_coordinates"]
    )
    expected_short_json = [entry.to_json() for entry in expected_short_coordinates]

    control_seed = _silver_dot(control, _silver_mat_vec(gram, seed_witness))
    control_swap = _silver_dot(control, _silver_mat_vec(gram, swap_witness))
    seed_isometry = _matrix_isometry(seed, gram)
    swap_isometry = _matrix_isometry(swap, gram)
    seed_fixed_future_timelike = (
        _silver_mat_vec(seed, seed_witness) == seed_witness
        and _silver_quadratic(seed_witness, gram).sign(0) < 0
        and control_seed.sign(0) < 0
    )
    swap_fixed_future_timelike = (
        _silver_mat_vec(swap, swap_witness) == swap_witness
        and _silver_quadratic(swap_witness, gram).sign(0) < 0
        and control_swap.sign(0) < 0
    )
    seed_preserves_future_cone = seed_isometry and seed_fixed_future_timelike
    swap_preserves_future_cone = swap_isometry and swap_fixed_future_timelike
    checks = {
        "input_has_literal_complete_prefix_through_6_plus_3sqrt2": (
            input_document.get("completed_through_control_coordinate")
            == Silver(6, 3).to_json()
            and input_document.get("next_admissible_control_coordinate")
            == Silver(5, 4).to_json()
            and archived_coordinates[: len(expected_short_coordinates)]
            == expected_short_coordinates
        ),
        "all_1291_roots_have_norm_two": all(
            _silver_quadratic(root, gram) == 2 for root in roots
        ),
        "field_form_rank_signatures_and_determinant_recomputed_from_gram": (
            len(form_pivots) == 10
            and all(form_pivots)
            and identity_signature == (9, 1, 0)
            and conjugate_signature == (10, 0, 0)
            and form_determinant == -2 - 2 * SQRT2
        ),
        "conjugate_positive_definiteness_implies_field_anisotropy": (
            field_anisotropic
        ),
        "height_zero_norm_two_root_system_is_exhaustive": (
            height_zero.passes
            and height_zero.reflection_closure_matches_exhaustion
            and height_zero.total_field_norm_two_roots == 82
            and roots[:9] == height_zero_walls
            and all(not root[0] for root in roots[:9])
            and all(root[0] for root in roots[9:])
        ),
        "seed_reversal_is_O_integral_isometric_involution": (
            all(entry.is_integral for row in seed for entry in row)
            and seed_isometry
            and _compose(seed, seed) == identity
        ),
        "seed_reversal_has_a_regular_fixed_future_timelike_point": (
            seed_fixed_future_timelike
            and all(pairing.sign(0) < 0 for pairing in seed_pairings)
        ),
        "seed_reversal_preserves_the_chosen_future_cone": (
            seed_preserves_future_cone
        ),
        "seed_separator_bound_is_inside_exhaustive_x_less_than_4_range": (
            seed_bound.sign(0) > 0
            and (Silver(16) - seed_bound).sign(0) > 0
            and expected_short_coordinates
            == (Silver(1), SQRT2, Silver(2), Silver(1) + SQRT2,
                Silver(2) + SQRT2, Silver(1) + 2 * SQRT2)
        ),
        "wall_swap_is_O_integral_isometric_involution": (
            _silver_dot(roots[4], root_images[18]) == -SQRT2
            and _silver_quadratic(difference, gram) == 4 + 2 * SQRT2
            and swap_coefficient == 1 - SQRT2 / 2
            and all(entry.is_integral for row in swap for entry in row)
            and swap_isometry
            and _compose(swap, swap) == identity
            and _silver_mat_vec(swap, roots[4]) == roots[18]
            and _silver_mat_vec(swap, roots[18]) == roots[4]
        ),
        "wall_swap_boundary_minimizer_is_fixed_with_advertised_zero_set": (
            _silver_mat_vec(swap, swap_boundary) == swap_boundary
            and _silver_quadratic(swap_boundary, gram) == -SQRT2 / 3
            and all(pairing.sign(0) <= 0 for pairing in swap_boundary_pairings)
            and tuple(
                index
                for index, pairing in enumerate(swap_boundary_pairings)
                if not pairing
            ) == (0, 1, 2, 3, 5, 6, 7)
            and _silver_mat_vec(swap, swap_strict) == swap_strict
            and all(pairing.sign(0) < 0 for pairing in swap_strict_pairings)
        ),
        "wall_swap_has_a_regular_fixed_future_timelike_point": (
            swap_fixed_future_timelike
            and all(pairing.sign(0) < 0 for pairing in swap_pairings)
        ),
        "wall_swap_preserves_the_chosen_future_cone": (
            swap_preserves_future_cone
        ),
        "wall_swap_separator_bound_is_inside_exhaustive_x_less_than_3_range": (
            swap_bound.sign(0) > 0
            and (Silver(9) - swap_bound).sign(0) > 0
            and archived_coordinates[: len(expected_swap_coordinates)]
            == expected_swap_coordinates
        ),
        "infinite_order_product_has_exact_trace_and_characteristic_polynomial": (
            product_trace == 10 + 4 * SQRT2
            and (product_trace - 10).sign(0) > 0
            and charpoly == expected_charpoly
        ),
    }

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_complete_norm_two_chamber_infinite_stabilizer",
        "verified": all(checks.values()),
        "input": {
            "manifest": str(INPUT.relative_to(ROOT)),
            "certificate_sha256": EXPECTED_INPUT_CERTIFICATE_SHA256,
            "ordered_roots_sha256": EXPECTED_INPUT_ROOTS_SHA256,
            "wall_count": len(roots),
            "completed_through_control_coordinate": Silver(6, 3).to_json(),
            "all_admissible_coordinates_below_4": expected_short_json,
        },
        "lattice": {
            "ring": "Z[sqrt(2)]",
            "rank": 10,
            "hyperbolic_dimension": 9,
            "field_gram_ldlt_pivots": [pivot.to_json() for pivot in form_pivots],
            "field_determinant": form_determinant.to_json(),
            "field_rank": len(form_pivots),
            "identity_inertia_positive_negative_zero": list(identity_signature),
            "conjugate_inertia_positive_negative_zero": list(conjugate_signature),
            "anisotropy_derivation": (
                "If 0 != z in Q(sqrt(2))^10 had h(z,z)=0, applying the "
                "nonidentity embedding would give a nonzero isotropic vector "
                "for the certified positive-definite conjugate form. Therefore "
                "h is anisotropic over Q(sqrt(2))."
            ),
            "control_norm": _silver_quadratic(control, gram).to_json(),
            "control_covector": [Silver(-1).to_json(), *[Silver(0).to_json() for _ in range(9)]],
        },
        "height_zero": {
            "finite_type": height_zero.finite_type,
            "trace_short_vector_counts": list(height_zero.trace_short_vector_counts),
            "field_norm_two_counts": list(height_zero.field_norm_two_counts),
            "total_field_norm_two_roots": height_zero.total_field_norm_two_roots,
            "reflection_closure_count": height_zero.reflection_closure_count,
            "reflection_closure_matches_exhaustion": (
                height_zero.reflection_closure_matches_exhaustion
            ),
            "simple_wall_indices": list(range(9)),
            "role_in_completeness_proof": (
                "The control point lies on these nine mirrors. Exhaustion of all "
                "82 height-zero roots proves that the strict negative direction "
                "of either witness enters the chosen finite Weyl sector; a later "
                "first-crossed mirror must therefore have x>0."
            ),
        },
        "seed_C5_reversal": {
            "basis_wall_indices": list(range(10)),
            "basis_wall_permutation": [1, 0, 9, 3, 4, 5, 6, 7, 8, 2],
            "matrix_sha256": _matrix_digest(seed),
            "trace": sum((seed[i][i] for i in range(10)), SILVER_ZERO).to_json(),
            "fixed_witness": _encoded_vector(seed_witness),
            "fixed_witness_norm": _silver_quadratic(seed_witness, gram).to_json(),
            "control_pairing": control_seed.to_json(),
            "largest_archived_wall_pairing": max(seed_pairings).to_json(),
            "largest_pairing_wall_indices": [
                index for index, value in enumerate(seed_pairings) if value == max(seed_pairings)
            ],
            "separator_control_coordinate_bound_squared": seed_bound.to_json(),
            "separator_bound_less_than": 4,
            "future_cone_preserved": seed_preserves_future_cone,
            "separator_hypothesis": "(fixed_witness,r)>=0, including equality",
            "complete_chamber_argument": (
                "Height zero is exhaustive. Starting into its chosen Weyl sector, "
                "the first mirror on a path from the control point to the fixed "
                "witness would be a Vinberg simple wall. The equality-inclusive "
                "separator lemma bounds that wall by x<4. Every admissible positive "
                "x<4 was exhausted in the frozen Vinberg chain, and the witness is "
                "strict on all 1291 accepted walls. Thus it is regular in the complete "
                "chamber. The isometry fixes it and preserves the future cone, so it "
                "stabilizes that complete chamber."
            ),
        },
        "wall_4_18_swap": {
            "wall_indices": [4, 18],
            "wall_pairing": _silver_dot(roots[4], root_images[18]).to_json(),
            "difference_norm": _silver_quadratic(difference, gram).to_json(),
            "reflection_coefficient": swap_coefficient.to_json(),
            "matrix_sha256": _matrix_digest(swap),
            "boundary_minimizer": _encoded_vector(swap_boundary),
            "boundary_minimizer_norm": _silver_quadratic(swap_boundary, gram).to_json(),
            "boundary_minimizer_zero_wall_indices": [0, 1, 2, 3, 5, 6, 7],
            "fixed_witness_construction": "midpoint(boundary_minimizer,strict_fixed_point)",
            "fixed_witness": _encoded_vector(swap_witness),
            "fixed_witness_norm": _silver_quadratic(swap_witness, gram).to_json(),
            "control_pairing": control_swap.to_json(),
            "largest_archived_wall_pairing": max(swap_pairings).to_json(),
            "largest_pairing_wall_indices": [
                index for index, value in enumerate(swap_pairings) if value == max(swap_pairings)
            ],
            "separator_control_coordinate_bound_squared": swap_bound.to_json(),
            "separator_bound_less_than": 3,
            "future_cone_preserved": swap_preserves_future_cone,
            "separator_hypothesis": "(fixed_witness,r)>=0, including equality",
            "complete_chamber_argument": (
                "After the exhaustive height-zero sector, any first-crossed mirror "
                "would be a positive-height Vinberg simple wall. The equality-inclusive "
                "separator lemma gives x<3, entirely inside the exhaustive Vinberg "
                "range. The midpoint is strict on every archived wall, hence regular "
                "in the complete chamber. The swap fixes it and preserves the future "
                "cone, so it stabilizes the complete chamber rather than merely the "
                "1291-wall outer cone."
            ),
        },
        "infinite_order_product": {
            "order": "seed_C5_reversal_times_wall_4_18_swap",
            "matrix_sha256": _matrix_digest(product),
            "trace": product_trace.to_json(),
            "trace_minus_dimension": (product_trace - 10).to_json(),
            "characteristic_polynomial_factorization": (
                "(t-1)^7*(t+1)*(t^2-(4+4sqrt(2))*t+1)"
            ),
            "characteristic_polynomial_coefficients": [entry.to_json() for entry in charpoly],
            "proof": (
                "A finite-order real 10 by 10 matrix has all eigenvalues on the unit "
                "circle, hence absolute trace at most 10. Here the identity-embedding "
                "trace is 10+4sqrt(2)>10."
            ),
        },
        "arithmetic_conclusion": {
            "anisotropic_over_Q_sqrt2_derived_from_conjugate_definiteness": (
                field_anisotropic
            ),
            "arithmetic_group": "Gamma=O^+(L,h)",
            "arithmetic_quotient": "compact (uniform lattice)",
            "reflection_group": "W=<s_r : r^2=2>",
            "normality_reason": "Gamma permutes the complete norm-two root set",
            "stabilizer_quotient": "Gamma/W is isomorphic to Stab_Gamma(C)",
            "complete_chamber_stabilizer": "infinite",
            "reflection_group_index": "infinite",
            "complete_chamber_hyperbolic_volume": "infinite",
            "standard_inputs": (
                "Borel-Harish-Chandra plus the anisotropy compactness criterion; "
                "Vinberg's standard chamber/simple-transitivity lemmas for reflection groups."
            ),
            "citation_dependent_inputs_not_reproved_by_this_certificate": [
                "Borel-Harish-Chandra: the integral orthogonal group is an arithmetic lattice.",
                "Godement compactness criterion: field anisotropy makes that lattice uniform.",
                "Vinberg reflection theory: the locally finite norm-two arrangement has W acting simply transitively on chambers.",
                "For W normal in Gamma, Gamma=W semidirect Stab_Gamma(C), and covering volume multiplies by [Gamma:W].",
            ],
            "computational_certificate_scope": (
                "The verified flag covers the exact lattice, root, symmetry, witness, "
                "separator, height-zero, and infinite-order computations. It invokes, "
                "but does not computationally reprove, the listed arithmetic and "
                "reflection-group theorems."
            ),
            "scope": (
                "This decides the complete intrinsic norm-two silver chamber negatively. "
                "It does not classify other reflective norms or prove a general H9 "
                "non-existence theorem."
            ),
        },
        "checks": checks,
    }
    payload["certificate_sha256"] = _canonical_digest(payload)
    return payload


def verify_frozen() -> dict[str, object]:
    expected = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    actual = build()
    if actual != expected:
        raise AssertionError("the frozen silver infinite-stabilizer certificate changed")
    if not actual["verified"]:
        raise AssertionError("the silver infinite-stabilizer checks do not all pass")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = verify_frozen() if args.check else build()
    print(json.dumps(document, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
