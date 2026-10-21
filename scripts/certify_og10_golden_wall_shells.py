#!/usr/bin/env python3
"""Exact certificate for the complete OG10 wall shells on the golden H10 slice.

The projective OG10 carrier has Neron--Severi lattice ``N`` of rank 22,
primitively embedded in the Beauville lattice ``Lambda``.  This certificate
proves two finite statements which make its wall search tractable:

1. ``N`` contains no primitive class of ambient divisibility three, so the
   OG10 wall types ``(-6,3)`` and ``(-24,3)`` do not occur; and
2. after restriction to the golden eigenspace, every remaining ``(-2,1)``
   or ``(-4,1)`` wall has field norm in one of three unit-square orbits,
   represented by ``2``, ``4``, and ``6-2*phi``.

The Markman/Mongardi--Onorati chamber-realization theorem is an external
geometric input and is stated in the accompanying research note.  This file
certifies the exact lattice and number-field calculations used there.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import isqrt
from typing import Any, Iterable

from sympy import Matrix, ZZ, eye
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5
from scripts import certify_og10_golden_h10_carrier as carrier


def _matrix_list(matrix: Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rref_mod_prime(
    matrix: Matrix, prime: int
) -> tuple[list[list[int]], list[int]]:
    """Return reduced row-echelon form and pivot columns over F_prime."""

    rows = [
        [int(matrix[row, column]) % prime for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(matrix.cols):
        source = next(
            (
                row
                for row in range(pivot_row, matrix.rows)
                if rows[row][column] % prime
            ),
            None,
        )
        if source is None:
            continue
        rows[pivot_row], rows[source] = rows[source], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [
            (entry * inverse) % prime for entry in rows[pivot_row]
        ]
        for row in range(matrix.rows):
            if row == pivot_row or not rows[row][column]:
                continue
            multiplier = rows[row][column]
            rows[row] = [
                (rows[row][index] - multiplier * rows[pivot_row][index])
                % prime
                for index in range(matrix.cols)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == matrix.rows:
            break
    return rows, pivot_columns


def _nullspace_mod_prime(matrix: Matrix, prime: int) -> list[list[int]]:
    reduced, pivots = _rref_mod_prime(matrix, prime)
    free_columns = [
        column for column in range(matrix.cols) if column not in pivots
    ]
    result: list[list[int]] = []
    for free in free_columns:
        vector = [0 for _ in range(matrix.cols)]
        vector[free] = 1
        for row, pivot in reversed(list(enumerate(pivots))):
            vector[pivot] = (
                -sum(
                    reduced[row][column] * vector[column]
                    for column in free_columns
                )
            ) % prime
        result.append(vector)
    return result


def _smith_nonunits(matrix: Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if abs(int(smith[index, index])) > 1
    ]


def _ambient_data() -> tuple[Matrix, Matrix]:
    """Rebuild Lambda and the primitive coordinate embedding N -> Lambda."""

    direct_sum = Matrix.diag(carrier.small_n_gram(), carrier.complement_gram())
    basis = carrier.overlattice_basis()
    glued = carrier._integer_matrix(basis.T * direct_sum * basis)
    ambient = Matrix.diag(carrier.core_gram(), glued)

    n_small_embedding = carrier._integer_matrix(
        basis.inv() * Matrix.vstack(eye(4), Matrix.zeros(2, 4))
    )
    n_embedding = Matrix.zeros(24, 22)
    n_embedding[:18, :18] = eye(18)
    n_embedding[18:, 18:] = n_small_embedding
    return ambient, n_embedding


def _golden_text(value: Golden) -> str:
    if value.b == 0:
        return str(value.a)
    sign = "+" if value.b > 0 else "-"
    if value.a == 0:
        return f"{value.b}*phi"
    return f"{value.a}{sign}{abs(value.b)}*phi"


def _shell(square_magnitude: int) -> list[Golden]:
    """Enumerate all totally positive even q=a+b*phi with a+b=s."""

    # If a+b=s, total positivity says -s*phi < b < s/phi.  Since
    # 1 < phi < 2, the deliberately wider interval [-2s,2s] is exhaustive.
    result: list[Golden] = []
    for b in range(-2 * square_magnitude, 2 * square_magnitude + 1, 2):
        a = square_magnitude - b
        q = Golden(a, b)
        if q.sign(0) > 0 and q.sign(1) > 0:
            result.append(q)
    return result


def _unit_square_exponent(left: Golden, right: Golden) -> int | None:
    """Return k with right=left*phi^(2k) for the two certified wall shells.

    At the identity embedding the seven shell values lie between
    ``2*phi^-4`` and ``2*phi^2``.  Thus a unit-square ratio between two of
    them has exponent difference at most three.
    """

    for exponent in range(-3, 4):
        if right == left * PHI ** (2 * exponent):
            return exponent
    return None


def _unit_square_orbits(values: Iterable[Golden]) -> list[list[Golden]]:
    remaining = set(values)
    result: list[list[Golden]] = []
    while remaining:
        seed = min(remaining, key=lambda value: (value.norm(), value.a, value.b))
        orbit = sorted(
            (
                value
                for value in remaining
                if _unit_square_exponent(seed, value) is not None
            ),
            key=lambda value: (value.a + value.b, value.b),
        )
        result.append(orbit)
        remaining.difference_update(orbit)
    return result


def _fraction_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Fraction(numerator, denominator)


def _is_square_in_field(value: Golden) -> bool:
    """Test exactly whether a totally real field element is a square in K."""

    norm_root = _fraction_square_root(value.norm())
    if norm_root is None:
        return False
    # Write value=A+B*sqrt(5).  For z=x+y*sqrt(5), z^2=value and
    # epsilon=x^2-5y^2=+/-sqrt(N(value)).
    rational_part = value.a + value.b / 2
    radical_part = value.b / 2
    for epsilon in {norm_root, -norm_root}:
        x_squared = (rational_part + epsilon) / 2
        y_squared = (rational_part - epsilon) / 10
        x = _fraction_square_root(x_squared)
        y = _fraction_square_root(y_squared)
        if x is None or y is None:
            continue
        if radical_part == 0 or 2 * x * y == abs(radical_part):
            return True
    return False


def _field_quadratic(vector: list[Golden | int]) -> Golden:
    form = carrier.field_gram()
    return sum(
        Golden.coerce(vector[row])
        * form[row][column]
        * Golden.coerce(vector[column])
        for row in range(len(vector))
        for column in range(len(vector))
    )


def certify() -> dict[str, Any]:
    carrier_certificate = carrier.certify()
    ambient, n_embedding = _ambient_data()
    coordinate_embedding = Matrix.vstack(eye(22), Matrix.zeros(2, 22))

    radical = _nullspace_mod_prime(ambient, 3)
    pairings_with_n = ambient * n_embedding
    _, pairing_pivots = _rref_mod_prime(pairings_with_n, 3)
    radical_intersection_n = _nullspace_mod_prime(pairings_with_n, 3)

    expected_radical = [0] * 19 + [2, 1, 0, 0, 1]
    radical_generator = Matrix(radical[0])
    radical_pairing = ambient * radical_generator
    golden_polynomial_mod_3 = [
        (value * value - value - 1) % 3 for value in range(3)
    ]

    alpha = -PHI * INV_SQRT5
    trace_coefficients = [(alpha * Golden(1)).trace(), (alpha * PHI).trace()]
    inverse_phi = PHI - 1
    inverse_root_trace_coefficients = [
        (alpha * inverse_phi * Golden(1)).trace(),
        (alpha * inverse_phi * PHI).trace(),
    ]
    d_square = (alpha * Golden(2)).trace()
    tau_d_square = (alpha * Golden(2) * PHI * PHI).trace()
    d_tau_d = (alpha * Golden(2) * PHI).trace()
    e_square = tau_d_square + d_square - 2 * d_tau_d
    d_e = d_tau_d - d_square
    divisorial_pair_gram = Matrix([[d_square, d_e], [d_e, e_square]])
    reflection_d = Matrix.diag(-1, 1)
    reflection_e = Matrix.diag(1, -1)
    paired_reflection = reflection_d * reflection_e
    tau_on_pair = Matrix([[1, 1], [1, 0]])
    field = carrier.field_gram()
    field_integral = all(entry.is_integral for row in field for entry in row)
    field_even_diagonal = all(
        (field[index][index] / 2).is_integral for index in range(len(field))
    )

    shell_2 = _shell(2)
    shell_4 = _shell(4)
    all_shell_values = shell_2 + shell_4
    orbits = _unit_square_orbits(all_shell_values)

    orbit_records: list[dict[str, Any]] = []
    preferred_representatives = [Golden(2), Golden(4), Golden(6, -2)]
    for preferred in preferred_representatives:
        orbit = next(
            group
            for group in orbits
            if any(_unit_square_exponent(preferred, value) is not None for value in group)
        )
        orbit = sorted(
            orbit, key=lambda value: _unit_square_exponent(preferred, value)
        )
        orbit_records.append(
            {
                "representative": preferred.to_json(),
                "representative_text": _golden_text(preferred),
                "field_norm": int(preferred.norm()),
                "wall_shell_values": [value.to_json() for value in orbit],
                "wall_shell_text": [_golden_text(value) for value in orbit],
                "unit_square_exponents_from_representative": [
                    _unit_square_exponent(preferred, value) for value in orbit
                ],
            }
        )

    representative_vectors = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
    representative_vector_norms = [
        _field_quadratic(vector) for vector in representative_vectors
    ]

    checks = {
        "inherited_og10_carrier_certificate_verified": carrier_certificate[
            "verified"
        ],
        "n_is_first_22_ambient_coordinates": n_embedding == coordinate_embedding,
        "ambient_determinant_minus_3": int(ambient.det()) == -3,
        "ambient_smith_nonunit_3": _smith_nonunits(ambient) == [3],
        "ambient_mod3_radical_dimension_1": len(radical) == 1,
        "ambient_mod3_radical_generator": radical[0] == expected_radical,
        "radical_generator_pairs_in_multiples_of_3": all(
            int(entry) % 3 == 0 for entry in radical_pairing
        ),
        "n_pairing_map_has_full_rank_mod3": len(pairing_pivots) == 22,
        "ambient_radical_intersects_n_mod3_trivially": not radical_intersection_n,
        "golden_polynomial_has_no_root_mod3": all(golden_polynomial_mod_3),
        "trace_formula_B_of_q_is_minus_a_plus_b": trace_coefficients == [-1, -1],
        "two_trace_pairings_detect_field_orthogonality": (
            inverse_root_trace_coefficients == [0, -1]
        ),
        "field_lattice_integral_and_even": field_integral and field_even_diagonal,
        "square_minus_2_shell_complete": shell_2 == [Golden(4, -2), Golden(2)],
        "square_minus_4_shell_complete": shell_4
        == [
            Golden(10, -6),
            Golden(8, -4),
            Golden(6, -2),
            Golden(4),
            Golden(2, 2),
        ],
        "seven_raw_shell_values": len(all_shell_values) == 7,
        "three_unit_square_mirror_orbits": len(orbits) == 3,
        "preferred_orbit_representatives_cover_shells": {
            value for record in orbit_records for value in record["wall_shell_text"]
        }
        == {_golden_text(value) for value in all_shell_values},
        "three_representative_families_are_nonempty": (
            representative_vector_norms == preferred_representatives
        ),
        "three_mirror_families_are_pairwise_distinct": (
            not _is_square_in_field(Golden(4) / Golden(2))
            and not _is_square_in_field(Golden(6, -2) / Golden(2))
            and not _is_square_in_field(Golden(6, -2) / Golden(4))
        ),
        "extra_kahler_mirrors_are_not_integral_field_reflections": (
            not (Golden(2) / Golden(4)).is_integral
            and not (Golden(2) / Golden(6, -2)).is_integral
        ),
        "norm_2_root_has_orthogonal_minus_2_pair": (
            divisorial_pair_gram == Matrix.diag(-2, -2)
        ),
        "paired_divisorial_reflections_are_minus_identity_on_O_root_line": (
            paired_reflection == -eye(2)
        ),
        "paired_reflection_centralizes_golden_action": (
            paired_reflection * tau_on_pair == tau_on_pair * paired_reflection
        ),
    }

    return {
        "candidate_id": "og10_golden_complete_wall_shells",
        "inherited_carrier_candidate_id": carrier_certificate["candidate_id"],
        "ambient_mod_3": {
            "ambient_rank": 24,
            "ambient_determinant": int(ambient.det()),
            "ambient_smith_nonunits": _smith_nonunits(ambient),
            "radical_dimension": len(radical),
            "radical_generator": radical[0],
            "integral_pairing_of_generator": [int(entry) for entry in radical_pairing],
            "n_coordinate_embedding": "first_22_coordinates",
            "pairing_map_from_n_rank_mod_3": len(pairing_pivots),
            "radical_intersection_n_dimension": len(radical_intersection_n),
            "golden_polynomial_values_mod_3": golden_polynomial_mod_3,
            "conclusion": (
                "ambient_divisibility_3_for_D_in_N_implies_D_is_3_divisible_in_N"
            ),
            "self_adjoint_extension_obstruction": (
                "no_integral_self_adjoint_golden_endomorphism_extends_tau_to_"
                "the_full_OG10_Beauville_lattice"
            ),
        },
        "og10_wall_types_on_this_n": {
            "surviving": [
                {"square": -2, "ambient_divisibility": 1},
                {"square": -4, "ambient_divisibility": 1},
            ],
            "excluded": [
                {"square": -6, "ambient_divisibility": 3},
                {"square": -24, "ambient_divisibility": 3},
            ],
            "complete_numerical_wall_arrangement": (
                "all_square_minus_2_and_square_minus_4_hyperplanes_in_N"
            ),
        },
        "field_trace_formula": {
            "h_value": "a+b*phi",
            "beauville_square": "-(a+b)",
            "h_is_even_over_Z_phi": True,
            "mirror_meets_golden_positive_cone_iff": (
                "h_value_is_positive_at_the_identity_embedding"
            ),
            "conjugate_embedding_is_automatically_positive_for_nonzero_vectors": True,
        },
        "raw_totally_positive_even_shells": {
            "square_minus_2": [value.to_json() for value in shell_2],
            "square_minus_2_text": [_golden_text(value) for value in shell_2],
            "square_minus_4": [value.to_json() for value in shell_4],
            "square_minus_4_text": [_golden_text(value) for value in shell_4],
        },
        "restricted_mirror_unit_square_orbits": orbit_records,
        "complete_restricted_search_norms": [
            preferred.to_json() for preferred in preferred_representatives
        ],
        "complete_restricted_search_norms_text": [
            _golden_text(preferred) for preferred in preferred_representatives
        ],
        "explicit_representative_vectors": [
            {
                "coordinates_in_field_basis": vector,
                "h_norm": norm.to_json(),
                "h_norm_text": _golden_text(norm),
            }
            for vector, norm in zip(
                representative_vectors, representative_vector_norms, strict=True
            )
        ],
        "pairwise_distinct_mirror_family_square_ratios": {
            "norm_4_ratio_to_norm_2": "2_is_not_a_square_in_Q(sqrt(5))",
            "norm_6_minus_2phi_ratio_to_norm_2": (
                "3-phi_has_field_norm_5_so_is_not_a_square"
            ),
            "norm_6_minus_2phi_ratio_to_norm_4": (
                "(3-phi)/2_has_field_norm_5/4_so_is_not_a_square"
            ),
        },
        "extra_kahler_families_are_nonreflective_over_Z_phi": {
            "test_pairing_with_adjacent_E8_root": -1,
            "norm_4_reflection_coefficient": (
                Golden(2) / Golden(4)
            ).to_json(),
            "norm_6_minus_2phi_reflection_coefficient": (
                Golden(2) / Golden(6, -2)
            ).to_json(),
        },
        "divisorial_factorization_of_norm_2_reflection": {
            "classes": ["D=r", "E=tau_inverse_D=(tau-1)D"],
            "pair_gram": _matrix_list(divisorial_pair_gram),
            "tau_on_pair_basis_D_E": _matrix_list(tau_on_pair),
            "reflection_in_D_on_pair": _matrix_list(reflection_d),
            "reflection_in_E_on_pair": _matrix_list(reflection_e),
            "product_on_O_root_line": _matrix_list(paired_reflection),
            "product_on_H_orthogonal": "identity",
            "conclusion": (
                "golden_field_reflection_equals_product_of_two_commuting_"
                "stably_prime_exceptional_OG10_reflections"
            ),
            "group_membership": (
                "restriction_to_N_of_W_Exc_intersection_Cent_O(N)(tau)"
            ),
        },
        "geometric_external_input": {
            "statement": (
                "every_component_of_the_full_numerical_wall_arrangement_is_the_"
                "marked_Kahler_cone_of_a_same_period_projective_OG10_model"
            ),
            "sources": [
                "Markman_Lemma_5.11_and_Theorem_5.16_arXiv_1101.4606",
                "Mongardi_Onorati_Proposition_5.4_and_Theorem_5.5_arXiv_2010.12511",
            ],
            "requires_hodge_monodromy_remarking_for_arbitrary_chambers": True,
        },
        "checks": checks,
        "verified": all(checks.values()),
    }


def main() -> None:
    print(json.dumps(certify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
