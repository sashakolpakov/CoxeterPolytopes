#!/usr/bin/env python3
"""Audit the finite arithmetic front end of the maximal OG10 H10 compiler.

For a rank-22 Picard lattice N primitively embedded in

    Lambda_OG10 = 3U + 2E8(-1) + A2(-1),  |det Lambda_OG10| = 3,

the orthogonal complement is a positive even binary lattice T.  If k is the
order of the discriminant gluing graph, then

    det(T) = 3*k^2/|det(N)|.

Both projections of the graph are injective.  Hence k divides both
|A_N| and |A_T|.  These elementary restrictions leave finitely many k and
finitely many Gauss-reduced positive even binary forms.  This script freezes
that first, deliberately necessary-only, enumeration for the four exact
rank-22 lattices already present in the repository.

The accompanying theorem explains the remaining finite steps: discriminant
anti-isometries, OG10-lattice identification, the two oriented rank-two
periods, and a compact-chamber ellipsoid enumeration of every possible OG10
wall divisor meeting the input polytope.
"""

from __future__ import annotations

import hashlib
import json
from math import isqrt
from pathlib import Path
from typing import Any

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import PHI
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    trace_gram,
)
from scripts.certify_og10_golden_a2_coupled_candidate import (
    field_gram as a2_coupled_field_gram,
)
from scripts.certify_og10_golden_h10_carrier import expected_trace_gram
from scripts.certify_og10_golden_norm_six_embedding import field_and_trace_grams


ROOT = Path(__file__).resolve().parents[1]
D4_CERTIFICATE_PATH = (
    ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith_diagonal(gram: Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i]
    ]


def _prime_factors(value: int) -> list[int]:
    value = abs(value)
    result: list[int] = []
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            result.append(prime)
            while value % prime == 0:
                value //= prime
        prime += 1
    if value > 1:
        result.append(value)
    return result


def _p_length(smith: list[int], prime: int) -> int:
    return sum(entry % prime == 0 for entry in smith)


def _divisors(value: int) -> list[int]:
    return [candidate for candidate in range(1, value + 1) if value % candidate == 0]


def reduced_even_positive_binary_forms(determinant: int) -> list[list[int]]:
    """Return one GL(2,Z) representative [[2a,b],[b,2c]] per reduced form.

    We use 0 <= b <= a <= c.  Changing the sign of one basis vector identifies
    b and -b for integral-lattice (rather than proper-class) equivalence.
    """

    forms: list[list[int]] = []
    # Reduction gives determinant = 4ac-b^2 >= 3a^2.
    for a in range(1, isqrt(determinant // 3) + 2):
        for b in range(0, a + 1):
            numerator = determinant + b * b
            denominator = 4 * a
            if numerator % denominator:
                continue
            c = numerator // denominator
            if a > c:
                continue
            forms.append([a, b, c])
    return forms


def _binary_gram(abc: list[int]) -> Matrix:
    a, b, c = abc
    return Matrix([[2 * a, b], [b, 2 * c]])


def _vectors_of_norm(gram: Matrix, target: int) -> list[tuple[int, int]]:
    """Enumerate exactly all integral vectors of a fixed positive norm."""

    two_a = int(gram[0, 0])
    b = int(gram[0, 1])
    determinant = int(gram.det())
    # q(x,y)=((two_a*x+b*y)^2+determinant*y^2)/two_a.
    y_bound = isqrt((two_a * target) // determinant) + 1
    linear_bound = isqrt(two_a * target)
    result: list[tuple[int, int]] = []
    for y in range(-y_bound, y_bound + 1):
        x_bound = (linear_bound + abs(b * y)) // two_a + 1
        for x in range(-x_bound, x_bound + 1):
            vector = Matrix([x, y])
            if int((vector.T * gram * vector)[0]) == target:
                result.append((x, y))
    return result


def binary_forms_are_gl2_isometric(left: list[int], right: list[int]) -> bool:
    """Decide integral isometry of two positive binary forms exactly."""

    gram_left = _binary_gram(left)
    gram_right = _binary_gram(right)
    if gram_left.det() != gram_right.det():
        return False
    first_vectors = _vectors_of_norm(gram_left, int(gram_right[0, 0]))
    second_vectors = _vectors_of_norm(gram_left, int(gram_right[1, 1]))
    target_pairing = int(gram_right[0, 1])
    for first in first_vectors:
        first_column = Matrix(first)
        for second in second_vectors:
            second_column = Matrix(second)
            if int((first_column.T * gram_left * second_column)[0]) != target_pairing:
                continue
            change = Matrix.hstack(first_column, second_column)
            if abs(int(change.det())) == 1:
                return True
    return False


def binary_forms_are_pairwise_gl2_distinct(forms: list[list[int]]) -> bool:
    """Audit that a reduced-form list has no integral-isometry duplicate."""

    return all(
        not binary_forms_are_gl2_isometric(forms[left], forms[right])
        for left in range(len(forms))
        for right in range(left + 1, len(forms))
    )


def determinant_lanes(determinant_n: int) -> list[dict[str, Any]]:
    """Necessary determinant/index lanes for primitive N -> Lambda_OG10."""

    lanes: list[dict[str, Any]] = []
    for glue_order in _divisors(determinant_n):
        numerator = 3 * glue_order * glue_order
        if numerator % determinant_n:
            continue
        determinant_t = numerator // determinant_n
        # The graph projection to A_T is injective.
        if determinant_t % glue_order:
            continue
        # An even binary Gram [[2a,b],[b,2c]] has determinant 0 or 3 mod 4.
        if determinant_t % 4 not in (0, 3):
            continue
        forms = reduced_even_positive_binary_forms(determinant_t)
        if not forms:
            continue
        lanes.append(
            {
                "glue_order": glue_order,
                "complement_determinant": determinant_t,
                "reduced_even_binary_forms_abc": forms,
                "reduced_form_count": len(forms),
            }
        )
    return lanes


def _d4_gram() -> Matrix:
    payload = json.loads(D4_CERTIFICATE_PATH.read_text())
    return Matrix(
        [
            [int(entry) for entry in row]
            for row in payload["coxeter_input"]["basis_gram"]
        ]
    )


def _known_grams() -> dict[str, Matrix]:
    golden = Matrix(expected_trace_gram())
    a2_coupled = Matrix(
        integral_gram(
            trace_gram(a2_coupled_field_gram(), -PHI * INV_SQRT5)
        )
    )
    _, norm_six = field_and_trace_grams()
    return {
        "golden_det44": golden,
        "a2_coupled_det20": a2_coupled,
        "norm_six_det36": norm_six,
        "borcherds_d4_det4": _d4_gram(),
    }


def _case_payload(name: str, gram: Matrix) -> dict[str, Any]:
    determinant_signed = int(gram.det())
    determinant = abs(determinant_signed)
    smith = _smith_diagonal(gram)
    primes = sorted(set(_prime_factors(determinant)) | {3})
    p_lengths = {str(prime): _p_length(smith, prime) for prime in primes}
    length_gate = all(
        length <= (3 if prime == 3 else 2)
        for prime, length in ((p, _p_length(smith, p)) for p in primes)
    )
    return {
        "id": name,
        "rank": gram.rows,
        "determinant_signed": determinant_signed,
        "determinant_absolute": determinant,
        "smith_nonunits": [entry for entry in smith if entry != 1],
        "p_primary_lengths": p_lengths,
        "passes_og10_corank_two_length_gate": length_gate,
        "determinant_lanes": determinant_lanes(determinant),
        "three_coprime": determinant % 3 != 0,
        "divisibility_three_picard_walls_automatically_absent": determinant % 3
        != 0,
    }


def build_certificate() -> dict[str, Any]:
    cases = {
        name: _case_payload(name, gram) for name, gram in _known_grams().items()
    }
    expected_lanes = {
        "golden_det44": [(44, 132, 4)],
        "a2_coupled_det20": [(20, 60, 4)],
        "norm_six_det36": [(12, 12, 2), (36, 108, 5)],
        "borcherds_d4_det4": [(4, 12, 2)],
    }
    observed_lanes = {
        name: [
            (
                lane["glue_order"],
                lane["complement_determinant"],
                lane["reduced_form_count"],
            )
            for lane in case["determinant_lanes"]
        ]
        for name, case in cases.items()
    }
    checks = {
        "all_four_inputs_have_rank_22": all(case["rank"] == 22 for case in cases.values()),
        "all_four_pass_corank_two_length_gate": all(
            case["passes_og10_corank_two_length_gate"] for case in cases.values()
        ),
        "determinant_lane_census_matches": observed_lanes == expected_lanes,
        "all_reduced_binary_lists_are_pairwise_GL2_distinct": all(
            binary_forms_are_pairwise_gl2_distinct(
                lane["reduced_even_binary_forms_abc"]
            )
            for case in cases.values()
            for lane in case["determinant_lanes"]
        ),
        "three_coprime_cases_have_one_forced_full_glue_lane": all(
            len(case["determinant_lanes"]) == 1
            and case["determinant_lanes"][0]["glue_order"]
            == case["determinant_absolute"]
            for name, case in cases.items()
            if name != "norm_six_det36"
        ),
        "norm_six_has_partial_and_full_three_primary_lanes": observed_lanes[
            "norm_six_det36"
        ]
        == expected_lanes["norm_six_det36"],
        "coherent_full_h2_quadratic_bound_is_H9_on_OG10": 2 * (9 + 3) == 24
        and 2 * (10 + 3) > 24,
        "coherent_full_h2_quadratic_bound_is_H8_on_K3": 2 * (8 + 3) == 22
        and 2 * (9 + 3) > 22,
    }
    if not all(checks.values()):
        raise AssertionError(f"maximal OG10 compiler audit failed: {checks}")

    payload: dict[str, Any] = {
        "schema_version": 1,
        "certificate_kind": "maximal_quadratic_H10_OG10_finite_compiler_front_end",
        "ambient": {
            "lattice": "3U + 2E8(-1) + A2(-1)",
            "rank": 24,
            "signature": [3, 21],
            "determinant_absolute": 3,
            "mod_p_radical_dimension": {"p_not_3": 0, "p_3": 1},
        },
        "universal_filters": {
            "p_length": "ell_p(A_N)<=2 for p!=3 and ell_3(A_N)<=3",
            "determinant_formula": "det(T)=3*k^2/abs(det(N))",
            "glue_injectivity": "k divides abs(det(N)) and k divides det(T)",
            "binary_form_reduction": (
                "Gram [[2a,b],[b,2c]] with 0<=b<=a<=c; b and -b are "
                "identified by an orientation-reversing basis change, and "
                "the frozen lists are pairwise audited for GL(2,Z) isometry"
            ),
            "three_coprime_consequence": (
                "if 3 does not divide det(N), then k=abs(det(N)), "
                "det(T)=3*abs(det(N)), and no primitive Picard class has "
                "ambient divisibility 3"
            ),
            "wall_types": [
                [-2, 1],
                [-4, 1],
                [-6, 3],
                [-24, 3],
            ],
            "compact_hidden_wall_bound": (
                "Qhat_p(delta)=2*(p,delta)^2/(p,p)-(delta,delta) "
                "<= (-(delta,delta))*(2*cosh(R)^2-1)"
            ),
            "coherent_full_h2_bound": (
                "a unital self-adjoint degree-d totally-real Hodge action "
                "requires d*(n+1) divisor slots and at least 2*d "
                "transcendental slots, hence d*(n+3)<=b2"
            ),
        },
        "known_rank22_replays": cases,
        "known_selected_complements": {
            "golden_det44": {"abc": [1, 0, 33], "gram": [[2, 0], [0, 66]]},
            "a2_coupled_det20": {"abc": [1, 0, 15], "gram": [[2, 0], [0, 30]]},
            "norm_six_det36": {"abc": [1, 0, 3], "gram": [[2, 0], [0, 6]]},
            "borcherds_d4_det4": {"abc": [2, 2, 2], "gram": [[4, 2], [2, 4]]},
        },
        "scope": {
            "proved": (
                "The arithmetic complement list is finite for every fixed "
                "rank-22 trace lattice; the accompanying theorem makes the "
                "complete OG10 hidden-wall test finite for a compact input P."
            ),
            "not_proved": (
                "No universal carrier-success theorem for all compact arithmetic "
                "H10 polytopes, no global exhaustion of all trace lattices, and "
                "no n<=9 conclusion for NS-only carriers.  The n<=9 conclusion "
                "is proved only for the coherent full-H2 real-multiplication branch."
            ),
        },
        "checks": checks,
        "source_sha256": {
            "golden_carrier_certificate": _sha256(
                ROOT / "data/candidates/og10_golden_h10_carrier.json"
            ),
            "a2_coupled_certificate": _sha256(
                ROOT / "data/candidates/og10_golden_a2_coupled_candidate.json"
            ),
            "norm_six_embedding_certificate": _sha256(
                ROOT / "data/candidates/og10_golden_norm_six_embedding.json"
            ),
            "borcherds_d4_certificate": _sha256(D4_CERTIFICATE_PATH),
        },
    }
    payload["certificate_payload_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    payload = build_certificate()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
