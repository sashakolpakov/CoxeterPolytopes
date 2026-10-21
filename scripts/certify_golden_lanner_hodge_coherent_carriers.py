"""Certify rational full-H2 golden actions for the complete Lanner returns.

The input chambers and their exact same-polytope K3 eigen-nef returns are
frozen by ``certify_golden_lanner_complete_chambers.py``.  Here we keep the
Picard trace lattice and extend multiplication by ``phi`` rationally across
the orthogonal complement in the K3 lattice.

For a field Gram ``H`` of rank ``r`` and determinant ``d``, put ``m=11-r``
and use the diagonal complementary field form

    C(t) = <-phi*t, d/t, 1, ..., 1>,

where ``t`` is totally positive.  Its trace has signature ``(2,2m-2)`` and
the correct determinant square class.  Hasse--Minkowski invariants decide
whether ``Tr(alpha C(t))`` is the rational orthogonal complement of the
Picard trace lattice inside the K3 lattice.  This is an existence
certificate for a rational full-H2 action.  Clearing denominators gives an
integral action of some finite-index order, but not a certified action of the
maximal order.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5, trace_gram
from coxeter_polytopes.lattices import k3_lattice_gram


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/golden_lanner_exhaustive_census.json"
ALPHA = -PHI * INV_SQRT5

# These twelve twists are exact survivors.  The nontrivial entries were found
# by the bounded ``--search-bound 5`` replay and are now frozen here.
TWISTS: dict[str, Golden | None] = {
    "triangle_255": Golden(1),
    "triangle_335": Golden(2),
    "triangle_355": Golden(2, 1),
    "triangle_555": Golden(2, 1),
    "chain_353": Golden(1),
    "chain_535": Golden(1),
    "star_5_3_3": Golden(1),
    "cycle4_one5": Golden(2, 1),
    "cycle4_opposite5": Golden(3),
    "chain_5333": Golden(1),
    "chain_5335": Golden(1),
    "fork5_terminal5": Golden(1),
}


def _diagonalize(gram: Sequence[Sequence[int | Fraction]]) -> list[Fraction]:
    """Diagonalize a nonsingular rational symmetric form by congruence."""

    current = [[Fraction(entry) for entry in row] for row in gram]
    diagonal: list[Fraction] = []
    while current:
        size = len(current)
        pivot = next((i for i in range(size) if current[i][i]), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(size)
                    for j in range(i + 1, size)
                    if current[i][j]
                ),
                None,
            )
            if pair is None:
                raise ArithmeticError("the rational form is singular")
            i, j = pair
            columns: list[list[Fraction]] = []
            first = [Fraction(0) for _ in range(size)]
            first[i] = first[j] = Fraction(1)
            columns.append(first)
            for index in range(size):
                if index == i:
                    continue
                vector = [Fraction(0) for _ in range(size)]
                vector[index] = Fraction(1)
                columns.append(vector)
            current = [
                [
                    sum(
                        (
                            columns[a][u]
                            * current[u][v]
                            * columns[b][v]
                            for u in range(size)
                            for v in range(size)
                        ),
                        start=Fraction(0),
                    )
                    for b in range(size)
                ]
                for a in range(size)
            ]
            pivot = 0
        if pivot:
            order = [pivot] + [i for i in range(size) if i != pivot]
            current = [[current[i][j] for j in order] for i in order]
        value = current[0][0]
        diagonal.append(value)
        current = [
            [
                current[i][j] - current[i][0] * current[0][j] / value
                for j in range(1, size)
            ]
            for i in range(1, size)
        ]
    return diagonal


def _valuation_and_unit(value: Fraction, prime: int) -> tuple[int, Fraction]:
    numerator = value.numerator
    denominator = value.denominator
    valuation = 0
    while numerator % prime == 0:
        numerator //= prime
        valuation += 1
    while denominator % prime == 0:
        denominator //= prime
        valuation -= 1
    return valuation, Fraction(numerator, denominator)


def _unit_mod(unit: Fraction, modulus: int) -> int:
    return unit.numerator * pow(unit.denominator, -1, modulus) % modulus


def _hilbert(left: Fraction, right: Fraction, prime: int | str) -> int:
    """Return the rational Hilbert symbol, including ``prime='infinity'``."""

    if prime == "infinity":
        return -1 if left < 0 and right < 0 else 1
    if not isinstance(prime, int):
        raise TypeError("a finite prime must be an integer")
    a, unit_left = _valuation_and_unit(left, prime)
    b, unit_right = _valuation_and_unit(right, prime)
    if prime == 2:
        u = _unit_mod(unit_left, 8)
        v = _unit_mod(unit_right, 8)
        exponent = (
            ((u - 1) // 2) * ((v - 1) // 2)
            + a * ((v * v - 1) // 8)
            + b * ((u * u - 1) // 8)
        ) % 2
        return -1 if exponent else 1

    sign = -1 if (a * b * ((prime - 1) // 2)) % 2 else 1

    def legendre(unit: Fraction) -> int:
        residue = pow(_unit_mod(unit, prime), (prime - 1) // 2, prime)
        return -1 if residue == prime - 1 else 1

    if b % 2:
        sign *= legendre(unit_left)
    if a % 2:
        sign *= legendre(unit_right)
    return sign


def _prime_factors(integer: int) -> set[int]:
    value = abs(integer)
    factors: set[int] = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.add(value)
    return factors


def rational_invariants(
    gram: Sequence[Sequence[int | Fraction]],
) -> dict[str, object]:
    diagonal = _diagonalize(gram)
    determinant = math.prod(diagonal, start=Fraction(1))
    primes = {2}
    for entry in diagonal:
        primes.update(_prime_factors(entry.numerator))
        primes.update(_prime_factors(entry.denominator))
    places: list[int | str] = sorted(primes) + ["infinity"]
    hasse = {
        str(place): math.prod(
            (
                _hilbert(diagonal[i], diagonal[j], place)
                for i in range(len(diagonal))
                for j in range(i + 1, len(diagonal))
            ),
            start=1,
        )
        for place in places
    }
    # At every omitted odd prime all diagonal coefficients are units, so all
    # pairwise Hilbert symbols are +1.  Hilbert reciprocity is therefore an
    # exact internal check that the finite place list is complete.
    reciprocity_product = math.prod(hasse.values(), start=1)
    if reciprocity_product != 1:
        raise ArithmeticError("Hilbert reciprocity failed")
    return {
        "rank": len(diagonal),
        "signature": [
            sum(entry > 0 for entry in diagonal),
            sum(entry < 0 for entry in diagonal),
        ],
        "determinant": str(determinant),
        "diagonal": [str(entry) for entry in diagonal],
        "hasse": hasse,
        "all_nontrivial_places_enumerated": True,
        "hilbert_reciprocity_product": reciprocity_product,
    }


def _is_square(value: Fraction) -> bool:
    return (
        value >= 0
        and math.isqrt(value.numerator) ** 2 == value.numerator
        and math.isqrt(value.denominator) ** 2 == value.denominator
    )


def _rational_trace_gram(field_gram: Sequence[Sequence[Golden]]) -> list[list[Fraction]]:
    result: list[list[Fraction]] = []
    for row in trace_gram(field_gram, ALPHA):
        rational_row: list[Fraction] = []
        for entry in row:
            if isinstance(entry, Golden):
                if not entry.is_rational:
                    raise ArithmeticError("trace Gram retained an irrational entry")
                rational_row.append(entry.a)
            else:
                rational_row.append(Fraction(entry))
        result.append(rational_row)
    return result


def complement_field_gram(
    field_determinant: Golden, rank: int, twist: Golden
) -> list[list[Golden]]:
    if twist.sign(0) <= 0 or twist.sign(1) <= 0:
        raise ArithmeticError("the complement twist is not totally positive")
    entries = [-PHI * twist, field_determinant / twist]
    entries.extend(Golden(1) for _ in range(rank - 2))
    gram = [[Golden(0) for _ in range(rank)] for _ in range(rank)]
    for i, entry in enumerate(entries):
        gram[i][i] = entry
    return gram


def _orthogonal_sum(
    left: Sequence[Sequence[int | Fraction]],
    right: Sequence[Sequence[int | Fraction]],
) -> list[list[Fraction]]:
    a = len(left)
    b = len(right)
    return [
        [Fraction(left[i][j]) if i < a and j < a else Fraction(0) for j in range(a + b)]
        for i in range(a)
    ] + [
        [Fraction(0) for _ in range(a)] + [Fraction(right[i][j]) for j in range(b)]
        for i in range(b)
    ]


def _golden_action(rank: int) -> list[list[int]]:
    """Multiplication by phi in the interleaved basis ``e_i, phi*e_i``."""

    action = [[0 for _ in range(2 * rank)] for _ in range(2 * rank)]
    for index in range(rank):
        e = 2 * index
        phi_e = e + 1
        action[phi_e][e] = 1
        action[e][phi_e] = 1
        action[phi_e][phi_e] = 1
    return action


def _matmul(
    left: Sequence[Sequence[int | Fraction]],
    right: Sequence[Sequence[int | Fraction]],
) -> list[list[Fraction]]:
    return [
        [
            sum(
                (Fraction(left[i][k]) * Fraction(right[k][j]) for k in range(len(right))),
                start=Fraction(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _transpose(
    matrix: Sequence[Sequence[int | Fraction]],
) -> list[list[int | Fraction]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def _identity(rank: int) -> list[list[int]]:
    return [[int(i == j) for j in range(rank)] for i in range(rank)]


def _matrix_add(
    left: Sequence[Sequence[int | Fraction]],
    right: Sequence[Sequence[int | Fraction]],
) -> list[list[Fraction]]:
    return [
        [Fraction(left[i][j]) + Fraction(right[i][j]) for j in range(len(left))]
        for i in range(len(left))
    ]


def _same_rational_form(left: dict[str, object], right: dict[str, object]) -> bool:
    if left["rank"] != right["rank"] or left["signature"] != right["signature"]:
        return False
    determinant_ratio = Fraction(str(left["determinant"])) / Fraction(
        str(right["determinant"])
    )
    if not _is_square(determinant_ratio):
        return False
    left_hasse = left["hasse"]
    right_hasse = right["hasse"]
    assert isinstance(left_hasse, dict) and isinstance(right_hasse, dict)
    return all(
        left_hasse.get(place, 1) == right_hasse.get(place, 1)
        for place in set(left_hasse) | set(right_hasse)
    )


def _search_twists(bound: int) -> Iterable[Golden]:
    seen: set[Golden] = set()
    for denominator in range(1, bound + 1):
        for a in range(-bound, bound + 1):
            for b in range(-bound, bound + 1):
                value = Golden(Fraction(a, denominator), Fraction(b, denominator))
                if value in seen or not value:
                    continue
                seen.add(value)
                if value.sign(0) > 0 and value.sign(1) > 0:
                    yield value


def analyze(search_bound: int = 0) -> dict[str, object]:
    source = json.loads(SOURCE.read_text())
    if source.get("all_checks_pass") is not True or len(source["examples"]) != 12:
        raise ArithmeticError("the upstream exhaustive same-P certificate is not frozen")
    if ALPHA.sign(0) >= 0 or ALPHA.sign(1) >= 0:
        raise ArithmeticError("the trace scalar must be negative at both real places")
    k3 = rational_invariants(k3_lattice_gram())
    rows: list[dict[str, object]] = []
    for example in source["examples"]:
        seed_id = example["seed_id"]
        package = example["trace_k3_package"]
        picard = package["gram_interleaved_basis_ei_phi_ei"]
        field_rank = int(example["dimension"]) + 1
        field_determinant = Golden.from_json(
            example.get("field_gram_determinant", example.get("field_determinant"))
        )
        if field_determinant.sign(0) >= 0 or field_determinant.sign(1) <= 0:
            raise ArithmeticError(f"wrong determinant signs for {seed_id}")
        complement_rank = 11 - field_rank
        if package["rank"] != 2 * field_rank:
            raise ArithmeticError(f"wrong Picard restriction rank for {seed_id}")
        if package["k3_complement_rank"] != 2 * complement_rank:
            raise ArithmeticError(f"wrong complement restriction rank for {seed_id}")

        picard_action = _golden_action(field_rank)
        picard_identity = _identity(2 * field_rank)
        picard_action_matches_source = (
            picard_action == package["golden_action_matrix"]
        )
        picard_action_minimal_polynomial = _matmul(
            picard_action, picard_action
        ) == _matrix_add(picard_action, picard_identity)
        picard_action_self_adjoint = _matmul(
            _transpose(picard_action), picard
        ) == _matmul(picard, picard_action)
        picard_invariants = rational_invariants(picard)
        if picard_invariants["signature"] != [1, 2 * field_rank - 1]:
            raise ArithmeticError(f"wrong Picard signature for {seed_id}")

        frozen_twist = TWISTS[seed_id]
        discovery_twist: Golden | None = None
        discovery_tested = 0
        if search_bound:
            for candidate in _search_twists(search_bound):
                discovery_tested += 1
                complement = _rational_trace_gram(
                    complement_field_gram(field_determinant, complement_rank, candidate)
                )
                total = rational_invariants(_orthogonal_sum(picard, complement))
                if _same_rational_form(total, k3):
                    discovery_twist = candidate
                    break
            if discovery_twist != frozen_twist:
                raise ArithmeticError(
                    f"bounded discovery did not recover frozen twist for {seed_id}"
                )
        twist = frozen_twist
        if twist is None:
            rows.append(
                {
                    "seed_id": seed_id,
                    "dimension": example["dimension"],
                    "same_polytope_input_certified": example["all_checks_pass"],
                    "picard_rank": package["rank"],
                    "complement_rank": package["k3_complement_rank"],
                    "status": "open_outside_the_certified_twist_family",
                    "bounded_twist_search": search_bound,
                }
            )
            continue
        complement_field = complement_field_gram(
            field_determinant, complement_rank, twist
        )
        complement = _rational_trace_gram(complement_field)
        action = _golden_action(complement_rank)
        identity = _identity(2 * complement_rank)
        action_square = _matmul(action, action)
        action_plus_identity = _matrix_add(action, identity)
        self_adjoint = _matmul(_transpose(action), complement) == _matmul(
            complement, action
        )
        sqrt5_action = [
            [2 * action[i][j] - identity[i][j] for j in range(len(action))]
            for i in range(len(action))
        ]
        five_identity = [
            [5 * identity[i][j] for j in range(len(identity))]
            for i in range(len(identity))
        ]
        sqrt5_action_square = _matmul(sqrt5_action, sqrt5_action)
        sqrt5_similarity = _matmul(
            _matmul(_transpose(sqrt5_action), complement), sqrt5_action
        ) == [
            [5 * Fraction(entry) for entry in row]
            for row in complement
        ]
        field_diagonal = [
            complement_field[i][i] for i in range(complement_rank)
        ]
        identity_eigenspace_signature = [
            sum(value.sign(0) < 0 for value in field_diagonal),
            sum(value.sign(0) > 0 for value in field_diagonal),
        ]
        conjugate_eigenspace_signature = [
            sum(value.sign(1) < 0 for value in field_diagonal),
            sum(value.sign(1) > 0 for value in field_diagonal),
        ]
        complement_invariants = rational_invariants(complement)
        expected_complement_signature = [2, 2 * complement_rank - 2]
        if complement_invariants["signature"] != expected_complement_signature:
            raise ArithmeticError(f"wrong complement trace signature for {seed_id}")

        total_gram = _orthogonal_sum(picard, complement)
        total_invariants = rational_invariants(total_gram)
        pass_equivalence = _same_rational_form(total_invariants, k3)
        full_action = _orthogonal_sum(picard_action, action)
        full_identity = _identity(22)
        full_action_minimal_polynomial = _matmul(
            full_action, full_action
        ) == _matrix_add(full_action, full_identity)
        full_action_self_adjoint = _matmul(
            _transpose(full_action), total_gram
        ) == _matmul(total_gram, full_action)
        all_action_gates = all(
            (
                picard_action_matches_source,
                picard_action_minimal_polynomial,
                picard_action_self_adjoint,
                action_square == action_plus_identity,
                self_adjoint,
                sqrt5_action_square == five_identity,
                sqrt5_similarity,
                full_action_minimal_polynomial,
                full_action_self_adjoint,
            )
        )
        period_plane_available = identity_eigenspace_signature == [
            2,
            complement_rank - 2,
        ] and conjugate_eigenspace_signature == [0, complement_rank]
        exact_ns_period_locus_nonempty = (
            pass_equivalence and all_action_gates and period_plane_available
        )
        certified = exact_ns_period_locus_nonempty and example["all_checks_pass"]
        rows.append(
            {
                "seed_id": seed_id,
                "dimension": example["dimension"],
                "same_polytope_input_certified": example["all_checks_pass"],
                "picard_rank": package["rank"],
                "complement_rank": package["k3_complement_rank"],
                "complement_K_rank": complement_rank,
                "twist": twist.to_json(),
                "bounded_discovery_search": search_bound,
                "bounded_discovery_candidates_tested": discovery_tested,
                "bounded_discovery_twist": (
                    discovery_twist.to_json() if discovery_twist is not None else None
                ),
                "bounded_discovery_matches_frozen": (
                    discovery_twist == frozen_twist if search_bound else None
                ),
                "field_determinant": field_determinant.to_json(),
                "trace_scalar_totally_negative": True,
                "field_determinant_signs": [-1, 1],
                "complement_field_diagonal": [
                    complement_field[i][i].to_json() for i in range(complement_rank)
                ],
                "picard_rational_invariants": picard_invariants,
                "complement_rational_invariants": complement_invariants,
                "picard_action_matches_upstream": picard_action_matches_source,
                "picard_action_minimal_polynomial": (
                    picard_action_minimal_polynomial
                ),
                "picard_action_self_adjoint": picard_action_self_adjoint,
                "golden_action_minimal_polynomial": action_square
                == action_plus_identity,
                "golden_action_self_adjoint": self_adjoint,
                "sqrt5_Hodge_operator_square_is_5I": (
                    sqrt5_action_square == five_identity
                ),
                "sqrt5_Hodge_similarity_multiplier": 5 if sqrt5_similarity else None,
                "kuga_satake_isogeny_follows": sqrt5_similarity,
                "X_times_X_algebraic_correspondence_certified": False,
                "kuga_satake_isogeny_used_as_carrier_gate": False,
                "full_H2_action_unital_and_faithful": full_action_minimal_polynomial,
                "full_H2_action_self_adjoint": full_action_self_adjoint,
                "distinguished_eigenspace_signature": identity_eigenspace_signature,
                "conjugate_eigenspace_signature": conjugate_eigenspace_signature,
                "total_rational_invariants": total_invariants,
                "rationally_isometric_to_K3_lattice": pass_equivalence,
                "witt_cancellation_for_every_fixed_primitive_embedding": (
                    pass_equivalence
                ),
                "period_plane_available": period_plane_available,
                "single_real_eigenspace_has_no_nonzero_Q_vectors": True,
                "exact_NS_Hodge_period_locus_nonempty": (
                    exact_ns_period_locus_nonempty
                ),
                "some_finite_index_integral_order_exists": True,
                "maximal_order_integrality_certified": False,
                "explicit_order_conductor_computed": False,
                "status": (
                    "certified_rational_full_H2_Hodge_coherent_same_P_return"
                    if certified
                    else "failed_exact_Hasse_Minkowski_gate"
                ),
            }
        )
    positive = [row for row in rows if str(row["status"]).startswith("certified")]
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "golden_lanner_rational_full_H2_coherent_return",
        "source_same_P_certificate": str(SOURCE.relative_to(ROOT)),
        "scope": (
            "rational Q(sqrt(5)) Hodge action on full H2; denominator clearing "
            "gives some finite-index integral order, but maximal-order integrality, "
            "an explicit conductor, and an algebraic correspondence are not asserted"
        ),
        "method": (
            "exact rational diagonalization, determinant square class, signatures, "
            "and all nontrivial rational Hasse invariants"
        ),
        "k3_rational_invariants": k3,
        "positive_count": len(positive),
        "rows": rows,
        "all_certified_rows_pass": all(
            row["status"]
            == "certified_rational_full_H2_Hodge_coherent_same_P_return"
            for row in rows
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search-bound", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(analyze(args.search_bound), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
