#!/usr/bin/env python3
"""Reverse-order audit of the silver Vinberg shell ``x=2+2sqrt(2)``.

The production recursion fixes the nine finite-simple-root weights from left
to right.  This replay reverses that order, rebuilds every leading Schur
bound, and reconstructs and classifies all exact solutions independently of
the production traversal and pruning transcript.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from typing import Any, Sequence

from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    Vector,
    _enumeration_data,
    _silver_inverse,
    _silver_mat_vec,
    _silver_pair,
    _silver_quadratic,
    run_silver_k3_vinberg_three_plus_sqrt2_extension,
    run_silver_k3_vinberg_two_plus_two_sqrt2_extension,
)


EXPECTED_REPLAY_SHA256 = "a31ecd05c1d816380ba20f78db78b352d0e617e9e9e2eae6749abc17d9a524f3"


def _jsonable(value: Any) -> Any:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if isinstance(value, tuple):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, list):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(entry) for key, entry in value.items()}
    return value


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _vector_key(vector: Sequence[Silver]) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((entry.a, entry.b) for entry in vector)


@dataclass(frozen=True, slots=True)
class ReverseSilverTwoPlusTwoShellAudit:
    coordinate_permutation: tuple[int, ...]
    alphabet: tuple[Silver, ...]
    visited_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    exact_weight_solutions: int
    nonintegral_reconstructions: int
    rejected_root_count: int
    accepted_roots: tuple[Vector, ...]
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    all_reconstructed_roots_sha256: str
    rejection_witnesses_sha256: str
    production_accepted_root_set_matches: bool
    production_rejection_witness_set_matches: bool
    replay_sha256: str

    @property
    def passes(self) -> bool:
        return (
            self.coordinate_permutation == tuple(reversed(range(9)))
            and len(self.alphabet) == 16
            and self.visited_by_depth
            == (1, 16, 256, 2336, 8816, 6832, 2592, 14160, 37264, 67264)
            and self.schur_pruned_by_depth
            == (0, 0, 110, 1785, 8389, 6670, 1707, 11831, 33060, 65602)
            and self.exact_weight_solutions == 277
            and self.nonintegral_reconstructions == 129
            and self.rejected_root_count == 92
            and len(self.accepted_roots) == 56
            and self.same_height_negative_pair_count == 1482
            and self.same_height_zero_pair_count == 58
            and self.production_accepted_root_set_matches
            and self.production_rejection_witness_set_matches
            and (
                EXPECTED_REPLAY_SHA256 == "TO_BE_FROZEN"
                or self.replay_sha256 == EXPECTED_REPLAY_SHA256
            )
        )


@lru_cache(maxsize=1)
def audit_reverse_silver_two_plus_two_shell(
) -> ReverseSilverTwoPlusTwoShellAudit:
    """Enumerate the shell with the finite weights fixed in reverse order."""

    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, _ = (
        _enumeration_data()
    )
    coordinate = 2 + 2 * SQRT2
    target = 2 + (1 + SQRT2) * coordinate * coordinate

    # The exact two-embedding square bounds are 11 and 2.  Hence |a|,|b|<=13
    # is an exhaustive coefficient box for z=a+b*sqrt(2).
    alphabet = tuple(
        sorted(
            {
                value
                for a in range(-13, 14)
                for b in range(-13, 14)
                if (value := Silver(a, b)).sign(0) >= 0
                and (2 * target - value * value).sign(0) >= 0
                and (2 * target - value * value).sign(1) >= 0
            }
        )
    )
    permutation = tuple(reversed(range(9)))
    finite_simple_gram = _silver_inverse(weight_form)
    permuted_gram = tuple(
        tuple(finite_simple_gram[i][j] for j in permutation)
        for i in permutation
    )
    prefixes = tuple(
        _silver_inverse(
            tuple(
                tuple(permuted_gram[i][j] for j in range(size))
                for i in range(size)
            )
        )
        for size in range(1, 10)
    )

    visited = [0 for _ in range(10)]
    pruned = [0 for _ in range(10)]
    weight_solutions: list[Vector] = []

    def extend(prefix: Vector) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            remainder = target - _silver_quadratic(prefix, prefixes[depth - 1])
            if remainder.sign(0) < 0 or remainder.sign(1) < 0:
                pruned[depth] += 1
                return
        if depth == 9:
            weight = [SILVER_ZERO for _ in range(9)]
            for permuted_index, entry in enumerate(prefix):
                weight[permutation[permuted_index]] = entry
            candidate = tuple(weight)
            if _silver_quadratic(candidate, weight_form) == target:
                weight_solutions.append(candidate)
            return
        for entry in alphabet:
            extend((*prefix, entry))

    extend(())

    prior = run_silver_k3_vinberg_three_plus_sqrt2_extension()
    reconstructed: list[Vector] = []
    nonintegral = 0
    accepted: list[Vector] = []
    rejected: list[tuple[Vector, tuple[int, ...]]] = []
    for weight in weight_solutions:
        finite_coordinates = _silver_mat_vec(
            finite_inverse,
            tuple(
                entry - cross[index] * coordinate
                for index, entry in enumerate(
                    _silver_mat_vec(inverse_simple_transpose, weight)
                )
            ),
        )
        if any(not entry.is_integral for entry in finite_coordinates):
            nonintegral += 1
            continue
        root = (coordinate, *finite_coordinates)
        if _silver_quadratic(root, gram) != 2:
            raise AssertionError("the reverse replay reconstructed a non-root")
        reconstructed.append(root)
        blockers = tuple(
            index
            for index, wall in enumerate(prior.roots)
            if _silver_pair(root, gram, wall).sign(0) > 0
        )
        if blockers:
            rejected.append((root, blockers))
        else:
            accepted.append(root)

    accepted_sorted = tuple(sorted(accepted, key=_vector_key))
    rejected_sorted = tuple(
        sorted(rejected, key=lambda item: (_vector_key(item[0]), item[1]))
    )
    all_roots_sorted = tuple(sorted(reconstructed, key=_vector_key))
    distinct_pairings = tuple(
        _silver_pair(accepted_sorted[i], gram, accepted_sorted[j])
        for i in range(len(accepted_sorted))
        for j in range(i)
    )
    production = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
    production_rejected = tuple(
        sorted(
            production.new_level.previous_wall_rejection_witnesses,
            key=lambda item: (_vector_key(item[0]), item[1]),
        )
    )

    payload = {
        "coordinate_permutation": permutation,
        "alphabet": alphabet,
        "visited_by_depth": tuple(visited),
        "schur_pruned_by_depth": tuple(pruned),
        "exact_weight_solutions": len(weight_solutions),
        "nonintegral_reconstructions": nonintegral,
        "rejected_root_count": len(rejected),
        "accepted_roots": accepted_sorted,
        "same_height_negative_pair_count": sum(
            pairing.sign(0) < 0 for pairing in distinct_pairings
        ),
        "same_height_zero_pair_count": sum(
            pairing.sign(0) == 0 for pairing in distinct_pairings
        ),
        "all_reconstructed_roots_sha256": _canonical_sha256(
            _jsonable(all_roots_sorted)
        ),
        "rejection_witnesses_sha256": _canonical_sha256(
            _jsonable(rejected_sorted)
        ),
        "production_accepted_root_set_matches": accepted_sorted
        == tuple(sorted(production.new_level.roots, key=_vector_key)),
        "production_rejection_witness_set_matches": rejected_sorted
        == production_rejected,
    }
    replay_sha256 = _canonical_sha256(_jsonable(payload))
    audit = ReverseSilverTwoPlusTwoShellAudit(
        **payload,
        replay_sha256=replay_sha256,
    )
    if not audit.passes:
        raise AssertionError("the reverse x=2+2sqrt(2) silver-shell replay failed")
    return audit


def main() -> None:
    audit = audit_reverse_silver_two_plus_two_shell()
    print(
        json.dumps(
            {
                "accepted_root_count": len(audit.accepted_roots),
                "all_reconstructed_roots_sha256": (
                    audit.all_reconstructed_roots_sha256
                ),
                "exact_weight_solutions": audit.exact_weight_solutions,
                "rejected_root_count": audit.rejected_root_count,
                "rejection_witnesses_sha256": audit.rejection_witnesses_sha256,
                "replay_sha256": audit.replay_sha256,
                "schur_pruned_by_depth": audit.schur_pruned_by_depth,
                "visited_by_depth": audit.visited_by_depth,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
