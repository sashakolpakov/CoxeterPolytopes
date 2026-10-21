#!/usr/bin/env python3
"""Block-factorized audit of the silver shell ``x=3+3sqrt(2)``.

Production uses a nine-level leading-Schur recursion.  This replay instead
uses the exact orthogonal decomposition ``B4 + F4 + A1``: it exhausts each
component box directly, groups the middle block by exact field norm, and
joins component vectors whose norms sum to the target.  It therefore shares
neither the production coordinate order nor its pruning tree.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from typing import Any, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver, silver_inertia
from coxeter_polytopes.search.silver_k3_vinberg import (
    Vector,
    _enumeration_data,
    _silver_inverse,
    _silver_mat_vec,
    _silver_pair,
    _silver_quadratic,
    admissible_silver_control_coordinates,
    run_silver_k3_vinberg_four_plus_two_sqrt2_extension,
    run_silver_k3_vinberg_three_plus_three_sqrt2_extension,
)


EXPECTED_REPLAY_SHA256 = "a25c0a841e2ce964f34c6ed0b0a34e585e3f24f9d80d715a760a66c13e8907f4"


def _jsonable(value: Any) -> Any:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Inertia):
        return {
            "positive": value.positive,
            "negative": value.negative,
            "zero": value.zero,
        }
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
class BlockSilverThreePlusThreeShellAudit:
    block_partition: tuple[tuple[int, ...], ...]
    alphabet: tuple[Silver, ...]
    component_box_counts: tuple[int, ...]
    component_bounded_counts: tuple[int, ...]
    norm_join_probe_count: int
    exact_weight_solutions: int
    nonintegral_reconstructions: int
    rejected_root_count: int
    accepted_roots: tuple[Vector, ...]
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    all_prefix_pairings_nonpositive: bool
    prefix_wall_count: int
    prefix_identity_inertia: Inertia
    prefix_conjugate_inertia: Inertia
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    next_admissible_control_coordinate: Silver
    all_reconstructed_roots_sha256: str
    rejection_witnesses_sha256: str
    production_accepted_root_set_matches: bool
    production_rejection_witness_set_matches: bool
    replay_sha256: str

    @property
    def passes(self) -> bool:
        return (
            self.block_partition == ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
            and len(self.alphabet) == 19
            and self.component_box_counts == (19**4, 19**4, 19)
            and self.component_bounded_counts == (166, 38, 19)
            and self.norm_join_probe_count == 3154
            and self.exact_weight_solutions == 154
            and self.nonintegral_reconstructions == 0
            and self.rejected_root_count == 95
            and len(self.accepted_roots) == 59
            and self.same_height_negative_pair_count == 1689
            and self.same_height_zero_pair_count == 22
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 59 * 58 // 2
            and self.all_prefix_pairings_nonpositive
            and self.prefix_wall_count == 333
            and self.prefix_identity_inertia == Inertia(9, 1, 323)
            and self.prefix_conjugate_inertia == Inertia(10, 0, 323)
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and self.next_admissible_control_coordinate == 5 + 2 * SQRT2
            and self.production_accepted_root_set_matches
            and self.production_rejection_witness_set_matches
            and (
                EXPECTED_REPLAY_SHA256 == "TO_BE_FROZEN"
                or self.replay_sha256 == EXPECTED_REPLAY_SHA256
            )
        )


@lru_cache(maxsize=1)
def audit_block_silver_three_plus_three_shell(
) -> BlockSilverThreePlusThreeShellAudit:
    """Enumerate and join the three orthogonal finite components exactly."""

    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, _ = (
        _enumeration_data()
    )
    coordinate = 3 + 3 * SQRT2
    target = 2 + (1 + SQRT2) * coordinate * coordinate

    # The exact square-root uppers at the two real embeddings are 17 and 2.
    # Hence |a|,|b|<=19 is an exhaustive box for z=a+b*sqrt(2).
    alphabet = tuple(
        sorted(
            {
                value
                for a in range(-19, 20)
                for b in range(-19, 20)
                if (value := Silver(a, b)).sign(0) >= 0
                and (2 * target - value * value).sign(0) >= 0
                and (2 * target - value * value).sign(1) >= 0
            }
        )
    )
    blocks = ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
    block_of = {
        coordinate_index: block_index
        for block_index, block in enumerate(blocks)
        for coordinate_index in block
    }
    if any(
        weight_form[i][j]
        for i in range(9)
        for j in range(9)
        if block_of[i] != block_of[j]
    ):
        raise AssertionError("the finite weight form lost its orthogonal blocks")

    bounded: list[tuple[tuple[Vector, Silver], ...]] = []
    for block in blocks:
        block_form = tuple(
            tuple(weight_form[i][j] for j in block) for i in block
        )
        candidates: list[tuple[Vector, Silver]] = []
        for vector in product(alphabet, repeat=len(block)):
            norm = _silver_quadratic(vector, block_form)
            remainder = target - norm
            if remainder.sign(0) >= 0 and remainder.sign(1) >= 0:
                candidates.append((vector, norm))
        bounded.append(tuple(candidates))

    middle_by_norm: dict[Silver, list[Vector]] = {}
    for vector, norm in bounded[1]:
        middle_by_norm.setdefault(norm, []).append(vector)
    weight_solutions: list[Vector] = []
    for left, left_norm in bounded[0]:
        for tail, tail_norm in bounded[2]:
            required = target - left_norm - tail_norm
            for middle in middle_by_norm.get(required, ()):
                weight = (*left, *middle, *tail)
                if _silver_quadratic(weight, weight_form) != target:
                    raise AssertionError("a block join has the wrong full norm")
                weight_solutions.append(weight)
    if len(set(weight_solutions)) != len(weight_solutions):
        raise AssertionError("the block join produced duplicate weights")

    prior = run_silver_k3_vinberg_four_plus_two_sqrt2_extension()
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
            raise AssertionError("the block replay reconstructed a non-root")
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
    combined_roots = (*prior.roots, *accepted_sorted)
    spanning_gram = tuple(
        tuple(
            _silver_pair(combined_roots[i], gram, combined_roots[j])
            for j in range(10)
        )
        for i in range(10)
    )
    spanning_identity = silver_inertia(spanning_gram, embedding=0)
    spanning_conjugate = silver_inertia(spanning_gram, embedding=1)
    if spanning_identity != Inertia(9, 1, 0) or spanning_conjugate != Inertia(
        10, 0, 0
    ):
        raise AssertionError("the independent prefix lost its rank-ten witness")

    edge = prior.edge_indices
    edge_roots = tuple(combined_roots[node] for node in edge)
    edge_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in edge_roots)
        for left in edge_roots
    )
    edge_inverse = _silver_inverse(edge_gram)
    independently_audited_endpoints = tuple(
        index
        for index, root in enumerate(combined_roots)
        if index not in edge
        and (
            residual := Silver(2)
            - _silver_quadratic(
                tuple(_silver_pair(root, gram, wall) for wall in edge_roots),
                edge_inverse,
            )
        ).sign(0)
        > 0
        and residual.sign(1) > 0
    )
    next_coordinates = tuple(
        candidate
        for candidate in admissible_silver_control_coordinates(10)
        if (candidate - coordinate).sign(0) > 0
    )
    if not next_coordinates:
        raise AssertionError("the independent control-coordinate box is too small")

    production = run_silver_k3_vinberg_three_plus_three_sqrt2_extension()
    production_rejected = tuple(
        sorted(
            production.new_level.previous_wall_rejection_witnesses,
            key=lambda item: (_vector_key(item[0]), item[1]),
        )
    )
    payload = {
        "block_partition": blocks,
        "alphabet": alphabet,
        "component_box_counts": tuple(
            len(alphabet) ** len(block) for block in blocks
        ),
        "component_bounded_counts": tuple(len(vectors) for vectors in bounded),
        "norm_join_probe_count": len(bounded[0]) * len(bounded[2]),
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
        # Accepted roots reached this list only after their complete blocker
        # tuples against every prior wall were found empty above.
        "all_prefix_pairings_nonpositive": all(
            pairing.sign(0) <= 0 for pairing in distinct_pairings
        ),
        "prefix_wall_count": len(combined_roots),
        "prefix_identity_inertia": Inertia(
            spanning_identity.positive,
            spanning_identity.negative,
            len(combined_roots) - 10,
        ),
        "prefix_conjugate_inertia": Inertia(
            spanning_conjugate.positive,
            spanning_conjugate.negative,
            len(combined_roots) - 10,
        ),
        "edge_indices": edge,
        "finite_endpoint_root_indices": independently_audited_endpoints,
        "new_finite_endpoint_root_indices": tuple(
            index
            for index in independently_audited_endpoints
            if index >= len(prior.roots)
        ),
        "next_admissible_control_coordinate": next_coordinates[0],
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
    audit = BlockSilverThreePlusThreeShellAudit(
        **payload,
        replay_sha256=replay_sha256,
    )
    if not audit.passes:
        raise AssertionError("the block x=3+3sqrt(2) silver-shell replay failed")
    return audit


def main() -> None:
    audit = audit_block_silver_three_plus_three_shell()
    print(
        json.dumps(
            _jsonable(
                {
                    "accepted_root_count": len(audit.accepted_roots),
                    "all_reconstructed_roots_sha256": (
                        audit.all_reconstructed_roots_sha256
                    ),
                    "component_bounded_counts": audit.component_bounded_counts,
                    "exact_weight_solutions": audit.exact_weight_solutions,
                    "finite_endpoint_root_indices": (
                        audit.finite_endpoint_root_indices
                    ),
                    "next_admissible_control_coordinate": (
                        audit.next_admissible_control_coordinate
                    ),
                    "prefix_conjugate_inertia": audit.prefix_conjugate_inertia,
                    "prefix_identity_inertia": audit.prefix_identity_inertia,
                    "rejected_root_count": audit.rejected_root_count,
                    "rejection_witnesses_sha256": (
                        audit.rejection_witnesses_sha256
                    ),
                    "replay_sha256": audit.replay_sha256,
                }
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
