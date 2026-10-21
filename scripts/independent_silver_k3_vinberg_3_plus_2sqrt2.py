#!/usr/bin/env python3
"""Block-factorized audit of the silver shell ``x=3+2sqrt(2)``.

Production uses a nine-level leading-Schur recursion.  This replay instead
exploits the exact orthogonal decomposition ``B4 + F4 + A1``: it exhausts the
three component boxes directly, groups the middle block by exact field norm,
and joins component vectors whose norms sum to the target.  It therefore has
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

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    Vector,
    _enumeration_data,
    _silver_mat_vec,
    _silver_pair,
    _silver_quadratic,
    run_silver_k3_vinberg_three_plus_two_sqrt2_extension,
    run_silver_k3_vinberg_two_plus_two_sqrt2_extension,
)


EXPECTED_REPLAY_SHA256 = "f45bcd9e02db07db2fb0c7a0bb9fe0363b7490803d45da5bc697e8fc8ae65600"


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
class BlockSilverThreePlusTwoShellAudit:
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
    all_reconstructed_roots_sha256: str
    rejection_witnesses_sha256: str
    production_accepted_root_set_matches: bool
    production_rejection_witness_set_matches: bool
    replay_sha256: str

    @property
    def passes(self) -> bool:
        return (
            self.block_partition == ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
            and len(self.alphabet) == 18
            and self.component_box_counts == (18**4, 18**4, 18)
            and self.component_bounded_counts == (154, 38, 18)
            and self.norm_join_probe_count == 2772
            and self.exact_weight_solutions == 138
            and self.nonintegral_reconstructions == 0
            and self.rejected_root_count == 78
            and len(self.accepted_roots) == 60
            and self.same_height_negative_pair_count == 1702
            and self.same_height_zero_pair_count == 68
            and self.production_accepted_root_set_matches
            and self.production_rejection_witness_set_matches
            and self.replay_sha256 == EXPECTED_REPLAY_SHA256
        )


@lru_cache(maxsize=1)
def audit_block_silver_three_plus_two_shell(
) -> BlockSilverThreePlusTwoShellAudit:
    """Enumerate and join the orthogonal finite components exactly."""

    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, _ = (
        _enumeration_data()
    )
    coordinate = 3 + 2 * SQRT2
    target = 2 + (1 + SQRT2) * coordinate * coordinate

    # The two exact square bounds are 13 and 2.  Thus |a|,|b|<=15 is an
    # exhaustive coefficient box for every z=a+b*sqrt(2).
    alphabet = tuple(
        sorted(
            {
                value
                for a in range(-15, 16)
                for b in range(-15, 16)
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

    prior = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
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
    production = run_silver_k3_vinberg_three_plus_two_sqrt2_extension()
    production_rejected = tuple(
        sorted(
            production.new_level.previous_wall_rejection_witnesses,
            key=lambda item: (_vector_key(item[0]), item[1]),
        )
    )

    payload = {
        "block_partition": blocks,
        "alphabet": alphabet,
        "component_box_counts": tuple(len(alphabet) ** len(block) for block in blocks),
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
    audit = BlockSilverThreePlusTwoShellAudit(
        **payload,
        replay_sha256=replay_sha256,
    )
    if not audit.passes:
        raise AssertionError("the block x=3+2sqrt(2) silver-shell replay failed")
    return audit


def main() -> None:
    audit = audit_block_silver_three_plus_two_shell()
    print(
        json.dumps(
            {
                "accepted_root_count": len(audit.accepted_roots),
                "all_reconstructed_roots_sha256": (
                    audit.all_reconstructed_roots_sha256
                ),
                "component_bounded_counts": audit.component_bounded_counts,
                "exact_weight_solutions": audit.exact_weight_solutions,
                "rejected_root_count": audit.rejected_root_count,
                "rejection_witnesses_sha256": audit.rejection_witnesses_sha256,
                "replay_sha256": audit.replay_sha256,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
