#!/usr/bin/env python3
"""Block-factorized audit of the silver shell ``x=4+3sqrt(2)``.

Production uses a nine-level leading-Schur recursion.  This replay instead
exhausts the exact orthogonal ``B4 + F4 + A1`` component boxes and joins them
by exact field norm, so it shares neither production order nor pruning tree.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
import sys
from typing import Any, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ZERO,
    SQRT2,
    Silver,
    silver_inertia,
)
from coxeter_polytopes.search.silver_k3_vinberg import (
    Vector,
    _enumeration_data,
    _silver_inverse,
    _silver_dot,
    _silver_mat_vec,
    _silver_pair,
    _silver_quadratic,
    admissible_silver_control_coordinates,
    run_silver_k3_vinberg_five_plus_two_sqrt2_extension,
    run_silver_k3_vinberg_four_plus_three_sqrt2_extension,
)


EXPECTED_REPLAY_SHA256 = "87313bb8f9911ba9a6102c622475cfe09a6940685cafe32f42ba57f0ab449adf"


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


def _phase(message: str) -> None:
    print(f"[independent x=4+3sqrt(2)] {message}", file=sys.stderr, flush=True)


def _direct_bounded_block(
    alphabet: tuple[Silver, ...],
    block_form: tuple[tuple[Silver, ...], ...],
    target: Silver,
) -> tuple[tuple[Vector, Silver], ...]:
    """Exhaust a raw component box with tabulated exact monomials.

    This is still a direct Cartesian-product census.  Tabulation merely avoids
    repeating the same field multiplications during certificate replays; it
    neither prunes prefixes nor imports production's Schur recursion.
    """

    size = len(block_form)
    diagonal = tuple(
        {entry: block_form[i][i] * entry * entry for entry in alphabet}
        for i in range(size)
    )
    cross = {
        (i, j): {
            (left, right): 2 * block_form[i][j] * left * right
            for left in alphabet
            for right in alphabet
        }
        for i in range(size)
        for j in range(i)
    }
    bounded: list[tuple[Vector, Silver]] = []
    for vector in product(alphabet, repeat=size):
        norm = SILVER_ZERO
        for i, entry in enumerate(vector):
            norm += diagonal[i][entry]
            for j in range(i):
                norm += cross[i, j][entry, vector[j]]
        remainder = target - norm
        if remainder.sign(0) >= 0 and remainder.sign(1) >= 0:
            bounded.append((vector, norm))
    return tuple(bounded)


@dataclass(frozen=True, slots=True)
class BlockSilverFourPlusThreeShellAudit:
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
    native_finite_endpoint_root_indices: tuple[int, ...]
    native_new_finite_endpoint_root_indices: tuple[int, ...]
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
            and len(self.alphabet) == 26
            and self.component_box_counts == (26**4, 26**4, 26)
            and self.component_bounded_counts == (488, 115, 26)
            and self.norm_join_probe_count == 12688
            and self.exact_weight_solutions == 2507
            and self.nonintegral_reconstructions == 1544
            and self.rejected_root_count == 570
            and len(self.accepted_roots) == 393
            and self.same_height_negative_pair_count == 76453
            and self.same_height_zero_pair_count == 575
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 393 * 392 // 2
            and self.all_prefix_pairings_nonpositive
            and self.prefix_wall_count == 726
            and self.prefix_identity_inertia == Inertia(9, 1, 716)
            and self.prefix_conjugate_inertia == Inertia(10, 0, 716)
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33, 622)
            and self.new_finite_endpoint_root_indices == (622,)
            and self.native_finite_endpoint_root_indices == (33, 443)
            and self.native_new_finite_endpoint_root_indices == (443,)
            and self.next_admissible_control_coordinate == 5 + 3 * SQRT2
            and self.all_reconstructed_roots_sha256
            == "60cf05c9a8618f1154503190b42d23f98b791aad59a53ac35899aaae47699e82"
            and self.rejection_witnesses_sha256
            == "351afd9d5fb1d1f383511d283d7f0552d5044a28d929cb0648f4c1721a474ce7"
            and self.production_accepted_root_set_matches
            and self.production_rejection_witness_set_matches
            and (
                EXPECTED_REPLAY_SHA256 == "TO_BE_FROZEN"
                or self.replay_sha256 == EXPECTED_REPLAY_SHA256
            )
        )


@lru_cache(maxsize=1)
def audit_block_silver_four_plus_three_shell(
) -> BlockSilverFourPlusThreeShellAudit:
    """Enumerate and join the three orthogonal finite components exactly."""

    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, _ = (
        _enumeration_data()
    )
    coordinate = 4 + 3 * SQRT2
    target = 2 + (1 + SQRT2) * coordinate * coordinate

    # Exact integral square-root uppers are 19 and 2, hence the exhaustive
    # coefficient box |a|,|b|<=21 for z=a+b*sqrt(2).
    alphabet = tuple(
        sorted(
            {
                value
                for a in range(-21, 22)
                for b in range(-21, 22)
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
        bounded.append(_direct_bounded_block(alphabet, block_form, target))
    _phase(
        "component boxes exhausted; bounded counts="
        + repr(tuple(len(vectors) for vectors in bounded))
    )

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
    _phase(f"norm join complete; exact weights={len(weight_solutions)}")

    prior = run_silver_k3_vinberg_five_plus_two_sqrt2_extension()
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
        gram_root = _silver_mat_vec(gram, root)
        blockers = tuple(
            index
            for index, wall in enumerate(prior.roots)
            if _silver_dot(wall, gram_root).sign(0) > 0
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
    _phase(
        "integral reconstruction/blockers complete; "
        f"nonintegral={nonintegral}, rejected={len(rejected)}, "
        f"accepted={len(accepted_sorted)}"
    )
    accepted_gram_images = tuple(
        _silver_mat_vec(gram, root) for root in accepted_sorted
    )
    distinct_pairings = tuple(
        _silver_dot(accepted_sorted[i], accepted_gram_images[j])
        for i in range(len(accepted_sorted))
        for j in range(i)
    )
    _phase(f"same-height census complete; pairings={len(distinct_pairings)}")
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
    independently_audited_native_endpoints = tuple(
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

    production = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    _phase("production extension complete; comparing exact sets")
    production_accepted_sorted = tuple(
        sorted(production.new_level.roots, key=_vector_key)
    )
    production_accepted_root_set_matches = (
        accepted_sorted == production_accepted_sorted
    )
    if not production_accepted_root_set_matches:
        raise AssertionError(
            "the independent accepted-root set differs from production"
        )
    production_index_by_root = {
        root: index for index, root in enumerate(production.roots)
    }
    if len(production_index_by_root) != len(production.roots):
        raise AssertionError("the production Vinberg prefix contains duplicates")
    # The independent block join deliberately has its own lexicographic root
    # order.  Endpoint roots are invariants; translate those vectors into the
    # certified production ordering only after exact set equality is known.
    production_order_endpoints = tuple(
        production_index_by_root[combined_roots[index]]
        for index in independently_audited_native_endpoints
    )
    production_order_new_endpoints = tuple(
        index
        for index in production_order_endpoints
        if index >= len(prior.roots)
    )
    production_rejected = tuple(
        sorted(
            production.new_level.previous_wall_rejection_witnesses,
            key=lambda item: (_vector_key(item[0]), item[1]),
        )
    )
    native_replay_payload = {
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
        "finite_endpoint_root_indices": independently_audited_native_endpoints,
        "new_finite_endpoint_root_indices": tuple(
            index
            for index in independently_audited_native_endpoints
            if index >= len(prior.roots)
        ),
        "next_admissible_control_coordinate": next_coordinates[0],
        "all_reconstructed_roots_sha256": _canonical_sha256(
            _jsonable(all_roots_sorted)
        ),
        "rejection_witnesses_sha256": _canonical_sha256(
            _jsonable(rejected_sorted)
        ),
        "production_accepted_root_set_matches": (
            production_accepted_root_set_matches
        ),
        "production_rejection_witness_set_matches": rejected_sorted
        == production_rejected,
    }
    # Preserve the frozen replay digest over the replay's native ordering.
    replay_sha256 = _canonical_sha256(_jsonable(native_replay_payload))
    payload = dict(native_replay_payload)
    payload["finite_endpoint_root_indices"] = production_order_endpoints
    payload["new_finite_endpoint_root_indices"] = production_order_new_endpoints
    audit = BlockSilverFourPlusThreeShellAudit(
        **payload,
        native_finite_endpoint_root_indices=(
            independently_audited_native_endpoints
        ),
        native_new_finite_endpoint_root_indices=tuple(
            index
            for index in independently_audited_native_endpoints
            if index >= len(prior.roots)
        ),
        replay_sha256=replay_sha256,
    )
    if not audit.passes:
        raise AssertionError(
            "the block x=4+3sqrt(2) silver-shell replay failed: "
            + json.dumps(
                _jsonable(
                    {
                        "component_bounded_counts": audit.component_bounded_counts,
                        "norm_join_probe_count": audit.norm_join_probe_count,
                        "exact_weight_solutions": audit.exact_weight_solutions,
                        "nonintegral_reconstructions": (
                            audit.nonintegral_reconstructions
                        ),
                        "rejected_root_count": audit.rejected_root_count,
                        "accepted_root_count": len(audit.accepted_roots),
                        "same_height_negative_pair_count": (
                            audit.same_height_negative_pair_count
                        ),
                        "same_height_zero_pair_count": (
                            audit.same_height_zero_pair_count
                        ),
                        "prefix_wall_count": audit.prefix_wall_count,
                        "finite_endpoint_root_indices": (
                            audit.finite_endpoint_root_indices
                        ),
                        "all_reconstructed_roots_sha256": (
                            audit.all_reconstructed_roots_sha256
                        ),
                        "rejection_witnesses_sha256": (
                            audit.rejection_witnesses_sha256
                        ),
                        "production_accepted_root_set_matches": (
                            audit.production_accepted_root_set_matches
                        ),
                        "production_rejection_witness_set_matches": (
                            audit.production_rejection_witness_set_matches
                        ),
                        "replay_sha256": audit.replay_sha256,
                    }
                ),
                sort_keys=True,
            )
        )
    return audit


def main() -> None:
    audit = audit_block_silver_four_plus_three_shell()
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
                    "nonintegral_reconstructions": (
                        audit.nonintegral_reconstructions
                    ),
                    "prefix_conjugate_inertia": audit.prefix_conjugate_inertia,
                    "prefix_identity_inertia": audit.prefix_identity_inertia,
                    "rejected_root_count": audit.rejected_root_count,
                    "rejection_witnesses_sha256": (
                        audit.rejection_witnesses_sha256
                    ),
                    "replay_sha256": audit.replay_sha256,
                    "same_height_negative_pair_count": (
                        audit.same_height_negative_pair_count
                    ),
                    "same_height_zero_pair_count": (
                        audit.same_height_zero_pair_count
                    ),
                }
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
