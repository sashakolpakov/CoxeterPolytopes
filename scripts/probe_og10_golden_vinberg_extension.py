#!/usr/bin/env python3
"""Extend and diagnose the archived exact OG10 norm-two Vinberg prefix.

The bound-eight artifact ends at ``x=3+3*phi``.  This script resumes the
same exhaustive height enumerator, preserving the archived ordering, and
can freeze one or both of the remaining control levels below identity bound
ten.  A hyperideal witness is searched locally at the archived ray after
the new walls are inserted and is then checked against every recorded wall.

This is still only a bounded norm-two prefix: neither later norm-two roots
nor the other OG10 wall-divisor shells are excluded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_vinberg import (
    og10_field_gram,
    run_og10_golden_vinberg_extension,
)
from coxeter_polytopes.verify.chamber import (
    CertificateError,
    _canonical_ray,
    _dot,
    _full_rank_kernel,
    _mat_vec,
    verify_hyperideal_extreme_ray,
)


ROOT = Path(__file__).resolve().parents[1]
BOUND_EIGHT_ARTIFACT = ROOT / "data/candidates/og10_golden_vinberg_bound8.json"
FIRST_EXTENSION_ARTIFACT = (
    ROOT / "data/candidates/og10_golden_vinberg_through_2_plus_4phi.json"
)
IDENTITY_BOUND = 10
DIMENSION = 10


def _vector(vector: Iterable[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _roots(value: object) -> tuple[tuple[Golden, ...], ...]:
    if not isinstance(value, list):
        raise TypeError("serialized roots must be a list")
    return tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in value
    )


def _digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _relation(inner: Golden) -> str:
    if inner == 0:
        return "m2"
    if inner == -1:
        return "m3"
    if inner == -PHI:
        return "m5"
    if (inner + 2).sign(0) <= 0:
        return "ultraparallel"
    raise ArithmeticError(f"unrecognized acute norm-two product {inner}")


def _local_hyperideal_witness(
    roots: Sequence[Sequence[Golden]],
    old_active: Sequence[int],
    old_ray: Sequence[Golden],
    first_new_index: int,
) -> tuple[tuple[int, ...], tuple[Golden, ...], Golden, tuple[Golden, ...]] | None:
    """Search rays created by cutting the archived hyperideal ray.

    If a new wall cuts an old extreme ray, every new ray adjacent to it is
    cut out by that wall and nine independent walls through the old ray.
    New walls containing the old ray are included in that local active pool.
    Every candidate is nevertheless checked against *all* walls exactly.
    """

    gram = og10_field_gram()
    functionals = [_mat_vec(gram, root) for root in roots]
    old_pairings = tuple(_pair(old_ray, gram, root) for root in roots)
    cutting = tuple(
        index
        for index in range(first_new_index, len(roots))
        if old_pairings[index].sign(0) > 0
    )
    containing = tuple(
        index
        for index in range(first_new_index, len(roots))
        if not old_pairings[index]
    )
    local_pool = tuple(dict.fromkeys((*old_active, *containing)))

    if not cutting:
        canonical = _canonical_ray(old_ray)
        canonical_pairings = tuple(_pair(canonical, gram, root) for root in roots)
        if all(value.sign(0) <= 0 for value in canonical_pairings):
            norm = _pair(canonical, gram, canonical)
            actual_active = tuple(
                index for index, value in enumerate(canonical_pairings) if not value
            )
            if norm.sign(0) > 0:
                return actual_active, canonical, norm, canonical_pairings

    witnesses: dict[
        tuple[Golden, ...],
        tuple[tuple[int, ...], tuple[Golden, ...], Golden, tuple[Golden, ...]],
    ] = {}
    for cutting_index in cutting:
        for inherited in combinations(local_pool, DIMENSION - 1):
            active_seed = (cutting_index, *inherited)
            ray = _full_rank_kernel([functionals[index] for index in active_seed])
            if ray is None:
                continue
            pairings = tuple(_dot(ray, functional) for functional in functionals)
            if all(value.sign(0) <= 0 for value in pairings):
                oriented = ray
                oriented_pairings = pairings
            elif all(value.sign(0) >= 0 for value in pairings):
                oriented = [-entry for entry in ray]
                oriented_pairings = tuple(-value for value in pairings)
            else:
                continue
            canonical = _canonical_ray(oriented)
            norm = _pair(canonical, gram, canonical)
            if norm.sign(0) <= 0:
                continue
            canonical_pairings = tuple(_pair(canonical, gram, root) for root in roots)
            actual_active = tuple(
                index for index, value in enumerate(canonical_pairings) if not value
            )
            witnesses[canonical] = (
                actual_active,
                canonical,
                norm,
                canonical_pairings,
            )
    if not witnesses:
        return None
    # Prefer a simple ray accepted directly by the independent verifier, then
    # choose deterministically by active-facet tuple.
    ordered = sorted(
        witnesses.values(),
        key=lambda item: (len(item[0]) != DIMENSION, item[0]),
    )
    return ordered[0]


def build(max_new_levels: int = 1) -> dict[str, object]:
    if max_new_levels not in {1, 2}:
        raise ValueError("max_new_levels must be one or two")
    archived = json.loads(BOUND_EIGHT_ARTIFACT.read_text())
    archived_roots = _roots(archived["simple_roots"])
    cone_count = archived["height_zero"]["simple_root_count"]
    start_after = Golden.from_json(archived["candidate_control_coordinates"][-1])
    if max_new_levels == 1:
        previous_roots = archived_roots
        previous_witness_payload = archived["hyperideal_witness"]
        resume_after = start_after
        extension = run_og10_golden_vinberg_extension(
            archived_roots[cone_count:],
            start_after=resume_after,
            identity_bound=IDENTITY_BOUND,
            max_new_levels=1,
        )
        all_candidate_coordinates = extension.candidate_coordinates
        all_levels_prefix: list[dict[str, object]] = []
        all_new_roots_prefix: tuple[tuple[Golden, ...], ...] = ()
        earlier_visited = 0
    else:
        if not FIRST_EXTENSION_ARTIFACT.exists():
            raise FileNotFoundError(
                "generate the one-level extension artifact before resuming level two"
            )
        previous_payload = json.loads(FIRST_EXTENSION_ARTIFACT.read_text())
        if previous_payload.get("verified") is not True:
            raise ValueError("the one-level extension checkpoint is not verified")
        previous_roots = _roots(previous_payload["simple_roots"])
        if len(previous_roots) != previous_payload.get("facet_count"):
            raise ValueError("the one-level checkpoint has the wrong facet count")
        if _digest([_vector(root) for root in previous_roots]) != previous_payload.get(
            "roots_sha256"
        ):
            raise ValueError("the one-level checkpoint root digest does not match")
        previous_witness_payload = previous_payload["hyperideal_witness"]
        checkpoint_witness = verify_hyperideal_extreme_ray(
            previous_roots,
            og10_field_gram(),
            previous_witness_payload["active_facets"],
            dimension=DIMENSION,
        )
        if _digest(
            [value.to_json() for value in checkpoint_witness.wall_pairings]
        ) != previous_witness_payload.get("wall_pairings_sha256"):
            raise ValueError("the one-level checkpoint witness digest does not match")
        resume_after = Golden.from_json(
            previous_payload["completed_through_control_coordinate"]
        )
        extension = run_og10_golden_vinberg_extension(
            previous_roots[cone_count:],
            start_after=resume_after,
            identity_bound=IDENTITY_BOUND,
            max_new_levels=1,
        )
        all_candidate_coordinates = (
            *tuple(
                Golden.from_json(value)
                for value in previous_payload["candidate_control_coordinates"]
            ),
            *extension.candidate_coordinates,
        )
        all_levels_prefix = list(previous_payload["levels"])
        all_new_roots_prefix = previous_roots[len(archived_roots) :]
        earlier_visited = int(previous_payload["extension_visited_node_count"])

    immediate_new_roots = extension.new_roots
    new_roots = (*all_new_roots_prefix, *immediate_new_roots)
    roots = (*previous_roots, *immediate_new_roots)
    roots_json = [_vector(root) for root in roots]
    gram = og10_field_gram()

    old_witness_payload = previous_witness_payload
    old_ray = tuple(
        Golden.from_json(entry) for entry in old_witness_payload["ray"]
    )
    new_pairings_at_old_ray = tuple(
        _pair(old_ray, gram, root) for root in immediate_new_roots
    )
    old_ray_signs = Counter(
        "positive" if value.sign(0) > 0 else "negative" if value.sign(0) < 0 else "zero"
        for value in new_pairings_at_old_ray
    )

    local = _local_hyperideal_witness(
        roots,
        old_witness_payload["active_facets"],
        old_ray,
        len(previous_roots),
    )
    if local is None:
        witness_payload: dict[str, object] | None = None
        witness_verified = False
    else:
        active, ray, norm, pairings = local
        witness_payload = {
            "active_facets": list(active),
            "ray": _vector(ray),
            "norm": norm.to_json(),
            "wall_pairings_sha256": _digest(
                [value.to_json() for value in pairings]
            ),
        }
        try:
            independently_verified = verify_hyperideal_extreme_ray(
                roots, gram, active, dimension=DIMENSION
            )
        except CertificateError:
            witness_verified = False
        else:
            witness_verified = (
                independently_verified.ray == ray
                and independently_verified.norm == norm
                and independently_verified.wall_pairings == pairings
            )

    root_functionals = tuple(_mat_vec(gram, root) for root in roots)
    relation_counts = Counter(
        _relation(_dot(left, root_functionals[right_index]))
        for left_index, left in enumerate(roots)
        for right_index in range(left_index + 1, len(roots))
    )
    current_levels = [
        {
            "control_coordinate": level.control_coordinate.to_json(),
            "height": level.height.to_json(),
            "target": level.target.to_json(),
            "root_count": len(level.roots),
            "visited_by_depth": list(level.visited_by_depth),
            "parity_pruned": level.parity_pruned,
            "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
            "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
            "leaf_norm_mismatches": level.leaf_norm_mismatches,
        }
        for level in extension.levels
    ]
    levels = [*all_levels_prefix, *current_levels]
    expected_first_level = (
        len(levels) >= 1
        and Golden.from_json(levels[0]["control_coordinate"]) == Golden(2, 4)
        and levels[0]["root_count"] == 150
    )
    checks = {
        "archived_bound_eight_prefix_verified": archived["verified"] is True,
        "candidate_levels_complete_in_stated_interval": (
            [str(value) for value in all_candidate_coordinates]
            == (["2 + 4*phi"] if max_new_levels == 1 else ["2 + 4*phi", "3 + 4*phi"])
        ),
        "first_extension_level_has_150_roots": expected_first_level,
        "roots_distinct": len(roots) == len(set(roots)),
        "all_integral_norm_two": all(
            all(entry.is_integral for entry in root) and _pair(root, gram, root) == 2
            for root in roots
        ),
        "all_pairs_are_Coxeter_or_ultraparallel": sum(relation_counts.values())
        == len(roots) * (len(roots) - 1) // 2,
        "immediate_new_wall_sign_partition_complete": sum(old_ray_signs.values())
        == len(immediate_new_roots),
        "resulting_hyperideal_ray_found": witness_payload is not None,
        "resulting_hyperideal_ray_independently_verified": witness_verified,
    }
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "bounded_exact_vinberg_extension_with_hyperideal_witness",
        "candidate_id": "og10_golden_h10_carrier",
        "hyperbolic_dimension": DIMENSION,
        "candidate_generation_identity_bound": IDENTITY_BOUND,
        "extension_start_after": start_after.to_json(),
        "completed_through_control_coordinate": all_candidate_coordinates[-1].to_json(),
        "candidate_control_coordinates": [
            value.to_json() for value in all_candidate_coordinates
        ],
        "levels": levels,
        "extension_visited_node_count": earlier_visited + sum(
            sum(level.visited_by_depth) for level in extension.levels
        ),
        "archived_facet_count": len(archived_roots),
        "immediate_previous_facet_count": len(previous_roots),
        "new_facet_count": len(new_roots),
        "immediate_new_facet_count": len(immediate_new_roots),
        "facet_count": len(roots),
        "previous_hyperideal_ray_immediate_new_wall_signs": dict(
            sorted(old_ray_signs.items())
        ),
        "previous_hyperideal_ray_status": (
            "cut"
            if old_ray_signs["positive"]
            else "survives_with_new_active_walls"
            if old_ray_signs["zero"]
            else "survives_strictly"
        ),
        "pair_relation_counts": dict(sorted(relation_counts.items())),
        "roots_sha256": _digest(roots_json),
        "simple_roots": roots_json,
        "hyperideal_witness": witness_payload,
        "checks": checks,
        "verified": all(checks.values()),
        "conclusion": (
            "the_recorded_norm_two_outer_cone_is_not_finite_volume;_later_"
            "norm_two_roots_may_cut_the_witness;_other_OG10_wall_types_are_"
            "separate_Kahler_subdivision_data"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-new-levels", type=int, choices=(1, 2), default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    payload = build(args.max_new_levels)
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.summary_only:
        summary = {key: value for key, value in payload.items() if key != "simple_roots"}
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
