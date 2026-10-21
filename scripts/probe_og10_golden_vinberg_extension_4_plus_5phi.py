#!/usr/bin/env python3
"""Freeze the exact OG10 norm-two Vinberg prefix through ``4+5*phi``.

The archived bound-ten prefix ends at ``x=3+4*phi`` with 457 walls.  The
next four admissible control levels returned with generation bound 13 are

    2+5*phi, 4+4*phi, 3+5*phi, 4+5*phi.

This script exhausts those four levels one at a time, archives all recursion
statistics, and propagates an exact spacelike extreme-ray witness after each
cut.  It deliberately stops before the fifth candidate ``3+6*phi``, whose
identity value is still below 13.  The result is therefore a bounded prefix
theorem through ``4+5*phi``, not a theorem through all ``x<13``.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
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
from scripts.probe_og10_golden_vinberg_extension import (
    _digest,
    _local_hyperideal_witness,
    _relation,
    _roots,
    _vector,
)
from scripts import probe_og10_golden_chamber_folding as folding


ROOT = Path(__file__).resolve().parents[1]
BOUND_TEN = ROOT / "data/candidates/og10_golden_vinberg_through_3_plus_4phi.json"
IDENTITY_BOUND = 13
DIMENSION = 10
EXPECTED_COORDINATES = (
    Golden(2, 5),
    Golden(4, 4),
    Golden(3, 5),
    Golden(4, 5),
)
EXPECTED_COUNTS = (71, 13, 371, 232)


def _level_payload(level) -> dict[str, object]:
    return {
        "control_coordinate": level.control_coordinate.to_json(),
        "height": level.height.to_json(),
        "target": level.target.to_json(),
        "root_count": len(level.roots),
        "roots_sha256": _digest([_vector(root) for root in level.roots]),
        "visited_by_depth": list(level.visited_by_depth),
        "visited_node_count": sum(level.visited_by_depth),
        "parity_pruned": level.parity_pruned,
        "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
        "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
        "leaf_norm_mismatches": level.leaf_norm_mismatches,
    }


def _witness_payload(
    roots,
    old_witness: dict[str, object],
    first_new_index: int,
) -> tuple[dict[str, object] | None, Counter[str]]:
    gram = og10_field_gram()
    old_ray = tuple(Golden.from_json(entry) for entry in old_witness["ray"])
    new_pairings = tuple(
        _pair(old_ray, gram, root) for root in roots[first_new_index:]
    )
    signs = Counter(
        "positive"
        if value.sign() > 0
        else "negative"
        if value.sign() < 0
        else "zero"
        for value in new_pairings
    )
    functionals = tuple(_mat_vec(gram, root) for root in roots)
    cutting = tuple(
        index
        for index in range(first_new_index, len(roots))
        if _pair(old_ray, gram, roots[index]).sign() > 0
    )
    local = None
    # Generic truncation of a ten-facet extreme ray: one cutting wall and
    # nine of the old active walls.  Trying this small exact set first avoids
    # a combinatorial explosion when many new walls also contain the old ray.
    for cutting_index in cutting:
        for inherited in combinations(old_witness["active_facets"], DIMENSION - 1):
            candidate = _full_rank_kernel(
                [
                    functionals[index]
                    for index in (cutting_index, *inherited)
                ]
            )
            if candidate is None:
                continue
            candidate_pairings = tuple(
                _dot(candidate, functional) for functional in functionals
            )
            if all(value.sign() <= 0 for value in candidate_pairings):
                oriented = candidate
            elif all(value.sign() >= 0 for value in candidate_pairings):
                oriented = [-entry for entry in candidate]
            else:
                continue
            ray = _canonical_ray(oriented)
            pairings = tuple(_pair(ray, gram, root) for root in roots)
            norm = _pair(ray, gram, ray)
            if norm.sign() <= 0:
                continue
            active = tuple(i for i, value in enumerate(pairings) if not value)
            local = (active, ray, norm, pairings)
            break
        if local is not None:
            break
    if local is None:
        # Degenerate fallback: include new walls containing the old ray, as
        # in the independently used bound-ten reconstruction.
        local = _local_hyperideal_witness(
            roots,
            old_witness["active_facets"],
            old_ray,
            first_new_index,
        )
    if local is None:
        return None, signs
    active, ray, norm, pairings = local
    payload = {
        "active_facets": list(active),
        "ray": _vector(ray),
        "norm": norm.to_json(),
        "wall_pairings_sha256": _digest(
            [value.to_json() for value in pairings]
        ),
    }
    try:
        independent = verify_hyperideal_extreme_ray(
            roots, gram, active, dimension=DIMENSION
        )
    except CertificateError:
        payload["independently_verified"] = False
    else:
        payload["independently_verified"] = (
            independent.ray == ray
            and independent.norm == norm
            and independent.wall_pairings == pairings
        )
    return payload, signs


def build(*, progress: bool = False) -> dict[str, object]:
    archived = json.loads(BOUND_TEN.read_text())
    archived_roots = _roots(archived["simple_roots"])
    if len(archived_roots) != 457 or archived.get("verified") is not True:
        raise ValueError("the bound-ten checkpoint is missing or unverified")
    if _digest(archived["simple_roots"]) != archived.get("roots_sha256"):
        raise ValueError("the bound-ten root digest does not match")
    cone_count = 10
    roots = list(archived_roots)
    nonzero_roots = list(archived_roots[cone_count:])
    start_after = Golden.from_json(archived["completed_through_control_coordinate"])
    witness = archived["hyperideal_witness"]
    levels: list[dict[str, object]] = []
    witness_history: list[dict[str, object]] = []

    if progress:
        print(
            "enumerating four levels through 4+5*phi",
            file=sys.stderr,
            flush=True,
        )
    extension = run_og10_golden_vinberg_extension(
        nonzero_roots,
        start_after=start_after,
        identity_bound=IDENTITY_BOUND,
        max_new_levels=4,
    )
    if extension.candidate_coordinates != EXPECTED_COORDINATES:
        raise AssertionError("unexpected next four admissible control coordinates")
    if progress:
        print("height enumeration complete; reconstructing rays", file=sys.stderr, flush=True)

    for level, expected_coordinate, expected_count in zip(
        extension.levels, EXPECTED_COORDINATES, EXPECTED_COUNTS, strict=True
    ):
        if len(level.roots) != expected_count:
            raise AssertionError(
                f"unexpected root count at {expected_coordinate}: {len(level.roots)}"
            )
        first_new_index = len(roots)
        roots.extend(level.roots)
        level_record = _level_payload(level)
        new_witness, signs = _witness_payload(roots, witness, first_new_index)
        level_record["previous_ray_new_wall_signs"] = dict(sorted(signs.items()))
        level_record["previous_ray_status"] = (
            "cut"
            if signs["positive"]
            else "survives_with_new_active_walls"
            if signs["zero"]
            else "survives_strictly"
        )
        if new_witness is None:
            raise ArithmeticError(
                f"no local hyperideal replacement found at {expected_coordinate}"
            )
        witness = new_witness
        witness_history.append(
            {
                "completed_through_control_coordinate": (
                    expected_coordinate.to_json()
                ),
                "facet_count": len(roots),
                "witness": witness,
            }
        )
        levels.append(level_record)
        if progress:
            print(
                f"completed {expected_coordinate}: {len(level.roots)} roots, "
                f"{level_record['visited_node_count']} nodes, "
                f"ray signs {dict(signs)}",
                file=sys.stderr,
                flush=True,
            )

    roots_tuple = tuple(roots)
    roots_json = [_vector(root) for root in roots_tuple]
    gram = og10_field_gram()
    functionals = tuple(_mat_vec(gram, root) for root in roots_tuple)
    relation_counts = Counter(
        _relation(_dot(left, functionals[right_index]))
        for left_index, left in enumerate(roots_tuple)
        for right_index in range(left_index + 1, len(roots_tuple))
    )
    final_witness = verify_hyperideal_extreme_ray(
        roots_tuple,
        gram,
        witness["active_facets"],
        dimension=DIMENSION,
    )
    final_witness_digest = _digest(
        [value.to_json() for value in final_witness.wall_pairings]
    )

    pell = folding._pell_seed()
    h = folding._compose(
        folding._root_reflection_matrix(
            roots_tuple[0], folding._covector(roots_tuple[0], gram)
        ),
        pell,
    )
    root_lookup = set(roots_tuple)
    transported_new = tuple(
        dict.fromkeys(
            image
            for root in roots_tuple
            for image in (folding._apply(h, root),)
            if image not in root_lookup
        )
    )
    transported_closed = (*roots_tuple, *transported_new)
    transported_lookup = set(transported_closed)
    transported_covectors = tuple(
        folding._covector(root, gram) for root in transported_closed
    )
    transported_pairwise_acute = all(
        folding._dot_covector(left, transported_covectors[j]).sign() <= 0
        for i, left in enumerate(transported_closed)
        for j in range(i + 1, len(transported_closed))
    )
    transported_levels = Counter(str(root[0]) for root in transported_new)
    transported_new_json = [_vector(root) for root in transported_new]
    transported_closed_json = [_vector(root) for root in transported_closed]

    checks = {
        "bound_ten_checkpoint_verified": archived["verified"] is True,
        "candidate_levels_are_exactly_the_next_four_after_bound_ten": tuple(
            Golden.from_json(level["control_coordinate"]) for level in levels
        )
        == EXPECTED_COORDINATES,
        "level_root_counts_are_71_13_371_232": tuple(
            level["root_count"] for level in levels
        )
        == EXPECTED_COUNTS,
        "facet_count_is_1144": len(roots_tuple) == 1144,
        "roots_are_distinct_integral_norm_two": (
            len(roots_tuple) == len(set(roots_tuple))
            and all(
                all(entry.is_integral for entry in root)
                and _pair(root, gram, root) == 2
                for root in roots_tuple
            )
        ),
        "all_pairs_are_Coxeter_or_ultraparallel": sum(relation_counts.values())
        == len(roots_tuple) * (len(roots_tuple) - 1) // 2,
        "every_intermediate_hyperideal_ray_independently_verified": all(
            record["witness"]["independently_verified"]
            for record in witness_history
        ),
        "final_hyperideal_ray_independently_verified": (
            final_witness.ray
            == tuple(Golden.from_json(entry) for entry in witness["ray"])
            and final_witness.norm == Golden.from_json(witness["norm"])
            and final_witness.norm.sign() > 0
            and final_witness_digest == witness["wall_pairings_sha256"]
        ),
        "h_transports_373_roots_beyond_the_1144_wall_prefix": (
            len(transported_new) == 373
        ),
        "h_transport_closure_has_1517_pairwise_acute_roots": (
            len(transported_closed) == 1517 and transported_pairwise_acute
        ),
        "h_transport_closure_is_involution_invariant": all(
            folding._apply(h, root) in transported_lookup
            for root in transported_closed
        ),
        "h_transport_level_counts_are_106_171_5_91": transported_levels
        == {
            "3 + 6*phi": 106,
            "4 + 6*phi": 171,
            "5 + 6*phi": 5,
            "4 + 7*phi": 91,
        },
    }
    return {
        "schema_version": "0.1.0",
        "certificate_kind": (
            "bounded_exact_OG10_norm_two_Vinberg_prefix_through_4_plus_5phi"
        ),
        "candidate_id": "og10_golden_h10_carrier",
        "hyperbolic_dimension": DIMENSION,
        "candidate_generation_identity_bound": IDENTITY_BOUND,
        "extension_start_after": archived[
            "completed_through_control_coordinate"
        ],
        "completed_through_control_coordinate": EXPECTED_COORDINATES[-1].to_json(),
        "next_admissible_control_coordinate": Golden(3, 6).to_json(),
        "candidate_control_coordinates": [
            value.to_json() for value in EXPECTED_COORDINATES
        ],
        "levels": levels,
        "extension_visited_node_count": sum(
            level["visited_node_count"] for level in levels
        ),
        "archived_facet_count": len(archived_roots),
        "new_facet_count": len(roots_tuple) - len(archived_roots),
        "facet_count": len(roots_tuple),
        "pair_relation_counts": dict(sorted(relation_counts.items())),
        "roots_sha256": _digest(roots_json),
        "simple_roots": roots_json,
        "hyperideal_witness_history": witness_history,
        "hyperideal_witness": witness,
        "complete_chamber_involution_transport": {
            "archived_root_count": len(roots_tuple),
            "images_retained_in_archived_prefix": (
                len(roots_tuple) - len(transported_new)
            ),
            "images_outside_archived_prefix": len(transported_new),
            "involution_closed_root_count": len(transported_closed),
            "new_control_level_counts": dict(
                sorted(transported_levels.items(), key=lambda item: item[0])
            ),
            "new_roots_sha256": _digest(transported_new_json),
            "closed_root_set_sha256": _digest(transported_closed_json),
            "all_1517_pairwise_acute": transported_pairwise_acute,
            "scope": (
                "the transported roots are genuine complete-chamber walls, "
                "but they do not exhaust any newly reached control level"
            ),
        },
        "checks": checks,
        "verified": all(checks.values()),
        "scope": (
            "the norm-two Vinberg search is exhaustive only through control "
            "coordinate 4+5*phi; candidate-generation bound 13 also contains "
            "the deliberately uncomputed next level 3+6*phi, which may cut "
            "the witness"
        ),
        "conclusion": (
            "the_1144_wall_through_4_plus_5phi_norm_two_outer_cone_has_an_exact_"
            "spacelike_extreme_ray_and_is_not_finite_volume;_this_is_a_"
            "bounded_prefix_theorem_only"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    payload = build(progress=args.progress)
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.summary_only:
        summary = {key: value for key, value in payload.items() if key != "simple_roots"}
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
