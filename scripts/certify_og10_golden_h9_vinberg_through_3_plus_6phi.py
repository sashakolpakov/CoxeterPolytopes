#!/usr/bin/env python3
"""Freeze the direct H9 Vinberg continuation through ``3+6*phi``.

The source lattice is ``<-2*phi> + E8 + <2>``.  This certificate first
replays its Vinberg search from height zero and requires an exact ordered
match with all 61 fixed roots extracted from the independently archived
ambient OG10 prefix.  It then exhausts the unique remaining admissible
control coordinate below identity bound 13, namely ``3+6*phi``.

The existing 61-wall hyperideal ray is checked against every newly accepted
wall and then independently reverified in the enlarged cone.  No claim is
made about any later height.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_h9_vinberg import (
    _positive_definite_enumeration_data,
    admissible_h9_control_coordinates,
    h9_cone_roots,
    h9_field_gram,
    run_og10_golden_h9_vinberg_extension,
)
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray


ROOT = Path(__file__).resolve().parents[1]
AMBIENT_PREFIX = (
    ROOT / "data/candidates/og10_golden_vinberg_through_4_plus_5phi.json"
)
FIXED_SLICE = ROOT / "data/candidates/og10_golden_h9_fixed_slice.json"
IDENTITY_BOUND = 13
START_AFTER = Golden(4, 5)
COMPLETED_THROUGH = Golden(3, 6)
NEXT_ADMISSIBLE = Golden(4, 6)
EXPECTED_CANDIDATES = (
    Golden(1),
    PHI,
    Golden(1, 1),
    Golden(0, 2),
    Golden(1, 2),
    Golden(2, 2),
    Golden(1, 3),
    Golden(2, 3),
    Golden(3, 3),
    Golden(2, 4),
    Golden(3, 4),
    Golden(2, 5),
    Golden(3, 5),
    Golden(4, 5),
    Golden(3, 6),
)
EXPECTED_CALIBRATION_COUNTS = (2, 2, 2, 0, 4, 3, 2, 7, 0, 6, 7, 1, 13, 3)
EXPECTED_NEXT_ROOT_COUNT = 10
EXPECTED_NEXT_VISITED_NODE_COUNT = 156_899
EXPECTED_G_FORCED_INDICES = (1, 4, 5, 8)

Vector = tuple[Golden, ...]


def _encoded_vector(vector: Sequence[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode(
        "ascii"
    )
    return hashlib.sha256(encoded).hexdigest()


def _load_ambient_fixed_roots() -> tuple[dict[str, object], tuple[Vector, ...]]:
    payload = json.loads(AMBIENT_PREFIX.read_text())
    ambient = tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )
    fixed = tuple(
        (root[0], *root[2:])
        for root in ambient
        if not root[1]
    )
    return payload, fixed


def _level_payload(level) -> dict[str, object]:
    roots_json = [_encoded_vector(root) for root in level.roots]
    return {
        "control_coordinate": level.control_coordinate.to_json(),
        "height": level.height.to_json(),
        "target": level.target.to_json(),
        "root_count": len(level.roots),
        "roots_sha256": _digest(roots_json),
        "alphabet": {
            "cardinality": len(level.alphabet.values),
            "identity_upper": level.alphabet.identity_upper,
            "conjugate_upper": level.alphabet.conjugate_upper,
            "a_bound": level.alphabet.a_bound,
            "b_bound": level.alphabet.b_bound,
        },
        "visited_by_depth": list(level.visited_by_depth),
        "visited_node_count": sum(level.visited_by_depth),
        "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
        "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
        "leaf_norm_mismatches": level.leaf_norm_mismatches,
    }


def _apply_g(root: Sequence[Golden]) -> Vector:
    x, *middle, u = root
    return (PHI * x - u, *middle, PHI * x - PHI * u)


def build() -> dict[str, object]:
    ambient_payload, archived_roots = _load_ambient_fixed_roots()
    fixed_payload = json.loads(FIXED_SLICE.read_text())
    gram = h9_field_gram()
    cone = h9_cone_roots()

    direct = run_og10_golden_h9_vinberg_extension(
        (),
        start_after=Golden(0),
        identity_bound=IDENTITY_BOUND,
    )
    calibration_levels = direct.levels[:-1]
    next_level = direct.levels[-1]
    replayed_61 = (
        *cone,
        *(root for level in calibration_levels for root in level.roots),
    )
    all_roots = (*replayed_61, *next_level.roots)
    roots_json = [_encoded_vector(root) for root in all_roots]

    archived_digest = _digest([_encoded_vector(root) for root in archived_roots])
    fixed_digest = fixed_payload["known_wall_orbits"]["fixed_roots_sha256"]

    witness_payload = fixed_payload["known_61_wall_outer_cone"]
    old_ray = tuple(
        Golden.from_json(entry)
        for entry in witness_payload["hyperideal_witness_ray"]
    )
    old_active = tuple(
        witness_payload["hyperideal_witness_active_facets"]
    )
    new_pairings = tuple(
        _pair(old_ray, gram, root) for root in next_level.roots
    )
    witness = verify_hyperideal_extreme_ray(
        all_roots,
        gram,
        old_active,
        dimension=9,
    )

    g_images_outside_61 = tuple(
        dict.fromkeys(
            image
            for root in replayed_61
            for image in (_apply_g(root),)
            if image not in set(replayed_61)
        )
    )
    g_forced_indices = tuple(
        i for i, root in enumerate(next_level.roots) if root in g_images_outside_61
    )
    g_forced_source_indices = tuple(
        next(
            i for i, root in enumerate(replayed_61) if _apply_g(root) == forced
        )
        for forced in (next_level.roots[i] for i in g_forced_indices)
    )

    _, schur = _positive_definite_enumeration_data()
    conditional_diagonals = tuple(
        matrix[depth][depth] for depth, matrix in enumerate(schur)
    )
    pair_signs = Counter(
        "positive"
        if value.sign() > 0
        else "negative"
        if value.sign() < 0
        else "zero"
        for value in new_pairings
    )
    expected_next_pairings = (
        Golden(-22, -32) / 13,
        Golden(-42, -42) / 13,
        Golden(-10, 6) / 13,
        Golden(-22, -12) / 13,
        Golden(-32, -32) / 13,
        Golden(-32, -22) / 13,
        Golden(-16) / 13,
        Golden(-16, -6) / 13,
        Golden(-22, -22) / 13,
        Golden(-18, -12) / 13,
    )

    checks = {
        "ambient_checkpoint_is_verified_1144_wall_prefix": (
            ambient_payload.get("verified") is True
            and len(ambient_payload["simple_roots"]) == 1144
        ),
        "fixed_slice_checkpoint_is_verified": fixed_payload.get("verified") is True,
        "ambient_fixed_extraction_has_61_roots": len(archived_roots) == 61,
        "ambient_fixed_digest_matches_fixed_slice_checkpoint": (
            archived_digest == fixed_digest
        ),
        "candidate_list_through_identity_bound_13_is_exact": (
            direct.candidate_coordinates == EXPECTED_CANDIDATES
            and admissible_h9_control_coordinates(IDENTITY_BOUND)
            == EXPECTED_CANDIDATES
        ),
        "calibration_level_counts_match_61_wall_archive": (
            tuple(len(level.roots) for level in calibration_levels)
            == EXPECTED_CALIBRATION_COUNTS
        ),
        "direct_replay_reproduces_all_61_roots_in_identical_order": (
            replayed_61 == archived_roots
        ),
        "coordinate_alphabet_is_exhaustive_for_doubled_form": all(
            value >= 1 for value in conditional_diagonals
        ),
        "unique_level_after_4_plus_5phi_below_bound_13_is_3_plus_6phi": (
            tuple(
                value
                for value in direct.candidate_coordinates
                if (value - START_AFTER).sign() > 0
            )
            == (COMPLETED_THROUGH,)
        ),
        "next_level_has_10_roots_and_156899_nodes": (
            next_level.control_coordinate == COMPLETED_THROUGH
            and len(next_level.roots) == EXPECTED_NEXT_ROOT_COUNT
            and sum(next_level.visited_by_depth)
            == EXPECTED_NEXT_VISITED_NODE_COUNT
        ),
        "all_71_roots_are_distinct_integral_norm_two": (
            len(all_roots) == len(set(all_roots)) == 71
            and all(
                all(entry.is_integral for entry in root)
                and _pair(root, gram, root) == 2
                for root in all_roots
            )
        ),
        "all_71_roots_are_pairwise_acute": all(
            _pair(left, gram, right).sign() <= 0
            for i, left in enumerate(all_roots)
            for right in all_roots[i + 1 :]
        ),
        "all_10_new_walls_are_strict_on_old_hyperideal_ray": (
            pair_signs == {"negative": 10}
            and new_pairings == expected_next_pairings
        ),
        "unchanged_hyperideal_ray_is_independently_verified_for_71_walls": (
            witness.ray == old_ray
            and witness.norm
            == Golden.from_json(witness_payload["hyperideal_witness_norm"])
            and witness.norm.sign() > 0
            and witness.active_facets == old_active
        ),
        "four_next_level_roots_are_forced_by_g_transport": (
            len(g_images_outside_61) == 12
            and g_forced_indices == EXPECTED_G_FORCED_INDICES
            and g_forced_source_indices == (37, 51, 59, 60)
        ),
        "next_admissible_coordinate_is_4_plus_6phi": (
            next(
                value
                for value in admissible_h9_control_coordinates(14)
                if (value - COMPLETED_THROUGH).sign() > 0
            )
            == NEXT_ADMISSIBLE
        ),
    }

    return {
        "schema_version": "0.1.0",
        "certificate_kind": (
            "bounded_exact_direct_OG10_golden_H9_norm_two_Vinberg_prefix_"
            "through_3_plus_6phi"
        ),
        "candidate_id": "og10_golden_h9_fixed_lattice",
        "hyperbolic_dimension": 9,
        "fixed_lattice": {
            "decomposition": "<-2*phi> + E8 + <2>",
            "coordinate_order": ["x", "E8^8", "u"],
            "rank_over_Z_phi": 10,
        },
        "candidate_generation_identity_bound": IDENTITY_BOUND,
        "conjugate_control_coordinate_square_bound": PHI.to_json(),
        "completed_through_control_coordinate": COMPLETED_THROUGH.to_json(),
        "next_admissible_control_coordinate": NEXT_ADMISSIBLE.to_json(),
        "candidate_control_coordinates": [
            value.to_json() for value in direct.candidate_coordinates
        ],
        "positive_form_conditional_diagonal_coefficients": [
            str(value) for value in conditional_diagonals
        ],
        "calibration": {
            "method": (
                "from-scratch direct rank-10 field-lattice replay from the "
                "E8+A1 height-zero cone"
            ),
            "ambient_fixed_checkpoint_root_count": len(archived_roots),
            "ambient_fixed_checkpoint_sha256": archived_digest,
            "replayed_root_count": len(replayed_61),
            "replayed_roots_sha256": _digest(
                [_encoded_vector(root) for root in replayed_61]
            ),
            "ordered_root_match": replayed_61 == archived_roots,
            "levels": [_level_payload(level) for level in calibration_levels],
        },
        "new_level": {
            **_level_payload(next_level),
            "roots": [_encoded_vector(root) for root in next_level.roots],
            "old_hyperideal_ray_pairings": [
                value.to_json() for value in new_pairings
            ],
            "old_hyperideal_ray_pairing_signs": dict(sorted(pair_signs.items())),
        },
        "pell_fold_g_cross_check": {
            "g_row_major_x_u_block": [
                [PHI.to_json(), Golden(-1).to_json()],
                [PHI.to_json(), (-PHI).to_json()],
            ],
            "images_of_61_outside_61_count": len(g_images_outside_61),
            "next_level_forced_root_indices": list(g_forced_indices),
            "source_root_indices_in_61": list(g_forced_source_indices),
            "forced_roots": [
                _encoded_vector(next_level.roots[i]) for i in g_forced_indices
            ],
            "forced_root_old_ray_pairings": [
                new_pairings[i].to_json() for i in g_forced_indices
            ],
        },
        "facet_count": len(all_roots),
        "roots_sha256": _digest(roots_json),
        "simple_roots": roots_json,
        "hyperideal_witness": {
            "active_facets": list(witness.active_facets),
            "ray": _encoded_vector(witness.ray),
            "norm": witness.norm.to_json(),
            "wall_pairings_sha256": _digest(
                [value.to_json() for value in witness.wall_pairings]
            ),
            "unchanged_from_61_wall_checkpoint": True,
        },
        "checks": checks,
        "verified": all(checks.values()),
        "scope": (
            "The direct norm-two Vinberg search is exhaustive only through "
            "control coordinate 3+6*phi (equivalently through identity bound "
            "13).  Later simple roots may cut the displayed hyperideal ray."
        ),
        "conclusion": (
            "The exhaustive 3+6*phi level contributes ten walls.  All ten are "
            "strict on the previous spacelike extreme ray, so the 71-wall "
            "outer cone still has that exact hyperideal escape and is not of "
            "finite volume."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.summary_only:
        print(
            json.dumps(
                {key: value for key, value in payload.items() if key != "simple_roots"},
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
