"""Resume an archived exact H11 Vinberg prefix without recomputing it."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.h11_vinberg import run_h11_vinberg_extension


DEFAULT_PREFIX = (
    Path(__file__).parents[1]
    / "data"
    / "candidates"
    / "h11_vinberg_prefix_bound20.json"
)


def _golden(value: Golden) -> dict[str, str]:
    return value.to_json()


def _vector(value: tuple[Golden, ...]) -> list[dict[str, str]]:
    return [_golden(entry) for entry in value]


def _decode_vector(value: Any) -> tuple[Golden, ...]:
    if not isinstance(value, list):
        raise ValueError("a serialized H11 root must be an array")
    root = tuple(Golden.from_json(entry) for entry in value)
    if len(root) != 10:
        raise ValueError("a serialized H11 root must have ten coordinates")
    return root


def _level_record(level: Any) -> dict[str, Any]:
    return {
        "control_coordinate": _golden(level.control_coordinate),
        "height": _golden(level.height),
        "positive_definite_target": _golden(level.target),
        "coordinate_alphabet": {
            "size": len(level.alphabet.values),
            "identity_upper": level.alphabet.identity_upper,
            "conjugate_upper": level.alphabet.conjugate_upper,
            "a_bound": level.alphabet.a_bound,
            "b_bound": level.alphabet.b_bound,
        },
        "visited_by_depth": list(level.visited_by_depth),
        "parity_pruned": level.parity_pruned,
        "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
        "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
        "leaf_norm_mismatches": level.leaf_norm_mismatches,
        "solution_count": len(level.roots),
        "roots": [_vector(root) for root in level.roots],
    }


def extend_archived_prefix(
    prefix_path: Path,
    *,
    identity_bound: int,
    max_new_levels: int | None,
) -> dict[str, Any]:
    raw = prefix_path.read_bytes()
    repository = Path(__file__).parents[1].resolve()
    try:
        displayed_prefix_path = str(prefix_path.resolve().relative_to(repository))
    except ValueError:
        displayed_prefix_path = str(prefix_path)
    prefix = json.loads(raw)
    levels = prefix.get("levels")
    if not isinstance(levels, list) or not levels:
        raise ValueError("the base transcript has no completed positive-height levels")
    last = Golden.from_json(levels[-1]["control_coordinate"])
    previous = tuple(
        _decode_vector(root)
        for level in levels
        for root in level.get("roots", ())
    )
    declared_walls = prefix.get("output", {}).get("wall_count")
    height_zero_walls = len(prefix.get("height_zero", {}).get("simple_roots", ()))
    if declared_walls != height_zero_walls + len(previous):
        raise ValueError("the base transcript wall count is inconsistent")

    extension = run_h11_vinberg_extension(
        previous,
        start_after=last,
        identity_bound=identity_bound,
        max_new_levels=max_new_levels,
    )
    completed_levels = [_level_record(level) for level in extension.levels]
    new_wall_count = len(extension.new_roots)
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "exact_norm_two_vinberg_height_extension",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "lattice": "[[-2*phi,1],[1,2]] + E8",
        "base_prefix": {
            "path": displayed_prefix_path,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "identity_bound": prefix["control_coordinate_scan"]["identity_bound"],
            "last_processed_control_coordinate": _golden(last),
            "completed_level_count": len(levels),
            "positive_height_wall_count": len(previous),
            "total_wall_count": declared_walls,
        },
        "extension_scan": {
            "requested_identity_bound": identity_bound,
            "max_new_levels": max_new_levels,
            "candidates": [_golden(x) for x in extension.candidate_coordinates],
            "levels": completed_levels,
        },
        "output": {
            "completed_new_level_count": len(extension.levels),
            "new_wall_count": new_wall_count,
            "total_wall_count": declared_walls + new_wall_count,
            "visited_nodes": sum(
                sum(level.visited_by_depth) for level in extension.levels
            ),
        },
        "status": {
            "kind": "resumed_bounded_prefix",
            "last_processed_control_coordinate": (
                _golden(extension.candidate_coordinates[-1])
                if extension.candidate_coordinates
                else _golden(last)
            ),
            "completion_claim": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--identity-bound", type=int, default=21)
    parser.add_argument("--max-new-levels", type=int)
    parser.add_argument(
        "--summary",
        action="store_true",
        help="print provenance, counters, and roots but omit per-depth arrays",
    )
    args = parser.parse_args()
    transcript = extend_archived_prefix(
        args.prefix,
        identity_bound=args.identity_bound,
        max_new_levels=args.max_new_levels,
    )
    if args.summary:
        transcript["extension_scan"]["levels"] = [
            {
                "control_coordinate": level["control_coordinate"],
                "solution_count": level["solution_count"],
                "visited_nodes": sum(level["visited_by_depth"]),
                "roots": level["roots"],
            }
            for level in transcript["extension_scan"]["levels"]
        ]
    print(json.dumps(transcript, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
