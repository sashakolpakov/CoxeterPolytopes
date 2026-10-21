#!/usr/bin/env python3
"""Certify the silver Vinberg shell ``x=6+3sqrt(2)`` exactly."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_six_plus_three_sqrt2_extension,
)
from scripts.independent_silver_k3_vinberg_6_plus_3sqrt2 import (
    EXPECTED_REPLAY_SHA256,
    build_census,
)


PRIOR_CERTIFICATE_SHA256 = (
    "2caccab1462719095839350612a2278f52e01d66ac24ea25cd0c4748f5d10f8f"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "01c8c722b0db00f757d67789d82674f2e9724ec2e81e4b7e720870e1700a67f5"
)
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.json"
)


def _phase(message: str) -> None:
    print(f"[x=6+3sqrt(2) certificate] {message}", file=sys.stderr, flush=True)


def _jsonable(value: Any) -> Any:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if isinstance(value, Inertia):
        return {
            "positive": value.positive,
            "negative": value.negative,
            "zero": value.zero,
        }
    if isinstance(value, tuple):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, list):
        return [_jsonable(entry) for entry in value]
    if isinstance(value, Mapping):
        return {str(key): _jsonable(entry) for key, entry in value.items()}
    return value


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _vector_key(
    vector: Sequence[Silver],
) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((entry.a, entry.b) for entry in vector)


def _decode_roots(rows: object) -> tuple[tuple[Silver, ...], ...]:
    if not isinstance(rows, list):
        raise ValueError("the independent accepted-root transcript is malformed")
    roots = tuple(
        tuple(
            Silver.from_json(entry)
            for entry in row
            if isinstance(entry, Mapping)
        )
        for row in rows
        if isinstance(row, list)
    )
    if len(roots) != len(rows) or any(len(root) != 10 for root in roots):
        raise ValueError("an independent accepted root is malformed")
    return roots


def build_manifest() -> dict[str, object]:
    _phase("building the exact production extension")
    extension = run_silver_k3_vinberg_six_plus_three_sqrt2_extension()
    level = extension.new_level
    _phase(
        f"production extension complete; walls={len(extension.roots)}, "
        f"new roots={len(level.roots)}"
    )
    _phase("running the independent block-factorized replay")
    replay = build_census()
    _phase("independent replay complete; comparing exact root and blocker sets")

    independent_roots = _decode_roots(replay["accepted_roots"])
    production_roots_sorted = tuple(sorted(level.roots, key=_vector_key))
    production_rejections_sorted = tuple(
        sorted(
            level.previous_wall_rejection_witnesses,
            key=lambda item: (_vector_key(item[0]), item[1]),
        )
    )
    production_rejection_sha256 = _canonical_sha256(
        _jsonable(production_rejections_sorted)
    )
    reconstructed_sorted = tuple(
        sorted(
            (
                *(root for root, _blockers in production_rejections_sorted),
                *level.roots,
            ),
            key=_vector_key,
        )
    )
    accepted_set_matches = independent_roots == production_roots_sorted
    rejection_set_matches = (
        production_rejection_sha256 == replay["rejection_witnesses_sha256"]
    )
    reconstructed_set_matches = _canonical_sha256(
        _jsonable(reconstructed_sorted)
    ) == replay["all_reconstructed_roots_sha256"]
    if not (accepted_set_matches and rejection_set_matches and reconstructed_set_matches):
        raise AssertionError("production and independent x=6+3sqrt(2) sets differ")

    first_new = extension.prior_wall_count
    new_indices = tuple(range(first_new, len(extension.roots)))
    prior_pairings = tuple(
        tuple(extension.wall_gram[index][old] for old in range(first_new))
        for index in new_indices
    )
    same_height_pairings = tuple(
        extension.wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    blocker_size_histogram = tuple(
        sorted(
            Counter(
                len(blockers)
                for _root, blockers in level.previous_wall_rejection_witnesses
            ).items()
        )
    )

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_six_plus_three_shell",
        "claim": (
            "At x=6+3sqrt(2), exact two-embedding enumeration gives 697 "
            "weights; 425 reconstruct nonintegrally, 199 integral roots are "
            "blocked, and 73 roots enter the current norm-two prefix."
        ),
        "terminology": (
            "The 73 roots are simple only relative to the exhausted lower "
            "norm-two prefix and complete same-height cohort. The 1,291-wall "
            "list is not claimed to be a complete chamber."
        ),
        "prior_certificate_sha256": PRIOR_CERTIFICATE_SHA256,
        "ordered_prior_prefix": {
            "wall_count": first_new,
            "certificate_sha256": PRIOR_CERTIFICATE_SHA256,
            "reproducer": (
                "run_silver_k3_vinberg_four_plus_four_sqrt2_extension().roots"
            ),
            "ordered_roots_sha256": _canonical_sha256(
                _jsonable(extension.roots[:first_new])
            ),
        },
        "ordered_prefix_root_count": len(extension.roots),
        "ordered_prefix_root_order": (
            "1,218-root certified lower-height prefix followed by the "
            "73-root production weight-recursion order"
        ),
        "ordered_prefix_roots": _jsonable(extension.roots),
        "ordered_prefix_roots_sha256": _canonical_sha256(
            _jsonable(extension.roots)
        ),
        "completed_through_control_coordinate": _jsonable(
            extension.completed_through_control_coordinate
        ),
        "next_admissible_control_coordinate": _jsonable(
            extension.next_admissible_control_coordinate
        ),
        "candidate_coordinates": _jsonable(extension.candidate_coordinates),
        "new_level": {
            "control_coordinate": _jsonable(level.control_coordinate),
            "height": _jsonable(level.height),
            "target": _jsonable(level.target),
            "alphabet": _jsonable(asdict(level.alphabet)),
            "visited_by_depth": list(level.visited_by_depth),
            "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
            "exact_weight_solutions": level.exact_weight_solutions,
            "nonintegral_reconstructions": level.nonintegral_reconstructions,
            "integral_reconstructions": (
                level.exact_weight_solutions - level.nonintegral_reconstructions
            ),
            "previous_wall_rejections": level.previous_wall_rejections,
            "blocker_set_size_histogram": _jsonable(blocker_size_histogram),
            "production_rejection_witnesses": _jsonable(
                production_rejections_sorted
            ),
            "production_rejection_witness_set_sha256": (
                production_rejection_sha256
            ),
            "accepted_root_indices": list(new_indices),
            "accepted_root_count": len(level.roots),
            "accepted_root_order": "production_weight_recursion_order",
            "accepted_roots": _jsonable(level.roots),
            "accepted_roots_sha256": _canonical_sha256(_jsonable(level.roots)),
            "accepted_root_set_sha256": _canonical_sha256(
                _jsonable(production_roots_sorted)
            ),
            "pairings_to_prior_walls_sha256": _canonical_sha256(
                _jsonable(prior_pairings)
            ),
            "all_pairings_to_prior_walls_nonpositive": all(
                value.sign(0) <= 0 for row in prior_pairings for value in row
            ),
            "same_height_pairings_sha256": _canonical_sha256(
                _jsonable(same_height_pairings)
            ),
            "same_height_pair_count": len(same_height_pairings),
            "same_height_negative_pair_count": (
                extension.same_height_negative_pair_count
            ),
            "same_height_zero_pair_count": extension.same_height_zero_pair_count,
            "all_distinct_same_height_pairings_nonpositive": all(
                pairing.sign(0) <= 0 for pairing in same_height_pairings
            ),
        },
        "independent_block_factorized_replay": {
            key: replay[key]
            for key in (
                "block_partition",
                "component_box_counts",
                "component_bounded_counts",
                "norm_join_probe_count",
                "exact_weight_solutions",
                "nonintegral_reconstructions",
                "integral_reconstructions",
                "rejected_root_count",
                "blocker_set_size_histogram",
                "accepted_root_count",
                "same_height_negative_pair_count",
                "same_height_zero_pair_count",
                "prefix_wall_count",
                "next_admissible_control_coordinate",
                "all_reconstructed_roots_sha256",
                "rejection_witnesses_sha256",
                "accepted_roots_sha256",
                "census_sha256",
                "stopping_status",
            )
        },
        "independent_production_comparison": {
            "accepted_root_set_matches": accepted_set_matches,
            "rejection_witness_set_matches": rejection_set_matches,
            "all_reconstructed_root_set_matches": reconstructed_set_matches,
            "expected_replay_sha256": EXPECTED_REPLAY_SHA256,
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
        "global_diagram_audit": "not_run_for_the_incomplete_1291_wall_prefix",
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "Every norm-two Vinberg candidate at x=6+3sqrt(2) was exhausted exactly.",
                "The independent block replay reproduces all accepted roots and complete blocker sets.",
                "All 73 accepted roots enter simultaneously; all same-height pairings are negative.",
                "The complete ordered 1,291-root prefix and its exact two-embedding inertia are archived.",
                "The next admissible control coordinate is x=5+4sqrt(2).",
            ],
            "not_proved": [
                "Simple-root completeness or any global Vinberg stopping criterion.",
                "Compactness, finite covolume, or infinite covolume of the eventual chamber.",
                "Exhaustion of primitive reflective vectors of norms other than two.",
            ],
        },
    }
    manifest["certificate_sha256"] = _canonical_sha256(manifest)
    return manifest


def _summary(manifest: Mapping[str, object]) -> dict[str, object]:
    level = manifest["new_level"]
    replay = manifest["independent_block_factorized_replay"]
    if not all(isinstance(value, Mapping) for value in (level, replay)):
        raise ValueError("the x=6+3sqrt(2) manifest summary is malformed")
    return {
        "certificate_sha256": manifest["certificate_sha256"],
        "exact_weight_solutions": level["exact_weight_solutions"],
        "accepted_root_count": level["accepted_root_count"],
        "wall_count": manifest["wall_count"],
        "independent_replay_sha256": replay["census_sha256"],
        "wall_identity_inertia": manifest["wall_identity_inertia"],
        "wall_conjugate_inertia": manifest["wall_conjugate_inertia"],
        "next_admissible_control_coordinate": manifest[
            "next_admissible_control_coordinate"
        ],
        "stopping_status": manifest["stopping_status"],
    }


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    expected = build_manifest()
    if manifest != expected:
        raise ValueError(
            "the supplied x=6+3sqrt(2) silver K3 manifest does not reproduce"
        )
    digest = manifest.get("certificate_sha256")
    if (
        not isinstance(digest, str)
        or (
            EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
            and digest != EXPECTED_CERTIFICATE_SHA256
        )
    ):
        raise ValueError("the x=6+3sqrt(2) silver K3 certificate is not frozen")
    return _summary(manifest)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    recomputed = args.recompute or args.write
    manifest = (
        build_manifest()
        if recomputed
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    if args.write:
        ARCHIVE.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if recomputed:
        if (
            EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
            and manifest.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
        ):
            raise ValueError("the recomputed x=6+3sqrt(2) certificate changed")
        result = _summary(manifest)
    else:
        result = verify_manifest(manifest)
    print(json.dumps(manifest if args.full else result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
