#!/usr/bin/env python3
"""Certify the empty silver Vinberg shell ``x=1+2sqrt(2)`` exactly."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    run_silver_k3_vinberg_one_plus_two_sqrt2_extension,
)


PRIOR_CERTIFICATE_SHA256 = (
    "9538f7d9845b0f31088280438be3843603f9dbbca7da2da450b9cf4c5e314b01"
)
EXPECTED_CERTIFICATE_SHA256 = "e644e073dacd4c691821409920c723f8c2edebc2ee7c9d751cb4b2059d0950c2"
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/candidates/silver_k3_norm2_vinberg_extension_1_plus_2sqrt2.json"
)


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


def build_manifest() -> dict[str, object]:
    extension = run_silver_k3_vinberg_one_plus_two_sqrt2_extension()
    level = extension.new_level
    rejection_witnesses = [
        {
            "root": _jsonable(root),
            "blocking_wall_indices": list(blockers),
        }
        for root, blockers in level.previous_wall_rejection_witnesses
    ]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_empty_shell",
        "claim": (
            "The next admissible control level x=1+2sqrt(2) has exactly four "
            "integral norm-equation solutions. Every solution is cut off by "
            "an earlier simple wall, so this shell adds no simple root and "
            "does not supply the targeted edge's missing endpoint."
        ),
        "prior_edge_extension_certificate_sha256": PRIOR_CERTIFICATE_SHA256,
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
            "previous_wall_rejections": level.previous_wall_rejections,
            "previous_wall_rejection_witnesses": rejection_witnesses,
            "accepted_root_count": len(level.roots),
            "accepted_roots_sha256": _canonical_sha256(_jsonable(level.roots)),
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
        "target_edge_audit": {
            "edge_indices": list(extension.edge_indices),
            "known_finite_endpoint_root_indices": list(
                extension.finite_endpoint_root_indices
            ),
            "new_finite_endpoint_root_indices": list(
                extension.new_finite_endpoint_root_indices
            ),
            "ideal_endpoint_audit_status": extension.ideal_endpoint_audit_status,
            "second_endpoint_status": extension.second_endpoint_status,
        },
        "endpoint_completion_criterion": {
            "finite": (
                "A second accepted simple root r for which the principal "
                "Gram matrix on E union {r} is positive definite of rank nine."
            ),
            "ideal": (
                "A full componentwise-affine zero-neighbourhood of total "
                "parabolic rank eight incident to E, checked only against "
                "the complete simple-root list."
            ),
            "terminal": (
                "After simple-root completeness, E must have exactly two "
                "distinct finite-or-ideal endpoints; every other elliptic "
                "rank-eight edge and every parabolic neighbourhood must pass "
                "the same literal Vinberg audit."
            ),
        },
        "global_diagram_audit": "not_run_for_the_incomplete_34_wall_prefix",
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "Every norm-two Vinberg candidate at x=1+2sqrt(2) was exhausted exactly.",
                "All four integral norm-equation solutions are blocked by previously accepted walls.",
                "No new simple wall, hence no new finite or ideal endpoint data, arises at this shell.",
                "The next admissible control coordinate is x=3+sqrt(2).",
            ],
            "not_proved": [
                "A second finite or ideal endpoint for the targeted edge.",
                "A full parabolic-neighbourhood audit for the incomplete wall prefix.",
                "Any global Vinberg stopping criterion or simple-root completeness.",
                "Compactness, finite covolume, or infinite covolume of the eventual chamber.",
                "Exhaustion of primitive reflective vectors of norms other than two.",
            ],
        },
    }
    manifest["certificate_sha256"] = _canonical_sha256(manifest)
    return manifest


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    expected = build_manifest()
    if manifest != expected:
        raise ValueError(
            "the supplied x=1+2sqrt(2) silver K3 manifest does not reproduce"
        )
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the x=1+2sqrt(2) silver K3 digest is missing")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("the x=1+2sqrt(2) silver K3 certificate is not frozen")
    new_level = manifest["new_level"]
    target = manifest["target_edge_audit"]
    if not isinstance(new_level, Mapping) or not isinstance(target, Mapping):
        raise ValueError("the x=1+2sqrt(2) summary is malformed")
    return {
        "certificate_sha256": digest,
        "exact_weight_solutions": new_level["exact_weight_solutions"],
        "accepted_root_count": new_level["accepted_root_count"],
        "wall_count": manifest["wall_count"],
        "known_finite_endpoint_root_indices": target[
            "known_finite_endpoint_root_indices"
        ],
        "next_admissible_control_coordinate": manifest[
            "next_admissible_control_coordinate"
        ],
        "stopping_status": manifest["stopping_status"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    manifest = (
        build_manifest()
        if args.recompute
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    result = verify_manifest(manifest)
    if args.full:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
