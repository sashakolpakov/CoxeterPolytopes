#!/usr/bin/env python3
"""Recompute the silver ``x=2+sqrt(2)`` shell and targeted edge endpoint."""

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
    run_silver_k3_vinberg_edge_extension,
)


PRIOR_CERTIFICATE_SHA256 = (
    "f77614309c85b85366c829aacdfd8f03a028866b293ca5325b8646983f72cf81"
)
EXPECTED_CERTIFICATE_SHA256 = "9538f7d9845b0f31088280438be3843603f9dbbca7da2da450b9cf4c5e314b01"
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/candidates/silver_k3_norm2_vinberg_edge_extension.json"


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
    extension = run_silver_k3_vinberg_edge_extension()
    level = extension.new_level
    endpoint = extension.edge_endpoint
    pattern_codes = [
        "".join(
            "0" if value == Silver(0) else "3" if value == Silver(-1) else "4"
            for value in pattern
        )
        for pattern in endpoint.allowable_pairing_patterns
    ]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_targeted_edge_extension",
        "claim": (
            "The next admissible level x=2+sqrt(2) adds exactly 18 simple "
            "roots. Exactly one of them supplies a finite endpoint to the "
            "previously zero-endpoint A2+F4+B2 rank-eight edge."
        ),
        "prior_extension_certificate_sha256": PRIOR_CERTIFICATE_SHA256,
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
            "accepted_root_count": len(level.roots),
            "accepted_roots_sha256": _canonical_sha256(_jsonable(level.roots)),
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
        "target_edge_endpoint_audit": {
            "edge_indices": list(endpoint.edge_indices),
            "finite_type": endpoint.finite_type,
            "edge_gram_sha256": _canonical_sha256(_jsonable(endpoint.edge_gram)),
            "edge_determinant": _jsonable(endpoint.edge_determinant),
            "pairing_alphabet": _jsonable(endpoint.pairing_alphabet),
            "pairing_pattern_code_legend": {
                "0": "0",
                "3": "-1",
                "4": "-sqrt(2)",
            },
            "allowable_pairing_pattern_codes": pattern_codes,
            "new_endpoint_root_indices": list(endpoint.new_endpoint_root_indices),
            "endpoint_pairing": _jsonable(endpoint.endpoint_pairing),
            "endpoint_gram_determinant": _jsonable(
                endpoint.endpoint_gram_determinant
            ),
            "endpoint_identity_inertia": _jsonable(
                endpoint.endpoint_identity_inertia
            ),
            "endpoint_conjugate_inertia": _jsonable(
                endpoint.endpoint_conjugate_inertia
            ),
        },
        "endpoint_root": _jsonable(
            extension.roots[endpoint.new_endpoint_root_indices[0]]
        ),
        "global_diagram_audit": "not_run_for_the_34_wall_incomplete_prefix",
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "All norm-two Vinberg candidates through x=2+sqrt(2) were exhausted exactly.",
                "The A2+F4+B2 edge admits only 11 pairing patterns for a finite endpoint root.",
                "Exactly one new root in this shell extends that edge elliptically, with rank-nine determinant four.",
            ],
            "not_proved": [
                "A second finite or ideal endpoint for the targeted edge.",
                "Any global Vinberg stopping criterion for the 34-wall prefix.",
                "Completeness, compactness, or finite covolume of the eventual chamber.",
                "Exhaustion of primitive reflective vectors of norms other than two.",
            ],
        },
    }
    manifest["certificate_sha256"] = _canonical_sha256(manifest)
    return manifest


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    expected = build_manifest()
    if manifest != expected:
        raise ValueError("the supplied silver K3 edge-extension manifest does not reproduce")
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the silver K3 edge-extension digest is missing")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("the silver K3 edge-extension certificate is not frozen")
    return {
        "certificate_sha256": digest,
        "new_wall_count": manifest["new_level"]["accepted_root_count"],
        "wall_count": manifest["wall_count"],
        "allowable_edge_pairing_patterns": len(
            manifest["target_edge_endpoint_audit"]["allowable_pairing_pattern_codes"]
        ),
        "new_endpoint_root_indices": manifest["target_edge_endpoint_audit"][
            "new_endpoint_root_indices"
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
