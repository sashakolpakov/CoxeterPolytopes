#!/usr/bin/env python3
"""Recompute the first silver K3 Vinberg shell beyond ``x=2``."""

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
    run_silver_k3_vinberg_extension,
)


PRIOR_CERTIFICATE_SHA256 = (
    "43364d42b68527be4192f905e7c905ee255e4d88569eb1e25e7af9b7518a4dc6"
)
EXPECTED_CERTIFICATE_SHA256 = "f77614309c85b85366c829aacdfd8f03a028866b293ca5325b8646983f72cf81"
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_1_plus_sqrt2.json"


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
    extension = run_silver_k3_vinberg_extension()
    level = extension.new_level
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_first_extension",
        "claim": (
            "The first admissible norm-two Vinberg level beyond x=2 on the "
            "Alexeev-aligned silver K3 carrier is x=1+sqrt(2). It adds exactly "
            "four walls. The resulting 16-wall prefix does not satisfy a "
            "Vinberg stopping criterion."
        ),
        "prior_bound_two_certificate_sha256": PRIOR_CERTIFICATE_SHA256,
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "identity_bound": extension.identity_bound,
        "candidate_coordinates": _jsonable(extension.candidate_coordinates),
        "completed_through_control_coordinate": _jsonable(
            extension.completed_through_control_coordinate
        ),
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
            "accepted_roots": _jsonable(level.roots),
        },
        "prior_wall_count": extension.prior_wall_count,
        "wall_count": len(extension.roots),
        "wall_gram_sha256": _canonical_sha256(_jsonable(extension.wall_gram)),
        "wall_identity_inertia": _jsonable(extension.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(extension.wall_conjugate_inertia),
        "prefix_diagram_diagnostic": _jsonable(asdict(extension.diagram)),
        "diagram_pruning_lemma": (
            "For norm-two acute walls, every PSD principal submatrix has all "
            "pairings >=-2 and every PD principal submatrix has all pairings >-2; "
            "exact inertia is still checked after this necessary filter."
        ),
        "stopping_status": extension.stopping_status,
        "scope": {
            "proved": [
                "Every admissible norm-two Vinberg level through x<=3 was exhausted exactly.",
                "The only new control coordinate after x=2 is x=1+sqrt(2).",
                "That level has exactly four accepted simple roots.",
                "The resulting 16-wall prefix remains acute and has no connected affine subdiagram.",
            ],
            "not_proved": [
                "The complete norm-two simple-root set.",
                "Compactness or finite covolume of the eventual chamber.",
                "Exhaustion of primitive reflective vectors of norms other than two.",
            ],
        },
    }
    manifest["certificate_sha256"] = _canonical_sha256(manifest)
    return manifest


def verify_manifest(manifest: Mapping[str, object]) -> dict[str, object]:
    expected = build_manifest()
    if manifest != expected:
        raise ValueError("the supplied silver K3 extension manifest does not reproduce")
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the silver K3 extension digest is missing")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("the silver K3 extension certificate is not frozen")
    return {
        "certificate_sha256": digest,
        "new_wall_count": len(manifest["new_level"]["accepted_roots"]),
        "wall_count": manifest["wall_count"],
        "open_edge_count": manifest["prefix_diagram_diagnostic"]["open_edge_count"],
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
