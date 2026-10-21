#!/usr/bin/env python3
"""Recompute the exact norm-two silver K3 Vinberg prefix through ``x=2``."""

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
    run_silver_k3_vinberg_prefix,
)


EXPECTED_CERTIFICATE_SHA256 = "43364d42b68527be4192f905e7c905ee255e4d88569eb1e25e7af9b7518a4dc6"
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/candidates/silver_k3_norm2_vinberg_prefix.json"


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
    prefix = run_silver_k3_vinberg_prefix()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_intrinsic_norm_two_vinberg_prefix",
        "claim": (
            "The Alexeev-aligned silver K3 carrier has an exact norm-two "
            "Vinberg prefix through 0<x<=2 with 9 height-zero and 3 "
            "positive-height walls. The prefix does not satisfy a Vinberg "
            "stopping criterion, so chamber completeness and volume remain "
            "undecidable."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "field_basis_gram": _jsonable(prefix.field_gram),
        "field_determinant": _jsonable(prefix.field_determinant),
        "field_identity_inertia": _jsonable(prefix.field_identity_inertia),
        "field_conjugate_inertia": _jsonable(prefix.field_conjugate_inertia),
        "control_vector": _jsonable(prefix.control_vector),
        "control_norm": _jsonable(prefix.control_norm),
        "control_covector": _jsonable(prefix.control_covector),
        "height_formula": "delta=(1+sqrt(2))*x^2/2",
        "height_zero": _jsonable(asdict(prefix.height_zero)),
        "cone_roots": _jsonable(prefix.cone_roots),
        "identity_bound": prefix.identity_bound,
        "candidate_coordinates": _jsonable(prefix.candidate_coordinates),
        "levels": [
            {
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
            }
            for level in prefix.levels
        ],
        "wall_count": len(prefix.roots),
        "roots_in_vinberg_order": _jsonable(prefix.roots),
        "wall_gram": _jsonable(prefix.wall_gram),
        "wall_identity_inertia": _jsonable(prefix.wall_identity_inertia),
        "wall_conjugate_inertia": _jsonable(prefix.wall_conjugate_inertia),
        "prefix_diagram_diagnostic": _jsonable(asdict(prefix.diagram)),
        "saturation_divisor_wall": _jsonable(asdict(prefix.saturation_wall)),
        "stopping_status": prefix.stopping_status,
        "scope": {
            "proved": [
                "All norm-two Vinberg candidates with 0<x<=2 were exhausted exactly.",
                "The height-zero norm-two system is exactly B4+F4+A1.",
                (
                    "The saturation (-2)-divisor of field norm 3-sqrt(2) "
                    "is not an integral O-reflection wall."
                ),
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
        raise ValueError("the supplied silver K3 prefix manifest does not reproduce")
    digest = manifest.get("certificate_sha256")
    if not isinstance(digest, str):
        raise ValueError("the silver K3 prefix digest is missing")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("the silver K3 prefix certificate is not frozen")
    return {
        "certificate_sha256": digest,
        "height_zero_root_count": manifest["height_zero"]["total_field_norm_two_roots"],
        "wall_count": manifest["wall_count"],
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
