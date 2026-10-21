#!/usr/bin/env python3
"""Freeze the exact silver C5 saturation and rank-20 K3 carrier theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.geometric.silver_k3_carrier import (
    silver_finite_completion_audits,
    silver_k3_carrier_theorem,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/candidates/silver_k3_carrier.json"
SOURCE_FILES = (
    ROOT / "src/coxeter_polytopes/geometric/silver_k3_carrier.py",
    ROOT / "tests/test_silver_k3_carrier.py",
)
EXPECTED_CERTIFICATE_SHA256 = "f26310dc2c9475653f8d64fdabe7ba15abe312953621f5f757fc7fb344dec0b9"


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_sha256(matrix: Sequence[Sequence[object]]) -> str:
    def encode(value: object) -> object:
        if isinstance(value, Silver):
            return value.to_json()
        if isinstance(value, Fraction):
            return str(value)
        return value

    return _canonical_sha256([[encode(value) for value in row] for row in matrix])


def _silver(value: Silver) -> dict[str, str]:
    return value.to_json()


def _inertia(value) -> list[int]:
    return [value.positive, value.negative, value.zero]


def _fraction(value: Fraction) -> str:
    return str(value)


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def build() -> dict[str, object]:
    theorem = silver_k3_carrier_theorem()
    if not theorem.passes:
        raise AssertionError("the silver K3 carrier theorem failed")
    seed = theorem.seed
    embedding = theorem.k3_embedding

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_c5_saturation_and_rank20_k3_carrier",
        "claim": (
            "The one-4 silver C5 Lanner seed has a unique nonzero even "
            "index-two trace overlattice stable under sqrt(2).  It is "
            "E8(-1)+U; adding B5 or F4+A1 gives two presentations of the "
            "(20,2,1) K3 Picard lattice, with an explicit primitive K3 "
            "embedding and complement <2>+<2>."
        ),
        "c5_saturation": {
            "raw_trace_type": "E8(-1) + U(2)",
            "raw_trace_determinant": seed.raw_trace_determinant,
            "raw_trace_inertia": _inertia(seed.raw_trace_inertia),
            "raw_smith_invariants": list(seed.raw_smith_invariants),
            "raw_trace_gram_sha256": _matrix_sha256(seed.raw_trace_gram),
            "raw_sqrt2_action_sha256": _matrix_sha256(seed.raw_sqrt2_action),
            "order_two_classes": [
                {
                    "support": [index for index, bit in enumerate(item.bits) if bit],
                    "square": _fraction(item.square),
                    "even": item.is_even,
                    "sqrt2_residue_support": [
                        index for index, bit in enumerate(item.sqrt2_residue) if bit
                    ],
                    "sqrt2_stable": item.sqrt2_stable,
                }
                for item in seed.order_two_classes
            ],
            "unique_stable_nonzero_even_support": [
                index for index, bit in enumerate(seed.saturation_bits) if bit
            ],
            "saturation_index": 2,
            "saturated_trace_type": seed.standard_trace_type,
            "saturated_trace_determinant": seed.saturated_trace_determinant,
            "saturated_trace_inertia": _inertia(seed.saturated_trace_inertia),
            "saturated_smith_invariants": list(seed.saturated_smith_invariants),
            "trace_isometry_determinant": seed.trace_isometry_determinant,
            "trace_isometry_verified": seed.trace_isometry_verified,
            "saturated_trace_gram_sha256": _matrix_sha256(
                seed.saturated_trace_gram
            ),
            "saturated_sqrt2_action_sha256": _matrix_sha256(
                seed.saturated_sqrt2_action
            ),
            "trace_isometry_sha256": _matrix_sha256(seed.trace_isometry),
            "sqrt2_action_integral": True,
            "sqrt2_action_self_adjoint": True,
            "sqrt2_action_square": 2,
            "extra_saturation_root": {
                "trace_coordinates": [str(value) for value in seed.saturation_vector],
                "trace_square": seed.extra_root_square,
                "trace_sqrt2_pairing": seed.extra_root_sqrt2_pairing,
                "sqrt2_root_coordinates": list(seed.sqrt2_extra_root),
                "sqrt2_root_square": seed.sqrt2_extra_root_square,
                "field_norm": _silver(seed.extra_root_field_norm),
                "field_norm_identity_sign": seed.extra_root_field_norm.sign(0),
                "field_norm_conjugate_sign": seed.extra_root_field_norm.sign(1),
            },
        },
        "finite_completions": [
            {
                "name": item.name,
                "field_rank": len(item.field_gram),
                "field_determinant": _silver(item.field_determinant),
                "field_identity_inertia": _inertia(item.field_identity_inertia),
                "field_conjugate_inertia": _inertia(item.field_conjugate_inertia),
                "trace_type": item.trace_type,
                "trace_determinant": item.trace_determinant,
                "trace_inertia": _inertia(item.trace_inertia),
                "trace_smith_invariants": list(item.trace_smith_invariants),
                "discriminant_class_squares_mod_2": [
                    _fraction(value) for value in item.discriminant_class_squares
                ],
                "coparity": item.coparity,
                "trace_isometry_determinant": item.trace_isometry_determinant,
                "trace_isometry_verified": item.trace_isometry_verified,
                "sqrt2_action_self_adjoint": item.sqrt2_action_self_adjoint,
                "sqrt2_action_square": 2,
                "field_gram_sha256": _matrix_sha256(item.field_gram),
                "trace_gram_sha256": _matrix_sha256(item.trace_gram),
                "trace_isometry_sha256": _matrix_sha256(item.trace_isometry),
            }
            for item in silver_finite_completion_audits()
        ],
        "carrier_models": [
            {
                "completion": item.name,
                "field_rank": len(item.field_gram),
                "field_determinant": _silver(item.field_determinant),
                "field_identity_inertia": _inertia(item.field_identity_inertia),
                "field_conjugate_inertia": _inertia(item.field_conjugate_inertia),
                "saturation_index": item.saturation_index,
                "trace_type": item.trace_type,
                "trace_determinant": item.trace_determinant,
                "trace_inertia": _inertia(item.trace_inertia),
                "trace_smith_invariants": list(item.trace_smith_invariants),
                "two_elementary_invariants": list(item.two_elementary_invariants),
                "trace_isometry_determinant": item.trace_isometry_determinant,
                "trace_isometry_verified": item.trace_isometry_verified,
                "sqrt2_action_integral": True,
                "sqrt2_action_self_adjoint": item.sqrt2_action_self_adjoint,
                "sqrt2_action_square": 2,
                "field_gram_sha256": _matrix_sha256(item.field_gram),
                "raw_trace_gram_sha256": _matrix_sha256(item.raw_trace_gram),
                "saturated_trace_gram_sha256": _matrix_sha256(
                    item.saturated_trace_gram
                ),
                "standard_trace_gram_sha256": _matrix_sha256(
                    item.standard_trace_gram
                ),
                "trace_isometry_sha256": _matrix_sha256(item.trace_isometry),
                "sqrt2_action_sha256": _matrix_sha256(item.sqrt2_action),
                "standard_sqrt2_action_sha256": _matrix_sha256(
                    item.standard_sqrt2_action
                ),
            }
            for item in theorem.models
        ],
        "k3_embedding": {
            "neron_severi_type": theorem.alexeev_lattice_type,
            "neron_severi_invariants": list(theorem.alexeev_invariants),
            "transcendental_gram": [list(row) for row in embedding.transcendental_gram],
            "discriminant_anti_isometry_squares": [
                [_fraction(left), _fraction(right)]
                for left, right in embedding.discriminant_anti_isometry_squares
            ],
            "gluing_index": embedding.gluing_index,
            "gluing_vectors": [
                [str(value) for value in vector] for vector in embedding.gluing_vectors
            ],
            "ambient_type": embedding.ambient_type,
            "ambient_determinant": embedding.ambient_determinant,
            "ambient_inertia": _inertia(embedding.ambient_inertia),
            "ambient_even": embedding.ambient_even,
            "ambient_unimodular": embedding.ambient_unimodular,
            "primitive_embedding_verified": embedding.primitive_embedding_verified,
            "orthogonal_complement_verified": embedding.orthogonal_complement_verified,
            "complement_gram": [list(row) for row in embedding.complement_gram],
            "ambient_gram_sha256": _matrix_sha256(embedding.ambient_gram),
            "embedding_matrix_sha256": _matrix_sha256(embedding.embedding_matrix),
            "complement_basis_sha256": _matrix_sha256(embedding.complement_basis),
        },
        "alexeev_identification": {
            "arxiv": "2209.09110v4",
            "ancillary_source_sha256": (
                "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
            ),
            "published_invariants": [20, 2, 1],
            "published_lattice_type": "U + E8(-1)^2 + A1(-1)^2",
            "classification_reason": (
                "An indefinite even 2-elementary lattice is determined by "
                "its signature and (r,a,delta); both carrier presentations "
                "have (20,2,1)."
            ),
        },
        "source_files": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": _file_sha256(path),
            }
            for path in SOURCE_FILES
        ],
        "scope": {
            "proved": (
                "A diagram-designed, sqrt(2)-stable rank-20 Picard-lattice "
                "carrier with an explicit primitive K3 embedding."
            ),
            "not_proved": [
                "Effectivity or a common nef orientation for the proposed roots.",
                "Closure under all (-2)-walls after saturation.",
                "An acute complete Coxeter chamber or a finite-volume H9 polytope.",
                "Compactness of any chamber cut out by these seed roots.",
            ],
            "wall_closure_required": theorem.wall_closure_required,
        },
    }

    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("silver K3 carrier certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("silver K3 carrier certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("silver K3 carrier certificate is not frozen")
    return {
        "seed_saturation": manifest["c5_saturation"]["saturated_trace_type"],
        "carrier_types": [row["trace_type"] for row in manifest["carrier_models"]],
        "k3_complement": manifest["k3_embedding"]["complement_gram"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--print-certificate", action="store_true")
    args = parser.parse_args(argv)
    if args.print_certificate:
        print(json.dumps(build(), indent=2, sort_keys=True))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
