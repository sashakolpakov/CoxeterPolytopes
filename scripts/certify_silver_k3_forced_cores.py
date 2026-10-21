#!/usr/bin/env python3
"""Freeze the silver K3 finite-core and Lanner-census theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.exact.silver import (
    SILVER_K3_TRACE_SCALAR,
    SQRT2,
    Silver,
)
from coxeter_polytopes.geometric.silver_finite import (
    silver_finite_trace_lift,
    silver_lanner_census,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/candidates/silver_k3_forced_cores.json"
SOURCE_FILES = (
    ROOT / "src/coxeter_polytopes/exact/silver.py",
    ROOT / "src/coxeter_polytopes/geometric/silver_finite.py",
    ROOT / "tests/test_silver_finite_trace_lifts.py",
)
EXPECTED_CERTIFICATE_SHA256 = "8f814d5db8a8073191c708879fe701cb489e69999122b66910e300995b2ee1e0"


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


def _silver(value: Silver) -> dict[str, str]:
    return value.to_json()


def _silver_matrix(matrix) -> list[list[dict[str, str]]]:
    return [[_silver(Silver.coerce(entry)) for entry in row] for row in matrix]


def _matrix_sha256(matrix, *, silver: bool = False) -> str:
    value = _silver_matrix(matrix) if silver else [list(row) for row in matrix]
    return _canonical_sha256(value)


def _inertia(value) -> list[int]:
    return [value.positive, value.negative, value.zero]


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def build() -> dict[str, object]:
    finite = []
    for field_type in ("I2(4)", "B3", "F4"):
        lift = silver_finite_trace_lift(field_type)
        if not lift.passes:
            raise AssertionError(f"the {field_type} silver trace lift failed")
        finite.append(
            {
                "field_type": field_type,
                "trace_type": lift.trace_type,
                "field_rank": lift.field_rank,
                "trace_rank": lift.trace_rank,
                "field_determinant": _silver(lift.field_determinant),
                "trace_determinant": lift.trace_determinant,
                "field_identity_inertia": _inertia(lift.field_identity_inertia),
                "field_conjugate_inertia": _inertia(lift.field_conjugate_inertia),
                "trace_inertia": _inertia(lift.trace_inertia),
                "isometry_determinant": lift.isometry_determinant,
                "isometry_verified": lift.isometry_verified,
                "field_gram_sha256": _matrix_sha256(lift.field_gram, silver=True),
                "trace_gram_sha256": _matrix_sha256(lift.trace_gram),
                "trace_isometry_sha256": _matrix_sha256(lift.trace_isometry),
            }
        )

    census = silver_lanner_census()
    if not census.passes:
        raise AssertionError("the silver Lanner census failed")
    lanners = []
    for block in census.blocks:
        lanners.append(
            {
                "order": block.order,
                "canonical_pair_code": list(block.canonical_label_key),
                "field_determinant": _silver(block.field_determinant),
                "field_identity_inertia": _inertia(block.field_identity_inertia),
                "field_conjugate_inertia": _inertia(block.field_conjugate_inertia),
                "proper_identity_inertias": [
                    _inertia(value) for value in block.proper_identity_inertias
                ],
                "trace_type": block.trace_type,
                "trace_determinant": block.trace_determinant,
                "trace_inertia": _inertia(block.trace_inertia),
                "trace_isometry_determinant": block.trace_isometry_determinant,
                "trace_isometry_verified": block.trace_isometry_verified,
                "curve_inertia": _inertia(block.curve_inertia),
                "fiber_square": block.fiber_square,
                "fiber_sqrt2_pairing": block.fiber_sqrt2_pairing,
                "sqrt2_fiber_square": block.sqrt2_fiber_square,
                "primitive_root_trace_u2_gram": [
                    list(row) for row in block.primitive_u2_gram
                ],
                "field_gram_sha256": _matrix_sha256(block.field_gram, silver=True),
                "trace_gram_sha256": _matrix_sha256(block.trace_gram),
                "trace_isometry_sha256": _matrix_sha256(block.trace_isometry),
            }
        )

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_forced_core_and_lanner_census",
        "claim": (
            "The norm-two Q(sqrt(2)) K3 track has finite trace lifts "
            "B3->D6(-1), F4->E8(-1), and exactly two arithmetic Lanner "
            "types in orders four/five: the one-4 cycles."
        ),
        "trace_scalar": _silver(SILVER_K3_TRACE_SCALAR),
        "trace_scalar_trace": str(SILVER_K3_TRACE_SCALAR.trace()),
        "trace_scalar_identity_sign": SILVER_K3_TRACE_SCALAR.sign(0),
        "trace_scalar_conjugate_sign": SILVER_K3_TRACE_SCALAR.sign(1),
        "codifferent_integral_witness": _silver(2 * SQRT2 * SILVER_K3_TRACE_SCALAR),
        "finite_trace_lifts": finite,
        "lanner_enumeration_stages": [
            {
                "order": stage.order,
                "graph_masks": stage.graph_masks,
                "graph_prefilter_survivors": stage.graph_prefilter_survivors,
                "decorated_candidates": stage.decorated_candidates,
                "labelled_exact_survivors": stage.labelled_exact_survivors,
                "isomorphism_types": stage.isomorphism_types,
            }
            for stage in census.stages
        ],
        "lanner_blocks": lanners,
        "source_files": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": _file_sha256(path),
            }
            for path in SOURCE_FILES
        ],
        "scope": {
            "proved": (
                "An exhaustive arithmetic seed classification and explicit "
                "integral trace-lattice isometries for the norm-two silver branch."
            ),
            "not_proved": [
                "Primitive embedding of a completed rank-20 trace lattice in the K3 lattice.",
                "Effectivity or common nef orientation of every root.",
                "Wall completeness, compactness, or a new H9 Coxeter polytope.",
                "Any exclusion for arbitrary diagrams mixing labels four and five.",
            ],
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("silver K3 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("silver K3 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("silver K3 certificate is not frozen")
    return {
        "finite_trace_types": [
            row["trace_type"] for row in manifest["finite_trace_lifts"]
        ],
        "lanner_orders": [row["order"] for row in manifest["lanner_blocks"]],
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
