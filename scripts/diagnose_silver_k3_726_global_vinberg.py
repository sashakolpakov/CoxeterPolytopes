#!/usr/bin/env python3
"""Run the exact global Vinberg prefix diagnostic on the 726 silver walls.

The input is the frozen ``x=4+3sqrt(2)`` manifest.  Unlike the shell search,
this script does not enumerate roots: it consumes the literal ordered wall
list, rebuilds its Gram matrix from the fixed rank-ten lattice, and applies
the sparse exact rank-eight edge/rank-nine finite endpoint census together
with the complete connected-affine/full-zero-neighbourhood ideal-endpoint
audit.

A failed prefix audit is only a search diagnostic.  A successful compact or
finite-volume audit can be joined to the exact height transcript by Vinberg
stopping; the finite-volume branch allows only ideal endpoints whose entire
zero-neighbourhood has the required componentwise affine rank.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Mapping, Sequence

from coxeter_polytopes.exact.silver import SILVER_ZERO, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    silver_k3_alexeev_field_gram,
)
from coxeter_polytopes.verify.silver_vinberg_volume import (
    SilverSparseVinbergAudit,
    audit_silver_vinberg_wall_gram,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_4_plus_3sqrt2.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_726_global_diagnostic.json"
)


def _note(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{stamp}] {message}", file=sys.stderr, flush=True)


def _decode_roots(document: Mapping[str, object]) -> tuple[tuple[Silver, ...], ...]:
    rows = document.get("ordered_prefix_roots")
    if not isinstance(rows, list):
        raise ValueError("manifest has no literal ordered_prefix_roots array")
    roots: list[tuple[Silver, ...]] = []
    for row_index, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != 10:
            raise ValueError(f"ordered_prefix_roots[{row_index}] is not rank ten")
        try:
            roots.append(
                tuple(
                    Silver.from_json(entry)
                    for entry in row
                    if isinstance(entry, Mapping)
                )
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid silver root at ordered_prefix_roots[{row_index}]"
            ) from exc
        if len(roots[-1]) != 10:
            raise ValueError(
                f"ordered_prefix_roots[{row_index}] has a malformed entry"
            )
    declared = document.get("wall_count")
    if not isinstance(declared, int) or len(roots) != declared:
        raise ValueError("literal root count does not match wall_count")
    return tuple(roots)


def _mat_vec(
    matrix: Sequence[Sequence[Silver]], vector: Sequence[Silver]
) -> tuple[Silver, ...]:
    return tuple(
        sum(
            (
                matrix[row][column] * vector[column]
                for column in range(len(vector))
                if matrix[row][column]
            ),
            SILVER_ZERO,
        )
        for row in range(len(matrix))
    )


def _dot(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    return sum(
        (left[index] * right[index] for index in range(len(left))),
        SILVER_ZERO,
    )


def build_wall_gram(
    roots: Sequence[Sequence[Silver]],
) -> tuple[tuple[Silver, ...], ...]:
    lattice = silver_k3_alexeev_field_gram()
    size = len(roots)
    work = [[SILVER_ZERO for _ in range(size)] for _ in range(size)]
    for row, root in enumerate(roots):
        image = _mat_vec(lattice, root)
        for column in range(row + 1):
            value = _dot(roots[column], image)
            work[row][column] = value
            work[column][row] = value
        if (row + 1) % 50 == 0 or row + 1 == size:
            _note(f"wall Gram rows {row + 1}/{size}")
    return tuple(tuple(row) for row in work)


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _audit_record(
    audit: SilverSparseVinbergAudit,
    *,
    input_path: Path,
    input_certificate: str,
) -> dict[str, object]:
    finite_endpoint_histogram = dict(audit.finite_endpoint_histogram)
    ideal_endpoint_histogram = dict(audit.ideal_endpoint_histogram)
    total_endpoint_histogram = dict(audit.finite_plus_ideal_endpoint_histogram)
    record: dict[str, object] = {
        "schema_version": "0.2.0-exploratory",
        "kind": "silver_k3_726_global_vinberg_prefix_diagnostic",
        "input_manifest": str(input_path.relative_to(ROOT)),
        "input_certificate_sha256": input_certificate,
        "dimension": audit.dimension,
        "facet_count": audit.facet_count,
        "backtracking_positive_prefix_visits_by_order_1_through_8": list(
            audit.backtracking_positive_prefix_visits
        ),
        "elliptic_rank_eight_edge_count": len(audit.elliptic_edges),
        "elliptic_rank_nine_finite_vertex_count": len(audit.finite_vertices),
        "finite_endpoint_histogram": {
            str(key): value for key, value in finite_endpoint_histogram.items()
        },
        "ideal_vertex_count": len(audit.ideal_vertices),
        "ideal_vertices": [
            {
                "facets": list(endpoint.facets),
                "affine_components": [
                    list(component) for component in endpoint.affine_components
                ],
            }
            for endpoint in audit.ideal_vertices
        ],
        "ideal_endpoint_histogram": {
            str(key): value for key, value in ideal_endpoint_histogram.items()
        },
        "finite_plus_ideal_endpoint_histogram": {
            str(key): value for key, value in total_endpoint_histogram.items()
        },
        "open_edge_count": len(audit.open_edges),
        "open_edge_sample": [list(edge) for edge in audit.open_edges[:100]],
        "finite_volume_open_edge_count": len(audit.finite_volume_open_edges),
        "finite_volume_open_edge_sample": [
            list(edge) for edge in audit.finite_volume_open_edges[:100]
        ],
        "connected_affine_subdiagram_count": len(
            audit.connected_affine_subdiagrams
        ),
        "connected_affine_subdiagrams": [
            list(vertices) for vertices in audit.connected_affine_subdiagrams
        ],
        "parabolic_neighborhood_audit_count": len(
            audit.parabolic_neighborhood_audits
        ),
        "bad_parabolic_neighborhood_count": len(
            audit.bad_parabolic_neighborhoods
        ),
        "parabolic_neighborhood_audits": [
            {
                "affine_component": list(row.affine_component),
                "full_zero_neighborhood": list(row.zero_neighborhood),
                "neighborhood_components": [
                    list(component) for component in row.neighborhood_components
                ],
                "is_rank_dimension_minus_one_parabolic": (
                    row.is_rank_dimension_minus_one_parabolic
                ),
            }
            for row in audit.parabolic_neighborhood_audits
        ],
        "invalid_finite_angle_pair_count": len(audit.invalid_finite_angle_pairs),
        "invalid_finite_angle_pair_sample": [
            [left, right, value.to_json()]
            for left, right, value in audit.invalid_finite_angle_pairs[:100]
        ],
        "compact_stopping_available": audit.compact_stopping_available,
        "finite_volume_stopping_available": (
            audit.finite_volume_stopping_available
        ),
        "literal_transcript": {
            "elliptic_edges": [list(edge) for edge in audit.elliptic_edges],
            "finite_vertices": [list(vertex) for vertex in audit.finite_vertices],
            "finite_endpoint_indices_by_edge": [
                list(row) for row in audit.finite_endpoint_indices_by_edge
            ],
            "ideal_endpoint_indices_by_edge": [
                list(row) for row in audit.ideal_endpoint_indices_by_edge
            ],
        },
        "scope": {
            "compact_success_meaning": (
                "When joined to the independently verified exact Vinberg height "
                "prefix, compact_stopping_available proves current-wall "
                "completeness and compactness."
            ),
            "finite_volume_success_meaning": (
                "When joined to the independently verified exact Vinberg height "
                "prefix, finite_volume_stopping_available proves current-wall "
                "completeness and finite covolume; ideal endpoints are counted "
                "only after Vinberg's full zero-neighbourhood test."
            ),
            "failure_meaning": (
                "Open edges or contaminated affine neighbourhoods in an "
                "incomplete prefix are undecided; later walls may cut off the "
                "corresponding escape."
            ),
            "other_reflective_norms_exhausted": False,
        },
    }
    record["diagnostic_sha256"] = _canonical_digest(record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-walls", type=int)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    document = json.loads(args.input.read_text(encoding="utf-8"))
    input_certificate = document.get("certificate_sha256")
    if not isinstance(input_certificate, str):
        raise ValueError("input certificate has no certificate_sha256")
    roots = _decode_roots(document)
    if args.max_walls is not None:
        if args.max_walls < 10:
            raise ValueError("--max-walls must be at least ten")
        roots = roots[: args.max_walls]
    _note(f"decoded {len(roots)} literal roots")
    wall_gram = build_wall_gram(roots)
    pairing_histogram = Counter(
        str(wall_gram[left][right])
        for left in range(len(wall_gram))
        for right in range(left)
    )
    _note(
        "wall Gram complete; distinct off-diagonal values="
        f"{len(pairing_histogram)}"
    )
    audit = audit_silver_vinberg_wall_gram(wall_gram, 9, progress=_note)
    record = _audit_record(
        audit,
        input_path=args.input,
        input_certificate=input_certificate,
    )
    record["off_diagonal_pairing_histogram"] = dict(
        sorted(pairing_histogram.items())
    )
    # Add the histogram before freezing the diagnostic digest.
    record.pop("diagnostic_sha256")
    record["diagnostic_sha256"] = _canonical_digest(record)
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
