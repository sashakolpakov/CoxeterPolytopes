"""Run Vinberg's exact volume gate on the complete high-rank K3 nef diagrams.

The rank-15 and rank-19 cases have hundreds of thousands of candidate edge
subdiagrams, so a complete rebuild is intentionally a standalone certificate
job rather than part of the ordinary high-rank Gram audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from coxeter_polytopes.verify.vinberg_volume import verify_vinberg_volume_criterion
from scripts.verify_k3_high_rank_nef_chambers import (
    DATA_PATH,
    build_gram,
    load_artifact,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_PATH = (
    ROOT / "data" / "known_examples" / "k3_high_rank_nef_vinberg_volume.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_audit_digest(result: Any) -> str:
    """Hash every exact subdiagram/extension output in a stable JSON order."""

    payload = {
        "finite_vertices": result.finite_vertices,
        "ideal_vertices": [
            {
                "facets": endpoint.facets,
                "affine_components": endpoint.affine_components,
            }
            for endpoint in result.ideal_vertices
        ],
        "edge_audits": [
            {
                "facets": audit.facets,
                "finite_endpoints": audit.finite_endpoints,
                "ideal_endpoints": audit.ideal_endpoints,
            }
            for audit in result.edge_audits
        ],
        "parabolic_neighborhood_audits": [
            {
                "affine_component": audit.affine_component,
                "zero_neighborhood": audit.zero_neighborhood,
                "neighborhood_components": audit.neighborhood_components,
                "is_rank_dimension_minus_one_parabolic": (
                    audit.is_rank_dimension_minus_one_parabolic
                ),
            }
            for audit in result.parabolic_neighborhood_audits
        ],
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def analyze_case(case: dict[str, Any]) -> dict[str, Any]:
    """Apply the complete-list gate to one published K3 curve diagram."""

    curve_gram = build_gram(case)
    wall_gram = [
        [-int(curve_gram[left, right]) for right in range(curve_gram.cols)]
        for left in range(curve_gram.rows)
    ]
    rho = int(case["rho"])
    result = verify_vinberg_volume_criterion(
        wall_gram,
        rho - 1,
        # Roulleau/Kondo's archived configuration is the complete set of
        # irreducible (-2)-curve walls, not a bounded root-search prefix.
        simple_roots_complete=True,
    )
    return {
        "rho": rho,
        "hyperbolic_dimension": rho - 1,
        "facet_count": curve_gram.rows,
        "simple_roots_complete": result.simple_roots_complete,
        "completeness_evidence": "published_complete_irreducible_curve_configuration",
        "decision": result.decision,
        "elliptic_edge_count": result.elliptic_edge_count,
        "finite_vertex_count": result.finite_vertex_count,
        "ideal_vertex_count": result.ideal_vertex_count,
        "connected_affine_audit_count": len(result.parabolic_neighborhood_audits),
        "bad_full_parabolic_neighborhood_count": sum(
            not audit.is_rank_dimension_minus_one_parabolic
            for audit in result.parabolic_neighborhood_audits
        ),
        "finite_volume_bad_edge_count": len(result.finite_volume_bad_edges),
        "canonical_full_audit_sha256": _canonical_audit_digest(result),
        "all_checks_pass": (
            result.finite_volume is True
            and result.compact is False
            and not result.finite_volume_bad_edges
            and all(
                audit.is_rank_dimension_minus_one_parabolic
                for audit in result.parabolic_neighborhood_audits
            )
        ),
    }


def build_certificate(rhos: set[int] | None = None) -> dict[str, Any]:
    artifact = load_artifact()
    selected = [
        case
        for case in artifact["cases"]
        if rhos is None or int(case["rho"]) in rhos
    ]
    return {
        "schema_version": 1,
        "artifact_id": "k3_high_rank_nef_vinberg_volume",
        "source_artifact": str(DATA_PATH.relative_to(ROOT)),
        "source_sha256": _sha256(DATA_PATH),
        "criterion": (
            "exact_elliptic_edge_endpoint_and_full_parabolic_zero_neighborhood"
        ),
        "audit_digest_schema": (
            "canonical_json_v1_finite_vertices_ideal_components_edges_parabolic_audits"
        ),
        "checker": "coxeter_polytopes.verify.vinberg_volume",
        "cases": [analyze_case(case) for case in selected],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rho",
        type=int,
        action="append",
        help="recompute only this Picard rank (repeatable)",
    )
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="compare a full recomputation with the archived certificate",
    )
    arguments = parser.parse_args()
    certificate = build_certificate(set(arguments.rho) if arguments.rho else None)
    if arguments.check_archive:
        if arguments.rho:
            raise SystemExit("--check-archive requires a full five-case rebuild")
        archived = json.loads(ARCHIVE_PATH.read_text())
        if certificate != archived:
            raise SystemExit("recomputed certificate differs from archive")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
