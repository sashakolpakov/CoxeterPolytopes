#!/usr/bin/env python3
"""Verify the frozen exact finite-volume diagnostic for the 726-wall prefix."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTIC = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_726_global_diagnostic.json"
)
INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_4_plus_3sqrt2.json"
)
DIAGNOSE_SOURCE = ROOT / "scripts/diagnose_silver_k3_726_global_vinberg.py"
VERIFIER_SOURCE = ROOT / "src/coxeter_polytopes/verify/silver_vinberg_volume.py"
PINNED_FILES = {
    VERIFIER_SOURCE: "ced3bfcff9f74c50fd789d0df0b69519bf122b46313d539b84c2aa829703c2c3",
    DIAGNOSE_SOURCE: "1da7608b3d467f82c9ed19d30067c52742f4ed8ec09b5f4e09af4356816154b2",
    INPUT: "006fd85d81adc9aa6f88b272c56926e323b8038918798a019759e915a2910018",
}
EXPECTED_PREFIX_VISITS = (719, 6_429, 23_154, 46_687, 61_224, 57_746, 43_123, 25_764)
EXPECTED_FINITE_HISTOGRAM = {0: 2_332, 1: 14_428, 2: 9_004}


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"silver 726 diagnostic {name} is malformed")
    return value


def _integer_histogram(value: object, name: str) -> dict[int, int]:
    mapping = _mapping(value, name)
    return {int(key): int(count) for key, count in mapping.items()}


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "silver_k3_726_global_vinberg_prefix_diagnostic":
        raise ValueError("wrong silver 726 diagnostic kind")
    if document.get("schema_version") != "0.2.0-exploratory":
        raise ValueError("silver 726 diagnostic schema changed")
    for path, expected in PINNED_FILES.items():
        if _file_sha256(path) != expected:
            raise ValueError(f"pinned silver audit source changed: {path.relative_to(ROOT)}")

    frozen_digest = document.get("diagnostic_sha256")
    payload = dict(document)
    payload.pop("diagnostic_sha256", None)
    if not isinstance(frozen_digest, str) or _canonical_digest(payload) != frozen_digest:
        raise ValueError("silver 726 diagnostic payload changed")

    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    if (
        document.get("input_manifest") != str(INPUT.relative_to(ROOT))
        or document.get("input_certificate_sha256")
        != input_document.get("certificate_sha256")
        or input_document.get("wall_count") != 726
        or len(input_document.get("ordered_prefix_roots", ())) != 726
    ):
        raise ValueError("silver 726 input provenance changed")
    if document.get("dimension") != 9 or document.get("facet_count") != 726:
        raise ValueError("silver 726 ambient dimensions changed")
    if tuple(document.get("backtracking_positive_prefix_visits_by_order_1_through_8", ())) != EXPECTED_PREFIX_VISITS:
        raise ValueError("silver 726 positive backtracking transcript changed")

    transcript = _mapping(document.get("literal_transcript"), "literal transcript")
    raw_edges = transcript.get("elliptic_edges")
    raw_vertices = transcript.get("finite_vertices")
    raw_finite_rows = transcript.get("finite_endpoint_indices_by_edge")
    raw_ideal_rows = transcript.get("ideal_endpoint_indices_by_edge")
    if not all(isinstance(value, list) for value in (raw_edges, raw_vertices, raw_finite_rows, raw_ideal_rows)):
        raise ValueError("silver 726 literal transcript arrays are malformed")
    edges = tuple(tuple(map(int, row)) for row in raw_edges)  # type: ignore[arg-type]
    vertices = tuple(tuple(map(int, row)) for row in raw_vertices)  # type: ignore[arg-type]
    finite_rows = tuple(tuple(map(int, row)) for row in raw_finite_rows)  # type: ignore[arg-type]
    ideal_rows = tuple(tuple(map(int, row)) for row in raw_ideal_rows)  # type: ignore[arg-type]
    if (
        len(edges) != 25_764
        or len(vertices) != 3_604
        or len(finite_rows) != len(edges)
        or len(ideal_rows) != len(edges)
        or document.get("elliptic_rank_eight_edge_count") != len(edges)
        or document.get("elliptic_rank_nine_finite_vertex_count") != len(vertices)
    ):
        raise ValueError("silver 726 edge/vertex transcript sizes changed")
    if (
        len(set(edges)) != len(edges)
        or any(tuple(sorted(set(edge))) != edge or len(edge) != 8 for edge in edges)
        or len(set(vertices)) != len(vertices)
        or any(tuple(sorted(set(vertex))) != vertex or len(vertex) != 9 for vertex in vertices)
        or any(node < 0 or node >= 726 for row in (*edges, *vertices) for node in row)
    ):
        raise ValueError("silver 726 edge/vertex transcript is not normalized")

    vertex_set = set(vertices)
    for edge, endpoint_nodes in zip(edges, finite_rows, strict=True):
        if tuple(sorted(set(endpoint_nodes))) != endpoint_nodes:
            raise ValueError("silver 726 finite endpoint row is not normalized")
        for node in endpoint_nodes:
            if node in edge or tuple(sorted((*edge, node))) not in vertex_set:
                raise ValueError("silver 726 finite endpoint incidence changed")
    finite_histogram = dict(sorted(Counter(map(len, finite_rows)).items()))
    if (
        finite_histogram != EXPECTED_FINITE_HISTOGRAM
        or _integer_histogram(document.get("finite_endpoint_histogram"), "finite histogram")
        != finite_histogram
    ):
        raise ValueError("silver 726 finite endpoint histogram changed")

    raw_ideal_vertices = document.get("ideal_vertices")
    if not isinstance(raw_ideal_vertices, list):
        raise ValueError("silver 726 ideal vertices are malformed")
    ideal_vertices = tuple(
        (
            tuple(map(int, row["facets"])),
            tuple(tuple(map(int, component)) for component in row["affine_components"]),
        )
        for row in raw_ideal_vertices
    )
    for edge, endpoint_indices in zip(edges, ideal_rows, strict=True):
        if tuple(sorted(set(endpoint_indices))) != endpoint_indices:
            raise ValueError("silver 726 ideal endpoint row is not normalized")
        edge_set = set(edge)
        for endpoint_index in endpoint_indices:
            if endpoint_index < 0 or endpoint_index >= len(ideal_vertices):
                raise ValueError("silver 726 ideal endpoint index is out of range")
            facets, components = ideal_vertices[endpoint_index]
            if not edge_set.issubset(facets) or any(
                len(set(component) - edge_set) != 1 for component in components
            ):
                raise ValueError("silver 726 ideal endpoint incidence changed")
    ideal_histogram = dict(sorted(Counter(map(len, ideal_rows)).items()))
    total_histogram = dict(
        sorted(
            Counter(
                len(finite) + len(ideal)
                for finite, ideal in zip(finite_rows, ideal_rows, strict=True)
            ).items()
        )
    )
    if (
        document.get("ideal_vertex_count") != len(ideal_vertices)
        or len(ideal_vertices) != 0
        or ideal_histogram != {0: 25_764}
        or _integer_histogram(document.get("ideal_endpoint_histogram"), "ideal histogram")
        != ideal_histogram
        or _integer_histogram(document.get("finite_plus_ideal_endpoint_histogram"), "total histogram")
        != total_histogram
        or total_histogram != finite_histogram
    ):
        raise ValueError("silver 726 ideal/total endpoint transcript changed")

    affine = document.get("connected_affine_subdiagrams")
    neighborhoods = document.get("parabolic_neighborhood_audits")
    if (
        affine != []
        or neighborhoods != []
        or document.get("connected_affine_subdiagram_count") != 0
        or document.get("parabolic_neighborhood_audit_count") != 0
        or document.get("bad_parabolic_neighborhood_count") != 0
    ):
        raise ValueError("silver 726 affine/full-neighborhood census changed")
    if (
        document.get("invalid_finite_angle_pair_count") != 0
        or document.get("invalid_finite_angle_pair_sample") != []
    ):
        raise ValueError("silver 726 finite-angle validation changed")

    finite_volume_open = tuple(
        edge
        for edge, finite, ideal in zip(edges, finite_rows, ideal_rows, strict=True)
        if len(finite) + len(ideal) != 2
    )
    if (
        len(finite_volume_open) != 16_760
        or document.get("finite_volume_open_edge_count") != len(finite_volume_open)
        or document.get("open_edge_count") != len(finite_volume_open)
        or document.get("finite_volume_open_edge_sample")
        != [list(edge) for edge in finite_volume_open[:100]]
        or document.get("open_edge_sample")
        != [list(edge) for edge in finite_volume_open[:100]]
    ):
        raise ValueError("silver 726 open-edge transcript changed")
    if (
        document.get("compact_stopping_available") is not False
        or document.get("finite_volume_stopping_available") is not False
    ):
        raise ValueError("silver 726 negative stopping result changed")

    pair_histogram = _mapping(
        document.get("off_diagonal_pairing_histogram"), "pairing histogram"
    )
    if sum(map(int, pair_histogram.values())) != comb(726, 2):
        raise ValueError("silver 726 off-diagonal histogram does not sum")
    scope = _mapping(document.get("scope"), "scope")
    if scope.get("other_reflective_norms_exhausted") is not False:
        raise ValueError("silver 726 reflective-norm caveat disappeared")

    return {
        "facet_count": 726,
        "elliptic_edge_count": 25_764,
        "finite_vertex_count": 3_604,
        "ideal_vertex_count": 0,
        "finite_plus_ideal_endpoint_histogram": EXPECTED_FINITE_HISTOGRAM,
        "finite_volume_open_edge_count": 16_760,
        "connected_affine_subdiagram_count": 0,
        "full_neighborhood_audit_complete": True,
        "compact_stopping_available": False,
        "finite_volume_stopping_available": False,
        "final_chamber_infinite_volume_proved": False,
        "other_reflective_norms_exhausted": False,
    }


def fresh_replay() -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, str(DIAGNOSE_SOURCE), "--input", str(INPUT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    rebuilt = json.loads(completed.stdout)
    frozen = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    if rebuilt != frozen:
        raise AssertionError("fresh silver 726 diagnostic differs from archive")
    return {"fresh_source_replay_matches_archive": True}


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args(argv)
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    result = verify(document)
    if args.fresh:
        result["fresh_replay"] = fresh_replay()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
