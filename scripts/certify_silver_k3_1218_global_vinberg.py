#!/usr/bin/env python3
"""Verify the frozen exact global diagnostic for the 1,218-wall prefix."""

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
    ROOT / "data/candidates/silver_k3_norm2_vinberg_1218_global_diagnostic.json"
)
INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_4_plus_4sqrt2.json"
)
DIAGNOSE_SOURCE = ROOT / "scripts/diagnose_silver_k3_1218_global_vinberg.py"
DIAGNOSE_HELPER_SOURCE = ROOT / "scripts/diagnose_silver_k3_726_global_vinberg.py"
VERIFIER_SOURCE = ROOT / "src/coxeter_polytopes/verify/silver_vinberg_volume.py"
PINNED_FILES = {
    VERIFIER_SOURCE: "ced3bfcff9f74c50fd789d0df0b69519bf122b46313d539b84c2aa829703c2c3",
    DIAGNOSE_SOURCE: "5f68318d909367036d3b13996b8bc4756efa6b6dd973f090a8950be5f51f5480",
    DIAGNOSE_HELPER_SOURCE: (
        "1da7608b3d467f82c9ed19d30067c52742f4ed8ec09b5f4e09af4356816154b2"
    ),
    INPUT: "4d6d23e8baa75251aefff796bc8600a485823db3ddf4a1a8c0024f74340e7ab8",
}
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "2caccab1462719095839350612a2278f52e01d66ac24ea25cd0c4748f5d10f8f"
)
EXPECTED_DIAGNOSTIC_SHA256 = (
    "cfb9272463b66ba8295a57b79a80bea4a311f64b6873a987e3098d4a5f1b0dea"
)
EXPECTED_PAIRING_HISTOGRAM_SHA256 = (
    "29989ed7befef001484c5d68731b208fbac47538d0e801d48760200e0c7ca99d"
)
EXPECTED_PREFIX_VISITS = (
    1211,
    11054,
    39647,
    78884,
    101295,
    94467,
    68937,
    42297,
)
EXPECTED_AFFINE_GROWTH = (
    731828,
    17502,
    30365,
    34816,
    25215,
    15458,
    7310,
    1949,
    6,
    0,
)
EXPECTED_FINITE_HISTOGRAM = {0: 3865, 1: 23989, 2: 14443}


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"silver 1218 diagnostic {name} is malformed")
    return value


def _integer_histogram(value: object, name: str) -> dict[int, int]:
    mapping = _mapping(value, name)
    return {int(key): int(count) for key, count in mapping.items()}


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "silver_k3_1218_global_vinberg_prefix_diagnostic":
        raise ValueError("wrong silver 1218 diagnostic kind")
    if document.get("schema_version") != "1.0.0":
        raise ValueError("silver 1218 diagnostic schema changed")
    for path, expected in PINNED_FILES.items():
        if _file_sha256(path) != expected:
            raise ValueError(
                f"pinned silver 1218 audit source changed: {path.relative_to(ROOT)}"
            )

    frozen_digest = document.get("diagnostic_sha256")
    payload = dict(document)
    payload.pop("diagnostic_sha256", None)
    if (
        frozen_digest != EXPECTED_DIAGNOSTIC_SHA256
        or _canonical_digest(payload) != frozen_digest
    ):
        raise ValueError("silver 1218 diagnostic payload changed")

    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    if (
        document.get("input_manifest") != str(INPUT.relative_to(ROOT))
        or document.get("input_certificate_sha256")
        != EXPECTED_INPUT_CERTIFICATE_SHA256
        or input_document.get("certificate_sha256")
        != EXPECTED_INPUT_CERTIFICATE_SHA256
        or input_document.get("wall_count") != 1218
        or len(input_document.get("ordered_prefix_roots", ())) != 1218
    ):
        raise ValueError("silver 1218 input provenance changed")
    if document.get("dimension") != 9 or document.get("facet_count") != 1218:
        raise ValueError("silver 1218 ambient dimensions changed")
    if tuple(
        document.get(
            "backtracking_positive_prefix_visits_by_order_1_through_8", ()
        )
    ) != EXPECTED_PREFIX_VISITS:
        raise ValueError("silver 1218 production-order backtracking changed")

    transcript = _mapping(document.get("literal_transcript"), "literal transcript")
    raw_edges = transcript.get("elliptic_edges")
    raw_vertices = transcript.get("finite_vertices")
    raw_finite_rows = transcript.get("finite_endpoint_indices_by_edge")
    raw_ideal_rows = transcript.get("ideal_endpoint_indices_by_edge")
    if not all(
        isinstance(value, list)
        for value in (raw_edges, raw_vertices, raw_finite_rows, raw_ideal_rows)
    ):
        raise ValueError("silver 1218 literal transcript arrays are malformed")
    edges = tuple(tuple(map(int, row)) for row in raw_edges)  # type: ignore[arg-type]
    vertices = tuple(tuple(map(int, row)) for row in raw_vertices)  # type: ignore[arg-type]
    finite_rows = tuple(tuple(map(int, row)) for row in raw_finite_rows)  # type: ignore[arg-type]
    ideal_rows = tuple(tuple(map(int, row)) for row in raw_ideal_rows)  # type: ignore[arg-type]
    if (
        len(edges) != 42297
        or len(vertices) != 5875
        or len(finite_rows) != len(edges)
        or len(ideal_rows) != len(edges)
        or document.get("elliptic_rank_eight_edge_count") != len(edges)
        or document.get("elliptic_rank_nine_finite_vertex_count") != len(vertices)
    ):
        raise ValueError("silver 1218 edge/vertex transcript sizes changed")
    if (
        len(set(edges)) != len(edges)
        or any(tuple(sorted(set(edge))) != edge or len(edge) != 8 for edge in edges)
        or len(set(vertices)) != len(vertices)
        or any(
            tuple(sorted(set(vertex))) != vertex or len(vertex) != 9
            for vertex in vertices
        )
        or any(node < 0 or node >= 1218 for row in (*edges, *vertices) for node in row)
    ):
        raise ValueError("silver 1218 edge/vertex transcript is not normalized")

    vertex_set = set(vertices)
    for edge, endpoint_nodes in zip(edges, finite_rows, strict=True):
        if tuple(sorted(set(endpoint_nodes))) != endpoint_nodes:
            raise ValueError("silver 1218 finite endpoint row is not normalized")
        for node in endpoint_nodes:
            if node in edge or tuple(sorted((*edge, node))) not in vertex_set:
                raise ValueError("silver 1218 finite endpoint incidence changed")
    finite_histogram = dict(sorted(Counter(map(len, finite_rows)).items()))
    if (
        finite_histogram != EXPECTED_FINITE_HISTOGRAM
        or _integer_histogram(
            document.get("finite_endpoint_histogram"), "finite histogram"
        )
        != finite_histogram
    ):
        raise ValueError("silver 1218 finite endpoint histogram changed")

    raw_ideal_vertices = document.get("ideal_vertices")
    if not isinstance(raw_ideal_vertices, list):
        raise ValueError("silver 1218 ideal vertices are malformed")
    if (
        raw_ideal_vertices != []
        or document.get("ideal_vertex_count") != 0
        or any(ideal_rows)
        or _integer_histogram(
            document.get("ideal_endpoint_histogram"), "ideal histogram"
        )
        != {0: 42297}
        or _integer_histogram(
            document.get("finite_plus_ideal_endpoint_histogram"),
            "combined endpoint histogram",
        )
        != finite_histogram
    ):
        raise ValueError("silver 1218 ideal/combined endpoint census changed")

    affine = document.get("connected_affine_subdiagrams")
    neighborhoods = document.get("parabolic_neighborhood_audits")
    if (
        tuple(document.get("connected_affine_growth_next_by_order", ()))
        != EXPECTED_AFFINE_GROWTH
        or affine != []
        or neighborhoods != []
        or document.get("connected_affine_subdiagram_count") != 0
        or document.get("parabolic_neighborhood_audit_count") != 0
        or document.get("bad_parabolic_neighborhood_count") != 0
    ):
        raise ValueError("silver 1218 affine/full-neighborhood census changed")
    if (
        document.get("invalid_finite_angle_pair_count") != 0
        or document.get("invalid_finite_angle_pair_sample") != []
    ):
        raise ValueError("silver 1218 finite-angle validation changed")

    finite_volume_open = tuple(
        edge
        for edge, finite, ideal in zip(edges, finite_rows, ideal_rows, strict=True)
        if len(finite) + len(ideal) != 2
    )
    if (
        len(finite_volume_open) != 27854
        or document.get("finite_volume_open_edge_count") != len(finite_volume_open)
        or document.get("open_edge_count") != len(finite_volume_open)
        or document.get("finite_volume_open_edge_sample")
        != [list(edge) for edge in finite_volume_open[:100]]
        or document.get("open_edge_sample")
        != [list(edge) for edge in finite_volume_open[:100]]
    ):
        raise ValueError("silver 1218 open-edge transcript changed")
    if (
        document.get("compact_stopping_available") is not False
        or document.get("finite_volume_stopping_available") is not False
    ):
        raise ValueError("silver 1218 negative stopping result changed")

    pair_histogram = dict(
        _mapping(
            document.get("off_diagonal_pairing_histogram"),
            "pairing histogram",
        )
    )
    if (
        len(pair_histogram) != 173
        or sum(map(int, pair_histogram.values())) != comb(1218, 2)
        or _canonical_digest(pair_histogram)
        != EXPECTED_PAIRING_HISTOGRAM_SHA256
        or document.get("off_diagonal_pairing_histogram_sha256")
        != EXPECTED_PAIRING_HISTOGRAM_SHA256
    ):
        raise ValueError("silver 1218 off-diagonal pairing histogram changed")
    scope = _mapping(document.get("scope"), "scope")
    if scope.get("other_reflective_norms_exhausted") is not False:
        raise ValueError("silver 1218 reflective-norm caveat disappeared")

    return {
        "diagnostic_sha256": EXPECTED_DIAGNOSTIC_SHA256,
        "pairing_histogram_sha256": EXPECTED_PAIRING_HISTOGRAM_SHA256,
        "facet_count": 1218,
        "elliptic_edge_count": 42297,
        "finite_vertex_count": 5875,
        "ideal_vertex_count": 0,
        "finite_plus_ideal_endpoint_histogram": EXPECTED_FINITE_HISTOGRAM,
        "finite_volume_open_edge_count": 27854,
        "connected_affine_subdiagram_count": 0,
        "full_neighborhood_audit_complete": True,
        "invalid_finite_angle_pair_count": 0,
        "compact_stopping_available": False,
        "finite_volume_stopping_available": False,
        "final_chamber_infinite_volume_proved": False,
        "other_reflective_norms_exhausted": False,
    }


def fresh_replay() -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, str(DIAGNOSE_SOURCE), "--full"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    rebuilt = json.loads(completed.stdout)
    frozen = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    if rebuilt != frozen:
        raise AssertionError("fresh silver 1218 diagnostic differs from archive")
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
