#!/usr/bin/env python3
"""Certify exact ordinary-Coxeter separation of JK-H12 stars and ABBB paths.

The two frozen five-vertex archives have one overlapping facet header: the
6,117 K1,4 rows and all 36,247 ABBB P5 rows with 779 facets.  Their production
archives recorded a zero SHA-256 screen on the exact degree-profile invariant,
but deliberately did not promote it to a direct exact comparison.  This
certifier rebuilds every relevant boundary and compares the underlying bytes.

Hashes are used only to pin inputs and to check that each rebuilt byte string
is the one represented in its source archive.  Cross-shape separation itself
is the literal Python ``bytes`` membership test ``degree_key in star_keys``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_abbb_census as abbb
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_small_subtree_census as seed


STAR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py"
STAR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
P5_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py"
P5_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
ABBB_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py"
ABBB_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_star_abbb_exact_separation.json"

SCHEMA = "coxeter.jacobian_kummer.h12_star_abbb_exact_separation.v1"
STATUS_DATE = "2026-08-08"
FACETS = 779
STAR_779_COUNT = 6_117
ABBB_COUNT = 36_247
COMPLETE_STAR_COUNT = 16_872
EXACT_UNION_COUNT = COMPLETE_STAR_COUNT + ABBB_COUNT
DEGREE_KEY_BYTES = FACETS * 5 * 2

PINNED_SHA256 = {
    "five_vertex_star_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "five_vertex_star_archive": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "P5_B_le_2_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_B_le_2_archive": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
    "P5_ABBB_builder": "6217c4902b8a66bbb4cef57e4a2edc188da723d6ceffcc2cad82cd8d172eb18e",
    "P5_ABBB_archive": "b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _pinned_payloads() -> tuple[dict[str, Any], dict[str, Any]]:
    paths = {
        "five_vertex_star_builder": STAR_SCRIPT,
        "five_vertex_star_archive": STAR_ARCHIVE,
        "P5_B_le_2_builder": P5_SCRIPT,
        "P5_B_le_2_archive": P5_ARCHIVE,
        "P5_ABBB_builder": ABBB_SCRIPT,
        "P5_ABBB_archive": ABBB_ARCHIVE,
    }
    actual = {name: _sha256(path) for name, path in paths.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen star/ABBB input changed: {actual}")

    # These private production gates transitively pin the seed, action,
    # boundary engine, ordinary serializer, and predecessor archives.
    star._pinned_predecessors()
    abbb._frozen_gate()

    star_payload = json.loads(STAR_ARCHIVE.read_text(encoding="utf-8"))
    abbb_payload = json.loads(ABBB_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        star_payload["all_checks_pass"]
        and all(star_payload["checks"].values())
        and star_payload["theorem"]["ordinary_Coxeter_type_count"]
        == COMPLETE_STAR_COUNT
        and star_payload["theorem"]["pairwise_nonisomorphic"]
        and star_payload["theorem"]["type_facet_histogram"][str(FACETS)]
        == STAR_779_COUNT
        and abbb_payload["all_checks_pass"]
        and all(abbb_payload["checks"].values())
        and abbb_payload["theorem"]["ordinary_Coxeter_type_count"] == ABBB_COUNT
        and abbb_payload["theorem"]["pairwise_nonisomorphic_within_ABBB_shard"]
        and abbb_payload["theorem"]["type_facet_histogram"]
        == {str(FACETS): ABBB_COUNT}
    ):
        raise ArithmeticError("a frozen star/ABBB theorem is not green")
    return star_payload, abbb_payload


def _stock_colors(row: dict[str, Any]) -> tuple[int, int, int, int]:
    colors = tuple(
        int(color) - star.STOCK_START - 1 for color in row["edge_colors_one_based"]
    )
    if not (
        len(colors) == 4
        and all(0 <= color < star.STOCK_SIZE for color in colors)
    ):
        raise ArithmeticError("an archived five-vertex color word left the stock")
    return colors  # type: ignore[return-value]


def _stream_digest_update(digest: Any, degree_key: bytes) -> None:
    digest.update(len(degree_key).to_bytes(4, "big"))
    digest.update(degree_key)


def _rebuild_degree_key(
    engine: Any, row: dict[str, Any], *, expected_shape: str
) -> bytes:
    boundary, boundary_data = engine.boundary(_stock_colors(row))
    _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
    if not (
        row["tree_shape"] == expected_shape
        and row["facet_count"] == FACETS
        and boundary_data["facet_count"] == FACETS
        and boundary_data["marked_boundary_root_sha256"]
        == row["marked_boundary_root_sha256"]
        and ordinary_data["ordinary_degree_profile_multiset_sha256"]
        == row["ordinary_degree_profile_multiset_sha256"]
        == hashlib.sha256(degree_key).hexdigest()
        and ordinary_data["ordinary_Coxeter_relation_histogram"]
        == row["ordinary_Coxeter_relation_histogram"]
        and len(degree_key) == DEGREE_KEY_BYTES
    ):
        raise ArithmeticError("an exact archived degree-profile replay changed")
    return degree_key


def _progress(label: str, completed: int, total: int, started: float) -> None:
    print(
        f"{label}_progress={completed}/{total} "
        f"elapsed_seconds={time.perf_counter() - started:.1f}",
        flush=True,
    )


def build(*, verbose: bool = False) -> dict[str, Any]:
    star_payload, abbb_payload = _pinned_payloads()
    _carrier, gram, roots, interior, _bounded = seed._load_inputs()
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)

    star_rows = tuple(
        row
        for row in star_payload["orbit_representatives"]
        if int(row["facet_count"]) == FACETS
    )
    abbb_rows = tuple(abbb_payload["orbit_representatives"])
    if not (
        len(star_rows) == STAR_779_COUNT
        and len(abbb_rows) == ABBB_COUNT
        and all(row["edge_class_word"] == "ABBB" for row in abbb_rows)
    ):
        raise ArithmeticError("the archived cross-shape denominator changed")

    star_degree_keys: set[bytes] = set()
    star_stream_digest = hashlib.sha256()
    started = time.perf_counter()
    for index, row in enumerate(star_rows, start=1):
        degree_key = _rebuild_degree_key(star_engine, row, expected_shape="K1,4")
        star_degree_keys.add(degree_key)
        _stream_digest_update(star_stream_digest, degree_key)
        if verbose and (index % 500 == 0 or index == len(star_rows)):
            _progress("star_779", index, len(star_rows), started)
    cross_matches = 0
    abbb_stream_digest = hashlib.sha256()
    started = time.perf_counter()
    for index, row in enumerate(abbb_rows, start=1):
        degree_key = _rebuild_degree_key(p5_engine, row, expected_shape="P5")
        _stream_digest_update(abbb_stream_digest, degree_key)
        # This is the theorem-producing comparison: exact bytes, not hashes.
        if degree_key in star_degree_keys:
            cross_matches += 1
        if verbose and (index % 500 == 0 or index == len(abbb_rows)):
            _progress("ABBB", index, len(abbb_rows), started)
    checks = {
        "all_frozen_builder_and_archive_SHA256_values_match": True,
        "frozen_star_and_ABBB_theorems_are_green": True,
        "exactly_6117_star_779_rows_rebuilt": len(star_rows) == STAR_779_COUNT,
        "all_36247_ABBB_rows_rebuilt": len(abbb_rows) == ABBB_COUNT,
        "every_rebuilt_exact_byte_key_matches_its_archived_SHA256": True,
        "cross_shape_comparison_uses_literal_degree_profile_bytes": True,
        "star_779_and_ABBB_exact_degree_profile_key_sets_are_disjoint": (
            cross_matches == 0
        ),
        "complete_star_plus_ABBB_union_count_is_53119": EXACT_UNION_COUNT == 53_119,
        "claim_is_project_scoped_and_worldwide_priority_remains_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"star/ABBB exact separation failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": dict(PINNED_SHA256),
        "exact_invariant": {
            "name": "sorted ordinary-label degree-profile bytes",
            "ordinary_relation_order": list(star.LABEL_NAMES),
            "profile_shape": [FACETS, 5],
            "entry_encoding": "little-endian unsigned 16-bit integer",
            "byte_length_per_polytope": DEGREE_KEY_BYTES,
            "isomorphism_invariance": (
                "an ordinary Coxeter-matrix isomorphism permutes facet profiles, "
                "so their lexicographically sorted byte serialization is unchanged"
            ),
            "cross_comparison": "literal Python bytes equality/membership",
            "hash_role": (
                "input pinning, per-row archive linkage, and compact stream digest only; "
                "not the cross-shape separation test"
            ),
        },
        "replay": {
            "star_779_row_count": len(star_rows),
            "star_779_unique_exact_degree_profile_key_count": len(star_degree_keys),
            "ABBB_row_count": len(abbb_rows),
            "cross_shape_exact_byte_match_count": cross_matches,
            "star_779_archive_order_exact_key_stream_sha256": (
                star_stream_digest.hexdigest()
            ),
            "ABBB_archive_order_exact_key_stream_sha256": (
                abbb_stream_digest.hexdigest()
            ),
        },
        "theorem": {
            "scope": (
                "complete K1,4 five-vertex layer versus complete ABBB P5 shard "
                "inside the frozen JK-H12 W61--W172 stock"
            ),
            "overlapping_facet_header": FACETS,
            "star_779_type_count": STAR_779_COUNT,
            "ABBB_type_count": ABBB_COUNT,
            "exact_cross_shape_ordinary_Coxeter_isomorphism_count": 0,
            "star_779_plus_ABBB_pairwise_nonisomorphic_type_count": (
                STAR_779_COUNT + ABBB_COUNT
            ),
            "complete_star_type_count": COMPLETE_STAR_COUNT,
            "complete_star_plus_ABBB_pairwise_nonisomorphic_type_count": (
                EXACT_UNION_COUNT
            ),
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiple": 5,
        },
        "carrier_and_novelty_boundary": {
            "direction": "SURFACE_TO_POLYTOPE carrier-preserving redoubling census",
            "same_weak_full_Picard_K3_carrier": (
                "U(2)+D5(-1)^2+<-4>"
            ),
            "full_nef_effectivity_or_surface_action_claim": False,
            "raw_POLYTOPE_TO_SURFACE_recognizer_claim": False,
            "project_cross_shape_accounting_promoted": True,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
            "full_P5_or_all_five_vertex_shape_classification_claim": False,
        },
        "checks": checks,
    }


def _probe(limit: int) -> None:
    if limit <= 0:
        raise SystemExit("--probe must be positive")
    star_payload, abbb_payload = _pinned_payloads()
    _carrier, gram, roots, interior, _bounded = seed._load_inputs()
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    star_rows: Iterable[dict[str, Any]] = (
        row
        for row in star_payload["orbit_representatives"]
        if int(row["facet_count"]) == FACETS
    )
    for index, row in enumerate(star_rows):
        if index >= limit:
            break
        _rebuild_degree_key(star_engine, row, expected_shape="K1,4")
    for row in abbb_payload["orbit_representatives"][:limit]:
        _rebuild_degree_key(p5_engine, row, expected_shape="P5")
    print(f"probe_rows_per_shape={limit}")
    print("probe=PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.probe is not None:
        _probe(arguments.probe)
        return

    payload = build(verbose=True)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"archive is missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact star/ABBB rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(
        "cross_shape_exact_byte_match_count="
        f"{payload['replay']['cross_shape_exact_byte_match_count']}"
    )
    print(
        "complete_star_plus_ABBB_types="
        f"{payload['theorem']['complete_star_plus_ABBB_pairwise_nonisomorphic_type_count']}"
    )
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
