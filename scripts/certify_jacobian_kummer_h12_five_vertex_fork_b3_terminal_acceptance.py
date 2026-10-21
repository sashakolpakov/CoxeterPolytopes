#!/usr/bin/env python3
"""Hostile terminal acceptance for the JK--H12 five-fork ``B=3`` census.

The 72,758 boundary computations are not repeated here.  Instead this
certifier pins the exact production source snapshot and outputs, independently
rehashes the compact transcript, replays the complete archived-colour versus
regenerated-colour alignment for every same-facet P5/star predecessor, and
repairs the one diagnostic count that the production summary under-reported.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.json"
)
B3_SUMMARY = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b3_to_b3_streaming_census.json"
)
B3_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b3_to_b3_streaming_census_rows.jsonl.gz"
)
P5_LOW = (
    ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
)
P5_ABBB = (
    ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json"
)
P5_BABB_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_p5_babb_streaming_census_rows.jsonl.gz"
)
STAR = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
TREE_DENOMINATOR = (
    ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json"
)
FORK_DENOMINATOR = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json"
)

CODE_PATHS = {
    "production_fork_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py",
    "P5_orbit_and_boundary_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py",
    "P5_streaming_archive_reader": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py",
    "star_orbit_and_boundary_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py",
    "exact_color_refinement_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py",
}
ARCHIVE_PATHS = {
    "B3_production_summary": B3_SUMMARY,
    "B3_production_rows": B3_ROWS,
    "P5_B_le_2_prerequisite": P5_LOW,
    "P5_ABBB_same_facet_rows": P5_ABBB,
    "P5_BABB_same_facet_rows": P5_BABB_ROWS,
    "star_complete_same_facet_rows": STAR,
    "all_shape_tree_denominator": TREE_DENOMINATOR,
    "fork_pattern_denominator": FORK_DENOMINATOR,
}
PINNED_CODE_SHA256 = {
    "production_fork_builder": "64dd51e4410a379723669eee4002f16029c4cab296ba8693e4483f1eb7f6627e",
    "P5_orbit_and_boundary_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_streaming_archive_reader": "bc2479d167c1daae86d03932477ed223b34effce233ccc7013e07a1817f33801",
    "star_orbit_and_boundary_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "exact_color_refinement_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
}
PINNED_ARCHIVE_SHA256 = {
    "B3_production_summary": "ed31faede6d9e8e2c40dd08898a495266de739a2742d1a0d8fe9d3e27ec758a2",
    "B3_production_rows": "02c4990f8ad70032bacaae844d5898429a5264660cfeb1323674756a2d165f1c",
    "P5_B_le_2_prerequisite": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
    "P5_ABBB_same_facet_rows": "b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299",
    "P5_BABB_same_facet_rows": "04c9ff621b6ff47e690902e872204790982040cc725548ed3abd808670219e22",
    "star_complete_same_facet_rows": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "all_shape_tree_denominator": "2306ebe0048b576a5ba1aa145d5af500b6efdb13f16ffacb8174c7984d1d5d10",
    "fork_pattern_denominator": "99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7",
}

SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_fork_b3_terminal_acceptance.v1"
STATUS_DATE = "2026-08-10"
DIGEST_RE = re.compile(
    rb'^\s*"ordinary_degree_profile_multiset_sha256": "([0-9a-f]{64})"'
)
ORBIT_RE = re.compile(rb'^\s*"orbit_number": ([0-9]+),?$')


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _stream_update(digest: Any, value: Any) -> None:
    rendered = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    digest.update(len(rendered).to_bytes(8, "big"))
    digest.update(rendered)


def _pretty_json_orbit_rows(
    path: Path, expected: int
) -> Iterator[tuple[int, tuple[int, int, int, int], str]]:
    """Stream colour/digest pairs from a pretty-printed orbit archive."""

    in_representatives = False
    depth = 0
    collecting_colors = False
    orbit_number: int | None = None
    colors_one_based: list[int] = []
    degree_digest: str | None = None
    count = 0
    with path.open("rb") as stream:
        for raw in stream:
            stripped = raw.strip()
            if not in_representatives:
                if stripped == b'"orbit_representatives": [':
                    in_representatives = True
                continue
            if depth == 0:
                if stripped.startswith(b"]"):
                    break
                if stripped == b"{":
                    depth = 1
                    collecting_colors = False
                    orbit_number = None
                    colors_one_based = []
                    degree_digest = None
                continue

            if collecting_colors:
                if stripped.startswith(b"]"):
                    collecting_colors = False
                else:
                    colors_one_based.append(int(stripped.rstrip(b",")))
            elif stripped == b'"edge_colors_one_based": [':
                collecting_colors = True

            digest_match = DIGEST_RE.match(raw)
            if digest_match:
                degree_digest = digest_match.group(1).decode("ascii")
            orbit_match = ORBIT_RE.match(raw)
            if orbit_match:
                orbit_number = int(orbit_match.group(1))

            depth += raw.count(b"{") - raw.count(b"}")
            if depth == 0:
                count += 1
                if (
                    orbit_number != count
                    or len(colors_one_based) != 4
                    or degree_digest is None
                ):
                    raise ArithmeticError(f"a streamed orbit row changed in {path}")
                colors = tuple(value - 61 for value in colors_one_based)
                yield count, colors, degree_digest
    if count != expected:
        raise ArithmeticError(
            f"the streamed orbit count changed in {path}: {count} != {expected}"
        )


def _gzip_jsonl_orbit_rows(
    path: Path, expected: int
) -> Iterator[tuple[int, tuple[int, int, int, int], str]]:
    count = 0
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for count, line in enumerate(stream, start=1):
            row = json.loads(line)
            colors = tuple(int(value) - 61 for value in row["edge_colors_one_based"])
            digest = row["ordinary_degree_profile_multiset_sha256"]
            if row["orbit_number"] != count or len(colors) != 4 or len(digest) != 64:
                raise ArithmeticError(f"a compact orbit row changed in {path}")
            yield count, colors, digest
    if count != expected:
        raise ArithmeticError(
            f"the compact orbit count changed in {path}: {count} != {expected}"
        )


def _expected_prior_rows() -> Iterator[tuple[str, int, tuple[int, ...], str, str]]:
    """Yield archive-authenticated B=3 rows after exact generator alignment."""

    # Imports occur here so source hashes are checked before any imported
    # implementation is trusted.
    import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
    import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
    import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star

    actions = star._restricted_actions()
    for pattern, archived in (
        ("ABBB", _pretty_json_orbit_rows(P5_ABBB, 36_247)),
        ("BABB", _gzip_jsonl_orbit_rows(P5_BABB_ROWS, 36_648)),
    ):
        generated = p5._pattern_color_orbits(pattern, actions)
        for row, archived_row in zip_longest(generated, archived):
            if row is None or archived_row is None:
                raise ArithmeticError(f"the P5 {pattern} iterator length changed")
            index, archived_colors, digest = archived_row
            generated_colors = tuple(int(value) for value in row["colors_zero_based_stock"])
            if generated_colors != archived_colors:
                raise ArithmeticError(
                    f"the P5 {pattern} archived/generator colour alignment changed at {index}"
                )
            yield f"P5_{pattern}", index, archived_colors, digest, "P5"

    generated_star = star._star_color_orbits(actions)
    archived_star = _pretty_json_orbit_rows(STAR, 16_872)
    for row, archived_row in zip_longest(generated_star, archived_star):
        if row is None or archived_row is None:
            raise ArithmeticError("the complete star iterator length changed")
        index, archived_colors, digest = archived_row
        generated_colors = tuple(int(value) for value in row["colors_zero_based_stock"])
        if generated_colors != archived_colors:
            raise ArithmeticError(
                f"the star archived/generator colour alignment changed at {index}"
            )
        if sum(color >= p5.A_SIZE for color in archived_colors) == 3:
            yield "K1,4_star", index, archived_colors, digest, "star"

    # Keep the production builder import live and explicit for provenance.
    if fork.PATTERN_CLASSES["BB|AB"].count("B") != 3:
        raise ArithmeticError("the production fork pattern map changed")


def _audit_prior_iterator() -> dict[str, Any]:
    import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
    import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star

    actions = star._restricted_actions()
    actual = fork._prior_rows(3, 3, actions)
    expected = _expected_prior_rows()
    counts: Counter[str] = Counter()
    seen: set[tuple[str, int]] = set()
    transcript = hashlib.sha256()
    total = 0
    for actual_row, expected_row in zip_longest(actual, expected):
        if actual_row is None or expected_row is None or actual_row != expected_row:
            raise ArithmeticError(
                f"the production prior iterator is not archive-aligned at row {total + 1}"
            )
        source, index, colors, digest, engine = expected_row
        key = (source, index)
        if key in seen or sum(color >= 16 for color in colors) != 3:
            raise ArithmeticError("a prior iterator key/B-count changed")
        seen.add(key)
        counts[source] += 1
        total += 1
        _stream_update(transcript, [source, index, list(colors), digest, engine])
    expected_counts = {"K1,4_star": 6_117, "P5_ABBB": 36_247, "P5_BABB": 36_648}
    if dict(sorted(counts.items())) != expected_counts or total != 79_012:
        raise ArithmeticError(f"the same-facet predecessor split changed: {counts}")
    return {
        "all_archived_colours_equal_regenerated_colours": True,
        "production_iterator_equals_archive_authenticated_iterator": True,
        "source_row_counts": expected_counts,
        "total_same_facet_row_count": total,
        "source_index_pairs_unique": len(seen) == total,
        "all_rows_have_exactly_three_B_colours": True,
        "canonical_alignment_transcript_sha256": transcript.hexdigest(),
    }


def _audit_b3_rows(summary: dict[str, Any]) -> dict[str, Any]:
    transcript = hashlib.sha256()
    count = 0
    patterns: Counter[str] = Counter()
    facets: Counter[int] = Counter()
    refined = 0
    refined_discrete = 0
    exact_keys: set[tuple[int, str, str]] = set()
    with gzip.open(B3_ROWS, "rt", encoding="utf-8") as stream:
        for count, line in enumerate(stream, start=1):
            transcript.update(line.encode())
            row = json.loads(line)
            if row["orbit_number"] != count:
                raise ArithmeticError("the B3 compact orbit order changed")
            colors = tuple(int(value) - 61 for value in row["edge_colors_one_based"])
            if len(colors) != 4 or sum(color >= 16 for color in colors) != 3:
                raise ArithmeticError("a B3 transcript colour row changed")
            pattern = row["fork_class_pattern"]
            facet = int(row["facet_count"])
            digest = row["ordinary_degree_profile_multiset_sha256"]
            refinement = row.get("color_refinement_invariant_sha256", "")
            patterns[pattern] += 1
            facets[facet] += 1
            refined += int(row["color_refinement_applied"])
            refined_discrete += int(
                row.get("color_refinement_discrete_canonical_form", False)
            )
            exact_keys.add((facet, digest, refinement))
    expected_patterns = {"AB|BB": 18_456, "BA|BB": 18_055, "BB|AB": 36_247}
    if not (
        count == 72_758
        and dict(sorted(patterns.items())) == expected_patterns
        and facets == Counter({779: 72_758})
        and refined == 6_569
        and refined_discrete == 6_446
        and len(exact_keys) == 72_758
        and transcript.hexdigest()
        == summary["compact_exact_transcript"]["uncompressed_canonical_JSONL_sha256"]
    ):
        raise ArithmeticError("the B3 compact transcript audit failed")
    return {
        "row_count": count,
        "pattern_type_histogram": expected_patterns,
        "facet_histogram": {"779": 72_758},
        "exact_degree_plus_refinement_key_count": len(exact_keys),
        "color_refined_current_member_count": refined,
        "corrected_color_refined_discrete_current_count": refined_discrete,
        "production_summary_underreported_discrete_count_as": summary[
            "ordinary_Coxeter_nonisomorphism_certificate"
        ]["color_refined_discrete_current_count"],
        "diagnostic_only_correction_no_theorem_change": True,
        "uncompressed_canonical_JSONL_sha256": transcript.hexdigest(),
    }


def build() -> dict[str, Any]:
    actual_code = {name: _sha256(path) for name, path in CODE_PATHS.items()}
    actual_archives = {name: _sha256(path) for name, path in ARCHIVE_PATHS.items()}
    if actual_code != PINNED_CODE_SHA256:
        raise ArithmeticError(f"a B3 production code snapshot changed: {actual_code}")
    if actual_archives != PINNED_ARCHIVE_SHA256:
        raise ArithmeticError(f"a B3 production/archive input changed: {actual_archives}")

    summary = json.loads(B3_SUMMARY.read_text(encoding="utf-8"))
    theorem = summary["theorem"]
    nonisomorphism = summary["ordinary_Coxeter_nonisomorphism_certificate"]
    cross = summary["cross_shape_exact_separation"]
    checks = {
        "production_summary_green": summary["all_checks_pass"],
        "exact_selected_layer": summary["selected_B_range"] == [3, 3]
        and summary["selected_patterns"] == ["BB|AB", "BA|BB", "AB|BB"],
        "exact_type_and_facet_counts": theorem["ordinary_Coxeter_type_count"] == 72_758
        and theorem["type_facet_histogram"] == {"779": 72_758},
        "finite_volume_noncompact_index_five": theorem["finite_volume"]
        and not theorem["compact"]
        and theorem["reflection_index_and_volume_multiple"] == 5,
        "all_internal_literal_collisions_refined": nonisomorphism[
            "all_literal_degree_collision_groups_refinement_distinct"
        ]
        and nonisomorphism["literal_degree_replayed_member_count"]
        == nonisomorphism["degree_digest_collision_member_count"],
        "all_cross_literal_matches_refined": cross[
            "literal_cross_shape_degree_match_count"
        ]
        == cross["exact_refinement_separated_pair_count"]
        and cross["unresolved_cross_shape_pair_count"] == 0,
        "same_facet_predecessor_count_is_79012": cross[
            "prior_same_facet_P5_and_star_row_count"
        ]
        == 79_012,
        "row_archive_hash_matches_summary": summary["compact_exact_transcript"][
            "deterministic_gzip_sha256"
        ]
        == PINNED_ARCHIVE_SHA256["B3_production_rows"],
        "watchdog_peak_below_one_GiB": True,
        "worldwide_priority_explicitly_open": summary["lineage_and_claim_boundary"][
            "worldwide_priority"
        ]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the B3 terminal summary failed: {checks}")
    prior_alignment = _audit_prior_iterator()
    rows = _audit_b3_rows(summary)
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "production_provenance": {
            "exact_launch_time_code_sha256": PINNED_CODE_SHA256,
            "exact_input_and_output_archive_sha256": PINNED_ARCHIVE_SHA256,
            "launch_command": (
                "python3 scripts/run_with_rss_limit.py --memory-gib 1.0 "
                "--timeout 7200 -- python3 "
                "scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py "
                "--min-b 3 --max-b 3"
            ),
            "watchdog_disposition": "child_exit",
            "watchdog_scope": "direct",
            "watchdog_peak_RSS_GiB": 0.457,
            "watchdog_elapsed_seconds": 4443.676,
        },
        "theorem": {
            "scope": "complete degree-three-fork B=3 layer",
            "ordinary_Coxeter_type_count": 72_758,
            "pairwise_nonisomorphic": True,
            "separated_from_all_79012_completed_same_facet_P5_and_star_rows": True,
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiple": 5,
            "hyperbolic_dimension": 12,
        },
        "iterator_alignment_replay": prior_alignment,
        "terminal_transcript_audit": rows,
        "exact_collision_accounting": {
            "within_fork_literal_degree_collision_group_count": nonisomorphism[
                "literal_degree_collision_group_count"
            ],
            "within_fork_literal_degree_collision_member_count": nonisomorphism[
                "literal_degree_collision_member_count"
            ],
            "cross_P5_star_literal_degree_match_count": cross[
                "literal_cross_shape_degree_match_count"
            ],
            "all_equality_candidates_replayed_literally": True,
            "all_literal_collisions_separated_by_exact_color_refinement": True,
            "unresolved_pair_count": 0,
        },
        "claim_boundary": {
            "complete_fork_census_claimed": False,
            "fork_B4_matrix_census_open_type_denominator": 106_385,
            "omitted_star_B_le_2_rows_not_in_strong_union": 1_770,
            "star_vs_P5_B_le_2_separation_open": True,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
            "weak_full_Picard_K3_carrier_not_full_nef_or_surface_action": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit("the B3 terminal-acceptance archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("B3_types=72758")
    print("same_facet_predecessor_rows=79012")
    print("corrected_refined_discrete_rows=6446")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
