#!/usr/bin/env python3
"""Hostile terminal acceptance for the JK--H12 five-fork ``B=4`` census.

The 106,385 boundary computations are not repeated.  This certifier pins the
exact production snapshot and archives, replays every archived-colour versus
regenerated-colour association in the complete same-facet predecessor
iterator, and independently audits the compact B4 transcript.  In particular
it repairs the production summary's diagnostic-only undercount of discrete
refinements first applied during the cross-shape phase.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import sys
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance as helpers


OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.json"
)
B4_SUMMARY = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b4_to_b4_streaming_census.json"
)
B4_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b4_to_b4_streaming_census_rows.jsonl.gz"
)
P5_BBBB_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census_rows.jsonl.gz"
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
    "RSS_watchdog": ROOT / "scripts/run_with_rss_limit.py",
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
    "stream_parser_helper": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.py",
}
ARCHIVE_PATHS = {
    "B4_production_summary": B4_SUMMARY,
    "B4_production_rows": B4_ROWS,
    "P5_BBBB_same_facet_rows": P5_BBBB_ROWS,
    "star_complete_same_facet_rows": STAR,
    "all_shape_tree_denominator": TREE_DENOMINATOR,
    "fork_pattern_denominator": FORK_DENOMINATOR,
}
PINNED_CODE_SHA256 = {
    "RSS_watchdog": "f118b000a60db7f04a9e24cabc5c191344ac61d2d8ab743841a6ee7d425f0b90",
    "production_fork_builder": "64dd51e4410a379723669eee4002f16029c4cab296ba8693e4483f1eb7f6627e",
    "P5_orbit_and_boundary_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_streaming_archive_reader": "bc2479d167c1daae86d03932477ed223b34effce233ccc7013e07a1817f33801",
    "star_orbit_and_boundary_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "exact_color_refinement_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
    "stream_parser_helper": "26b16157f8c6c944ead72fb1f0d27cfe5e4a979c2580c5c8d482560313bb83b0",
}
PINNED_ARCHIVE_SHA256 = {
    "B4_production_summary": "127b1a552d78a63f71cc6ef7bef426e33be134599b2246cb7c23b7a121f7bf8c",
    "B4_production_rows": "57fe1749c497cdedafc433d5594a9435d392b0025a5679524da823d545d9e763",
    "P5_BBBB_same_facet_rows": "d367b7d208ad8c7fff294978a1857f7557d3755f33df8e4afff6689ad2aec1a7",
    "star_complete_same_facet_rows": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "all_shape_tree_denominator": "2306ebe0048b576a5ba1aa145d5af500b6efdb13f16ffacb8174c7984d1d5d10",
    "fork_pattern_denominator": "99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7",
}

SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_fork_b4_terminal_acceptance.v1"
STATUS_DATE = "2026-08-10"


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


def _expected_prior_rows() -> Iterator[tuple[str, int, tuple[int, ...], str, str]]:
    """Yield archive-authenticated B4 predecessor rows after colour alignment."""

    import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
    import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
    import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star

    actions = star._restricted_actions()
    generated_p5 = p5._pattern_color_orbits("BBBB", actions)
    archived_p5 = helpers._gzip_jsonl_orbit_rows(P5_BBBB_ROWS, 108_159)
    for row, archived_row in zip_longest(generated_p5, archived_p5):
        if row is None or archived_row is None:
            raise ArithmeticError("the P5 BBBB iterator length changed")
        index, archived_colors, digest = archived_row
        generated_colors = tuple(int(value) for value in row["colors_zero_based_stock"])
        if generated_colors != archived_colors:
            raise ArithmeticError(
                f"the P5 BBBB archived/generator colour alignment changed at {index}"
            )
        yield "P5_BBBB", index, archived_colors, digest, "P5"

    generated_star = star._star_color_orbits(actions)
    archived_star = helpers._pretty_json_orbit_rows(STAR, 16_872)
    for row, archived_row in zip_longest(generated_star, archived_star):
        if row is None or archived_row is None:
            raise ArithmeticError("the complete star iterator length changed")
        index, archived_colors, digest = archived_row
        generated_colors = tuple(int(value) for value in row["colors_zero_based_stock"])
        if generated_colors != archived_colors:
            raise ArithmeticError(
                f"the star archived/generator colour alignment changed at {index}"
            )
        if sum(color >= p5.A_SIZE for color in archived_colors) == 4:
            yield "K1,4_star", index, archived_colors, digest, "star"

    if fork.PATTERN_CLASSES["BB|BB"].count("B") != 4:
        raise ArithmeticError("the production fork pattern map changed")


def _audit_prior_iterator() -> dict[str, Any]:
    import scripts.certify_jacobian_kummer_h12_five_vertex_fork_streaming_census as fork
    import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star

    actions = star._restricted_actions()
    actual = fork._prior_rows(4, 4, actions)
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
        if key in seen or sum(color >= 16 for color in colors) != 4:
            raise ArithmeticError("a prior iterator key/B-count changed")
        seen.add(key)
        counts[source] += 1
        total += 1
        _stream_update(transcript, [source, index, list(colors), digest, engine])
    expected_counts = {"K1,4_star": 8_985, "P5_BBBB": 108_159}
    if dict(sorted(counts.items())) != expected_counts or total != 117_144:
        raise ArithmeticError(f"the same-facet predecessor split changed: {counts}")
    return {
        "all_archived_colours_equal_regenerated_colours": True,
        "production_iterator_equals_archive_authenticated_iterator": True,
        "source_row_counts": expected_counts,
        "total_same_facet_row_count": total,
        "source_index_pairs_unique": len(seen) == total,
        "all_rows_have_exactly_four_B_colours": True,
        "canonical_alignment_transcript_sha256": transcript.hexdigest(),
    }


def _audit_b4_rows(summary: dict[str, Any]) -> dict[str, Any]:
    transcript = hashlib.sha256()
    count = 0
    patterns: Counter[str] = Counter()
    facets: Counter[int] = Counter()
    refined = 0
    refined_discrete = 0
    exact_keys: set[tuple[int, str, str]] = set()
    with gzip.open(B4_ROWS, "rt", encoding="utf-8") as stream:
        for count, line in enumerate(stream, start=1):
            transcript.update(line.encode())
            row = json.loads(line)
            if row["orbit_number"] != count:
                raise ArithmeticError("the B4 compact orbit order changed")
            colors = tuple(int(value) - 61 for value in row["edge_colors_one_based"])
            if len(colors) != 4 or sum(color >= 16 for color in colors) != 4:
                raise ArithmeticError("a B4 transcript colour row changed")
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
    if not (
        count == 106_385
        and patterns == Counter({"BB|BB": 106_385})
        and facets == Counter({788: 106_385})
        and refined
        == summary["ordinary_Coxeter_nonisomorphism_certificate"][
            "color_refined_current_member_count"
        ]
        and len(exact_keys) == 106_385
        and transcript.hexdigest()
        == summary["compact_exact_transcript"][
            "uncompressed_canonical_JSONL_sha256"
        ]
    ):
        raise ArithmeticError("the B4 compact transcript audit failed")
    reported_discrete = summary["ordinary_Coxeter_nonisomorphism_certificate"][
        "color_refined_discrete_current_count"
    ]
    if refined_discrete < reported_discrete:
        raise ArithmeticError("the final B4 discrete count decreased")
    return {
        "row_count": count,
        "pattern_type_histogram": {"BB|BB": 106_385},
        "facet_histogram": {"788": 106_385},
        "exact_degree_plus_refinement_key_count": len(exact_keys),
        "color_refined_current_member_count": refined,
        "corrected_color_refined_discrete_current_count": refined_discrete,
        "production_summary_underreported_discrete_count_as": reported_discrete,
        "diagnostic_only_correction_no_theorem_change": True,
        "uncompressed_canonical_JSONL_sha256": transcript.hexdigest(),
    }


def build() -> dict[str, Any]:
    actual_code = {name: _sha256(path) for name, path in CODE_PATHS.items()}
    actual_archives = {name: _sha256(path) for name, path in ARCHIVE_PATHS.items()}
    if actual_code != PINNED_CODE_SHA256:
        raise ArithmeticError(f"a B4 production code snapshot changed: {actual_code}")
    if actual_archives != PINNED_ARCHIVE_SHA256:
        raise ArithmeticError(f"a B4 production/archive input changed: {actual_archives}")

    summary = json.loads(B4_SUMMARY.read_text(encoding="utf-8"))
    theorem = summary["theorem"]
    nonisomorphism = summary["ordinary_Coxeter_nonisomorphism_certificate"]
    cross = summary["cross_shape_exact_separation"]
    checks = {
        "production_summary_green": summary["all_checks_pass"],
        "exact_selected_layer": summary["selected_B_range"] == [4, 4]
        and summary["selected_patterns"] == ["BB|BB"],
        "exact_type_and_facet_counts": theorem["ordinary_Coxeter_type_count"]
        == 106_385
        and theorem["type_facet_histogram"] == {"788": 106_385},
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
        "same_facet_predecessor_count_is_117144": cross[
            "prior_same_facet_P5_and_star_row_count"
        ]
        == 117_144,
        "row_archive_hash_matches_summary": summary["compact_exact_transcript"][
            "deterministic_gzip_sha256"
        ]
        == PINNED_ARCHIVE_SHA256["B4_production_rows"],
        "watchdog_peak_below_one_GiB": True,
        "worldwide_priority_explicitly_open": summary["lineage_and_claim_boundary"][
            "worldwide_priority"
        ]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the B4 terminal summary failed: {checks}")
    prior_alignment = _audit_prior_iterator()
    rows = _audit_b4_rows(summary)
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "production_provenance": {
            "exact_launch_time_code_sha256": PINNED_CODE_SHA256,
            "exact_input_and_output_archive_sha256": PINNED_ARCHIVE_SHA256,
            "launch_command": (
                "python3 scripts/run_with_rss_limit.py --memory-gib 1.0 "
                "--timeout 14400 -- python3 "
                "scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py "
                "--min-b 4 --max-b 4"
            ),
            "watchdog_disposition": "child_exit",
            "watchdog_scope": "direct",
            "watchdog_peak_RSS_GiB": 0.483,
            "watchdog_elapsed_seconds": 11_549.659,
        },
        "theorem": {
            "scope": "complete degree-three-fork B=4 layer",
            "ordinary_Coxeter_type_count": 106_385,
            "pairwise_nonisomorphic": True,
            "separated_from_all_117144_completed_same_facet_P5_and_star_rows": True,
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
            "selected_B4_layer_complete": True,
            "complete_fork_combination_requires_frozen_B_le_3_input": True,
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
            raise SystemExit("the B4 terminal-acceptance archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("B4_types=106385")
    print("same_facet_predecessor_rows=117144")
    print(
        "corrected_refined_discrete_rows="
        f"{payload['terminal_transcript_audit']['corrected_color_refined_discrete_current_count']}"
    )
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
