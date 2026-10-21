#!/usr/bin/env python3
"""Certify the complete ABBB five-vertex P5 shard in the JK-H12 factory.

This is the next immutable binary class-word shard after the frozen B<=2
P5 theorem.  It imports that final exact engine without changing its source,
enumerates all 36,247 ABBB seed-and-reversal orbits, and rebuilds every
779-facet ordinary Coxeter boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as frozen


FROZEN_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py"
FROZEN_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json"
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_p5_abbb_census.v1"
STATUS_DATE = "2026-08-08"
PATTERNS = ("ABBB",)

PINNED_FROZEN_SHA256 = {
    "P5_B_le_2_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_B_le_2_archive": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _frozen_gate() -> dict[str, Any]:
    actual = {
        "P5_B_le_2_builder": _sha256(FROZEN_SCRIPT),
        "P5_B_le_2_archive": _sha256(FROZEN_ARCHIVE),
    }
    if actual != PINNED_FROZEN_SHA256:
        raise ArithmeticError(f"the frozen B<=2 P5 theorem changed: {actual}")
    payload = json.loads(FROZEN_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        payload["all_checks_pass"]
        and all(payload["checks"].values())
        and payload["theorem"]["ordinary_Coxeter_type_count"] == 20_693
        and payload["theorem"]["pairwise_nonisomorphic_within_selected_layer"]
        and payload["theorem"]["type_facet_histogram"]
        == {"752": 142, "761": 2_091, "770": 18_460}
        and payload["theorem"]["finite_volume"]
        and not payload["theorem"]["compact"]
    ):
        raise ArithmeticError("the frozen B<=2 P5 theorem is not green")
    return payload


def build(*, verbose: bool = False) -> dict[str, Any]:
    previous = _frozen_gate()
    payload = frozen.build(patterns=PATTERNS, verbose=verbose)
    if not (
        payload["all_checks_pass"]
        and payload["coloring_denominator"]["selected_binary_class_word_shards"]
        == ["ABBB"]
        and payload["coloring_denominator"]["selected_raw_unrooted_proper_colorings"]
        == 13_862_400
        and payload["coloring_denominator"]["selected_seed_and_reversal_orbit_count"]
        == 36_247
        and payload["theorem"]["ordinary_Coxeter_type_count"] == 36_247
        and payload["theorem"]["pairwise_nonisomorphic_within_selected_layer"]
        and payload["theorem"]["type_facet_histogram"] == {"779": 36_247}
        and payload["pilot_recovery"]["eligible_B_0_1_2_pilot_count"] == 1
    ):
        raise ArithmeticError("the ABBB production tranche did not close")

    combined_count = (
        previous["theorem"]["ordinary_Coxeter_type_count"]
        + payload["theorem"]["ordinary_Coxeter_type_count"]
    )
    remaining_count = (
        frozen.EXPECTED_PATTERN_ORBITS["BABB"]
        + frozen.EXPECTED_PATTERN_ORBITS["BBBB"]
    )
    combined_raw = (
        previous["coloring_denominator"]["selected_raw_unrooted_proper_colorings"]
        + payload["coloring_denominator"]["selected_raw_unrooted_proper_colorings"]
    )
    completed_patterns = [*frozen.PRODUCTION_PATTERNS, "ABBB"]
    combined_facet_histogram = {
        **previous["theorem"]["type_facet_histogram"],
        "779": 36_247,
    }
    wrapper_checks = {
        "frozen_B_le_2_builder_and_archive_match_pinned_SHA256": True,
        "complete_ABBB_explicit_and_Burnside_denominator_is_36247": True,
        "all_36247_ABBB_boundaries_rebuilt": True,
        "all_36247_ABBB_ordinary_matrices_are_pairwise_nonisomorphic": True,
        "every_ABBB_type_has_779_facets": True,
        "ABBB_is_exactly_disjoint_from_frozen_B_le_2_by_facet_count": not (
            {779} & {int(value) for value in previous["theorem"]["type_facet_histogram"]}
        ),
        "combined_completed_P5_count_is_56940": combined_count == 56_940,
        "remaining_BABB_BBBB_denominator_is_144807": remaining_count == 144_807,
        "predecessor_B_le_2_scope_pinned_and_green": True,
        "production_scope_is_exactly_ABBB": payload["coloring_denominator"][
            "selected_binary_class_word_shards"
        ]
        == ["ABBB"],
        "complete_P5_census_not_claimed": remaining_count > 0,
        "worldwide_priority_explicitly_open": True,
        "cross_star_hash_screen_not_promoted_to_exact_theorem": True,
    }
    if not all(wrapper_checks.values()):
        raise ArithmeticError(f"the ABBB wrapper checks failed: {wrapper_checks}")

    payload["schema"] = SCHEMA
    payload["status_date"] = STATUS_DATE
    payload["pinned_inputs"] = {
        **payload["pinned_inputs"],
        **PINNED_FROZEN_SHA256,
    }
    payload["theorem"]["scope"] = "complete ABBB five-vertex P5 class-word shard"
    payload["theorem"]["pairwise_nonisomorphic_within_ABBB_shard"] = payload[
        "theorem"
    ].pop("pairwise_nonisomorphic_within_selected_layer")
    payload["coloring_denominator"][
        "complement_of_selected_shard_patterns"
    ] = payload["coloring_denominator"].pop("remaining_disjoint_patterns")
    payload["coloring_denominator"][
        "complement_of_selected_shard_orbit_count"
    ] = payload["coloring_denominator"].pop("remaining_disjoint_orbit_count")
    payload["coloring_denominator"]["previously_frozen_binary_class_word_shards"] = list(
        frozen.PRODUCTION_PATTERNS
    )
    payload["coloring_denominator"]["globally_completed_binary_class_word_shards"] = (
        completed_patterns
    )
    payload["coloring_denominator"]["globally_completed_raw_unrooted_colorings"] = (
        combined_raw
    )
    payload["coloring_denominator"]["globally_completed_seed_and_reversal_orbits"] = (
        combined_count
    )
    payload["coloring_denominator"]["globally_remaining_binary_class_word_shards"] = [
        "BABB",
        "BBBB",
    ]
    payload["coloring_denominator"]["globally_remaining_seed_and_reversal_orbits"] = (
        remaining_count
    )
    payload["combined_exact_P5_progress"] = {
        "frozen_B_le_2_type_count": 20_693,
        "new_ABBB_type_count": 36_247,
        "combined_completed_type_count": combined_count,
        "combined_completed_type_facet_histogram": combined_facet_histogram,
        "exact_cross_tranche_nonisomorphism_reason": (
            "ABBB has 779 facets, whereas frozen B<=2 has 752, 761, or 770"
        ),
        "four_of_five_original_P5_pilots_recovered": True,
        "combined_positive_increment_beyond_all_five_original_P5_pilots": 56_936,
        "remaining_explicit_BBBB_pilot_has_788_facets": True,
    }
    payload["pilot_recovery"][
        "eligible_ABBB_B3_pilot_count"
    ] = payload["pilot_recovery"].pop("eligible_B_0_1_2_pilot_count")
    payload["pilot_recovery"][
        "ABBB_positive_increment_beyond_original_B3_pilot"
    ] = payload["pilot_recovery"].pop(
        "project_positive_increment_over_the_five_path_pilot_within_this_layer"
    )
    payload["resumable_full_census_boundary"] = {
        "all_ten_binary_reversal_shards": list(frozen.ALL_PATTERNS),
        "globally_completed_shards": completed_patterns,
        "remaining_shards": ["BABB", "BBBB"],
        "remaining_orbit_count": remaining_count,
        "next_single_shard_CLI": (
            "frozen P5 engine with --pattern BABB --output PATH"
        ),
        "final_single_shard_CLI": (
            "frozen P5 engine with --pattern BBBB --output PATH"
        ),
        "shards_are_pairwise_disjoint_because_the_seed_action_preserves_A_and_B": True,
        "ABBB_vs_frozen_B_le_2_exact_merge_closed_by_facet_count": True,
        "BABB_and_BBBB_matrix_censuses_remain_open": True,
    }
    payload["lineage_and_claim_boundary"] = {
        "frozen_B_le_2_complete_type_count": 20_693,
        "complete_ABBB_type_count": 36_247,
        "combined_exact_completed_P5_type_count": combined_count,
        "ABBB_positive_increment_beyond_original_B3_pilot": 36_246,
        "combined_positive_increment_beyond_all_five_original_P5_pilots": 56_936,
        "complete_P5_census_claimed": False,
        "remaining_BABB_BBBB_orbit_denominator": remaining_count,
        "degree_three_fork_complete_census_claimed": False,
        "all_five_vertex_shapes_complete_census_claimed": False,
        "cross_star_pure_exact_nonisomorphism_claimed": False,
        "worldwide_priority": "open",
        "weak_carrier_not_full_nef": True,
    }
    payload["checks"].pop("production_scope_is_exactly_B_at_most_2")
    payload["checks"].pop("full_P5_census_not_claimed_for_production_scope")
    payload["checks"] = {**payload["checks"], **wrapper_checks}
    if not all(payload["checks"].values()):
        raise ArithmeticError("the completed ABBB payload is not wholly green")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.probe is not None:
        _frozen_gate()
        frozen._probe(PATTERNS, arguments.probe)
        return
    payload = build(verbose=True)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"archive is missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact ABBB rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(f"ABBB_types={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(
        "combined_completed_P5_types="
        f"{payload['combined_exact_P5_progress']['combined_completed_type_count']}"
    )
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
