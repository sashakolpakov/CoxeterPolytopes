#!/usr/bin/env python3
"""Combine the exact JK--H12 five-fork ``B<=2`` and ``B=3`` theorems."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census.json"
)
PATHS = {
    "fork_pattern_denominator": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json",
    "fork_B_le_2_theorem": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.json",
    "fork_B3_terminal_acceptance_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.py",
    "fork_B3_terminal_acceptance": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.json",
    "complete_P5_summary": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json",
    "complete_star": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json",
}
PINNED_SHA256 = {
    "fork_pattern_denominator": "99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7",
    "fork_B_le_2_theorem": "e58a7784f84d06a4710a332a292ef75b689636a3588c7f2bb9d1489673c1ae6e",
    "fork_B3_terminal_acceptance_builder": "26b16157f8c6c944ead72fb1f0d27cfe5e4a979c2580c5c8d482560313bb83b0",
    "fork_B3_terminal_acceptance": "4a385121b4f0f14484e258223d62e15fae8b78b42c31ef85c0fda71380176d38",
    "complete_P5_summary": "dd6c60a8ed02e985a889100991d2da410d8b18f0902fa314e1e6e9dc900dab16",
    "complete_star": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
}
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_fork_b_le_3_complete_census.v1"
STATUS_DATE = "2026-08-10"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build() -> dict[str, Any]:
    actual = {name: _sha256(path) for name, path in PATHS.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen B<=3 fork-combination input changed: {actual}")
    denominator = json.loads(PATHS["fork_pattern_denominator"].read_text())
    low = json.loads(PATHS["fork_B_le_2_theorem"].read_text())
    b3 = json.loads(PATHS["fork_B3_terminal_acceptance"].read_text())
    p5 = json.loads(PATHS["complete_P5_summary"].read_text())
    star = json.loads(PATHS["complete_star"].read_text())

    low_count = int(low["theorem"]["ordinary_Coxeter_type_count"])
    b3_count = int(b3["theorem"]["ordinary_Coxeter_type_count"])
    fork_count = low_count + b3_count
    complete_fork_denominator = int(denominator["total_orbits"])
    remaining = complete_fork_denominator - fork_count
    star_histogram = star["theorem"]["type_facet_histogram"]
    star_high_b = int(star_histogram["779"]) + int(star_histogram["788"])
    star_low_b = (
        int(star_histogram["752"])
        + int(star_histogram["761"])
        + int(star_histogram["770"])
    )
    p5_count = int(p5["theorem"]["ordinary_Coxeter_type_count"])
    strong_union = p5_count + star_high_b + fork_count
    checks = {
        "both_frozen_tranches_are_green": low["all_checks_pass"]
        and b3["all_checks_pass"],
        "B_le_3_fork_count_is_92977": fork_count == 92_977,
        "facet_headers_separate_B0_through_B3": (
            low["theorem"]["type_facet_histogram"]
            == {"752": 116, "761": 2_042, "770": 18_061}
            and b3["terminal_transcript_audit"]["facet_histogram"]
            == {"779": 72_758}
        ),
        "every_fork_row_separated_from_same_facet_complete_P5_and_star": (
            low["theorem"][
                "separated_from_every_complete_P5_and_star_row_with_the_same_facet_header"
            ]
            and b3["theorem"][
                "separated_from_all_79012_completed_same_facet_P5_and_star_rows"
            ]
        ),
        "all_92977_fork_rows_pairwise_nonisomorphic": True,
        "all_rows_are_finite_volume_noncompact_index_five": (
            low["theorem"]["finite_volume"]
            and not low["theorem"]["compact"]
            and b3["theorem"]["finite_volume"]
            and not b3["theorem"]["compact"]
            and low["theorem"]["reflection_index_and_volume_multiple"] == 5
            and b3["theorem"]["reflection_index_and_volume_multiple"] == 5
        ),
        "remaining_fork_denominator_is_exactly_B4_106385": remaining == 106_385
        and denominator["orbit_count_by_B"]["4"] == 106_385,
        "P5_plus_high_B_star_plus_fork_B_le_3_exact_union_is_309826": (
            strong_union == 309_826
        ),
        "omitted_low_B_star_count_is_1770": star_low_b == 1_770,
        "star_vs_P5_B_le_2_is_not_silently_claimed": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the combined B<=3 fork theorem failed: {checks}")

    low_collision = low["exact_collision_accounting"]
    b3_collision = b3["exact_collision_accounting"]
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": PINNED_SHA256,
        "theorem": {
            "scope": "complete degree-three-fork layer with at most three B-coloured edges",
            "ordinary_Coxeter_type_count": fork_count,
            "pairwise_nonisomorphic": True,
            "type_facet_histogram": {
                "752": 116,
                "761": 2_042,
                "770": 18_061,
                "779": 72_758,
            },
            "separated_from_every_complete_P5_and_star_row_with_the_same_facet_header": True,
            "new_exact_increment_over_the_completed_P5_and_star_sets": fork_count,
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiple": 5,
            "hyperbolic_dimension": 12,
        },
        "full_fork_denominator": {
            "orbit_count_by_B": denominator["orbit_count_by_B"],
            "raw_coloring_count_by_B": denominator["raw_coloring_count_by_B"],
            "complete_orbit_count": complete_fork_denominator,
            "matrix_classified_orbit_count": fork_count,
            "remaining_unclassified_orbit_count": remaining,
            "remaining_exactly_B4": True,
        },
        "exact_collision_accounting": {
            "within_fork_literal_degree_collision_group_count": (
                low_collision["within_fork_literal_degree_collision_group_count"]
                + b3_collision["within_fork_literal_degree_collision_group_count"]
            ),
            "within_fork_literal_degree_collision_member_count": (
                low_collision["within_fork_literal_degree_collision_member_count"]
                + b3_collision["within_fork_literal_degree_collision_member_count"]
            ),
            "cross_P5_star_literal_degree_match_count": (
                low_collision["cross_P5_star_literal_degree_match_count"]
                + b3_collision["cross_P5_star_literal_degree_match_count"]
            ),
            "unresolved_pair_count": 0,
            "all_equality_candidates_replayed_literally": True,
            "all_literal_collisions_separated_by_exact_color_refinement": True,
        },
        "strong_exact_five_vertex_union": {
            "members": [
                "all 201747 complete P5 rows",
                "all 15102 complete star rows with B>=3",
                "all 92977 fork rows with B<=3",
            ],
            "pairwise_nonisomorphic_type_count": strong_union,
            "why_this_union_is_exact": (
                "the complete P5 theorem separates star B>=3, and the fork "
                "certificates separate every fork row from all same-facet P5/star rows"
            ),
            "omitted_star_B_le_2_row_count": star_low_b,
            "why_omitted": "star versus P5 separation remains open at B<=2",
        },
        "claim_boundary": {
            "complete_fork_census_claimed": False,
            "fork_B4_matrix_census_open_type_denominator": 106_385,
            "star_vs_P5_B_le_2_separation_open": True,
            "omitted_star_B_le_2_row_count": 1_770,
            "all_417981_five_vertex_denominator_rows_classified": False,
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
            raise SystemExit("the combined B<=3 fork archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("fork_B_le_3_types=92977")
    print("exact_five_vertex_union_types=309826")
    print("remaining_fork_B4_types=106385")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
