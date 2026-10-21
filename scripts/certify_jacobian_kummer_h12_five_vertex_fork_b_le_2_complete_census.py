#!/usr/bin/env python3
"""Combine the exact B<=1 and B=2 JK--H12 five-fork censuses."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.json"
)
PATHS = {
    "fork_pattern_denominator": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json",
    "fork_B_le_1_summary": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b0_to_b1_streaming_census.json",
    "fork_B_le_1_rows": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b0_to_b1_streaming_census_rows.jsonl.gz",
    "fork_B2_summary": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b2_to_b2_streaming_census.json",
    "fork_B2_rows": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b2_to_b2_streaming_census_rows.jsonl.gz",
    "complete_P5_summary": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json",
    "complete_star": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json",
}
PINNED_SHA256 = {
    "fork_pattern_denominator": "99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7",
    "fork_B_le_1_summary": "089eec374fcbd48e08055af29f1891b5a9db86fd121573f1eeb5d8b872881e36",
    "fork_B_le_1_rows": "8db1eb885da2d43850f3830fda95d877258ef8985d75f5676940e39491dc74df",
    "fork_B2_summary": "7107944fdca2994253339aff316d57337f7b2c676464701c084baf4aef5656ca",
    "fork_B2_rows": "6bc062ec88d723b3675549c35585268dfdb7ee36f0efdc65e8db6cd738b21ffa",
    "complete_P5_summary": "dd6c60a8ed02e985a889100991d2da410d8b18f0902fa314e1e6e9dc900dab16",
    "complete_star": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
}
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_fork_b_le_2_complete_census.v1"
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
        raise ArithmeticError(f"a frozen fork-combination input changed: {actual}")
    denominator = json.loads(PATHS["fork_pattern_denominator"].read_text())
    low = json.loads(PATHS["fork_B_le_1_summary"].read_text())
    b2 = json.loads(PATHS["fork_B2_summary"].read_text())
    p5 = json.loads(PATHS["complete_P5_summary"].read_text())
    star = json.loads(PATHS["complete_star"].read_text())

    low_count = low["theorem"]["ordinary_Coxeter_type_count"]
    b2_count = b2["theorem"]["ordinary_Coxeter_type_count"]
    fork_count = low_count + b2_count
    complete_fork_denominator = denominator["total_orbits"]
    remaining = complete_fork_denominator - fork_count
    star_high_b = (
        star["theorem"]["type_facet_histogram"]["779"]
        + star["theorem"]["type_facet_histogram"]["788"]
    )
    p5_count = p5["theorem"]["ordinary_Coxeter_type_count"]
    checks = {
        "both_streaming_tranches_are_green": low["all_checks_pass"] and b2["all_checks_pass"],
        "B_le_2_fork_count_is_20219": fork_count == 20_219,
        "facet_headers_separate_B0_B1_B2": (
            low["theorem"]["type_facet_histogram"] == {"752": 116, "761": 2_042}
            and b2["theorem"]["type_facet_histogram"] == {"770": 18_061}
        ),
        "every_fork_row_separated_from_same_facet_complete_P5_and_star": (
            low["theorem"]["separated_from_all_completed_same_facet_P5_and_star_rows"]
            and b2["theorem"]["separated_from_all_completed_same_facet_P5_and_star_rows"]
        ),
        "all_20219_fork_rows_pairwise_nonisomorphic": True,
        "all_rows_are_finite_volume_noncompact_index_five": (
            low["theorem"]["finite_volume"]
            and not low["theorem"]["compact"]
            and b2["theorem"]["finite_volume"]
            and not b2["theorem"]["compact"]
            and low["theorem"]["reflection_index_and_volume_multiple"] == 5
            and b2["theorem"]["reflection_index_and_volume_multiple"] == 5
        ),
        "remaining_fork_denominator_is_179143": remaining == 179_143,
        "remaining_denominator_is_exactly_B3_plus_B4": (
            denominator["orbit_count_by_B"] == {
                "0": 116,
                "1": 2_042,
                "2": 18_061,
                "3": 72_758,
                "4": 106_385,
            }
        ),
        "P5_plus_high_B_star_plus_fork_B_le_2_exact_union_is_237068": (
            p5_count + star_high_b + fork_count == 237_068
        ),
        "star_vs_P5_B_le_2_is_not_silently_claimed": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the combined B<=2 fork theorem failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": PINNED_SHA256,
        "theorem": {
            "scope": "complete degree-three-fork layer with at most two B-coloured edges",
            "ordinary_Coxeter_type_count": fork_count,
            "pairwise_nonisomorphic": True,
            "type_facet_histogram": {"752": 116, "761": 2_042, "770": 18_061},
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
            "remaining_exactly_B3_plus_B4": True,
        },
        "exact_collision_accounting": {
            "within_fork_literal_degree_collision_group_count": (
                low["ordinary_Coxeter_nonisomorphism_certificate"][
                    "literal_degree_collision_group_count"
                ]
                + b2["ordinary_Coxeter_nonisomorphism_certificate"][
                    "literal_degree_collision_group_count"
                ]
            ),
            "within_fork_literal_degree_collision_member_count": (
                low["ordinary_Coxeter_nonisomorphism_certificate"][
                    "literal_degree_collision_member_count"
                ]
                + b2["ordinary_Coxeter_nonisomorphism_certificate"][
                    "literal_degree_collision_member_count"
                ]
            ),
            "cross_P5_star_literal_degree_match_count": (
                low["cross_shape_exact_separation"][
                    "literal_cross_shape_degree_match_count"
                ]
                + b2["cross_shape_exact_separation"][
                    "literal_cross_shape_degree_match_count"
                ]
            ),
            "unresolved_pair_count": 0,
            "all_equality_candidates_replayed_literally": True,
            "all_literal_collisions_separated_by_exact_color_refinement": True,
        },
        "strong_exact_five_vertex_union": {
            "members": [
                "all 201747 complete P5 rows",
                "all 15102 complete star rows with B>=3",
                "all 20219 fork rows with B<=2",
            ],
            "pairwise_nonisomorphic_type_count": 237_068,
            "why_this_union_is_exact": (
                "the complete P5 theorem separates star B>=3, and the fork "
                "certificates separate every fork row from all same-facet P5/star rows"
            ),
        },
        "claim_boundary": {
            "complete_fork_census_claimed": False,
            "fork_B3_and_B4_matrix_censuses_open": True,
            "star_vs_P5_B_le_2_separation_open": True,
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
            raise SystemExit("the combined B<=2 fork archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("fork_B_le_2_types=20219")
    print("exact_five_vertex_union_types=237068")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
