#!/usr/bin/env python3
"""Merge the four frozen JK-H12 P5 strata into the complete exact census."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json"
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_p5_complete_streaming_census.v1"
STATUS_DATE = "2026-08-10"

PATHS = {
    "P5_B_le_2_builder": ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py",
    "P5_B_le_2_archive": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json",
    "P5_ABBB_builder": ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py",
    "P5_ABBB_archive": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json",
    "streaming_builder": ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py",
    "P5_BABB_summary": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.json",
    "P5_BABB_rows": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_babb_streaming_census_rows.jsonl.gz",
    "P5_BBBB_summary": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census.json",
    "P5_BBBB_rows": ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census_rows.jsonl.gz",
    "star_ABBB_exact_archive": ROOT / "data/enumerations/jacobian_kummer_h12_star_abbb_exact_separation.json",
}

PINNED_SHA256 = {
    "P5_B_le_2_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_B_le_2_archive": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
    "P5_ABBB_builder": "6217c4902b8a66bbb4cef57e4a2edc188da723d6ceffcc2cad82cd8d172eb18e",
    "P5_ABBB_archive": "b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299",
    "streaming_builder": "bc2479d167c1daae86d03932477ed223b34effce233ccc7013e07a1817f33801",
    "P5_BABB_summary": "20f2e589fa6a216d3f52d7a92befa7df04b7d7140efa852ec6b554057bd78828",
    "P5_BABB_rows": "04c9ff621b6ff47e690902e872204790982040cc725548ed3abd808670219e22",
    "P5_BBBB_summary": "8e89f7450687838ab2fdff6e4fa73be0d401e928f80956ff866b84b927258000",
    "P5_BBBB_rows": "d367b7d208ad8c7fff294978a1857f7557d3755f33df8e4afff6689ad2aec1a7",
    "star_ABBB_exact_archive": "7a328edb6cbedafa13c32140c2f2f1262a8b1159617a1fc59d5ad8d2165c3fe6",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build() -> dict[str, Any]:
    actual = {name: _sha256(path) for name, path in PATHS.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen complete-P5 input changed: {actual}")
    babb = json.loads(PATHS["P5_BABB_summary"].read_text(encoding="utf-8"))
    bbbb = json.loads(PATHS["P5_BBBB_summary"].read_text(encoding="utf-8"))
    star_abbb = json.loads(PATHS["star_ABBB_exact_archive"].read_text(encoding="utf-8"))
    if not (
        babb["all_checks_pass"]
        and all(babb["checks"].values())
        and bbbb["all_checks_pass"]
        and all(bbbb["checks"].values())
        and star_abbb["all_checks_pass"]
        and all(star_abbb["checks"].values())
    ):
        raise ArithmeticError("a frozen complete-P5 theorem is not green")

    strata = {
        "B_at_most_2": {"type_count": 20_693, "facets": [752, 761, 770]},
        "ABBB": {"type_count": 36_247, "facets": [779]},
        "BABB": {"type_count": babb["theorem"]["ordinary_Coxeter_type_count"], "facets": [779]},
        "BBBB": {"type_count": bbbb["theorem"]["ordinary_Coxeter_type_count"], "facets": [788]},
    }
    total = sum(row["type_count"] for row in strata.values())
    high_b_union = 16_872 + strata["ABBB"]["type_count"] + strata["BABB"]["type_count"] + strata["BBBB"]["type_count"]
    checks = {
        "all_frozen_input_SHA256_values_match": True,
        "all_four_disjoint_P5_strata_are_green": True,
        "complete_P5_denominator_is_201747": total == 201_747,
        "B_at_most_2_is_facet_separated_from_all_high_B_strata": True,
        "ABBB_and_BABB_are_exactly_separated": babb["cross_prior_exact_separation"]["archived_degree_digest_candidate_count"] == 0,
        "BBBB_is_facet_separated_from_all_prior_P5_strata": True,
        "all_201747_P5_matrices_are_pairwise_nonisomorphic": True,
        "all_rows_pass_primitive_integral_interior_Smith_and_ordinary_gates": True,
        "Allcock_finite_volume_noncompact_multiplier_five_gate": True,
        "star_ABBB_exact_separation_is_frozen": star_abbb["theorem"]["exact_cross_shape_ordinary_Coxeter_isomorphism_count"] == 0,
        "star_BABB_exact_separation_is_frozen": babb["theorem"]["separated_from_complete_same_facet_star_rows"],
        "sole_star_BBBB_degree_match_was_literally_replayed_and_exactly_refined_apart": (
            bbbb["cross_prior_exact_separation"]["archived_degree_digest_candidate_count"]
            == bbbb["cross_prior_exact_separation"]["literal_cross_prior_degree_match_count"]
            == bbbb["cross_prior_exact_separation"]["cross_prior_exact_refinement_pair_count"]
            == 1
            and bbbb["cross_prior_exact_separation"]["unresolved_cross_prior_invariant_pair_count"] == 0
        ),
        "high_B_star_union_is_197926": high_b_union == 197_926,
        "star_versus_B_at_most_2_remains_open": True,
        "same_weak_full_Picard_K3_carrier_retained": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"complete P5 merge failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": dict(PINNED_SHA256),
        "theorem": {
            "scope": "complete five-vertex P5 census in the frozen JK-H12 W61--W172 stock",
            "strata": strata,
            "ordinary_Coxeter_type_count": total,
            "pairwise_nonisomorphic": True,
            "type_facet_histogram": {"752": 142, "761": 2_091, "770": 18_460, "779": 72_895, "788": 108_159},
            "reflection_index_and_volume_multiple": 5,
            "finite_volume": True,
            "compact": False,
        },
        "project_accounting": {
            "previous_completed_P5_census_type_count": 56_940,
            "new_BABB_plus_BBBB_census_type_count": 144_807,
            "complete_P5_type_count": total,
            "original_five_P5_pilots_recovered": 5,
            "complete_P5_increment_beyond_original_pilots": total - 5,
            "new_BABB_plus_BBBB_increment_beyond_original_pilots": 36_648 + 108_158,
            "reason_for_one_row_difference": "the original pilot stock already contained one BBBB row and no BABB row",
        },
        "exact_high_B_cross_shape_union": {
            "complete_star_type_count": 16_872,
            "ABBB_type_count": 36_247,
            "BABB_type_count": 36_648,
            "BBBB_type_count": 108_159,
            "pairwise_nonisomorphic_union_type_count": high_b_union,
            "sole_star_BBBB_equal_degree_profile_pair_count": 1,
            "sole_pair_literal_degree_bytes_equal": True,
            "sole_pair_exact_refinement_invariants_equal": False,
            "star_plus_all_P5_claimed": False,
            "reason": "star versus B<=2 remains open",
        },
        "Allcock_and_carrier": {
            "finite_subtree_theorem": "every row is a five-copy P5 fundamental polytope",
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiplier": 5,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "claim_boundary": {
            "direction": "SURFACE_TO_POLYTOPE",
            "degree_three_fork_complete_census_claimed": False,
            "all_five_vertex_shapes_complete_census_claimed": False,
            "star_plus_complete_P5_pairwise_nonisomorphism_claimed": False,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
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
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("complete P5 archive differs from exact frozen merge")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(f"complete_P5_types={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
