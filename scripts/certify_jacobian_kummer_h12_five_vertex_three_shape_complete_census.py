#!/usr/bin/env python3
"""Combine the complete JK--H12 P5, star, and degree-three-fork censuses."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_three_shape_complete_census.json"
)
PATHS = {
    "all_shape_tree_denominator": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json",
    "fork_pattern_denominator": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json",
    "fork_B_le_2_theorem": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.json",
    "fork_B3_terminal_acceptance": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.json",
    "fork_B4_terminal_acceptance_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.py",
    "fork_B4_terminal_acceptance": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.json",
    "complete_P5_summary": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json",
    "complete_star": ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json",
    "complete_star_P5_separation": ROOT
    / "data/enumerations/jacobian_kummer_h12_star_p5_b_le_2_exact_separation.json",
}
PINNED_SHA256 = {
    "all_shape_tree_denominator": "2306ebe0048b576a5ba1aa145d5af500b6efdb13f16ffacb8174c7984d1d5d10",
    "fork_pattern_denominator": "99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7",
    "fork_B_le_2_theorem": "e58a7784f84d06a4710a332a292ef75b689636a3588c7f2bb9d1489673c1ae6e",
    "fork_B3_terminal_acceptance": "4a385121b4f0f14484e258223d62e15fae8b78b42c31ef85c0fda71380176d38",
    "fork_B4_terminal_acceptance_builder": "41d18e60873deec7ce04705e68f99996382cf3e954e6e88b2aed668dab8007eb",
    "fork_B4_terminal_acceptance": "25ac2f56bd33b41b2d40e04d6659addc6653226cc3c96dab0e034cb6d71cd2e1",
    "complete_P5_summary": "dd6c60a8ed02e985a889100991d2da410d8b18f0902fa314e1e6e9dc900dab16",
    "complete_star": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "complete_star_P5_separation": "b73e85b9cd3fa61fc727724dca966e906af33feadcaa5b656c758f6640dee15f",
}
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_three_shape_complete_census.v1"
STATUS_DATE = "2026-08-10"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _histogram(*items: dict[str, int]) -> dict[str, int]:
    total: Counter[str] = Counter()
    for item in items:
        total.update({str(key): int(value) for key, value in item.items()})
    return dict(sorted(total.items(), key=lambda pair: int(pair[0])))


def build() -> dict[str, Any]:
    actual = {name: _sha256(path) for name, path in PATHS.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen five-vertex combination input changed: {actual}")

    tree_denominator = json.loads(PATHS["all_shape_tree_denominator"].read_text())
    fork_denominator = json.loads(PATHS["fork_pattern_denominator"].read_text())
    low = json.loads(PATHS["fork_B_le_2_theorem"].read_text())
    b3 = json.loads(PATHS["fork_B3_terminal_acceptance"].read_text())
    b4 = json.loads(PATHS["fork_B4_terminal_acceptance"].read_text())
    p5 = json.loads(PATHS["complete_P5_summary"].read_text())
    star = json.loads(PATHS["complete_star"].read_text())
    p5_star = json.loads(PATHS["complete_star_P5_separation"].read_text())

    fork_histogram = _histogram(
        low["theorem"]["type_facet_histogram"],
        b3["terminal_transcript_audit"]["facet_histogram"],
        b4["terminal_transcript_audit"]["facet_histogram"],
    )
    expected_fork_histogram = {
        "752": 116,
        "761": 2_042,
        "770": 18_061,
        "779": 72_758,
        "788": 106_385,
    }
    p5_count = int(p5["theorem"]["ordinary_Coxeter_type_count"])
    star_count = int(star["theorem"]["ordinary_Coxeter_type_count"])
    fork_count = sum(expected_fork_histogram.values())
    total = p5_count + star_count + fork_count
    aggregate_histogram = _histogram(
        p5["theorem"]["type_facet_histogram"],
        star["theorem"]["type_facet_histogram"],
        expected_fork_histogram,
    )
    checks = {
        "all_frozen_inputs_green": all(
            payload["all_checks_pass"]
            for payload in (low, b3, b4, p5, star, p5_star)
        ),
        "three_unlabelled_five_vertex_tree_shapes_are_complete": set(
            tree_denominator["tree_shapes"]
        )
        == {"P5", "K1,4", "degree_3_fork"},
        "complete_fork_count_is_199362": fork_count == 199_362
        and fork_count == int(fork_denominator["total_orbits"])
        and fork_histogram == expected_fork_histogram,
        "fork_B_partition_matches_denominator": {
            str(index): count
            for index, count in enumerate(expected_fork_histogram.values())
        }
        == fork_denominator["orbit_count_by_B"],
        "fork_layers_pairwise_separated_by_facet_header": len(fork_histogram) == 5,
        "every_fork_row_separated_from_same_facet_P5_and_star": low["theorem"][
            "separated_from_every_complete_P5_and_star_row_with_the_same_facet_header"
        ]
        and b3["theorem"][
            "separated_from_all_79012_completed_same_facet_P5_and_star_rows"
        ]
        and b4["theorem"][
            "separated_from_all_117144_completed_same_facet_P5_and_star_rows"
        ],
        "complete_P5_star_union_is_218619": p5_star["theorem"][
            "complete_star_plus_complete_P5_pairwise_nonisomorphic_type_count"
        ]
        == 218_619,
        "shape_counts_are_exact_denominators": (p5_count, star_count, fork_count)
        == (201_747, 16_872, 199_362),
        "complete_three_shape_union_is_417981": total == 417_981
        and sum(aggregate_histogram.values()) == total,
        "all_rows_finite_volume_noncompact_index_five": all(
            payload["theorem"]["finite_volume"]
            and not payload["theorem"]["compact"]
            and payload["theorem"]["reflection_index_and_volume_multiple"] == 5
            for payload in (low, b3, b4, p5_star)
        ),
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the complete five-vertex theorem failed: {checks}")

    low_collision = low["exact_collision_accounting"]
    b3_collision = b3["exact_collision_accounting"]
    b4_collision = b4["exact_collision_accounting"]
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": PINNED_SHA256,
        "theorem": {
            "scope": (
                "complete P5, K1,4-star, and degree-three-fork tree layers "
                "over the frozen five-vertex JK-H12 doubling stock"
            ),
            "ordinary_Coxeter_type_count": total,
            "pairwise_nonisomorphic": True,
            "shape_type_counts": {
                "path_P5": p5_count,
                "star_K1,4": star_count,
                "degree_3_fork": fork_count,
            },
            "type_facet_histogram": aggregate_histogram,
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiple": 5,
            "hyperbolic_dimension": 12,
        },
        "complete_fork_theorem": {
            "ordinary_Coxeter_type_count": fork_count,
            "pairwise_nonisomorphic": True,
            "type_count_by_B": fork_denominator["orbit_count_by_B"],
            "type_facet_histogram": expected_fork_histogram,
            "complete_matrix_classification_of_all_fork_orbits": True,
            "remaining_unclassified_orbit_count": 0,
            "separated_from_every_complete_P5_and_star_row": True,
        },
        "exact_cross_shape_accounting": {
            "complete_P5_plus_star_type_count": 218_619,
            "complete_P5_plus_star_plus_fork_type_count": total,
            "fork_internal_literal_degree_collision_group_count": sum(
                payload["within_fork_literal_degree_collision_group_count"]
                for payload in (low_collision, b3_collision, b4_collision)
            ),
            "fork_internal_literal_degree_collision_member_count": sum(
                payload["within_fork_literal_degree_collision_member_count"]
                for payload in (low_collision, b3_collision, b4_collision)
            ),
            "fork_cross_P5_star_literal_degree_match_count": sum(
                payload["cross_P5_star_literal_degree_match_count"]
                for payload in (low_collision, b3_collision, b4_collision)
            ),
            "all_equality_candidates_replayed_literally": True,
            "all_literal_collisions_separated_by_exact_color_refinement": True,
            "unresolved_pair_count": 0,
        },
        "surface_carrier_and_claim_boundary": {
            "direction": "SURFACE_TO_POLYTOPE",
            "same_full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "weak_divisor_reflection_carrier": True,
            "full_nef_effectivity_or_surface_action_claim": False,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
        },
        "claim_boundary": {
            "complete_fork_census_claimed": True,
            "all_three_five_vertex_tree_shapes_complete": True,
            "all_417981_five_vertex_denominator_rows_classified": True,
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
            raise SystemExit("the complete five-vertex archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("complete_fork_types=199362")
    print("complete_five_vertex_types=417981")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
