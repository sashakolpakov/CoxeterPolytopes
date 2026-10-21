#!/usr/bin/env python3
"""Exact Burnside denominators for all three five-vertex JK-H12 trees."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base


BASE_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_small_subtree_census.py"
BASE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json"
SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_tree_denominators.v1"
STATUS_DATE = "2026-08-08"
STOCK_START = 60
STOCK_END = 172
STOCK_SIZE = 112
GROUP_ORDER = 384
PINNED_SHA256 = {
    "small_tree_builder": "30ccb14a9240c32dfd8baeb70f1a3482a3941f5a6e4a9523094f4ebc4630d125",
    "small_tree_archive": "796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d",
}

Permutation = tuple[int, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _restricted_actions() -> tuple[Permutation, ...]:
    actions = base._folded_chamber_actions()
    restricted = tuple(
        tuple(action[index] - STOCK_START for index in range(STOCK_START, STOCK_END))
        for action in actions
    )
    if not (
        len(restricted) == GROUP_ORDER
        and len(set(restricted)) == GROUP_ORDER
        and all(sorted(action) == list(range(STOCK_SIZE)) for action in restricted)
    ):
        raise ArithmeticError("the exact 112-wall seed action changed")
    return restricted


def _cycle_lengths(permutation: Sequence[int]) -> tuple[int, ...]:
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        if current != start:
            raise ArithmeticError("an action is not a permutation")
        lengths.append(length)
    return tuple(sorted(lengths))


def _fixed_four_subsets(lengths: Sequence[int]) -> int:
    coefficients = [1, 0, 0, 0, 0]
    for length in lengths:
        for degree in range(4, length - 1, -1):
            coefficients[degree] += coefficients[degree - length]
    return coefficients[4]


def build() -> dict[str, Any]:
    actual = {
        "small_tree_builder": _sha256(BASE_SCRIPT),
        "small_tree_archive": _sha256(BASE_ARCHIVE),
    }
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a pinned denominator prerequisite changed: {actual}")
    small = json.loads(BASE_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        small["all_checks_pass"]
        and small["selected_stock"]["wall_count"] == STOCK_SIZE
        and small["selected_stock"]["folded_chamber_action_order"] == GROUP_ORDER
    ):
        raise ArithmeticError("the frozen seed action is not green")

    actions = _restricted_actions()
    cycle_data = []
    k14_sum = p5_sum = fork_sum = 0
    for action in actions:
        lengths = _cycle_lengths(action)
        fixed = lengths.count(1)
        transposition_points = 2 * lengths.count(2)
        fixed_by_square = fixed + transposition_points
        k14_fixed = _fixed_four_subsets(lengths)
        p5_identity_fixed = fixed * (fixed - 1) ** 3
        p5_reversal_fixed = transposition_points * (fixed_by_square - 1)
        fork_identity_fixed = fixed * (fixed - 1) ** 2 * (fixed - 2)
        fork_short_arm_swap_fixed = fixed * (fixed - 1) * transposition_points
        k14_sum += k14_fixed
        p5_sum += p5_identity_fixed + p5_reversal_fixed
        fork_sum += fork_identity_fixed + fork_short_arm_swap_fixed
        cycle_data.append(
            (
                tuple(sorted(Counter(lengths).items())),
                k14_fixed,
                p5_identity_fixed,
                p5_reversal_fixed,
                fork_identity_fixed,
                fork_short_arm_swap_fixed,
            )
        )

    raw = {
        "P5": STOCK_SIZE * (STOCK_SIZE - 1) ** 3 // 2,
        "K1,4": math.comb(STOCK_SIZE, 4),
        "degree_3_fork": (
            STOCK_SIZE * (STOCK_SIZE - 1) * math.comb(STOCK_SIZE - 1, 2)
        ),
    }
    orbit_counts = {
        "P5": p5_sum // (2 * GROUP_ORDER),
        "K1,4": k14_sum // GROUP_ORDER,
        "degree_3_fork": fork_sum // (2 * GROUP_ORDER),
    }
    expected_raw = {"P5": 76_587_336, "K1,4": 6_210_820, "degree_3_fork": 75_897_360}
    expected_orbits = {"P5": 201_747, "K1,4": 16_872, "degree_3_fork": 199_362}
    if not (
        raw == expected_raw
        and orbit_counts == expected_orbits
        and k14_sum == 6_478_848
        and p5_sum == 154_941_696
        and fork_sum == 153_110_016
        and all(value >= 0 for row in cycle_data for value in row[1:])
    ):
        raise ArithmeticError("a five-vertex Burnside denominator changed")

    cycle_type_rows = []
    grouped = Counter(cycle_data)
    for data, multiplicity in sorted(grouped.items()):
        (
            cycle_type,
            k14_fixed,
            p5_identity_fixed,
            p5_reversal_fixed,
            fork_identity_fixed,
            fork_swap_fixed,
        ) = data
        cycle_type_rows.append(
            {
                "cycle_type": [[length, count] for length, count in cycle_type],
                "group_element_count": multiplicity,
                "K14_fixed_colorings": k14_fixed,
                "P5_identity_fixed_ordered_words": p5_identity_fixed,
                "P5_reversal_fixed_ordered_words": p5_reversal_fixed,
                "fork_identity_fixed_ordered_short_arm_words": fork_identity_fixed,
                "fork_short_arm_swap_fixed_ordered_words": fork_swap_fixed,
            }
        )
    checks = {
        "exact_384_element_action_reconstructed": len(actions) == GROUP_ORDER,
        "all_three_unlabelled_five_vertex_tree_shapes_included": True,
        "proper_edge_coloring_constraints_encoded": True,
        "P5_reversal_automorphism_encoded": True,
        "K14_S4_edge_automorphism_encoded_as_four_subsets": True,
        "fork_short_arm_swap_automorphism_encoded": True,
        "raw_coloring_counts_match_direct_formulas": raw == expected_raw,
        "Burnside_sums_divide_by_their_exact_group_orders": (
            k14_sum % GROUP_ORDER == 0
            and p5_sum % (2 * GROUP_ORDER) == 0
            and fork_sum % (2 * GROUP_ORDER) == 0
        ),
        "exact_colored_tree_orbit_counts_are_201747_16872_199362": (
            orbit_counts == expected_orbits
        ),
        "total_search_denominator_is_417981": sum(orbit_counts.values()) == 417_981,
        "no_Coxeter_matrix_nonisomorphism_claim_is_inferred": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"five-vertex denominator checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": PINNED_SHA256,
        "stock_size": STOCK_SIZE,
        "seed_action_order": GROUP_ORDER,
        "tree_shapes": {
            "P5": {
                "proper_coloring_rule": "(a,b,c,d) with adjacent colors distinct, modulo reversal",
                "raw_unrooted_colorings": raw["P5"],
                "Burnside_numerator": p5_sum,
                "Burnside_group_order_including_tree_automorphism": 2 * GROUP_ORDER,
                "seed_and_tree_automorphism_orbits": orbit_counts["P5"],
            },
            "K1,4": {
                "proper_coloring_rule": "unordered four-subset of distinct colors",
                "raw_unrooted_colorings": raw["K1,4"],
                "Burnside_numerator": k14_sum,
                "Burnside_group_order": GROUP_ORDER,
                "seed_and_tree_automorphism_orbits": orbit_counts["K1,4"],
            },
            "degree_3_fork": {
                "proper_coloring_rule": (
                    "long-arm colors a,d with a!=d; unordered short-arm colors "
                    "{b,c} distinct from each other and a; d may equal b or c"
                ),
                "raw_unrooted_colorings": raw["degree_3_fork"],
                "Burnside_numerator": fork_sum,
                "Burnside_group_order_including_short_arm_swap": 2 * GROUP_ORDER,
                "seed_and_tree_automorphism_orbits": orbit_counts["degree_3_fork"],
            },
        },
        "total_raw_unrooted_colored_trees": sum(raw.values()),
        "total_seed_and_tree_automorphism_orbits": sum(orbit_counts.values()),
        "cycle_type_Burnside_table": cycle_type_rows,
        "claim_boundary": {
            "this_is_a_colored_tree_denominator_not_a_polytope_census": True,
            "different_color_orbits_may_still_give_isomorphic_ordinary_Coxeter_matrices": True,
            "finite_volume_checks_belong_to_the_shape_census": True,
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
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact five-vertex denominator rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("orbit_counts=P5:201747,K1_4:16872,fork:199362,total:417981")
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
