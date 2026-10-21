#!/usr/bin/env python3
"""Lossless full-source generator-lex-slice census for Record 14 ``a=0``."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from scripts.explore_h9_14_record14_a1_symmetry_z3 import (
    _swap_action,
    enumerate_orbit,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    BASE_G3,
    BASE_G4,
    EDGES,
    FULL,
    N,
    hits_all,
    independent,
    macaulay_raise,
)


def build_a0_problem() -> dict[str, object]:
    all_fives = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    lower = []
    upper = []
    for vertex in range(N):
        before = BASE_G4[vertex]
        lower.append(max(0, before - macaulay_raise(BASE_G3[vertex], 3)))
        upper.append(before - 1)
    candidates = tuple(
        face
        for face in all_fives
        if all(
            not (face & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(N)
        )
    )
    base = tuple(EDGES)
    transversals = tuple(
        value for value in range(1, FULL + 1) if hits_all(value, base)
    )
    return {
        "a": 0,
        "orbit_index": 0,
        "fours": (),
        "shadow": (),
        "addition_count": 8,
        "scheduled_facets": 110,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small": tuple(
            value for value in transversals if value.bit_count() <= 4
        ),
        "fives": tuple(
            value for value in transversals if value.bit_count() == 5
        ),
    }


def source_generators() -> tuple[tuple[int, ...], ...]:
    isolate_swaps = tuple(_swap_action(left, left + 1) for left in range(3, 11))
    return (*isolate_swaps, _swap_action(0, 1))


def enumerate_a0(max_models: int | None = None) -> dict[str, object]:
    result = enumerate_orbit(
        0,
        max_models,
        problem_override=build_a0_problem(),
        generator_override=source_generators(),
        result_kind="h9_14_record14_a0_full_source_generator_lex_slice_z3_census",
    )
    result.pop("source_positive_A4_orbit_index", None)
    result["source_branch_count"] = 1
    result["source_symmetry_statement"] = (
        "the adjacent isolate swaps generate S9 and the active twin swap "
        "generates S2, the full automorphism group of the Record-14 edge shadow"
    )
    return result


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-models", type=int)
    args = parser.parse_args(argv)
    print(json.dumps(enumerate_a0(args.max_models), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
