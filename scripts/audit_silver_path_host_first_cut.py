#!/usr/bin/env python3
"""Run the complete exact H22 silver-path first-cut transcript."""

from coxeter_polytopes.search.silver_path_host import (
    CONTROL_COUNT,
    INTEGRAL_CANDIDATE_COUNT,
    accepted_first_cut_prefix,
    enumerate_first_cut_prefix,
    first_cut_roots,
    silver_path_first_cut_audit,
)


def main() -> None:
    controls, candidates = enumerate_first_cut_prefix()
    accepted = accepted_first_cut_prefix()
    if len(controls) != CONTROL_COUNT:
        raise AssertionError(f"control count changed: {len(controls)}")
    if len(candidates) != INTEGRAL_CANDIDATE_COUNT:
        raise AssertionError(f"candidate count changed: {len(candidates)}")
    expected = ((1,) + (0,) * 22, *first_cut_roots())
    actual = tuple(candidate.vector for candidate in accepted)
    # Compare e0 without importing a second constructor spelling.
    if actual[1:] != expected[1:] or not all(
        actual[0][i] == expected[0][i] for i in range(23)
    ):
        raise AssertionError("accepted chronology changed")
    print(silver_path_first_cut_audit())
    print("full exact first-cut transcript: PASS")


if __name__ == "__main__":
    main()
