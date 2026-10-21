from __future__ import annotations

from functools import cache
import json

from scripts import certify_jacobian_kummer_h12_path_cube_family as cert


@cache
def payload() -> dict:
    return cert.build()


def test_all_eight_configurations_and_both_phases() -> None:
    result = payload()
    assert result["all_checks_pass"]
    assert result["complete_configuration_count"] == 8
    expected = [886 * n + 55 for n in range(1, 9)]
    for row in result["direct_N_1_through_8_census"]:
        assert row["phase_0_counts_N_1_through_8"] == expected
        assert row["phase_1_counts_N_1_through_8"] == expected


def test_canonical_exact_boundaries_and_cube_compatibility() -> None:
    result = payload()
    canonical = result["canonical_configuration"]
    assert canonical["N_1_matches_frozen_H12_941_cube"]
    assert canonical["N_2_matches_frozen_H12_1827_cube"]
    assert [row["facet_count"] for row in canonical["boundaries"]] == [
        941,
        1827,
        2713,
        3599,
        4485,
        5371,
        6257,
        7143,
    ]
    assert all(
        row["full_linear_rank"] == 13
        and row["all_boundary_roots_primitive"]
        and row["all_boundary_reflections_integral"]
        for row in canonical["boundaries"]
    )


def test_frozen_archive_matches_fresh_build() -> None:
    archived = json.loads(cert.OUTPUT.read_text(encoding="utf-8"))
    assert archived == payload()
