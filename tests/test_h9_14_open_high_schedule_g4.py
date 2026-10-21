from __future__ import annotations

import json

from scripts.certify_h9_14_open_high_schedule_g4 import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_EXCLUSIONS,
    OUTPUT,
    generate,
    verify,
)
from scripts.independent_h9_14_open_high_schedule_g4 import replay


def test_independent_strict_link_g4_cutoffs() -> None:
    audit = replay()
    assert audit["total_schedule_branches_excluded"] == 14
    assert {
        row["record_index"]: tuple(row["strict_link_g4_excluded_schedules"])
        for row in audit["records"]
    } == EXPECTED_EXCLUSIONS
    for record in audit["records"]:
        sums = [row["summed_vertex_link_g4"] for row in record["schedules"]]
        assert all(first - second == 6 for first, second in zip(sums, sums[1:]))
        assert sums[-2:] == [12, 6]


def test_frozen_strict_link_g4_certificate_rebuilds() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived == generate()
    summary = verify(archived)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["total_schedule_branches_excluded"] == 14
    assert archived["theorem"]["verified"] is True


def test_strict_link_g4_scope_is_not_overpromoted() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived["scope"] == {
        "compact_only": True,
        "record_indices": [0, 1, 3, 4, 5, 7, 9],
        "listed_schedule_branches_excluded": True,
        "whole_graph_shadows_excluded": False,
        "remaining_schedules_excluded": False,
        "finite_volume_noncompact_chambers_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    assert "generalized lower-bound equality theorem" in archived["proof"][
        "strict_local_bound"
    ]
