from __future__ import annotations

import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_open_selected_coxeter_no_go import (
    TARGET_RECORDS,
    audit_all_open_selected_nerves,
)
from scripts.certify_h9_14_open_selected_coxeter_no_go import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    generate,
    verify,
)
from scripts.independent_h9_14_open_selected_coxeter_audit import replay


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_DECISIVE_SUPPORTS = (
    (0, 7, 9, 2, 1),
    (0, 4, 1, 5, 6, 8),
    (0, 7, 4, 5, 8),
    (0, 3, 10, 4),
    (0, 4, 5),
    (0, 7, 11, 1, 8),
    (0, 1, 5),
)


@pytest.fixture(scope="session")
def production_audit():
    return audit_all_open_selected_nerves()


@pytest.fixture(scope="session")
def hostile_replay():
    return replay()


def test_all_seven_selected_nerves_reach_exact_zero(production_audit) -> None:
    assert production_audit.verified
    assert production_audit.all_seven_selected_nerves_excluded
    assert tuple(audit.nerve.record_index for audit in production_audit.audits) == (
        TARGET_RECORDS
    )
    assert tuple(
        tuple(stage.support_index for stage in audit.stages)
        for audit in production_audit.audits
    ) == EXPECTED_DECISIVE_SUPPORTS
    assert [audit.stages[-1].passing_state_count for audit in production_audit.audits] == [
        0
    ] * 7


def test_hostile_replay_matches_every_stage_and_state_hash(
    production_audit, hostile_replay
) -> None:
    source = (
        ROOT / "scripts/independent_h9_14_open_selected_coxeter_audit.py"
    ).read_text(encoding="utf-8")
    assert "h9_14_open_selected_coxeter_no_go" not in source
    assert hostile_replay["all_seven_selected_nerves_excluded"] is True
    assert hostile_replay["lanner_table_counts"] == {"4": 108, "5": 420}
    for production, independent in zip(
        production_audit.audits, hostile_replay["record_gates"], strict=True
    ):
        assert [stage.to_record() for stage in production.stages] == independent[
            "stages"
        ]


def test_frozen_selected_nerve_certificate_rebuilds() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived == generate()
    summary = verify(archived)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["record_indices"] == list(TARGET_RECORDS)
    assert summary["final_survivors"] == [0] * 7
    assert archived["scope"] == {
        "record_indices": list(TARGET_RECORDS),
        "selected_nerves_only": True,
        "whole_graph_shadows_compact_coxeter_excluded": False,
        "alternative_completions_excluded": False,
        "compact_only": True,
        "finite_volume_noncompact_chambers_excluded": False,
        "global_lorentzian_gram_values_enumerated": False,
    }


def test_archive_states_lossless_perron_frobenius_scope() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    proof = archived["proof"]["perron_frobenius_step"]
    assert "G=I-A" in proof
    assert "rho(A)<1" in proof
    assert "rho(A')>=rho(A)" in proof
    assert archived["scope"]["selected_nerves_only"] is True
    assert archived["scope"]["whole_graph_shadows_compact_coxeter_excluded"] is False
