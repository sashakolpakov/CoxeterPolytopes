from __future__ import annotations

import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record17_all_coxeter_labels import (
    EXPECTED_ORBIT_SIGNATURES,
    EXPECTED_STAGE_COUNTS,
    EXPECTED_STAGE_DIGESTS,
    audit_record17_all_coxeter_labels,
)
from scripts.certify_h9_14_record17_all_coxeter_labels import (
    EXPECTED_CERTIFICATE_SHA256,
    INPUT,
    OUTPUT,
    generate,
    verify,
)
from scripts.independent_h9_14_record17_all_coxeter_audit import (
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record17_all_coxeter_audit.py"
)
HOSTILE_DIGEST = (
    "e19c45bb7134690b3649b9476c69829072076a777335c124216aeee5b0e3b609"
)


@pytest.fixture(scope="session")
def all_orbit_audit():
    return audit_record17_all_coxeter_labels(INPUT)


@pytest.fixture(scope="session")
def hostile_replay():
    return replay(INPUT)


def test_every_residual_orbit_has_an_exact_zero_survivor_gate(
    all_orbit_audit,
) -> None:
    assert all_orbit_audit.verified
    assert all_orbit_audit.whole_shadow_compact_coxeter_no_go
    assert tuple(orbit.nerve.signature for orbit in all_orbit_audit.orbits) == (
        EXPECTED_ORBIT_SIGNATURES
    )
    assert tuple(
        tuple(stage.count_signature for stage in orbit.stages)
        for orbit in all_orbit_audit.orbits
    ) == EXPECTED_STAGE_COUNTS
    assert tuple(
        tuple(stage.passing_state_sha256 for stage in orbit.stages)
        for orbit in all_orbit_audit.orbits
    ) == EXPECTED_STAGE_DIGESTS
    assert [orbit.stages[-1].passing_state_count for orbit in all_orbit_audit.orbits] == [
        0
    ] * 8


def test_hostile_replay_is_no_import_and_matches_every_state_hash(
    all_orbit_audit, hostile_replay
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert hostile_replay["all_eight_orbits_excluded"] is True
    assert hostile_replay["summary_sha256"] == HOSTILE_DIGEST
    for production, independent in zip(
        all_orbit_audit.orbits, hostile_replay["orbit_gates"], strict=True
    ):
        assert [stage.to_record() for stage in production.stages] == independent[
            "stages"
        ]


def test_frozen_certificate_rebuilds_with_whole_shadow_scope() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived == generate()
    summary = verify(archived)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["whole_graph_shadow_compact_coxeter_excluded"] is True
    assert summary["final_survivors"] == [0] * 8
    assert archived["scope"] == {
        "record_index": 17,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "selected_nerve_only": False,
        "all_eight_residual_incidence_orbits_covered": True,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "pl_sphericity_of_abstract_residuals_claimed": False,
        "convex_polytopality_of_abstract_residuals_claimed": False,
        "global_lorentzian_gram_values_enumerated": False,
    }


def test_archive_states_the_exact_perron_frobenius_implication() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    proof = archived["proof"]["perron_frobenius_step"]
    assert "G=I-A" in proof
    assert "rho(A)<1" in proof
    assert "rho(A')>=rho(A)" in proof
    assert archived["proof"]["why_no_global_gram_search_is_needed"]
