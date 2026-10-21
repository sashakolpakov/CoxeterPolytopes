from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record12_all_coxeter_labels import (
    audit_record12_all_coxeter_labels,
    load_record12_residual_nerves,
)
from scripts.certify_h9_14_record12_all_coxeter_labels import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_HOSTILE_SUMMARY_SHA256,
    INPUT,
    OUTPUT,
    generate,
    verify,
)
from scripts.independent_h9_14_record12_all_coxeter_audit import replay
from scripts.independent_h9_14_record12_coxeter_probe import (
    EXPECTED_ALL_SUMMARY_SHA256 as EXPECTED_SECONDARY_SUMMARY_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record12_all_coxeter_audit.py"
)
SECONDARY_SOURCE = (
    ROOT / "scripts/independent_h9_14_record12_coxeter_probe.py"
)
EXPECTED_SUPPORT_HISTOGRAM = ((2, 49), (3, 33), (4, 152), (5, 129), (6, 8))
EXPECTED_TRANSCRIPT_SHA256 = (
    "29e7904bd2af75c692f24b0a4c502426c5316a8b2d1926ca4b07d3d7711f810a"
)
EXPECTED_ARCHIVE_SHA256 = (
    "fe6c0ee7f494a67ecac88b51aed301a4fdfe26fbddbf68034916f7b7aecc97eb"
)


@pytest.fixture(scope="session")
def all_orbit_audit():
    return audit_record12_all_coxeter_labels(INPUT)


@pytest.fixture(scope="session")
def hostile_replay():
    return replay(INPUT)


def test_all_371_residual_orbits_reconstruct_and_reach_zero(
    all_orbit_audit,
) -> None:
    nerves = load_record12_residual_nerves(INPUT)
    assert len(nerves) == 371
    assert Counter(nerve.missing_four_count for nerve in nerves) == {
        0: 48,
        1: 95,
        2: 144,
        3: 63,
        4: 21,
    }
    assert Counter(len(nerve.rank_nine_vertices) for nerve in nerves) == {
        122: 21,
        124: 63,
        126: 144,
        128: 95,
        130: 48,
    }
    assert all_orbit_audit.verified
    assert all_orbit_audit.whole_shadow_compact_coxeter_no_go
    assert all_orbit_audit.support_count_histogram == EXPECTED_SUPPORT_HISTOGRAM
    assert all_orbit_audit.orbit_transcript_sha256 == EXPECTED_TRANSCRIPT_SHA256
    assert sum(len(orbit.stages) for orbit in all_orbit_audit.orbits) == 1498
    assert max(
        stage.passing_state_count
        for orbit in all_orbit_audit.orbits
        for stage in orbit.stages
    ) == 864
    assert [orbit.stages[-1].passing_state_count for orbit in all_orbit_audit.orbits] == [
        0
    ] * 371


def test_hostile_no_import_replay_matches_every_stage(
    all_orbit_audit, hostile_replay
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    secondary_source = SECONDARY_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert "from coxeter_polytopes" not in secondary_source
    assert "import coxeter_polytopes" not in secondary_source
    assert hostile_replay["all_371_orbits_excluded"] is True
    assert hostile_replay["summary_sha256"] == EXPECTED_HOSTILE_SUMMARY_SHA256
    assert hostile_replay["support_count_histogram"] == {
        str(key): value for key, value in EXPECTED_SUPPORT_HISTOGRAM
    }
    assert len(hostile_replay["orbit_gates"]) == 371
    for production, independent in zip(
        all_orbit_audit.orbits, hostile_replay["orbit_gates"], strict=True
    ):
        assert production.nerve.orbit_index == independent["orbit_index"]
        assert list(production.decisive_support_indices) == independent[
            "decisive_support_indices"
        ]
        assert [stage.to_record() for stage in production.stages] == independent[
            "stages"
        ]


def test_frozen_certificate_rebuilds_with_whole_shadow_scope() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_SHA256
    assert archived == generate()
    summary = verify(archived)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["hostile_summary_sha256"] == EXPECTED_HOSTILE_SUMMARY_SHA256
    assert summary["secondary_summary_sha256"] == EXPECTED_SECONDARY_SUMMARY_SHA256
    assert summary["residual_orbits"] == 371
    assert summary["whole_graph_shadow_compact_coxeter_excluded"] is True
    assert archived["secondary_independent_replay"] == {
        "method": (
            "standard-library only; corrected matrix reconstruction, two-way "
            "actual rank-nine face equality, independent 9+5 Lanner closure, "
            "independent finite Coxeter classifier, no label-state symmetry "
            "quotient, and deterministic touching-first support joins"
        ),
        "input_file_sha256": (
            "faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05"
        ),
        "input_certificate_sha256": (
            "50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4"
        ),
        "lanner_tables": {
            "rank_four_placements": 108,
            "rank_four_isomorphism_types": 9,
            "rank_five_placements": 420,
            "rank_five_isomorphism_types": 5,
        },
        "orbit_count": 371,
        "schedule_histogram": {"0": 48, "1": 95, "2": 144, "3": 63, "4": 21},
        "support_count_histogram": {"2": 187, "3": 164, "4": 15, "5": 5},
        "stage_count": 951,
        "maximum_transient_state_count": 3481,
        "all_371_orbits_excluded": True,
        "full_replay_summary_sha256": EXPECTED_SECONDARY_SUMMARY_SHA256,
    }
    assert archived["scope"] == {
        "record_index": 12,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "selected_nerve_only": False,
        "all_371_residual_incidence_orbits_covered": True,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "pl_sphericity_of_abstract_residuals_claimed": False,
        "convex_polytopality_of_abstract_residuals_claimed": False,
        "global_lorentzian_gram_values_enumerated": False,
    }


def test_certificate_records_global_label_bound_and_lossless_pruning() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    bound = archived["proof"]["global_finite_label_bound"]
    pruning = archived["proof"]["perron_frobenius_step"]
    assert "Lemma 10.10" in bound
    assert "m in {2,3,4,5}" in bound
    assert archived["theorem"]["finite_label_alphabet"] == [2, 3, 4, 5]
    assert "G=I-A" in pruning
    assert "rho(A)<1" in pruning
    assert "rho(A')>=rho(A)" in pruning
    assert archived["proof"]["why_no_global_gram_search_is_needed"]
