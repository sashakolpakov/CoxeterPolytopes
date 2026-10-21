from __future__ import annotations

import json

from scripts.certify_h9_14_record1_a1_orbit0_target_orbit import (
    EXPECTED_CERTIFICATE_DIGEST,
    OUTPUT,
    fresh_record,
    verify,
)
from scripts.independent_h9_14_record1_a1_orbit0_target_orbit_audit import (
    replay as independent_replay,
)


def test_frozen_record1_target_orbit_theorem_and_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 1,
        "schedule_a": 1,
        "source_orbit_index": 0,
        "labelled_target_orbit_size": 60_480,
        "target_orbit_compact_label_no_go": True,
        "quotient_incidence_census_complete": False,
        "whole_source_orbit_compact_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }
    assert document["cutoff_correction"] == {
        "complete_quotient_incidence_census_available": False,
        "minimum_literal_pure_target_count": 60_480,
        "reason": (
            "The exploratory command used --max-models 10, preserved no UNSAT "
            "transcript, and the single frozen target already has 60,480 "
            "distinct source-stabilizer images."
        ),
        "retracted_statement": "exactly nine literal pure targets",
    }


def test_production_and_hostile_orbit_replays_agree_exactly() -> None:
    fresh = fresh_record()
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert fresh == frozen
    independent = independent_replay()
    orbit = frozen["orbit_evidence"]
    assert independent == orbit["independent_replay"]
    assert independent["source_stabilizer_structure"] == "S2 x S3 x S7"
    assert independent["source_stabilizer_order"] == 60_480
    assert independent["representative_target_stabilizer_order"] == 1
    assert independent["labelled_target_orbit_size"] == 60_480
    assert independent["representative_minimal_transversal_size_histogram"] == {
        "5": 158
    }
    assert independent["representative_ridge_degree_histogram"] == {
        "0": 2292,
        "2": 711,
    }
    assert independent["representative_compact_label_stage_counts"] == [
        108,
        84,
        8,
        0,
    ]
    assert independent["whole_labelled_target_orbit_compact_label_no_go"]
    assert not independent["quotient_incidence_census_complete"]
    assert not independent["whole_source_orbit_compact_no_go"]
