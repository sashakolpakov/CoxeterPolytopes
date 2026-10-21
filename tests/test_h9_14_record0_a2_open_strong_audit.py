import json
from copy import deepcopy

from scripts.certify_h9_14_record0_a2_open_strong_audit import (
    EXPECTED_SOURCE29_STRICT_FAILURE_DIGEST,
    OPEN_SOURCES,
    OUTPUT,
    SOURCE29_SHADOW_WITNESS,
    verify,
)


def test_frozen_open_strong_audit() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "mathematical_digest": document["mathematical_digest"],
        "strong_audit_open_source_orbits": list(OPEN_SOURCES),
        "frontier_changed": False,
        "source29_scheduled_shadow_witness": True,
        "source29_strict_link_failure_count": 48,
        "new_hyperbolic_coxeter_polytope": False,
    }


def test_frontier_and_nonclaims_remain_strict() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert document["frontier"]["partition_identity"] == "8 + 24 + 4 = 36"
    assert document["frontier"]["strong_audit_open_source_orbits"] == list(
        OPEN_SOURCES
    )
    assert document["frontier"]["frontier_changed"] is False
    assert document["theorem"] == {
        "new_compact_incidence_exclusion": False,
        "new_compact_incidence_sat_witness": False,
        "new_hyperbolic_coxeter_polytope": False,
        "finite_volume_noncompact_checked": False,
    }


def test_source29_shadow_witness_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    witness = document["positive_shadow_witnesses"][0]
    assert witness["source_orbit_index"] == 29
    assert witness["additional_A5_masks"] == list(SOURCE29_SHADOW_WITNESS)
    assert witness["scheduled_missing_face_counts"] == {
        "6": 147,
        "7": 420,
        "8": 629,
    }
    assert witness["strict_link_failure_count"] == 48
    assert witness["strict_link_failure_histogram_by_face_size"] == {
        "1": 0,
        "2": 2,
        "3": 46,
    }
    assert (
        witness["strict_link_failure_digest"]
        == EXPECTED_SOURCE29_STRICT_FAILURE_DIGEST
    )
    assert witness["minimal_transversal_order_histogram"] == {
        "4": 10,
        "5": 90,
        "6": 27,
    }
    assert witness["full_compact_incidence_model_pass"] is False


def test_timeout_telemetry_is_noncanonical() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    changed = deepcopy(document)
    changed["run_metadata"]["full_primal_no_lex"][0][
        "wall_time_seconds"
    ] += 1000
    assert verify(changed)["mathematical_digest"] == document[
        "mathematical_digest"
    ]
