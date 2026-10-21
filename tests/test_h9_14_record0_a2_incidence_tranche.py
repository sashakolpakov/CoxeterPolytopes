import json
from copy import deepcopy

from scripts.certify_h9_14_record0_a2_incidence_tranche import (
    EXCLUDED_SOURCE_ORBITS,
    EXPECTED_SOURCE_REDUCTION,
    INCIDENCE_SAT_SOURCE_ORBITS,
    OUTPUT,
    SOLVER_OPEN_SOURCE_ORBITS,
    normalize_document,
    verify,
)
from scripts.explore_h9_14_record0_a2_dual_facets import (
    SOURCE_A4_FAMILIES,
    source_representatives,
    source_stabilizer_generators,
)
from scripts.independent_h9_14_record0_a2_tranche_z3 import SOURCES


def test_frozen_record0_a2_incidence_tranche() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(document)
    assert summary == {
        "mathematical_digest": document["mathematical_digest"],
        "exact_source_orbit_count": 36,
        "excluded_source_orbit_indices": list(EXCLUDED_SOURCE_ORBITS),
        "two_solver_closed_source_orbit_count": 8,
        "remaining_source_orbit_count": 28,
        "whole_record0_a2_schedule_excluded": False,
        "whole_record0_compact_shadow_excluded": False,
        "finite_volume_noncompact_checked": False,
    }
    assert normalize_document(document) == document


def test_source_partition_and_independent_literals() -> None:
    representatives = source_representatives()
    assert len(representatives) == 36
    assert sum(
        int(size) * count
        for size, count in EXPECTED_SOURCE_REDUCTION[
            "source_orbit_size_histogram"
        ].items()
    ) == EXPECTED_SOURCE_REDUCTION[
        "strict_aggregate_viable_labelled_count"
    ]
    for source in EXCLUDED_SOURCE_ORBITS:
        assert tuple(tuple(face) for face in SOURCE_A4_FAMILIES[source]) == SOURCES[source]
        actions, order = source_stabilizer_generators(representatives[source])
        assert actions
        assert order > 0
        for action in actions:
            assert all(action[action[vertex]] == vertex for vertex in range(14))


def test_strict_nonclaims() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    theorem = document["theorem"]
    coverage = document["coverage"]
    assert theorem["excluded_compact_source_orbits"] == list(
        EXCLUDED_SOURCE_ORBITS
    )
    assert theorem["whole_record0_a2_schedule_excluded"] is False
    assert theorem["whole_record0_compact_shadow_excluded"] is False
    assert theorem["finite_volume_noncompact_checked"] is False
    assert coverage["remaining_source_orbit_count"] == 28
    assert document["scout_status_ledger"] == {
        "two_solver_unsat_source_orbits": list(EXCLUDED_SOURCE_ORBITS),
        "explicit_incidence_sat_source_orbits": list(
            INCIDENCE_SAT_SOURCE_ORBITS
        ),
        "solver_open_source_orbits": list(SOLVER_OPEN_SOURCE_ORBITS),
        "partition_identity": "8 + 24 + 4 = 36",
        "sat_scope": (
            "An explicit solution of the compact-incidence overcover was "
            "found; this is not a Coxeter-label or Gram survivor."
        ),
        "open_scope": "Recorded timeouts are not exclusions.",
    }


def test_solver_telemetry_is_noncanonical() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    changed = deepcopy(document)
    changed["run_metadata"]["per_source"][0][
        "production_wall_time_seconds"
    ] += 123.0
    assert verify(changed)["mathematical_digest"] == document[
        "mathematical_digest"
    ]
