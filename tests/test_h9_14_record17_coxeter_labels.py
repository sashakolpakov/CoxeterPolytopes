import json

from coxeter_polytopes.combinatorics.h9_14_record17_coxeter_labels import (
    COMMON_CORE,
    OUTER_FACETS,
    audit_record17_coxeter_labels,
    record17_vertex_facets,
)
from scripts.certify_h9_14_record17_coxeter_labels import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record17_coxeter_labels_audit import (
    compare_archive,
    replay,
)


def test_record17_common_core_lanner_gate_is_locally_impossible() -> None:
    audit = audit_record17_coxeter_labels()
    assert audit.verified
    assert audit.vertex_count == 110
    assert audit.gate.outer_facets == OUTER_FACETS == (0, 1, 2)
    assert audit.gate.common_core == COMMON_CORE == (3, 4, 5, 6)
    assert audit.gate.join_counts == (420, 1_332, 4_620)
    assert [stage.passing_state_count for stage in audit.gate.stages] == [
        4_620,
        80,
        0,
    ]
    assert audit.unrestricted_compact_coxeter_no_go


def test_record17_decisive_cluster_has_no_dotted_pair() -> None:
    audit = audit_record17_coxeter_labels()
    vertices = set(record17_vertex_facets())
    assert len(vertices) == 110
    assert audit.dotted_pairs_in_cluster_union == ()
    assert all(
        any(set(restriction) <= set(vertex) for vertex in vertices)
        for restriction in audit.gate.actual_vertex_restrictions
    )


def test_frozen_record17_compact_coxeter_no_go_replays() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 17,
        "rank_nine_vertices": 110,
        "joined_lanner_states": 4_620,
        "rank_five_survivors": 80,
        "rank_six_survivors": 0,
        "compact_coxeter_realization_exists": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["scope"]["selected_nerve_only"] is True
    assert document["scope"]["whole_graph_shadow_excluded"] is False
    assert document["scope"]["finite_volume_noncompact_chamber_excluded"] is False


def test_independent_record17_label_replay_agrees_with_archive() -> None:
    result = replay()
    compare_archive(result)
    assert result["rank_six_failure_multiplicity_histogram"] == {
        "1": 72,
        "2": 8,
    }
    assert result["unrestricted_compact_coxeter_no_go"] is True
