from collections import Counter
import json

from coxeter_polytopes.combinatorics.h9_14_record32_a0 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    audit_record32_a0_common_four_core,
    audit_record32_a0_incidence,
    record32_a0_minimal_transversals,
)
from scripts.certify_h9_14_record32_a0 import OUTPUT, verify
from scripts.independent_h9_14_record32_a0_audit import replay


def test_record32_a0_blocker_is_pure_with_114_facets() -> None:
    blocker = record32_a0_minimal_transversals()
    assert Counter(mask.bit_count() for mask in blocker) == {5: 114}
    audit = audit_record32_a0_incidence(TARGET_GRAPH6)
    assert audit.verified is True
    assert audit.facet_complement_sha256 == (
        "f1d397014fbda68430f8dfc0bc557138984dc79088bc055df891f03a7a20cdf2"
    )
    assert audit.packed_face_table_sha256 == (
        "aa97de9f03cb3f07c56588c5abf33e2597e73fba155b9735b35e4c423771d513"
    )


def test_record32_a0_common_four_core_is_unrestricted_no_go() -> None:
    audit = audit_record32_a0_common_four_core()
    assert audit.join_counts == (420, 1332, 4620)
    assert tuple(
        (stage.order, len(stage.restrictions), stage.passing_state_count)
        for stage in audit.stages
    ) == ((4, 13, 4620), (5, 12, 96), (6, 4, 0))
    assert audit.sparse_shadow_passing_state_count == 0
    assert audit.unrestricted_external_label_no_go is True


def test_frozen_record32_a0_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 114,
        "face_vector": list(TARGET_F_VECTOR),
        "rank_six_common_core_survivors": 0,
        "unrestricted_compact_coxeter_no_go": True,
        "certificate_sha256": (
            "14c3f9f5a425eb6934d7f53a7a7f071467eba7cf2e4a64b7cf9a194d3490d0bb"
        ),
    }


def test_record32_a0_hostile_replay() -> None:
    result = replay()
    assert result["minimal_transversal_size_histogram"] == {"5": 114}
    assert result["face_vector"] == list(TARGET_F_VECTOR)
    assert result["ridge_degree_histogram"] == {"2": 513}
    assert result["dual_component_count"] == 1
    assert result["all_links_mod2_homology_spheres"] is True
    assert result["common_four_core"]["join_counts"] == [420, 1332, 4620]
    assert result["common_four_core"]["staged_survivor_counts"][-1] == [6, 4, 0]
    assert result["common_four_core"]["unrestricted_external_label_no_go"] is True
    assert result["ordinary_polytopality_claimed"] is False
