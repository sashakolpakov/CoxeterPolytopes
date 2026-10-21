import json

from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    audit_record18_incidence,
    audit_record18_structure,
)
from scripts.certify_h9_14_record18 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record18_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


def test_record18_complete_blocker_and_structural_reduction() -> None:
    incidence = audit_record18_incidence(TARGET_GRAPH6)
    structure = audit_record18_structure()
    assert incidence.verified is True
    assert incidence.minimal_transversal_size_histogram == ((5, 110),)
    assert incidence.face_table_true_count == 5379
    assert (
        incidence.facet_complement_sha256
        == "969da679ad87ca730e822ceb7c5a4eee7f14b461062cf71543293c50c9872afc"
    )
    assert (
        incidence.packed_face_table_sha256
        == "620c83d073c4a141ef2ec5cf9333326cdc8d7c7e94d5e48e647b19b030943eaf"
    )
    assert structure.verified is True
    assert (structure.selected_y_count, structure.selected_z_count, structure.selected_w_count) == (60, 25, 25)


def test_frozen_record18_candidate_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "record_index": 18,
        "status": "abstract_audited_nerve_candidate",
        "facet_count": 110,
        "face_vector": list(TARGET_F_VECTOR),
        "gale_sign_count": 1001,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }


def test_record18_hostile_replay() -> None:
    summary = replay()
    assert summary["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert summary["minimal_transversal_size_histogram"] == {"5": 110}
    assert summary["ridge_degree_histogram"] == {"2": 495}
    assert summary["dual_component_count"] == 1
    assert summary["all_links_mod2_homology_spheres"] is True
    assert summary["positive_circuit_count"] == 110
    assert summary["gale_audit_passes"] is True
    assert summary["status"] == "abstract_audited_nerve_candidate"
