import json

from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    audit_record43_incidence,
)
from scripts.certify_h9_14_record43 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record43_audit import replay


def test_record43_complete_blocker_is_pure_with_ninety_facets() -> None:
    audit = audit_record43_incidence(TARGET_GRAPH6)
    assert audit.verified is True
    assert audit.minimal_transversal_size_histogram == ((5, 90),)
    assert len(audit.facet_complement_masks) == TARGET_F_VECTOR[-1] == 90
    assert audit.face_table_true_count == 4743
    assert (
        audit.packed_face_table_sha256
        == "32cb562b346841aa4c4a850a911a04da96d907cf6503ed6947f4960bd0212b51"
    )


def test_frozen_record43_incidence_and_gale_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 90,
        "face_vector": list(TARGET_F_VECTOR),
        "gale_sign_count": 1001,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }


def test_record43_hostile_replay() -> None:
    summary = replay()
    assert summary["minimal_transversal_size_histogram"] == {"5": 90}
    assert summary["ridge_degree_histogram"] == {"2": 405}
    assert summary["dual_component_count"] == 1
    assert summary["all_links_mod2_homology_spheres"] is True
    assert summary["positive_circuit_count"] == 90
    assert summary["gale_audit_passes"] is True
