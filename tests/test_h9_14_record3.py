import json

from coxeter_polytopes.combinatorics.h9_14_record3 import (
    BASE_F_VECTOR,
    TARGET_F_VECTOR,
    audit_record3_polytopality,
    target_facet_complements,
)
from scripts.certify_h9_14_record3 import EXPECTED_CERTIFICATE_SHA256, OUTPUT, verify


def test_record3_exact_rational_stellar_multiwedge() -> None:
    audit = audit_record3_polytopality()
    assert audit.verified is True
    assert audit.stellar_facet_counts == (20, 28, 36)
    assert audit.radial_formulas_verified == (True, True)
    assert audit.coordinate_support_facet_count == 36
    assert audit.exact_coordinate_support_match is True
    assert audit.base_face_vector == BASE_F_VECTOR
    assert audit.expanded_minimal_nonfaces == audit.target_minimal_nonfaces
    assert audit.blocker_order_histogram == ((5, 162),)
    assert len(target_facet_complements()) == 162


def test_frozen_record3_ordinary_polytopality_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "facet_count": 162,
        "face_vector": list(TARGET_F_VECTOR),
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    incidence = manifest["global_incidence_audit"]
    assert incidence["passes_selection_model"] is True
    assert incidence["passes_connected_homology_audit"] is True
    assert manifest["polytopality_audit"]["ordinary_rational_polytopal"] is True
