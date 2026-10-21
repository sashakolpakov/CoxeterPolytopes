import json

from coxeter_polytopes.combinatorics.h9_14_record5 import (
    BASE_F_VECTOR,
    TARGET_F_VECTOR,
    audit_record5_polytopality,
    target_facet_complements,
)
from scripts.certify_h9_14_record5 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)


def test_record5_exact_rational_stellar_multiwedge() -> None:
    audit = audit_record5_polytopality()
    assert audit.verified is True
    assert audit.base_face_vector == BASE_F_VECTOR
    assert (
        audit.cyclic_facet_count,
        audit.first_stellar_facet_count,
        audit.second_stellar_facet_count,
    ) == (20, 28, 34)
    assert audit.first_radial_formula_verified is True
    assert audit.second_radial_formula_verified is True
    assert audit.coordinate_support_facet_count == 34
    assert audit.exact_coordinate_support_match is True
    assert audit.exact_multiwedge_nonface_match is True
    assert audit.blocker_order_histogram == ((5, 136),)
    assert len(target_facet_complements()) == TARGET_F_VECTOR[-1]


def test_frozen_record5_ordinary_polytopality_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "facet_count": 136,
        "face_vector": list(TARGET_F_VECTOR),
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    incidence = manifest["global_incidence_audit"]
    assert incidence["matching_schedule_index"] == 3
    assert incidence["missing_four_count"] == 3
    assert incidence["passes_selection_model"] is True
    assert incidence["passes_connected_homology_audit"] is True
    assert manifest["polytopality_audit"]["ordinary_rational_polytopal"] is True
