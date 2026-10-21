import json

from coxeter_polytopes.combinatorics.h9_14_record18_polytopality import (
    BASE_F_VECTOR,
    audit_record18_polytopality,
)
from scripts.certify_h9_14_record18_polytopality import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)


def test_record18_exact_rational_base_and_inverse_multiwedge() -> None:
    audit = audit_record18_polytopality()
    assert audit.verified is True
    assert audit.base_face_vector == BASE_F_VECTOR
    assert audit.cyclic_facet_count == 40
    assert audit.first_stellar_facet_count == 52
    assert audit.second_stellar_facet_count == 62
    assert audit.first_radial_formula_verified is True
    assert audit.second_radial_formula_verified is True
    assert audit.coordinate_support_facet_count == 62
    assert audit.exact_coordinate_support_match is True
    assert audit.exact_target_base_match is True
    assert audit.exact_multiwedge_nonface_match is True


def test_frozen_record18_polytopality_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "base_facet_count": 62,
        "base_face_vector": list(BASE_F_VECTOR),
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
