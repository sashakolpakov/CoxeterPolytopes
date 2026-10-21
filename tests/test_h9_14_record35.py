import json

from coxeter_polytopes.combinatorics.h9_14_record35 import (
    TARGET_F_VECTOR,
    audit_record35_base_realization,
    audit_record35_incidence,
    audit_record35_multiwedge,
)
from scripts.certify_h9_14_record35 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)


def test_record35_complete_blocker_and_face_table() -> None:
    audit = audit_record35_incidence()
    assert audit.verified is True
    assert audit.minimal_transversal_size_histogram == ((5, 102),)
    assert len(audit.facet_complement_masks) == TARGET_F_VECTOR[-1]
    assert audit.facet_complement_sha256 == (
        "2e9b1664ed27d6f19310748efd4717cd96201ceb4ae94bf845faa9fb6cff1466"
    )
    assert audit.face_table_true_count == 5163
    assert audit.packed_face_table_sha256 == (
        "168a3653a6fd2db909153fc32169f3eb6ab3a3e99847b874d3b61d8f1fc7ba96"
    )


def test_record35_integral_base_and_inverse_multiwedge() -> None:
    base = audit_record35_base_realization()
    assert base.verified is True
    assert len(base.support_facets) == 19
    assert base.five_point_determinant_count == 126
    assert base.zero_five_point_determinant_count == 0
    assert base.minimum_absolute_determinant == 9
    assert base.minimum_facet_side_absolute_determinant == 21

    multiwedge = audit_record35_multiwedge()
    assert multiwedge.verified is True
    assert multiwedge.wedge_count == 5
    assert multiwedge.lifted_polytope_dimension == 9
    assert multiwedge.lifted_facet_count == 14


def test_frozen_record35_polytopality_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 102,
        "face_vector": list(TARGET_F_VECTOR),
        "base_facet_count": 19,
        "minimum_facet_side_absolute_determinant": 21,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
