import json

from coxeter_polytopes.combinatorics.h9_14_record14 import (
    BASE_GALE_MATRIX,
    BASE_POSITIVE_DEPENDENCE,
    TARGET_F_VECTOR,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    audit_record14_polytopality,
    record14_facet_complements,
)
from scripts.certify_h9_14_record14 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record14_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


def test_record14_clutter_has_a_pure_106_member_blocker() -> None:
    assert tuple(map(len, TARGET_MISSING_EDGES)) == (2,) * 5
    assert tuple(map(len, TARGET_MISSING_FOURS)) == (4,) * 2
    assert tuple(map(len, TARGET_MISSING_FIVES)) == (5,) * 4
    facets = record14_facet_complements()
    assert len(facets) == 106
    assert all(value.bit_count() == 5 for value in facets)


def test_record14_exact_gale_multiwedge_is_ordinary_polytopal() -> None:
    audit = audit_record14_polytopality()
    assert audit.verified is True
    assert audit.face_vector == TARGET_F_VECTOR
    assert audit.ridge_degree_histogram == ((2, 477),)
    assert audit.base_f_vector == (9, 28, 38, 19)
    assert audit.base_nonzero_minor_count == 126
    assert audit.base_positive_circuit_count == 19
    assert audit.base_positive_circuits_match_facets is True
    assert audit.base_positive_dependence_verified is True
    assert audit.clone_block_sizes == (3, 2, 3)
    assert audit.target_dimension == 9
    assert len(BASE_GALE_MATRIX) == 4
    assert len(BASE_POSITIVE_DEPENDENCE) == 9


def test_record14_frozen_certificate_and_hostile_replay() -> None:
    hostile = independent_audit()
    assert hostile["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile["ordinary_polytopal"] is True
    assert hostile["positive_kernel_product"] == [0, 0, 0, 0]
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_index": 14,
        "status": "certified_ordinary_polytopal_completion",
        "facet_count": 106,
        "face_vector": list(TARGET_F_VECTOR),
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert manifest["vinberg_scope"]["finite_volume_criterion_checked"] is False
    assert manifest["construction"]["exact_minimal_nonface_match"] is True
