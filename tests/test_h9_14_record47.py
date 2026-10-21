from collections import Counter
import json

from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_FIVES,
    audit_record47_incidence,
    record47_minimal_transversals,
)
from scripts.certify_h9_14_record47 import OUTPUT, verify
from scripts.independent_h9_14_record47_audit import replay


def test_record47_complete_blocker_is_pure_with_ninety_facets() -> None:
    blocker = record47_minimal_transversals()
    assert Counter(mask.bit_count() for mask in blocker) == {5: 90}
    audit = audit_record47_incidence(TARGET_GRAPH6)
    assert audit.verified is True
    assert audit.missing_fives == TARGET_MISSING_FIVES
    assert audit.facet_complement_sha256 == (
        "8306d9c6c39a90aae519bd5ae587743b3737371f9efbe01ebd3757105b0b56ed"
    )
    assert audit.packed_face_table_sha256 == (
        "4d07f78cba8551c9d0fe9fc4d6041f3df8ae35c5017e9942e10005686c628aba"
    )


def test_frozen_record47_incidence_and_gale_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 90,
        "face_vector": list(TARGET_F_VECTOR),
        "gale_sign_count": 1001,
        "certificate_sha256": (
            "960ef6564c6b8013caf0989c95a46cc8401c846d32d47c66c057091340301546"
        ),
    }
    assert manifest["gale_sign_sha256"] == (
        "b8a9b3c5a8b4712e908892a53e75934e49db2e6a2fa6856c1428a38c8d0749b2"
    )
    assert manifest["gale_audit"]["positive_circuit_count"] == 90
    assert manifest["gale_audit"]["independent_sign_audit_passes"] is True


def test_record47_hostile_replay() -> None:
    result = replay()
    assert result["minimal_transversal_size_histogram"] == {"5": 90}
    assert result["face_vector"] == list(TARGET_F_VECTOR)
    assert result["ridge_degree_histogram"] == {"2": 405}
    assert result["dual_component_count"] == 1
    assert result["all_links_mod2_homology_spheres"] is True
    assert result["gale_sign_count"] == 1001
    assert result["positive_circuit_count"] == 90
    assert result["gale_audit_passes"] is True
