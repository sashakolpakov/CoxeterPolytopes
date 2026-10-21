from collections import Counter

from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_FIVES,
    audit_record44_incidence,
    record44_minimal_transversals,
)
from scripts.certify_h9_14_record44 import OUTPUT, verify
from scripts.independent_h9_14_record44_audit import replay

import json


def test_record44_complete_blocker_is_pure_with_ninety_facets() -> None:
    blocker = record44_minimal_transversals()
    assert Counter(mask.bit_count() for mask in blocker) == {5: 90}
    audit = audit_record44_incidence(TARGET_GRAPH6)
    assert audit.verified is True
    assert audit.missing_fives == TARGET_MISSING_FIVES
    assert audit.facet_complement_sha256 == (
        "a706233135e61171db1da19172f4cb5234191f0a0a703a3388b4a4d9a9877917"
    )
    assert audit.packed_face_table_sha256 == (
        "449e514bd9f090a3c52dba38db833014e696c7d9cb20753a6b8f4f7390cb6ad6"
    )


def test_frozen_record44_incidence_and_gale_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 90,
        "face_vector": list(TARGET_F_VECTOR),
        "gale_sign_count": 1001,
        "certificate_sha256": (
            "e66e493d0c68f7ae28ffe03e63f96209276a8d83a5ceb27c498f35b42c4a1867"
        ),
    }
    assert manifest["gale_sign_sha256"] == (
        "2e7c05277d8f34686093fc27b5c3a170be321edefad31c8eff586943aaa39444"
    )
    assert manifest["gale_audit"]["positive_circuit_count"] == 90
    assert manifest["gale_audit"]["independent_sign_audit_passes"] is True


def test_record44_hostile_replay() -> None:
    result = replay()
    assert result["minimal_transversal_size_histogram"] == {"5": 90}
    assert result["face_vector"] == list(TARGET_F_VECTOR)
    assert result["ridge_degree_histogram"] == {"2": 405}
    assert result["dual_component_count"] == 1
    assert result["all_links_mod2_homology_spheres"] is True
    assert result["gale_sign_count"] == 1001
    assert result["positive_circuit_count"] == 90
    assert result["gale_audit_passes"] is True
