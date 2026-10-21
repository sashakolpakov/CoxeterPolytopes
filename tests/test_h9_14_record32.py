from collections import Counter
import json

from coxeter_polytopes.combinatorics.h9_14_record32 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    audit_record32_incidence,
    audit_record32_multiwedge,
    record32_minimal_transversals,
)
from scripts.certify_h9_14_record32 import OUTPUT, verify
from scripts.independent_h9_14_record32_audit import replay


def test_record32_complete_blocker_is_pure_with_110_facets() -> None:
    blocker = record32_minimal_transversals()
    assert Counter(mask.bit_count() for mask in blocker) == {5: 110}
    audit = audit_record32_incidence(TARGET_GRAPH6)
    assert audit.verified is True
    assert audit.missing_fours == TARGET_MISSING_FOURS
    assert audit.missing_fives == TARGET_MISSING_FIVES
    assert audit.facet_complement_sha256 == (
        "271b5bc535b7ec43295a29d375a820b8288c91a59f31fc005468eb3e2263954e"
    )
    assert audit.packed_face_table_sha256 == (
        "ac1ef3f84874cbbdeab6e5f8a2082f03b5b0309636eaaea9459c1575d5504f66"
    )


def test_record32_has_exact_convex_multiwedge_factorization() -> None:
    audit = audit_record32_multiwedge()
    assert audit.verified is True
    assert audit.base_face_vector == (7, 15, 10)
    assert audit.target_polytope_dimension == 9
    assert audit.target_vertex_count == 14


def test_frozen_record32_polytopality_and_gale_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "facet_count": 110,
        "face_vector": list(TARGET_F_VECTOR),
        "ordinary_polytopal": True,
        "gale_sign_count": 1001,
        "certificate_sha256": (
            "94a65380b2d56d4787a3de28d5cca5dede0540cc2ec0ff23b2a3d65b4cfbd318"
        ),
    }
    assert manifest["gale_sign_sha256"] == (
        "f9972dda2d51d44c0fc7d86d7220c8c5a4f4f8111260c94b602e06d4eff69222"
    )
    assert manifest["gale_audit"]["positive_circuit_count"] == 110
    assert manifest["gale_audit"]["independent_sign_audit_passes"] is True


def test_record32_hostile_replay() -> None:
    result = replay()
    assert result["minimal_transversal_size_histogram"] == {"5": 110}
    assert result["face_vector"] == list(TARGET_F_VECTOR)
    assert result["ridge_degree_histogram"] == {"2": 495}
    assert result["dual_component_count"] == 1
    assert result["all_links_mod2_homology_spheres"] is True
    assert result["base_face_vector"] == [7, 15, 10]
    assert result["ordinary_polytopal"] is True
    assert result["gale_sign_count"] == 1001
    assert result["positive_circuit_count"] == 110
    assert result["gale_audit_passes"] is True
