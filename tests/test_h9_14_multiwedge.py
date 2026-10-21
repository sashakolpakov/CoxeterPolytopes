import json

from coxeter_polytopes.combinatorics.h9_14_multiwedge import audit_all
from scripts.certify_h9_14_multiwedge import OUTPUT, verify
from scripts.independent_h9_14_multiwedge_audit import replay


def test_records_43_44_47_have_exact_multiwedge_factorizations() -> None:
    audits = audit_all()
    assert [audit.record_index for audit in audits] == [43, 44, 47]
    assert all(audit.verified for audit in audits)
    assert [audit.base_face_vector for audit in audits] == [
        (7, 15, 10),
        (6, 6),
        (7, 15, 10),
    ]
    assert [audit.wedge_count for audit in audits] == [2, 3, 2]
    assert [audit.clone_vertices for audit in audits] == [
        (8, 9, 10),
        (7, 8, 9, 10),
        (8, 9, 10),
    ]
    assert all(audit.exact_nonface_match for audit in audits)
    assert all(audit.multiwedge_factor_dimension == 5 for audit in audits)
    assert all(audit.product_dimension == 9 for audit in audits)
    assert all(audit.product_facet_count == 14 for audit in audits)


def test_frozen_multiwedge_polytopality_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_indices": [43, 44, 47],
        "ordinary_polytopal_count": 3,
        "certificate_sha256": (
            "d5efd409cf5d5a5a62e9784f8c9e71eebb5a08ee20e8ea7f1d05e894d305a183"
        ),
    }
    assert manifest["updated_theorem_partition"] == {
        "independently_certified_impossible": 119,
        "certified_explicit_ordinary_polytopal": 5,
        "theorem_undecided": 18,
        "total_graph_shadows": 142,
    }
    assert manifest["algebro_geometric_corollary"]["picard_rank"] == 5
    assert all(audit["ordinary_polytopal"] for audit in manifest["audits"])


def test_independent_exact_coordinate_and_face_polynomial_replay() -> None:
    result = replay()
    assert result["ordinary_polytopality_follows_from_facet_wedge_lemma"] is True
    assert [record["record_index"] for record in result["records"]] == [43, 44, 47]
    assert [record["base_facet_count"] for record in result["records"]] == [10, 6, 10]
    assert all(record["exact_nonface_expansion"] for record in result["records"])
    assert all(
        record["residual_face_polynomial"] == [1, 9, 30, 50, 45, 18]
        for record in result["records"]
    )
    assert all(
        record["full_face_polynomial"]
        == [1, 14, 85, 300, 690, 1088, 1190, 880, 405, 90]
        for record in result["records"]
    )
