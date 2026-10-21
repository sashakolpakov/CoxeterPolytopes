from __future__ import annotations

from scripts.certify_p6_10_exact_known_example import analyze, load_artifact


def test_p6_10_exact_gram_scaling_and_lineage_certificate() -> None:
    report = analyze()
    assert report["all_checks_pass"] is True
    assert report["rank"] == 7
    assert report["kernel_dimension"] == 3
    assert report["inertia_positive_negative_zero"] == [6, 1, 3]
    assert report["mixed_diagonal"] == [2, 2, 2, 2, 2, 2, 4, 4, 2, 2]
    assert report["compactness_status"] == "accepted_known_compact_example"


def test_certificate_does_not_import_disputed_uniqueness() -> None:
    payload = load_artifact()
    boundary = payload["classification_boundary"]
    provenance = payload["compactness_provenance"]
    assert boundary["arXiv_2608_02894v1_uniqueness"] == (
        "unaccepted_open_due_to_invalid_facet_filter"
    )
    assert boundary["worldwide_new_abstract_polytope_claim"] is False
    assert provenance["local_certificate_scope"] == (
        "Gram, algebraic scaling, and lineage consistency only"
    )
    assert any("vertex enumeration" in item for item in provenance["not_claimed_by_this_artifact"])
