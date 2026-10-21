import json

import pytest

from coxeter_polytopes.combinatorics.h9_14_sparse_coxeter_labels import (
    audit_all_sparse_h9_14_labels,
    audit_sparse_h9_14_labels,
)
from scripts.certify_h9_14_sparse_coxeter_labels import OUTPUT, verify
from scripts.certify_h9_14_ordinary_coxeter_label_gate import (
    OUTPUT as UNIFIED_OUTPUT,
    verify as verify_unified,
)


def test_complete_sparse_label_censuses_and_symmetry_orbits() -> None:
    audits = audit_all_sparse_h9_14_labels()
    assert [audit.record_index for audit in audits] == [43, 44, 47]
    assert [audit.isolated_simplex_lanner_face for audit in audits] == [
        (3, 4, 5, 6, 7),
        (2, 3, 4, 5, 6),
        (3, 4, 5, 6, 7),
    ]
    assert [audit.cluster_join_counts for audit in audits] == [
        (420, 1332, 4212),
        (420, 1332, 4620),
        (420, 1332, 222),
    ]
    assert [audit.passing_cluster_assignment_count for audit in audits] == [
        30,
        144,
        18,
    ]
    assert [audit.full_labelled_assignment_count for audit in audits] == [
        12_600,
        60_480,
        7_560,
    ]
    assert [audit.symmetry_orbit_count for audit in audits] == [15, 25, 5]
    assert [audit.symmetry_orbit_size_histogram for audit in audits] == [
        ((360, 3), (720, 8), (1440, 4)),
        ((1440, 12), (2880, 11), (5760, 2)),
        ((1080, 3), (2160, 2)),
    ]
    assert all(audit.verified for audit in audits)


def test_sparse_orbit_representatives_are_exportable_complete_label_maps() -> None:
    for audit in audit_all_sparse_h9_14_labels():
        representatives = audit.orbit_representative_label_maps()
        assert len(representatives) == audit.symmetry_orbit_count
        assert all(len(labels) == 85 for labels in representatives)
        assert len({
            tuple(sorted(labels.items())) for labels in representatives
        }) == audit.symmetry_orbit_count
        assert audit.representative_all_ninety_vertices_finite is True
        assert audit.representative_all_four_missing_fives_lanner is True
        assert audit.orthogonal_negative_index_lower_bound == 2
        assert len(audit.minimal_nonface_components) == 2
        assert audit.unrestricted_minimum_cross_component_nonright_edges == 1


def test_sparse_census_rejects_an_unsupported_record() -> None:
    with pytest.raises(ValueError, match="records 43, 44, and 47"):
        audit_sparse_h9_14_labels(30)


def test_frozen_sparse_label_certificate_replays() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_indices": [43, 44, 47],
        "full_labelled_solutions": 80_640,
        "symmetry_orbits": 45,
        "certificate_sha256": manifest["certificate_sha256"],
    }
    assert manifest["scope"]["discrete_sparse_label_status"] == "sat"
    assert manifest["scope"]["global_sparse_gram_status"] == "unsat"


def test_unified_five_record_label_gate_replays() -> None:
    manifest = json.loads(UNIFIED_OUTPUT.read_text(encoding="utf-8"))
    assert verify_unified(manifest) == {
        "record_indices": [30, 42, 43, 44, 47],
        "sparse_labelled_solutions": 46_540_944,
        "full_nerve_symmetry_orbits": 6_210,
        "certificate_sha256": (
            "9c3943ea538733020bc09f34c0150cdf6a7d66af743d5643a32ac6033057f4ea"
        ),
    }
    assert all(
        record["local_discrete_status"] == "sat"
        and record["sparse_global_gram_status"] == "unsat"
        for record in manifest["records"]
    )
    assert [
        record["unrestricted_minimum_cross_component_nonright_edges"]
        for record in manifest["records"]
    ] == [4, 2, 1, 1, 1]
