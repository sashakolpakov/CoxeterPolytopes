import json
import shutil

import pytest
import sympy as sp

from scripts.audit_k3_rank15_order3_h12_scharlau_kirschmer import (
    CERTIFICATE_PATH,
    _canonical_sha256,
    build_certificate,
)


def _archived() -> dict:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_lineage_certificate_is_self_authenticating() -> None:
    payload = _archived()
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256")
    assert digest == _canonical_sha256(body)
    assert all(payload["checks"].values())


def test_exact_integral_lattice_identity() -> None:
    payload = _archived()
    source = sp.Matrix(payload["source_freeze"]["gram12_entry_13x13"])
    candidate_payload = json.loads(
        (
            CERTIFICATE_PATH.parents[1]
            / "candidates/k3_rho15_h12_cyclic_fold_audit.json"
        ).read_text()
    )
    candidate = sp.Matrix(candidate_payload["matrices"]["fixed_lattice_gram_13"])
    isometry = sp.Matrix(payload["exact_lattice_identity"]["integral_isometry_P"])

    assert abs(isometry.det()) == 1
    assert isometry.T * (-candidate) * isometry == source
    assert payload["lattice_invariants"]["genus_symbol"] == (
        "II_{12,1}([2^{1}]_{1} 3^{-1})"
    )


def test_every_face_number_and_cusp_matches() -> None:
    comparison = _archived()["complete_face_and_cusp_comparison"]
    expected = [96, 623, 2054, 4433, 6807, 7680, 6444, 4017, 1835, 596, 130, 17, 1]

    assert comparison["all_face_dimensions_and_cusps_equal"] is True
    assert comparison["source_cusps"] == comparison["candidate_cusps"] == 3
    assert comparison["source_recovered_full_f_vector_including_ideal_vertices"] == expected
    assert comparison["candidate_full_f_vector"] == expected


def test_all_reflective_roots_and_identity_hierarchy() -> None:
    payload = _archived()
    roots = payload["full_reflective_root_audit"]
    hierarchy = payload["identity_hierarchy"]

    assert roots["classification"]["allowed_half_norm_and_divisibility"] == [
        [1, 1],
        [1, 2],
        [3, 3],
    ]
    assert roots["height_zero"]["reflective_root_pairs_by_negative_norm"] == {
        "-2": 91,
        "-6": 1,
    }
    assert roots["height_zero"]["bad_mixed_sign_root_count"] == 0
    assert [
        level["expected_candidate_facets"]
        for level in roots["positive_heights"]["level_audits"]
    ] == [[12], [14], [1], [], [], [16, 17]]
    assert all(
        level["extra_root_query"] == "unsat"
        for level in roots["positive_heights"]["level_audits"]
    )
    assert hierarchy["rung_1_same_genus_and_integral_lattice"]["proved"]
    assert hierarchy["rung_2_same_complete_face_and_cusp_invariants"]["proved"]
    assert not hierarchy["rung_3_explicit_weighted_diagram_facet_permutation"]["proved"]
    assert hierarchy["rung_4_full_reflection_group_identity"]["proved"]
    assert payload["novelty_conclusion"]["new_abstract_hyperbolic_coxeter_polytope"] is False


@pytest.mark.skipif(shutil.which("gp") is None, reason="PARI/GP is required")
def test_full_lineage_certificate_recomputes() -> None:
    assert build_certificate() == _archived()
