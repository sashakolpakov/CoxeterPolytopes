import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.k3 import K3_WALL_NORMS
from coxeter_polytopes.known_examples.bugaenko import BUGAENKO_VERTICES
from scripts.bugaenko_k3_eigen_nef import (
    analyze,
    field_form_is_even,
    mirrors_meeting_closed_chamber,
    spherical_mirrors_at_vertex,
    surviving_even_k3_shells,
    trace_roots_meeting_closed_chamber,
    verify_field_reflection_factorization,
    verify_ns_phi_endomorphism,
)


ARCHIVE = (
    Path(__file__).parents[1]
    / "data"
    / "trace_lattices"
    / "bugaenko_k3_eigen_nef.json"
)


def test_evenness_reduces_five_K3_shells_to_one_mirror_arrangement() -> None:
    assert field_form_is_even() is True
    assert surviving_even_k3_shells() == (Golden(4, -2), Golden(2))
    assert PHI * PHI * Golden(4, -2) == 2
    assert (PHI - 1) ** 2 * 2 == Golden(4, -2)
    factorization = verify_field_reflection_factorization()
    assert all(factorization.values())


def test_every_mirror_touching_the_closed_chamber_is_enumerated_locally() -> None:
    local_counts = [
        len(spherical_mirrors_at_vertex(vertex)) for vertex in BUGAENKO_VERTICES
    ]
    mirrors = mirrors_meeting_closed_chamber()

    assert len(mirrors) == 320
    assert sum(local_counts) == 1856
    assert {count: local_counts.count(count) for count in set(local_counts)} == {
        20: 2,
        21: 2,
        23: 2,
        24: 2,
        26: 8,
        28: 2,
        32: 4,
        35: 2,
        41: 2,
        43: 2,
        64: 2,
        67: 2,
        68: 2,
        70: 2,
        76: 2,
        120: 3,
    }


def test_all_signed_K3_vectors_in_the_five_shells_are_accounted_for() -> None:
    roots = trace_roots_meeting_closed_chamber()
    assert set(roots) == set(K3_WALL_NORMS)
    assert {str(shell): len(roots[shell]) for shell in K3_WALL_NORMS} == {
        "5 - 3*phi": 0,
        "4 - 2*phi": 640,
        "3 - 1*phi": 0,
        "2": 640,
        "1 + 1*phi": 0,
    }
    assert len(set(roots[Golden(2)])) == 640
    assert len(set(roots[Golden(4, -2)])) == 640


def test_NS_structure_is_self_adjoint_endomorphism_not_isometry() -> None:
    checks = verify_ns_phi_endomorphism()
    assert checks == {
        "phi_polynomial_verified": True,
        "phi_self_adjoint_verified": True,
        "phi_is_quadratic_form_isometry": False,
    }


def test_archived_eigen_nef_completion_is_reproducible() -> None:
    result = analyze()
    assert result == json.loads(ARCHIVE.read_text())
    assert result["closed_chamber"] == {
        "distinct_k3_mirrors_meeting": 320,
        "mirror_representative_sha256": (
            "815c0c4cad8e4ab9da32ef491c104ac71a781d4f2b4300f988ecb484506da69e"
        ),
        "facet_mirrors": 11,
        "nonfacet_mirrors_touching_lower_faces": 309,
        "mirrors_cutting_interior": 0,
        "contact_dimension_distribution": {
            "0": 88,
            "1": 43,
            "2": 26,
            "3": 17,
            "4": 96,
            "5": 25,
            "6": 14,
            "7": 11,
        },
    }
    realization = result["eigen_nef_realization"]
    assert realization["extra_nonreflection_facets"] == 0
    assert realization["restricted_chamber_equals_bugaenko_chamber"] is True
    assert realization["explicit_K3_equation_or_Weyl_word_archived"] is False
