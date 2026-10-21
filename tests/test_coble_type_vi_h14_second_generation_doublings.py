from __future__ import annotations

from collections import Counter
import json

import pytest

from scripts.certify_coble_type_vi_h14_second_generation_doublings import (
    OUTPUT,
    PINNED,
    SEED_RAYS,
    SINGLE,
    SOURCE,
    _sha256,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_and_pinned_lineage(certificate: dict) -> None:
    assert certificate["all_checks_pass"] is True
    assert certificate["schema"].endswith("second_generation_doublings.v1")
    assert {
        "single": _sha256(SINGLE),
        "source": _sha256(SOURCE),
        "seed_rays": _sha256(SEED_RAYS),
    } == PINNED
    assert all(certificate["checks"].values())


def test_seventeen_intrinsic_orbits_are_exhaustive(certificate: dict) -> None:
    classification = certificate["classification"]
    assert classification["parent_intrinsic_automorphism_orders"] == {
        "H14/100": 16,
        "H14/101": 12,
    }
    assert classification["intrinsic_doubling_orbit_counts"] == {
        "H14/100": 7,
        "H14/101": 10,
    }
    cases = certificate["cases"]
    assert len(cases) == 17
    assert Counter(case["parent"] for case in cases) == {
        "H14/100": 7,
        "H14/101": 10,
    }
    for parent in ("H14/100", "H14/101"):
        orbits = [
            set(case["doubling_orbit_one_based"])
            for case in cases
            if case["parent"] == parent
        ]
        assert sum(map(len, orbits)) == 68
        assert len(set().union(*orbits)) == 68
        assert all(left.isdisjoint(right) for i, left in enumerate(orbits) for right in orbits[i + 1 :])


def test_pairwise_abstract_nonisometry_scope(certificate: dict) -> None:
    cases = certificate["cases"]
    hashes = [
        case["abstract_isometry_invariant_sha256"]
        for case in cases
    ]
    assert len(set(hashes)) == 17
    assert "uses no lattice root-length decoration" in certificate[
        "classification"
    ]["pairwise_nonisometry_proof"]
    assert certificate["novelty"] == {
        "new_to_project_abstract_polytopes": 17,
        "worldwide_priority": "open",
        "reason": (
            "Allcock already gives infinite H14 families; no comprehensive "
            "publication-level isometry census was available"
        ),
    }


@pytest.mark.parametrize(
    ("facets", "multiplicity", "finite", "ideal", "edges", "endpoint_types"),
    [
        (180, 3, 107_264, 3_925, 940_516, {"0,2": 37_134, "1,1": 305_068, "2,0": 598_314}),
        (181, 4, 108_076, 3_975, 948_264, {"0,2": 37_626, "1,1": 308_212, "2,0": 602_426}),
        (182, 4, 108_888, 4_025, 956_012, {"0,2": 38_118, "1,1": 311_356, "2,0": 606_538}),
        (183, 6, 109_700, 4_075, 963_760, {"0,2": 38_610, "1,1": 314_500, "2,0": 610_650}),
    ],
)
def test_complete_ray_and_literal_vinberg_profiles(
    certificate: dict,
    facets: int,
    multiplicity: int,
    finite: int,
    ideal: int,
    edges: int,
    endpoint_types: dict[str, int],
) -> None:
    assert sum(case["facet_count"] == facets for case in certificate["cases"]) == multiplicity
    profile = certificate["audit_profiles_by_facet_count"][f"F{facets}"]
    rays = profile["ray_census"]
    literal = profile["literal_vinberg_endpoint_cusp_audit"]
    assert rays["finite"] == finite
    assert rays["ideal"] == ideal
    assert rays["hyperideal"] == 0
    assert literal["elliptic_rank_13_edge_count"] == edges
    assert literal["endpoint_count_histogram"] == {"2": edges}
    assert literal["endpoint_type_histogram"] == endpoint_types
    assert literal["bad_edge_count"] == 0
    assert literal["invalid_full_parabolic_neighborhood_count"] == 0
    assert literal["full_zero_neighborhood_failure_count"] == 0
    assert literal["decision"] == "finite_volume_noncompact"
    assert "independently complete union/ray transport" in literal[
        "extreme_ray_dependency"
    ]
    assert "no-ray-import" in literal["extreme_ray_dependency"]


def test_divisor_carrier_irredundancy_and_redoubling(certificate: dict) -> None:
    classification = certificate["classification"]
    assert classification["ambient_lattice_inertia"] == [1, 14, 0]
    assert classification["every_arithmetic_and_physical_wall_span"] == {
        "rank": 15,
        "smith_diagonal": [1] * 15,
    }
    for case in certificate["cases"]:
        assert case["wall_gram_inertia_by_full_span"] == [
            1,
            14,
            case["facet_count"] - 15,
        ]
        assert case["minimum_facet_witness_slack"] == "2"
        assert case["all_checks_pass"] is True
        profile = certificate["audit_profiles_by_facet_count"][
            case["audit_profile"]
        ]
        assert profile["redoubling"]["doubling_wall_count"] == 134
        assert profile["redoubling"][
            "parallel_pairs_among_doubling_walls"
        ] > 0


def test_union_ray_lemma_and_algebraic_claim_boundary(certificate: dict) -> None:
    lemma = certificate["polyhedral_union_ray_transport_lemma"]
    assert "No other extreme rays occur" in lemma["statement"]
    assert "s_h(a)-a" in lemma["proof"]
    assert "causal" in lemma["causal_type_consequence"] or "type" in lemma["causal_type_consequence"]
    scope = certificate["algebraic_scope"]
    assert "primitive integral physical NS divisor class" in scope["positive"]
    assert "effectivity" in scope["not_claimed"]
    assert "nef-cone equality" in scope["not_claimed"]
    assert "surface-automorphism" in scope["not_claimed"]
