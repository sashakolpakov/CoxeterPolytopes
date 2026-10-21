from __future__ import annotations

import hashlib
import json

import pytest

from scripts.certify_coble_type_vi_h14_single_doublings import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    return build_certificate()


def test_archive_is_exact_replay(certificate: dict) -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert certificate == archived
    assert certificate["all_checks_pass"] is True


def test_two_exact_weighted_gram_orbits(certificate: dict) -> None:
    seed = certificate["seed"]
    assert seed["weighted_gram_automorphism_group_order"] == 120
    assert sorted(map(len, seed["doubling_wall_orbits_one_based"])) == [15, 20]
    assert set(seed["doubling_walls_one_based"]) == set(range(26, 61))
    assert set(certificate["positive_doubles"]) == {"H14/100", "H14/101"}


@pytest.mark.parametrize(
    ("name", "facets", "finite", "ideal", "orbit", "mirror_types"),
    [
        ("H14/100", 100, 54_318, 1_993, 15, {"-2,1": 18, "-2,2": 68, "-4,2": 14}),
        ("H14/101", 101, 55_130, 2_043, 20, {"-2,1": 20, "-2,2": 68, "-4,2": 13}),
    ],
)
def test_complete_positive_chambers(
    certificate: dict,
    name: str,
    facets: int,
    finite: int,
    ideal: int,
    orbit: int,
    mirror_types: dict[str, int],
) -> None:
    record = certificate["positive_doubles"][name]
    construction = record["construction"]
    assert record["seed_automorphism_orbit_size"] == orbit
    assert record["reflection_subgroup_index"] == 2
    assert record["finite_volume"] is True and record["compact"] is False
    assert record["ray_census"]["finite"] == finite
    assert record["ray_census"]["ideal"] == ideal
    assert record["ray_census"]["hyperideal"] == 0
    assert construction["facet_count"] == facets
    assert construction["wall_gram_inertia_by_full_span"] == [1, 14, facets - 15]
    assert construction["wall_span_index_one"] is True
    assert construction["all_roots_primitive_and_reflections_integral"] is True
    assert construction["all_facets_irredundant"] is True
    assert construction["invalid_finite_angle_count"] == 0
    assert construction["primitive_square_divisibility_histogram"] == mirror_types

    physical = construction["physical_NS_divisor_carrier"]
    assert physical["NS_identification"] == "I_(1,14)"
    assert physical["glued_seed_divisor_square"] == -1
    assert physical["glued_divisor_reflection_integral_on_NS"] is True
    assert physical["wall_span_rank"] == 15
    assert physical["wall_span_index_one"] is True
    assert physical["normalized_gram_equals_arithmetic_normalized_gram_entrywise"] is True
    assert len(physical["wall_vectors_in_NS_basis"]) == facets
    assert all(len(vector) == 15 for vector in physical["wall_vectors_in_NS_basis"])

    vectors = construction["wall_vectors_in_arithmetic_basis"]
    raw_gram = construction["wall_gram"]
    normalized = construction["exact_normalized_squared_gram"]
    assert len(vectors) == len(raw_gram) == len(normalized) == facets
    assert all(len(row) == 15 for row in vectors)
    assert all(len(row) == facets for row in raw_gram)
    assert all(len(row) == facets for row in normalized)
    digest = hashlib.sha256(
        json.dumps(normalized, separators=(",", ":")).encode()
    ).hexdigest()
    assert digest == construction["exact_normalized_squared_gram_sha256"]


def test_exact_redoublable_pair_and_infinite_family(certificate: dict) -> None:
    expected = {"H14/100": 316, "H14/101": 324}
    for name, pairs in expected.items():
        frontier = certificate["positive_doubles"][name]["construction"][
            "redoubling_frontier"
        ]
        assert frontier["doubling_wall_count"] == 68
        assert frontier["parallel_pairs_among_doubling_walls"] == pairs
        pair = frontier["explicit_redoublable_pair"]
        assert pair["squares"] == [-2, -2]
        assert pair["raw_inner_product"] == 2
        assert pair["normalized_squared_ratio"] == "1"
        assert pair["sum_isotropic"] is True
        assert pair["classification"] == "parallel_asymptotic_disjoint_in_H14"
        family = frontier["allcock_infinite_family"]
        assert family["index_in_this_double_reflection_group"] == "N"
        assert family["index_in_parent_H14_60_reflection_group"] == "2*N"
        assert family["volume_multiple_of_parent_H14_60"] == "2*N"
        assert family["facet_formula"] == "not claimed without a separate collision recurrence"


def test_type_mii_is_exact_negative_factory_control(certificate: dict) -> None:
    control = certificate["type_MII_H17_negative_doubling_control"]
    assert control["doubling_wall_count"] == 0
    assert control["m3_neighbors_per_wall_histogram"] == {"2": 16, "4": 8, "8": 24}
    assert all(control["checks"].values())


def test_lineage_and_surface_scope_are_strict(certificate: dict) -> None:
    assert certificate["lineage"]["worldwide_priority"] == "open"
    assert "not nef cones" in certificate["claim"]["scope"]
    assert "not asserted effective" in certificate["divisor_carrier_scope"]["negative"]
    for header in ("H14/100", "H14/101"):
        counts = certificate["lineage"]["bounded_corpus"]["headers"][header]
        assert counts == {"ordinary": 0, "doubling": 0}
