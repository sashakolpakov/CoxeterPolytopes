from collections import Counter

from scripts.search_15_nodal_quartic_chamber_folds import build_artifacts


def test_official_15_nodal_stock_is_parsed_and_acted_on_exactly() -> None:
    stock, _ = build_artifacts()
    assert stock["lattice"]["rank"] == 16
    assert stock["lattice"]["determinant"] == -128
    assert stock["lattice"]["inertia"] == [1, 15, 0]
    assert stock["lattice"]["ample_square"] == "17/2"
    assert stock["wall_summary"] == {
        "count": 314,
        "kind_histogram": {"inner": 264, "outer": 50},
        "published_orbit_histogram": {
            "1": 10,
            "2": 15,
            "3": 15,
            "4": 10,
            "5": 6,
            "6": 45,
            "7": 6,
            "8": 15,
            "9": 120,
            "10": 72,
        },
        "square_histogram": {
            "-2": 25,
            "-3/2": 6,
            "-1": 66,
            "-3/4": 120,
            "-1/2": 25,
            "-1/4": 72,
        },
    }

    group = stock["chamber_isometry_group"]
    assert group["order"] == 720
    assert group["published_identification"] == "S6"
    assert group["element_order_histogram"] == {
        "1": 1,
        "2": 75,
        "3": 80,
        "4": 180,
        "5": 144,
        "6": 240,
    }
    assert len(group["elements"]) == 720
    assert all(group["checks"].values())
    assert not group["surface_automorphism_realization_used_as_search_gate"]


def test_prime_order_fold_census_and_high_dimensional_exhaustion() -> None:
    _, census = build_artifacts()
    screens = census["subgroups"]
    assert [screen["s6_cycle_type"] for screen in screens] == [
        [2, 1, 1, 1, 1],
        [2, 2, 1, 1],
        [2, 2, 2],
        [3, 1, 1, 1],
        [3, 3],
        [5, 1],
    ]
    assert [screen["hyperbolic_dimension"] for screen in screens] == [
        11,
        9,
        9,
        7,
        5,
        3,
    ]
    assert [screen["wall_orbit_count"] for screen in screens] == [
        178,
        176,
        172,
        114,
        112,
        66,
    ]
    assert [screen["facet_count"] for screen in screens] == [128, 68, 66, 31, 21, 9]
    assert [screen["non_coxeter_pair_count"] for screen in screens] == [
        82,
        143,
        0,
        15,
        0,
        0,
    ]
    assert all(all(screen["checks"].values()) for screen in screens)

    high = [screen for screen in screens if screen["hyperbolic_dimension"] > 9]
    assert len(high) == 1
    assert high[0]["s6_cycle_type"] == [2, 1, 1, 1, 1]
    assert high[0]["non_coxeter_ratio_histogram"] == {
        "1/3": 44,
        "1/9": 6,
        "2/3": 32,
    }
    assert not high[0]["all_facet_pairs_are_coxeter"]
    assert census["dimension_above_nine_exhaustion_theorem"]["conclusion"].startswith(
        "the transposition C2 screen is exhaustive"
    )
    assert census["summary"]["dimension_above_nine_coxeter_count"] == 0


def test_volume_claims_are_separated_from_the_coxeter_gate() -> None:
    _, census = build_artifacts()
    screens = census["subgroups"]
    coxeter = [screen for screen in screens if screen["all_facet_pairs_are_coxeter"]]
    assert Counter(screen["hyperbolic_dimension"] for screen in coxeter) == {
        9: 1,
        5: 1,
        3: 1,
    }

    h9, h5, h3 = coxeter
    assert h9["vinberg"]["status"] == "deferred_lower_dimensional"
    assert h9["vinberg"]["finite_volume"] is None
    assert h5["vinberg"]["status"] == "deferred_lower_dimensional"
    assert h5["vinberg"]["finite_volume"] is None

    assert h3["vinberg"] == {
        "bad_edge_count": 0,
        "compact": False,
        "decision": "finite_volume_noncompact",
        "edge_endpoint_count_histogram": {"2": 19},
        "elliptic_edge_count": 19,
        "finite_vertex_count": 10,
        "finite_volume": True,
        "ideal_vertex_count": 2,
        "status": "completed_exactly",
    }
