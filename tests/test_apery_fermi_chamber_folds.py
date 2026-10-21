from collections import Counter

from scripts.search_apery_fermi_chamber_folds import build_artifacts


def test_apery_fermi_stock_is_parsed_exactly() -> None:
    stock, _ = build_artifacts()
    assert stock["lattice"]["rank"] == 19
    assert stock["lattice"]["determinant"] == 12
    assert stock["lattice"]["inertia"] == [1, 18, 0]
    assert stock["lattice"]["ample_square"] == "32"
    assert stock["wall_summary"] == {
        "count": 80,
        "kind_histogram": {"inner": 56, "outer": 24},
        "square_histogram": {
            "-2": 24,
            "-4/3": 4,
            "-1": 8,
            "-3/4": 24,
            "-1/3": 4,
            "-1/12": 16,
        },
    }
    stabilizer = stock["weighted_wall_symmetry_group"]
    assert stabilizer["order"] == 32
    assert stabilizer["element_order_histogram"] == {
        "1": 1,
        "2": 15,
        "4": 8,
        "8": 8,
    }
    assert len({tuple(permutation) for permutation in stabilizer["wall_permutations"]}) == 32
    assert stabilizer[
        "equals_full_abstract_weighted_wall_graph_automorphism_group"
    ]
    assert stock["geometric_chamber_subgroup"] == {
        "index_in_weighted_wall_symmetry_group": 2,
        "order": 16,
        "published_source_assignment": "AutXD0",
        "used_as_a_search_gate": False,
    }


def test_all_chamber_stabilizer_subgroups_are_screened() -> None:
    _, census = build_artifacts()
    assert census["subgroup_count"] == 58
    assert census["subgroup_order_histogram"] == {
        "1": 1,
        "2": 15,
        "4": 19,
        "8": 15,
        "16": 7,
        "32": 1,
    }
    assert census["summary"]["all_exact_reduction_checks_pass"]
    assert census["summary"]["genuine_h2_through_h8_coxeter_count"] == 33
    assert (
        census["summary"]["genuine_h2_through_h8_vinberg_finite_volume_count"]
        == 33
    )
    assert census["summary"]["degenerate_h1_coxeter_compatible_count"] == 8
    assert census["summary"]["vinberg_was_not_run_in_h1"]

    nontrivial_high = [
        subgroup
        for subgroup in census["subgroups"]
        if subgroup["order"] > 1 and subgroup["hyperbolic_dimension"] > 9
    ]
    assert len(nontrivial_high) == 11
    assert Counter(
        subgroup["hyperbolic_dimension"] for subgroup in nontrivial_high
    ) == {11: 6, 10: 5}
    assert Counter(subgroup["facet_count"] for subgroup in nontrivial_high) == {
        22: 6,
        28: 5,
    }
    assert sorted(subgroup["non_coxeter_pair_count"] for subgroup in nontrivial_high) == [
        2,
        2,
        2,
        2,
        4,
        8,
        8,
        10,
        10,
        10,
        10,
    ]
    assert not any(subgroup["all_facet_pairs_are_coxeter"] for subgroup in nontrivial_high)

    coxeter_dimensions = Counter(
        subgroup["hyperbolic_dimension"]
        for subgroup in census["subgroups"]
        if subgroup["all_facet_pairs_are_coxeter"]
    )
    assert coxeter_dimensions[8] == 4
    assert all(
        subgroup["facet_count"] == 14
        for subgroup in census["subgroups"]
        if subgroup["hyperbolic_dimension"] == 8
    )
