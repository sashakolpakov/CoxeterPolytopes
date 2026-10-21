from functools import lru_cache
import json

from scripts.certify_k3_rank18_19_all_subgroup_completion import OUTPUT, build


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


@lru_cache(maxsize=1)
def _fresh() -> dict:
    return build()


def test_rank18_19_all_subgroup_archive_replays_exactly() -> None:
    archived = _archive()
    assert _fresh() == archived
    assert archived["all_checks_pass"] is True
    assert archived["summary"]["weighted_group_orders"] == {"18": 2, "19": 6}
    assert archived["summary"]["element_conjugacy_class_counts"] == {
        "18": 2,
        "19": 3,
    }
    assert archived["summary"]["subgroup_conjugacy_class_counts"] == {
        "18": 2,
        "19": 4,
    }


def test_only_high_class_is_the_frozen_rank19_transposition() -> None:
    archived = _archive()
    queue = archived["proper_fixed_rank_at_least_11_queue"]
    assert len(queue) == 1
    assert (
        queue[0]["rho"],
        queue[0]["group_type"],
        queue[0]["conjugacy_class_size"],
        queue[0]["fixed_rank"],
        queue[0]["hyperbolic_dimension"],
    ) == (19, "C2", 3, 11, 10)
    assert archived["summary"]["new_noninvolution_high_classes"] == 0
    assert archived["summary"]["new_hyperbolic_Coxeter_polytopes"] == 0
    assert archived["summary"]["new_algebraic_realizations"] == 0


def test_survivor_has_strict_cone_and_literal_vinberg_certificates() -> None:
    survivor = _archive()["survivor_geometry"]
    assert all(survivor["checks"].values())
    assert len(survivor["strict_facet_audits"]) == 13
    assert all(
        row["irredundant_with_strict_relative_interior"]
        for row in survivor["strict_facet_audits"]
    )
    assert survivor["coxeter_edge_counts"] == {
        "3": 9,
        "4": 2,
        "infinity": 2,
    }
    assert (
        survivor["exact_cone"]["finite_vertex_count"],
        survivor["exact_cone"]["ideal_vertex_count"],
        survivor["exact_cone"]["hyperideal_vertex_count"],
    ) == (34, 4, 0)
    literal = survivor["literal_coxiter_prop_1_10"]
    assert literal["spherical_rank_n_minus_one_count"] == 208
    assert literal["extension_histogram"] == {"2": 208}
    assert literal["endpoint_type_histogram"] == {
        "0,2": 5,
        "1,1": 66,
        "2,0": 137,
    }
    assert literal["bad_edges"] == []
    assert survivor["finite_volume"] is True
    assert survivor["compact"] is False


def test_fixed_lattice_and_exact_deduplication() -> None:
    archived = _archive()
    fixed = archived["survivor_geometry"]["fixed_lattice"]
    assert (fixed["rank"], fixed["determinant"], fixed["inertia"]) == (
        11,
        512,
        [1, 10, 0],
    )
    assert fixed["smith_diagonal_absolute"] == [1, 1, *([2] * 9)]
    assert fixed["facet_mirrors_generate_saturated_fixed_lattice"] is True
    assert all(
        row["reflection_integral"]
        for row in fixed["mirror_records"]
        if row["role"] == "facet"
    )

    dedup = archived["exact_deduplication"]
    assert dedup["frozen_published_rank19_fold"][
        "entrywise_exact_facet_gram_match"
    ] is True
    assert len(dedup["three_realized_S3_conjugates"]) == 3
    assert all(
        row["entrywise_exact_facet_gram_match"]
        for row in dedup["three_realized_S3_conjugates"]
    )
    assert dedup["classical_f2_H10_13"]["isomorphic"] is False
    corpus = dedup["pinned_coxiter_H10_13"]
    assert len(corpus["ordinary_H10_13_records"]) == 4
    assert len(corpus["doubling_H10_13_records"]) == 4
    assert corpus["ordinary_exact_matches"] == []
    assert corpus["doubling_exact_matches"] == []
