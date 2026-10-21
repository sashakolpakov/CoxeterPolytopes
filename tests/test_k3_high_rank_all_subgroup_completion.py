import json

from scripts.certify_k3_high_rank_all_subgroup_completion import OUTPUT, build


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_full_rank15_17_subgroup_archive_replays_exactly() -> None:
    archived = _archive()
    assert build() == archived
    assert archived["all_checks_pass"] is True
    assert archived["summary"]["weighted_group_orders"] == {
        "15": 12,
        "16": 4,
        "17": 2,
    }
    assert archived["summary"]["element_conjugacy_class_counts"] == {
        "15": 6,
        "16": 4,
        "17": 2,
    }
    assert archived["summary"]["subgroup_conjugacy_class_counts"] == {
        "15": 10,
        "16": 5,
        "17": 2,
    }


def test_only_new_high_subgroup_is_same_cone_S3() -> None:
    archived = _archive()
    queue = archived["proper_fixed_rank_at_least_11_queue"]
    assert len(queue) == 4
    assert sorted(
        (row["rho"], row["group_type"], row["fixed_rank"])
        for row in queue
    ) == [
        (15, "C2", 14),
        (15, "C3", 13),
        (15, "S3", 13),
        (16, "C2", 15),
    ]
    s3 = next(row for row in queue if row["group_type"] == "S3")
    assert s3["geometric_status"] == (
        "new_subgroup_class_but_same_fixed_cone_as_C3_triality"
    )
    assert archived["summary"]["new_high_subgroup_classes_with_new_fixed_cone"] == 0
    assert archived["summary"]["new_hyperbolic_Coxeter_polytopes"] == 0


def test_C3_S3_fixed_cone_inherits_strict_and_literal_certificates() -> None:
    equality = _archive()["rank15_C3_S3_same_fixed_cone_certificate"]
    assert all(equality["checks"].values())
    assert len(equality["wall_orbits_one_based"]) == 17
    assert equality["transferred_geometry"] == {
        "hyperbolic_dimension": 12,
        "facet_count": 17,
        "finite_vertex_count": 93,
        "ideal_vertex_count": 3,
        "hyperideal_vertex_count": 0,
        "elliptic_edge_count": 623,
        "endpoint_type_histogram": {
            "0,2": 3,
            "1,1": 124,
            "2,0": 496,
        },
        "volume": "finite_volume_noncompact",
    }
