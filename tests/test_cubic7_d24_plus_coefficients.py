import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7, THETA7
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    DELTA7,
    P7_D24_SEED,
    TWO_P7_D24_SEED,
    TWO_THETA_D24_SEED,
    d24_coefficient_pair,
    enumerate_coefficient_control_levels,
    enumerate_coefficient_level,
    unit_coefficient_isometry_audit,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/cubic7_d24_plus_coefficient_screen_h100.json"
)
P7_H300_ARTIFACT = ROOT / "data/enumerations/cubic7_p7_d24_plus_h300.json"


def test_unit_squareclass_is_deduplicated_and_norm_classes_are_complete() -> None:
    assert unit_coefficient_isometry_audit()
    assert TWO_THETA_D24_SEED.canonical_mirror_norms == (Cubic7(2),)
    assert P7_D24_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )
    assert TWO_P7_D24_SEED.canonical_mirror_norms == (Cubic7(2), 2 * DELTA7)
    assert all(
        seed.admissible
        for seed in (TWO_THETA_D24_SEED, P7_D24_SEED, TWO_P7_D24_SEED)
    )


def test_three_first_packets_are_exact_affine_d24_walls() -> None:
    expected = {
        TWO_THETA_D24_SEED: (
            10,
            2 * THETA7 + THETA7**2,
            Cubic7(2),
            1 + THETA7,
        ),
        P7_D24_SEED: (
            50,
            7 + 22 * THETA7 + 10 * THETA7**2,
            Cubic7(1),
            1 + 2 * THETA7 + THETA7**2,
        ),
        TWO_P7_D24_SEED: (
            20,
            2 + 7 * THETA7 + 3 * THETA7**2,
            Cubic7(2),
            2 * THETA7 + THETA7**2,
        ),
    }
    for seed, (bound, height, norm, spatial_value) in expected.items():
        levels = enumerate_coefficient_control_levels(seed, bound)
        level = next(item for item in levels if item.height == height)
        roots = enumerate_coefficient_level(seed, level)
        assert level.height == height
        assert len(roots) == 1
        root = roots[0]
        assert root.norm == norm and root.coset == "D"
        assert root.vector[1:3] == (spatial_value, spatial_value)
        assert all(not value for value in root.vector[3:])
        assert d24_coefficient_pair(seed, root.vector, root.vector) == norm


def test_frozen_h100_coefficient_ranking_and_ray_hashes() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["best_exact_survivor"] == P7_D24_SEED.seed_id
    records = {record["seed_id"]: record for record in payload["seeds"]}
    expected = {
        TWO_THETA_D24_SEED.seed_id: (
            5,
            10,
            [(25, 25, 2, 0, 23), (27, 113, 4, 0, 109)],
            "03aa42cdc2dd3b4df2126ef3e08c40ae791c2e4251b326aa0137cb36c20b9fae",
        ),
        P7_D24_SEED.seed_id: (
            31,
            1,
            [(25, 25, 2, 0, 23)],
            "d6feb8d1ddfac35443db8332bcfc06f55e93046e63905ee989b78944fbb8cc15",
        ),
        TWO_P7_D24_SEED.seed_id: (
            8,
            8,
            [
                (25, 25, 2, 0, 23),
                (27, 113, 4, 0, 109),
                (32, 404, 6, 0, 398),
            ],
            "9956beaa6cde790a7e31d03ffab68d9f76944cd74e22454e9df2d0c647f72732",
        ),
    }
    for seed_id, (level_count, accepted_count, snapshots, digest) in expected.items():
        record = records[seed_id]
        assert record["control_level_count"] == level_count
        assert record["accepted_root_count"] == accepted_count
        assert record["spinor_accepted_count"] == 0
        assert [
            (
                item["wall_count"],
                item["extreme_ray_count"],
                item["finite_vertex_count"],
                item["cusp_count"],
                item["hyperideal_ray_count"],
            )
            for item in record["packet_snapshots"]
        ] == snapshots
        assert record["geometry_final"]["ray_stream_sha256"] == digest
        assert not record["geometry_final"]["finite_volume"]


def test_frozen_p7_h300_prefix_is_complete_and_still_open() -> None:
    payload = json.loads(P7_H300_ARTIFACT.read_text())
    assert payload["seed_id"] == P7_D24_SEED.seed_id
    assert payload["height_bound"] == 300
    assert payload["canonical_mirror_norms"] == [
        {
            "display": str(norm),
            "exact": norm.to_json(),
            "signs": [norm.sign(place) for place in range(3)],
        }
        for norm in (Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7)
    ]
    levels = enumerate_coefficient_control_levels(P7_D24_SEED, 300)
    assert len(levels) == payload["control_level_count"] == 51
    assert sum(
        packet["candidate_count"] for packet in payload["height_packets"]
    ) == 1
    assert sum(
        packet["accepted_count"] for packet in payload["height_packets"]
    ) == payload["accepted_root_count"] == 1
    assert payload["spinor_accepted_count"] == 0
    assert payload["same_height_audit"]["passed"]
    assert payload["frontier"][
        "all_controls_at_height_at_most_300_exhausted"
    ]
    assert payload["geometry"] == {
        "all_active_sets_rank_24": True,
        "compact": False,
        "cusp_count": 0,
        "engine": (
            "exact double description; adjacency by common-active rank 23"
        ),
        "extreme_ray_count": 25,
        "finite_vertex_count": 2,
        "finite_volume": False,
        "hyperideal_ray_count": 23,
        "ray_stream_sha256": (
            "d6feb8d1ddfac35443db8332bcfc06f55e93046e63905ee989b78944fbb8cc15"
        ),
        "wall_count": 25,
    }
