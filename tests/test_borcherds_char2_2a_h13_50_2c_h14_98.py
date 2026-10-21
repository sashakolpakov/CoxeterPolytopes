import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_construction.json"
)
VINBERG = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_vinberg.json"
)


def _cases(path: Path) -> dict[str, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    return {record["id"]: record for record in payload["cases"]}


def test_h13_h14_exact_construction_archive() -> None:
    cases = _cases(CONSTRUCTION)
    assert set(cases) == {"SS2-2A-H13-50", "SS2-2C-H14-98"}
    h13 = cases["SS2-2A-H13-50"]
    assert h13["fixed_lattice"]["determinant"] == -1024
    assert h13["fixed_lattice"]["smith_diagonal_absolute"] == [1] * 4 + [2] * 10
    assert h13["restricted_chamber"]["primitive_facet_type_histogram"] == {
        "orbit1_square-2_div1": 10,
        "orbit1_square-4_div2": 24,
        "orbit2_square-4_div2": 16,
    }
    assert h13["restricted_chamber"]["facet_gram_sha256"] == (
        "7aaea6ed54d34abd3e6d205c1ba132b243bfddcc0ad2db455a1a45a28a21bd51"
    )

    h14 = cases["SS2-2C-H14-98"]
    assert h14["fixed_lattice"]["determinant"] == 256
    assert h14["fixed_lattice"]["smith_diagonal_absolute"] == (
        [1] * 8 + [2] * 6 + [4]
    )
    assert h14["restricted_chamber"]["primitive_facet_type_histogram"] == {
        "orbit1_square-2_div1": 14,
        "orbit1_square-4_div2": 14,
        "orbit2_square-4_div2": 14,
        "orbit2_square-4_div4": 56,
    }
    assert h14["restricted_chamber"]["facet_gram_sha256"] == (
        "257c0444e7774232672520096ef5639074d83ef2e1a58d625ab2435007f2f330"
    )


def test_h13_literal_vinberg_archive() -> None:
    h13 = _cases(VINBERG)["SS2-2A-H13-50"]
    assert h13["finite_volume"] is True
    assert h13["compact"] is False
    assert h13["elliptic_edge_count"] == 391614
    assert h13["finite_vertex_count"] == 45954
    assert h13["ideal_vertex_count"] == 2598
    assert h13["connected_affine_component_count"] == 6210
    assert h13["edge_endpoint_type_histogram"] == {
        "0,2": 23340,
        "1,1": 139146,
        "2,0": 229128,
    }
    assert h13["bad_edge_count"] == 0
    assert h13["bad_parabolic_neighborhood_count"] == 0


def test_h14_literal_vinberg_and_alexeev_counts() -> None:
    h14 = _cases(VINBERG)["SS2-2C-H14-98"]
    assert h14["finite_volume"] is True
    assert h14["compact"] is False
    assert h14["finite_special_subset_counts_by_rank"] == [
        98,
        1666,
        14616,
        83559,
        342888,
        1056573,
        2490674,
        4501714,
        6177500,
        6299076,
        4598006,
        2249968,
        646786,
        75898,
    ]
    assert h14["ideal_vertex_count"] == 2114
    assert h14["connected_affine_component_count"] == 5838
    assert h14["edge_endpoint_type_histogram"] == {
        "0,2": 19628,
        "1,1": 191744,
        "2,0": 435414,
    }
    assert h14["bad_edge_count"] == 0
    assert h14["bad_parabolic_neighborhood_count"] == 0
