import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LINEAGE = ROOT / "data/enumerations/borcherds_char2_2d_h10_47_lineage.json"


def _load() -> dict:
    return json.loads(LINEAGE.read_text(encoding="utf-8"))


def test_explicit_f2_10_integral_isometry() -> None:
    record = _load()
    lattice = record["lattice_isometry"]
    assert record["all_checks_pass"] is True
    assert lattice["target_model"] == "f_2^10 = <2> + <-1>^10"
    assert lattice["target_gram"] == [
        [2 if row == 0 and column == 0 else
         -1 if row == column and row > 0 else 0
         for column in range(11)]
        for row in range(11)
    ]
    assert lattice["basis_change_determinant"] == -1
    assert lattice["half_gram_smith_diagonal_absolute"] == [1] * 10 + [2]
    assert all(lattice["checks"].values())


def test_SK_entry_2_and_maximal_chamber() -> None:
    record = _load()
    source = record["scharlau_kirschmer_source"]
    maximal = record["maximal_chamber_certificate"]
    assert source["entry_one_based"] == 2
    assert source["list_entry"] == [
        2,
        "I_{10,1}([1^{10}2^{1}]_{1})",
        [13, 76, 264, 606, 966, 1092, 876, 488, 175, 29],
        3,
    ]
    assert maximal["facet_count"] == 13
    assert maximal["finite_vertex_count"] == 29
    assert maximal["ideal_vertex_count"] == 3
    assert maximal["vinberg_stopping_fired"] is True


def test_complete_solomon_index_and_volume() -> None:
    record = _load()
    solomon = record["solomon_index"]
    subgroup = solomon["subgroup_47_walls"]
    full = solomon["full_f2_10_group_13_walls"]
    assert subgroup["finite_special_subset_count_by_rank"] == [
        1, 47, 540, 3230, 12215, 31212,
        54582, 64050, 47760, 19910, 3136,
    ]
    assert subgroup["connected_finite_counts_by_rank_1_through_10"] == [
        47, 50, 120, 240, 300, 360, 420, 420, 240, 120,
    ]
    assert subgroup["rational_euler_characteristic"] == "-527/120960"
    assert full["finite_special_subset_count_by_rank"] == [
        1, 13, 76, 264, 606, 966, 1092, 876, 488, 175, 29,
    ]
    assert full["rational_euler_characteristic"] == "-31/11147673600"
    assert solomon["index"] == 1_566_720
    assert solomon["index_prime_factorization"] == {
        "2": 11, "3": 2, "5": 1, "17": 1
    }
    assert solomon["hyperbolic_covolume"] == "pi^5 * 527/3572100"
    assert solomon["full_face_vector"] == [
        3352, 19910, 47760, 64050, 54582, 31212,
        12215, 3230, 540, 47, 1,
    ]
    assert all(solomon["checks"].values())


def test_lineage_decision_keeps_literal_priority_open() -> None:
    decision = _load()["decision"]
    assert decision == {
        "ambient_lattice_and_full_reflection_group": (
            "known: Vinberg f_2^10; Scharlau--Kirschmer entry 2"
        ),
        "finite_volume_noncompact": True,
        "literal_47_wall_subgroup_chamber": (
            "not located in the bounded ordinary/doubling CoxIter audit; "
            "worldwide priority open"
        ),
        "reflection_subgroup_index": 1_566_720,
    }
