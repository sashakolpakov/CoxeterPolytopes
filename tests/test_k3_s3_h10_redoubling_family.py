from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.certify_k3_s3_h10_redoubling_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/k3_s3_h10_redoubling_family.json"


@pytest.fixture(scope="module")
def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def fresh() -> dict:
    return build_certificate(include_lineage=True)


def test_archived_payload_hash_and_global_checks(archived: dict) -> None:
    payload = dict(archived)
    expected = payload.pop("payload_sha256")
    actual = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert actual == expected
    assert archived["all_checks_pass"]
    assert all(archived["checks"].values())


def test_parallel_doubling_normalization_is_not_confused(archived: dict) -> None:
    warning = archived["seed"]["parallel_normalization_warning"]
    assert archived["seed"]["doubling_walls_one_based"] == [10, 13]
    assert warning == {
        "conclusion": "parallel, not ultraparallel",
        "norm_product": 4,
        "normalized_squared_ratio": "1",
        "raw_inner_product": 2,
        "raw_inner_product_squared": 4,
    }


def test_arbitrary_path_theorem_and_parity_fingerprints(archived: dict) -> None:
    rows = archived["formula_replay_N1_through_N12"]
    assert [row["N"] for row in rows] == list(range(1, 13))
    for row in rows:
        expected_same = row["N"] % 2 == 1
        assert row["A_B_relation_isomorphic"] is expected_same
        assert row["A_B_lattice_decorated_Gram_isomorphic"] is expected_same
        for orientation in row["orientations"]:
            assert orientation["facets"] == row["N"] + 12
            assert orientation["doubling_wall_count"] == 2


def test_first_two_positive_doubles_have_exact_finite_volume(archived: dict) -> None:
    audits = archived["exact_N2_volume_audits"]
    assert (
        audits["N2_A"]["finite_vertices"],
        audits["N2_A"]["cusps"],
        audits["N2_A"]["elliptic_rank9_edges"],
    ) == (41, 5, 252)
    assert (
        audits["N2_B"]["finite_vertices"],
        audits["N2_B"]["cusps"],
        audits["N2_B"]["elliptic_rank9_edges"],
    ) == (52, 6, 318)
    assert all(all(audit["checks"].values()) for audit in audits.values())


def test_exact_volume_and_index_scale_for_every_N(archived: dict) -> None:
    invariants = archived["theorem"]["exact_orbifold_invariants"]
    assert invariants["seed_rational_Euler_characteristic"] == "-527/33443020800"
    assert invariants["P_N_rational_Euler_characteristic"] == "-N*527/33443020800"
    assert invariants["P_N_covolume"] == "N*pi^5*527/987614208000"
    for row in archived["representative_polytopes"]:
        assert row["reflection_subgroup_index"] == row["copy_count_N"]
        assert row["volume_multiple_of_seed"] == row["copy_count_N"]


def test_bounded_lineage_uses_full_labelled_graphs(archived: dict) -> None:
    lineage = archived["bounded_lineage"]
    assert all(lineage["checks"].values())
    assert lineage["corpus_inventory"]["combined_H10_records"] == 49
    nonseed = [
        row
        for row in lineage["histogram_collision_controls"]
        if row["candidate"] != "N1_A"
    ]
    assert len(nonseed) == 1
    assert nonseed[0]["candidate"] == "N2_B"
    assert nonseed[0]["same_global_histogram"]
    assert not nonseed[0]["same_node_profile_multiset"]
    assert not nonseed[0]["exact_labelled_graph_isomorphic"]


def test_fresh_exact_replay_matches_archive(fresh: dict, archived: dict) -> None:
    assert fresh == archived
