import json
import os
from pathlib import Path

import pytest

from scripts.explore_shimada_disc11_ambient_carrier import analyze_ambient_carrier


DISC11 = Path(
    os.environ.get(
        "SHIMADA_DISC11",
        "/tmp/shimada_threediscs/threediscs/disc11.txt",
    )
)
RADIUS5 = (
    Path(__file__).parents[1]
    / "data"
    / "candidates"
    / "shimada_disc11_ambient_carrier_radius5.json"
)
RADIUS7 = (
    Path(__file__).parents[1]
    / "data"
    / "candidates"
    / "shimada_disc11_ambient_carrier_radius7.json"
)
RADIUS8 = (
    Path(__file__).parents[1]
    / "data"
    / "candidates"
    / "shimada_disc11_ambient_carrier_radius8.json"
)


def test_archived_exact_ambient_carrier_radius5_certificate() -> None:
    result = json.loads(RADIUS5.read_text())

    assert result["verified"] is True
    assert result["complete_ambient_carrier_radius"] == 5
    assert result["complete_layer_counts"] == {
        "0": 1,
        "1": 4,
        "2": 7,
        "3": 12,
        "4": 17,
        "5": 23,
    }
    assert result["states_in_complete_carrier_ball"] == 64
    assert result["queued_radius6_frontier_states"] == 28
    assert result["auxiliary_ambient_wall_crossings_tested"] == 1408
    centralizer = result["centralizer_decision"]
    assert centralizer["only_identity_and_A_in_complete_carrier_ball"] is True
    assert centralizer[
        "any_centralizer_element_other_than_identity_or_A_requires_ambient_carrier_distance_at_least"
    ] == 6
    assert centralizer["known_V4_elements_not_seen"] == ["B_minus", "B_plus"]


def test_archived_exact_ambient_carrier_radius7_certificate() -> None:
    result = json.loads(RADIUS7.read_text())

    assert result["verified"] is True
    assert result["complete_ambient_carrier_radius"] == 7
    assert result["complete_layer_counts"] == {
        "0": 1,
        "1": 4,
        "2": 7,
        "3": 12,
        "4": 17,
        "5": 23,
        "6": 28,
        "7": 36,
    }
    assert result["states_in_complete_carrier_ball"] == 128
    assert result["queued_radius8_frontier_states"] == 47
    assert result["auxiliary_ambient_wall_crossings_tested"] == 3172
    assert result["state_records_sha256"] == (
        "47e4aaca8745af730cba3ad64ac649875439c286953f5d6eb2b3b4315194fa09"
    )
    centralizer = result["centralizer_decision"]
    assert centralizer["only_identity_and_A_in_complete_carrier_ball"] is True
    assert centralizer[
        "any_centralizer_element_other_than_identity_or_A_requires_ambient_carrier_distance_at_least"
    ] == 8
    assert centralizer["known_V4_elements_not_seen"] == ["B_minus", "B_plus"]
    assert result["D35_regression"][
        "both_branches_occur_in_the_radius7_BFS"
    ] is True


def test_archived_exact_ambient_carrier_radius8_certificate() -> None:
    result = json.loads(RADIUS8.read_text())

    assert result["verified"] is True
    assert result["complete_ambient_carrier_radius"] == 8
    assert result["complete_layer_counts"] == {
        "0": 1,
        "1": 4,
        "2": 7,
        "3": 12,
        "4": 17,
        "5": 23,
        "6": 28,
        "7": 36,
        "8": 47,
    }
    assert result["states_in_complete_carrier_ball"] == 175
    assert result["queued_radius9_frontier_states"] == 65
    assert result["auxiliary_ambient_wall_crossings_tested"] == 4572
    assert result["retained_crossings"] == 874
    assert result["rejected_crossings"] == 3698
    assert result["root_boundary_wall_members_not_crossed"] == 3486
    assert result["state_records_sha256"] == (
        "3b4e5379a3bac0ce8535fc24aebff3c7b73bb3a4df13f204b92ee5dcbc58062f"
    )
    assert result["queue_records_sha256"] == (
        "ce10f96ca06b190ba4d8dd86103e8bbf1976ee1a415e3f7e0ebb4cb83f3b5adf"
    )
    assert result["edge_records_sha256"] == (
        "18ac41c6802a5761a80cb198f7877ca6aa51e6addcfd13c6df3cfd15f1124b08"
    )
    assert result["checkpoint_payload_sha256"] == (
        "08089ebf8498f87b4674a454f9e072aa620ca00ee571f73df21b63bd3a6347ad"
    )
    assert result["checkpoint_file_sha256"] == (
        "04ad5e761c6ee7e579264aeb7909e807fa471fba5cdd0dd5347ed4dfedf3a58f"
    )
    assert result["checkpoint_file_size_bytes"] == 123929116
    centralizer = result["centralizer_decision"]
    assert centralizer["only_identity_and_A_in_complete_carrier_ball"] is True
    assert centralizer[
        "any_centralizer_element_other_than_identity_or_A_requires_ambient_carrier_distance_at_least"
    ] == 9
    assert centralizer[
        "missing_known_V4_elements_require_ambient_carrier_distance_at_least"
    ] == 9
    assert centralizer["known_V4_elements_not_seen"] == ["B_minus", "B_plus"]
    assert centralizer["D68_return_records_sha256"] == (
        "f7017aca8c10636741179a308121bd8067391b7c1fb33f2e72caf5be1631fe30"
    )
    assert result["D35_regression"][
        "both_branches_occur_in_the_radius8_BFS"
    ] is True


@pytest.mark.skipif(
    not DISC11.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_exact_ambient_carrier_base_layer(tmp_path: Path) -> None:
    checkpoint = tmp_path / "ambient-carrier.json"
    result = analyze_ambient_carrier(
        DISC11,
        max_depth=0,
        checkpoint_path=checkpoint,
        resume=False,
        checkpoint_every=1,
    )

    assert result["verified_exact_closed_section_decisions"] is True
    assert result["processed_carrier_states"] == 1
    assert result["queued_carrier_states"] == 4
    assert result["auxiliary_ambient_wall_crossings_tested"] == 11
    assert result["retained_crossings"] == 4
    assert result["rejected_crossings"] == 7
    assert result["root_boundary_wall_members_not_crossed"] == 17
    assert {item["representative"] for item in result["queue"]} == {
        10,
        69,
        86,
        181,
    }
    assert result["edge_records_sha256"] == (
        "86a941f8616f9771a8d4f583177dc21e071634d73d1941382fe5c9cb7dff7b92"
    )

    centralizer = result["centralizer_decision_ball"]
    assert centralizer["complete_ambient_carrier_radius"] == 0
    assert centralizer["complete_layer_counts"] == {"0": 1}
    assert centralizer["states_in_complete_carrier_ball"] == 1
    assert centralizer["D68_return_candidates_tested"] == 2
    assert centralizer["centralizing_automorphisms_in_complete_carrier_ball"] == 1
    assert centralizer["known_V4_elements_seen"] == ["identity"]
    assert centralizer["only_identity_and_A_in_complete_carrier_ball"] is False
    assert centralizer[
        "centralizing_set_in_complete_carrier_ball_is_exactly_known_V4"
    ] is False
    assert centralizer["carrier_queue_closed"] is False
    assert centralizer["decision_status"] == (
        "pending_larger_exact_ambient_carrier_ball"
    )

    resumed = analyze_ambient_carrier(
        DISC11,
        max_depth=0,
        checkpoint_path=checkpoint,
        checkpoint_every=1,
    )
    assert resumed["checkpoint"]["resumed"] is True
    assert resumed["checkpoint"]["initial_processed_states"] == 1
    assert resumed["edge_records_sha256"] == result["edge_records_sha256"]
    assert resumed["queue_records_sha256"] == result["queue_records_sha256"]
