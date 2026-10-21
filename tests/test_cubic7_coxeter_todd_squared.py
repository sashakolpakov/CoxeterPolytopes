import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    SUPPORTED_K12_SEEDS,
    TWO_THETA_K12_SEED,
    screen_cubic7_k12_seed,
    two_theta_k12_first_two_packets,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h100.json"


def test_finite_support_manifest_is_exact_and_distinct():
    assert len(SUPPORTED_K12_SEEDS) == 19
    assert len({seed.seed_id for seed in SUPPORTED_K12_SEEDS}) == 19
    assert [seed.radical_norm for seed in SUPPORTED_K12_SEEDS] == [
        1,
        7,
        8,
        13,
        13,
        13,
        27,
        56,
        91,
        91,
        91,
        104,
        104,
        104,
        189,
        216,
        351,
        351,
        351,
    ]
    assert all(seed.coefficient.sign(0) > 0 for seed in SUPPORTED_K12_SEEDS)
    assert all(
        seed.coefficient.sign(place) < 0
        for seed in SUPPORTED_K12_SEEDS
        for place in (1, 2)
    )


def test_h100_necessary_screen_has_exactly_three_surviving_lanes():
    screens = tuple(screen_cubic7_k12_seed(seed, 100) for seed in SUPPORTED_K12_SEEDS)
    nonempty = tuple(screen for screen in screens if screen.nonempty)
    assert [(screen.seed.support, len(screen.nonempty)) for screen in nonempty] == [
        (("p2",), 1),
        (("p2", "p7"), 1),
        (("p2", "p3"), 1),
    ]
    assert sum(len(screen.nonempty) for screen in screens) == 3
    assert sum(
        len(types) for screen in screens for _level, types in screen.nonempty
    ) == 3


def test_selected_lane_first_two_packets_are_realized_counted_and_acute():
    assert TWO_THETA_K12_SEED.coefficient == Cubic7(0, 2)
    packets = two_theta_k12_first_two_packets()
    assert packets.first_level.height == Cubic7(5, 14, 6)
    assert packets.first_wall_count == 1_512
    assert packets.first_same_height_maximum == Cubic7(-4, -14, -6)
    assert packets.second_level.height == Cubic7(16, 45, 20)
    assert packets.minimal_vectors_per_factor == 756
    assert packets.minimal_neighbors_inner_product_minus_one == 160
    assert packets.second_wall_count == 241_920
    assert packets.first_second_maximum_bound == Cubic7(-3, -7, -3)
    assert packets.second_same_height_maximum_bound == Cubic7(-1)


def test_frozen_artifact_records_bounded_negative_conclusion():
    payload = json.loads(ARTIFACT.read_text())
    assert payload["schema"] == "coxeter.cubic7_k12_squared_h100.v1"
    assert payload["local_arithmetic"]["seed_count"] == 19
    assert payload["screen_summary"]["empty_necessary_seed_count"] == 16
    assert payload["screen_summary"]["nonempty_necessary_seed_count"] == 3
    selected = payload["selected_lane"]
    assert selected["first_packet"]["wall_count"] == 1_512
    assert selected["second_packet"]["wall_count"] == 241_920
    assert selected["minimal_shell_audit"]["neighbors_at_inner_product_minus_one"] == 160
    assert payload["conclusion"]["new_compact_example"] is False
    assert payload["conclusion"]["new_finite_volume_example"] is False
    assert "later walls remain open" in payload["conclusion"]["scope"]
