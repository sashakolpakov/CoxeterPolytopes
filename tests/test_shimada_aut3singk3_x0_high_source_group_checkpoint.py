from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x0_high_source_group_checkpoint.json"
)
CHECKPOINT_SHA256 = "6773b94c9650d27deb5fcec70275cc9cef001fb35ae3b9942c0c37e6f1972931"


def _payload() -> dict:
    return json.loads(CHECKPOINT.read_text(encoding="utf-8"))


def test_x0_checkpoint_is_frozen() -> None:
    assert hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest() == CHECKPOINT_SHA256
    payload = _payload()
    assert payload["schema"] == (
        "coxeter.shimada.aut3singk3.x0.high_source_group_checkpoint.v1"
    )
    assert payload["all_checks_pass"] is True


def test_x0_complete_source_and_high_source_group_census() -> None:
    payload = _payload()
    stock = payload["complete_source_stock"]
    assert stock["neron_severi_rank"] == 20
    assert stock["neron_severi_determinant"] == -36
    assert stock["wall_count"] == 3278
    assert stock["wall_orbit_sizes"] == [
        60, 40, 180, 10, 144, 144, 240, 360, 180, 240, 240, 720, 720
    ]
    assert stock["genuine_ample_stabilizer"] == "Aut(X_0,a_0), order 1440"

    census = payload["high_source_group_census"]
    assert census["fixed_rank_histogram"] == {
        "2": 144,
        "4": 584,
        "6": 450,
        "8": 170,
        "10": 46,
        "12": 45,
        "20": 1,
    }
    assert census["unique_high_class"] == {
        "cardinality": 45,
        "centralizer_order": 32,
        "fixed_rank": 12,
        "hyperbolic_dimension": 11,
        "representative_zero_based": 3,
    }
    assert census["high_closed_V4_count"] == 30
    assert census["high_closed_V4_common_fixed_rank_histogram"] == {"8": 30}
    assert census["only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup"]


def test_x0_downstream_gates_remain_explicitly_open() -> None:
    pending = _payload()["pending_gates"]
    assert pending == {
        "Coxeter_pair_angle_gate": "pending",
        "Vinberg_volume_gate": "not authorized before a Coxeter survivor",
        "exact_irredundancy": "pending",
        "full_abstract_oriented_weighted_wall_group": "pending",
        "restriction_of_all_3278_walls": "pending",
    }

