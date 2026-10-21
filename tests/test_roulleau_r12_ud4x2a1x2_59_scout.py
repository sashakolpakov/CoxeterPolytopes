import json
from pathlib import Path

from scripts.scout_roulleau_r12_ud4x2a1x2_59 import coordinate_model
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_scout.json"


def test_formula_coordinate_reconstruction() -> None:
    basis_gram, walls, wall_gram = coordinate_model()
    assert basis_gram.det() == -64
    assert exact_inertia(basis_gram) == (1, 11, 0)
    assert len(walls) == 59
    assert wall_gram.rank() == 12
    assert exact_inertia(wall_gram) == (1, 11, 47)
    assert all(wall_gram[index, index] == -2 for index in range(59))


def test_frozen_weighted_group_denominator() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert all(payload["coordinate_checks"].values())
    assert payload["weighted_group"]["order"] == 2880
    assert payload["weighted_group"]["element_order_histogram"] == {
        "1": 1,
        "2": 303,
        "3": 80,
        "4": 720,
        "5": 144,
        "6": 1200,
        "10": 432,
    }
    assert all(payload["weighted_group"]["checks"].values())
    assert len(payload["weighted_group"]["elements_zero_based"]) == 2880
