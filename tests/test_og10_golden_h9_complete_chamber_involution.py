import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from scripts.certify_og10_golden_h9_complete_chamber_involution import build


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/candidates/og10_golden_h9_complete_chamber_involution.json"
)


def test_frozen_h9_complete_chamber_involution_is_reproducible() -> None:
    result = build()
    assert result == json.loads(ARTIFACT.read_text())
    assert result["verified"] is True
    assert all(result["checks"].values())


def test_pell_and_inverse_fold_to_the_same_complete_chamber_reflection() -> None:
    result = build()
    fold = result["pell_fold"]
    assert fold["pell_violated_wall_index_in_61"] == 9
    assert fold["inverse_violated_wall_index_in_61"] == 8
    assert fold["reflection_normal_square"] == Golden(4, 4).to_json()

    membership = result["complete_chamber_membership"]
    assert membership["archived_fixed_wall_count"] == 61
    assert membership["largest_archived_wall_pairing"] == {
        "a": "-1/100",
        "b": "0",
    }
    assert membership["possible_violator_x_square_bound"] == {
        "a": "202727309617/200773372141",
        "b": "-1021288932/200773372141",
    }
    assert membership["next_uncomputed_control_coordinate_square"] == {
        "a": "45",
        "b": "72",
    }


def test_precise_ambient_extension_scope() -> None:
    result = build()
    obstruction = result["golden_H10_extension_obstruction"]
    assert obstruction["both_lifts_are_O_nonintegral"] == [True, True]
    assert obstruction["both_lifts_have_y_to_u_coefficient"] == "-1/2"

    extension = result["full_OG10_lattice_extension"]
    assert extension["verified"] is True
    assert all(extension["checks"].values())
    assert extension["positive_complement_characteristic_coefficients"] == [
        1, -2, 0, 2, -1
    ]
    assert extension["positive_complement_one_eigenspace_dimension"] == 1
    assert "does not preserve the original period plane" in extension["scope"]
