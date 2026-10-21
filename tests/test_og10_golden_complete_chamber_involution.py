import json
from pathlib import Path

from scripts.probe_og10_golden_chamber_folding import build


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/candidates/og10_golden_complete_chamber_involution.json"
)


def test_complete_chamber_involution_and_monodromy_extension() -> None:
    frozen = json.loads(ARTIFACT.read_text())
    result = build()

    assert result["verified"] is True
    assert frozen["verified"] is True
    assert result["conclusion"] == frozen["conclusion"]

    checks = result["checks"]
    assert checks["h_is_nonidentity"] is True
    assert checks["h_has_order_two"] is True
    assert checks["h_is_integral_reflection_in_norm_6_plus_2phi_vector"] is True
    assert checks["h_endpoint_global_violator_bound_below_100"] is True
    assert checks["h_retains_274_and_exports_183_of_457_roots"] is True
    assert checks["h_closure_has_640_roots"] is True

    involution = result["complete_chamber_involution"]
    assert involution["field_matrix_sha256"] == frozen[
        "complete_chamber_involution"
    ]["field_matrix_sha256"]
    assert involution["integral_reflection_normal_square"] == {
        "a": "6",
        "b": "2",
    }
    assert involution["global_membership_x_square_bound"] == frozen[
        "global_membership"
    ]["folded_h_endpoint_possible_violator_x_square_bound"]

    walls = result["chamber_involution_wall_transport"]
    frozen_walls = frozen["wall_transport"]
    assert walls["new_roots_sha256"] == frozen_walls["new_roots_sha256"]
    assert walls["closed_root_set_sha256"] == frozen_walls[
        "closed_root_set_sha256"
    ]
    assert walls["involution_closed_root_count"] == 640
    assert walls["images_outside_archived_prefix"] == 183

    monodromy = result["monodromy_and_birational_scope"]
    frozen_monodromy = frozen["monodromy_extension"]
    assert monodromy["verified"] is True
    assert monodromy["gluing_generator_action"] == [
        "n1_to_n1",
        "n2_to_minus_n2",
    ]
    assert monodromy["trace_action_sha256"] == frozen_monodromy[
        "trace_action_sha256"
    ]
    assert monodromy["ambient_extension_sha256"] == frozen_monodromy[
        "ambient_extension_sha256"
    ]

    # This theorem is finite: it must never be silently relabelled as the
    # infinite-index obstruction sought by the dual search.
    assert result["complete_chamber_involution"]["order"] == 2
    assert "no_infinite_order" in result["conclusion"]
