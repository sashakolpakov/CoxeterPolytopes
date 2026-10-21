from __future__ import annotations

import json

from coxeter_polytopes.exact.golden import Golden
from scripts.certify_bugaenko_h8_weyl_vector_obstruction import OUTPUT, analyze


def _archive() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_replay_matches_frozen_certificate() -> None:
    payload = analyze()
    assert payload == _archive()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["certificate_sha256"] == (
        "b58a500597702a42fb739f03109bd30016ea6ea625df8803335720ab3dde2b61"
    )


def test_first_nine_walls_force_a_candidate_rejected_by_the_last_two() -> None:
    candidate = _archive()["unique_nine_wall_candidate"]
    assert candidate["pivot_root_coordinate_determinant"] == {"a": "-3", "b": "2"}
    assert candidate["failing_wall_indices_zero_based"] == [9, 10]
    assert candidate["pairings_with_all_eleven_walls"][:9] == [
        {"a": "-1", "b": "0"}
    ] * 9
    assert candidate["pairings_with_all_eleven_walls"][9:] == [
        {"a": "3", "b": "-8"}
    ] * 2
    residual = Golden.from_json(candidate["common_nonzero_residual"])
    assert residual == Golden(4, -8)
    assert residual.norm() == -80


def test_inconsistency_is_over_the_field_not_an_integrality_failure() -> None:
    payload = _archive()
    assert payload["system"]["coefficient_rank"] == 9
    assert payload["system"]["augmented_rank"] == 10
    assert payload["conclusion"] == {
        "generalized_weyl_vector_exists_over_Q_phi": False,
        "integral_weyl_vector_exists_over_Z_phi": False,
        "reason": (
            "the first nine wall equations have a unique solution, while each "
            "remaining wall has residual 4-8*phi of field norm -80"
        ),
        "scope": (
            "direct linear-algebra obstruction for the marked Bugaenko H8 "
            "chamber; not an application or extension of the integral "
            "Z-lattice classification of Gritsenko--Nikulin"
        ),
    }
