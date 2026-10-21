import importlib.util
from pathlib import Path

import pytest

from scripts.search_shimada_disc11_centralizer_words4 import search_length_four


DISC11 = Path("/tmp/shimada_threediscs/threediscs/disc11.txt")
HAS_FLINT = importlib.util.find_spec("flint") is not None


@pytest.mark.skipif(
    not DISC11.exists() or not HAS_FLINT,
    reason="external Shimada disc11 atlas or optional python-flint absent",
)
def test_only_identity_and_known_involution_centralize_at_length_four() -> None:
    result = search_length_four(DISC11)
    exact = result["exact_length_four"]

    assert result["verified"] is True
    assert exact["ordered_words_exhausted"] == 346_083_947_521
    assert exact["half_words_per_side"] == 588_289
    assert exact["modular_candidate_quadruples"] == 1_139_377
    assert exact["exact_integral_commutator_checks"] == 1_139_377
    assert exact["false_modular_candidates"] == 0
    assert exact["centralizing_ordered_words"] == 1_139_377
    assert exact["identity_word_multiplicity"] == 1_138_938
    assert exact["nonidentity_word_multiplicity"] == 439
    assert exact["distinct_centralizing_actions"] == 2
    assert exact["word_quadruples_sha256"] == (
        "012f1d1e0314b1a70925cafa7e5e633971d86d55ff27d51393fe92b675df352c"
    )

    actions = {
        action["column_action_sha256"]: action["ordered_word_multiplicity"]
        for action in exact["distinct_actions"]
    }
    assert actions == {
        "8817db0b345cc277b1b2e649062c98d9a923e6a5159c08ff8b075a3947d0d259": 1_138_938,
        "6b4bc88ee93385bcb4a4c55828924860aa90158f84a8469ffc238dbe4f33e4eb": 439,
    }
    nonidentity = [
        action for action in exact["distinct_actions"] if not action["is_identity"]
    ]
    assert len(nonidentity) == 1
    assert nonidentity[0]["is_known_D68_D27_involution"] is True
