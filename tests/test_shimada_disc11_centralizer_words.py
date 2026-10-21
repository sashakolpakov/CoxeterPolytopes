from pathlib import Path

import pytest

from scripts.search_shimada_disc11_centralizer_words import search_length_two


DISC11 = Path("/tmp/shimada_threediscs/threediscs/disc11.txt")


@pytest.mark.skipif(not DISC11.exists(), reason="external Shimada disc11 atlas absent")
def test_only_the_known_involution_centralizes_through_word_length_three() -> None:
    result = search_length_two(DISC11, max_word_length=3)

    assert result["verified"] is True
    assert result["individual_centralizing_generators"] == 0
    assert result["centralizing_ordered_length_two_words"] == 682
    assert result["nonidentity_centralizing_length_two_actions"] == 0
    length_three = result["exhaustive_length_three"]
    assert length_three["centralizing_ordered_words"] == 4989
    assert length_three["distinct_centralizing_actions"] == 2
    assert {
        action["characteristic_polynomial"]
        for action in length_three["distinct_actions"]
    } == {
        "(lambda - 1)**20",
        "(lambda - 1)**18*(lambda + 1)**2",
    }
