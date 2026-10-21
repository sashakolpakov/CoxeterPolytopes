from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d181 import analyze_d181


DISC11 = Path("/tmp/shimada_threediscs/threediscs/disc11.txt")


@pytest.mark.skipif(not DISC11.exists(), reason="external Shimada disc11 atlas absent")
def test_physical_D181_restricted_section_and_internal_star() -> None:
    result = analyze_d181(DISC11)

    assert result["verified"] is True
    assert result["physical_state"] == {
        "reached_from": 68,
        "crossed_D68_wall": [14, 0],
        "representative": 181,
        "primitive_wall_count": 26,
        "wall_orbit_count": 16,
        "root_wall_orbits": list(range(10)),
    }
    assert result["positive_proportional_restrictions"]["class_count"] == 17

    section = result["irredundant_section"]
    assert section["facet_count"] == 11
    assert section["root_facets"] == 9
    assert section["internal_facets"] == 2
    assert section["internal_facet_classes"] == [[[10, 0]], [[11, 0]]]
    assert len(section["positive_timelike_face_witnesses"]) == 11
    assert len(section["positive_redundancy_dependences"]) == 6
    assert all(
        item["verified"] for item in section["positive_redundancy_dependences"]
    )
    assert section["generic_exact_halfspace_kernel_verified"] is True

    crossings = result["one_star_farther"]
    assert [item["to_physical_representative"] for item in crossings] == [30, 68]
    assert all(item["verified"] for item in crossings)
    assert crossings[0]["returns_to_original_physical_D68"] is False
    assert crossings[1]["returns_to_original_physical_D68"] is True
