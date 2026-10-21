import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d30 import analyze_d30


DISC11_PATH = Path(
    os.environ.get(
        "SHIMADA_DISC11",
        "/tmp/shimada_threediscs/threediscs/disc11.txt",
    )
)


@pytest.mark.skipif(
    not DISC11_PATH.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_physical_D30_two_path_identity_and_restricted_facets() -> None:
    result = analyze_d30(DISC11_PATH)

    assert result["verified"] is True
    state = result["physical_state"]
    assert state["representative"] == 30
    assert state["D68_D181_D30_equals_D68_D86_D30"] is True
    assert state["chamber_stabilizer_order"] == 1
    assert state["fingerprints_verified"] is True

    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 34
    assert section["positive_proportional_classes"] == 26
    assert section["irredundant_facets"] == 12
    assert section["root_facets"] == 8
    assert section["internal_facets"] == 4
    assert section["internal_facet_classes"] == [
        [[19, 0]], [[20, 0]], [[22, 0]], [[26, 0]]
    ]
    assert len(section["facet_witnesses"]) == 12
    assert len(section["redundancy_dependences"]) == 14
    assert all(item["verified"] for item in section["redundancy_dependences"])
    assert section["generic_exact_halfspace_kernel_verified"] is True


@pytest.mark.skipif(
    not DISC11_PATH.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_physical_D30_root_filter_and_internal_first_star() -> None:
    result = analyze_d30(DISC11_PATH)

    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 18
    assert roots["integral_field_roots"] == 17
    assert roots["shell_distribution"] == {"2": 8, "4 - 2*phi": 9}
    assert roots["field_span_rank"] == 9
    assert roots["non_field_root_walls"] == [[10, 0]]
    assert roots["field_root_nonfacet_walls"] == [[13, 0]]
    assert roots["root_walls_on_section_facets"] == 16

    star = result["internal_first_star"]
    assert star["all_crossings_verified"] is True
    assert star["representative_targets"] == [16, 86, 61, 181]
    by_orbit = {item["D30_orbit"]: item for item in star["crossings"]}
    assert by_orbit[20]["returns_to_same_physical_D86"] is True
    assert by_orbit[26]["returns_to_same_physical_D181"] is True
