import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d27 import analyze_d27


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
def test_physical_D27_exact_restricted_section() -> None:
    result = analyze_d27(DISC11_PATH)

    assert result["verified"] is True
    assert result["physical_state"]["representative"] == 27
    assert result["physical_state"]["fingerprints_verified"] is True
    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 33
    assert section["positive_proportional_classes"] == 24
    assert section["irredundant_facets"] == 12
    assert section["root_facets"] == 8
    assert section["internal_facets"] == 4
    assert section["partition_verified"] is True
    assert section["witnesses_verified"] is True
    assert section["dependences_verified"] is True
    assert section["generic_exact_halfspace_kernel_verified"] is True


@pytest.mark.skipif(
    not DISC11_PATH.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_physical_D27_roots_simultaneous_crossing_and_diamond() -> None:
    result = analyze_d27(DISC11_PATH)

    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 19
    assert roots["integral_field_roots"] == 18
    assert roots["shell_distribution"] == {"2": 9, "4 - 2*phi": 9}
    assert roots["field_span_rank"] == 9
    assert roots["non_field_root_walls"] == [[14, 0]]
    assert roots["all_root_facet_members_are_field_roots"] is True

    star = result["internal_first_star"]
    assert star["verified"] is True
    assert star["simple_representative_targets"] == [86, 14, 61]
    assert star["wall_20_returns_to_same_physical_D86"] is True
    assert star["D27_and_D55_reach_same_physical_D14"] is True
    assert star["simultaneous_paths_same_physical_endpoint"] is True
    assert star["simultaneous_endpoint_is_distinct_physical_D27"] is True
    assert [item["representative_sequence"] for item in star["simultaneous_paths"]] == [
        [27, 11, 11, 74, 27],
        [27, 74, 11, 11, 27],
    ]
