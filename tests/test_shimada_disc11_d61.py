import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d61 import analyze_d61


DISC11_PATH = Path(
    os.environ.get("SHIMADA_DISC11", "/tmp/shimada_threediscs/threediscs/disc11.txt")
)


@pytest.mark.skipif(not DISC11_PATH.is_file(), reason="Shimada disc11.txt is not vendored")
def test_common_physical_D61_and_exact_section() -> None:
    result = analyze_d61(DISC11_PATH)
    assert result["verified"] is True
    state = result["physical_state"]
    assert state["D30_wall_22_and_D27_wall_23_same_row_isometry"] is True
    assert state["D30_wall_22_and_D27_wall_23_same_physical_wall_set"] is True
    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 32
    assert section["positive_proportional_classes"] == 25
    assert section["irredundant_facets"] == 11
    assert section["root_facets"] == 7
    assert section["internal_facets"] == 4
    assert section["generic_exact_halfspace_kernel_verified"] is True


@pytest.mark.skipif(not DISC11_PATH.is_file(), reason="Shimada disc11.txt is not vendored")
def test_D61_root_filter_and_physical_first_star() -> None:
    result = analyze_d61(DISC11_PATH)
    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 18
    assert roots["integral_field_roots"] == 17
    assert roots["shell_distribution"] == {"2": 8, "4 - 2*phi": 9}
    assert roots["field_span_rank"] == 10
    assert roots["non_field_root_walls"] == [[16, 0]]
    star = result["internal_first_star"]
    assert star["representative_targets"] == [84, 30, 27, 37]
    assert star["returns_to_same_physical_D30"] is True
    assert star["returns_to_same_physical_D27"] is True
    assert star["return_walls_are_opposite_entry_walls"] is True
    assert star["new_physical_neighbors"] == [84, 37]
