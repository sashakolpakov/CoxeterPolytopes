import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d37 import analyze_d37


DISC11_PATH = Path(
    os.environ.get("SHIMADA_DISC11", "/tmp/shimada_threediscs/threediscs/disc11.txt")
)


@pytest.mark.skipif(not DISC11_PATH.is_file(), reason="Shimada disc11.txt is not vendored")
def test_exact_root_facet_branch_selection() -> None:
    result = analyze_d37(DISC11_PATH)
    selection = result["branch_selection"]
    assert selection["verified"] is True
    assert selection["D84_root_screen"]["verified"] is True
    assert selection["D84_root_screen"]["root_facets"] == 6
    assert selection["D37_root_facets"] == 9
    assert selection["selected"] == 37


@pytest.mark.skipif(not DISC11_PATH.is_file(), reason="Shimada disc11.txt is not vendored")
def test_physical_D37_exact_section_roots_and_crossings() -> None:
    result = analyze_d37(DISC11_PATH)
    assert result["verified"] is True
    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 31
    assert section["positive_proportional_classes"] == 22
    assert section["irredundant_facets"] == 11
    assert section["root_facets"] == 9
    assert section["internal_facets"] == 2
    assert section["generic_exact_halfspace_kernel_verified"] is True
    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 20
    assert roots["integral_field_roots"] == 19
    assert roots["field_span_rank"] == 10
    assert roots["non_field_root_walls"] == [[11, 0]]
    assert roots["field_root_nonfacet_walls"] == [[10, 0]]
    star = result["internal_first_star"]
    assert star["representative_targets"] == [61, 22]
    assert star["returns_to_same_physical_D61"] is True
    assert star["return_wall_is_opposite_D61_entry_wall"] is True
    assert star["new_physical_neighbor"] == 22
