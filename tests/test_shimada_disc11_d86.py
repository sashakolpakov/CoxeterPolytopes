import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d86 import analyze_d86


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
def test_physical_D86_restricted_facet_certificate() -> None:
    result = analyze_d86(DISC11_PATH)

    assert result["verified"] is True
    assert result["physical_state"]["representative"] == 86
    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 35
    assert section["positive_proportional_classes"] == 27
    assert section["irredundant_facets"] == 11
    assert section["root_facets"] == 7
    assert section["internal_facets"] == 4
    assert section["partition_verified"] is True
    assert section["witnesses_verified"] is True
    assert section["dependences_verified"] is True
    assert section["generic_exact_halfspace_kernel_verified"] is True


@pytest.mark.skipif(
    not DISC11_PATH.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_physical_D86_root_filter_and_internal_first_star() -> None:
    result = analyze_d86(DISC11_PATH)

    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 17
    assert roots["integral_field_roots"] == 16
    assert roots["shell_distribution"] == {"2": 8, "4 - 2*phi": 8}
    assert roots["field_span_rank"] == 9
    assert roots["non_field_root_walls"] == [[15, 0]]

    star = result["internal_first_star"]
    assert star["all_crossings_verified"] is True
    assert star["representative_targets"] == [55, 30, 27, 68]
    assert star["orbit_26_is_opposite_D68_entry_wall"] is True
