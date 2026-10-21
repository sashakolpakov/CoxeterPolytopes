import os
from pathlib import Path

import pytest

from scripts.analyze_shimada_disc11_d55 import analyze_d55


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
def test_physical_D55_exact_restricted_section() -> None:
    result = analyze_d55(DISC11_PATH)

    assert result["verified"] is True
    assert result["physical_state"]["representative"] == 55
    assert result["physical_state"]["fingerprints_verified"] is True
    section = result["restricted_section"]
    assert section["ambient_primitive_walls"] == 41
    assert section["positive_proportional_classes"] == 33
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
def test_physical_D55_roots_and_physical_returns() -> None:
    result = analyze_d55(DISC11_PATH)

    roots = result["root_filter"]
    assert roots["ambient_minus_two_walls"] == 18
    assert roots["integral_field_roots"] == 17
    assert roots["shell_distribution"] == {"2": 8, "4 - 2*phi": 9}
    assert roots["field_span_rank"] == 10
    assert roots["non_field_root_walls"] == [[3, 1]]
    assert roots["all_root_facet_members_are_field_roots"] is True

    star = result["internal_first_star"]
    assert star["all_crossings_verified"] is True
    assert star["representative_targets"] == [20, 14, 16, 86]
    assert star["crossing_physical_wall_sets_pairwise_distinct"] is True
    assert star["wall_16_member_1_is_opposite_D86_entry_wall"] is True
    assert star["wall_16_member_1_returns_to_same_physical_D86"] is True
    assert star["crossings"][-1]["physical_matches_among_known_states"] == ["D86"]
