import os
from pathlib import Path

import pytest

from scripts.certify_shimada_disc11_d35_facet_subdivision import (
    certify_d35_facet_subdivision,
)
from scripts.explore_shimada_disc11_ambient_carrier import (
    certify_d35_carrier_regression,
)


DISC11 = Path(
    os.environ.get(
        "SHIMADA_DISC11",
        "/tmp/shimada_threediscs/threediscs/disc11.txt",
    )
)


@pytest.mark.skipif(
    not DISC11.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_d35_restricted_facet_has_two_distinct_opposite_induced_cells() -> None:
    certificate = certify_d35_facet_subdivision(DISC11)

    assert certificate["verified"] is True
    assert certificate["state_representative"] == 35
    assert certificate["physical_state_sha256"] == (
        "69782e666dacd063973e126a43a9cdbacc0d2124e152efc053c3403688653c6d"
    )
    assert certificate["wall_class"] == [[26, 0]]
    assert certificate["distinct_opposite_physical_states"] == 2

    crossings = certificate["witness_crossings"]
    assert [record["gallery_representatives"] for record in crossings] == [
        [35, 162, 137, 122],
        [35, 162, 137],
    ]
    assert {
        record["target_physical_state_sha256"] for record in crossings
    } == {
        "9fa7a82af850497441dd366fcead7b2987662b9b59e6b54ee306d306bc72e7ec",
        "6546aa34711a8fb1f334b132526f4b6100acc4c753ac7e297862a43e25d51bda",
    }
    assert all(
        record["full_exact_relative_interior_witness"] for record in crossings
    )


@pytest.mark.skipif(
    not DISC11.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_exact_ambient_carrier_retains_both_d35_subdivision_branches() -> None:
    certificate = certify_d35_carrier_regression(DISC11)

    assert certificate["verified"] is True
    assert certificate["distinct_opposite_carrier_targets"] == 2
    assert [record["representatives"] for record in certificate["galleries"]] == [
        [35, 162, 137, 122],
        [35, 162, 137],
    ]
    assert all(
        record["all_states_retained"] for record in certificate["galleries"]
    )
    assert [
        [crossing["ambient_wall"] for crossing in record["ambient_crossings"]]
        for record in certificate["galleries"]
    ] == [
        [[26, 0], [22, 0], [24, 0]],
        [[26, 0], [22, 0]],
    ]
