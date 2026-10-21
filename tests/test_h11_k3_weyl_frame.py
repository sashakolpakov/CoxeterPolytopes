from pathlib import Path

import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.fold_h11_seed_via_shimada import fold_h11_seed
from scripts.frame_h11_vinberg_prefix_on_k3 import frame_vinberg_prefix


DISC11 = Path("/tmp/shimada_threediscs/threediscs/disc11.txt")
EXTENSION21 = Path("data/candidates/h11_vinberg_extension_bound21.json")


@pytest.mark.skipif(not DISC11.exists(), reason="external Shimada disc11 atlas absent")
def test_resumed_vinberg_prefix_has_an_exact_k3_weyl_frame() -> None:
    result = frame_vinberg_prefix(DISC11, extension_path=EXTENSION21)

    assert result["verified"] is True
    assert result["archived_simple_walls"] == 97
    assert result["distinct_transported_curve_wall_rays"] == 97
    assert result["weyl_field_reflections"] == 55
    assert result["auxiliary_shimada_crossings"] == 5
    assert result["terminal_physical_representative"] == 181
    assert result["independent_local_wall_rays"] == 11
    assert result["independent_local_rays_matched_in_prefix"] == 11
    assert result["prefix_hyperideal_witness"]["norm"] == {
        "a": "81/28",
        "b": "-25/14",
    }
    assert (
        result["prefix_hyperideal_witness"][
            "wall_pairings_have_one_consistent_orientation"
        ]
        is True
    )


@pytest.mark.skipif(not DISC11.exists(), reason="external Shimada disc11 atlas absent")
def test_h11_pell_seed_folds_to_the_known_enhanced_k3_involution() -> None:
    odd = fold_h11_seed(DISC11, seed_power=1)
    even = fold_h11_seed(DISC11, seed_power=2)

    assert odd["terminal_same_chamber"] is True
    assert odd["root_reflections_used"] == 1
    assert odd["infinite_order_trace_power"] is None
    odd_matrix = tuple(
        tuple(Golden.from_json(entry) for entry in row)
        for row in odd["folded_row_action"]
    )
    expected = tuple(
        tuple(
            (1 - PHI)
            if (i, j) == (1, 0)
            else Golden(-1)
            if (i, j) == (1, 1)
            else Golden(int(i == j))
            for j in range(10)
        )
        for i in range(10)
    )
    assert odd_matrix == expected

    assert even["terminal_same_chamber"] is True
    assert even["root_reflections_used"] == 2
    assert even["infinite_order_trace_power"] is None
    even_matrix = tuple(
        tuple(Golden.from_json(entry) for entry in row)
        for row in even["folded_row_action"]
    )
    assert even_matrix == tuple(
        tuple(Golden(int(i == j)) for j in range(10)) for i in range(10)
    )
