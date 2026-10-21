import json
from pathlib import Path

import pytest

from scripts.certify_rational_p18_p19_mmp import build_certificate


@pytest.fixture(scope="module")
def certificate():
    return build_certificate()


def test_p18_p19_surface_mmp_canonical_markings_are_exact_diagram_orbits(certificate):
    frozen = json.loads(
        Path("data/enumerations/rational_p18_p19_mmp.json").read_text()
    )
    assert certificate == frozen

    expected = {
        "18": {
            "facets": 37,
            "endpoints": 24,
            "automorphisms": 24,
            "state_counts": [
                1,
                15,
                6,
                12,
                12,
                24,
                24,
                24,
                24,
                48,
                24,
                24,
                24,
                24,
                24,
                24,
                24,
                48,
                24,
            ],
        },
        "19": {
            "facets": 50,
            "endpoints": 120,
            "automorphisms": 120,
            "state_counts": [
                1,
                5,
                75,
                30,
                60,
                60,
                120,
                120,
                120,
                120,
                240,
                120,
                120,
                120,
                120,
                120,
                120,
                120,
                240,
                120,
            ],
        },
    }
    for n, summary in expected.items():
        dimension = certificate["dimensions"][n]
        assert dimension["facet_count"] == summary["facets"]
        assert dimension["state_counts_by_descending_rank"] == summary["state_counts"]
        assert dimension["rank_two_endpoint_types"] == {
            "F1": summary["endpoints"],
            "F2": summary["endpoints"],
        }
        assert dimension["rank_one_P2_endpoint_count"] == summary["endpoints"]
        assert dimension["canonical_class_count"] == summary["endpoints"]
        assert (
            dimension["weighted_diagram_automorphism_order"]
            == summary["automorphisms"]
        )
        assert dimension["p2_and_f2_reconstructions_agree"]
        assert dimension["canonical_classes_equal_standard_diagram_orbit"]
        assert all(dimension["checks"].values())


def test_p18_dead_numerical_minus_one_branches_and_p19_first_orbit_are_frozen(
    certificate,
):
    p18 = certificate["dimensions"]["18"]
    p19 = certificate["dimensions"]["19"]

    assert p18["initial_square_minus_one_facet_count"] == 15
    assert p18["first_contraction_endpoint_survival_counts"] == {
        "18": 8,
        "25": 8,
        "26": 8,
    }
    assert p18["contraction_level_census"][1][
        "states_without_square_minus_one_image_count"
    ] == 12

    assert p19["initial_square_minus_one_facet_count"] == 5
    assert p19["first_contraction_endpoint_survival_counts"] == {
        "19": 24,
        "27": 24,
        "31": 24,
        "32": 24,
        "35": 24,
    }
    assert p19["contraction_level_census"][1][
        "states_without_square_minus_one_image_count"
    ] == 0
