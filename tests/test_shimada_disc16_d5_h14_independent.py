from __future__ import annotations

from collections import Counter
import json

from scripts.independent_certify_shimada_disc16_d5_h14 import (
    SUMMARY_OUTPUT as CONE_CERTIFICATE,
    load_ray_archive,
)
from scripts.independent_screen_shimada_disc16_d5_h14 import (
    OUTPUT as SCREEN_CERTIFICATE,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_primary_d5_reconstruction_and_intrinsic_representative() -> None:
    screen = _load(SCREEN_CERTIFICATE)

    assert screen["source"]["sha256"] == (
        "a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6"
    )
    assert screen["ambient_lattice"]["rank"] == 20
    assert screen["ambient_lattice"]["determinant"] == -16
    assert screen["ambient_lattice"]["inertia"] == [1, 19, 0]
    assert screen["chamber"]["index"] == 5
    assert screen["chamber"]["published_stabilizer_order"] == 20
    assert screen["chamber"]["published_wall_count"] == 165
    assert screen["chamber"]["ample_square"] == "35/2"
    assert screen["representative"]["selection"] == (
        "unique element of the full archived D5 stabilizer with fixed rank 15"
    )
    assert screen["representative"]["order"] == 2
    assert screen["representative"]["trace"] == 10
    assert screen["representative"]["fixed_rank"] == 15
    assert screen["representative"]["hyperbolic_dimension"] == 14


def test_exact_h14_facets_and_complete_pair_gate() -> None:
    screen = _load(SCREEN_CERTIFICATE)

    assert screen["wall_orbit_count"] == 100
    assert screen["wall_orbit_size_histogram"] == {"1": 35, "2": 65}
    assert screen["orbit_gram_inertia"] == [1, 14, 85]
    assert screen["facet_count"] == 60
    assert screen["redundant_orbit_count"] == 40
    assert screen["facet_normal_square_histogram"] == {
        "-1/2": 20,
        "-2": 30,
        "-4": 10,
    }
    assert screen["facet_pair_count"] == 1_770
    assert screen["coxeter_edge_counts"] == {
        "2": 865,
        "3": 35,
        "4": 70,
        "infinity": 310,
        "ultraparallel": 490,
    }
    assert screen["non_coxeter_pair_count"] == 0
    assert screen["status"]["all_facet_pairs_are_coxeter"] is True
    assert screen["status"]["finite_volume_claimed"] is False


def test_independent_exact_cone_proves_finite_volume_noncompact() -> None:
    cone = _load(CONE_CERTIFICATE)
    exact = cone["exact_double_description"]

    assert cone["cone"]["hyperbolic_dimension"] == 14
    assert cone["cone"]["linear_dimension"] == 15
    assert cone["cone"]["facet_count"] == 60
    assert cone["cone"]["inequality_rank"] == 15
    assert cone["cone"]["lineality_dimension"] == 0
    assert exact["extreme_ray_count"] == 28_872
    assert exact["finite_vertex_count"] == 27_845
    assert exact["ideal_vertex_count"] == 1_027
    assert exact["hyperideal_vertex_count"] == 0
    assert exact["finite_vertex_active_facet_count_histogram"] == {
        "14": 27_845
    }
    assert exact["ideal_vertex_active_facet_count_histogram"] == {
        "15": 432,
        "16": 540,
        "17": 55,
    }
    assert cone["decision"]["finite_volume"] is True
    assert cone["decision"]["compact"] is False
    assert cone["decision"]["classification"] == "finite_volume_noncompact"
    assert cone["checks"]["no_hyperideal_extreme_ray"] is True

    header, rays = load_ray_archive()
    assert header["schema"] == (
        "shimada_disc16_d5_h14_independent_exact_cone_rays_v1"
    )
    assert header["ray_count"] == len(rays) == 28_872
    assert Counter(len(active) for _, active in rays) == {
        14: 27_845,
        15: 432,
        16: 540,
        17: 55,
    }
    assert len({active for _, active in rays}) == 28_872
