from __future__ import annotations

import json

from scripts.certify_alexeev_a14_8_1_c01_h9_14 import (
    INEQUALITY,
    OUTPUT,
    RAYS,
    build,
)


def test_frozen_h9_14_archive_rebuilds_exactly() -> None:
    payload, inequality, rays = build()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert INEQUALITY.read_bytes() == inequality
    assert RAYS.read_bytes() == rays

    decision = payload["decision"]
    assert decision["classification"] == (
        "finite_volume_noncompact_hyperbolic_coxeter_polytope"
    )
    assert decision["hyperbolic_dimension"] == 9
    assert decision["facet_count"] == 14
    assert decision["finite_vertex_count"] == 112
    assert decision["ideal_vertex_count"] == 15
    assert decision["hyperideal_vertex_count"] == 0

    orientation = payload["time_orientation"]
    assert orientation == {
        "all_extreme_rays_same_time_cone": True,
        "minimum_strict_interior_pairing": "21",
        "maximum_strict_interior_pairing": "226",
    }

    section = payload["section"]
    assert section["exact_reduction"]["active_hyperbolic_facet_count"] == 14
    assert section["coxeter_gate"]["all_pairs_coxeter"] is True
    assert section["restricted_mirrors"][
        "all_active_reflections_integral_on_saturated_fixed_lattice"
    ] is True

    vinberg = payload["literal_vinberg"]
    assert vinberg["literal_vinberg_passes"] is True
    assert vinberg["elliptic_edge_count"] == 707
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
