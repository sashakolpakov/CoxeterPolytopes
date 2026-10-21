from __future__ import annotations

import json

from scripts.certify_silver_k3_726_global_vinberg import DIAGNOSTIC, verify


def test_frozen_726_wall_finite_volume_prefix_diagnostic() -> None:
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))

    assert verify(document) == {
        "facet_count": 726,
        "elliptic_edge_count": 25_764,
        "finite_vertex_count": 3_604,
        "ideal_vertex_count": 0,
        "finite_plus_ideal_endpoint_histogram": {
            0: 2_332,
            1: 14_428,
            2: 9_004,
        },
        "finite_volume_open_edge_count": 16_760,
        "connected_affine_subdiagram_count": 0,
        "full_neighborhood_audit_complete": True,
        "compact_stopping_available": False,
        "finite_volume_stopping_available": False,
        "final_chamber_infinite_volume_proved": False,
        "other_reflective_norms_exhausted": False,
    }
