from __future__ import annotations

import json

import pytest

from scripts.certify_silver_k3_1291_global_vinberg import (
    DIAGNOSTIC,
    EXPECTED_DIAGNOSTIC_SHA256,
    EXPECTED_PAIRING_HISTOGRAM_SHA256,
    verify,
)


def test_frozen_1291_wall_global_vinberg_diagnostic() -> None:
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))

    assert verify(document) == {
        "diagnostic_sha256": EXPECTED_DIAGNOSTIC_SHA256,
        "pairing_histogram_sha256": EXPECTED_PAIRING_HISTOGRAM_SHA256,
        "facet_count": 1291,
        "elliptic_edge_count": 43891,
        "finite_vertex_count": 6091,
        "ideal_vertex_count": 0,
        "finite_plus_ideal_endpoint_histogram": {0: 4065, 1: 24833, 2: 14993},
        "finite_volume_open_edge_count": 28898,
        "connected_affine_subdiagram_count": 0,
        "full_neighborhood_audit_complete": True,
        "invalid_finite_angle_pair_count": 0,
        "compact_stopping_available": False,
        "finite_volume_stopping_available": False,
        "final_chamber_infinite_volume_proved": False,
        "other_reflective_norms_exhausted": False,
    }


def test_frozen_1291_wall_global_diagnostic_rejects_mutation() -> None:
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    document["finite_volume_open_edge_count"] = 28897

    with pytest.raises(ValueError, match="payload changed"):
        verify(document)
