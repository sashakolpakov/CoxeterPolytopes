from __future__ import annotations

import json

import pytest

from scripts.certify_silver_k3_1076_global_vinberg import (
    DIAGNOSTIC,
    EXPECTED_DIAGNOSTIC_SHA256,
    EXPECTED_PAIRING_HISTOGRAM_SHA256,
    verify,
)


def test_frozen_1076_wall_global_vinberg_diagnostic() -> None:
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))

    assert verify(document) == {
        "diagnostic_sha256": EXPECTED_DIAGNOSTIC_SHA256,
        "pairing_histogram_sha256": EXPECTED_PAIRING_HISTOGRAM_SHA256,
        "facet_count": 1076,
        "elliptic_edge_count": 38564,
        "finite_vertex_count": 5352,
        "ideal_vertex_count": 0,
        "finite_plus_ideal_endpoint_histogram": {0: 3554, 1: 21852, 2: 13158},
        "finite_volume_open_edge_count": 25406,
        "connected_affine_subdiagram_count": 0,
        "full_neighborhood_audit_complete": True,
        "invalid_finite_angle_pair_count": 0,
        "compact_stopping_available": False,
        "finite_volume_stopping_available": False,
        "final_chamber_infinite_volume_proved": False,
        "other_reflective_norms_exhausted": False,
    }


def test_frozen_1076_wall_global_diagnostic_rejects_mutation() -> None:
    document = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    document["finite_volume_open_edge_count"] = 25405

    with pytest.raises(ValueError, match="payload changed"):
        verify(document)
